#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <time.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cusparse_v2.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Prototypes for BLAS and LAPACK functions.                                                                              ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern "C"
{
  // Lapack routines.
  double dlamch_(const char   &cmach);

  void   dsyev_ (const char   &jobz,
                 const char   &uplo,
                 const int    &n,
                 double       *a,
                 const int    &lda,
                 double       *w,
                 double       *work,
                 int          &lwork,
                 int          &info);

  void   dsyevx_(const char   &jobz,
                 const char   &range,
                 const char   &uplo,
                 const int    &n,
                 double       *a,
                 const int    &lda,
                 const double &vl,
                 const double &vu,
                 const int    &il,
                 const int    &iu,
                 const double &abstol,
                 int          &m,
                 double       *w,
                 double       *z,
                 const int    &ldz,
                 double       *work,
                 const int    &lwork,
                 int          *iwork,
                 int          *ifail,
                 int          &info);

  // Prototype for private function in "gugatime.c".
  void gettime(double *tuser, double *tsys, double *twall);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Device functions.                                                                                                      ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__global__ void gpuInitBasisVectors(double    *B_d,
                                    int        ldb,
                                    const int *jrefconf_d,
                                    int        ni,
                                    int        nj)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < ni  &&  j < nj)
  {
    if (i == jrefconf_d[j] - 1)  B_d[ldb * j + i] = 1.0;
    else                         B_d[ldb * j + i] = 0.0;
  }
}


__global__ void gpuAddBasisVectors(double *B_d,
                                   int     ldb,
                                   double *Q_d,
                                   int     ldq,
                                   double *Hii_d,
                                   double *eig_d,
                                   int    *iroot_d,
                                   double  qtol,
                                   int     ni,
                                   int     nj)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < ni  &&  j < nj)
  {
    int    k         = iroot_d[j];
    double x         = eig_d[k] - Hii_d[i];

    // If x is too small, replace it by
    // +qtol or -qtol depending on its sign.
    //
    if (fabs(x) < qtol)
      x = copysign(qtol, x);

    B_d[ldb * j + i] = Q_d[ldq * k + i] / x;
  }
}


__global__ void gpuMultByDiag1(double       *A_d,
                               int           lda,
                               const double *b_d,
                               int           ni,
                               int           nj)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < ni  &&  j < nj)
  {
    A_d[lda * j + i] *= b_d[j];
  }
}


__global__ void gpuMultByDiag2(double       *C_d,
                               int           ldc,
                               const double *A_d,
                               int           lda,
                               const double *b_d,
                               int           ni,
                               int           nj)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < ni  &&  j < nj)
  {
    C_d[ldc * j + i] = A_d[lda * j + i] * b_d[j];
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Error handling functions.                                                                                              ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline void Error(cudaError_t err, const char *str)
{
  if (err != cudaSuccess)
  {
    printf(" %s(): %s\n", str, cudaGetErrorString(err));
    exit(1);
  }
}


inline void Error(cublasStatus_t err, const char *str)
{
  if (err != CUBLAS_STATUS_SUCCESS)
  {
    printf(" %s() returned error code %d.\n", str, err);
    exit(1);
  }
}


inline void Error(cusparseStatus_t err, const char *str)
{
  if (err != CUSPARSE_STATUS_SUCCESS)
  {
    printf(" %s() returned error code %d.\n", str, err);
    exit(1);
  }
}


inline void Sync(void)
{
  Error(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Class definitions.                                                                                                     ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Namespace is necessary because gpudavliu.cu
// and cpudavliu.cpp share the same names.

namespace gpu
{

// Forward declarations.
struct DeviceMatBase;
struct DeviceTrans;
struct ColumnVec;
struct ScaledCol;
struct ProdMatVec;
struct ProdTraVec;
struct DeviceCSR;
struct HostMatBase;


struct BlaHandleWrapper
{
  cublasHandle_t handle;

  BlaHandleWrapper(void)
  {
    Error(cublasCreate(&handle), "cublasCreate");
  }

 ~BlaHandleWrapper(void)
  {
    Error(cublasDestroy(handle), "cublasDestroy");
  }

  operator cublasHandle_t (void)
  {
    return handle;
  }
};


struct SpaHandleWrapper
{
  cusparseHandle_t handle;

  SpaHandleWrapper(void)
  {
    Error(cusparseCreate(&handle), "cusparseCreate");
  }

 ~SpaHandleWrapper(void)
  {
    Error(cusparseDestroy(handle), "cusparseDestroy");
  }

  operator cusparseHandle_t (void)
  {
    return handle;
  }
};


struct MatDescrWrapper
{
  cusparseMatDescr_t descrA;

  MatDescrWrapper(cusparseMatrixType_t matrixType,
                  cusparseFillMode_t   fillMode,
                  cusparseDiagType_t   diagType,
                  cusparseIndexBase_t  indexBase)
  {
    Error(cusparseCreateMatDescr(&descrA),             "cusparseCreateMatDescr");
    Error(cusparseSetMatType     (descrA, matrixType), "cusparseSetMatType");
    Error(cusparseSetMatFillMode (descrA, fillMode),   "cusparseSetMatFillMode");
    Error(cusparseSetMatDiagType (descrA, diagType),   "cusparseSetMatDiagType");
    Error(cusparseSetMatIndexBase(descrA, indexBase),  "cusparseSetMatIndexBase");
  }


 ~MatDescrWrapper(void)
  {
    Error(cusparseDestroyMatDescr(descrA), "cusparseDestroyMatDescr");
  }


  operator cusparseMatDescr_t (void) const
  {
    return descrA;
  }
};


struct Counters
{
  long isync;
  long imemcpy;
  long idcopy;
  long idnrm2s;  // cublasDnrm2(), ein Aufruf vor Sync().
  long idnrm2p;  // cublasDnrm2(), mehrere Aufrufe vor Sync().
  long idscal;
  long idgemv;
  long idgemm;
  long idcsrmm;
  long ieigen;
  long iother;


  void Init(void)
  {
    isync   = 0L;
    imemcpy = 0L;
    idcopy  = 0L;
    idnrm2s = 0L;
    idnrm2p = 0L;
    idscal  = 0L;
    idgemv  = 0L;
    idgemm  = 0L;
    idcsrmm = 0L;
    ieigen  = 0L;
    iother  = 0L;
  }


  long Total(void)
  {
    return isync + imemcpy + idcopy + idnrm2s + idnrm2p + idscal + idgemv + idgemm + idcsrmm + ieigen + iother;
  };
};


struct Timing
{
  static Counters ncalls;
  static Counters twall;

  static void Init(void)
  {
    ncalls.Init();
    twall.Init();
  }


  static long GetTime(clockid_t clk_id)
  {
    timespec ts;

    if (clock_gettime(clk_id, &ts))
    {
      printf("clock_gettime(): %s\n", strerror(errno));
      exit(1);
    }

    return 1000000000L * ts.tv_sec + ts.tv_nsec;
  }


  static void Print(void)
  {
    long MiB = (ncalls.imemcpy + 524288L) / 1048576L;
    printf("\n Wall clock time statistics of function calls:\n");
    printf(" init:  %6ld calls,%12.4le s\n", ncalls.isync,   1e-9 * twall.isync);
    printf(" memcpy:%6ld MiB,  %12.4le s\n", MiB,            1e-9 * twall.imemcpy);
    printf(" dcopy: %6ld calls,%12.4le s\n", ncalls.idcopy,  1e-9 * twall.idcopy);
    printf(" dnrm2s:%6ld calls,%12.4le s\n", ncalls.idnrm2s, 1e-9 * twall.idnrm2s);
    printf(" dnrm2p:%6ld calls,%12.4le s\n", ncalls.idnrm2p, 1e-9 * twall.idnrm2p);
    printf(" dscal: %6ld calls,%12.4le s\n", ncalls.idscal,  1e-9 * twall.idscal);
    printf(" dgemv: %6ld calls,%12.4le s\n", ncalls.idgemv,  1e-9 * twall.idgemv);
    printf(" dgemm: %6ld calls,%12.4le s\n", ncalls.idgemm,  1e-9 * twall.idgemm);
    printf(" dcsrmm:%6ld calls,%12.4le s\n", ncalls.idcsrmm, 1e-9 * twall.idcsrmm);
    printf(" eigen: %6ld calls,%12.4le s\n", ncalls.ieigen,  1e-9 * twall.ieigen);
    printf(" other: %6ld calls,%12.4le s\n", ncalls.iother,  1e-9 * twall.iother);
    printf(" Total:            %14.4le s\n",                 1e-9 * twall.Total());
  }
};


template <class T> struct DeviceVecBase
{
  T      *v;
  size_t  n;

  DeviceVecBase(size_t m, T *b_h):  v(b_h), n(m)  {}
  DeviceVecBase(void)  {}
 ~DeviceVecBase(void)  {}


  operator T * (void) const
  {
    return v;
  }


  T &operator [] (int i) const
  {
    return v[i];
  }


  const DeviceVecBase &mult(const DeviceMatBase &A, const double *b, double alpha, double beta) const;
  const DeviceVecBase &mult(const DeviceTrans   &A, const double *b, double alpha, double beta) const;
  const DeviceVecBase &operator  = (const ProdMatVec &prod) const;
  const DeviceVecBase &operator  = (const ProdTraVec &prod) const;
  const DeviceVecBase &operator += (const ProdMatVec &prod) const;
  const DeviceVecBase &operator += (const ProdTraVec &prod) const;
  const DeviceVecBase &operator -= (const ProdMatVec &prod) const;
  const DeviceVecBase &operator -= (const ProdTraVec &prod) const;


  void copyFromHost(const T *b, size_t n) const
  {
    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cudaMemcpy(v, b, n * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy");
    Timing::twall.imemcpy  += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.imemcpy += n * (long)sizeof(T);
  }
};


template <class T> struct DeviceVecMem: public DeviceVecBase<T>
{
  DeviceVecMem(size_t m, int prtlevel)
  {
    T           *tp;
    size_t       size  = m * sizeof(T);
    cudaError_t  err   = cudaMalloc(&tp, size);
    DeviceVecBase<T>::v = tp;
    DeviceVecBase<T>::n = m;

    if (err != cudaSuccess)
    {
      printf(" Failed to allocate %ld bytes of global memory on CUDA device.\n", size);
      printf("cudaMalloc(): %s\n", cudaGetErrorString(err));
      exit(1);
    }

    if (prtlevel >= 2)
      printf(" Allocated %10ld bytes of GPU memory.\n", size);
  }


 ~DeviceVecMem(void)
  {
    Error(cudaFree(DeviceVecBase<T>::v), "cudaFree");
  }
};


template <class T> struct DeviceVec: public DeviceVecMem<T>
{
  DeviceVec(size_t m,               int prtlevel):  DeviceVecMem<T>(m, prtlevel)  {}


  DeviceVec(size_t m, const T *b_h, int prtlevel):  DeviceVecMem<T>(m, prtlevel)
  {
    this->copyFromHost(b_h, m);
  }


 ~DeviceVec(void)  {}
};


template <class T> struct HostVec
{
  T *v;

  HostVec(size_t n)
  {
    if (!(v = new T[n]))
    {
      printf(" Failed to allocate %ld bytes of CPU memory.\n", n * sizeof(T));
      exit(1);
    }
  }


 ~HostVec(void)
  {
    delete [] v;
  }


  operator T * (void) const
  {
    return v;
  }


  T &operator [] (int i) const
  {
    return v[i];
  }
};


struct ColumnVec
{
  const DeviceMatBase &M;
  int                  j;

  ColumnVec(const DeviceMatBase &M, int j);
 ~ColumnVec(void)  {}

  operator double * (void) const;

  const ColumnVec &mult(const DeviceMatBase &A, const double *b, double alpha, double beta) const;
  const ColumnVec &operator  = (const ScaledCol  &scol) const;
  const ColumnVec &operator  = (const ProdMatVec &prod) const;
  const ColumnVec &operator += (const ProdMatVec &prod) const;
  const ColumnVec &operator -= (const ProdMatVec &prod) const;
  const ColumnVec &operator *= (double b) const;
  const ScaledCol  operator *  (double b) const;
};


struct ScaledCol
{
  const ColumnVec &a;
  double           b;

  ScaledCol(const ColumnVec &a, double b):  a(a), b(b)  {}
 ~ScaledCol(void)  {}
};


struct ProdMatVec
{
  const DeviceMatBase &A;
  const double        *b;

  ProdMatVec(const DeviceMatBase &A, const double *b):  A(A), b(b)  {}
 ~ProdMatVec(void)  {}
};


struct ProdTraVec
{
  const DeviceTrans &A;
  const double      *b;

  ProdTraVec(const DeviceTrans &A, const double *b):  A(A), b(b)  {}
 ~ProdTraVec(void)  {}
};


struct ProdMatMat
{
  const DeviceMatBase &A;
  const DeviceMatBase &B;

  ProdMatMat(const DeviceMatBase &A, const DeviceMatBase &B):  A(A), B(B)  {}
 ~ProdMatMat(void)  {}
};


struct ProdTraMat
{
  const DeviceTrans   &A;
  const DeviceMatBase &B;

  ProdTraMat(const DeviceTrans &A, const DeviceMatBase &B):  A(A), B(B)  {}
 ~ProdTraMat(void)  {}
};


struct ProdSpaDen
{
  const DeviceCSR     &A;
  const DeviceMatBase &B;

  ProdSpaDen(const DeviceCSR &A, const DeviceMatBase &B):  A(A), B(B)  {}
 ~ProdSpaDen(void)  {}
};


struct DeviceTrans
{
  cublasHandle_t  handle;
  double         *M;
  int             ld;
  int             ni;  // Number of rows of the transposed matrix.
  int             nj;  // Number of columns of the transposed matrix.

  DeviceTrans(cublasHandle_t h, double *M, int ld, int mi, int mj):  handle(h), M(M), ld(ld), ni(mi), nj(mj)  {}

  operator double * (void) const
  {
    return M;
  }

  ProdTraVec operator * (const ColumnVec             &b) const;
  ProdTraVec operator * (const DeviceVecBase<double> &b) const;
  ProdTraMat operator * (const DeviceMatBase         &B) const;
};


struct DeviceSym
{
  double *M;
  int     ld;
  int     ni;  // Number of rows (and columns) of the symmetric matrix.

  DeviceSym(double *M, int ld, int mi):  M(M), ld(ld), ni(mi)  {}
  DeviceSym(int mi):                     ni(mi)                {}

  const DeviceSym &operator = (const ProdTraMat &prod) const;
};


struct DeviceCSR
{
  cusparseHandle_t  handle;
  int               nrows;
  int               nnz;
  DeviceVec<double> a;
  DeviceVec<int>    ia;
  DeviceVec<int>    ja;

  DeviceCSR(cusparseHandle_t h, int n, const double *Hij, const int *ifirst, const int *icol, int prtlevel):
    handle(h), nrows(n), nnz(ifirst[n]-1), a(nnz, Hij, prtlevel), ia(n+1, ifirst, prtlevel), ja(nnz, icol, prtlevel)  {}

 ~DeviceCSR(void)  {}

  ProdSpaDen operator * (const DeviceMatBase &B) const
  {
    return ProdSpaDen(*this, B);
  }
};


struct DeviceMatBase
{
  cublasHandle_t  handle;
  double         *M;
  int             ld;
  int             ni;
  int             nj;

  DeviceMatBase(cublasHandle_t h, double *M, int ld, int mi, int mj):  handle(h), M(M), ld(ld), ni(mi), nj(mj)  {}
  DeviceMatBase(cublasHandle_t h, int mi, int mj):                     handle(h), ni(mi), nj(mj)                {}


  operator double * (void) const
  {
    return M;
  }


  ColumnVec colvec(int j) const
  {
    return ColumnVec(*this, j);
  }


  DeviceMatBase submat(int mi, int mj) const
  {
    if (mi > ni  ||  mj > nj)
    {
      printf(" Submatrix too large.\n");
      exit(1);
    }

    return DeviceMatBase(handle, M, ld, mi, mj);
  }


  DeviceMatBase submat(int mi, int mj, int i, int j) const
  {
    if (mi + i > ni  ||  mj + j > nj)
    {
      printf(" Submatrix too large.\n");
      exit(1);
    }

    return DeviceMatBase(handle, M+ld*j+i, ld, mi, mj);
  }


  DeviceSym subsym(int mi) const
  {
    if (mi > ni  ||  mi > nj)
    {
      printf(" Submatrix too large.\n");
      exit(1);
    }

    return DeviceSym(M, ld, mi);
  }


  DeviceSym subsym(int mi, int i) const
  {
    if (mi + i > ni  ||  mi + i > nj)
    {
      printf(" Submatrix too large.\n");
      exit(1);
    }

    return DeviceSym(M+(ld+1)*i, ld, mi);
  }


  double &operator () (int i, int j) const
  {
    if (i >= ni  ||  j >= nj)
    {
      printf(" Subscript out of range.\n");
      exit(1);
    }

    return M[ld*j+i];
  }


  DeviceTrans operator ~ (void) const
  {
    return DeviceTrans(handle, M, ld, nj, ni);
  }


  ProdMatVec operator * (const ColumnVec &b) const
  {
    if (nj != b.M.ni)
    {
      printf(" Incompatible matrix dimensions during M * v.\n");
      printf("   M: %6d x %6d\n", ni, nj);
      printf("   v:       %6d\n", b.M.ni);
      exit(1);
    }

    return ProdMatVec(*this, b);
  }


  ProdMatVec operator * (const DeviceVecBase<double> &b) const
  {
    return ProdMatVec(*this, b);
  }


  ProdMatMat operator * (const DeviceMatBase &B) const
  {
    return ProdMatMat(*this, B);
  }


  const DeviceMatBase &mult(const DeviceMatBase &A, const DeviceMatBase &B, double alpha, double beta) const
  {
    if (ni != A.ni  ||  nj != B.nj  ||  A.nj != B.ni)
    {
      printf(" Incompatible matrix dimensions during C = A * B.\n");
      printf("   A: %6d x %6d\n", A.ni, A.nj);
      printf("   B: %6d x %6d\n", B.ni, B.nj);
      printf("   C: %6d x %6d\n", ni,   nj);
      exit(1);
    }

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, ni, nj, A.nj,
                      &alpha, A, A.ld, B, B.ld, &beta, M, ld), "cublasDgemm");
    Sync();
    Timing::twall.idgemm += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.idgemm++;
    return *this;
  }


  const DeviceMatBase &mult(const DeviceTrans &A, const DeviceMatBase &B, double alpha, double beta) const
  {
    if (ni != A.ni  ||  nj != B.nj  ||  A.nj != B.ni)
    {
      printf(" Incompatible matrix dimensions during C = A**T * B.\n");
      printf("   A**T: %6d x %6d\n", A.ni, A.nj);
      printf("   B:    %6d x %6d\n", B.ni, B.nj);
      printf("   C:    %6d x %6d\n", ni,   nj);
      exit(1);
    }

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, ni, nj, A.nj,
                      &alpha, A, A.ld, B, B.ld, &beta, M, ld), "cublasDgemm");
    Sync();
    Timing::twall.idgemm += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.idgemm++;
    return *this;
  }


  const DeviceMatBase &operator = (const ProdMatMat &prod) const
  {
    return mult(prod.A, prod.B, 1.0, 0.0);
  }


  const DeviceMatBase &operator = (const ProdTraMat &prod) const
  {
    return mult(prod.A, prod.B, 1.0, 0.0);
  }


  const DeviceMatBase &operator = (const ProdSpaDen &prod) const
  {
    if (ni != prod.A.nrows  ||  nj != prod.B.nj  ||  prod.A.nrows != prod.B.ni)
    {
      printf(" Incompatible matrix dimensions during C = CSR * B.\n");
      printf("   CSR: %6d x %6d\n", prod.A.nrows, prod.A.nrows);
      printf("   B:   %6d x %6d\n", prod.B.ni,    prod.B.nj);
      printf("   C:   %6d x %6d\n", ni,           nj);
      exit(1);
    }

    MatDescrWrapper descrA(CUSPARSE_MATRIX_TYPE_SYMMETRIC,
                           CUSPARSE_FILL_MODE_UPPER,
                           CUSPARSE_DIAG_TYPE_NON_UNIT,
                           CUSPARSE_INDEX_BASE_ONE);
    double alpha = 1.0;
    double beta  = 0.0;

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cusparseDcsrmm(prod.A.handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                         prod.A.nrows, prod.B.nj, prod.A.nrows, prod.A.nnz,
                         &alpha, descrA, prod.A.a, prod.A.ia, prod.A.ja,
                         prod.B, prod.B.ld, &beta, *this, ld), "cusparseDcsrmm");
    Sync();
    Timing::twall.idcsrmm += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.idcsrmm++;
    return *this;
  }


  const DeviceMatBase &operator += (const ProdMatMat &prod) const
  {
    return mult(prod.A, prod.B, 1.0, 1.0);
  }


  const DeviceMatBase &operator += (const ProdTraMat &prod) const
  {
    return mult(prod.A, prod.B, 1.0, 1.0);
  }


  const DeviceMatBase &operator -= (const ProdMatMat &prod) const
  {
    return mult(prod.A, prod.B, -1.0, 1.0);
  }


  const DeviceMatBase &operator -= (const ProdTraMat &prod) const
  {
    return mult(prod.A, prod.B, -1.0, 1.0);
  }


  const DeviceMatBase &multByDiag(double *b_d) const
  {
    const dim3 blockDim(256, 1);
    const dim3 gridDim((ni + blockDim.x - 1) / blockDim.x, (nj + blockDim.y - 1) / blockDim.y);

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    gpuMultByDiag1<<<gridDim, blockDim>>>(M, ld, b_d, ni, nj);
    Sync();
    Timing::twall.iother += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.iother++;

    return *this;
  }


  const DeviceMatBase &multByDiag(const DeviceMatBase &A, double *b_d) const
  {
    const dim3 blockDim(256, 1);
    const dim3 gridDim((ni + blockDim.x - 1) / blockDim.x, (nj + blockDim.y - 1) / blockDim.y);

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    gpuMultByDiag2<<<gridDim, blockDim>>>(M, ld, A, A.ld, b_d, ni, nj);
    Sync();
    Timing::twall.iother += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.iother++;

    return *this;
  }


  void colnorms(double *b) const
  {
    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    for (int j=0; j<nj; ++j)
      Error(cublasDnrm2(handle, ni, M+ld*j, 1, b+j), "cublasDnrm2");
    Sync();
    Timing::twall.idnrm2p  += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.idnrm2p += nj;
  }


  void colnorm(int j, double &result) const
  {
    if (j >= nj)
    {
      printf(" Column index out of range.\n");
      exit(1);
    }

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cublasDnrm2(handle, ni, M+ld*j, 1, &result), "cublasDnrm2");
    Sync();
    Timing::twall.idnrm2s += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.idnrm2s++;
  }


  void swap(DeviceMatBase &B)
  {
    // Swap handle.
    cublasHandle_t h = B.handle;
    B.handle         = handle;
    handle           = h;

    // Swap pointer.
    double        *x = B.M;
    B.M              = M;
    M                = x;

    // Swap leading dimension.
    int            i = B.ld;
    B.ld             = ld;
    ld               = i;

    // Swap number of rows.
    i                = B.ni;
    B.ni             = ni;
    ni               = i;

    // Swap number of columns.
    i                = B.nj;
    B.nj             = nj;
    nj               = i;
  }


  void init(const int *jrefconf_d) const
  {
    const dim3 blockDim(256, 1);
    const dim3 gridDim((ni + blockDim.x - 1) / blockDim.x, (nj + blockDim.y - 1) / blockDim.y);

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    gpuInitBasisVectors<<<gridDim, blockDim>>>(M, ld, jrefconf_d, ni, nj);
    Sync();
    Timing::twall.iother += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.iother++;
  }


  void copyToHost(double *B, int ldb) const
  {
    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cudaMemcpy2D(B, ldb*sizeof(double),
                       M, ld *sizeof(double),
                       ni*sizeof(double), nj,
                       cudaMemcpyDeviceToHost), "cudaMemcpy2D");
    Timing::twall.imemcpy  += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.imemcpy += (long)ni * (long)nj * (long)sizeof(double);
  }


  void copyFromHost(const double *B, int ldb) const
  {
    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cudaMemcpy2D(M, ld *sizeof(double),
                       B, ldb*sizeof(double),
                       ni*sizeof(double), nj,
                       cudaMemcpyHostToDevice), "cudaMemcpy2D");
    Timing::twall.imemcpy  += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.imemcpy += (long)ni * (long)nj * (long)sizeof(double);
  }


  void copyFromHost(const HostMatBase &B) const;
  void diag(HostVec<double> &eig_h) const;
};


struct DeviceMat: public DeviceMatBase
{
  DeviceMat(cublasHandle_t h, int ni, int nj, int prtlevel):  DeviceMatBase(h, ni, nj)
  {
    ld               = ((ni + 15) / 16) * 16;
    size_t      size = ld * nj * sizeof(double);
    cudaError_t err  = cudaMalloc(&M, size);

    if (err != cudaSuccess)
    {
      printf(" Failed to allocate %ld bytes of global memory on CUDA device.\n", size);
      printf(" cudaMalloc(): %s\n", cudaGetErrorString(err));
      exit(1);
    }

    if (prtlevel >= 2)
      printf(" Allocated %10ld bytes of GPU memory.\n", size);
  }

 ~DeviceMat(void)
  {
    Error(cudaFree(M), "cudaFree");
  }
};


struct HostMatBase
{
  double *M;
  int     ld;
  int     ni;
  int     nj;

  HostMatBase(double *M, int ld, int mi, int mj):  M(M), ld(ld), ni(mi), nj(mj)  {}
  HostMatBase(int mi, int mj):                     ni(mi), nj(mj)                {}


  operator double * (void) const
  {
    return M;
  }


  double &operator () (int i, int j) const
  {
    return M[j*ld+i];
  }


  HostMatBase submat(int mi, int mj) const
  {
    if (mi > ni  ||  mj > nj)
    {
      printf(" Submatrix too large.\n");
      exit(1);
    }

    return HostMatBase(M, ld, mi, mj);
  }


  void copyFromDevice(DeviceMatBase &B) const
  {
    if (ni > B.ni  ||  nj > B.nj)
    {
      printf(" Incompatible matrix dimensions during copying.\n");
      exit(1);
    }

    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    Error(cudaMemcpy2D(M, ld  *sizeof(double),
                       B, B.ld*sizeof(double),
                       ni*sizeof(double), nj,
                       cudaMemcpyDeviceToHost), "cudaMemcpy2D");
    Timing::twall.imemcpy  += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.imemcpy += (long)ni * (long)nj * (long)sizeof(double);
  }


  void diag(const HostVec<double> &eig_h) const
  {
    if (ni != nj)
    {
      printf(" Square matrix expected.\n");
      exit(1);
    }

    int    lwork = -1;
    int    info  =  0;
    double work1;
    dsyev_('V', 'L', ni, M, ld, eig_h, &work1, lwork, info);

    if (info != 0)
    {
      printf(" Error from DSYEV() when querying workspace size.\n");
      exit(1);
    }

    lwork = work1;
    HostVec<double> work(lwork);
    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    dsyev_('V', 'L', ni, M, ld, eig_h, work, lwork, info);
    Timing::twall.ieigen += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.ieigen++;

    if (info != 0)
    {
      printf(" Error from DSYEV() when diagonalizing.\n");
      exit(1);
    }
  }


  void eigen(const HostMatBase &U_h, const HostVec<double> &eig_h) const
  {
    if (ni != nj)
    {
      printf(" Square matrix expected.\n");
      exit(1);
    }

    if (U_h.ni != ni  ||  U_h.nj > ni)
    {
      printf(" Incompatible matrix dimensions.\n");
      exit(1);
    }

    HostVec<int> iwork(ni * 5);
    HostVec<int> ifail(ni);
    double       abstol = 2.0 * dlamch_('S');
    double       work1;
    int          lwork  = -1;
    int          info   =  0;
    int          m;

    dsyevx_('V', 'I', 'L', ni, M, ld, 0.0, 0.0, 1, U_h.nj, abstol, m, eig_h, U_h, U_h.ld, &work1, lwork, iwork, ifail, info);

    if (info != 0)
    {
      printf(" Error from DSYEVX() when querying workspace size.\n");
      exit(1);
    }

    lwork = work1;
    HostVec<double> work(lwork);
    long t1 = Timing::GetTime(CLOCK_MONOTONIC);
    dsyevx_('V', 'I', 'L', ni, M, ld, 0.0, 0.0, 1, U_h.nj, abstol, m, eig_h, U_h, U_h.ld, work, lwork, iwork, ifail, info);
    Timing::twall.ieigen += Timing::GetTime(CLOCK_MONOTONIC) - t1;
    Timing::ncalls.ieigen++;

    if (info != 0  ||  m < U_h.nj)
    {
      printf(" Error from DSYEVX() when diagonalizing.\n");
      exit(1);
    }
  }


  void Print(void)
  {
    for (int k=0; k<nj; k+=6)
    {
      for (int j=k; j<k+6; ++j)
	printf("%20d", j+1);

      printf("\n");

      for (int i=k; i<ni; ++i)
      {
	int l = k+6;
	if (nj  < l)  l = nj;
        if (i+1 < l)  l = i+1;

	printf("%8d", i+1);

	for (int j=k; j<l; ++j)
	  printf("%20.15lf", (*this)(i,j));

	printf("\n");
      }

      printf("\n");
    }
  }
};


struct HostMat: public HostMatBase
{
  HostMat(int ni, int nj):  HostMatBase(ni, nj)
  {
    ld           = ni;
    size_t nelem = ld * nj;

    if (!(M = new double[nelem]))
    {
      printf(" Failed to allocate %ld bytes of CPU memory.\n", nelem * sizeof(double));
      exit(1);
    }
  }

 ~HostMat(void)
  {
    delete [] M;
  }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Member and friend functions.                                                                                           ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <> inline const DeviceVecBase<double> &
  DeviceVecBase<double>::mult(const DeviceMatBase &A, const double *b, double alpha, double beta) const
{
  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDgemv(A.handle, CUBLAS_OP_N, A.ni, A.nj, &alpha, A, A.ld, b, 1, &beta, *this, 1), "cublasDgemv");
  Sync();
  Timing::twall.idgemv += Timing::GetTime(CLOCK_MONOTONIC) - t1;
  Timing::ncalls.idgemv++;
  return *this;
}


template <> inline const DeviceVecBase<double> &
  DeviceVecBase<double>::mult(const DeviceTrans &A, const double *b, double alpha, double beta) const
{
  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDgemv(A.handle, CUBLAS_OP_T, A.nj, A.ni, &alpha, A, A.ld, b, 1, &beta, *this, 1), "cublasDgemv");
  Sync();
  Timing::twall.idgemv += Timing::GetTime(CLOCK_MONOTONIC) - t1;
  Timing::ncalls.idgemv++;
  return *this;
}


template <> inline const DeviceVecBase<double> &DeviceVecBase<double>::operator = (const ProdMatVec &prod) const
{
  return mult(prod.A, prod.b, 1.0, 0.0);
}


template <> inline const DeviceVecBase<double> &DeviceVecBase<double>::operator = (const ProdTraVec &prod) const
{
  return mult(prod.A, prod.b, 1.0, 0.0);
}


template <> inline const DeviceVecBase<double> &DeviceVecBase<double>::operator += (const ProdMatVec &prod) const
{
  return mult(prod.A, prod.b, 1.0, 1.0);
}


template <> inline const DeviceVecBase<double> &DeviceVecBase<double>::operator += (const ProdTraVec &prod) const
{
  return mult(prod.A, prod.b, 1.0, 1.0);
}


template <> inline const DeviceVecBase<double> &DeviceVecBase<double>::operator -= (const ProdMatVec &prod) const
{
  return mult(prod.A, prod.b, -1.0, 1.0);
}


template <> inline const DeviceVecBase<double> &DeviceVecBase<double>::operator -= (const ProdTraVec &prod) const
{
  return mult(prod.A, prod.b, -1.0, 1.0);
}


inline ColumnVec::ColumnVec(const DeviceMatBase &M, int j):  M(M), j(j)
{
  if (j >= M.nj)
  {
    printf(" Column index out of range.");
    exit(1);
  }
}


inline ColumnVec::operator double * (void) const
{
  return M.M + M.ld * j;
}


inline const ColumnVec &ColumnVec::mult(const DeviceMatBase &A, const double *b, double alpha, double beta) const
{
  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDgemv(A.handle, CUBLAS_OP_N, A.ni, A.nj, &alpha, A, A.ld, b, 1, &beta, *this, 1), "cublasDgemv");
  Sync();
  Timing::twall.idgemv += Timing::GetTime(CLOCK_MONOTONIC) - t1;
  Timing::ncalls.idgemv++;
  return *this;
}


inline const ColumnVec &ColumnVec::operator = (const ScaledCol &scol) const
{
  if (M.ni != scol.a.M.ni)
  {
    printf(" Incompatible matrix dimensions during y = x * a.");
    exit(1);
  }

  // Copy and scale in two steps, there seems to be no combined function.
  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDcopy(M.handle, M.ni, scol.a, 1, *this, 1), "cublasDcopy");
  Sync();
  long t2 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDscal(M.handle, M.ni, &scol.b, *this, 1), "cublasDaxpy");
  Sync();
  Timing::twall.idscal += Timing::GetTime(CLOCK_MONOTONIC) - t2;
  Timing::twall.idcopy += t2 - t1;
  Timing::ncalls.idcopy++;
  Timing::ncalls.idscal++;
  return *this;
}


inline const ColumnVec &ColumnVec::operator = (const ProdMatVec &prod) const
{
  if (M.ni != prod.A.ni)
  {
    printf(" Incompatible matrix dimensions during c = A * b.");
    exit(1);
  }

  return mult(prod.A, prod.b, 1.0, 0.0);
}


inline const ColumnVec &ColumnVec::operator += (const ProdMatVec &prod) const
{
  if (M.ni != prod.A.ni)
  {
    printf(" Incompatible matrix dimensions during c = A * b.");
    exit(1);
  }

  return mult(prod.A, prod.b, 1.0, 1.0);
}


inline const ColumnVec &ColumnVec::operator -= (const ProdMatVec &prod) const
{
  if (M.ni != prod.A.ni)
  {
    printf(" Incompatible matrix dimensions during c = A * b.");
    exit(1);
  }

  return mult(prod.A, prod.b, -1.0, 1.0);
}


inline const ColumnVec &ColumnVec::operator *= (double b) const
{
  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDscal(M.handle, M.ni, &b, M, 1), "cublasDscal");
  Sync();
  Timing::twall.idscal += Timing::GetTime(CLOCK_MONOTONIC) - t1;
  Timing::ncalls.idscal++;
  return *this;
}


const ScaledCol ColumnVec::operator * (double b) const
{
  return ScaledCol(*this, b);
}


inline ProdTraVec DeviceTrans::operator * (const ColumnVec &b) const
{
  if (nj != b.M.ni)
  {
    printf(" Incompatible matrix dimensions during M**T * v.\n");
    printf("   M**T: %6d x %6d\n", ni, nj);
    printf("   v:          %6d\n", b.M.ni);
    exit(1);
  }

  return ProdTraVec(*this, b);
}


inline ProdTraVec DeviceTrans::operator * (const DeviceVecBase<double> &b) const
{
  return ProdTraVec(*this, b);
}


inline ProdTraMat DeviceTrans::operator * (const DeviceMatBase &B) const
{
  return ProdTraMat(*this, B);
}


const DeviceSym &DeviceSym::operator = (const ProdTraMat &prod) const
{
  if (ni != prod.A.ni  ||  ni != prod.B.nj  ||  prod.A.nj != prod.B.ni)
  {
    printf(" Incompatible matrix dimensions during Sym = A**T * B.\n");
    printf("   A**T: %6d x %6d\n", prod.A.ni, prod.A.nj);
    printf("   B:    %6d x %6d\n", prod.B.ni, prod.B.nj);
    printf("   Sym:  %6d x %6d\n", ni,        ni);
    exit(1);
  }

  double alpha = 1.0;
  double beta  = 0.0;
  /*
  Error(cublasDtrmm(prod.B.handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
                    prod.B.ni, prod.B.nj, &alpha, prod.A, prod.A.ld, prod.B, prod.B.ld, M, ld), "cublasDtrmm");
  */
  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  Error(cublasDgemm(prod.B.handle, CUBLAS_OP_T, CUBLAS_OP_N, ni, ni, prod.B.ni,
                    &alpha, prod.A, prod.A.ld, prod.B, prod.B.ld, &beta, M, ld), "cublasDgemm");
  Sync();
  Timing::twall.idgemm += Timing::GetTime(CLOCK_MONOTONIC) - t1;
  Timing::ncalls.idgemm++;
  return *this;
}


inline void DeviceMatBase::copyFromHost(const HostMatBase &B) const
{
  if (ni > B.ni  ||  nj > B.nj)
  {
    printf(" Incompatible matrix dimensions during copying.\n");
    exit(1);
  }

  copyFromHost(B, B.ld);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Host functions.                                                                                                        ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void AddBasisVectors(const DeviceMatBase     &B_d,
                     const DeviceMatBase     &Q_d,
                     const DeviceVec<double> &Hii_d,
                     const DeviceVec<double> &eig_d,
                     const DeviceVec<int>    &iroot_d,
                     double                   qtol)
{
  // Effectively, the following operation is performed:
  //
  // for (int i=0; i<B_d.ni; ++i)
  //   for (int j=0; j<B_d.nj; ++j)
  //   {
  //     int k    = iroot_d[j];
  //     double x = eig_d[k] - Hii[i];
  //
  //     if (fabs(x) < qtol)
  //       x = copysign(qtol, x);
  //
  //     B_d(i,j) = Q_d(i,k) / x;
  //   }

  const dim3 blockDim(256, 1);
  const dim3 gridDim((B_d.ni + blockDim.x - 1) / blockDim.x, (B_d.nj + blockDim.y - 1) / blockDim.y);

  long t1 = Timing::GetTime(CLOCK_MONOTONIC);
  gpuAddBasisVectors<<<gridDim, blockDim>>>(B_d, B_d.ld, Q_d, Q_d.ld, Hii_d, eig_d, iroot_d, qtol, B_d.ni, B_d.nj);
  Sync();
  Timing::twall.iother += Timing::GetTime(CLOCK_MONOTONIC) - t1;
  Timing::ncalls.iother++;
}


void Schmidt1(const DeviceMatBase &B_d,
              const DeviceMatBase &Q_d,
              const DeviceMatBase &A_d)
{
  // Schmidt orthogonalize the column vectors of Q_d against
  // the column vectors of B_d. A_d is used as scratch space.

  A_d  = ~B_d * Q_d;
  Q_d -=  B_d * A_d;
}


void Schmidt2(const DeviceMatBase         &B_d,
              const DeviceMatBase         &Q_d,
              const DeviceVecBase<double> &v_d,
              int                         &newdav,
              double                       qtol)
{
  // Schmidt orthogonalize the column vectors of Q among each other.
  // Vectors with a norm less than qtol after projecting out the
  // previous vectors are dropped. The final vectors are copied
  // into B. Their number is returned in newdav.

  int nci = B_d.ni;

  if (nci != Q_d.ni  ||  B_d.nj < Q_d.nj)
  {
    printf(" Incompatible matrix dimensions in Schmidt2().\n");
    printf("   B: %6d x %6d\n", B_d.ni, B_d.nj);
    printf("   Q: %6d x %6d\n", Q_d.ni, Q_d.nj);
    exit(1);
  }

  newdav = 0;

  for (int j=0; j<Q_d.nj; ++j)
  {
    if (newdav > 0)
    {
      // Project previous vectors out of current (j-th) vector.
      v_d            = ~B_d.submat(nci, newdav) * Q_d.colvec(j);
      Q_d.colvec(j) -=  B_d.submat(nci, newdav) * v_d;
    }

    double scalf;
    Q_d.colnorm(j, scalf);

    if (scalf > qtol)
    {
      // Normalize current vector and copy into B.
      scalf                = 1.0 / scalf;
      B_d.colvec(newdav++) = Q_d.colvec(j) * scalf;
    }
  }
}


void Schmidt3(const DeviceMatBase         &Q_d,
              const DeviceVecBase<double> &v_d)
{
  for (int j=0; j<Q_d.nj; ++j)
  {
    if (j > 0)
    {
      // Project previous vectors out of current (j-th) vector.
      v_d            = ~Q_d.submat(Q_d.ni, j) * Q_d.colvec(j);
      Q_d.colvec(j) -=  Q_d.submat(Q_d.ni, j) * v_d;
    }

    // Normalize current vector.
    double scalf;
    Q_d.colnorm(j, scalf);
    scalf          = 1.0 / scalf;
    Q_d.colvec(j) *= scalf;
  }
}


void DavLiu(double *Hii,
            double *Hij,
            int    *ifirst,
            int    *icol,
            int     nci,
            int    *jrefconf,
            int     inidav,
            double *E,
            double *C,
            int     ldc,
            int     nroots,
            int     mindav,
            int     maxdav,
            int     kitdav,
            double  qtol,
            int     prtlevel,
            int    &icall)
{
  Timing::Init();
  long              ibegin = Timing::GetTime(CLOCK_MONOTONIC);
  Sync();
  Timing::twall.isync      = Timing::GetTime(CLOCK_MONOTONIC) - ibegin;
  Timing::ncalls.isync++;
  BlaHandleWrapper  bla;
  SpaHandleWrapper  spa;
  DeviceVec<double> Hii_d(nci, Hii, prtlevel);
  DeviceCSR         Ham_d(spa, nci, Hij, ifirst, icol, prtlevel);
  DeviceMat         B_d(bla, nci, maxdav, prtlevel);
  DeviceMat         P_d(bla, nci, maxdav, prtlevel);
  DeviceMat         G_d(bla, maxdav, maxdav, prtlevel);
  DeviceMat         A_d(bla, maxdav, maxdav, prtlevel);
  HostMat           A_h(maxdav, maxdav);
  DeviceMat         U_d(bla, maxdav, maxdav, prtlevel);
  HostMat           U_h(maxdav, mindav);
  DeviceVec<double> v_d(maxdav, prtlevel);
  DeviceVec<double> eig_d(maxdav, prtlevel);
  HostVec<double>   eig_h(maxdav);
  DeviceVec<double> qnorm_d(nroots, prtlevel);
  HostVec<double>   qnorm_h(nroots);
  DeviceVec<int>    iroot_d(nroots, prtlevel);
  HostVec<int>      iroot_h(nroots);
  DeviceVec<int>    jrefconf_d(inidav, jrefconf, prtlevel);
  double            twall1, twall2, twall3;
  int               numdav = inidav;  // Current number of basis vectors.

  if (prtlevel >= 2)
  {
    // Print greeting message.
    printf("\n GPU version of the C++ implementation of the Davidson-Liu diagonalizer.\n\n");
    printf(" nci=%d, nroots=%d, inidav=%d, mindav=%d, maxdav=%d, kitdav=%d, qtol=%7.1le\n",
	   nci, nroots, inidav, mindav, maxdav, kitdav, qtol);
    double percent = 200.0 * (double)Ham_d.nnz / ((double)nci*(double)(nci+1));
    printf(" nnz=%ld (%4.2lf%)\n\n", Ham_d.nnz, percent);
    // fflush(stdout);
  }

  // Get initial wall clock time.
  gettime(0, 0, &twall1);
  twall2 = twall1;

  // Initial basis vectors.
  B_d.submat(nci, inidav).init(jrefconf_d);

  // Initialize projected Hamiltonian (G).
  P_d.submat(nci, numdav) =  Ham_d                   * B_d.submat(nci, numdav);
  G_d.subsym(numdav)      = ~B_d.submat(nci, numdav) * P_d.submat(nci, numdav);

  if (prtlevel >= 2)
  {
    // Print header for list of Davidson iterations.
    printf(" Iter.   Basis   To go   Worst      Eigenvalue           Norm      Wall clock   Total wall\n");
    //         1       30      5       5     3.123456789012345    0.1234E+12     1.23 s      123.45 s
    // fflush(stdout);
  }

  // Perform Davidson iterations.
  for (int iter=1;;)
  {
    // Temporary copy of G, will be destroyed by solving eigenvalue problem.
    A_h.submat(numdav, numdav).copyFromDevice(G_d);

    // Eigenvalues und eigenvectors of the submatrix.
    // A_h.submat(numdav, numdav).diag(eig_h);
    // A_d.submat(numdav, numdav).copyFromHost(A_h);
    // eig_d.copyFromHost(eig_h, numdav);

    // Unused columns of B and P will be used as scratch space.
    // If less than nroots columns are unused, a restart will be performed.
    if (numdav + nroots <= maxdav)
    {
      // Calculate nroots eigenvalues und eigenvectors of the submatrix.
      A_h.submat(numdav, numdav).eigen(U_h.submat(numdav, nroots), eig_h);
      A_d.submat(numdav, nroots).copyFromHost(U_h);
      eig_d.copyFromHost(eig_h, nroots);

      // Calculate Q matrix and column vector norms (appending Q to B).
      B_d.submat(nci, nroots, 0, numdav)  = P_d.submat(nci, numdav) * A_d.submat(numdav, nroots);
      U_d.submat(numdav, nroots).multByDiag(A_d.submat(numdav, nroots), eig_d);
      B_d.submat(nci, nroots, 0, numdav) -= B_d.submat(nci, numdav) * U_d.submat(numdav, nroots);
      B_d.submat(nci, nroots, 0, numdav).colnorms(qnorm_h);

      // Find roots that have not converged.
      int newdav = 0;  // Number of new basis vectors (one per unconverged root).

      for (int k=0; k<nroots; ++k)
        if (qnorm_h[k] > qtol)
          iroot_h[newdav++] = k;

      // Find root with the largest norm.
      int kbad = 0;  // Index of "worst" root (with the largest norm).

      for (int k=1; k<nroots; ++k)
        if (qnorm_h[k] > qnorm_h[kbad])
          kbad = k;

      // Get wall clock time.
      gettime(0, 0, &twall3);

      // Print information about iterations:
      if (prtlevel >= 2)
      {
        printf("%4d%9d%7d%8d%22.15lf%14.4le%9.2lf s%12.2lf s\n",
               iter, numdav, newdav, kbad+1, eig_h[kbad], qnorm_h[kbad],
               twall3-twall2, twall3-twall1);
        // fflush(stdout);
      }

      // Save current wall clock time:
      twall2 = twall3;

      // Terminate if all roots have converged.
      if (!newdav)
        break;

      // Error return if maximum number of iterations has been reached.
      if (iter >= kitdav)
      {
        printf(" Maximum number of Davidson iterations reached.\n");
        icall = -1;
        return;
      }

      // Calculate new basis vectors from Q matrix (appending Q' to P).
      iroot_d.copyFromHost(iroot_h, newdav);
      AddBasisVectors(P_d.submat(nci, newdav, 0, numdav), B_d.submat(nci, nroots, 0, numdav), Hii_d, eig_d, iroot_d, qtol);

      // Normalization factors of the column vectors of Q'.
      P_d.submat(nci, newdav, 0, numdav).colnorms(qnorm_h);
      for (int k=0; k<newdav; ++k)
        qnorm_h[k] = 1.0 / qnorm_h[k];

      // Normalize column vectors of Q' (in P).
      qnorm_d.copyFromHost(qnorm_h, newdav);
      P_d.submat(nci, newdav, 0, numdav).multByDiag(qnorm_d);

      // Schmidt orthogonalize (appending surviving vectors to B in the second step).
      Schmidt1(B_d.submat(nci, numdav),            P_d.submat(nci, newdav, 0, numdav), U_d.submat(numdav, newdav));
      Schmidt2(B_d.submat(nci, newdav, 0, numdav), P_d.submat(nci, newdav, 0, numdav), v_d, newdav, qtol);

      // Terminate if no new basis vectors have survived.
      if (newdav == 0)
        break;

      // Repeat orthogonalization (to remove numerical noise).
      Schmidt1(B_d.submat(nci, numdav),            B_d.submat(nci, newdav, 0, numdav), U_d.submat(numdav, newdav));
      Schmidt3(B_d.submat(nci, newdav, 0, numdav), v_d);

      // Update product matrix (P).
      P_d.submat(nci, newdav, 0, numdav)    = Ham_d * B_d.submat(nci, newdav, 0, numdav);

      // Update projected Hamiltonian (G) and counters.
      int numold  = numdav;
      numdav     += newdav;
      G_d.submat(newdav, numdav, numold, 0) = ~B_d.submat(nci, newdav, 0, numold) * P_d.submat(nci, numdav);
      ++iter;
    }
    else
    {
      // Perform a restart.

      // Calculate mindav eigenvalues und eigenvectors of the submatrix.
      A_h.submat(numdav, numdav).eigen(U_h.submat(numdav, mindav), eig_h);
      A_d.submat(numdav, mindav).copyFromHost(U_h);
      eig_d.copyFromHost(eig_h, mindav);

      // Collapse basis vectors (in B) into P and swap with B.
      P_d.submat(nci, mindav) = B_d.submat(nci, numdav) * A_d.submat(numdav, mindav);
      P_d.swap(B_d);
      numdav = mindav;

      // Schmidt orthogonalize (twice).
      Schmidt3(B_d.submat(nci, mindav), v_d);
      Schmidt3(B_d.submat(nci, mindav), v_d);

      // Calculate product matrix.
      P_d.submat(nci, numdav) = Ham_d * B_d.submat(nci, numdav);

      // Update projected Hamiltonian (G).
      G_d.subsym(numdav)      = ~B_d.submat(nci, numdav) * P_d.submat(nci, numdav);

      if (prtlevel >= 2)
        printf(" The Davidson algorithm has been restarted.\n");
    }
  }

  // Collapse basis vectors (in B) into P.
  P_d.submat(nci, nroots) = B_d.submat(nci, numdav) * A_d.submat(numdav, nroots);
  P_d.submat(nci, nroots).copyToHost(C, ldc);

  // Copy eigenvalues.
  for (int k=0; k<nroots; ++k)
    E[k] = eig_h[k];

  // Print statistics.
  long iend = Timing::GetTime(CLOCK_MONOTONIC);
  Timing::Print();
  printf("\n Overall wall clock time: %4.2lf s\n", 1e-9 * (double)(iend - ibegin));

  // Flush output to avoid mixing of Fortran and C output.
  fflush(stdout);
}


// Static members.
Counters Timing::ncalls;
Counters Timing::twall;

} // End of namespace.


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Interface functions.                                                                                                   ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern "C"
{
  void gpudavliu_(double *Hii,
                  double *Hij,
                  int    *ifirst,
                  int    *icol,
                  int    &nci,
                  int    *jrefconf,
                  int    &inidav,
                  double *E,
                  double *C,
                  int    &ldc,
                  int    &nroots,
                  int    &mindav,
                  int    &maxdav,
                  int    &kitdav,
                  double &qtol,
                  int    &prtlevel,
                  int    &icall)
  {
    gpu::DavLiu(Hii, Hij, ifirst, icol, nci, jrefconf, inidav, E, C, ldc, nroots, mindav, maxdav, kitdav, qtol, prtlevel, icall);
  }
}
