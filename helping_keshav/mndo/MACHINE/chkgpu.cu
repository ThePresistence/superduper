#include <cstring>


extern "C"
{
  void  chkgpu (int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu);
  void  chkgpu_(int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu);
  void _chkgpu (int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu);
  void _chkgpu_(int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu);
}


void chkgpu(int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu)
{
  // Query number of CUDA devices.
  if (cudaGetDeviceCount(&numgpu) != cudaSuccess)
  {
    numgpu = 0;
    return;
  }

  // Limit number of CUDA devices.
  if (numgpu > maxgpu)
    numgpu = maxgpu;

  // Initialize device names.
  memset(gpunam, ' ', 80 * numgpu);

  // Query properties of each CUDA device.
  for (int igpu=0; igpu<numgpu; ++igpu)
  {
    cudaDeviceProp prop;

    if (cudaGetDeviceProperties(&prop, igpu) == cudaSuccess)
    {
      int len      = strlen(prop.name);
      memcpy(gpunam + 80 * igpu, prop.name, (len < 80) ? len : 80);
      nvcapa[igpu] = 10 * prop.major + prop.minor;
      mibglo[igpu] = static_cast<int>(static_cast<long long>(prop.totalGlobalMem + 524287LL) / 1048576LL);
    }
    else
    {
      // Querying device properties failed,
      // limit number of CUDA devices.
      numgpu = igpu;
      return;
    }
  }
}


void _chkgpu(int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu)
{
  chkgpu(numgpu, gpunam, nvcapa, mibglo, maxgpu);
}


void chkgpu_(int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu)
{
  chkgpu(numgpu, gpunam, nvcapa, mibglo, maxgpu);
}


void _chkgpu_(int &numgpu, char *gpunam, int *nvcapa, int *mibglo, int maxgpu)
{
  chkgpu(numgpu, gpunam, nvcapa, mibglo, maxgpu);
}
