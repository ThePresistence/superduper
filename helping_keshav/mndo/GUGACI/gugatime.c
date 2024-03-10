#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/* Get CPU and wall clock time. */

void gettime(double *tuser, double *tsys, double *twall)
{
  struct rusage   ru;
  struct timeval  tv;
  struct timezone tz;

  if (tuser  &&  tsys)
  {
    if (!getrusage(RUSAGE_SELF, &ru))
    {
      *tuser = ru.ru_utime.tv_sec + 1e-6 * ru.ru_utime.tv_usec;
      *tsys  = ru.ru_stime.tv_sec + 1e-6 * ru.ru_stime.tv_usec;
    }
    else *tuser = *tsys = 0.;
  }

  if (twall)
  {
    tz.tz_minuteswest = 0;
    tz.tz_dsttime     = 0;

    if (!gettimeofday(&tv, &tz))
    {
      *twall = tv.tv_sec + 1e-6 * tv.tv_usec;
    }
    else *twall = 0.;
  }
}


void gettime_(double *tuser, double *tsys, double *twall)
{
  gettime(tuser, tsys, twall);
}


void _gettime(double *tuser, double *tsys, double *twall)
{
  gettime(tuser, tsys, twall);
}


void _gettime_(double *tuser, double *tsys, double *twall)
{
  gettime(tuser, tsys, twall);
}


void gettime__(double *tuser, double *tsys, double *twall)
{
  gettime(tuser, tsys, twall);
}


void __gettime(double *tuser, double *tsys, double *twall)
{
  gettime(tuser, tsys, twall);
}


/* Calculate the difference between two pointers
   to obtain the leading dimension of a matrix.  */

int ldget(double *a, double *b)
{
  return b - a;
}


int ldget_(double *a, double *b)
{
  return b - a;
}


int _ldget(double *a, double *b)
{
  return b - a;
}


int _ldget_(double *a, double *b)
{
  return b - a;
}


int ldget__(double *a, double *b)
{
  return b - a;
}


int __ldget(double *a, double *b)
{
  return b - a;
}
