#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num () { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif

int main(int argc, char** argv)
{
  return 0;
}
