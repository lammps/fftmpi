// Error class

#ifndef FFT_ERROR_H
#define FFT_ERROR_H

#include <mpi.h>

namespace FFTMPI_NS {

class Error {
 public:
  Error(MPI_Comm);
  void all(const char *);
  void one(const char *);
  void warning(const char *);

 private:
  MPI_Comm world;
};

}

#endif
