// Memory class

#ifndef FFT_MEMORY_H
#define FFT_MEMORY_H

#include <stdint.h>

namespace FFTMPI_NS {

class Memory {
 public:
  Memory() {}
  void *smalloc(int64_t n);
  void *srealloc(void *, int64_t n);
  void sfree(void *);
};

}

#endif
