#ifndef PARALLELIZER2_H
#define PARALLELIZER2_H
#include "Concurrency.h"
#include "Vector.h"

#ifdef USE_PTHREADS
#include "ParallelizerPthread.h"
#else

#ifdef _OPENMP
#include "ParallelizerOpenMP.h"
#else
#include "ParallelizerSerial.h"
#endif

#endif

#endif // PARALLELIZER2_H
