#ifndef CODESECTION_PARAMS_H
#define CODESECTION_PARAMS_H
#include "AllocatorCpu.h"

namespace PsimagLite {

struct CodeSectionParams {

	CodeSectionParams(SizeType threads)
	    : npthreads(threads)
	    , npthreadsLevelTwo(1)
	    , setAffinities(false)
	    , stackSize(0)
	    , name("")
	{ }

	CodeSectionParams(SizeType threads,
	                  SizeType nthreadsLevel2,
	                  bool     a,
	                  size_t   s,
	                  String   sectionName = "")
	    : npthreads(threads)
	    , npthreadsLevelTwo(nthreadsLevel2)
	    , setAffinities(a)
	    , stackSize(s)
	    , name(sectionName)
	{ }

	SizeType npthreads;
	SizeType npthreadsLevelTwo;
	bool     setAffinities;
	size_t   stackSize;
	String   name;
};
} // namespace PsimagLite
#endif // CODESECTION_PARAMS_H
