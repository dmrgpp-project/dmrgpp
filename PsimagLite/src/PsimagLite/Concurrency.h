// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Concurrency.h
 *
 */
#ifndef CONCURRENCY_HEADER_H
#define CONCURRENCY_HEADER_H
#include "CodeSectionParams.h"
#include "FloatingPoint.h"
#include "LAPACK.h"
#include "LabelDisabled.h"
#include "NotMpi.h"
#include "Vector.h"
#include <cassert>
#include <stdexcept>
#include <sys/syscall.h>
#include <unistd.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif

namespace PsimagLite {

class Concurrency {

	using MpiDisabledType = LabelDisabled;

public:

	static SizeType          mode;
	static CodeSectionParams codeSectionParams;

#ifndef USE_PTHREADS

	using MutexType    = int;
	using PthreadtType = SizeType;

	static void mutexLock(MutexType*) { }

	static void mutexUnlock(MutexType*) { }

	static void mutexInit(MutexType*) { }

	static void mutexDestroy(MutexType*) { }

	static PthreadtType threadSelf() { return 0; }

#else

	using MutexType    = pthread_mutex_t;
	using PthreadtType = pthread_t;

	static void mutexInit(MutexType* mutex)
	{
		if (pthread_mutex_init(mutex, 0) != 0)
			std::cerr << "WARNING: mutexInit returned non zero\n";
	}

	static void mutexDestroy(MutexType* mutex)
	{
		if (pthread_mutex_destroy(mutex) != 0)
			std::cerr << "WARNING: mutexDestroy returned non zero\n";
	}

	static void mutexLock(MutexType* mutex)
	{
		if (pthread_mutex_lock(mutex) != 0)
			std::cerr << "WARNING: mutexLock returned non zero\n";
	}

	static void mutexUnlock(MutexType* mutex)
	{
		if (pthread_mutex_unlock(mutex) != 0)
			std::cerr << "WARNING: mutexUnlock returned non zero\n";
	}

	// DON'T EVEN think of using this as a thread id
	// This ISN'T 0, 1, 2, ...
	// it's a rather large number
	// Only assured thing is its uniqueness per PROCESS
	static PthreadtType threadSelf() { return pthread_self(); }

#endif

	enum
	{
		SERIAL           = 0,
		PTHREADS         = 1,
		MPI              = 2,
		PTHREADS_AND_MPI = 3
	};

	static SizeType storageSize(SizeType npthreads)
	{
		switch (mode) {
		case SERIAL:
			assert(npthreads == 1);
		case PTHREADS:
		case PTHREADS_AND_MPI:
			return npthreads;
		case MPI:
			return 1;
		}
		throw RuntimeError("storageSize: wrong mode\n");
	}

	static SizeType storageIndex(SizeType threadNum)
	{
		switch (mode) {
		case SERIAL:
			assert(threadNum == 0);
		case PTHREADS:
		case PTHREADS_AND_MPI:
			return threadNum;
		case MPI:
			return 0;
		}
		throw RuntimeError("storageIndex: wrong mode\n");
	}

	Concurrency(int* argc, char*** argv, size_t nthreads)
	{
		FloatingPoint::enableExcept();
		codeSectionParams.npthreads         = nthreads;
		codeSectionParams.npthreadsLevelTwo = 1;
		mode                                = 0;
#ifdef USE_PTHREADS
		mode |= 1;
		if (!psimag::LAPACK::isThreadSafe())
			std::cerr << "WARNING: You LAPACK might not be thread safe\n";
#else
		if (nthreads != 1)
			throw RuntimeError("nthreads>1 but no USE_PTHREADS "
			                   "support compiled\n");
#endif

		MPI::init(argc, argv);
		MPI::info(std::cout);
		MPI::version(std::cout);

		if (MPI::hasMpi())
			mode |= 2;
	}

	~Concurrency() { MPI::finalize(); }

	static bool root(MPI::CommType comm = MPI::COMM_WORLD)
	{
		return (MPI::commRank(comm) == 0);
	}

	static SizeType nprocs(MPI::CommType comm = MPI::COMM_WORLD) { return MPI::commSize(comm); }

	static SizeType rank(MPI::CommType comm = MPI::COMM_WORLD) { return MPI::commRank(comm); }

	static bool hasMpi() { return (mode & MPI); }

	static bool hasPthreads() { return (mode & PTHREADS); }

	static void mpiDisable(String label)
	{
		if (!hasMpi())
			return;
		mpiDisabled_.disable(label);
	}

	static void
	mpiDisableIfNeeded(SizeType& mpiRank, SizeType& blockSize, String label, SizeType total)
	{
		if (!hasMpi())
			return;
		if (!mpiDisabled_(label))
			return;
		mpiRank   = 0;
		blockSize = total;
		if (!hasPthreads())
			return;
		String str(__FILE__);
		str += " mpiDisableIfNeeded label = " + label + "\n";
		throw RuntimeError(str);
	}

	static bool isMpiDisabled(String label)
	{
		if (!hasMpi())
			return false;
		return mpiDisabled_(label);
	}

	static void setOptions(const CodeSectionParams& cs)
	{
		codeSectionParams = cs;
		if (codeSectionParams.npthreads == 1 && codeSectionParams.npthreadsLevelTwo == 1)
			return;

#ifndef USE_PTHREADS
		PsimagLite::String message1(__FILE__);
		message1 += " FATAL: You are requesting nthreads > 0 but you ";
		message1 += "did not compile with USE_PTHREADS enabled\n";
		message1 += " Either set Threads=1 in the input file (you won't ";
		message1 += "have threads though) or\n";
		message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile ";
		message1 += "and recompile\n";
		throw PsimagLite::RuntimeError(message1.c_str());
#else
		std::cout << "Concurrency::npthreads=" << codeSectionParams.npthreads << "\n";
		std::cout << "Concurrency::npthreads2=" << codeSectionParams.npthreadsLevelTwo
		          << "\n";
		std::cout << "Concurrency::setAffinitiesDefault=" << codeSectionParams.setAffinities
		          << "\n";
#endif
	}

private:

	static MpiDisabledType mpiDisabled_;
};

} // namespace PsimagLite

/*@}*/
#endif
