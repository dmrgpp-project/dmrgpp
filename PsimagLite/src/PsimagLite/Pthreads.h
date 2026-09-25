// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Pthreads .h
 *
 *  A C++ pthreads class that implements the Concurrency interface
 *
 */
#ifndef PTHREADS_HEADER_H
#define PTHREADS_HEADER_H

#include "AllocatorCpu.h"
#include <iostream>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#ifdef _GNU_SOURCE
#include <cstring>
#include <errno.h>
#endif

template <typename PthreadFunctionHolderType> struct PthreadFunctionStruct {
	PthreadFunctionStruct()
	    : pfh(0)
	    , threadNum(0)
	    , blockSize(0)
	    , total(0)
	    , mutex(0)
	    , cpu(0)
	{ }

	PthreadFunctionHolderType* pfh;
	int                        threadNum;
	SizeType                   blockSize;
	SizeType                   total;
	pthread_mutex_t*           mutex;
	SizeType                   cpu;
};

template <typename PthreadFunctionHolderType> void* thread_function_wrapper(void* dummyPtr)
{
	PthreadFunctionStruct<PthreadFunctionHolderType>* pfs
	    = (PthreadFunctionStruct<PthreadFunctionHolderType>*)dummyPtr;

	PthreadFunctionHolderType* pfh = pfs->pfh;

	int s = 0;
#ifdef __linux___
	s = sched_getcpu();
#endif
	if (s >= 0)
		pfs->cpu = s;

	pfh->thread_function_(pfs->threadNum, pfs->blockSize, pfs->total, pfs->mutex);

	int retval = 0;
	pthread_exit(static_cast<void*>(&retval));
	return 0;
}

namespace PsimagLite {
template <typename PthreadFunctionHolderType> class Pthreads {

public:

	Pthreads(SizeType npthreads, int = 0)
	    : nthreads_(npthreads)
	    , cores_(1)
	{
		std::cerr << "Pthreads is deprecated, please use PthreadsNg\n";
		int cores = sysconf(_SC_NPROCESSORS_ONLN);
		cores_    = (cores > 0) ? cores : 1;
	}

	void loopCreate(SizeType total, PthreadFunctionHolderType& pfh)
	{
		PthreadFunctionStruct<PthreadFunctionHolderType>* pfs;
		pfs = new PthreadFunctionStruct<PthreadFunctionHolderType>[nthreads_];
		pthread_mutex_init(&(mutex_), NULL);
		pthread_t*       thread_id = new pthread_t[nthreads_];
		pthread_attr_t** attr      = new pthread_attr_t*[nthreads_];

		for (SizeType j = 0; j < nthreads_; j++) {
			pfs[j].threadNum = j;
			pfs[j].pfh       = &pfh;
			pfs[j].total     = total;
			pfs[j].blockSize = total / nthreads_;
			if (total % nthreads_ != 0)
				pfs[j].blockSize++;
			pfs[j].mutex = &mutex_;

			attr[j] = new pthread_attr_t;
			int ret = pthread_attr_init(attr[j]);
			checkForError(ret);

			setAffinity(attr[j], j, cores_);

			ret = pthread_create(&thread_id[j],
			                     attr[j],
			                     thread_function_wrapper<PthreadFunctionHolderType>,
			                     &pfs[j]);
			checkForError(ret);
		}

		for (SizeType j = 0; j < nthreads_; ++j)
			pthread_join(thread_id[j], 0);
		for (SizeType j = 0; j < nthreads_; ++j) {
			int ret = pthread_attr_destroy(attr[j]);
			checkForError(ret);
			delete attr[j];
			attr[j] = 0;
		}

		delete[] attr;

#ifndef NDEBUG
#ifdef __linux__
		for (SizeType j = 0; j < nthreads_; j++) {
			std::cout << "Pthreads: Pthread number " << j << " runs on core number ";
			std::cout << pfs[j].cpu << "\n";
		}
#endif
#endif

		pthread_mutex_destroy(&mutex_);
		delete[] thread_id;
		delete[] pfs;
	}

	String name() const { return "pthreads"; }

	SizeType threads() const { return nthreads_; }

	SizeType mpiProcs() const { return 1; }

private:

	void setAffinity(std::thread& thread, SizeType threadNum, SizeType cores) const { }

	void checkForError(int ret) const
	{
		if (ret == 0)
			return;
#ifdef _GNU_SOURCE
		std::cerr << "Pthreads ERROR: " << strerror(ret) << "\n";
#endif
	}

	SizeType        nthreads_;
	SizeType        cores_;
	pthread_mutex_t mutex_;
}; // Pthreads class

} // namespace PsimagLite

/*@}*/
#endif
