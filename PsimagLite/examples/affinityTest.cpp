// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.]

#include <PsimagLite/Concurrency.h>
#include <PsimagLite/Parallelizer.h>
#include <cstdlib>
#include <iostream>

class MyHelper {

	typedef PsimagLite::Concurrency            ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	MyHelper(SizeType ntasks, SizeType nthreads)
	    : ntasks_(ntasks)
	    , x_(nthreads, 0)
	{ }

	SizeType tasks() const { return ntasks_; }

	int result() const { return x_[0]; }

	void doTask(SizeType taskNumber, SizeType threadNum) { x_[threadNum] += taskNumber; }

	void sync()
	{
		for (SizeType i = 1; i < x_.size(); ++i)
			x_[0] += x_[i];
	}

private:

	SizeType       ntasks_;
	VectorSizeType x_;
}; // class MyHelper

int main(int argc, char* argv[])
{
	typedef PsimagLite::Concurrency ConcurrencyType;

	if (argc != 3) {
		std::cout << "USAGE: " << argv[0] << " nthreads ntasks\n";
		return 1;
	}

	SizeType nthreads = atoi(argv[1]);
	SizeType ntasks   = atoi(argv[2]);

	ConcurrencyType concurrency(&argc, &argv, nthreads);

	typedef MyHelper                             HelperType;
	typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;

	PsimagLite::CodeSectionParams csp(nthreads, 1, true, 0);
	ParallelizerType              threadObject(csp);

	HelperType helper(ntasks, nthreads);

	std::cout << "Using " << threadObject.name();
	std::cout << " with " << nthreads << " threads.\n";
	threadObject.loopCreate(helper);
	helper.sync();
	std::cout << "Sum of all tasks= " << helper.result() << "\n";
}
