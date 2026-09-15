// SPDX-FileCopyrightText: Copyright (c) 2009-2015-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]


#ifndef DISKSTACK_NG_H
#define DISKSTACK_NG_H

// All these includes are in PsimagLite
#include "BasisTraits.hh"
#include <PsimagLite/Io/IoNg.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/Stack.h>
#include <exception>

// A disk stack, similar to std::stack but stores in disk not in memory
namespace Dmrg {
template <typename DataType> class DiskStack {

	using IoInType  = typename PsimagLite::IoNg::In;
	using IoOutType = typename PsimagLite::IoNg::Out;

public:

	DiskStack(const PsimagLite::String filename,
	          bool                     needsToRead,
	          PsimagLite::String       label,
	          const BasisTraits&       basisTraits)
	    : ioOut_((needsToRead) ? 0 : new IoOutType(filename, PsimagLite::IoNg::ACC_RDW))
	    , ioIn_((needsToRead) ? new IoInType(filename) : 0)
	    , label_("DiskStack" + label)
	    , basisTraits_(basisTraits)
	    , total_(0)
	    , progress_("DiskStack")
	    , dt_(0)
	{
		if (!needsToRead) {
			ioOut_->createGroup(label_);
			ioOut_->write(total_, label_ + "/Size");
			return;
		}

		ioIn_->read(total_, label_ + "/Size");
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Read from file " + filename + " succeeded";
		progress_.printline(msgg, std::cout);
	}

	~DiskStack()
	{
		delete dt_;
		dt_ = 0;
		delete ioIn_;
		ioIn_ = 0;
		delete ioOut_;
		ioOut_ = 0;
	}

	void flush()
	{
		assert(ioOut_);
		ioOut_->flush();
	}

	bool inDisk() const { return true; }

	void push(const DataType& d)
	{
		assert(ioOut_);

		try {
			d.write(*ioOut_,
			        label_ + "/" + ttos(total_),
			        IoOutType::Serializer::NO_OVERWRITE,
			        DataType::SaveEnum::ALL);
		} catch (...) {
			d.write(*ioOut_,
			        label_ + "/" + ttos(total_),
			        IoOutType::Serializer::ALLOW_OVERWRITE,
			        DataType::SaveEnum::ALL);
		}

		++total_;

		ioOut_->write(total_, label_ + "/Size", IoOutType::Serializer::ALLOW_OVERWRITE);
	}

	void pop()
	{
		if (total_ == 0)
			err("Can't pop; the stack is empty!\n");

		--total_;

		if (!ioOut_)
			return;

		ioOut_->write(total_, label_ + "/Size", IoOutType::Serializer::ALLOW_OVERWRITE);
	}

	void restore(SizeType total)
	{
		total_ = total;
		if (!ioOut_)
			return;

		ioOut_->write(total_, label_ + "/Size", IoOutType::Serializer::ALLOW_OVERWRITE);
	}

	const DataType& top() const
	{
		if (!ioIn_)
			err("DiskStack::top() called with ioIn_ as nullptr\n");

		assert(total_ > 0);
		delete dt_;
		dt_ = 0;
		dt_ = new DataType(*ioIn_, label_ + "/" + ttos(total_ - 1), basisTraits_);
		return *dt_;
	}

	SizeType size() const { return total_; }

private:

	DiskStack(const DiskStack&);

	DiskStack& operator=(const DiskStack&);

	IoOutType*                    ioOut_;
	IoInType*                     ioIn_;
	PsimagLite::String            label_;
	const BasisTraits&            basisTraits_;
	int                           total_;
	PsimagLite::ProgressIndicator progress_;
	mutable DataType*             dt_;
}; // class DiskStack

} // namespace Dmrg

#endif
