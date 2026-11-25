// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file SpinSquaredHelper.h
 *
 *  Helper class for SpinSquared
 *
 */
#ifndef SPIN_SQ_HELPER_H
#define SPIN_SQ_HELPER_H

namespace Dmrg {
template <typename FieldType_, typename Word_> class SpinSquaredHelper {

public:

	using FieldType = FieldType_;
	using Word      = Word_;

	SpinSquaredHelper()
	    : data_(0)
	    , ketSaved_(0)
	{ }

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres, SizeType, PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "SpinSquaredHelper";

		const char* start = reinterpret_cast<const char*>(this);
		const char* end   = reinterpret_cast<const char*>(&ketSaved_);
		SizeType    total = mres.memResolv(&data_, end - start, str + " data");

		total += mres.memResolv(&ketSaved_, sizeof(*this) - total, str + " ketSaved");

		return total;
	}

	void operator()(const Word& ket, const Word& bra, const FieldType& value)
	{
		if (ket != bra) {
			throw PsimagLite::RuntimeError("SpinSquaredHelper::operator(): ket!=bra\n");
		}
		if (ket != ketSaved_)
			data_ = 0;
		ketSaved_ = ket;
		data_ += value;
	}

	//! receives m, returns (2*j,m+j)
	std::pair<SizeType, SizeType> getJmPair(const FieldType& m) const
	{
		SizeType j      = getJvalue();
		SizeType mtilde = getMvalue(m, j);
		return std::pair<SizeType, SizeType>(j, mtilde);
	}

	void clear() { data_ = 0; }

	void write(PsimagLite::String, PsimagLite::IoNg::Out::Serializer&) const { }

private:

	int perfectSquareOrCrash(const FieldType& t) const
	{
		FieldType r  = sqrt(t);
		int       ri = int(r);
		if (ri != r)
			PsimagLite::RuntimeError("SpinSquaredHelper:: sqrt(1+4d) not an integer\n");

		return ri;
	}

	SizeType getJvalue() const
	{
		if (data_ < 0)
			PsimagLite::RuntimeError("SpinSquaredHelper::getJvalue(): d<0\n");
		int tmp = perfectSquareOrCrash(1.0 + 4.0 * data_);
		if (tmp < 1)
			PsimagLite::RuntimeError("SpinSquaredHelper::getJvalue(): sqrt(1+4d)<1\n");
		SizeType ret = tmp - 1;
		return ret;
	}

	SizeType getMvalue(const FieldType& m, SizeType j) const
	{
		FieldType tmp = m + j * 0.5;
		if (tmp < 0)
			PsimagLite::RuntimeError("SpinSquaredHelper::getMvalue(): j+m <0\n");
		SizeType ret = SizeType(tmp);
		if (ret != tmp)
			PsimagLite::RuntimeError(
			    "SpinSquaredHelper::getMvalue: j+m not SizeType\n");
		return ret;
	}

	// serializr start class SpinSquaredHelper
	// serializr normal data_
	FieldType data_;
	// serializr normal ketSaved_
	Word ketSaved_;

}; // class SpinSquaredHelper
} // namespace Dmrg

/*@}*/
#endif
