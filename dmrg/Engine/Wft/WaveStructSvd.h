// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file WaveStructSvd.h
 *
 *  DOC NEEDED FIXME (This file should go in Wft/ directory perhaps)
 */
#ifndef WAVE_STRUCT_SVD_H
#define WAVE_STRUCT_SVD_H
#include "BasisTraits.hh"
#include "ProgramGlobals.h"
#include <PsimagLite/Io/IoNg.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename LeftRightSuperType_> struct WaveStructSvd {

	using LeftRightSuperType      = LeftRightSuperType_;
	using BasisWithOperatorsType  = typename LeftRightSuperType::BasisWithOperatorsType;
	using BlockDiagonalMatrixType = typename BasisWithOperatorsType::BlockDiagonalMatrixType;
	using OperatorType            = typename BasisWithOperatorsType::OperatorType;
	using SparseMatrixType        = typename OperatorType::StorageType;
	using SparseElementType       = typename SparseMatrixType::value_type;
	using BasisType               = typename BasisWithOperatorsType::BasisType;
	using QnType                  = typename BasisType::QnType;
	using VectorSizeType          = PsimagLite::Vector<SizeType>::Type;
	using RealType                = typename BasisWithOperatorsType::RealType;
	using VectorRealType          = typename PsimagLite::Vector<RealType>::Type;
	using VectorVectorRealType    = typename PsimagLite::Vector<VectorRealType>::Type;
	using MatrixType              = typename PsimagLite::Matrix<SparseElementType>;
	using VectorMatrixType        = typename PsimagLite::Vector<MatrixType>::Type;
	using VectorQnType            = typename BasisWithOperatorsType::VectorQnType;
	using IoInType                = typename PsimagLite::IoNg::In;

	enum class SaveEnum
	{
		ALL
	};

	WaveStructSvd() { }

	WaveStructSvd(IoInType& io, PsimagLite::String label, const BasisTraits&)
	{
		read(label, io.serializer());
	}

	WaveStructSvd(const BlockDiagonalMatrixType& u,
	              const VectorMatrixType&        vts,
	              const VectorVectorRealType&    s,
	              const VectorQnType&            qns)
	    : u_(u)
	    , vts_(vts)
	    , // Not yet used, will be used by WftAccelSvd
	    s_(s)
	    , // Not yet used, will be used by WftAccelSvd
	    qns_(qns)
	{ }

	const BlockDiagonalMatrixType& u() const { return u_; }

	const VectorMatrixType& vts() const { return vts_; }

	const VectorVectorRealType& s() const { return s_; }

	const VectorQnType& qns() const { return qns_; }

	void read(PsimagLite::IoNg::In& io, PsimagLite::String prefix) { io.read(*this, prefix); }

	void read(PsimagLite::String prefix, PsimagLite::IoNgSerializer& io)
	{
		u_.read(prefix + "/u", io);
		io.read(vts_, prefix + "/vts");
		io.read(s_, prefix + "/s");
		QnType::readVector(qns_, prefix + "/qns", io);
	}

	void write(PsimagLite::IoNg::Out&                io,
	           PsimagLite::String                    prefix,
	           PsimagLite::IoNgSerializer::WriteMode writeMode,
	           SaveEnum) const
	{
		if (writeMode != PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
			io.createGroup(prefix);
		io.write(u_, prefix + "/u", writeMode);
		io.write(vts_, prefix + "/vts", writeMode);
		io.write(s_, prefix + "/s", writeMode);
		io.write(qns_, prefix + "/qns", writeMode);
	}

	void write(PsimagLite::String prefix, PsimagLite::IoNgSerializer& io) const
	{
		io.createGroup(prefix);
		u_.write(prefix + "/u", io);
		io.write(prefix + "/vts", vts_);
		io.write(prefix + "/s", s_);
		io.write(prefix + "/qns", qns_);
	}

private:

	BlockDiagonalMatrixType u_;
	VectorMatrixType        vts_;
	VectorVectorRealType    s_;
	VectorQnType            qns_;
}; // struct WaveStructSvd

} // namespace Dmrg

/*@}*/
#endif // WAVE_STRUCT_SVD_H
