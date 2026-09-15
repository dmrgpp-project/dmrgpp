// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetQuantumElectrons.h
 *
 *
 */
#ifndef TargetQuantumElectrons_H
#define TargetQuantumElectrons_H
#include "AlgebraicStringToNumber.h"
#include "ProgramGlobals.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
//! Hubbard Model Parameters
template <typename RealType, typename QnType> class TargetQuantumElectrons {

public:

	using VectorSizeType = PsimagLite::Vector<SizeType>::Type;
	using PairSizeType   = typename QnType::PairSizeType;
	using VectorQnType   = typename QnType::VectorQnType;

	template <typename IoInputType>
	TargetQuantumElectrons(IoInputType& io)
	    : totalNumberOfSites_(0)
	    , isSu2_(false)
	{
		io.readline(totalNumberOfSites_, "TotalNumberOfSites=");
		int tmp = 0;
		try {
			io.readline(tmp, "UseSu2Symmetry=");
		} catch (std::exception&) { }

		isSu2_ = (tmp > 0);

		bool     hasNqns = false;
		SizeType nqns    = 0;
		try {
			io.readline(nqns, "NumberOfTargetQns=");
			hasNqns = true;
		} catch (std::exception&) { }

		if (!hasNqns) {
			readOneTarget(io, "");
			return;
		}

		for (SizeType i = 0; i < nqns; ++i)
			readOneTarget(io, ttos(i));
	}

	SizeType sizeOfOther() const
	{
		const SizeType n = vqn_.size();
		if (n == 0)
			return 0;
		const SizeType answer = vqn_[0].other.size();
		for (SizeType i = 1; i < n; ++i) {
			if (vqn_[i].other.size() == answer)
				continue;
			err("sizeOfOther must be the same for all target qns\n");
		}

		return answer;
	}

	SizeType size() const { return vqn_.size(); }

	const QnType& qn(SizeType ind) const
	{
		assert(ind < vqn_.size());
		return vqn_[ind];
	}

	void updateQuantumSector(VectorQnType&                 quantumSector,
	                         SizeType                      sites,
	                         ProgramGlobals::DirectionEnum direction,
	                         SizeType                      step,
	                         const VectorQnType&           adjustQuantumNumbers) const
	{
		const SizeType maxSites = totalNumberOfSites_;

		if (direction == ProgramGlobals::DirectionEnum::INFINITE && sites < maxSites
		    && adjustQuantumNumbers.size() > step) {
			if (quantumSector.size() != 1)
				err("adjustQuantumNumbers only with single target\n");
			quantumSector[0] = adjustQuantumNumbers[step];
			return;
		} else {
			quantumSector = vqn_;
		}

		const SizeType n = quantumSector.size();

		for (SizeType i = 0; i < n; ++i)
			quantumSector[i].scale(sites, totalNumberOfSites_, direction, isSu2_);
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/TargetQuantumElectrons";
		io.createGroup(label);
		io.write(label + "/TotalNumberOfSites", totalNumberOfSites_);
		io.write(label + "/isSu2", isSu2_);
		vqn_.write(label + "/qn", io);
	}

private:

	TargetQuantumElectrons(const TargetQuantumElectrons&);

	TargetQuantumElectrons& operator=(const TargetQuantumElectrons&);

	template <typename IoInputType>
	void readOneTarget(IoInputType& io, const PsimagLite::String label)
	{
		QnType         qn(QnType::zero());
		VectorSizeType qnOther;
		const bool     allowUpDown = true;

		PsimagLite::String msg("TargetQuantumElectrons: ");
		bool               hasTwiceJ = false;
		try {
			io.readline(qn.jmPair.first, "TargetSpinTimesTwo" + label + "=");
			hasTwiceJ = true;
		} catch (std::exception&) { }

		SizeType ready = 0;
		if (allowUpDown) {
			SizeType electronsUp   = 0;
			SizeType electronsDown = 0;
			try {
				electronsUp
				    = readNumberOrExpression(io, "TargetElectronsUp" + label + "=");
				electronsDown = readNumberOrExpression(
				    io, "TargetElectronsDown" + label + "=");
				SizeType tmp    = electronsUp + electronsDown;
				qn.oddElectrons = (tmp & 1);
				qnOther.push_back(tmp);
				qnOther.push_back(electronsUp);
				ready = 2;
			} catch (std::exception&) { }
		}

		try {
			SizeType tmp
			    = readNumberOrExpression(io, "TargetElectronsTotal" + label + "=");
			qn.oddElectrons = (tmp & 1);
			qnOther.push_back(tmp);
			ready++;
		} catch (std::exception&) { }

		try {
			SizeType szPlusConst
			    = readNumberOrExpression(io, "TargetSzPlusConst" + label + "=");
			qnOther.push_back(szPlusConst);
		} catch (std::exception&) { }

		if (ready == 3) {
			msg += "Provide either up/down or total/sz but not both.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		try {
			std::string parity;
			io.readline(parity, "TargetFermionicParity=");
			if (ready != 0) {
				err("TargetFermionicParity cannot be used if you provide "
				    "TargetElectrons*\n");
			}

			if (parity != "even" && parity != "odd") {
				err("TargetFermionicParity must be either even or odd, not "
				    + parity + "\n");
			}

			qn.oddElectrons = (parity == "odd");
			std::cout << "Using TargetFermionicParity=" << qn.oddElectrons << "\n";
			std::cerr << "Using TargetFermionicParity=" << qn.oddElectrons << "\n";
		} catch (std::exception&) { }

		bool flag = false;
		try {
			readNumberOrExpression(io, "TargetExtra" + label + "=");
			flag = true;
		} catch (std::exception&) { }

		if (flag)
			err("Instead of TargetExtra" + label + "= please use a vector\n");

		try {
			VectorSizeType extra;
			io.read(extra, "TargetExtra" + label);
			for (SizeType i = 0; i < extra.size(); ++i)
				qnOther.push_back(extra[i]);
		} catch (std::exception&) { }

		qn.other.fromStdVector(qnOther);

		if (isSu2_ && !hasTwiceJ) {
			msg += "Please provide TargetSpinTimesTwo when running with SU(2).\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (isSu2_)
			qn.oddElectrons = (totalNumberOfSites_ & 1);

		vqn_.push_back(qn);
	}

	template <typename IoInputType>
	SizeType readNumberOrExpression(IoInputType& io, PsimagLite::String fullLabel)
	{
		using AlgebraicStringToNumberType = AlgebraicStringToNumber<RealType>;

		PsimagLite::String val;
		io.readline(val, fullLabel);

		const PsimagLite::String    msg = "Target number for " + fullLabel;
		AlgebraicStringToNumberType algebraicStringToNumber(msg, totalNumberOfSites_);

		SizeType p = algebraicStringToNumber.procLength(val);

		std::cout << fullLabel << p << "\n";
		return p;
	}

	SizeType     totalNumberOfSites_;
	bool         isSu2_;
	VectorQnType vqn_;
};
} // namespace Dmrg

/*@}*/
#endif
