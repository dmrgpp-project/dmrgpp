// SPDX-FileCopyrightText: Copyright (c) 2008-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ObserverHelper.h
 *
 *  A class to read and serve precomputed data to the observer
 *
 */
#ifndef PRECOMPUTED_H
#define PRECOMPUTED_H
#include "DmrgSerializer.h"
#include "ProgramGlobals.h"
#include "TimeSerializer.h"
#include "VectorWithOffsets.h" // to include norm
#include <PsimagLite/GetBraOrKet.h>
#include <PsimagLite/Io/HDF5DisableExceptionPrinting.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/SparseVector.h>

#include <memory>

namespace Dmrg {

template <typename IoInputType_,
          typename MatrixType_,
          typename VectorType_,
          typename VectorWithOffsetType_,
          typename LeftRightSuperType>
class ObserverHelper {

public:

	using IoInputType                = IoInputType_;
	using MatrixType                 = MatrixType_;
	using VectorType                 = VectorType_;
	using VectorWithOffsetType       = VectorWithOffsetType_;
	using IndexType                  = SizeType;
	using FieldType                  = typename VectorType::value_type;
	using BasisWithOperatorsType     = typename LeftRightSuperType::BasisWithOperatorsType;
	using RealType                   = typename BasisWithOperatorsType::RealType;
	using SparseMatrixType           = typename BasisWithOperatorsType::SparseMatrixType;
	using TimeSerializerType         = TimeSerializer<VectorWithOffsetType>;
	using BasisType                  = typename BasisWithOperatorsType::BasisType;
	using OperatorType               = typename BasisWithOperatorsType::OperatorType;
	using DmrgSerializerType         = DmrgSerializer<LeftRightSuperType, VectorWithOffsetType>;
	using BlockDiagonalMatrixType    = typename DmrgSerializerType::BlockDiagonalMatrixType;
	using FermionSignType            = typename DmrgSerializerType::FermionSignType;
	using VectorSizeType             = PsimagLite::Vector<SizeType>::Type;
	using VectorShortIntType         = PsimagLite::Vector<short int>::Type;
	using GetBraOrKetType            = PsimagLite::GetBraOrKet;
	using PairLeftRightSuperSizeType = std::pair<std::unique_ptr<LeftRightSuperType>, SizeType>;

	enum class SaveEnum
	{
		YES,
		NO
	};

	ObserverHelper(IoInputType& io,
	               SizeType     start,
	               SizeType     nf,
	               SizeType     trail,
	               bool         withLegacyBugs,
	               bool         readOnDemand)
	    : io_(io)
	    , withLegacyBugs_(withLegacyBugs)
	    , readOnDemand_(readOnDemand)
	    , progress_("ObserverHelper")
	    , noMoreData_(false)
	    , numberOfSites_(0)
	    , lrsStorage_(PairLeftRightSuperSizeType(nullptr, 0))
	{
		bool hasConcurrency
		    = (PsimagLite::Concurrency::codeSectionParams.npthreads > 1
		       || PsimagLite::Concurrency::codeSectionParams.npthreadsLevelTwo > 1);
		if (hasConcurrency && readOnDemand_)
			err(std::string("ReadOnDemand does not support threading. ")
			    + "Set Threads=1 in input, or use -S 1 in command line.\n");

		typename BasisWithOperatorsType::VectorBoolType odds;
		io_.read(odds, "OddElectronsOneSite");
		SizeType n = odds.size();
		signsOneSite_.resize(n);
		for (SizeType i = 0; i < n; ++i)
			signsOneSite_[i] = (odds[i]) ? -1 : 1;

		if (readOnDemand_) {
			std::cout << "ObserverHelper: observeReadOnDemand is ON\n";
			std::cerr << "ObserverHelper: observeReadOnDemand is ON\n";
		}

		if (nf > 0)
			if (!init(start, start + nf, SaveEnum::YES))
				return;

		if (trail > 0)
			if (!init(start, start + trail, SaveEnum::NO))
				return;
	}

	~ObserverHelper() = default;

	const SizeType& numberOfSites() const { return numberOfSites_; }

	bool endOfData() const { return noMoreData_; }

	void transform(SparseMatrixType& ret, const SparseMatrixType& O2, SizeType ind) const
	{
		checkIndex(ind);

		if (!readOnDemand_)
			return dSerializerV_[ind]->transform(ret, O2);

		const PsimagLite::String prefix = "Serializer/" + ttos(ind);
		BlockDiagonalMatrixType  transformStorage(io_, prefix + "/transform", false);
		DmrgSerializerType::transform(ret, O2, transformStorage);
	}

	SizeType cols(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->cols();
	}

	SizeType rows(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->rows();
	}

	short int signsOneSite(SizeType site) const
	{
		assert(site < signsOneSite_.size());
		return signsOneSite_[site];
	}

	const FermionSignType& fermionicSignLeft(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->fermionicSignLeft();
	}

	const FermionSignType& fermionicSignRight(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->fermionicSignRight();
	}

	const LeftRightSuperType& leftRightSuper(SizeType ind) const
	{
		checkIndex(ind);

		if (readOnDemand_) {
			if (ind != lrsStorage_.second) {
				lrsStorage_.first.reset();
			}

			if (!lrsStorage_.first) {
				const PsimagLite::String prefix = "Serializer/" + ttos(ind);

				lrsStorage_.first = std::make_unique<LeftRightSuperType>(
				    io_, prefix, BasisTraits { true, true });
				lrsStorage_.second = ind;
			}

			return *lrsStorage_.first;
		}

		return dSerializerV_[ind]->leftRightSuper();
	}

	ProgramGlobals::DirectionEnum direction(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->direction();
	}

	const VectorWithOffsetType&
	psiConst(SizeType ind, SizeType sectorIndex, SizeType levelIndex) const
	{
		checkIndex(ind);

		return dSerializerV_[ind]->psiConst(sectorIndex, levelIndex);
	}

	RealType time(SizeType ind) const
	{
		if (timeSerializerV_.empty())
			return 0.0;
		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->time();
	}

	SizeType site(SizeType ind) const
	{
		if (timeSerializerV_.empty()) {
			checkIndex(ind);
			return this->siteInternal(this->leftRightSuper(ind), this->direction(ind));
		}

		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->site();
	}

	SizeType size() const { return dSerializerV_.size(); }

	const VectorWithOffsetType& getVectorFromBracketId(const PsimagLite::GetBraOrKet& braOrKet,
	                                                   SizeType index) const
	{
		if (braOrKet.isPvector()) {
			const SizeType pIndex = braOrKet.pIndex();
			return timeVector(pIndex, index);
		}

		return psiConst(index, braOrKet.sectorIndex(), braOrKet.levelIndex());
	}

	const VectorWithOffsetType& timeVector(SizeType braketId, SizeType ind) const
	{
		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->vector(braketId);
	}

	bool withLegacyBugs() const { return withLegacyBugs_; }

private:

	SizeType siteInternal(const LeftRightSuperType&     lrs,
	                      ProgramGlobals::DirectionEnum direction) const
	{
		return (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		    ? lrs.right().block()[0] - 1
		    : lrs.right().block()[0];
	}

	bool init(SizeType start, SizeType end, SaveEnum saveOrNot)
	{
		PsimagLite::String prefix = "Serializer";
		SizeType           total  = 0;
		io_.read(total, prefix + "/Size");
		if (start >= end || start >= total)
			return false;

		if (end > total) {
			end = total;
		}

		for (SizeType i = start; i < end; ++i) {

			auto dSerializer
			    = std::make_unique<DmrgSerializerType>(io_,
			                                           prefix + "/" + ttos(i),
			                                           false,
			                                           BasisTraits { true, true },
			                                           readOnDemand_);

			SizeType tmp = dSerializer->leftRightSuper().sites();
			if (tmp > 0 && numberOfSites_ == 0)
				numberOfSites_ = tmp;

			if (readOnDemand_)
				dSerializer->freeLrs();

			if (saveOrNot == SaveEnum::YES)
				dSerializerV_.push_back(std::move(dSerializer));

			try {
				HDF5DisableExceptionPrinting disable;
				PsimagLite::String           prefix("/TargetingCommon/" + ttos(i));
				auto ts = std::make_unique<TimeSerializerType>(io_, prefix);
				std::cerr << "Read TimeSerializer\n";
				std::cout << "Time = " << ts->time() << '\n';
				if (saveOrNot == SaveEnum::YES)
					timeSerializerV_.push_back(std::move(ts));
			} catch (...) { }

			std::cerr << __FILE__ << " read " << i << " out of " << (end - start)
			          << "\n";
			progress_.printMemoryUsage();
		}

		noMoreData_ = (end == total);
		return !dSerializerV_.empty();
	}

	static SizeType braketStringToNumber(const PsimagLite::String& str)
	{
		GetBraOrKetType ketOrBra(str);
		if (!ketOrBra.isPvector())
			return 0;

		return ketOrBra.levelIndex();
	}

	void checkIndex(SizeType ind) const
	{
		if (ind >= dSerializerV_.size())
			err("Index " + ttos(ind) + " greater or equal to "
			    + ttos(dSerializerV_.size()));

		if (dSerializerV_[ind])
			return;

		err("dSerializerV_ at index " + ttos(ind) + " point to 0x0\n");
	}

	IoInputType&                                                           io_;
	typename PsimagLite::Vector<std::unique_ptr<DmrgSerializerType>>::Type dSerializerV_;
	typename PsimagLite::Vector<std::unique_ptr<TimeSerializerType>>::Type timeSerializerV_;
	const bool                                                             withLegacyBugs_;
	const bool                                                             readOnDemand_;
	PsimagLite::ProgressIndicator                                          progress_;
	bool                                                                   noMoreData_;
	VectorShortIntType                                                     signsOneSite_;
	SizeType                                                               numberOfSites_;
	mutable PairLeftRightSuperSizeType                                     lrsStorage_;
}; // ObserverHelper
} // namespace Dmrg

/*@}*/
#endif
