// SPDX-FileCopyrightText: Copyright (c) 2008-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file OnePointCorrelations.h
 *
 *
 *  A class to perform post-processing calculation of one-point correlations
 *  <state1 | A_i | state2>
 * Note that site below is not necessarily the site, but there's a mapping
 * mapping(site) = i
 */

#ifndef ONE_POINT_H
#define ONE_POINT_H
#include "OffsetVector.hpp"
#include "VectorWithOffsets.h" // for operator*
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/GetBraOrKet.h>

namespace Dmrg {

template <typename ObserverHelperType, typename ModelType> class OnePointCorrelations {

	using MatrixType             = typename ObserverHelperType::MatrixType;
	using VectorType             = typename ObserverHelperType::VectorType;
	using VectorWithOffsetType   = typename ObserverHelperType::VectorWithOffsetType;
	using BasisWithOperatorsType = typename ObserverHelperType::BasisWithOperatorsType;
	using FieldType              = typename VectorType::value_type;
	using RealType               = typename BasisWithOperatorsType::RealType;
	using OffsetVectorType       = OffsetVector<FieldType>;
	using OffsetVectorBaseType   = OffsetVectorBase<FieldType>;

public:

	OnePointCorrelations(const ObserverHelperType& helper, const ModelType& model)
	    : helper_(helper)
	    , model_(model)
	{ }

	template <typename ApplyOperatorType>
	FieldType operator()(SizeType                                        ptr,
	                     const typename ApplyOperatorType::OperatorType& A,
	                     SizeType                                        site,
	                     typename ApplyOperatorType::BorderEnum          corner,
	                     const PsimagLite::GetBraOrKet&                  bra,
	                     const PsimagLite::GetBraOrKet&                  ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);

			return onePointInternal<ApplyOperatorType>(
			    A, site, src1, src2, corner, ptr);
		} catch (std::exception& e) {
			std::cerr << "CAUGHT: " << e.what();
			std::cerr << "WARNING: Observer::onePoint(...): Nothing here yet\n";
			return 0;
		}
	}

	template <typename ApplyOperatorType>
	FieldType hookForZero(SizeType                                        ptr,
	                      const typename ApplyOperatorType::OperatorType& A,
	                      SizeType                                        splitSize,
	                      const PsimagLite::GetBraOrKet&                  bra,
	                      const PsimagLite::GetBraOrKet&                  ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);

			return onePointInternalHookForZero<ApplyOperatorType>(
			    A, splitSize, src1, src2, ptr);
		} catch (std::exception& e) {
			std::cerr << "CAUGHT: " << e.what();
			std::cerr << "WARNING: Observer::onePoint(...): Nothing here yet\n";
			return 0;
		}
	}

private:

	template <typename ApplyOperatorType>
	FieldType onePointInternal(const typename ApplyOperatorType::OperatorType& A,
	                           SizeType                                        site,
	                           const VectorWithOffsetType&                     src1,
	                           const VectorWithOffsetType&                     src2,
	                           typename ApplyOperatorType::BorderEnum          corner,
	                           SizeType                                        ptr) const
	{
		if (src1.sectors() == 0 || src2.sectors() == 0)
			return 0.0;

		SizeType             splitSize = model_.hilbertSize(site);
		ApplyOperatorType    applyOpLocal1(helper_.leftRightSuper(ptr), true);
		VectorWithOffsetType dest;
		applyOpLocal1(dest,
		              src1,
		              A,
		              helper_.fermionicSignLeft(ptr),
		              splitSize,
		              helper_.direction(ptr),
		              corner);

		FieldType                   sum = static_cast<FieldType>(0.0);
		const VectorWithOffsetType& v1  = dest;
		const VectorWithOffsetType& v2  = src2;
		OffsetVectorType            v1proxy;
		const auto&                 v1opt = v1proxy.makeOffsetVector(v1);
		OffsetVectorType            v2proxy;
		const auto&                 v2opt = v2proxy.makeOffsetVector(v2);

		for (SizeType ii = 0; ii < v1.sectors(); ii++) {
			SizeType i = v1.sector(ii);
			for (SizeType jj = 0; jj < v1.sectors(); jj++) {
				SizeType j = v2.sector(jj);
				if (i != j)
					continue;
				SizeType offset = v1.offset(i);
				for (SizeType k = 0; k < v1.effectiveSize(i); k++)
					sum += v1opt.slowAccess(k + offset)
					    * PsimagLite::conj(v2opt.slowAccess(k + offset));
			}
		}

		return sum;
	}

	template <typename ApplyOperatorType>
	FieldType onePointInternalHookForZero(const typename ApplyOperatorType::OperatorType& A,
	                                      SizeType                    splitSize,
	                                      const VectorWithOffsetType& src1,
	                                      const VectorWithOffsetType& src2,
	                                      SizeType                    ptr) const
	{

		ApplyOperatorType    applyOpLocal1(helper_.leftRightSuper(ptr), true);
		VectorWithOffsetType dest;
		applyOpLocal1.hookForZero(dest,
		                          src1,
		                          A,
		                          splitSize,
		                          // helper_.fermionicSignLeft(ptr),
		                          helper_.direction(ptr));

		FieldType                   sum = static_cast<FieldType>(0.0);
		const VectorWithOffsetType& v1  = dest;
		const VectorWithOffsetType& v2  = src2;

		for (SizeType ii = 0; ii < v1.sectors(); ii++) {
			SizeType i = v1.sector(ii);
			for (SizeType jj = 0; jj < v2.sectors(); jj++) {
				SizeType j = v2.sector(jj);
				if (i != j)
					continue;
				assert(v1.offset(i) == v2.offset(j));
				for (SizeType k = 0; k < v1.effectiveSize(i); k++)
					sum += v1.fastAccess(i, k)
					    * PsimagLite::conj(v2.fastAccess(j, k));
			}
		}

		return sum;
	}

	const ObserverHelperType& helper_;
	const ModelType&          model_;
}; // class OnePointCorrelations
} // namespace Dmrg

/*@}*/
#endif // ONE_POINT_H
