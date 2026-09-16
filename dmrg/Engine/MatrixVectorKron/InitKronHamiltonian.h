// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file InitKronHamiltonian.h
 *
 *
 */
#ifndef INITKRON_HAMILTONIAN_H
#define INITKRON_HAMILTONIAN_H
#include "InitKronBase.h"
#include "OpsForLink.hh"
#include "ProgramGlobals.h"
#include <PsimagLite/Profiling.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename ModelType_>
class InitKronHamiltonian : public InitKronBase<typename ModelType_::LeftRightSuperType> {

	using VectorBoolType = typename PsimagLite::Vector<bool>::Type;

public:

	using ModelType                  = ModelType_;
	using HamiltonianConnectionType  = typename ModelType::HamiltonianConnectionType;
	using ModelHelperType            = typename ModelType::ModelHelperType;
	using LeftRightSuperType         = typename ModelHelperType::LeftRightSuperType;
	using OperatorStorageType        = typename ModelHelperType::OperatorStorageType;
	using BasisType                  = typename LeftRightSuperType::BasisType;
	using BaseType                   = InitKronBase<LeftRightSuperType>;
	using SparseMatrixType           = typename ModelHelperType::SparseMatrixType;
	using LinkType                   = typename HamiltonianConnectionType::LinkType;
	using RealType                   = typename ModelHelperType::RealType;
	using ComplexOrRealType          = typename SparseMatrixType::value_type;
	using ArrayOfMatStructType       = typename BaseType::ArrayOfMatStructType;
	using GenIjPatchType             = typename ArrayOfMatStructType::GenIjPatchType;
	using VectorArrayOfMatStructType = typename PsimagLite::Vector<ArrayOfMatStructType*>::Type;
	using VectorType                 = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using VectorSizeType             = typename ArrayOfMatStructType::VectorSizeType;
	using SuperBlockType             = typename LeftRightSuperType::SuperBlockType;
	using ParametersDmrgSolverType   = typename ModelType_::ParametersType;

	InitKronHamiltonian(const ModelType&                     model,
	                    const HamiltonianConnectionType&     hc,
	                    const typename ModelHelperType::Aux& aux)
	    : BaseType(hc.modelHelper().leftRightSuper(),
	               aux.m(),
	               hc.modelHelper().quantumNumber(aux.m()),
	               model.params().denseSparseThreshold,
	               !model.params().options.isSet("KronNoUseLowerPart")
	                   && !model.params().options.isSet("BatchedGemm"))
	    , model_(model)
	    , hc_(hc)
	    , vstart_(BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size() + 1)
	    , offsetForPatches_(BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size() + 1)
	{
		addHlAndHr();

		{
			PsimagLite::Profiling profiling("convertXcYcArrays", std::cout);

			convertXcYcArrays();
		}

		BaseType::setUpVstart(vstart_, BaseType::NEW);
		assert(vstart_.size() > 0);
		SizeType nsize = vstart_[vstart_.size() - 1];
		assert(nsize > 0);
		yin_.resize(nsize, 0.0);
		xout_.resize(nsize, 0.0);
		BaseType::computeOffsets(offsetForPatches_, BaseType::NEW);
	}

	bool isWft() const { return false; }

	bool loadBalance() const { return model_.params().options.isSet("KronLoadBalance"); }

	SizeType gemmRnb() const { return model_.params().gemmRnb; }

	SizeType nthreads2() const { return model_.params().nthreads2; }

	// -------------------
	// copy vin(:) to yin(:)
	// -------------------
	void copyIn(const VectorType& vout, const VectorType& vin)
	{
		VectorType& xout = xout_;
		VectorType& yin  = yin_;

		const VectorSizeType& permInverse
		    = BaseType::lrs(BaseType::NEW).super().permutationInverse();
		const SparseMatrixType& leftH
		    = BaseType::lrs(BaseType::NEW).left().hamiltonian().getCRS();
		SizeType nl = leftH.rows();

		SizeType offset       = BaseType::offset(BaseType::NEW);
		SizeType npatches     = BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size();
		const BasisType& left = BaseType::lrs(BaseType::NEW).left();
		const BasisType& right = BaseType::lrs(BaseType::NEW).right();

		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {

			SizeType igroup
			    = BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup
			    = BaseType::patch(BaseType::NEW, GenIjPatchType::RIGHT)[ipatch];

			assert(left.partition(igroup + 1) >= left.partition(igroup));
			SizeType sizeLeft = left.partition(igroup + 1) - left.partition(igroup);

			assert(right.partition(jgroup + 1) >= right.partition(jgroup));
			SizeType sizeRight = right.partition(jgroup + 1) - right.partition(jgroup);

			SizeType left_offset  = left.partition(igroup);
			SizeType right_offset = right.partition(jgroup);

			for (SizeType ileft = 0; ileft < sizeLeft; ++ileft) {
				for (SizeType iright = 0; iright < sizeRight; ++iright) {

					SizeType i = ileft + left_offset;
					SizeType j = iright + right_offset;

					SizeType ij = i + j * nl;

					assert(i < nl);
					assert(j < BaseType::lrs(BaseType::NEW)
					               .right()
					               .hamiltonian()
					               .rows());

					assert(ij < permInverse.size());

					SizeType r = permInverse[ij];
					assert(
					    !((r < offset)
					      || (r >= (offset + BaseType::size(BaseType::NEW)))));

					SizeType ip
					    = vstart_[ipatch] + (iright + ileft * sizeRight);
					assert(ip < yin.size());

					assert((r >= offset) && ((r - offset) < vin.size()));
					yin[ip]  = vin[r - offset];
					xout[ip] = vout[r - offset];
				}
			}
		}
	}

	// -------------------
	// copy xout(:) to vout(:)
	// -------------------
	void copyOut(VectorType& vout) const { BaseType::copyOut(vout, xout_, vstart_); }

	const VectorType& yin() const { return yin_; }

	VectorType& xout() { return xout_; }

	const SizeType& offsetForPatches(typename BaseType::WhatBasisEnum, SizeType ind) const
	{
		assert(ind < offsetForPatches_.size());
		return offsetForPatches_[ind];
	}

	const ParametersDmrgSolverType& params() const { return model_.params(); }

private:

	void addHlAndHr()
	{
		const RealType             value = 1.0;
		const OperatorStorageType& aL
		    = hc_.modelHelper().leftRightSuper().left().hamiltonian();
		const OperatorStorageType& aR
		    = hc_.modelHelper().leftRightSuper().right().hamiltonian();
		identityL_.makeDiagonal(aL.rows(), value);
		identityR_.makeDiagonal(aR.rows(), value);

		BaseType::addOneConnection(
		    aL, identityR_, value, ProgramGlobals::FermionOrBosonEnum::BOSON);
		BaseType::addOneConnection(
		    identityL_, aR, value, ProgramGlobals::FermionOrBosonEnum::BOSON);
	}

	void convertXcYcArrays()
	{
		SizeType total = hc_.tasks();

		OpsForLink<LeftRightSuperType> opsForLink = hc_.opsForLink();
		for (SizeType ix = 0; ix < total; ++ix) {

			opsForLink.setPointer(ix);
			const OperatorStorageType& A     = opsForLink.A();
			const OperatorStorageType& B     = opsForLink.B();
			const LinkType&            link2 = opsForLink.link();

			if (link2.type == ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM) {
				LinkType link3 = link2;
				link3.type     = ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON;
				if (link3.fermionOrBoson
				    == ProgramGlobals::FermionOrBosonEnum::FERMION)
					link3.value *= -1.0;

				BaseType::addOneConnection(B, A, link3.value, link2.fermionOrBoson);
				continue;
			}

			BaseType::addOneConnection(A, B, link2.value, link2.fermionOrBoson);
		}
	}

	InitKronHamiltonian(const InitKronHamiltonian&);

	InitKronHamiltonian& operator=(const InitKronHamiltonian&);

	const ModelType&                 model_;
	const HamiltonianConnectionType& hc_;
	OperatorStorageType              identityL_;
	OperatorStorageType              identityR_;
	VectorSizeType                   vstart_;
	VectorType                       yin_;
	VectorType                       xout_;
	VectorSizeType                   offsetForPatches_;
};
} // namespace Dmrg

#endif // INITKRON_HAMILTONIAN_H
