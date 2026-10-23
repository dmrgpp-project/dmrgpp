#ifndef PURESTATE_HPP
#define PURESTATE_HPP

#include "LeftRightSuper.h"
#include "Basis.h"
#include "BasisWithOperators.h"
#include "BlockDiagonalMatrix.h"
#include "Wft/WaveFunctionTransfFactory.h"
#include "OneSiteSpaces.hh"
#include "VectorWithOffsets.h"
#include "DmrgOptions.h"
#include "ModelBase.h"
#include "ModelHelperLocal.h"
#include "ParametersDmrgSolver.h"
#include "SuperGeometry.h"
#include "ProgramGlobals.h"
#include "InputCheck.h"

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/PackIndices.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/Matrix.h>

#include <vector>
#include <cassert>

namespace Dmrg {

template<typename ComplexOrRealType>
class PureState {
public:

	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorType = std::vector<ComplexOrRealType>;
	using CrsMatrixType = PsimagLite::CrsMatrix<ComplexOrRealType>;
	using BasisType = Basis<CrsMatrixType>;
	using BasisWithOperatorsType = BasisWithOperators<BasisType>;
	using LrsType = LeftRightSuper<BasisWithOperatorsType, BasisType>;
	using MatrixType              = PsimagLite::Matrix<ComplexOrRealType>;
	using BlockDiagMatrixType = BlockDiagonalMatrix<MatrixType>;
	using VectorPtrBlockDiagMatrixType = std::vector<BlockDiagMatrixType*>;
	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType>;
	using InputType = PsimagLite::InputNg<InputCheck>::Readable;
	using OptionsType           = DmrgOptions<InputType>;
	using ModelHelperType = ModelHelperLocal<LrsType>;
	using ParametersDmrgType = ParametersDmrgSolver<RealType, InputType, Qn>;
	using ProgramGlobalsType = ProgramGlobals;
	using SuperGeometryType = SuperGeometry<ComplexOrRealType, InputType, ProgramGlobalsType>;
	using ModelType = ModelBase<ModelHelperType,
	                                               ParametersDmrgType,
	                                               InputType,
	                                               SuperGeometryType>;
	using OneSiteSpacesType = OneSiteSpaces<ModelType>;
	using WaveFunctionTransfType  = WaveFunctionTransfFactory<LrsType,
	                                                               VectorWithOffsetsType,
	                                                               OptionsType,
	                                                               OneSiteSpacesType>;

	struct MettsPrev {
		SizeType       fixed = 0;
		std::vector<SizeType> perm_inverse;
	};

	PureState(const WaveFunctionTransfType& wft,
	          const std::vector<SizeType>& pure_indices,
	          SizeType hilbert)
	    : wft_(wft)
	    , pure_indices_(pure_indices)
	    , hilbert_(hilbert)
	    , progress_("PureState")
	{}

	~PureState()
	{
		SizeType n = garbage_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete garbage_[i];
			garbage_[i] = nullptr;
		}
	}

	// direction here is INFINITE
	void getNewPures(VectorType& dest, SizeType i0, SizeType site0, SizeType site1, const LrsType& lrs)
	{
		assert(site0 < pure_indices_.size());
		assert(site1 < pure_indices_.size());

		SizeType alpha_fixed = pure_indices_[site0];
		SizeType beta_fixed = pure_indices_[site1];

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "New pures for site=" << site0 << " is " << alpha_fixed;
		msg << " and for site=" << site1 << " is " << beta_fixed;
		progress_.printline(msgg, std::cerr);

		// Pure for System
		{
			const BlockDiagMatrixType& transform_sys = getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM);
			VectorType new_vec1(transform_sys.rows()); // <--- FIXME: avoid allocation?

			getNewPure(new_vec1,
			           pure_vectors_.first,
			           ProgramGlobals::SysOrEnvEnum::SYSTEM,
			           hilbert_,
			           lrs.left(),
			           transform_sys,
			           site0);
			pure_vectors_.first = std::move(new_vec1);
		}

		// Pure for Environ
		{
			const BlockDiagMatrixType& transform_env = getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON);
			VectorType new_vec2(transform_env.rows()); // <--- FIXME: avoid allocation?

			getNewPure(new_vec2,
			           pure_vectors_.second,
			           ProgramGlobals::SysOrEnvEnum::ENVIRON,
			           hilbert_,
			           lrs.right(),
			           transform_env,
			           site1);
			pure_vectors_.second = std::move(new_vec2);
		}

		// Pure for Superblock = Environ x System
		getFullVector(dest, i0, lrs);

		// Prepare for next step
		sys_previous_.fixed               = alpha_fixed;
		sys_previous_.perm_inverse  = lrs.left().permutationInverse();
		env_previous_.fixed              = beta_fixed;
		env_previous_.perm_inverse = lrs.right().permutationInverse();
	}

private:

	const BlockDiagMatrixType& getTransform(ProgramGlobals::SysOrEnvEnum sysOrEnv)
	{
		const SizeType stackSize = wft_.size(sysOrEnv);

		if (stackSize > 0)
			return wft_.getTransform(sysOrEnv);

		BlockDiagMatrixType* m = new BlockDiagMatrixType;
		garbage_.push_back(m);

		return *m;
	}

	       // do the tensor product of pure_vectors_.first * pure_vectors_.second
	void getFullVector(VectorType& v, SizeType m, const LrsType& lrs) const
	{
		int offset = lrs.super().partition(m);
		int total  = lrs.super().partition(m + 1) - offset;

		PsimagLite::PackIndices pack(lrs.left().size());
		v.resize(total);
		assert(PsimagLite::norm(pure_vectors_.first) > 1e-6);
		assert(PsimagLite::norm(pure_vectors_.second) > 1e-6);
		for (int i = 0; i < total; i++) {
			SizeType alpha = 0;
			SizeType beta = 0;
			pack.unpack(alpha, beta, lrs.super().permutation(i + offset));
			v[i] = pure_vectors_.first[alpha] * pure_vectors_.second[beta];
		}
	}

	void getNewPure(VectorType&                  newVector,
	                VectorType&                  oldVector,
	                const ProgramGlobals::SysOrEnvEnum direction,
	                SizeType                           alphaFixed,
	                const BasisWithOperatorsType&      basis,
	                const BlockDiagMatrixType&     transform,
	                SizeType site)
	{
		if (oldVector.size() == 0)
			setInitialPure(oldVector, site);

		VectorType tmpVector;
		if (transform.rows() == 0) {
			tmpVector = oldVector;
			assert(PsimagLite::norm(tmpVector) > 1e-6);
		} else {
			MatrixType transform1;
			transform.toDense(transform1);
			delayedTransform(tmpVector, oldVector, direction, transform1);
			assert(PsimagLite::norm(tmpVector) > 1e-6);
		}

		SizeType       ns = tmpVector.size();
		SizeType newSize
		    = (transform.cols() == 0) ? (ns * ns) : transform.cols() * hilbert_;

		newVector.resize(newSize); // <--- note that new vector alway allocated here
		std::fill(newVector.begin(), newVector.end(), 0);

		for (SizeType alpha = 0; alpha < ns; ++alpha) {
			SizeType gamma   = (direction == ProgramGlobals::SysOrEnvEnum::SYSTEM)
			? basis.permutationInverse(alpha + alphaFixed * ns)
			: basis.permutationInverse(alphaFixed + alpha * hilbert_);
			newVector[gamma] = tmpVector[alpha];
		}

		PsimagLite::OstringStream                     msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2 << "Old size of pure is " << ns << " norm=" << PsimagLite::norm(tmpVector);
		progress_.printline(msgg2, std::cerr);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "New size of pure is " << newSize << " norm=" << PsimagLite::norm(newVector);
		progress_.printline(msgg, std::cerr);
		assert(PsimagLite::norm(newVector) > 1e-6);
	}

	void delayedTransform(VectorType&                  newVector,
	                      VectorType&                  oldVector,
	                      const ProgramGlobals::SysOrEnvEnum direction,
	                      const MatrixType&                  transform)
	{
		assert(oldVector.size() == transform.rows());

		SizeType ne = hilbert_;

		const std::vector<SizeType>& permutationInverse
		    = (direction == ProgramGlobals::SysOrEnvEnum::SYSTEM)
		    ? sys_previous_.perm_inverse
		    : env_previous_.perm_inverse;
		SizeType nsPrev = permutationInverse.size() / ne;

		assert(newVector.size() == 0);
		newVector.resize(transform.cols());

		// newVector = oldVector * transform;
		for (SizeType gamma = 0; gamma < newVector.size(); ++gamma) {
			newVector[gamma] = 0;
			for (SizeType alpha = 0; alpha < nsPrev; ++alpha) {
				SizeType noPermIndex
				    = (direction == ProgramGlobals::SysOrEnvEnum::SYSTEM)
				    ? alpha + sys_previous_.fixed * nsPrev
				    : env_previous_.fixed + alpha * ne;

				SizeType gammaPrime = permutationInverse[noPermIndex];

				assert(gammaPrime < transform.rows());
				newVector[gamma]
				    += transform(gammaPrime, gamma) * oldVector[gammaPrime];
			}
		}
	}

	void setInitialPure(VectorType& oldVector, SizeType site)
	{
//		site = (site == 1) ? 0 : site + 1; // <--- FIXME???

		assert(site < pure_indices_.size());
		SizeType alphaFixed = pure_indices_[site];

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "New pures for site "<<site<<" ";
		msg << " is " << alphaFixed;
		progress_.printline(msgg, std::cerr);

		oldVector.resize(hilbert_);
		assert(alphaFixed < oldVector.size());
		for (SizeType i = 0; i < oldVector.size(); i++) {
			oldVector[i] = (i == alphaFixed) ? 1 : 0;
		}

		assert(PsimagLite::norm(oldVector) > 1e-6);
	}

	const WaveFunctionTransfType& wft_;
	std::vector<SizeType> pure_indices_;
	SizeType hilbert_;
	PsimagLite::ProgressIndicator progress_;
	std::pair<VectorType, VectorType> pure_vectors_;
	MettsPrev sys_previous_;
	MettsPrev env_previous_;
	VectorPtrBlockDiagMatrixType                 garbage_;

};
}
#endif // PURESTATE_HPP
