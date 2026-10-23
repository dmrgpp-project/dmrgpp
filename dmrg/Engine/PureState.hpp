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
	PureState(const WaveFunctionTransfType& wft,
	          SizeType hilbert)
	    : wft_(wft)
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
	void getNewPures(VectorWithOffsetsType& dest, SizeType site0, SizeType site1, const LrsType& lrs)
	{
		SizeType alpha_fixed = 0; // FIXME
		SizeType beta_fixed = 0; // FIXME

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "New pures for site=" << site0 << " is " << alpha_fixed;
		msg << " and for site=" << site1 << " is " << beta_fixed;
		progress_.printline(msgg, std::cerr);

		       // Pure for system
		{
			const BlockDiagMatrixType& transform_sys = getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM);
			VectorType new_vec1(transform_sys.rows());

			getNewPure(new_vec1,
			           pure_vectors_.first,
			           ProgramGlobals::SysOrEnvEnum::SYSTEM,
			           hilbert_,
			           lrs.left(),
			           transform_sys,
			           site0);
			pure_vectors_.first = std::move(new_vec1);
		}

		{
			const BlockDiagMatrixType& transform_env = getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON);
			VectorType new_vec2(transform_env.rows());

			getNewPure(new_vec2,
			           pure_vectors_.second,
			           ProgramGlobals::SysOrEnvEnum::ENVIRON,
			           hilbert_,
			           lrs.right(),
			           transform_env,
			           site1);
			pure_vectors_.second = std::move(new_vec2);
		}

		setFromInfinite(dest, lrs);
		assert(norm(dest) > 1e-6);

		// systemPrev_.fixed               = alphaFixedVolume;
		// systemPrev_.permutationInverse  = lrs_.left().permutationInverse();
		// environPrev_.fixed              = betaFixedVolume;
		// environPrev_.permutationInverse = lrs_.right().permutationInverse();
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

	// The pure_vectors.{first,second} are regular std::vector;
	// they need to be converted into vector with offsets;
	// this function does it.
	void setFromInfinite(VectorWithOffsetsType& phi, const LrsType& lrs) const
	{
		phi.populateSectors(lrs.super());
		for (SizeType ii = 0; ii < phi.sectors(); ii++) {
			SizeType         i0 = phi.sector(ii);
			VectorType v;
			getFullVector(v, i0, lrs);
			RealType tmpNorm = PsimagLite::norm(v);
			if (std::fabs(tmpNorm - 1.0) < 1e-6) {
				const Qn& j = lrs.super().qnEx(i0);
				std::cerr << "setFromInfinite: qns= ";
				std::cerr << j << "\n";
			}

			phi.setDataInSector(v, i0);
		}

		phi.collapseSectors();
		assert(norm(phi) > 1e-6);
	}

	void getFullVector(VectorType& v, SizeType m, const LrsType& lrs) const
	{
		int offset = lrs.super().partition(m);
		int total  = lrs.super().partition(m + 1) - offset;

		PsimagLite::PackIndices pack(lrs.left().size());
		v.resize(total);
		assert(PsimagLite::norm(pure_vectors_.first) > 1e-6);
		assert(PsimagLite::norm(pure_vectors_.second) > 1e-6);
		for (int i = 0; i < total; i++) {
			SizeType alpha, beta;
			pack.unpack(alpha, beta, lrs.super().permutation(i + offset));
			v[i] = pure_vectors_.first[alpha] * pure_vectors_.second[beta];
		}
	}


	const WaveFunctionTransfType& wft_;
	SizeType hilbert_;
	PsimagLite::ProgressIndicator progress_;
	std::pair<VectorType, VectorType> pure_vectors_;
	VectorPtrBlockDiagMatrixType                 garbage_;

};
}
#endif // PURESTATE_HPP
