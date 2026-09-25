// SPDX-FileCopyrightText: Copyright (c) 2009-2014-2020, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file LeftRightSuper.h
 *
 *  A class that contains the left block or system, the
 *  right block or environ, and the superblock
 */
#ifndef LEFT_RIGHT_SUPER_H
#define LEFT_RIGHT_SUPER_H

#include "BasisTraits.hh"
#include "ProgramGlobals.h"
#include <PsimagLite/Io/IoNg.h>
#include <PsimagLite/ProgressIndicator.h>

namespace Dmrg {

template <typename BasisWithOperatorsType_, typename SuperBlockType_> class LeftRightSuper {

public:

	using SuperBlockType         = SuperBlockType_;
	using RealType               = typename SuperBlockType::RealType;
	using BasisWithOperatorsType = BasisWithOperatorsType_;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using OperatorsType          = typename BasisWithOperatorsType::OperatorsType;
	using OperatorType           = typename OperatorsType::OperatorType;
	using OperatorStorageType    = typename OperatorType::StorageType;
	using BlockType              = typename BasisType::BlockType;
	using ProgressIndicatorType  = PsimagLite::ProgressIndicator;
	using ThisType               = LeftRightSuper<BasisWithOperatorsType_, SuperBlockType>;
	using QnType                 = typename BasisType::QnType;

	template <typename IoInputter>
	LeftRightSuper(
	    IoInputter&        io,
	    PsimagLite::String prefix,
	    const BasisTraits& basisTraits,
	    typename PsimagLite::EnableIf<PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	    : progress_("LeftRightSuper")
	    , left_(0)
	    , right_(0)
	    , super_(0)
	    , refCounter_(0)
	{
		prefix += "/LRS";

		PsimagLite::String nameSuper;
		io.read(nameSuper, prefix + "/NameSuper");

		PsimagLite::String nameLeft;
		io.read(nameLeft, prefix + "/NameSystem");

		PsimagLite::String nameRight;
		io.read(nameRight, prefix + "/NameEnviron");

		BasisTraits basisTraits2   = basisTraits;
		basisTraits2.isObserveCode = true;
		super_ = new SuperBlockType(io, prefix + "/" + nameSuper, basisTraits);
		left_  = new BasisWithOperatorsType(io, prefix + "/" + nameLeft, basisTraits2);
		right_ = new BasisWithOperatorsType(io, prefix + "/" + nameRight, basisTraits2);
	}

	LeftRightSuper(const PsimagLite::String& slabel,
	               const PsimagLite::String& elabel,
	               const PsimagLite::String& selabel,
	               const BasisTraits&        basisTraits)
	    : progress_("LeftRightSuper")
	    , left_(0)
	    , right_(0)
	    , super_(0)
	    , refCounter_(0)
	{
		left_  = new BasisWithOperatorsType(slabel, basisTraits);
		right_ = new BasisWithOperatorsType(elabel, basisTraits);
		super_ = new SuperBlockType(selabel, basisTraits);
	}

	~LeftRightSuper()
	{
		if (refCounter_ > 0) {
			--refCounter_;
			return;
		}

		delete left_;
		left_ = 0;
		delete right_;
		right_ = 0;
		delete super_;
		super_ = 0;
	}

	LeftRightSuper(BasisWithOperatorsType& left,
	               BasisWithOperatorsType& right,
	               SuperBlockType&         super)
	    : progress_("LeftRightSuper")
	    , left_(&left)
	    , right_(&right)
	    , super_(&super)
	    , refCounter_(1)
	{ }

	LeftRightSuper(const ThisType& rls)
	    : progress_("LeftRightSuper")
	    , refCounter_(1)
	{
		left_  = rls.left_;
		right_ = rls.right_;
		super_ = rls.super_;
	}

	void dontCopyOperators(const ThisType& rls)
	{
		assert(left_);
		left_->dontCopyOperators(rls.left());
		assert(right_);
		right_->dontCopyOperators(rls.right());
		assert(super_);
		assert(rls.super_);
		*super_ = *rls.super_;
		if (refCounter_ > 0)
			--refCounter_;
	}

	template <typename SomeModelType>
	SizeType growLeftBlock(const SomeModelType&    model,
	                       BasisWithOperatorsType& pS,
	                       BlockType const&        X,
	                       RealType                time)
	{
		assert(left_);
		return grow(
		    *left_, model, pS, X, ProgramGlobals::DirectionEnum::EXPAND_SYSTEM, time);
	}

	template <typename SomeModelType>
	SizeType growRightBlock(const SomeModelType&    model,
	                        BasisWithOperatorsType& pE,
	                        BlockType const&        X,
	                        RealType                time)
	{
		assert(right_);
		return grow(
		    *right_, model, pE, X, ProgramGlobals::DirectionEnum::EXPAND_ENVIRON, time);
	}

	void printSizes(const PsimagLite::String& label, std::ostream& os) const
	{
		assert(left_);
		assert(right_);

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << label << ": left-block basis=" << left_->size();
		msg << ", right-block basis=" << right_->size();
		msg << " sites=" << left_->block().size() << "+";
		msg << right_->block().size();
		progress_.printline(msgg, os);
	}

	SizeType sites() const
	{
		assert(left_);
		assert(right_);
		return left_->block().size() + right_->block().size();
	}

	/*!PTEX_LABEL{setToProductLrs} */
	void setToProduct(SizeType initialSizeOfHashTable)
	{
		assert(left_);
		assert(right_);
		assert(super_);
		super_->setToProduct(*left_, *right_, initialSizeOfHashTable);
	}

	void write(PsimagLite::IoNg::Out&                    io,
	           PsimagLite::String                        prefix,
	           typename BasisWithOperatorsType::SaveEnum option,
	           bool                                      minimizeWrite) const
	{
		prefix += "/LRS";
		io.createGroup(prefix);

		io.write(super_->name(), prefix + "/NameSuper");
		io.write(left_->name(), prefix + "/NameSystem");
		io.write(right_->name(), prefix + "/NameEnviron");

		assert(left_);
		assert(right_);
		assert(super_);

		super_->write(io, PsimagLite::IoSerializer::NO_OVERWRITE, prefix, minimizeWrite);
		left_->write(io, PsimagLite::IoSerializer::NO_OVERWRITE, prefix, option);
		right_->write(io, PsimagLite::IoSerializer::NO_OVERWRITE, prefix, option);
	}

	const BasisWithOperatorsType& left() const
	{
		assert(left_);
		return *left_;
	}

	const BasisWithOperatorsType& right() const
	{
		assert(right_);
		return *right_;
	}

	BasisWithOperatorsType& leftNonConst()
	{
		assert(left_);
		return *left_;
	}

	BasisWithOperatorsType& rightNonConst()
	{
		assert(right_);
		return *right_;
	}

	const SuperBlockType& super() const
	{
		assert(super_);
		return *super_;
	}

	void left(const BasisWithOperatorsType& left)
	{
		if (refCounter_ > 0)
			err("LeftRightSuper::left(...): not the owner\n");
		assert(left_);
		*left_ = left; // deep copy
	}

	void right(const BasisWithOperatorsType& right)
	{
		if (refCounter_ > 0)
			err("LeftRightSuper::right(...): not the owner\n");
		assert(right_);
		*right_ = right; // deep copy
	}

	template <typename IoInputType>
	void
	read(IoInputType&       io,
	     PsimagLite::String prefix,
	     typename PsimagLite::EnableIf<PsimagLite::IsInputLike<IoInputType>::True, int>::Type
	     = 0)
	{
		prefix += "/LRS";

		PsimagLite::String nameSuper;
		io.read(nameSuper, prefix + "/NameSuper");
		PsimagLite::String nameSys;
		io.read(nameSys, prefix + "/NameSystem");
		PsimagLite::String nameEnviron;
		io.read(nameEnviron, prefix + "/NameEnviron");

		super_->read(io, prefix + "/" + nameSuper);
		left_->read(io, prefix + "/" + nameSys);
		right_->read(io, prefix + "/" + nameEnviron);
	}

private:

	//! add block X to basis pS and put the result in left_:
	/* PSIDOC LeftRightSuperGrow
	        Local operators are set for the basis in question with a call to
	        \cppClass{BasisWithOperators}'s member function \cppFunction{setOperators()}.
	        When adding sites to the system or environment the program does a
	        full outer product, i.e., it increases the size of all local operators.
	        This is performed by the call to \verb!setToProduct!
	        \verb!(pSprime,pS,Xbasis,dir,option)!
	        in the grow function, which actually calls \verb!pSprime.setToProduct!
	        \verb!(pS,xBasis,dir)!
	        This function also recalculates the Hamiltonian in the outer product
	        of (i) the previous system basis $pS$, and (ii) the basis $Xbasis$
	        corresponding to the site(s) that is (are) being added.
	        To do this, the Hamiltonian connection between the two parts
	        needs to be calculated and added, and this is done in the call to
	        \cppFunction{addHamiltonianConnection}.
	        Finally, the resulting dmrgBasis object for the outer product,
	        pSprime, is set to contain this full Hamiltonian with the call
	        to  \cppFunction{pSprime.setHamiltonian(matrix)}.
	        */
	template <typename SomeModelType>
	SizeType grow(BasisWithOperatorsType&       leftOrRight,
	              const SomeModelType&          model,
	              BasisWithOperatorsType&       pS,
	              const BlockType&              X,
	              ProgramGlobals::DirectionEnum dir,
	              RealType                      time)
	{
		BasisWithOperatorsType Xbasis("Xbasis", pS.traits());
		using LeftRightSuper2Type = LeftRightSuper<BasisWithOperatorsType, BasisType>;
		SizeType oneSiteTruncSize = Xbasis.setOneSite(X, model, time);

		assert(X.size() == 1);
		SizeType lastS = pS.block().size();
		assert(lastS > 0);
		model.superOpHelper().setToProduct(pS.block()[--lastS], X[0], dir);
		leftOrRight.setToProduct(pS, Xbasis, model.superOpHelper());

		SparseMatrixType matrix = leftOrRight.hamiltonian().getCRS();

		LeftRightSuper2Type* lrs;
		BasisType*           leftOrRightL = &leftOrRight;
		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			lrs = new LeftRightSuper2Type(pS, Xbasis, *leftOrRightL);
		} else {
			lrs = new LeftRightSuper2Type(Xbasis, pS, *leftOrRightL);
		}
		//! PTEX_LABEL{295}
		model.addHamiltonianConnection(matrix, *lrs, time);
		std::vector<OperatorType> nonLocalOps;
		model.fillNewNonLocals(nonLocalOps, *lrs, time);
		delete lrs;
		leftOrRight.setHamiltonian(matrix);
		leftOrRight.addNewNonLocals(nonLocalOps);

		return oneSiteTruncSize;
	}

	LeftRightSuper(LeftRightSuper&);

	LeftRightSuper& operator=(const LeftRightSuper&);

	ProgressIndicatorType   progress_;
	BasisWithOperatorsType* left_;
	BasisWithOperatorsType* right_;
	SuperBlockType*         super_;
	SizeType                refCounter_;
}; // class LeftRightSuper

} // namespace Dmrg

/*@}*/
#endif // LEFT_RIGHT_SUPER_H
