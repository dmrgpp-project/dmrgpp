// SPDX-FileCopyrightText: Copyright (c) 2009-2012-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file WaveFunctionTransfFactory.h
 *
 *  This class implements the wave function transformation factory,
 *  see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_FACTORY_H
#define WFT_FACTORY_H
#include "Utils.h"
#include "WaveFunctionTransfLocal.h"
#include "WaveStructCombined.h"
#include <PsimagLite/Io/IoSelector.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/Random48.h>

namespace Dmrg {
template <typename LeftRightSuperType_,
          typename VectorWithOffsetType_,
          typename OptionsType_,
          typename OneSiteSpacesType_>
class WaveFunctionTransfFactory {

	using IoType = PsimagLite::IoSelector;

public:

	using VectorSizeType             = PsimagLite::Vector<SizeType>::Type;
	using VectorStringType           = PsimagLite::Vector<PsimagLite::String>::Type;
	using LeftRightSuperType         = LeftRightSuperType_;
	using BasisWithOperatorsType     = typename LeftRightSuperType::BasisWithOperatorsType;
	using BlockDiagonalMatrixType    = typename BasisWithOperatorsType::BlockDiagonalMatrixType;
	using SparseMatrixType           = typename BasisWithOperatorsType::SparseMatrixType;
	using BasisType                  = typename BasisWithOperatorsType::BasisType;
	using SparseElementType          = typename SparseMatrixType::value_type;
	using VectorType                 = typename PsimagLite::Vector<SparseElementType>::Type;
	using RealType                   = typename BasisWithOperatorsType::RealType;
	using WaveStructCombinedType     = WaveStructCombined<LeftRightSuperType>;
	using VectorVectorRealType       = typename WaveStructCombinedType::VectorVectorRealType;
	using VectorMatrixType           = typename WaveStructCombinedType::VectorMatrixType;
	using VectorQnType               = typename WaveStructCombinedType::VectorQnType;
	using VectorWithOffsetType       = VectorWithOffsetType_;
	using OneSiteSpacesType          = OneSiteSpacesType_;
	using WaveFunctionTransfBaseType = WaveFunctionTransfBase<WaveStructCombinedType,
	                                                          VectorWithOffsetType,
	                                                          OptionsType_,
	                                                          OneSiteSpacesType_>;
	using WaveFunctionTransfLocalType = WaveFunctionTransfLocal<WaveStructCombinedType,
	                                                            VectorWithOffsetType,
	                                                            OptionsType_,
	                                                            OneSiteSpacesType_>;
	using WftOptionsType              = typename WaveFunctionTransfBaseType::WftOptionsType;
	using WaveStructSvdType           = typename WaveStructCombinedType::WaveStructSvdType;

	template <typename SomeParametersType>
	WaveFunctionTransfFactory(SomeParametersType& params)
	    : isEnabled_(!(params.options.isSet("nowft")))
	    , wftOptions_(ProgramGlobals::DirectionEnum::INFINITE,
	                  params.options,
	                  true,
	                  true,
	                  params.denseSparseThreshold,
	                  params.gemmRnb,
	                  params.nthreads2)
	    , progress_("WaveFunctionTransf")
	    , filenameIn_(params.checkpoint.filename())
	    , filenameOut_(params.filename)
	    , waveStructCombined_(
	          params.options.isSet("wftstacksondisk"),
	          params.filename,
	          { params.options.isSet("observe"), params.options.isSet("noSaveOperators") })
	    , wftImpl_(0)
	    , rng_(3433117)
	    , noLoad_(false)
	    , save_(!params.options.isSet("noSaveWft") && !params.options.isSet("minimizeDisk"))
	{
		if (!isEnabled_)
			return;

		bool b = (params.options.isSet("restart") || params.autoRestart);

		if (b) {
			if (params.options.isSet("noloadwft"))
				noLoad_ = true;
			else
				read();
		} else {
			if (params.options.isSet("noloadwft")) {
				PsimagLite::String str(
				    "Error: noloadwft needs restart or checkpoint\n");
				throw PsimagLite::RuntimeError(str.c_str());
			}
		}

		wftImpl_ = new WaveFunctionTransfLocalType(waveStructCombined_, wftOptions_);
	}

	~WaveFunctionTransfFactory()
	{
		if (!isEnabled_)
			return;

		if (save_) {
			IoType::Out ioOut(filenameOut_, IoType::ACC_RDW);
			write(ioOut);
		}

		delete wftImpl_;
	}

	void setStage(ProgramGlobals::DirectionEnum stage)
	{
		if (stage == wftOptions_.dir)
			return;
		wftOptions_.dir    = stage;
		wftOptions_.bounce = true;
	}

	void triggerOn()
	{
		bool allow = false;
		switch (wftOptions_.dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			allow = false;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			allow = true;
			break;
		}

		// FIXME: Must check the below change when using SU(2)!!
		// if (m<0) allow = false; // isEnabled_=false;

		if (noLoad_)
			allow = false;

		if (!isEnabled_ || !allow)
			return;

		waveStructCombined_.beforeWft(
		    wftOptions_.dir, wftOptions_.twoSiteDmrg, wftOptions_.bounce);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Window open, ready to transform vectors";
		progress_.printline(msgg, std::cout);
	}

	// FIXME: change name to transformVector
	void setInitialVector(VectorWithOffsetType&       dest,
	                      const VectorWithOffsetType& src,
	                      const LeftRightSuperType&   lrs,
	                      const OneSiteSpacesType&    oneSiteSpaces) const
	{
		bool allow = false;
		switch (wftOptions_.dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			allow = false;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
			allow = true;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			allow = true;
			break;
		}

		// FIXME: Must check the below change when using SU(2)!!
		// if (m<0) allow = false; // isEnabled_=false;

		if (noLoad_)
			allow = false;

		if (isEnabled_ && allow) {
			RealType eps = 1e-12;
			RealType x   = norm(src);
			bool     b   = (x < eps);
			if (b) {
				err(std::string("An important vector has norm=") + ttos(x)
				    + " which is too small\n");
			}

			createVector(dest, src, lrs, oneSiteSpaces);

		} else {
			createRandomVector(dest);
		}
	}

	void triggerOff(const LeftRightSuperType& lrs)
	{
		bool allow = false;
		switch (wftOptions_.dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			allow = false;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
			allow = true;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			allow = true;
			break;
		}

		// FIXME: Must check the below change when using SU(2)!!
		// if (m<0) allow = false; // isEnabled_=false;

		if (noLoad_)
			allow = false;

		if (!isEnabled_ || !allow)
			return;
		afterWft(lrs);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Window closed, no more transformations, please";
		progress_.printline(msgg, std::cout);
	}

	void createRandomVector(VectorWithOffsetType& y) const
	{
		for (SizeType jj = 0; jj < y.sectors(); jj++) {
			SizeType j = y.sector(jj);
			createRandomVector(y, j);
		}

		if (!isEnabled_)
			return; // don't make noise unless enabled
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Yes, I'm awake, but there's nothing heavy to do now";
		progress_.printline(msgg, std::cout);
	}

	void createRandomVector(VectorWithOffsetType& y, SizeType i0) const
	{
		SizeType                                  total = y.effectiveSize(i0);
		typename VectorWithOffsetType::value_type tmp;
		RealType                                  atmp = 0;
		for (SizeType i = 0; i < total; i++) {
			myRandomT(tmp);
			y.fastAccess(i0, i) = tmp;
			atmp += PsimagLite::real(tmp * PsimagLite::conj(tmp));
		}

		assert(fabs(atmp) > 1e-10);
		atmp = 1.0 / sqrt(atmp);
		for (SizeType i = 0; i < total; i++)
			y.fastAccess(i0, i) *= atmp;
	}

	void push(const BlockDiagonalMatrixType& transform,
	          ProgramGlobals::DirectionEnum  direction,
	          const LeftRightSuperType&      lrs,
	          const VectorMatrixType&        vts,
	          const VectorVectorRealType&    s,
	          const VectorQnType&            qns)
	{
		if (!isEnabled_)
			return;

		waveStructCombined_.push(transform, direction, vts, s, qns, wftOptions_.dir);
		waveStructCombined_.setLrs(lrs);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "OK, pushing option=" << ProgramGlobals::toString(direction);
		msg << " and stage=" << ProgramGlobals::toString(wftOptions_.dir);
		progress_.printline(msgg, std::cout);

		if (noLoad_) {
			SizeType center = computeCenter(lrs, direction);
			updateNoLoad(lrs, center);
		}
	}

	const BlockDiagonalMatrixType& getTransform(ProgramGlobals::SysOrEnvEnum dir) const
	{
		return waveStructCombined_.getTransform(dir);
	}

	SizeType size(ProgramGlobals::SysOrEnvEnum sysOrEnv) const
	{
		return waveStructCombined_.size(sysOrEnv);
	}

	bool isEnabled() const { return isEnabled_; }

	const WftOptionsType options() const { return wftOptions_; }

	void write(PsimagLite::IoSelector::Out& ioMain)
	{
		if (!isEnabled_)
			return;
		if (!save_)
			return;

		PsimagLite::String label = "Wft";
		writePartial(ioMain, label);
		waveStructCombined_.write(ioMain, label + "/WaveStructCombined");
	}

	void write(PsimagLite::IoSelector::Out& ioMain) const
	{
		if (!isEnabled_)
			return;
		if (!save_)
			return;

		PsimagLite::String label = "Wft";
		writePartial(ioMain, label);
		waveStructCombined_.write(ioMain, label + "/WaveStructCombined");
	}

	const BlockDiagonalMatrixType&
	multiPointGetTransform(SizeType ind, ProgramGlobals::DirectionEnum dir) const
	{
		return waveStructCombined_.multiPointGetTransform(ind, dir);
	}

private:

	void writePartial(PsimagLite::IoSelector::Out& ioMain, PsimagLite::String prefix) const
	{
		assert(isEnabled_);
		assert(save_);

		ioMain.createGroup(prefix);
		ioMain.write(isEnabled_, prefix + "/isEnabled");
		wftOptions_.write(ioMain, prefix + "/WftOptions");
	}

	void read()
	{
		if (!isEnabled_)
			throw PsimagLite::RuntimeError(
			    "WFT::read(...) called but wft is disabled\n");

		PsimagLite::IoSelector::In ioMain(filenameIn_);
		PsimagLite::String         label = "Wft";
		ioMain.read(isEnabled_, label + "/isEnabled");
		wftOptions_.read(ioMain, label + "/WftOptions");
		waveStructCombined_.read(ioMain, label + "/WaveStructCombined");
		ioMain.close();
	}

	void myRandomT(std::complex<RealType>& value) const
	{
		value = std::complex<RealType>(rng_() - 0.5, rng_() - 0.5);
	}

	void myRandomT(RealType& value) const { value = rng_() - 0.5; }

	void afterWft(const LeftRightSuperType& lrs)
	{
		waveStructCombined_.setLrs(lrs);
		waveStructCombined_.afterWft(wftOptions_.dir);
		wftOptions_.firstCall = false;
		wftOptions_.bounce    = false;
	}

	void createVector(VectorWithOffsetType&       psiDest,
	                  const VectorWithOffsetType& psiSrc,
	                  const LeftRightSuperType&   lrs,
	                  const OneSiteSpacesType&    oneSiteSpaces) const
	{

		RealType norm1 = norm(psiSrc);
		if (norm1 < 1e-5)
			std::cerr << "WFT Factory: norm1 = " << norm1 << " < 1e-5\n";

		wftImpl_->transformVector(psiDest, psiSrc, lrs, oneSiteSpaces);

		RealType                                      norm2 = norm(psiDest);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Transformation completed ";
		if (fabs(norm1 - norm2) > 1e-5) {
			msg << "WARNING: orig. norm= " << norm1 << " resulting norm= " << norm2;
		}

		if (norm2 < 1e-5)
			std::cerr << "WFT Factory: norm2 = " << norm2 << " < 1e-5\n";

		progress_.printline(msgg, std::cout);
	}

	SizeType computeCenter(const LeftRightSuperType&     lrs,
	                       ProgramGlobals::DirectionEnum direction) const
	{
		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			SizeType total = lrs.left().block().size();
			assert(total > 0);
			total--;
			return lrs.left().block()[total];
		}

		return lrs.right().block()[0];
	}

	void updateNoLoad(const LeftRightSuperType& lrs, SizeType center)
	{
		sitesSeen_.push_back(center);
		SizeType numberOfSites = lrs.super().block().size();
		if (checkSites(numberOfSites)) {
			noLoad_ = false;
			PsimagLite::OstringStream                     msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg << " now available";
			progress_.printline(msgg, std::cout);
		}
	}

	bool checkSites(SizeType numberOfSites) const
	{
		assert(numberOfSites > 0);
		for (SizeType i = 1; i < numberOfSites - 1; i++) {
			bool seen = (std::find(sitesSeen_.begin(), sitesSeen_.end(), i)
			             != sitesSeen_.end());
			if (!seen)
				return false;
		}

		return true;
	}

	bool                           isEnabled_;
	WftOptionsType                 wftOptions_;
	PsimagLite::ProgressIndicator  progress_;
	PsimagLite::String             filenameIn_;
	PsimagLite::String             filenameOut_;
	WaveStructCombinedType         waveStructCombined_;
	WaveFunctionTransfBaseType*    wftImpl_;
	PsimagLite::Random48<RealType> rng_;
	bool                           noLoad_;
	const bool                     save_;
	VectorSizeType                 sitesSeen_;
}; // class WaveFunctionTransformation
} // namespace Dmrg

/*@}*/
#endif
