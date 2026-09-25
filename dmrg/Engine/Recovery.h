// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file Recovery.h
 *
 *
 */

#ifndef DMRG_RECOVER_H
#define DMRG_RECOVER_H

#include "Checkpoint.h"
#include "OutputFileOrNot.h"
#include "ProgramGlobals.h"
#include <PsimagLite/Io/IoNg.h>
#include <PsimagLite/PredicateAwesome.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>
#include <dirent.h>
#include <fstream>
#include <sys/types.h>

namespace Dmrg {

template <typename CheckpointType, typename TargetingType> class Recovery {

	using ParametersType     = typename CheckpointType::ParametersType;
	using RecoveryStaticType = Recovery<ParametersType, int>;
	using ComplexOrRealType  = typename CheckpointType::ComplexOrRealType;

	struct OptionSpec {

		OptionSpec()
		    : keepFiles(false)
		    , maxFiles(10)
		{ }

		bool     keepFiles;
		SizeType maxFiles;
	};

	struct OpaqueRestart {

		OpaqueRestart()
		    : loopIndex(0)
		    , stepCurrent(0)
		{ }

		SizeType loopIndex;
		SizeType stepCurrent;
	};

public:

	using VectorStringType       = PsimagLite::Vector<PsimagLite::String>::Type;
	using BasisWithOperatorsType = typename CheckpointType::BasisWithOperatorsType;
	using WaveFunctionTransfType = typename CheckpointType::WaveFunctionTransfType;
	using IoType                 = typename CheckpointType::IoType;
	using VectorSizeType         = typename PsimagLite::Vector<SizeType>::Type;
	using DiskOrMemoryStackType  = typename CheckpointType::DiskOrMemoryStackType;
	using DiskStackType          = typename CheckpointType::DiskStackType;
	using VectorBlockType        = PsimagLite::Vector<VectorSizeType>::Type;
	using RealType               = typename CheckpointType::RealType;

	class SpecOptions {

	public:

		SpecOptions(bool& keepFiles, SizeType& maxFiles)
		    : keepFiles_(keepFiles)
		    , maxFiles_(maxFiles)
		{ }

		void operator()(PsimagLite::String str2)
		{
			if (str2.length() < 2)
				return;

			PsimagLite::String str = str2.substr(1, str2.length() - 1);

			if (str == "keep") {
				keepFiles_ = true;
				return;
			}

			if (str.length() < 3)
				dieWithError(str);

			if (str[0] == 'M' && str[1] == '=') {
				PsimagLite::String each = str.substr(2, str.length() - 2);
				maxFiles_               = atoi(each.c_str());
				std::cerr << "Recovery Max files= " << maxFiles_ << "\n";
				return;
			}
		}

	private:

		bool&     keepFiles_;
		SizeType& maxFiles_;
	};

	Recovery(const VectorBlockType&        siteIndices,
	         const CheckpointType&         checkpoint,
	         const WaveFunctionTransfType& wft,
	         const BasisWithOperatorsType& pS,
	         const BasisWithOperatorsType& pE)
	    : progress_("Recovery")
	    , predicateAwesome_(nullptr)
	    , siteIndices_(siteIndices)
	    , checkpoint_(checkpoint)
	    , wft_(wft)
	    , pS_(pS)
	    , pE_(pE)
	    , counter_(0)
	{
		procOptions();

		if (!checkpoint_.parameters().options.isSet("recoveryEnableRead"))
			return;

		if (!checkpoint_.parameters().autoRestart)
			return;

		readRecovery();

		VectorStringType parts;
		RecoveryStaticType::makeThreeParts(parts,
		                                   checkpoint_.parameters().checkpoint.filename());
		if (parts.size() == 3)
			counter_ = 1 + atoi(parts[1].c_str());
	}

	~Recovery()
	{
		for (SizeType i = 0; i < optionSpec_.maxFiles; ++i) {
			PsimagLite::String prefix(RecoveryStaticType::recoveryFilePrefix());
			prefix += ttos(i);
			PsimagLite::String savedName(prefix + checkpoint_.parameters().filename);
			if (optionSpec_.keepFiles)
				continue;
			unlink(savedName.c_str());
		}

		delete predicateAwesome_;
		predicateAwesome_ = nullptr;
	}

	SizeType indexOfFirstFiniteLoop() const { return opaqueRestart_.loopIndex; }

	SizeType stepCurrent(ProgramGlobals::DirectionEnum direction) const
	{
		return (checkpoint_.parameters().autoRestart) ? opaqueRestart_.stepCurrent
		                                              : nonRecoveryStepCurrent(direction);
	}

	bool byLoop(RealType loopIndex, RealType time, RealType loopLength) const
	{
		return predicateAwesome_->isTrue("%l", loopIndex, "%t", time, "%ll", loopLength);
	}

	void write(const TargetingType& psi,
	           SizeType             loopIndex,
	           SizeType             stepCurrent,
	           int                  lastSign,
	           OutputFileOrNot&     ioOutCurrent,
	           PsimagLite::String   inputBlob) const
	{
		PsimagLite::String prefix(RecoveryStaticType::recoveryFilePrefix());
		prefix += ttos(counter_++);
		PsimagLite::String savedName(prefix + checkpoint_.parameters().filename);
		ioOutCurrent.flush();

		// copyFile(savedName.c_str(), ioOutCurrent.filename());

		typename IoType::Out ioOut(savedName, IoType::ACC_TRUNC);

		PsimagLite::PsiBase64::Encode base64encode(inputBlob);
		ioOut.write(base64encode, "InputBase64Encoded");

		writeEnergies(ioOut, ioOutCurrent.filename());

		writeRecovery(ioOut, loopIndex, stepCurrent);

		checkpoint_.parameters().write("PARAMETERS", ioOut.serializer());

		// taken from end of finiteDmrgLoops
		checkpoint_.write(pS_, pE_, ioOut);
		ioOut.createGroup("FinalPsi");
		psi.write(siteIndices_[stepCurrent], ioOut, "FinalPsi");
		ioOut.write(lastSign, "LastLoopSign");
		ioOut.write(PsimagLite::IsComplexNumber<ComplexOrRealType>::True, "IsComplex");
		// wft dtor
		wft_.write(ioOut);

		ioOut.close();

		// checkpoint stacks
		checkpoint_.checkpointStacks(savedName);

		if (counter_ >= optionSpec_.maxFiles)
			counter_ = 0;
	}

private:

	/* PSIDOC RecoverySave
	  This is a comma-separated list options; whitespace isn't allowed.
	  Supported options are as follows.

	  no, which disables RecoverySave, and is assumed if RecoverySave= is
	  absent from the input file. In other words, RecoverySave is disabled by default.

	  keep, which keeps recovery files even if the run finishes normally.
	  By default, all recovery files are deleted if the run finishes normally.

	  M=n, where n is the maximum number of recovery files that will be saved, before
	  the oldest file is overwritten. Defaults to 10.

	  Any predicate awesome that can use the loop variable %l
	 */
	void procOptions()
	{
		PsimagLite::String str = checkpoint_.parameters().recoverySave;

		SpecOptions lambda(optionSpec_.keepFiles, optionSpec_.maxFiles);
		predicateAwesome_
		    = new PsimagLite::PredicateAwesome<SpecOptions>(str, ",", "&", &lambda);
	}

	static void dieWithError(PsimagLite::String str)
	{
		err("Syntax error for RecoverySave expression " + str + "\n");
	}

	static void copyFile(PsimagLite::String destName, PsimagLite::String sourceName)
	{
		std::ifstream source(sourceName.c_str(), std::ios::binary);
		std::ofstream dest(destName.c_str(), std::ios::binary);
		dest << source.rdbuf();
		source.close();
		dest.close();
	}

	void
	writeRecovery(typename IoType::Out& ioOut, SizeType loopIndex, SizeType stepCurrent) const
	{
		ioOut.createGroup("Recovery");

		ioOut.write(loopIndex, "Recovery/loopIndex");
		ioOut.write(stepCurrent, "Recovery/stepCurrent");
	}

	void readRecovery()
	{
		typename IoType::In ioIn2(checkpoint_.parameters().checkpoint.filename());

		ioIn2.read(opaqueRestart_.loopIndex, "Recovery/loopIndex");
		ioIn2.read(opaqueRestart_.stepCurrent, "Recovery/stepCurrent");
		ioIn2.close();
	}

	// set initial site to add to either system or environment:
	// this is a bit tricky and has been a source of endless bugs
	// basically we have pS on the left and pE on the right,
	// and we need to determine which site is to be added
	int nonRecoveryStepCurrent(ProgramGlobals::DirectionEnum direction) const
	{
		// all right, now we can get the actual site to add:
		SizeType       sitesPerBlock = checkpoint_.parameters().sitesPerBlock;
		VectorSizeType siteToAdd(sitesPerBlock);
		// left-most site of pE
		for (SizeType j = 0; j < sitesPerBlock; ++j)
			siteToAdd[j] = pE_.block()[j];

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) {
			// right-most site of pS
			for (SizeType j = 0; j < sitesPerBlock; ++j)
				siteToAdd[j] = pS_.block()[pS_.block().size() - 1 - j];
		}

		// now stepCurrent_ is such that sitesIndices_[stepCurrent_] = siteToAdd
		// so:
		int sc = PsimagLite::indexOrMinusOne(siteIndices_, siteToAdd);

		if (sc < 0)
			err("nonRecoveryStepCurrent(...): step current error\n");

		return sc; // phew!!, that's all folks, now bugs, go away!!
	}

	void writeEnergies(typename IoType::Out& ioOut, PsimagLite::String file) const
	{
		PsimagLite::String energyLabel
		    = checkpoint_.parameters().checkpoint.labelForEnergy();
		ioOut.flush();
		typename IoType::In ioIn(file);

		typename CheckpointType::VectorVectorRealType energies;
		CheckpointType::readEnergies(energies, energyLabel, ioIn);
		CheckpointType::writeEnergies(true, energyLabel, energies, ioOut);

		ioIn.close();
	}

	PsimagLite::ProgressIndicator              progress_;
	PsimagLite::PredicateAwesome<SpecOptions>* predicateAwesome_;
	OptionSpec                                 optionSpec_;
	OpaqueRestart                              opaqueRestart_;
	const VectorBlockType&                     siteIndices_;
	const CheckpointType&                      checkpoint_;
	const WaveFunctionTransfType&              wft_;
	const BasisWithOperatorsType&              pS_;
	const BasisWithOperatorsType&              pE_;
	mutable SizeType                           counter_;
}; // class Recovery

template <typename ParametersType> class Recovery<ParametersType, int> {

public:

	static PsimagLite::String recoveryFilePrefix() { return "Recovery"; }

	static void checkOptions(PsimagLite::String                          recoverySave,
	                         const typename ParametersType::OptionsType& options)
	{
		if (recoverySave == "" || recoverySave == "no")
			return;

		if (!options.isSet("minimizeDisk"))
			return;

		err("minimizeDisk cannot be used with RecoverySave\n");
	}

	// this function is called before the ctor
	static void autoRestart(ParametersType& params)
	{
		if (!params.options.isSet("recoveryEnableRead"))
			return;

		// params.filename must have been corrected already if necessary
		PsimagLite::String recoveryFile = getRecoveryFile(params.filename);
		if (recoveryFile == "")
			return;

		// *  add the line RestartFilename= pointing to the data file of the
		// run to be restarted.
		params.checkpoint.setFilename(recoveryFile);
		// params.checkRestart(params.filename, recoveryFile, params.options, "INTERNAL=");

		params.autoRestart = true;
	}

	static void makeThreeParts(std::vector<PsimagLite::String>& parts,
	                           PsimagLite::String               filename)
	{
		const PsimagLite::String prefix = recoveryFilePrefix();
		const SizeType           len    = prefix.length();
		if (filename.substr(0, len) != prefix)
			return;
		parts.push_back(prefix);

		PsimagLite::String buffer("");
		for (SizeType i = len; i < filename.length(); ++i) {
			if (isAdigit(filename[i]))
				buffer += filename[i];
			break;
		}

		if (buffer == "")
			return;

		parts.push_back(buffer);

		SizeType lastPartLen = filename.length() - buffer.length() - len;
		parts.push_back(filename.substr(len + buffer.length(), lastPartLen));
	}

private:

	static bool isAdigit(char c) { return (c > 47 && c < 58); }

	static void listFilesInDirectory(std::vector<PsimagLite::String>& files,
	                                 PsimagLite::String               path)
	{
		DIR*    dir = 0;
		dirent* ent = 0;
		if ((dir = opendir(path.c_str())) != 0) {
			while ((ent = readdir(dir)) != 0) {
				files.push_back(ent->d_name);
			}

			closedir(dir);
			return;
		}

		/* could not open directory */
		perror("");
	}

	static PsimagLite::String getRecoveryFile(PsimagLite::String filename)
	{
		const PsimagLite::String        prefix = recoveryFilePrefix();
		std::vector<PsimagLite::String> files;
		listFilesInDirectory(files, ".");

		if (files.size() == 0)
			return "";

		PsimagLite::String saved("");
		SizeType           max = 0;

		for (SizeType i = 0; i < files.size(); ++i) {
			std::vector<PsimagLite::String> parts;
			makeThreeParts(parts, files[i]);
			if (parts.size() != 3 || parts[0] != prefix || parts[2] != filename)
				continue;
			SizeType counter = atoi(parts[1].c_str());
			if (counter >= max && isValidFile(files[i])) {
				max   = counter;
				saved = files[i];
			}
		}

		return saved;
	}

	static bool isValidFile(PsimagLite::String file)
	{
		try {
			PsimagLite::IoNg::In ioIn(file);
			ioIn.close();
			return true;
		} catch (...) { }

		return false;
	}
};
} // namespace Dmrg
/*@}*/
#endif
