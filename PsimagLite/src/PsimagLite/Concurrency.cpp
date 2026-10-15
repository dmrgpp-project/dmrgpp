#include "Concurrency.h"

namespace PsimagLite {

SizeType                    Concurrency::mode = 0;
LabelDisabled               Concurrency::mpiDisabled_;
CodeSectionParams           Concurrency::codeSectionParams(1);
Map<String, SizeType>::Type Concurrency::codeSectionParamsOverrides_;

const Vector<String>::Type& Concurrency::codeSectionNames()
{
	static const Vector<String>::Type names = { "KronMatrix",
		                                    "ModelBase",
		                                    "Operators",
		                                    "DensityMatrixLocal",
		                                    "DensityMatrixSvd",
		                                    "TwoPointCorrelations",
		                                    "Observer",
		                                    "ObservableLibrary",
		                                    "WftAccelBlocks",
		                                    "BlockDiagWfCtor",
		                                    "BlockDiagWfTransform",
		                                    "WaveFunctionTransfLocal",
		                                    "WftAccelSvdOne",
		                                    "WftAccelSvdTwo" };

	return names;
}

bool Concurrency::isCodeSectionName(const String& name)
{
	const auto& names = codeSectionNames();
	return (std::find(names.begin(), names.end(), name) != names.end());
}

CodeSectionParams Concurrency::codeSectionParamsFor(const String& name)
{
	if (!isCodeSectionName(name))
		throw RuntimeError("Concurrency: unknown code section " + name + "\n");

	CodeSectionParams params = codeSectionParams;
	params.name              = name;
	const auto it            = codeSectionParamsOverrides_.find(name);
	if (it != codeSectionParamsOverrides_.end())
		params.npthreads = it->second;

	return params;
}

void Concurrency::setCodeSectionParamsOverride(const String& name, SizeType threads)
{
	if (!isCodeSectionName(name))
		throw RuntimeError("Concurrency: unknown code section " + name + "\n");

	checkThreadCount(threads, codeSectionParams.npthreadsLevelTwo);
	codeSectionParamsOverrides_[name] = threads;
}

void Concurrency::clearCodeSectionParamsOverrides() { codeSectionParamsOverrides_.clear(); }
} // namespace PsimagLite
