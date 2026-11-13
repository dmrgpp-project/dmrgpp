#ifndef QUASICANONICAL_H
#define QUASICANONICAL_H
#include "AST/ExpressionForAST.h"
#include "AST/PlusMinusMultiplyDivide.h"
#include "CmplxOrReal.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace PsimagLite {

// level one parens only
// TODO FIXME : Generalize to multiple levels
template <typename ComplexOrRealType> class QuasiCanonical {

public:

	using VectorStringType = Vector<String>::Type;
	using VectorBoolType   = Vector<bool>::Type;
	using VectorType       = typename Vector<ComplexOrRealType>::Type;
	using RealType         = typename Real<ComplexOrRealType>::Type;
	using ComplexType      = std::complex<RealType>;

	QuasiCanonical(String str)
	    : str_(str)
	{
		const SizeType len = str_.length();
		String         tempBuffer;
		bool           parensOpen = false;
		SizeType       bracesOpen = 0;
		char           prev       = '\0';
		for (SizeType i = 0; i < len; ++i) {
			const char c = str_[i];

			if (bracesOpen > 0) {
				mainBuffer_ += c;
				if (c == '{')
					++bracesOpen;
				else if (c == '}')
					--bracesOpen;
				continue;
			}

			if (c == '{') {
				mainBuffer_ += c;
				++bracesOpen;
				continue;
			}

			if (c == '}')
				throw RuntimeError("Unbalanced braces, closed\n");

			if (c == '-' && i > 0) {
				if (prev != '(' && prev != '+') {
					throw RuntimeError("The - sign must be preceded by "
					                   "nothing, parens, or +\n");
				}
			}

			prev = c;

			if (c == '(') {
				if (parensOpen)
					throw RuntimeError("Too many parens levels (one only "
					                   "supported for now)\n");
				parensOpen = true;
				continue;
			}

			if (c == ')') {
				if (!parensOpen)
					throw RuntimeError("Unbalanced parens, closed\n");
				parensOpen = false;
				mainBuffer_ += "@" + ttos(ats_.size()) + "@";
				ats_.push_back(tempBuffer);
				tempBuffer = "";
				continue;
			}

			if (!parensOpen) {
				mainBuffer_ += c;
			} else {
				tempBuffer += c;
			}
		}

		if (bracesOpen > 0)
			throw RuntimeError("Unbalanced braces, open\n");

		if (parensOpen)
			throw RuntimeError("Unbalanced parens, open\n");

		splitTopLevel(terms_, mainBuffer_, '+');

		const SizeType nscalars = ats_.size();
		cachedValues_.resize(nscalars);
		for (SizeType i = 0; i < nscalars; ++i) {
			cachedValues_[i] = simpleArithmetics(ats_[i]);
		}
	}

	SizeType numberOfTerms() const { return terms_.size(); }

	static void splitTopLevel(VectorStringType& terms, const String& str, char separator)
	{
		terms.clear();
		String   term;
		SizeType bracesOpen = 0;
		for (SizeType i = 0; i < str.length(); ++i) {
			const char c = str[i];
			if (c == '{') {
				++bracesOpen;
			} else if (c == '}') {
				if (bracesOpen == 0)
					throw RuntimeError("Unbalanced braces, closed\n");
				--bracesOpen;
			}

			if (c == separator && bracesOpen == 0) {
				terms.push_back(term);
				term = "";
				continue;
			}

			term += c;
		}

		if (bracesOpen > 0)
			throw RuntimeError("Unbalanced braces, open\n");

		terms.push_back(term);
	}

	const String& term(SizeType i) const
	{
		assert(i < terms_.size());
		return terms_[i];
	}

	int scalarIndex(String str) const
	{
		const SizeType len = str.length();
		if (len < 3)
			return -1;
		if (str[0] != '@' || str[len - 1] != '@')
			return -1;

		String   snumber = str.substr(1, len - 2);
		SizeType number  = atoi(snumber);
		if (number >= ats_.size())
			return -1;
		return number;
	}

	ComplexOrRealType scalarFromIndex(SizeType ind) const
	{
		assert(ind < cachedValues_.size());

		return cachedValues_[ind];
	}

	static bool isPureComplex(String t)
	{
		if (t == "i")
			return true;

		const SizeType n = t.length();
		if (n < 2)
			return false;
		String tmp = t.substr(0, n - 1);
		return isRealScalar(tmp);
	}

	static ComplexOrRealType pureComplex(String t)
	{
		static const bool isComplex = IsComplexNumber<ComplexOrRealType>::True;
		if (!isComplex)
			err("i = sqrt(-1) found in code path that is real\n");

		CpmlxOrReal<RealType, (isComplex) ? 1 : 0> cmplxOrReal(t);
		return cmplxOrReal.value();
	}

	static bool isRealScalar(String termStr)
	{
		const SizeType n = termStr.length();
		if (n == 0)
			err("CanonicalExpression: term must not be empty\n");

		for (SizeType i = 0; i < n; ++i) {
			char c       = termStr[i];
			bool isDigit = (c >= '0' && c <= '9');
			if (c == '.' || c == '-' || c == '+' || isDigit)
				continue;
			return false;
		}

		return true;
	}

private:

	ComplexOrRealType simpleArithmetics(String str)
	{
		replaceAll(str, "pi", ttos(M_PI));

		VectorStringType ve;
		split(ve, str, ":");

		using PrimitivesType = PlusMinusMultiplyDivide<ComplexOrRealType>;
		PrimitivesType                   primitives;
		ExpressionForAST<PrimitivesType> expresionForAST(ve, primitives);
		return expresionForAST.exec();
	}

	String           str_;
	String           mainBuffer_;
	VectorStringType ats_;
	VectorStringType terms_;
	VectorType       cachedValues_;
};
} // namespace PsimagLite
#endif // QUASICANONICAL_H
