#ifndef GSEXPRESSION_HPP
#define GSEXPRESSION_HPP
#include <PsimagLite/PsimagLite.h>

#include <string>
#include <vector>

namespace Dmrg {

class GsExpression {
public:

	static void setPureState(std::vector<SizeType>& indices, const std::string& expr)
	{
		if (expr.empty())
			return;

		const PsimagLite::String expression(expr.begin(), expr.end());
		const PsimagLite::String prefix("PureState[");
		if (expression.length() <= prefix.length()
		    || expression.substr(0, prefix.length()) != prefix
		    || expression[expression.length() - 1] != ']') {
			throw PsimagLite::RuntimeError("Expected PureState[...] for gs\n");
		}

		PsimagLite::Vector<PsimagLite::String>::Type tokens;
		PsimagLite::split(
		    tokens,
		    expression.substr(prefix.length(), expression.length() - prefix.length() - 1),
		    ",");

		std::vector<SizeType> pureIndices;
		pureIndices.reserve(tokens.size());
		for (const auto& token : tokens) {
			if (!PsimagLite::isAnInteger(token))
				throw PsimagLite::RuntimeError(
				    "PureState indices must be non-negative integers\n");
			pureIndices.push_back(PsimagLite::atoi(token));
		}

		indices.swap(pureIndices);
	}
};
}
#endif // GSEXPRESSION_HPP
