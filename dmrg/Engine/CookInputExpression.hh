//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   src/Engine/CookInputExpression.hh
 * \brief  Process input expressions
 * \note   Copyright (c) 2026 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef COOKINPUTEXPRESSION_HH
#define COOKINPUTEXPRESSION_HH
#include "InputCheck.h"
#include <PsimagLite/AST/ExpressionForAST.h>
#include <PsimagLite/AST/PlusMinusMultiplyDivide.h>
#include <PsimagLite/InputNg.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/PsimagLite.h>
#include <string>
#include <vector>

namespace Dmrg {
//===========================================================================//
/*!
 * \class CookInputExpression
 *
 * \brief Process input expressions including replacements and tables
 *
 * In the future, this class will be used in more places to allow for user
 * more complicated user input.
 *
 * \tparam ComplexOrRealType The type of numbers whether real or complex
 */
//===========================================================================//

template <typename ComplexOrRealType> class CookInputExpression {
public:

	using RealType             = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType          = PsimagLite::InputNg<Dmrg::InputCheck>;
	using InputNgValidatorType = InputNgType::Readable;

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Constructor
	 *
	 * It keeps a reference to the input parameters.
	 * If using Ainur the io could be const throughout. However, for the plain old input
	 * format, the input parameter reads needs also write access due to permitting
	 * duplicated labels. In the future we may disable support for the plain old input format.
	 *
	 * \param[in/out] input_params  The input object; sorry about the non-constness
	 */
	CookInputExpression(InputNgValidatorType& io)
	    : io_(io)
	{ }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Operator parens
	 *
	 * Convert an input string into a complex number, after replacements, and
	 * expression processing. We are using the AST class for expressions here.
	 * For example, the string "+:3:4" would return 7.
	 *
	 * \param[in] factor           The string coming from the input file entry
	 * \param[in] value_for_empty  The value to return if factor is empty
	 * \param[in] time             The time, whose meaning depends on the target being run
	 *
	 * \returns The real or complex value resulting from cooking the string factor
	 */
	ComplexOrRealType operator()(const std::string&       factor,
	                             const ComplexOrRealType& value_for_empty,
	                             const RealType&          time) const
	{
		if (factor.empty())
			return value_for_empty;

		SizeType number_of_sites = 0;
		io_.readline(number_of_sites, "TotalNumberOfSites=");

		std::string str = factor;
		PsimagLite::replaceAll(str, "%t", ttos(time));
		PsimagLite::replaceAll(str, "%n", ttos(number_of_sites));

		std::vector<std::string> ve;
		PsimagLite::split(ve, str, ":");

		for (SizeType i = 0; i < ve.size(); ++i) {
			ve[i] = replaceTables(ve[i]);
		}

		typedef PsimagLite::PlusMinusMultiplyDivide<ComplexOrRealType> PrimitivesType;
		PrimitivesType                                                 primitives;
		PsimagLite::ExpressionForAST<PrimitivesType> expresionForAST(ve, primitives);
		return expresionForAST.exec();
	}

private:

	//---------------------------------------------------------------------------//
	/*!
	 * \brief If present, evaluates a function with an argument
	 *
	 * See input9001.ain in the TestSuite for an example.
	 *
	 * Here other directives could be added in the future.
	 *
	 * \param[in] expr The expression to check for a function or table
	 *
	 * \returns The value of the function at the argument, or the expression unchanged
	 */
	std::string replaceTables(const std::string& expr) const
	{
		std::string label = "!readTable";
		if (expr.substr(0, label.size()) == label) {
			std::string str = expr.substr(label.size(), std::string::npos);

			// delete spaces if any
			str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());

			// delete parens if any
			SizeType length = str.length();
			if (length > 0 && str[0] == '(' && str[length - 1] == ')') {
				str.erase(0, 1);
				if (length > 1) {
					str.erase(length - 2, 1);
				}
			}

			// split on comma
			std::vector<std::string> args;
			PsimagLite::split(args, str, ",");
			if (args.size() != 2) {
				err("readTable expects two arguments\n");
			}

			PsimagLite::Matrix<RealType> matrix;
			io_.read(matrix, args[0]);
			RealType value = findValueFor(matrix, PsimagLite::atof(args[1]));
			return ttos(value);
		} else {
			return expr;
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Evaluate a function (given as a matrix) at a point
	 *
	 * See input9001.ain in the TestSuite for an example.
	 * Here other directives could be added in the future.
	 *
	 * \param[in] matrix The matrix that represents the function
	 * \param[in] t      The argument to the function (first column of the matrix)
	 *
	 * \returns The value of the function at t
	 */
	static RealType findValueFor(const PsimagLite::Matrix<RealType>& matrix, const RealType& t)
	{
		SizeType rows = matrix.rows();
		if (matrix.cols() != 2) {
			err("findValueFor(): not a table\n");
		}

		RealType eps = 1e-6;
		for (SizeType i = 0; i < rows; ++i) {
			if (std::fabs(matrix(i, 0) - t) < eps) {
				return matrix(i, 1);
			}
		}

		throw std::runtime_error("Value for time = " + ttos(t) + " not found in table\n");
	}

	// The input params object
	InputNgValidatorType& io_;
};
}
#endif // COOKINPUTEXPRESSION_HH
