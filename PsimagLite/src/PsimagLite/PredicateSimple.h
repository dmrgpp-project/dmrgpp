#ifndef PREDICATE_SIMPLE_H
#define PREDICATE_SIMPLE_H
#include "AST/ExpressionForAST.h"
#include "AST/PlusMinusMultiplyDivide.h"
#include "PsimagLite.h"
#include "Vector.h"

/* PSIDOC PredicateSimple
 PredicateSimple is of the form
 word operator word
 where operator is in {==, <, >, <=, >=, %%}
 with the usual meaning, and %% means divisible by.
 So l%%2 means that the simple predicate is true if l is divisible by 2.
 */
namespace PsimagLite {

class PredicateSimple {

public:

	using VectorStringType = Vector<String>::Type;

	PredicateSimple(String pred, String separator = ":")
	    : pred_(pred)
	    , separator_(separator)
	{
		SizeType length   = 0;
		size_t   location = String::npos;

		SizeType n = pred.length();
		if (n < 3)
			err("PredicateSimple: pred must have at least 3 "
			    "characters\n");

		for (SizeType i = 0; i < n - 1; ++i) {

			// order matters: test with LONGER FIRST
			String maybe_op = pred.substr(i, 2);
			if (operatorLength(maybe_op) == 2) {
				location = i;
				length   = 2;
				break;
			}

			maybe_op = pred.substr(i, 1);
			if (operatorLength(maybe_op) == 1) {
				location = i;
				length   = 1;
				break;
			}
		}

		if (length == 0)
			err("Could not find operator in predicate " + pred + "\n");

		lhs_ = pred.substr(0, location);
		op_  = pred.substr(location, length);
		rhs_ = pred.substr(location + length, n - location - length);
		if (lhs_ == "" || rhs_ == "")
			err("Left or right expression is empty\n");
	}

	template <typename SomeVectorType>
	bool isTrue(const VectorStringType& names, const SomeVectorType& vals)
	{
		using SomeValueType = typename SomeVectorType::value_type;
		SomeValueType lv    = getValue(lhs_, names, vals);
		SomeValueType rv    = getValue(rhs_, names, vals);
		return compareOnOp(lv, op_, rv);
	}

private:

	template <typename T> static bool compareOnOp(T lv, String op, T rv)
	{
		// {==, <, >, <=, >=, %%}
		// If you add something here, add it also to PredicateSimple.cpp
		if (op == "==")
			return (lv == rv);
		if (op == "!=")
			return (lv != rv);
		if (op == "<")
			return (lv < rv);
		if (op == ">")
			return (lv > rv);
		if (op == "<=")
			return (lv <= rv);
		if (op == ">=")
			return (lv >= rv);
		if (op == "%%")
			return divisibleBy<T>(lv, rv);
		throw RuntimeError("Unknown operator " + op + "\n");
	}

	template <typename T>
	static bool divisibleBy(SizeType lv,
	                        SizeType rv,
	                        typename std::enable_if<Loki::TypeTraits<T>::isArith, int*>::type
	                        = 0)
	{
		return ((lv % rv) == 0);
	}

	static SizeType operatorLength(String op)
	{
		const SizeType n = OPS.size();
		for (SizeType i = 0; i < n; ++i)
			if (op == OPS[i])
				return op.length();

		return 0;
	}

	template <typename SomeVectorType>
	typename SomeVectorType::value_type
	getValue(String hs, const VectorStringType& names, const SomeVectorType& vals)
	{
		String           numeric_hs = replaceVariables(hs, names, vals);
		VectorStringType tokens;
		split(tokens, numeric_hs, separator_);
		using PrimitivesType = PlusMinusMultiplyDivide<typename SomeVectorType::value_type>;
		PrimitivesType                   primitives;
		ExpressionForAST<PrimitivesType> expresion_for_ast(tokens, primitives);
		return expresion_for_ast.exec();
	}

	template <typename SomeVectorType>
	static String
	replaceVariables(String hs, const VectorStringType& names, const SomeVectorType& vals)
	{
		const SizeType n = names.size();
		assert(n == vals.size());
		String buffer = hs;
		for (SizeType i = 0; i < n; ++i) {
			buffer = replaceVariable(buffer, names[i], vals[i]);
		}

		return buffer;
	}

	template <typename T> static String replaceVariable(String hs, String name, T val)
	{
		const String   val_string  = ttos(val);
		const SizeType name_length = name.length();

		while (true) {
			size_t index = hs.find(name);
			if (index == String::npos)
				return hs;
			String part1 = (index == 0) ? "" : hs.substr(0, index);
			String part2
			    = hs.substr(index + name_length, hs.length() - name_length - index);
			hs = part1 + val_string + part2;
		}

		return hs;
	}

	static VectorStringType OPS;
	String                  pred_;
	String                  separator_;
	String                  lhs_;
	String                  op_;
	String                  rhs_;
};
} // namespace PsimagLite
#endif // PREDICATE_SIMPLE_H
