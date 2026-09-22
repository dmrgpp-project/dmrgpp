#include <PsimagLite/AST/ExpressionForAST.h>
#include <PsimagLite/AST/PlusMinusMultiplyDivide.h>
#include <PsimagLite/PredicateAwesome.h>
#include <PsimagLite/PsimagLite.h>

#ifdef USE_COMPLEX
typedef std::complex<double> ComplexOrRealType;
#else
typedef double ComplexOrRealType;
#endif

int main(int argc, char** argv)
{
	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	if (argc != 2) {
		std::cerr << "USAGE: " << argv[0] << " expression\n";
		return 1;
	}

	typedef PsimagLite::Vector<PsimagLite::String>::Type           VectorStringType;
	typedef PsimagLite::PlusMinusMultiplyDivide<ComplexOrRealType> PrimitivesType;
	PsimagLite::String                                             str(argv[1]);
	PsimagLite::replaceAll(str, "%t", "0.25");
	VectorStringType ve;
	PsimagLite::split(ve, str, ":");

	PrimitivesType                               primitives;
	PsimagLite::ExpressionForAST<PrimitivesType> expresion_for_ast(ve, primitives);

	std::cout << argv[1] << "\t" << expresion_for_ast.exec() << "\n";
}
