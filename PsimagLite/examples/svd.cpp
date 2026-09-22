#include <PsimagLite/Matrix.h>
#include <PsimagLite/Svd.h>
#include <PsimagLite/Vector.h>

int main()
{
	typedef PsimagLite::Matrix<double> MatrixType;

	MatrixType a(4, 2);
	a(0, 0) = 2;
	a(0, 1) = 4;
	a(1, 0) = 1;
	a(1, 1) = 3;

	std::cout << "A\n";
	std::cout << a;

	MatrixType m(a);

	PsimagLite::Vector<double>::Type s;
	MatrixType                       vt;

	PsimagLite::Svd<double> svd;
	svd('A', a, s, vt);

	std::cout << "U\n";
	std::cout << a;
	std::cout << "S\n";
	std::cout << s;

	std::cout << "\n\nFallback\n";
	PsimagLite::Svd<double> svd_fallback("gesvd");
	svd_fallback('A', m, s, vt);
	std::cout << "U\n";
	std::cout << m;
	std::cout << "S\n";
	std::cout << s;
}
