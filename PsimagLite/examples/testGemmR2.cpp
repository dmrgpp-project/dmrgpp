#include <PsimagLite/BLAS.h>
#include <PsimagLite/GemmR.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Parallelizer2.h>
#include <PsimagLite/Random48.h>

typedef double RealType;

template <typename SomeRngType>
void fillRandom(PsimagLite::Matrix<RealType>& m, RealType min, RealType max, SomeRngType& rng)
{
	const SizeType rows = m.rows();
	const SizeType cols = m.cols();
	for (SizeType i = 0; i < rows; ++i)
		for (SizeType j = 0; j < cols; ++j)
			m(i, j) = min + rng() * max;
}

bool equalMatrices(PsimagLite::Matrix<RealType>& a,
                   PsimagLite::Matrix<RealType>& b,
                   RealType                      tolerance)
{
	const SizeType rows = a.rows();
	const SizeType cols = b.cols();
	if (rows != b.rows() || cols != b.cols())
		return false;
	for (SizeType i = 0; i < rows; ++i)
		for (SizeType j = 0; j < cols; ++j)
			if (fabs(a(i, j) - b(i, j)) > tolerance)
				return false;

	return true;
}

int main(int argc, char** argv)
{
	if (argc < 2)
		throw PsimagLite::RuntimeError("USAGE: " + PsimagLite::String(argv[0])
		                               + " total nthreadsOuter nthreadsInner\n");

	const bool needs_printing = false;
	int const  nb             = 99;
	int        total          = atoi(argv[1]);
	int        nthreads_outer = atoi(argv[2]);
	int        nthreads_inner = atoi(argv[3]);

	PsimagLite::Concurrency concurrency(&argc, &argv, nthreads_inner);

	PsimagLite::Random48<RealType> rng(1234);

	auto lambda = [&rng, nthreads_inner](SizeType, SizeType)
	{
		PsimagLite::GemmR<RealType>  gemm_r(needs_printing, nb, nthreads_inner);
		SizeType                     lda = static_cast<SizeType>(rng() * 500) + 10;
		SizeType                     cda = lda;
		SizeType                     ldb = lda;
		SizeType                     cdb = lda;
		SizeType                     ldc = lda;
		SizeType                     cdc = lda;
		PsimagLite::Matrix<RealType> a(lda, cda);
		PsimagLite::Matrix<RealType> b(ldb, cdb);
		PsimagLite::Matrix<RealType> c(ldc, cdc);

		fillRandom(a, -10, 10, rng);
		fillRandom(b, -10, 10, rng);
		gemm_r(
		    'N', 'N', ldc, cdc, cda, 1.0, &a(0, 0), lda, &b(0, 0), ldb, 0.0, &c(0, 0), ldc);

		PsimagLite::Matrix<RealType> c2(ldc, cdc);
		psimag::BLAS::GEMM('N',
		                   'N',
		                   ldc,
		                   cdc,
		                   cda,
		                   1.0,
		                   &A(0, 0),
		                   lda,
		                   &B(0, 0),
		                   ldb,
		                   0.0,
		                   &C2(0, 0),
		                   ldc);
		if (!equalMatrices(c, c2, 1e-6))
			throw PsimagLite::RuntimeError("TEST FAILED\n");
	};

	PsimagLite::CodeSectionParams csp = PsimagLite::Concurrency::codeSectionParams;
	csp.npthreads                     = nthreads_outer;

	PsimagLite::Parallelizer2<> parallelizer2(csp);
	parallelizer2.parallelFor(0, total, lambda);
}
