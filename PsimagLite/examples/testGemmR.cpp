#include <PsimagLite/GemmR.h>
#include <PsimagLite/Matrix.h>

#include <iostream>

typedef std::complex<double> Zcomplex;

template <typename T> T makeVal(double const x, double const y) { return (x); }

template <> double makeVal<double>(double const x, double const /* y */) { return (x); }

template <> Zcomplex makeVal<Zcomplex>(double const x, double const y)
{
	Zcomplex z(x, y);
	return (z);
}

template <typename T>
int testGemmr(int const Mmax, int const Nmax, int const Kmax, int const nb, bool needsPrinting)
{
	int const idebug  = (needsPrinting) ? 1 : 0;
	int       nerrors = 0;

	char const trans_table[3] = { 'N', 'T', 'C' };

	T const alpha = make_val<T>(1.1, 2.1);
	T const beta  = make_val<T>(3.1, 4.1);

	PsimagLite::GemmR<T> gemm_r(
	    needsPrinting, nb, PsimagLite::Concurrency::codeSectionParams.npthreads);

	for (int k = 1; k <= Kmax; k += nb) {
		for (int n = 1; n <= Nmax; n += nb) {
			for (int m = 1; m <= Mmax; m += nb) {
				for (int itrans_b = 0; itrans_b < 3; itrans_b++) {
					for (int itrans_a = 0; itrans_a < 3; itrans_a++) {

						char const trans_a = trans_table[itrans_a];
						char const trans_b = trans_table[itrans_b];

						bool const is_trans_a
						    = (trans_a == 'T') || (trans_a == 't');
						bool const is_trans_b
						    = (trans_b == 'T') || (trans_b == 't');
						bool const is_conj_a
						    = (trans_a == 'C') || (trans_a == 'c');
						bool const is_conj_b
						    = (trans_b == 'C') || (trans_b == 'c');
						bool const is_notrans_a
						    = (!is_trans_a) && (!is_conj_a);
						bool const is_notrans_b
						    = (!is_trans_b) && (!is_conj_b);

						int const m_c = m;
						int const n_c = n;
						int const m_a = (is_notrans_a) ? m_c : k;
						int const n_a = (is_notrans_a) ? k : m_c;
						int const m_b = (is_notrans_b) ? k : n_c;
						int const n_b = (is_notrans_b) ? n_c : k;

						PsimagLite::Matrix<T> c(m_c, n_c);
						PsimagLite::Matrix<T> c_gemmr(m_c, n_c);
						PsimagLite::Matrix<T> a(m_a, n_a);
						PsimagLite::Matrix<T> b(m_b, n_b);

						int const ld_a = m_a;
						int const ld_b = m_b;
						int const ld_c = m_c;

						for (int j = 0; j < n_c; j++) {
							for (int i = 0; i < m_c; i++) {
								T cij = make_val<T>(
								    1.0 * (i + j) / (m_c + n_c),
								    1.0 * i * j / (m_c * n_c));
								c(i, j)       = cij;
								c_gemmr(i, j) = cij;
							}
						}

						for (int j = 0; j < n_a; j++) {
							for (int i = 0; i < m_a; i++) {
								T aij = make_val<T>(
								    -1.0 * (i + j + 1),
								    1.0 * (j - i + 1));
								a(i, j) = aij;
							}
						}

						for (int j = 0; j < n_b; j++) {
							for (int i = 0; i < m_b; i++) {
								T bij = make_val<T>(
								    1.0 * (i + j + 1) / (m_b * n_b),
								    -1.0 * (j - i + 1)
								        / (m_b * n_b));
								b(i, j) = bij;
							}
						}

						gemm_r(trans_a,
						       trans_b,
						       m,
						       n,
						       k,
						       alpha,
						       &(a(0, 0)),
						       ld_a,
						       &(b(0, 0)),
						       ld_b,
						       beta,
						       &(c_gemmr(0, 0)),
						       ld_c);

						psimag::BLAS::GEMM(transA,
						                   transB,
						                   m,
						                   n,
						                   k,
						                   alpha,
						                   &(A(0, 0)),
						                   ldA,
						                   &(B(0, 0)),
						                   ldB,
						                   beta,
						                   &(C(0, 0)),
						                   ldC);

						double max_err = 0;
						double c_norm  = 0;
						for (int j = 0; j < n_c; j++) {
							for (int i = 0; i < m_c; i++) {
								double const err = std::abs(
								    c(i, j) - c_gemmr(i, j));
								max_err = std::max(max_err, err);
								c_norm += std::abs(c(i, j));
							}
						}

						double const tol  = 0.0000001;
						bool const   isok = (max_err < tol);
						if (!isok) {
							nerrors++;
						}
						if ((!isok) || (idebug >= 1)) {
							std::cout << " transA " << trans_a
							          << " transB " << trans_b << " m "
							          << m << " n " << n << " k " << k
							          << " max_err " << max_err
							          << " c_norm " << c_norm << "\n";
						}
					}
				}
			}
		}
	}

	return (nerrors);
}

int main(int argc, char** argv)
{
	int const nmax          = 300;
	int const mmax          = 301;
	int const kmax          = 302;
	int       nerr_zcomplex = 0;

	if (argc < 2)
		throw PsimagLite::RuntimeError("USAGE: " + PsimagLite::String(argv[0])
		                               + " nthreads [nb] [debug]\n");

	int nthreads = atoi(argv[1]);

	int const  nb             = (argc >= 3) ? atoi(argv[2]) : 99;
	const bool needs_printing = (argc == 4) ? atoi(argv[3]) > 0 : false;

	PsimagLite::Concurrency concurrency(&argc, &argv, nthreads);

	int nerr_double = testGemmr<double>(mmax, nmax, kmax, nb, needs_printing);
	if (nerr_double == 0) {
		nerr_zcomplex = testGemmr<Zcomplex>(mmax, nmax, kmax, nb, needs_printing);
	}

	bool const all_passed = (nerr_double == 0) && (nerr_zcomplex == 0);
	if (all_passed) {
		std::cout << "ALL PASSED "
		          << "\n";
	} else {
		std::cout << " nerr_double = " << nerr_double
		          << " nerr_zcomplex = " << nerr_zcomplex << "\n";
	}
}
