#include "Matrix.h"

namespace PsimagLite {

void checkBlasStatus(int info, PsimagLite::String msg)
{
	if (info == 0)
		return;

	PsimagLite::String str = msg;
	str += " failed with info = " + ttos(info) + "\n";
	throw RuntimeError(str);
}

void diag(Matrix<double>& m, Vector<double>::Type& eigs, char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: dsyev_: NO LAPACK!\n");
#else
	char                 jobz = option;
	char                 uplo = 'U';
	int                  n    = m.rows();
	int                  lda  = m.cols();
	Vector<double>::Type work(3);
	int                  info;
	int                  lwork = -1;

	if (lda <= 0)
		throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	dsyev(&jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}

	const int nb = 256;
	lwork        = std::max(1 + static_cast<int>(work[0]), (nb + 2) * n);
	work.resize(lwork);
	// real work:
	dsyev(&jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
#endif
}

void diag(Matrix<std::complex<double>>& m, Vector<double>::Type& eigs, char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: zheev: NO LAPACK!\n");
#else
	char                               jobz = option;
	char                               uplo = 'U';
	int                                n    = m.rows();
	int                                lda  = m.cols();
	Vector<std::complex<double>>::Type work(3);
	Vector<double>::Type               rwork(3 * n);
	int                                info;
	int                                lwork = -1;

	eigs.resize(n);

	// query:
	zheev(
	    &jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &(rwork[0]), &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: zheev_: failed with info!=0.\n");
	}

	const int nb = 256;
	lwork        = std::max(1 + static_cast<int>(std::real(work[0])), (nb + 2) * n);
	work.resize(lwork);
	// real work:
	zheev(
	    &jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &(rwork[0]), &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: zheev: failed with info!=0.\n");
	}
#endif
}

void diag(Matrix<float>& m, Vector<float>::Type& eigs, char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: dsyev_: NO LAPACK!\n");
#else
	char                jobz = option;
	char                uplo = 'U';
	int                 n    = m.rows();
	int                 lda  = m.cols();
	Vector<float>::Type work(3);
	int                 info;
	int                 lwork = -1;

	if (lda <= 0)
		throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	ssyev(&jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}

	const int nb = 256;
	lwork        = std::max(1 + static_cast<int>(work[0]), (nb + 2) * n);
	work.resize(lwork);

	// real work:
	ssyev(&jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
#endif
}

void diag(Matrix<std::complex<float>>& m, Vector<float>::Type& eigs, char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: cheev: NO LAPACK!\n");
#else
	char                              jobz = option;
	char                              uplo = 'U';
	int                               n    = m.rows();
	int                               lda  = m.cols();
	Vector<std::complex<float>>::Type work(3);
	Vector<float>::Type               rwork(3 * n);
	int                               info, lwork = -1;

	eigs.resize(n);

	// query:
	cheev(
	    &jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &(rwork[0]), &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: cheev_: failed with info!=0.\n");
	}

	const int nb = 256;
	lwork        = std::max(1 + static_cast<int>(std::real(work[0])), (nb + 2) * n);
	work.resize(lwork);

	// real work:
	cheev(
	    &jobz, &uplo, &n, &(m(0, 0)), &lda, &(eigs[0]), &(work[0]), &lwork, &(rwork[0]), &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: cheev: failed with info!=0.\n");
	}
#endif
}

// complex zgeev version
void geev(char                               jobvl,
          char                               jobvr,
          Matrix<std::complex<double>>&      a,
          std::vector<std::complex<double>>& w,
          Matrix<std::complex<double>>&      vl,
          Matrix<std::complex<double>>&      vr)
{
	int                               n    = a.rows();
	int                               lda  = a.cols();
	int                               ldvl = vl.rows();
	int                               ldvr = vr.rows();
	int                               info = 0;
	std::vector<std::complex<double>> work(10, 0);
	std::vector<double>               rwork(2 * n + 1, 0);
	int                               lwork = -1;
	zgeev(&jobvl,
	      &jobvr,
	      &n,
	      &(a(0, 0)),
	      &lda,
	      &(w[0]),
	      &(vl(0, 0)),
	      &ldvl,
	      &(vr(0, 0)),
	      &ldvr,
	      &(work[0]),
	      &lwork,
	      &(rwork[0]),
	      &info);

	const int nb = 256;
	lwork        = std::max(1 + static_cast<int>(std::real(work[0])), (nb + 2) * n);
	work.resize(lwork);

	zgeev(&jobvl,
	      &jobvr,
	      &n,
	      &(a(0, 0)),
	      &lda,
	      &(w[0]),
	      &(vl(0, 0)),
	      &ldvl,
	      &(vr(0, 0)),
	      &ldvr,
	      &(work[0]),
	      &lwork,
	      &(rwork[0]),
	      &info);

	checkBlasStatus(info, "zgeev_");
}

// real dgeev version
// Note that w will be filled with the complex eigenvalues
// but dgeev splits them in real and imag so we have to postprocess
// which is done at the end of this function
void geev(char                               jobvl,
          char                               jobvr,
          Matrix<double>&                    a,
          std::vector<std::complex<double>>& w,
          Matrix<double>&                    vl,
          Matrix<double>&                    vr)
{
	int                 n    = a.rows();
	int                 lda  = a.cols();
	int                 ldvl = vl.rows();
	int                 ldvr = vr.rows();
	int                 info = 0;
	std::vector<double> wr(n);
	std::vector<double> wi(n);
	std::vector<double> work(10, 0);
	std::vector<double> rwork(2 * n + 1, 0);
	int                 lwork = -1;
	dgeev(&jobvl,
	      &jobvr,
	      &n,
	      &(a(0, 0)),
	      &lda,
	      &(wr[0]),
	      &(wi[0]),
	      &(vl(0, 0)),
	      &ldvl,
	      &(vr(0, 0)),
	      &ldvr,
	      &(work[0]),
	      &lwork,
	      &(rwork[0]),
	      &info);

	const int nb = 256;
	lwork        = std::max(1 + static_cast<int>(std::real(work[0])), (nb + 2) * n);
	work.resize(lwork);

	dgeev(&jobvl,
	      &jobvr,
	      &n,
	      &(a(0, 0)),
	      &lda,
	      &(wr[0]),
	      &(wi[0]),
	      &(vl(0, 0)),
	      &ldvl,
	      &(vr(0, 0)),
	      &ldvr,
	      &(work[0]),
	      &lwork,
	      &(rwork[0]),
	      &info);

	checkBlasStatus(info, "dgeev_");

	// fill eigenvalues into complex vector
	for (int i = 0; i < n; ++i) {
		w[i] = std::complex<double>(wr[i], wi[i]);
	}
}

} // namespace PsimagLite
