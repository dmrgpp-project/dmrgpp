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
	char                 jobz = option;
	char                 uplo = 'U';
	IntegerForLapackType n    = m.rows();
	IntegerForLapackType lda  = m.cols();
	Vector<double>::Type work(3);
	IntegerForLapackType info;
	IntegerForLapackType lwork = -1;

	if (lda <= 0)
		throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	dsyev_(&jobz, &uplo, &n, &(m(0, 0)), &lda, eigs.data(), work.data(), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}

	const IntegerForLapackType NB = 256;
	lwork = std::max(static_cast<IntegerForLapackType>(1 + static_cast<int>(work[0])),
	                 (NB + 2) * n);
	work.resize(lwork);
	// real work:
	dsyev_(&jobz, &uplo, &n, &(m(0, 0)), &lda, eigs.data(), work.data(), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
}

void diag(Matrix<std::complex<double>>& m, Vector<double>::Type& eigs, char option)
{
	char                               jobz = option;
	char                               uplo = 'U';
	IntegerForLapackType               n    = m.rows();
	IntegerForLapackType               lda  = m.cols();
	Vector<std::complex<double>>::Type work(3);
	Vector<double>::Type               rwork(3 * n);
	IntegerForLapackType               info;
	IntegerForLapackType               lwork = -1;

	eigs.resize(n);

	// query:
	zheev_(&jobz,
	       &uplo,
	       &n,
	       &(m(0, 0)),
	       &lda,
	       eigs.data(),
	       work.data(),
	       &lwork,
	       rwork.data(),
	       &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: zheev_: failed with info!=0.\n");
	}

	const IntegerForLapackType NB = 256;
	lwork
	    = std::max(static_cast<IntegerForLapackType>(1 + static_cast<int>(std::real(work[0]))),
	               (NB + 2) * n);
	work.resize(lwork);
	// real work:
	zheev_(&jobz,
	       &uplo,
	       &n,
	       &(m(0, 0)),
	       &lda,
	       eigs.data(),
	       work.data(),
	       &lwork,
	       rwork.data(),
	       &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: zheev: failed with info!=0.\n");
	}
}

void diag(Matrix<float>& m, Vector<float>::Type& eigs, char option)
{
	char                 jobz = option;
	char                 uplo = 'U';
	IntegerForLapackType n    = m.rows();
	IntegerForLapackType lda  = m.cols();
	Vector<float>::Type  work(3);
	IntegerForLapackType info;
	IntegerForLapackType lwork = -1;

	if (lda <= 0)
		throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	ssyev_(&jobz, &uplo, &n, &(m(0, 0)), &lda, eigs.data(), work.data(), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}

	const IntegerForLapackType NB = 256;
	lwork = std::max(static_cast<IntegerForLapackType>(1 + static_cast<int>(work[0])),
	                 (NB + 2) * n);
	work.resize(lwork);

	// real work:
	ssyev_(&jobz, &uplo, &n, &(m(0, 0)), &lda, eigs.data(), work.data(), &lwork, &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
}

void diag(Matrix<std::complex<float>>& m, Vector<float>::Type& eigs, char option)
{
	char                              jobz = option;
	char                              uplo = 'U';
	IntegerForLapackType              n    = m.rows();
	IntegerForLapackType              lda  = m.cols();
	Vector<std::complex<float>>::Type work(3);
	Vector<float>::Type               rwork(3 * n);
	IntegerForLapackType              info, lwork = -1;

	eigs.resize(n);

	// query:
	cheev_(&jobz,
	       &uplo,
	       &n,
	       &(m(0, 0)),
	       &lda,
	       eigs.data(),
	       work.data(),
	       &lwork,
	       rwork.data(),
	       &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: cheev_: failed with info!=0.\n");
	}

	const IntegerForLapackType NB = 256;
	lwork
	    = std::max(static_cast<IntegerForLapackType>(1 + static_cast<int>(std::real(work[0]))),
	               (NB + 2) * n);
	work.resize(lwork);

	// real work:
	cheev_(&jobz,
	       &uplo,
	       &n,
	       &(m(0, 0)),
	       &lda,
	       eigs.data(),
	       work.data(),
	       &lwork,
	       rwork.data(),
	       &info);
	if (info != 0) {
		std::cerr << "info=" << info << "\n";
		throw RuntimeError("diag: cheev: failed with info!=0.\n");
	}
}

// complex zgeev version
void geev(char                               jobvl,
          char                               jobvr,
          Matrix<std::complex<double>>&      a,
          std::vector<std::complex<double>>& w,
          Matrix<std::complex<double>>&      vl,
          Matrix<std::complex<double>>&      vr)
{
	IntegerForLapackType              n    = a.rows();
	IntegerForLapackType              lda  = a.cols();
	IntegerForLapackType              ldvl = vl.rows();
	IntegerForLapackType              ldvr = vr.rows();
	IntegerForLapackType              info = 0;
	std::vector<std::complex<double>> work(10, 0);
	std::vector<double>               rwork(2 * n + 1, 0);
	IntegerForLapackType              lwork = -1;
	zgeev_(&jobvl,
	       &jobvr,
	       &n,
	       &(a(0, 0)),
	       &lda,
	       w.data(),
	       &(vl(0, 0)),
	       &ldvl,
	       &(vr(0, 0)),
	       &ldvr,
	       work.data(),
	       &lwork,
	       rwork.data(),
	       &info);

	const IntegerForLapackType NB = 256;
	lwork
	    = std::max(static_cast<IntegerForLapackType>(1 + static_cast<int>(std::real(work[0]))),
	               (NB + 2) * n);
	work.resize(lwork);

	zgeev_(&jobvl,
	       &jobvr,
	       &n,
	       &(a(0, 0)),
	       &lda,
	       w.data(),
	       &(vl(0, 0)),
	       &ldvl,
	       &(vr(0, 0)),
	       &ldvr,
	       work.data(),
	       &lwork,
	       rwork.data(),
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
	IntegerForLapackType n    = a.rows();
	IntegerForLapackType lda  = a.cols();
	IntegerForLapackType ldvl = vl.rows();
	IntegerForLapackType ldvr = vr.rows();
	IntegerForLapackType info = 0;
	std::vector<double>  wr(n);
	std::vector<double>  wi(n);
	std::vector<double>  work(10, 0);
	IntegerForLapackType lwork = -1;
	dgeev_(&jobvl,
	       &jobvr,
	       &n,
	       &(a(0, 0)),
	       &lda,
	       wr.data(),
	       wi.data(),
	       &(vl(0, 0)),
	       &ldvl,
	       &(vr(0, 0)),
	       &ldvr,
	       work.data(),
	       &lwork,
	       &info);

	const IntegerForLapackType NB = 256;
	lwork
	    = std::max(static_cast<IntegerForLapackType>(1 + static_cast<int>(std::real(work[0]))),
	               (NB + 2) * n);
	work.resize(lwork);

	dgeev_(&jobvl,
	       &jobvr,
	       &n,
	       &(a(0, 0)),
	       &lda,
	       wr.data(),
	       wi.data(),
	       &(vl(0, 0)),
	       &ldvl,
	       &(vr(0, 0)),
	       &ldvr,
	       work.data(),
	       &lwork,
	       &info);

	checkBlasStatus(info, "dgeev_");

	// fill eigenvalues into complex vector
	for (int i = 0; i < n; ++i) {
		w[i] = std::complex<double>(wr[i], wi[i]);
	}
}

} // namespace PsimagLite
