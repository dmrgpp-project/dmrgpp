#ifndef SVD_H
#define SVD_H
#include "Matrix.h"
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

namespace PsimagLite {

template <typename ComplexOrRealType> class Svd {

public:

	using RealType = typename Real<ComplexOrRealType>::Type;

	Svd(String name = "gesdd")
	    : name_(name)
	{ }

	bool canTryAgain() const { return (name_ == "gesdd"); }

	String name() const { return name_; }

	void operator()(char                             jobz,
	                Matrix<ComplexOrRealType>&       a,
	                typename Vector<RealType>::Type& s,
	                Matrix<ComplexOrRealType>&       vt)
	{
		if (jobz != 'A' && jobz != 'S') {
			String msg("svd: jobz must be either A or S");
			String jobzString = " ";
			jobzString[0]     = jobz;
			throw RuntimeError(msg + ", not " + jobzString + "\n");
		}

		psimag::LAPACK::IntegerForLapackType m   = a.rows();
		psimag::LAPACK::IntegerForLapackType n   = a.cols();
		psimag::LAPACK::IntegerForLapackType lda = m;
		psimag::LAPACK::IntegerForLapackType min = (m < n) ? m : n;

		s.resize(min);
		psimag::LAPACK::IntegerForLapackType ldu  = m;
		psimag::LAPACK::IntegerForLapackType ucol = (jobz == 'A') ? m : min;
		Matrix<ComplexOrRealType>            u(ldu, ucol);
		psimag::LAPACK::IntegerForLapackType ldvt = (jobz == 'A') ? n : min;
		vt.resize(ldvt, n);
		psimag::LAPACK::IntegerForLapackType lrwork
		    = 2.0 * min * std::max(5 * min + 7, 2 * std::max(m, n) + 2 * min + 1);
		typename Vector<typename Real<ComplexOrRealType>::Type>::Type rwork(lrwork, 0.0);

		typename Vector<ComplexOrRealType>::Type           work(100, 0);
		psimag::LAPACK::IntegerForLapackType               info = 0;
		Vector<psimag::LAPACK::IntegerForLapackType>::Type iwork(8 * min, 0);

		// query optimal work
		psimag::LAPACK::IntegerForLapackType lwork = -1;
		mycall(&jobz,
		       &m,
		       &n,
		       &(a(0, 0)),
		       &lda,
		       &(s[0]),
		       &(u(0, 0)),
		       &ldu,
		       &(vt(0, 0)),
		       &ldvt,
		       &(work[0]),
		       &lwork,
		       &(rwork[0]),
		       iwork.data(),
		       &info);
		if (info != 0) {
			String str(__FILE__);
			str += " svd(...) failed at workspace size calculation ";
			str += "with info=" + ttos(info) + "\n";
			throw RuntimeError(str.c_str());
		}

		RealType lworkReal = PsimagLite::real(work[0]);
		lwork
		    = static_cast<psimag::LAPACK::IntegerForLapackType>(lworkReal) + (m + n) * 256;
		work.resize(lwork + 10);

		// real work:
		mycall(&jobz,
		       &m,
		       &n,
		       &(a(0, 0)),
		       &lda,
		       &(s[0]),
		       &(u(0, 0)),
		       &ldu,
		       &(vt(0, 0)),
		       &ldvt,
		       &(work[0]),
		       &lwork,
		       &(rwork[0]),
		       iwork.data(),
		       &info);
		if (info != 0) {
			if (info < 0)
				throw RuntimeError(String(__FILE__) + ": " + ttos(__LINE__)
				                   + " info= " + ttos(info));
			if (info > 0)
				std::cerr << "WARNING " << __FILE__ << ": " << __LINE__
				          << " info= " << info << "\n";
		}

		a = u;
	}

private:

	void mycall(char*                                 jobz,
	            psimag::LAPACK::IntegerForLapackType* m,
	            psimag::LAPACK::IntegerForLapackType* n,
	            ComplexOrRealType*                    a, // T*,
	            psimag::LAPACK::IntegerForLapackType* lda,
	            RealType*                             s,
	            ComplexOrRealType*                    u, // T*,
	            psimag::LAPACK::IntegerForLapackType* ldu,
	            ComplexOrRealType*                    vt, // T*,
	            psimag::LAPACK::IntegerForLapackType* ldvt,
	            ComplexOrRealType*                    work, // T*,
	            psimag::LAPACK::IntegerForLapackType* lwork,
	            RealType*                             rwork, // nothing
	            psimag::LAPACK::IntegerForLapackType* iwork,
	            psimag::LAPACK::IntegerForLapackType* info)
	{
		if (name_ == "gesdd") {
			psimag::LAPACK::GESDD(jobz,
			                      m,
			                      n,
			                      a,
			                      lda,
			                      s,
			                      u,
			                      ldu,
			                      vt,
			                      ldvt,
			                      work,
			                      lwork,
			                      rwork,
			                      iwork,
			                      info);
		} else if (name_ == "gesvd") {
			psimag::LAPACK::GESVD(jobz,
			                      jobz,
			                      m,
			                      n,
			                      a,
			                      lda,
			                      s,
			                      u,
			                      ldu,
			                      vt,
			                      ldvt,
			                      work,
			                      lwork,
			                      rwork,
			                      info);
		} else {
			throw PsimagLite::RuntimeError("Unknown backend " + name_ + "\n");
		}
	}

	Svd(const Svd&);

	Svd& operator=(const Svd&);

	String name_;
};

} // namespace PsimagLite
#endif // SVD_H
