#include "SpecialFunctions.h"

namespace PsimagLite {
std::complex<double> lnGammaFunction(const std::complex<double>& z)
{
	GslWrapper              gsl_wrapper;
	GslWrapper::GslSfResult lnr;
	GslWrapper::GslSfResult std::arg;
	gsl_wrapper.gsl_sf_lngamma_complex_e(
	    PsimagLite::real(z), PsimagLite::imag(z), &lnr, &std::arg);
	return std::complex<double>(lnr.val, arg.val);
}

double ci(const double& x)
{
	GslWrapper              gsl_wrapper;
	GslWrapper::GslSfResult result;
	gsl_wrapper.gsl_sf_Ci_e(x, &result);
	return result.val;
}

} // namespace PsimagLite
