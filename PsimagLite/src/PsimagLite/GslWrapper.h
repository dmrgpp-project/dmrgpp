// SPDX-FileCopyrightText: Copyright (c) 2011-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file GslWrapper.h
 *
 *  Wrapper for GSL functions and types
 */

#ifndef GSL_WRAPPER_H_
#define GSL_WRAPPER_H_

#include "Vector.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <stdexcept>

namespace PsimagLite {

class GslWrapper {
public:

	using gsl_integration_workspace = ::gsl_integration_workspace;
	using gsl_function              = ::gsl_function;
	using gsl_sf_result             = ::gsl_sf_result;

	void printError(int status) const
	{
		std::cerr << "GslWrapper: error: " << gsl_strerror(status) << "\n";
	}

	gsl_error_handler_t* gsl_set_error_handler(gsl_error_handler_t* new_handler) const
	{
		return ::gsl_set_error_handler(new_handler);
	}

	gsl_integration_workspace* gsl_integration_workspace_alloc(SizeType n) const
	{
		return ::gsl_integration_workspace_alloc(n);
	}

	void gsl_integration_workspace_free(gsl_integration_workspace* w) const
	{
		return ::gsl_integration_workspace_free(w);
	}

	int gsl_integration_qagi(gsl_function*              f,
	                         double                     epsabs,
	                         double                     epsrel,
	                         size_t                     limit,
	                         gsl_integration_workspace* workspace,
	                         double*                    result,
	                         double*                    abserr) const
	{
		return ::gsl_integration_qagi(f, epsabs, epsrel, limit, workspace, result, abserr);
	}

	int gsl_integration_qagiu(gsl_function*              f,
	                          double                     a,
	                          double                     epsabs,
	                          double                     epsrel,
	                          size_t                     limit,
	                          gsl_integration_workspace* workspace,
	                          double*                    result,
	                          double*                    abserr) const
	{
		return ::gsl_integration_qagiu(
		    f, a, epsabs, epsrel, limit, workspace, result, abserr);
	}

	int gsl_integration_qagp(const gsl_function*        f,
	                         double*                    pts,
	                         SizeType                   npts,
	                         double                     epsabs,
	                         double                     epsrel,
	                         SizeType                   limit,
	                         gsl_integration_workspace* workspace,
	                         double*                    result,
	                         double*                    abserr) const
	{
		return ::gsl_integration_qagp(
		    f, pts, npts, epsabs, epsrel, limit, workspace, result, abserr);
	}

	int gsl_integration_qag(const gsl_function*        f,
	                        double                     a,
	                        double                     b,
	                        double                     epsabs,
	                        double                     epsrel,
	                        size_t                     limit,
	                        int                        key,
	                        gsl_integration_workspace* workspace,
	                        double*                    result,
	                        double*                    abserr) const
	{
		return ::gsl_integration_qag(
		    f, a, b, epsabs, epsrel, limit, key, workspace, result, abserr);
	}

	int
	gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result* lnr, gsl_sf_result* arg) const
	{
		return ::gsl_sf_lngamma_complex_e(zr, zi, lnr, arg);
	}

	int gsl_sf_Ci_e(double x, gsl_sf_result* result) const { return ::gsl_sf_Ci_e(x, result); }

}; // class GslWrapper
} // namespace PsimagLite

/*@}*/
#endif // GSL_WRAPPER_H_
