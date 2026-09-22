//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for BLAS                                                     *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// * Richard Thigpen, ORNL, June 2003                                         *
// ****************************************************************************

#ifndef PSIMAG_BLAS
#define PSIMAG_BLAS
#include "AllocatorCpu.h"

#include <PsimagLite/KokkosGemm.h>

#include <complex>

/** \file BLAS.h
 *  \author Thomas C. Schulthess and Richard N. Thigpen
 */

/** \brief Namespace encapsulating all PsiMag tools.
 */
namespace psimag {

/** \brief Namespace for psimag wrappers of BLAS functions
 */
namespace BLAS {

#ifndef PSI_BLAS_64
	using IntegerForBlasType = int;
#else
	using IntegerForBlasType = long int;
#endif

	//===============================================================
	// MISSING STUFF (BY G.A.)
	// ==============================================================
	extern "C" double
	ddot(IntegerForBlasType*, double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	// ============================================================================
	// = Level 3 BLAS             GEMM
	// ============================================================================
	extern "C" void sgemm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dgemm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void cgemm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zgemm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	//*****************************************************************************
	//*                           SYMM
	//*****************************************************************************

	extern "C" void ssymm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*);

	extern "C" void dsymm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*);

	extern "C" void csymm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*);
	extern "C" void zsymm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*);

	//*****************************************************************************
	//*                           HEMM
	//*****************************************************************************

	extern "C" void chemm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zhemm(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	// ****************************************************************************
	// *                          SYRK
	// ****************************************************************************

	extern "C" void ssyrk(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dsyrk(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void csyrk(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zsyrk(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	// ****************************************************************************
	// *                          HERK
	// ****************************************************************************
	extern "C" void cherk(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zherk(char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ****************************************************************************
	// *                          SYR2K
	// ****************************************************************************
	extern "C" void ssyr2k(char*,
	                       char*,
	                       IntegerForBlasType*,
	                       IntegerForBlasType*,
	                       const float*,
	                       const float*,
	                       IntegerForBlasType*,
	                       const float*,
	                       IntegerForBlasType*,
	                       const float*,
	                       float*,
	                       IntegerForBlasType*);

	extern "C" void dsyr2k(char*,
	                       char*,
	                       IntegerForBlasType*,
	                       IntegerForBlasType*,
	                       const double*,
	                       const double*,
	                       IntegerForBlasType*,
	                       const double*,
	                       IntegerForBlasType*,
	                       const double*,
	                       double*,
	                       IntegerForBlasType*);

	extern "C" void csyr2k(char*,
	                       char*,
	                       IntegerForBlasType*,
	                       IntegerForBlasType*,
	                       const std::complex<float>*,
	                       const std::complex<float>*,
	                       IntegerForBlasType*,
	                       const std::complex<float>*,
	                       IntegerForBlasType*,
	                       const std::complex<float>*,
	                       std::complex<float>*,
	                       IntegerForBlasType*);

	extern "C" void zsyr2k(char*,
	                       char*,
	                       IntegerForBlasType*,
	                       IntegerForBlasType*,
	                       const std::complex<double>*,
	                       const std::complex<double>*,
	                       IntegerForBlasType*,
	                       const std::complex<double>*,
	                       IntegerForBlasType*,
	                       const std::complex<double>*,
	                       std::complex<double>*,
	                       IntegerForBlasType*);
	// ****************************************************************************
	// *                          HER2k
	// ****************************************************************************
	extern "C" void cher2k(char*,
	                       char*,
	                       IntegerForBlasType*,
	                       IntegerForBlasType*,
	                       const std::complex<float>*,
	                       const std::complex<float>*,
	                       IntegerForBlasType*,
	                       const std::complex<float>*,
	                       IntegerForBlasType*,
	                       const std::complex<float>*,
	                       std::complex<float>*,
	                       IntegerForBlasType*);

	extern "C" void zher2k(char*,
	                       char*,
	                       IntegerForBlasType*,
	                       IntegerForBlasType*,
	                       const std::complex<double>*,
	                       const std::complex<double>*,
	                       IntegerForBlasType*,
	                       const std::complex<double>*,
	                       IntegerForBlasType*,
	                       const std::complex<double>*,
	                       std::complex<double>*,
	                       IntegerForBlasType*);
	// ****************************************************************************
	// *                          TRMM
	// ****************************************************************************
	extern "C" void strmm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dtrmm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctrmm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztrmm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ****************************************************************************
	// *                          TRSM
	// ****************************************************************************
	extern "C" void strsm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dtrsm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctrsm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztrsm(char*,
	                      char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ****************************************************************************
	// *    Level 2 BLAS          GEMV
	// ****************************************************************************
	extern "C" void sgemv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dgemv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void cgemv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zgemv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ****************************************************************************
	// *                          GBMV
	// ****************************************************************************
	extern "C" void sgbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dgbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void cgbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zgbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	// ****************************************************************************
	// *                          HEMV
	// ****************************************************************************
	extern "C" void chemv(char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zhemv(char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ****************************************************************************
	// *                         HBMV
	// ****************************************************************************
	extern "C" void chbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zhbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ****************************************************************************
	// *                         HPMV
	// ****************************************************************************
	extern "C" void chpmv(char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zhpmv(char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         SYMV
	// ******************************************************************************
	extern "C" void ssymv(char*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dsymv(char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         SBMV
	// ******************************************************************************
	extern "C" void ssbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dsbmv(char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         SPMV
	// ******************************************************************************
	extern "C" void sspmv(char*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dspmv(char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         TRMV
	// ******************************************************************************
	extern "C" void strmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dtrmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctrmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztrmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	// ******************************************************************************
	// *                         TBMV
	// ******************************************************************************
	extern "C" void stbmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dtbmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctbmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztbmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         TPMV
	// ******************************************************************************
	extern "C" void
	stpmv(char*, char*, char*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dtpmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctpmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztpmv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         TRSV
	// ******************************************************************************
	extern "C" void strsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dtrsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctrsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztrsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         TBSV
	// ******************************************************************************
	extern "C" void stbsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dtbsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctbsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztbsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         TPSV
	// ******************************************************************************
	extern "C" void
	stpsv(char*, char*, char*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dtpsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ctpsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void ztpsv(char*,
	                      char*,
	                      char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         GER
	// ******************************************************************************
	extern "C" void sger(IntegerForBlasType*,
	                     IntegerForBlasType*,
	                     const float*,
	                     const float*,
	                     IntegerForBlasType*,
	                     const float*,
	                     IntegerForBlasType*,
	                     float*,
	                     IntegerForBlasType*);

	extern "C" void dger(IntegerForBlasType*,
	                     IntegerForBlasType*,
	                     const double*,
	                     const double*,
	                     IntegerForBlasType*,
	                     const double*,
	                     IntegerForBlasType*,
	                     double*,
	                     IntegerForBlasType*);
	// ******************************************************************************
	// *                         GERU
	// ******************************************************************************
	extern "C" void cgeru(IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zgeru(IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         GERC
	// ******************************************************************************
	extern "C" void cgerc(IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zgerc(IntegerForBlasType*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         HER
	// ******************************************************************************
	extern "C" void cher(char*,
	                     IntegerForBlasType*,
	                     const std::complex<float>*,
	                     const std::complex<float>*,
	                     IntegerForBlasType*,
	                     std::complex<float>*,
	                     IntegerForBlasType*);

	extern "C" void zher(char*,
	                     IntegerForBlasType*,
	                     const std::complex<double>*,
	                     const std::complex<double>*,
	                     IntegerForBlasType*,
	                     std::complex<double>*,
	                     IntegerForBlasType*);
	// ******************************************************************************
	// *                         HPR
	// ******************************************************************************
	extern "C" void chpr(char*,
	                     IntegerForBlasType*,
	                     const std::complex<float>*,
	                     const std::complex<float>*,
	                     IntegerForBlasType*,
	                     std::complex<float>*);

	extern "C" void zhpr(char*,
	                     IntegerForBlasType*,
	                     const std::complex<double>*,
	                     const std::complex<double>*,
	                     IntegerForBlasType*,
	                     std::complex<float>*);
	// ******************************************************************************
	// *                         HER2
	// ******************************************************************************
	extern "C" void cher2(char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zher2(char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         HPR2
	// ******************************************************************************
	extern "C" void chpr2(char*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*);

	extern "C" void zhpr2(char*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*);
	// ******************************************************************************
	// *                         SYR
	// ******************************************************************************
	extern "C" void ssyr(char*,
	                     IntegerForBlasType*,
	                     const float*,
	                     const float*,
	                     IntegerForBlasType*,
	                     float*,
	                     IntegerForBlasType*);

	extern "C" void dsyr(char*,
	                     IntegerForBlasType*,
	                     const double*,
	                     const double*,
	                     IntegerForBlasType*,
	                     double*,
	                     IntegerForBlasType*);
	// ******************************************************************************
	// *                         SPR
	// ******************************************************************************
	extern "C" void
	sspr(char*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, float*);

	extern "C" void dspr(char*,
	                     IntegerForBlasType*,
	                     const double*,
	                     const double*,
	                     IntegerForBlasType*,
	                     double*);
	// ******************************************************************************
	// *                         SYR2
	// ******************************************************************************
	extern "C" void ssyr2(char*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void dsyr2(char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);
	// ******************************************************************************
	// *                         SPR2
	// ******************************************************************************
	extern "C" void sspr2(char*,
	                      IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*);

	extern "C" void dspr2(char*,
	                      IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*);
	// ******************************************************************************
	// *Level 1 BLAS
	// ******************************************************************************

	extern "C" void srotg(float*, float*, float*, float*);
	extern "C" void drotg(double*, double*, double*, double*);

	extern "C" void srotmg(float*, float*, float*, float*, float*);
	extern "C" void drotmg(double*, double*, double*, double*, double*);

	extern "C" void srot(IntegerForBlasType*,
	                     float*,
	                     IntegerForBlasType*,
	                     float*,
	                     IntegerForBlasType*,
	                     const float*,
	                     const float*);
	extern "C" void drot(IntegerForBlasType*,
	                     double*,
	                     IntegerForBlasType*,
	                     double*,
	                     IntegerForBlasType*,
	                     const double*,
	                     const double*);

	extern "C" void srotm(IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*,
	                      const float*);
	extern "C" void drotm(IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*,
	                      const double*);

	extern "C" void
	sswap(IntegerForBlasType*, float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void
	dswap(IntegerForBlasType*, double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void cswap(IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zswap(IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	extern "C" void saxpy(IntegerForBlasType*,
	                      const float*,
	                      const float*,
	                      IntegerForBlasType*,
	                      float*,
	                      IntegerForBlasType*);

	extern "C" void daxpy(IntegerForBlasType*,
	                      const double*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void caxpy(IntegerForBlasType*,
	                      const std::complex<float>*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zaxpy(IntegerForBlasType*,
	                      const std::complex<double>*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	extern "C" void
	scopy(IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dcopy(IntegerForBlasType*,
	                      const double*,
	                      IntegerForBlasType*,
	                      double*,
	                      IntegerForBlasType*);

	extern "C" void ccopy(IntegerForBlasType*,
	                      const std::complex<float>*,
	                      IntegerForBlasType*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zcopy(IntegerForBlasType*,
	                      const std::complex<double>*,
	                      IntegerForBlasType*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	extern "C" void sscal(IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dscal(IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void cscal(IntegerForBlasType*,
	                      const std::complex<float>*,
	                      std::complex<float>*,
	                      IntegerForBlasType*);

	extern "C" void zscal(IntegerForBlasType*,
	                      const std::complex<double>*,
	                      std::complex<double>*,
	                      IntegerForBlasType*);

	// ============================================================================
	inline double dot(IntegerForBlasType n,
	                  double*            dx,
	                  IntegerForBlasType incx,
	                  double*            dy,
	                  IntegerForBlasType incy)
	{
		return ddot(&n, dx, &incx, dy, &incy);
	}

	// ============================================================================
	inline void gemm(char               c1,
	                 char               c2,
	                 IntegerForBlasType sX,
	                 IntegerForBlasType sY,
	                 IntegerForBlasType sZ,
	                 const float&       a,
	                 const float*       x,
	                 IntegerForBlasType sx,
	                 const float*       y,
	                 IntegerForBlasType sy,
	                 const float&       b,
	                 float*             z,
	                 IntegerForBlasType sz)
	{
		PsimagLite::kokkosGemm(c1, c2, sX, sY, sZ, a, x, sx, y, sy, b, z, sz);
	}

	inline void gemm(char               c1,
	                 char               c2,
	                 IntegerForBlasType sX,
	                 IntegerForBlasType sY,
	                 IntegerForBlasType sZ,
	                 const double&      a,
	                 const double*      x,
	                 IntegerForBlasType sx,
	                 const double*      y,
	                 IntegerForBlasType sy,
	                 const double&      b,
	                 double*            z,
	                 IntegerForBlasType sz)
	{
		PsimagLite::kokkosGemm(c1, c2, sX, sY, sZ, a, x, sx, y, sy, b, z, sz);
	}

	inline void gemm(char                       c1,
	                 char                       c2,
	                 IntegerForBlasType         sX,
	                 IntegerForBlasType         sY,
	                 IntegerForBlasType         sZ,
	                 const std::complex<float>& a,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         sx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType         sy,
	                 const std::complex<float>& b,
	                 std::complex<float>*       z,
	                 IntegerForBlasType         sz)
	{
		PsimagLite::kokkosGemm(c1, c2, sX, sY, sZ, a, x, sx, y, sy, b, z, sz);
	}

	inline void gemm(char                        c1,
	                 char                        c2,
	                 IntegerForBlasType          sX,
	                 IntegerForBlasType          sY,
	                 IntegerForBlasType          sZ,
	                 const std::complex<double>& a,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          sx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType          sy,
	                 const std::complex<double>& b,
	                 std::complex<double>*       z,
	                 IntegerForBlasType          sz)
	{
		/* When  TRANSA = 'N' or 'n' then
		        LDA must be at least  max( 1, m ), otherwise  LDA must be at
		least  max( 1, k ).*/

		if (c1 == 'N' || c1 == 'n') {
			if (sx < std::max(1, sX)) {
				throw PsimagLite::RuntimeError("GEMM lda < max(1, m)\n");
			}
		} else {
			if (sx < std::max(1, sZ)) {
				throw PsimagLite::RuntimeError("GEMM lda < max(1, k)\n");
			}
		}
		PsimagLite::kokkosGemm(c1, c2, sX, sY, sZ, a, x, sx, y, sy, b, z, sz);
	}

	// ***************************************************************************
	inline void symm(char               c1,
	                 char               c2,
	                 IntegerForBlasType sX,
	                 IntegerForBlasType sY,
	                 const float&       a,
	                 const float*       x,
	                 IntegerForBlasType sx,
	                 const float*       y,
	                 IntegerForBlasType sy,
	                 const float&       b,
	                 float*             z,
	                 IntegerForBlasType sz)
	{
		ssymm(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}

	inline void symm(char               c1,
	                 char               c2,
	                 IntegerForBlasType sX,
	                 IntegerForBlasType sY,
	                 const double&      a,
	                 const double*      x,
	                 IntegerForBlasType sx,
	                 const double*      y,
	                 IntegerForBlasType sy,
	                 const double&      b,
	                 double*            z,
	                 IntegerForBlasType sz)
	{
		dsymm(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}
	inline void symm(char                       c1,
	                 char                       c2,
	                 IntegerForBlasType         sX,
	                 IntegerForBlasType         sY,
	                 const std::complex<float>& a,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         sx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType         sy,
	                 const std::complex<float>& b,
	                 std::complex<float>*       z,
	                 IntegerForBlasType         sz)
	{
		csymm(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}
	inline void symm(char                        c1,
	                 char                        c2,
	                 IntegerForBlasType          sX,
	                 IntegerForBlasType          sY,
	                 const std::complex<double>& a,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          sx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType          sy,
	                 const std::complex<double>& b,
	                 std::complex<double>*       z,
	                 IntegerForBlasType          sz)
	{
		zsymm(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}
	// ---------------------------------------------------------------------------
	inline void hemm(char                       c1,
	                 char                       c2,
	                 IntegerForBlasType         sX,
	                 IntegerForBlasType         sY,
	                 const std::complex<float>& a,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         sx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType,
	                 const std::complex<float>& b,
	                 std::complex<float>*       z,
	                 IntegerForBlasType         sz)
	{
		chemm(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sx, &b, z, &sz);
	}
	inline void hemm(char                        c1,
	                 char                        c2,
	                 IntegerForBlasType          sX,
	                 IntegerForBlasType          sY,
	                 const std::complex<double>& a,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          sx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType,
	                 const std::complex<double>& b,
	                 std::complex<double>*       z,
	                 IntegerForBlasType          sz)
	{
		zhemm(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sx, &b, z, &sz);
	}
	// **************************************************************************
	inline void syrk(char               UPLO,
	                 char               TRANS,
	                 IntegerForBlasType N,
	                 IntegerForBlasType K,
	                 const float&       ALPHA,
	                 const float*       A,
	                 IntegerForBlasType LDA,
	                 const float&       BETA,
	                 float*             C,
	                 IntegerForBlasType LDC)
	{
		ssyrk(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void syrk(char               UPLO,
	                 char               TRANS,
	                 IntegerForBlasType N,
	                 IntegerForBlasType K,
	                 const double&      ALPHA,
	                 const double*      A,
	                 IntegerForBlasType LDA,
	                 const double&      BETA,
	                 double*            C,
	                 IntegerForBlasType LDC)
	{
		dsyrk(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void syrk(char                       UPLO,
	                 char                       TRANS,
	                 IntegerForBlasType         N,
	                 IntegerForBlasType         K,
	                 const std::complex<float>& ALPHA,
	                 const std::complex<float>* A,
	                 IntegerForBlasType         LDA,
	                 const std::complex<float>& BETA,
	                 std::complex<float>*       C,
	                 IntegerForBlasType         LDC)
	{
		csyrk(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void syrk(char                        UPLO,
	                 char                        TRANS,
	                 IntegerForBlasType          N,
	                 IntegerForBlasType          K,
	                 const std::complex<double>& ALPHA,
	                 const std::complex<double>* A,
	                 IntegerForBlasType          LDA,
	                 const std::complex<double>& BETA,
	                 std::complex<double>*       C,
	                 IntegerForBlasType          LDC)
	{
		zsyrk(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}

	// ***************************************************************************
	inline void herk(char                       UPLO,
	                 char                       TRANS,
	                 IntegerForBlasType         N,
	                 IntegerForBlasType         K,
	                 const std::complex<float>& ALPHA,
	                 const std::complex<float>* A,
	                 IntegerForBlasType         LDA,
	                 const std::complex<float>& BETA,
	                 std::complex<float>*       C,
	                 IntegerForBlasType         LDC)
	{
		cherk(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void herk(char                        UPLO,
	                 char                        TRANS,
	                 IntegerForBlasType          N,
	                 IntegerForBlasType          K,
	                 const std::complex<double>& ALPHA,
	                 const std::complex<double>* A,
	                 IntegerForBlasType          LDA,
	                 const std::complex<double>& BETA,
	                 std::complex<double>*       C,
	                 IntegerForBlasType          LDC)
	{
		zherk(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	// ***************************************************************************
	inline void syR2K(char               uplo,
	                  char               trans,
	                  IntegerForBlasType n,
	                  IntegerForBlasType k,
	                  const float&       alpha,
	                  const float*       A,
	                  IntegerForBlasType lda,
	                  const float*       B,
	                  IntegerForBlasType ldb,
	                  const float&       beta,
	                  float*             C,
	                  IntegerForBlasType ldc)
	{
		ssyr2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void syR2K(char               uplo,
	                  char               trans,
	                  IntegerForBlasType n,
	                  IntegerForBlasType k,
	                  const double&      alpha,
	                  const double*      A,
	                  IntegerForBlasType lda,
	                  const double*      B,
	                  IntegerForBlasType ldb,
	                  const double&      beta,
	                  double*            C,
	                  IntegerForBlasType ldc)
	{
		dsyr2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void syR2k(char                       uplo,
	                  char                       trans,
	                  IntegerForBlasType         n,
	                  IntegerForBlasType         k,
	                  const std::complex<float>& alpha,
	                  const std::complex<float>* A,
	                  IntegerForBlasType         lda,
	                  const std::complex<float>* B,
	                  IntegerForBlasType         ldb,
	                  const std::complex<float>& beta,
	                  std::complex<float>*       C,
	                  IntegerForBlasType         ldc)
	{
		csyr2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void syR2k(char                        uplo,
	                  char                        trans,
	                  IntegerForBlasType          n,
	                  IntegerForBlasType          k,
	                  const std::complex<double>& alpha,
	                  const std::complex<double>* A,
	                  IntegerForBlasType          lda,
	                  const std::complex<double>* B,
	                  IntegerForBlasType          ldb,
	                  const std::complex<double>& beta,
	                  std::complex<double>*       C,
	                  IntegerForBlasType          ldc)
	{
		zsyr2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	// ***************************************************************************
	inline void heR2k(char                       uplo,
	                  char                       trans,
	                  IntegerForBlasType         n,
	                  IntegerForBlasType         k,
	                  const std::complex<float>& alpha,
	                  const std::complex<float>* A,
	                  IntegerForBlasType         lda,
	                  const std::complex<float>* B,
	                  IntegerForBlasType         ldb,
	                  const std::complex<float>& beta,
	                  std::complex<float>*       C,
	                  IntegerForBlasType         ldc)
	{
		cher2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void heR2k(char                        uplo,
	                  char                        trans,
	                  IntegerForBlasType          n,
	                  IntegerForBlasType          k,
	                  const std::complex<double>& alpha,
	                  const std::complex<double>* A,
	                  IntegerForBlasType          lda,
	                  const std::complex<double>* B,
	                  IntegerForBlasType          ldb,
	                  const std::complex<double>& beta,
	                  std::complex<double>*       C,
	                  IntegerForBlasType          ldc)
	{
		zher2k(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	// ********************************************************************************
	inline void trmm(char               side,
	                 char               uplo,
	                 char               transa,
	                 char               diag,
	                 IntegerForBlasType m,
	                 IntegerForBlasType n,
	                 const float&       alpha,
	                 const float*       A,
	                 IntegerForBlasType lda,
	                 float*             B,
	                 IntegerForBlasType ldb)
	{
		strmm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void trmm(char               side,
	                 char               uplo,
	                 char               transa,
	                 char               diag,
	                 IntegerForBlasType m,
	                 IntegerForBlasType n,
	                 const double&      alpha,
	                 const double*      A,
	                 IntegerForBlasType lda,
	                 double*            B,
	                 IntegerForBlasType ldb)
	{
		dtrmm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void trmm(char                       side,
	                 char                       uplo,
	                 char                       transa,
	                 char                       diag,
	                 IntegerForBlasType         m,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* A,
	                 IntegerForBlasType         lda,
	                 std::complex<float>*       B,
	                 IntegerForBlasType         ldb)
	{
		ctrmm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void trmm(char                        side,
	                 char                        uplo,
	                 char                        transa,
	                 char                        diag,
	                 IntegerForBlasType          m,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* A,
	                 IntegerForBlasType          lda,
	                 std::complex<double>*       B,
	                 IntegerForBlasType          ldb)
	{
		ztrmm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	// ********************************************************************************
	inline void trsm(char               side,
	                 char               uplo,
	                 char               transa,
	                 char               diag,
	                 IntegerForBlasType m,
	                 IntegerForBlasType n,
	                 const float&       alpha,
	                 const float*       A,
	                 IntegerForBlasType lda,
	                 float*             B,
	                 IntegerForBlasType ldb)
	{
		strsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void trsm(char               side,
	                 char               uplo,
	                 char               transa,
	                 char               diag,
	                 IntegerForBlasType m,
	                 IntegerForBlasType n,
	                 const double&      alpha,
	                 const double*      A,
	                 IntegerForBlasType lda,
	                 double*            B,
	                 IntegerForBlasType ldb)
	{
		dtrsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void trsm(char                       side,
	                 char                       uplo,
	                 char                       transa,
	                 char                       diag,
	                 IntegerForBlasType         m,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* A,
	                 IntegerForBlasType         lda,
	                 std::complex<float>*       B,
	                 IntegerForBlasType         ldb)
	{
		ctrsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void trsm(char                        side,
	                 char                        uplo,
	                 char                        transa,
	                 char                        diag,
	                 IntegerForBlasType          m,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* A,
	                 IntegerForBlasType          lda,
	                 std::complex<double>*       B,
	                 IntegerForBlasType          ldb)
	{
		ztrsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	// ***************************************************************************

	inline void gemv(char               c,
	                 IntegerForBlasType M,
	                 IntegerForBlasType N,
	                 const float&       alpha,
	                 const float*       A,
	                 IntegerForBlasType ldA,
	                 const float*       x,
	                 IntegerForBlasType incX,
	                 const float&       beta,
	                 float*             y,
	                 IntegerForBlasType incY)
	{
		sgemv(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ----------------------------------------------------------------------------
	inline void gemv(char               c,
	                 IntegerForBlasType M,
	                 IntegerForBlasType N,
	                 const double&      alpha,
	                 const double*      A,
	                 IntegerForBlasType ldA,
	                 const double*      x,
	                 IntegerForBlasType incX,
	                 const double&      beta,
	                 double*            y,
	                 IntegerForBlasType incY)
	{
		dgemv(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ---------------------------------------------------------------------------
	inline void gemv(char                       c,
	                 IntegerForBlasType         M,
	                 IntegerForBlasType         N,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* A,
	                 IntegerForBlasType         ldA,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incX,
	                 const std::complex<float>& beta,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         incY)
	{
		cgemv(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ---------------------------------------------------------------------------
	inline void gemv(char                        c,
	                 IntegerForBlasType          M,
	                 IntegerForBlasType          N,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* A,
	                 IntegerForBlasType          ldA,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incX,
	                 const std::complex<double>& beta,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          incY)
	{
		zgemv(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ----------------------------------------------------------------------------
	inline void gbmv(char               trans,
	                 IntegerForBlasType m,
	                 IntegerForBlasType n,
	                 IntegerForBlasType kl,
	                 IntegerForBlasType ku,
	                 const float&       alpha,
	                 const float*       A,
	                 IntegerForBlasType lda,
	                 const float*       x,
	                 IntegerForBlasType incx,
	                 const float&       beta,
	                 float*             y,
	                 IntegerForBlasType incy)
	{
		sgbmv(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	inline void gbmv(char               trans,
	                 IntegerForBlasType m,
	                 IntegerForBlasType n,
	                 IntegerForBlasType kl,
	                 IntegerForBlasType ku,
	                 const double&      alpha,
	                 const double*      A,
	                 IntegerForBlasType lda,
	                 const double*      x,
	                 IntegerForBlasType incx,
	                 const double&      beta,
	                 double*            y,
	                 IntegerForBlasType incy)
	{
		dgbmv(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	inline void gbmv(char                       trans,
	                 IntegerForBlasType         m,
	                 IntegerForBlasType         n,
	                 IntegerForBlasType         kl,
	                 IntegerForBlasType         ku,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* A,
	                 IntegerForBlasType         lda,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>& beta,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         incy)
	{
		cgbmv(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	inline void gbmv(char                        trans,
	                 IntegerForBlasType          m,
	                 IntegerForBlasType          n,
	                 IntegerForBlasType          kl,
	                 IntegerForBlasType          ku,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* A,
	                 IntegerForBlasType          lda,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>& beta,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          incy)
	{
		zgbmv(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void hemv(char                       uplo,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* a,
	                 IntegerForBlasType         lda,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>& beta,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         incy)
	{
		chemv(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void hemv(char                        uplo,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* a,
	                 IntegerForBlasType          lda,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>& beta,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          incy)
	{
		zhemv(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// **************************************************************************
	inline void hbmv(char                       uplo,
	                 IntegerForBlasType         n,
	                 IntegerForBlasType         k,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* a,
	                 IntegerForBlasType         lda,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>& beta,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         incy)
	{
		chbmv(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void hbmv(char                        uplo,
	                 IntegerForBlasType          n,
	                 IntegerForBlasType          k,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* a,
	                 IntegerForBlasType          lda,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>& beta,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          incy)
	{
		zhbmv(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// ***************************************************************************
	inline void hpmv(char                       uplo,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* ap,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>& beta,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         incy)
	{
		chpmv(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	inline void hpmv(char                        uplo,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* ap,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>& beta,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          incy)
	{
		zhpmv(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	// ***************************************************************************
	inline void symv(char               uplo,
	                 IntegerForBlasType n,
	                 const float&       alpha,
	                 const float*       a,
	                 IntegerForBlasType lda,
	                 const float*       x,
	                 IntegerForBlasType incx,
	                 const float&       beta,
	                 float*             y,
	                 IntegerForBlasType incy)
	{
		ssymv(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void symv(char               uplo,
	                 IntegerForBlasType n,
	                 const double&      alpha,
	                 const double*      a,
	                 IntegerForBlasType lda,
	                 const double*      x,
	                 IntegerForBlasType incx,
	                 const double&      beta,
	                 double*            y,
	                 IntegerForBlasType incy)
	{
		dsymv(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void sbmv(char               uplo,
	                 IntegerForBlasType n,
	                 IntegerForBlasType k,
	                 const float&       alpha,
	                 const float*       a,
	                 IntegerForBlasType lda,
	                 const float*       x,
	                 IntegerForBlasType incx,
	                 const float&       beta,
	                 float*             y,
	                 IntegerForBlasType incy)
	{
		ssbmv(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void sbmv(char               uplo,
	                 IntegerForBlasType n,
	                 IntegerForBlasType k,
	                 const double&      alpha,
	                 const double*      a,
	                 IntegerForBlasType lda,
	                 const double*      x,
	                 IntegerForBlasType incx,
	                 const double&      beta,
	                 double*            y,
	                 IntegerForBlasType incy)
	{
		dsbmv(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void spmv(char               uplo,
	                 IntegerForBlasType n,
	                 const float&       alpha,
	                 const float*       ap,
	                 const float*       x,
	                 IntegerForBlasType incx,
	                 const float&       beta,
	                 float*             y,
	                 IntegerForBlasType incy)
	{
		sspmv(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	inline void spmv(char               uplo,
	                 IntegerForBlasType n,
	                 const double&      alpha,
	                 const double*      ap,
	                 const double*      x,
	                 IntegerForBlasType incx,
	                 const double&      beta,
	                 double*            y,
	                 IntegerForBlasType incy)
	{
		dspmv(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void trmv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const float*       a,
	                 IntegerForBlasType lda,
	                 float*             x,
	                 IntegerForBlasType incx)
	{
		strmv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void trmv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const double*      a,
	                 IntegerForBlasType lda,
	                 double*            x,
	                 IntegerForBlasType incx)
	{
		dtrmv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void trmv(char                       uplo,
	                 char                       trans,
	                 char                       diag,
	                 IntegerForBlasType         n,
	                 const std::complex<float>* a,
	                 IntegerForBlasType         lda,
	                 std::complex<float>*       x,
	                 IntegerForBlasType         incx)
	{
		ctrmv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void trmv(char                        uplo,
	                 char                        trans,
	                 char                        diag,
	                 IntegerForBlasType          n,
	                 const std::complex<double>* a,
	                 IntegerForBlasType          lda,
	                 std::complex<double>*       x,
	                 IntegerForBlasType          incx)
	{
		ztrmv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void tbmv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 IntegerForBlasType k,
	                 const float*       a,
	                 IntegerForBlasType lda,
	                 float*             x,
	                 IntegerForBlasType incx)
	{
		stbmv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void tbmv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 IntegerForBlasType k,
	                 const double*      a,
	                 IntegerForBlasType lda,
	                 double*            x,
	                 IntegerForBlasType incx)
	{
		dtbmv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void tbmv(char                       uplo,
	                 char                       trans,
	                 char                       diag,
	                 IntegerForBlasType         n,
	                 IntegerForBlasType         k,
	                 const std::complex<float>* a,
	                 IntegerForBlasType         lda,
	                 std::complex<float>*       x,
	                 IntegerForBlasType         incx)
	{
		ctbmv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void tbmv(char                        uplo,
	                 char                        trans,
	                 char                        diag,
	                 IntegerForBlasType          n,
	                 IntegerForBlasType          k,
	                 const std::complex<double>* a,
	                 IntegerForBlasType          lda,
	                 std::complex<double>*       x,
	                 IntegerForBlasType          incx)
	{
		ztbmv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void tpmv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const float*       ap,
	                 float*             x,
	                 IntegerForBlasType incx)
	{
		stpmv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void tpmv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const double*      ap,
	                 double*            x,
	                 IntegerForBlasType incx)
	{
		dtpmv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void tpmv(char                       uplo,
	                 char                       trans,
	                 char                       diag,
	                 IntegerForBlasType         n,
	                 const std::complex<float>* ap,
	                 std::complex<float>*       x,
	                 IntegerForBlasType         incx)
	{
		ctpmv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void tpmv(char                        uplo,
	                 char                        trans,
	                 char                        diag,
	                 IntegerForBlasType          n,
	                 const std::complex<double>* ap,
	                 std::complex<double>*       x,
	                 IntegerForBlasType          incx)
	{
		ztpmv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	// ****************************************************************************
	inline void trsv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const float*       a,
	                 IntegerForBlasType lda,
	                 float*             x,
	                 IntegerForBlasType incx)
	{
		strsv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void trsv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const double*      a,
	                 IntegerForBlasType lda,
	                 double*            x,
	                 IntegerForBlasType incx)
	{
		dtrsv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void trsv(char                       uplo,
	                 char                       trans,
	                 char                       diag,
	                 IntegerForBlasType         n,
	                 const std::complex<float>* a,
	                 IntegerForBlasType         lda,
	                 std::complex<float>*       x,
	                 IntegerForBlasType         incx)
	{
		ctrsv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void trsv(char                        uplo,
	                 char                        trans,
	                 char                        diag,
	                 IntegerForBlasType          n,
	                 const std::complex<double>* a,
	                 IntegerForBlasType          lda,
	                 std::complex<double>*       x,
	                 IntegerForBlasType          incx)
	{
		ztrsv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void tbsv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 IntegerForBlasType k,
	                 const float*       a,
	                 IntegerForBlasType lda,
	                 float*             x,
	                 IntegerForBlasType incx)
	{
		stbsv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void tbsv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 IntegerForBlasType k,
	                 const double*      a,
	                 IntegerForBlasType lda,
	                 double*            x,
	                 IntegerForBlasType incx)
	{
		dtbsv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void tbsv(char                       uplo,
	                 char                       trans,
	                 char                       diag,
	                 IntegerForBlasType         n,
	                 IntegerForBlasType         k,
	                 const std::complex<float>* a,
	                 IntegerForBlasType         lda,
	                 std::complex<float>*       x,
	                 IntegerForBlasType         incx)
	{
		ctbsv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void tbsv(char                        uplo,
	                 char                        trans,
	                 char                        diag,
	                 IntegerForBlasType          n,
	                 IntegerForBlasType          k,
	                 const std::complex<double>* a,
	                 IntegerForBlasType          lda,
	                 std::complex<double>*       x,
	                 IntegerForBlasType          incx)
	{
		ztbsv(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void tpsv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const float*       ap,
	                 float*             x,
	                 IntegerForBlasType incx)
	{
		stpsv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void tpsv(char               uplo,
	                 char               trans,
	                 char               diag,
	                 IntegerForBlasType n,
	                 const double*      ap,
	                 double*            x,
	                 IntegerForBlasType incx)
	{
		dtpsv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void tpsv(char                       uplo,
	                 char                       trans,
	                 char                       diag,
	                 IntegerForBlasType         n,
	                 const std::complex<float>* ap,
	                 std::complex<float>*       x,
	                 IntegerForBlasType         incx)
	{
		ctpsv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void tpsv(char                        uplo,
	                 char                        trans,
	                 char                        diag,
	                 IntegerForBlasType          n,
	                 const std::complex<double>* ap,
	                 std::complex<double>*       x,
	                 IntegerForBlasType          incx)
	{
		ztpsv(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	// ****************************************************************************
	inline void ger(IntegerForBlasType m,
	                IntegerForBlasType n,
	                const float&       alpha,
	                const float*       x,
	                IntegerForBlasType incx,
	                const float*       y,
	                IntegerForBlasType incy,
	                float*             a,
	                IntegerForBlasType lda)
	{
		sger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void ger(IntegerForBlasType m,
	                IntegerForBlasType n,
	                const double&      alpha,
	                const double*      x,
	                IntegerForBlasType incx,
	                const double*      y,
	                IntegerForBlasType incy,
	                double*            a,
	                IntegerForBlasType lda)
	{
		dger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// ****************************************************************************
	inline void geru(IntegerForBlasType         m,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType         incy,
	                 std::complex<float>*       a,
	                 IntegerForBlasType         lda)
	{
		cgeru(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void geru(IntegerForBlasType          m,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType          incy,
	                 std::complex<double>*       a,
	                 IntegerForBlasType          lda)
	{
		zgeru(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// ****************************************************************************
	inline void gerc(IntegerForBlasType         m,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType         incy,
	                 std::complex<float>*       a,
	                 IntegerForBlasType         lda)
	{
		cgerc(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void gerc(IntegerForBlasType          m,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType          incy,
	                 std::complex<double>*       a,
	                 IntegerForBlasType          lda)
	{
		zgerc(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// *****************************************************************************
	inline void her(char                       uplo,
	                IntegerForBlasType         n,
	                const std::complex<float>& alpha,
	                const std::complex<float>* x,
	                IntegerForBlasType         incx,
	                std::complex<float>*       a,
	                IntegerForBlasType         lda)
	{
		cher(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	inline void her(char                        uplo,
	                IntegerForBlasType          n,
	                const std::complex<double>& alpha,
	                const std::complex<double>* x,
	                IntegerForBlasType          incx,
	                std::complex<double>*       a,
	                IntegerForBlasType          lda)
	{
		zher(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	// *****************************************************************************
	inline void hpr(char                       uplo,
	                IntegerForBlasType         n,
	                const std::complex<float>& alpha,
	                const std::complex<float>* x,
	                IntegerForBlasType         incx,
	                std::complex<float>*       ap)
	{
		chpr(&uplo, &n, &alpha, x, &incx, ap);
	}
	inline void hpr(char                        uplo,
	                IntegerForBlasType          n,
	                const std::complex<double>& alpha,
	                const std::complex<double>* x,
	                IntegerForBlasType          incx,
	                std::complex<float>*        ap)
	{
		zhpr(&uplo, &n, &alpha, x, &incx, ap);
	}
	// *****************************************************************************
	inline void heR2(char                       uplo,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType         incy,
	                 std::complex<float>*       a,
	                 IntegerForBlasType         lda)
	{
		cher2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void heR2(char                        uplo,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType          incy,
	                 std::complex<double>*       a,
	                 IntegerForBlasType          lda)
	{
		zher2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// *****************************************************************************
	inline void hpR2(char                       uplo,
	                 IntegerForBlasType         n,
	                 const std::complex<float>& alpha,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         incx,
	                 const std::complex<float>* y,
	                 IntegerForBlasType         incy,
	                 std::complex<float>*       ap)
	{
		chpr2(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}
	inline void hpR2(char                        uplo,
	                 IntegerForBlasType          n,
	                 const std::complex<double>& alpha,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          incx,
	                 const std::complex<double>* y,
	                 IntegerForBlasType          incy,
	                 std::complex<double>*       ap)
	{
		zhpr2(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}
	// *****************************************************************************
	inline void syr(char               uplo,
	                IntegerForBlasType n,
	                const float&       alpha,
	                const float*       x,
	                IntegerForBlasType incx,
	                float*             a,
	                IntegerForBlasType lda)
	{
		ssyr(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	inline void syr(char               uplo,
	                IntegerForBlasType n,
	                const double&      alpha,
	                const double*      x,
	                IntegerForBlasType incx,
	                double*            a,
	                IntegerForBlasType lda)
	{
		dsyr(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	// ****************************************************************************
	inline void spr(char               uplo,
	                IntegerForBlasType n,
	                const float&       alpha,
	                const float*       x,
	                IntegerForBlasType incx,
	                float*             ap)
	{
		sspr(&uplo, &n, &alpha, x, &incx, ap);
	}
	inline void spr(char               uplo,
	                IntegerForBlasType n,
	                const double&      alpha,
	                const double*      x,
	                IntegerForBlasType incx,
	                double*            ap)
	{
		dspr(&uplo, &n, &alpha, x, &incx, ap);
	}
	// ****************************************************************************
	inline void syR2(char               uplo,
	                 IntegerForBlasType n,
	                 const float&       alpha,
	                 const float*       x,

	                 IntegerForBlasType incx,
	                 const float*       y,
	                 IntegerForBlasType incy,
	                 float*             a,
	                 IntegerForBlasType lda)
	{
		ssyr2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void syR2(char               uplo,
	                 IntegerForBlasType n,
	                 const double&      alpha,
	                 const double*      x,

	                 IntegerForBlasType incx,
	                 const double*      y,
	                 IntegerForBlasType incy,
	                 double*            a,
	                 IntegerForBlasType lda)
	{
		dsyr2(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// ****************************************************************************
	inline void spR2(char               uplo,
	                 IntegerForBlasType n,
	                 const float&       alpha,
	                 const float*       x,

	                 IntegerForBlasType incx,
	                 const float*       y,
	                 IntegerForBlasType incy,
	                 float*             ap)
	{
		sspr2(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}
	inline void spR2(char               uplo,
	                 IntegerForBlasType n,
	                 const double&      alpha,
	                 const double*      x,

	                 IntegerForBlasType incx,
	                 const double*      y,
	                 IntegerForBlasType incy,
	                 double*            ap)
	{
		dspr2(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}

	// ****************************************************************************

	inline void axpy(IntegerForBlasType size,
	                 const float&       a,
	                 const float*       x,
	                 IntegerForBlasType sx,
	                 float*             y,
	                 IntegerForBlasType sy)
	{
		saxpy(&size, &a, x, &sx, y, &sy);
	}

	inline void axpy(IntegerForBlasType size,
	                 const double&      a,
	                 const double*      x,
	                 IntegerForBlasType sx,
	                 double*            y,
	                 IntegerForBlasType sy)
	{
		daxpy(&size, &a, x, &sx, y, &sy);
	}

	inline void axpy(IntegerForBlasType         size,
	                 const std::complex<float>& a,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         sx,

	                 std::complex<float>* y,
	                 IntegerForBlasType   sy)
	{
		caxpy(&size, &a, x, &sx, y, &sy);
	}

	inline void axpy(IntegerForBlasType          size,
	                 const std::complex<double>& a,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          sx,

	                 std::complex<double>* y,
	                 IntegerForBlasType    sy)
	{
		zaxpy(&size, &a, x, &sx, y, &sy);
	}
	// ----------------------------------------------------------------------------
	inline void copy(IntegerForBlasType size,
	                 const float*       x,
	                 IntegerForBlasType sx,
	                 float*             y,
	                 IntegerForBlasType sy)
	{
		scopy(&size, x, &sx, y, &sy);
	}

	inline void copy(IntegerForBlasType size,
	                 const double*      x,
	                 IntegerForBlasType sx,
	                 double*            y,
	                 IntegerForBlasType sy)
	{
		dcopy(&size, x, &sx, y, &sy);
	}

	inline void copy(IntegerForBlasType         size,
	                 const std::complex<float>* x,
	                 IntegerForBlasType         sx,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         sy)
	{
		ccopy(&size, x, &sx, y, &sy);
	}

	inline void copy(IntegerForBlasType          size,
	                 const std::complex<double>* x,
	                 IntegerForBlasType          sx,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          sy)
	{
		zcopy(&size, x, &sx, y, &sy);
	}
	// ----------------------------------------------------------------------------
	inline void scal(IntegerForBlasType size, const float& a, float* y, IntegerForBlasType sy)
	{
		sscal(&size, &a, y, &sy);
	}

	inline void scal(IntegerForBlasType size, const double& a, double* y, IntegerForBlasType sy)
	{
		dscal(&size, &a, y, &sy);
	}

	inline void scal(IntegerForBlasType         size,
	                 const std::complex<float>& a,
	                 std::complex<float>*       y,
	                 IntegerForBlasType         sy)
	{
		cscal(&size, &a, y, &sy);
	}

	inline void scal(IntegerForBlasType          size,
	                 const std::complex<double>& a,
	                 std::complex<double>*       y,
	                 IntegerForBlasType          sy)
	{
		zscal(&size, &a, y, &sy);
	}
} /* namespace BLAS */
} /* namespace psimag */

#endif /* PSIMAG_BLAS */
