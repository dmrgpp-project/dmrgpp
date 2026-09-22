#ifndef UTIL_H
#define UTIL_H

#include "KronUtil.h"
#include <PsimagLite/Complex.h>
#include <PsimagLite/GemmR.h>
#include <PsimagLite/MatrixNonOwned.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>

template <typename ComplexOrRealType>
void estimateKronCost(const int          nrow_A,
                      const int          ncol_A,
                      const int          nnz_A,
                      const int          nrow_B,
                      const int          ncol_B,
                      const int          nnz_B,
                      ComplexOrRealType* p_kron_nnz,
                      ComplexOrRealType* p_kron_flops,
                      int*               p_imethod,
                      const typename PsimagLite::Real<ComplexOrRealType>::Type);

template <typename ComplexOrRealType>
void csrDenKronMultMethod(const int                                                   imethod,
                          const char                                                  transA,
                          const char                                                  transB,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&             a_,
                          const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                          SizeType                                                    offsetY,
                          typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                          SizeType                                                    offsetX,
                          PsimagLite::GemmR<ComplexOrRealType>&);

template <typename ComplexOrRealType>
bool csrIsEye(const PsimagLite::CrsMatrix<ComplexOrRealType>&);

template <typename ComplexOrRealType>
void csrTranspose(const int               nrow_A,
                  const int               ncol_A,
                  const int               arowptr[],
                  const int               acol[],
                  const ComplexOrRealType aval[],
                  int                     atrowptr[],
                  int                     atcol[],
                  ComplexOrRealType       atval[]);

template <typename ComplexOrRealType>
void csrKronMultMethod(const int  imethod,
                       const char transA,
                       const char transB,

                       const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>& b,

                       const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                       PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout);

template <typename ComplexOrRealType>
void csrMatmulPost(const char trans_A,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                   const int                                                  nrow_Y,
                   const int                                                  ncol_Y,
                   const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                   const int                                                  nrow_X,
                   const int                                                  ncol_X,
                   PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout);

template <typename ComplexOrRealType>
void csrMatmulPre(const char trans_A,

                  const PsimagLite::CrsMatrix<ComplexOrRealType>&,

                  const int                                                  nrow_Y,
                  const int                                                  ncol_Y,
                  const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,

                  const int                                      nrow_X,
                  const int                                      ncol_X,
                  PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout);

template <typename ComplexOrRealType>
void csrSubmatrix(const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                  const int                                       nrow_B,
                  const int                                       ncol_B,
                  const int                                       max_nnz,

                  const PsimagLite::Vector<int>::Type&      rindex,
                  const PsimagLite::Vector<int>::Type&      cindex,
                  PsimagLite::CrsMatrix<ComplexOrRealType>& b);

template <typename ComplexOrRealType>
void csrEye(const int nrow_B, const int ncol_B, PsimagLite::CrsMatrix<ComplexOrRealType>& b);

template <typename ComplexOrRealType>
void csrKronSubmatrix(const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                      const PsimagLite::CrsMatrix<ComplexOrRealType>& b,
                      const int                                       nrindex,
                      const int                                       ncindex,
                      const int                                       max_nnz,
                      const PsimagLite::Vector<int>::Type&            rindex,
                      const PsimagLite::Vector<int>::Type&            cindex,
                      PsimagLite::CrsMatrix<ComplexOrRealType>&       h);

int cscNnz(const int ncol_A, const PsimagLite::Vector<int>::Type& acolptr);

template <typename ComplexOrRealType>
void cscMatmulPre(const char                                                  trans_A,
                  const int                                                   nrow_A,
                  const int                                                   ncol_A,
                  const PsimagLite::Vector<int>::Type&                        acolptr,
                  const PsimagLite::Vector<int>::Type&                        arow,
                  const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                  const int                                                   nrow_Y,
                  const int                                                   ncol_Y,
                  const PsimagLite::Matrix<ComplexOrRealType>&                yin,
                  const int                                                   nrow_X,
                  const int                                                   ncol_X,
                  PsimagLite::Matrix<ComplexOrRealType>&                      xout);

template <typename ComplexOrRealType>
void cscMatmulPost(const char                                                  trans_A,
                   const int                                                   nrow_A,
                   const int                                                   ncol_A,
                   const PsimagLite::Vector<int>::Type&                        acolptr,
                   const PsimagLite::Vector<int>::Type&                        arow,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                   const int                                                   nrow_Y,
                   const int                                                   ncol_Y,
                   const PsimagLite::Matrix<ComplexOrRealType>&                yin,
                   const int                                                   nrow_X,
                   const int                                                   ncol_X,
                   PsimagLite::Matrix<ComplexOrRealType>&                      xout);

template <typename ComplexOrRealType>
void cscKronMultMethod(const int                                                   imethod,
                       const int                                                   nrow_A,
                       const int                                                   ncol_A,
                       const PsimagLite::Vector<int>::Type&                        acolptr,
                       const PsimagLite::Vector<int>::Type&                        arow,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                       const int                                                   nrow_B,
                       const int                                                   ncol_B,
                       const PsimagLite::Vector<int>::Type&                        bcolptr,
                       const PsimagLite::Vector<int>::Type&                        brow,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,
                       const PsimagLite::Matrix<ComplexOrRealType>&                yin,
                       PsimagLite::Matrix<ComplexOrRealType>&                      xout);

template <typename ComplexOrRealType>
void cscKronMult(const int                                                   nrow_A,
                 const int                                                   ncol_A,
                 const PsimagLite::Vector<int>::Type&                        acolptr,
                 const PsimagLite::Vector<int>::Type&                        arow,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                 const int                                                   nrow_B,
                 const int                                                   ncol_B,
                 const PsimagLite::Vector<int>::Type&                        bcolptr,
                 const PsimagLite::Vector<int>::Type&                        brow,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,
                 const PsimagLite::Matrix<ComplexOrRealType>&                yin,
                 PsimagLite::Matrix<ComplexOrRealType>&                      xout,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type);

template <typename ComplexOrRealType>
void coord2csr(const int               nrow_A,
               const int               ncol_A,
               const int               nnz,
               const int               ilist[],
               const int               jlist[],
               const ComplexOrRealType alist[],
               int                     arowptr[],
               int                     acol[],
               ComplexOrRealType       aval[]);

template <typename ComplexOrRealType>
void denCsrKronMultMethod(const int                                                   imethod,
                          const char                                                  transA,
                          const char                                                  transB,
                          const PsimagLite::Matrix<ComplexOrRealType>&                a_,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&             b,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                          SizeType                                                    offsetY,
                          typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                          SizeType                                                    offsetX,
                          PsimagLite::GemmR<ComplexOrRealType>&);

void denCopymat(const int nrow, const int ncol, const int asrc_[], int bdest_[]);

template <typename ComplexOrRealType>
void denZeros(const int nrow_A, const int ncol_A, PsimagLite::Matrix<ComplexOrRealType>& a_);

template <typename ComplexOrRealType>
void denTranspose(const int               nrow_A,
                  const int               ncol_A,
                  const ComplexOrRealType a_[],
                  ComplexOrRealType       at_[]);

template <typename ComplexOrRealType>
void denGenMatrix(const int                                                 nrow_A,
                  const int                                                 ncol_A,
                  const typename PsimagLite::Real<ComplexOrRealType>::Type& threshold,
                  PsimagLite::Matrix<ComplexOrRealType>&                    a_);

template <typename ComplexOrRealType>
void denMatmulPre(const char                                                 trans_A,
                  const int                                                  nrow_A,
                  const int                                                  ncol_A,
                  const PsimagLite::Matrix<ComplexOrRealType>&               a_,
                  const int                                                  nrow_Y,
                  const int                                                  ncol_Y,
                  const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                  const int                                                  nrow_X,
                  const int                                                  ncol_X,
                  PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout,
                  PsimagLite::GemmR<ComplexOrRealType>&);

template <typename ComplexOrRealType>
void denMatmulPost(const char                                                 trans_A,
                   const int                                                  nrow_A,
                   const int                                                  ncol_A,
                   const PsimagLite::Matrix<ComplexOrRealType>&               a_,
                   const int                                                  nrow_Y,
                   const int                                                  ncol_Y,
                   const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                   const int                                                  nrow_X,
                   const int                                                  ncol_X,
                   PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout,
                   PsimagLite::GemmR<ComplexOrRealType>&);

template <typename ComplexOrRealType>
void denKronSubmatrix(const int                                    nrow_A,
                      const int                                    ncol_A,
                      const PsimagLite::Matrix<ComplexOrRealType>& a_,
                      const int                                    nrow_B,
                      const int                                    ncol_B,
                      const PsimagLite::Matrix<ComplexOrRealType>& b_,
                      const int                                    nrindex,
                      const int                                    ncindex,
                      const PsimagLite::Vector<int>::Type&         rindex,
                      const PsimagLite::Vector<int>::Type&         cindex,
                      PsimagLite::Matrix<ComplexOrRealType>&       c_);

template <typename ComplexOrRealType>
void denKronMultMethod(const int                                                   imethod,
                       const char                                                  transA,
                       const char                                                  transB,
                       const PsimagLite::Matrix<ComplexOrRealType>&                a_,
                       const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                       SizeType                                                    offsetY,
                       typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                       SizeType                                                    offsetX,
                       PsimagLite::GemmR<ComplexOrRealType>&);

template <typename ComplexOrRealType> int denNnz(const PsimagLite::Matrix<ComplexOrRealType>&);

template <typename ComplexOrRealType> bool denIsEye(const PsimagLite::Matrix<ComplexOrRealType>&);

template <typename ComplexOrRealType> bool denIsZeros(const PsimagLite::Matrix<ComplexOrRealType>&);

template <typename ComplexOrRealType>
void denKronForm(const int                                    nrow_A,
                 const int                                    ncol_A,
                 const PsimagLite::Matrix<ComplexOrRealType>& a_,
                 const int                                    nrow_B,
                 const int                                    ncol_B,
                 const PsimagLite::Matrix<ComplexOrRealType>& b_,
                 PsimagLite::Matrix<ComplexOrRealType>&       c_);

template <typename ComplexOrRealType>
void denKronFormGeneral(const char                                   transA,
                        const char                                   transB,
                        const int                                    nrow_A,
                        const int                                    ncol_A,
                        const PsimagLite::Matrix<ComplexOrRealType>& a_,
                        const int                                    nrow_B,
                        const int                                    ncol_B,
                        const PsimagLite::Matrix<ComplexOrRealType>& b_,
                        PsimagLite::Matrix<ComplexOrRealType>&       c_);

template <typename ComplexOrRealType>
void denSubmatrix(const int                                    nrow_A,
                  const int                                    ncol_A,
                  const PsimagLite::Matrix<ComplexOrRealType>& a_,
                  const int                                    nrindex,
                  const int                                    ncindex,
                  const PsimagLite::Vector<int>::Type&         rindex,
                  const PsimagLite::Vector<int>::Type&         cindex,
                  PsimagLite::Matrix<ComplexOrRealType>&       c_);

template <typename ComplexOrRealType>
void denEye(const int nrow_A, const int ncol_A, PsimagLite::Matrix<ComplexOrRealType>& c_);

template <typename ComplexOrRealType>
bool csrIsZeros(const PsimagLite::CrsMatrix<ComplexOrRealType>& a)
{
	// ----------------------------------------------------
	// check whether a sparse matrix is the zero matrix
	// ----------------------------------------------------
	return isZero(a, 0.0);
}

template <typename ComplexOrRealType> int csrNnz(const PsimagLite::CrsMatrix<ComplexOrRealType>& a)
{
	return a.nonZeros();
}

#endif
