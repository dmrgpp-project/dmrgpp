#ifndef KRON_UTIL_HEADER_H
#define KRON_UTIL_HEADER_H
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/GemmR.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Vector.h>

#include <complex>

template <typename ComplexOrRealType>
void csrKronMult(const char                                                  transA,
                 const char                                                  transB,
                 const PsimagLite::CrsMatrix<ComplexOrRealType>&             a,
                 const PsimagLite::CrsMatrix<ComplexOrRealType>&             b,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                 SizeType                                                    offsetY,
                 typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                 SizeType                                                    offsetX,
                 const typename PsimagLite::Real<ComplexOrRealType>::Type);

//-----------------------------------------------------------------------------------

template <typename ComplexOrRealType>
void denCsrKronMult(const char                                   transA,
                    const char                                   transB,
                    const PsimagLite::Matrix<ComplexOrRealType>& a_,
                    const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                    SizeType                                                    offsetY,
                    typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                    SizeType                                                    offsetX,
                    const typename PsimagLite::Real<ComplexOrRealType>::Type,
                    PsimagLite::GemmR<ComplexOrRealType>&);

//-----------------------------------------------------------------------------------

template <typename ComplexOrRealType>
void denKronMult(const char                                                  transA,
                 const char                                                  transB,
                 const PsimagLite::Matrix<ComplexOrRealType>&                a_,
                 const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                 SizeType                                                    offsetY,
                 typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                 SizeType                                                    offsetX,
                 const typename PsimagLite::Real<ComplexOrRealType>::Type,
                 PsimagLite::GemmR<ComplexOrRealType>&);

//-----------------------------------------------------------------------------------

template <typename ComplexOrRealType>
void csrDenKronMult(const char transA,
                    const char transB,
                    const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                    const PsimagLite::Matrix<ComplexOrRealType>&                b_,
                    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                    SizeType                                                    offsetY,
                    typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                    SizeType                                                    offsetX,
                    const typename PsimagLite::Real<ComplexOrRealType>::Type,
                    PsimagLite::GemmR<ComplexOrRealType>&);
#endif
