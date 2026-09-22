#include "util.h"

template <typename ComplexOrRealType>
void denKronFormGeneral(const char                                   transA,
                        const char                                   transB,
                        const int                                    nrow_A,
                        const int                                    ncol_A,
                        const PsimagLite::Matrix<ComplexOrRealType>& a_,
                        const int                                    nrow_B,
                        const int                                    ncol_B,
                        const PsimagLite::Matrix<ComplexOrRealType>& b_,
                        PsimagLite::Matrix<ComplexOrRealType>&       c_)
{
	/*
	 * ---------------------------------------
	 * form C = kron(A,B),  where
	 * nrow_C = nrow_A * nrow_B
	 * ncol_C = ncol_A * ncol_B
	 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
	 * ---------------------------------------
	 */
	const int  idebug     = 0;
	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;

	const bool istrans_a      = (transA == 'T') || (transA == 't');
	const bool istrans_b      = (transB == 'T') || (transB == 't');
	const bool is_conjtrans_a = (transA == 'C') || (transA == 'c');
	const bool is_conjtrans_b = (transB == 'C') || (transB == 'c');

	int ia = 0;
	int ja = 0;
	int ib = 0;
	int jb = 0;

	int nrow_1 = (istrans_a || is_conjtrans_a) ? ncol_A : nrow_A;
	int nrow_2 = (istrans_b || is_conjtrans_b) ? ncol_B : nrow_B;

	int ncol_1 = (istrans_a || is_conjtrans_a) ? nrow_A : ncol_A;
	int ncol_2 = (istrans_b || is_conjtrans_b) ? nrow_B : ncol_B;

	int nrow_c = c_.rows();
	int ncol_c = c_.cols();

	assert(nrow_1 * nrow_2 == nrow_C);
	assert(ncol_1 * ncol_2 == ncol_C);

	if (idebug >= 1) {
		printf("den_kron_form_general: transA=%c, transB=%c, nrow_A=%d,ncol_A=%d\n",
		       transA,
		       transB,
		       nrow_A,
		       ncol_A);
		printf("nrow_B=%d, ncol_B=%d,    nrow_C=%d, ncol_C=%d\n",
		       nrow_B,
		       ncol_B,
		       nrow_c,
		       ncol_c);
		printf("(nrow_1,nrow_2) = (%d,%d), (ncol_1,ncol2) = (%d,%d) \n",
		       nrow_1,
		       nrow_2,
		       ncol_1,
		       ncol_2);
	};

	for (ja = 0; ja < ncol_A; ja++) {
		for (jb = 0; jb < ncol_B; jb++) {
			for (ia = 0; ia < nrow_A; ia++) {
				for (ib = 0; ib < nrow_B; ib++) {
					ComplexOrRealType aij = a_(ia, ja);
					ComplexOrRealType bij = b_(ib, jb);
					if (is_complex && is_conjtrans_a) {
						aij = PsimagLite::conj(aij);
					};
					if (is_complex && is_conjtrans_b) {
						bij = PsimagLite::conj(bij);
					};

					int iia = (istrans_a || is_conjtrans_a) ? ja : ia;
					int jja = (istrans_a || is_conjtrans_a) ? ia : ja;

					int iib = (istrans_b || is_conjtrans_b) ? jb : ib;
					int jjb = (istrans_b || is_conjtrans_b) ? ib : jb;

					// ------------------------------
					// note index for B varies faster
					// ------------------------------
					int ic = iib + iia * nrow_2;
					int jc = jjb + jja * ncol_2;

					if (idebug >= 2) {
						printf("(ia,ja)=(%d,%d), (ib,jb)=(%d,%d) ",
						       ia,
						       ja,
						       ib,
						       jb);
						printf("(iia,jja)=(%d,%d), (iib,jjb)=(%d,%d) ",
						       iia,
						       jja,
						       iib,
						       jjb);
						printf("(ic,jc)=(%d,%d) \n", ic, jc);
					};

					c_(ic, jc) = aij * bij;
				};
			};
		};
	};
}
