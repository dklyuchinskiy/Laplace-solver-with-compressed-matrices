#include "Header.h"

// all source files should contain templates of functions

int main()
{
	int n2 = N * N;
	int n3 = 3 * N;
	int nbl = N;
	int ldA = n2;

	double *A_f = alloc_arr(n2 * n2);
	double *L_f = alloc_arr(n2 * n2);
	double *U_f = alloc_arr(n2 * n2);

	A_f[0:n2*n2] = 0;
	U_f[0:n2*n2] = 0;
	L_f[0:n2*n2] = 0;

	int ldaf, lduf, ldlf;
	ldaf = lduf = ldlf = n2;

	int **ipiv_mat = (int**)malloc(nbl * sizeof(int*));
	for (int j = 0; j < nbl; j++)
		ipiv_mat[j] = (int*)malloc(nbl * sizeof(int));

	// Matrix with non-zero elements (dense)
	double *A = alloc_arr(n2 * n3);

	// Initialize dense matrix A
	Init_matrix(nbl, n2, n3, A, ldA);

	// Construct sparse matrix A
	construct_A(nbl, n2, n3, A, ldA, A_f, ldaf);

	// Factorize LU in dense matrix A
	sparse_lu(nbl, n2, n3, A, ldA, ipiv_mat);

	// Construct sparse matrix L
	construct_L(nbl, n2, n3, A, ldA, L_f, ldlf);

	// Construct sparse matrix U
	construct_U(nbl, n2, n3, A, ldA, U_f, lduf);

	// Test sparse residual || A - L * U ||
	test(nbl, n2, n3, A_f, ldaf, U_f, lduf, L_f, ldlf, ipiv_mat);


	free_arr(&A_f);
	free_arr(&L_f);
	free_arr(&U_f);
	free_arr(&A);

	for (int j = 0; j < nbl; j++)
		free(ipiv_mat[j]);

	free(ipiv_mat);

	system("pause");

	return 0;

}