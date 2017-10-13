#include "Header.h"

// all source files should contain templates of functions

int main()
{

#if (PROBLEM == 0)
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

#elif (PROBLEM == 1)
	int n = N;
	double *H = new double[n*n]; // init
	double *H1 = new double[n*n]; // compressed
	double *H2 = new double[n*n]; // recovered init

	H[0:n*n] = 0;
	H1[0:n*n] = 0;
	H2[0:n*n] = 0;

	int ldh = n;
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			H[i + ldh * j] = 1.0 / (i + j + 1);
			H1[i + ldh * j] = 1.0 / (i + j + 1);
		}
	print(n, n, H1, ldh);
	Test_SymRecCompress(n, H, H1, H2, ldh);

#endif
	system("pause");

	return 0;

}