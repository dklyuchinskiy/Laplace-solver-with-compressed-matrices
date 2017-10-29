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
	int n = 10;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		Test_SymRecCompress(n, eps, method, smallsize);

#elif (PROBLEM == 2)
	int n1 = 4; // number of point across the directions
	int n2 = 4;
	int n3 = 4;
	int n = n1 * n2; // size of blocks
	int NB = n3; // number of blocks
	int size = n * NB; // size of vector x and f: n1 * n2 * n3
	int smallsize = 400;
	double thresh = 1e-6; // stop level of algorithm by relative error
	int ItRef = 100; // Maximal number of iterations in refirement
	char bench[255] = "display"; // parameter into solver to show internal results
	int sparse_size = n + 2 * (n - 1) + 2 * (n - n1);

	// init
	double *D = alloc_arr(sparse_size * NB);
	double *B = alloc_arr(size - NB); // vector of diagonal elements
	double *x1 = alloc_arr(size);
	double *f = alloc_arr(size);

	// result obtained
	double *G = alloc_arr(1);
	double *x = alloc_arr(size);

	int ldd = 0;
	int ldb = 0;
	int ldg = 0;

	int success = 0;
	int itcount = 0;
	double RelRes = 0;

	// Generation matrix of coefficients, vector of solution (to compare with obtained) and vector of RHS
	GenMatrixandRHSandSolution(n1, n2, n3, D, ldd, B, ldb, x1, f);

	print_vec(size, x1, f, "vector x1 and f");

	printf("Solving %d x %d x %d Laplace equation\n", n1, n2, n3);
	printf("The system has %d diagonal blocks of size %d x %d\n", n3, n1*n2, n1*n2);
	printf("Compressed blocks method\n");
	printf("Parameters: thresh=%g, smallsize=%d \n", thresh, smallsize);

	// Calling the solver
	Block3DSPDSolveFast(n1, n2, n3, D, ldd, B, ldb, f, thresh, smallsize, ItRef, bench, G, ldg, x, success, RelRes, itcount);

	printf("success=%d, itcount=%d\n", success, itcount);
	printf("-----------------------------------\n");
	printf("Computing error\n");
	
	double normx = 0;
	double normy = 0;
	int ione = 1;


/*
#pragma omp parallel for
	for (int i = 0; i < size; i++)
		x[i] -= x1[i];

	normy = dnrm2(&size, x, &ione);
	normx = dnrm2(&size, x1, &ione);

	printf("relative error %lf\n", normy/normx);*/

#elif (PROBLEM == 3)
	int n = 10;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	// Test compress relative error of Hd and H2
	for (int n = 3; n <= 10; n++)
		Test_DiagMult(n, eps, method, smallsize);

#elif (PROBLEM == 4)
///	int n = 5;
	//int k = 6;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	// Test compress relative error of Y1(recovered) and Y (init)
	for (int n = 3; n <= 10; n++)
		for (int k = 1; k <= 10; k++)
			Test_RecMultL(n, k, eps, method, smallsize);
#elif (PROBLEM == 5)
	int m = 8;
	int n = 10;
	double eps = 1e-5;
	char method[255] = "SVD";
	for (int i = 0; i < TEST_IT; i++)
	{
		Test_LowRankApprox_InitA(m, n, eps, method);
	}
#elif (PROBLEM == 6)
	int n = 5;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;
	double alpha = 1.0;
	double beta = -1.0;

	for (int n = 3; n <= 10; n++)
		Test_add(n, alpha, beta, smallsize, eps, method);

#elif (PROBLEM == 7)
	int m = 3;
	int n = 4;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for (int m = 2; m <= 10; m++)
			Test_transpose(m, n, smallsize, eps, method);

#elif (PROBLEM == 8)

	//int n = 3;
	//int k = 1;
	double alpha = -1.3;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for (int k = 1; k <= 10; k++)
			Test_SymCompUpdate2(n, k, alpha, smallsize, eps, method);

#elif (PROBLEM == 9) // test for inversion of compressed matrix

	int n = 5;
	double eps = 1e-06;
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		Test_SymCompRecInv(n, smallsize, eps, method);

#endif
	system("pause");

	return 0;

}