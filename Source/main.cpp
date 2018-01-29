#include "definitions.h"
#include "templates.h"
#include "TestSuite.h"

// all source files should contain templates of functions

int main()
{
	TestAll();

	int n1 = 39; // number of point across the directions
	int n2 = 39;
	int n3 = 39;
	int n = n1 * n2; // size of blocks
	int NB = n3; // number of blocks
	int size = n * NB; // size of vector x and f: n1 * n2 * n3
	int smallsize = 300;
	double thresh = 1e-6; // stop level of algorithm by relative error
	int ItRef = 200; // Maximal number of iterations in refirement
	char bench[255] = "no"; // parameter into solver to show internal results
	int sparse_size = n + 2 * (n - 1) + 2 * (n - n1);
	int non_zeros_in_3diag = n + (n - 1) * 2 + (n - n1) * 2 - (n1 - 1) * 2;

	size_m x, y, z;

	x.n = n1;
	y.n = n2;
	z.n = n3;

	x.l = y.l = z.l = 40;
	x.h = x.l / (double)(x.n + 1);
	y.h = y.l / (double)(y.n + 1);
	z.h = z.l / (double)(z.n + 1);


	// the size of matrix A: n^3 * n^3 = n^6
	// here we have n^2 * n = n^3 
	// init

#ifndef ONLINE
	double *D = alloc_arr(size * n); // it's a matrix with size n^3 * n^2 = size * n
	double *B_mat = alloc_arr((size - n) * n); int ldb = size - n;
	int ldd = size;
	int ldg = size;
#else
	double *D = alloc_arr(n * n); // it's a matrix with size n^3 * n^2 = size * n
	double *B_mat = alloc_arr(n * n);

	int ldd = n;
	int ldb = n;
#endif

	double *B = alloc_arr(size - n); // vector of diagonal elementes
	double *x_orig = alloc_arr(size);
	double *x_sol = alloc_arr(size);
	double *f = alloc_arr(size);

#ifdef CSR_FORMAT
	// Memory for CSR matrix
	dcsr *Dcsr;
	int non_zeros_in_block3diag = (n + (n - 1) * 2 + (n - x.n) * 2 - (x.n - 1) * 2) * z.n + 2 * (size - n);
	Dcsr = (dcsr*)malloc(sizeof(dcsr));
	Dcsr->values = (double*)malloc(non_zeros_in_block3diag * sizeof(double));
	Dcsr->ia = (int*)malloc((size + 1) * sizeof(int));
	Dcsr->ja = (int*)malloc(non_zeros_in_block3diag * sizeof(int));
	Dcsr->ia[size] = non_zeros_in_block3diag + 1;
#endif

	int success = 0;
	int itcount = 0;
	double RelRes = 0;
	double norm = 0;
	int nthr;

#pragma omp parallel
{
	nthr = omp_get_num_threads();
	int ithr = omp_get_thread_num();
	if (ithr == 0) printf("Run in parallel on %d threads\n", nthr);
}
#if (TEST == 0)
	// Generation matrix of coefficients, vector of solution (to compare with obtained) and vector of RHS
	GenMatrixandRHSandSolution(n1, n2, n3, D, ldd, B, x_orig, f);

	printf("Solving %d x %d x %d Laplace equation\n", n1, n2, n3);
	printf("The system has %d diagonal blocks of size %d x %d\n", n3, n1*n2, n1*n2);
	printf("Compressed blocks method\n");
	printf("Parameters: thresh=%g, smallsize=%d \n", thresh, smallsize);

	// Calling the solver
	Block3DSPDSolveFast(n1, n2, n3, D, ldd, B, f, thresh, smallsize, ItRef, bench, G, ldg, x_sol, success, RelRes, itcount);

	printf("success=%d, itcount=%d\n", success, itcount);
	printf("-----------------------------------\n");
	
	//print_vec(size, x_orig, x_sol, "x_orig and x_sol");

	printf("Computing error ||x_{comp}-x_{orig}||/||x_{orig}||\n");
	rel_error(n, 1, x_sol, x_orig, size, thresh);
#elif (TEST == 1)
		
	printf("%lf %lf %lf\n", x.h, y.h, z.h);

#ifndef ONLINE
	GenMatrixandRHSandSolution2(x, y, z, D, ldd, B, x_orig, f, thresh);
#else
	GenRHSandSolution(x, y, z, x_orig, f);
	// Set vector B
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < z.n - 1; j++)
		for (int i = 0; i < n; i++)
			B[ind(j, n) + i] = 1.0 / (z.h * z.h);
#endif

#ifdef CSR_FORMAT
#ifndef ONLINE
	GenSparseMatrix(x, y, z, B_mat, ldb, D, ldd, B_mat, ldb, Dcsr);
#else
	GenSparseMatrix2(x, y, z, B_mat, n, D, n, B_mat, n, Dcsr);
#endif
	free_arr(&B_mat);
#endif

#ifdef ONLINE
	free_arr(&D);
#endif

#ifdef CSR_FORMAT
	//	Test_CompareColumnsOfMatrix(n1, n2, n3, D, ldd, B, Dcsr, thresh);
	Test_TransferBlock3Diag_to_CSR(n1, n2, n3, Dcsr, x_orig, f, thresh);
#endif

#if 0
	printf("non_zeros: %d\n", non_zeros_in_3diag);
	double *values = alloc_arr(non_zeros_in_3diag);
	int *ia = (int*)malloc(non_zeros_in_3diag * sizeof(int));
	int *ja = (int*)malloc(non_zeros_in_3diag * sizeof(int));
	map<vector<int>, double> CSR;
	CSR = dense_to_CSR(n, n, &D[0], ldd, ia, ja, values);
	print(n, n, &D[0], ldd, "D[0]");
	print_map(CSR);
	free(ia);
	free(ja);
	free(values);
	system("pause");
#else

//	print(n, n, &B_mat[0], ldb, "B[0]");
//	print(size, n, D, ldd, "D");

//	printf("all non_zero elements: %d\n", non_zeros_in_block3diag);
//	for (int i = 0; i < size + 1; i ++)
//	printf("%d: ia = %d value(ia) = %lf diff = %d\n", i, Dcsr->ia[i], Dcsr->values[Dcsr->ia[i] - 1], Dcsr->ia[i+1]- Dcsr->ia[i]);

//	print_vec(non_zeros_in_block3diag, Dcsr->ja, Dcsr->values, "ja and values");
//	print_map(CSR);

#endif
#if 1
	printf("non_zeros: %d\n", non_zeros_in_block3diag);
	printf("Solving %d x %d x %d Laplace equation\n", n1, n2, n3);
	printf("The system has %d diagonal blocks of size %d x %d\n", n3, n1*n2, n1*n2);
	printf("Compressed blocks method\n");
	printf("Parameters: thresh = %g, smallsize = %d \n", thresh, smallsize);

	// Calling the solver
	
#ifndef STRUCT
	double *G = alloc_arr(size * n);
	Block3DSPDSolveFast(n1, n2, n3, D, ldd, B, f, thresh, smallsize, ItRef, bench, G, ldg, x_sol, success, RelRes, itcount);
#else
	mnode **Gstr;
#ifndef ONLINE
	Block3DSPDSolveFastStruct(x, y, z, D, ldd, B, f, Dcsr, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
#else
	Block3DSPDSolveFastStruct(x, y, z, NULL, ldd, B, f, Dcsr, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
#endif

#endif
	printf("success = %d, itcount = %d\n", success, itcount);
	printf("-----------------------------------\n");

	//print_vec(size, x_orig, x_sol, "x_orig and x_sol");

	printf("Computing error ||x_{exact}-x_{comp}||/||x_{exact}||\n");
	norm = rel_error(n, 1, x_sol, x_orig, size, thresh);

	if (norm < thresh) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, thresh);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, thresh);

	printf("----------Trees information-----------\n");

	for (int i = z.n - 1; i >= 0; i--)
	{
		printf("For block %2d. Size: %d, MaxDepth: %d, Ranks: ", i, TreeSize(Gstr[i]), MaxDepth(Gstr[i]));
		PrintRanksInWidthList(Gstr[i]);
		printf("\n");
	}

	for (int i = z.n - 1; i >= 0; i--)
		FreeNodes(n, Gstr[i], smallsize);

	free(Gstr);

#endif
#endif

#ifndef ONLINE
	free_arr(&D);
	free_arr(&B);
#endif
	free_arr(&x_orig);
	free_arr(&x_sol);
	free_arr(&f);

	system("pause");

	return 0;

}