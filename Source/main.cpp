#include "Header.h"
#include "templates.h"

// all source files should contain templates of functions

int main()
{

#if (PROBLEM == 1)
	int n = 6;
	double eps = 1e-2;
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for(double eps = 1e-2; eps > 1e-8; eps /= 10)
		Test_SymRecCompress(n, eps, method, smallsize);
	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
			Test_SymRecCompressStruct(n, eps, method, smallsize);

#elif (PROBLEM == 2)
	int n1 = 39; // number of point across the directions
	int n2 = 39;
	int n3 = 39;
	int n = n1 * n2; // size of blocks
	int NB = n3; // number of blocks
	int size = n * NB; // size of vector x and f: n1 * n2 * n3
	int smallsize = 400;
	double thresh = 1e-6; // stop level of algorithm by relative error
	int ItRef = 200; // Maximal number of iterations in refirement
	char bench[255] = "display"; // parameter into solver to show internal results
	int sparse_size = n + 2 * (n - 1) + 2 * (n - n1);

	// the size of matrix A: n^3 * n^3 = n^6
	// here we have n^2 * n = n^3 
	// init
	double *D = alloc_arr(size * n); // it's a matrix with size n^3 * n^2 = size * n
	double *B = alloc_arr(size - n); // vector of diagonal elements
	double *x_orig = alloc_arr(size);
	double *f = alloc_arr(size);

	// result obtained
	double *G = alloc_arr(size * n);
	double *x_sol = alloc_arr(size);

	int ldd = size;
	int ldg = size;

	int success = 0;
	int itcount = 0;
	double RelRes = 0;
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
		
	size_m x, y, z;

	x.n = n1;
	y.n = n2;
	z.n = n3;

	x.l = y.l = z.l = 1;
	x.h = x.l / (double)(x.n + 1);
	y.h = y.l / (double)(y.n + 1);
	z.h = z.l / (double)(z.n + 1);

	printf("%lf %lf %lf\n", x.h, y.h, z.h);

	GenMatrixandRHSandSolution2(x, y, z, D, ldd, B, x_orig, f);

	printf("Solving %d x %d x %d Laplace equation\n", n1, n2, n3);
	printf("The system has %d diagonal blocks of size %d x %d\n", n3, n1*n2, n1*n2);
	printf("Compressed blocks method\n");
	printf("Parameters: thresh=%g, smallsize=%d \n", thresh, smallsize);

	// Calling the solver
//	Block3DSPDSolveFast(n1, n2, n3, D, ldd, B, f, thresh, smallsize, ItRef, bench, G, ldg, x_sol, success, RelRes, itcount);
	mnode **Gstr;
	Block3DSPDSolveFastStruct(n1, n2, n3, D, ldd, B, f, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
	printf("success=%d, itcount=%d\n", success, itcount);
	printf("-----------------------------------\n");

	//print_vec(size, x_orig, x_sol, "x_orig and x_sol");

	printf("Computing error ||x_{exact}-x_{comp}||/||x_{exact}||\n");
	rel_error(n, 1, x_sol, x_orig, size, thresh);


#endif

	free_arr(&D);
	free_arr(&B);
	free_arr(&x_orig);
	free_arr(&x_sol);
	free_arr(&f);

#elif (PROBLEM == 3)
	int n = 6;
	double eps = 1e-2;
	char method[255] = "SVD";
	int smallsize = 3;

	// Test compress relative error of Hd and H2
	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
			Test_DiagMult(n, eps, method, smallsize);

	printf("--------------------Structure-------------\n");
	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
			Test_DiagMultStruct(n, eps, method, smallsize);

#elif (PROBLEM == 4)
///	int n = 5;
	//int k = 6;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	// Test compress relative error of Y1(recovered) and Y (init)
	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int k = 1; k <= 10; k++)
				Test_RecMultL(n, k, eps, method, smallsize);

	printf("--------------------Structure-------------\n");
	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int k = 1; k <= 10; k++)
				Test_RecMultLStruct(n, k, eps, method, smallsize);
#elif (PROBLEM == 5)
	int m = 8;
	int n = 10;
	double eps = 1e-5;
	char method[255] = "SVD";

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int m = 3; m <= 10; m++)
			for (int n = 1; n <= 10; n++)
				Test_LowRankApprox(m, n, eps, method);

	cout << "----------------Structured approach---------------\n" << endl;
	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int m = 3; m <= 10; m++)
			for (int n = 1; n <= 10; n++)
				Test_LowRankApproxStruct(m, n, eps, method);

#elif (PROBLEM == 6)
	int n = 7;
	double eps = 1e-2;
	char method[255] = "SVD";
	int smallsize = 3;
	double alpha = 1.0;
	double beta = -1.0;

	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
			Test_Add(n, alpha, beta, smallsize, eps, method);
	printf("-----------------Structure----------------\n");
	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
			Test_AddStruct(n, alpha, beta, smallsize, eps, method);

#elif (PROBLEM == 7)
	int m = 3;
	int n = 4;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int m = 2; m <= 10; m++)
				Test_Transpose(m, n, smallsize, eps, method);

#elif (PROBLEM == 8)

	//int n = 3;
	//int k = 1;
	double alpha = -1.3;
	double eps = 1e-4;
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int k = 1; k <= 10; k++)
				Test_SymCompUpdate2(n, k, alpha, smallsize, eps, method);

	printf("----------------Structure----------------\n");
	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int k = 1; k <= 10; k++)
				Test_SymCompUpdate2Struct(n, k, alpha, smallsize, eps, method);

#elif (PROBLEM == 9) // test for inversion of compressed matrix

	int n = 5;
	double eps = 1e-06;
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			Test_SymCompRecInv(n, smallsize, eps, method);

	printf("----------------------Structure---------------\n");
	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			Test_SymCompRecInvStruct(n, smallsize, eps, method);

#elif (PROBLEM == 10)
	int n = 5;
	double eps = 1e-06;
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			Test_CopyStruct(n, eps, method, smallsize);

#elif (PROBLEM == 11)
	
	node *root1 = NULL;
	node *root2 = NULL;

	root1 = insert1(root1, 5);
	root2 = insert1(root2, 5);
	insert2(&root2, 7);
	insert2(&root2, 3);
	insert2(&root2, 2);
	insert2(&root2, 8);
	insert2(&root2, 4);
	insert2(&root2, 1);
	insert2(&root2, 6);
	insert2(&root2, 9);

//	cout << "val1 = " << root1->val << endl;
	cout << "val2 = " << root2->val << endl;

	printf("size = %d\n", TreeSize(root2));
	printf("maxDepth = %d\n", MaxDepth(root2));
	printf("minValue = %d\n", MinValue(root2));
	printf("maxValue = %d\n", MaxValue(root2));
	PrintInorder(root2);
	printf("\n");
	PrintPostorder(root2);
	printf("\n");
	printf("IsBST = %d\n", IsBST(root2));


/*
//	bool is_node;
	is_node = lookup(root2, 3);
	cout << is_node << endl;
	is_node = lookup(root2, 5);
	cout << is_node << endl;*/


#endif
	system("pause");

	return 0;

}