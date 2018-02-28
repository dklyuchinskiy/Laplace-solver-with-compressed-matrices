#include "definitions.h"
#include "templates.h"
#include "TestSuite.h"
#include "TestFramework.h"

/***************************************
Source file contains tests for testing
functionalites, implemented in
functions.cpp and BinaryTrees.cpp.

The interface is declared in TestSuite.h
****************************************/

void TestAll()
{
	TestRunner runner;
	//void(*pt_func)(int&) = NULL;
	//pt_func = &Shell_SymRecCompress;

	printf("***** TEST LIBRARY FUNCTIONS *******\n");

	runner.RunTests(Shell_LowRankApprox, Test_LowRankApprox, Test_LowRankApproxStruct, "Test_LowRankApprox");
	runner.RunTests(Shell_SymRecCompress, Test_SymRecCompress, Test_SymRecCompressStruct, "Test_SymRecCompress");
	runner.RunTests(Shell_DiagMult, Test_DiagMult, Test_DiagMultStruct, "Test_DiagMult");
	runner.RunTests(Shell_RecMultL, Test_RecMultL, Test_RecMultLStruct, "Test_RecMultL");
	runner.RunTests(Shell_Add, Test_Add, Test_AddStruct, "Test_Add");
	runner.RunTests(Shell_SymCompUpdate2, Test_SymCompUpdate2, Test_SymCompUpdate2Struct, "Test_SymCompUpdate2");
	runner.RunTests(Shell_SymCompRecInv, Test_SymCompRecInv, Test_SymCompRecInvStruct, "Test_SymCompRecInv");
	runner.RunTest(Shell_CopyStruct, Test_CopyStruct,  "Test_CopyStruct");
	runner.RunTest(Shell_LowRankApproxTranspose, Test_LowRankApproxTranspose, "Test_LowRankApproxTranspose");

	printf("********************\n");
	printf("ALL TESTS: %d\nPASSED: %d \nFAILED: %d\n", runner.GetAll(), runner.GetPassed(), runner.GetFailed());

	printf("***** THE END OF TESTING*******\n\n");

}

void Shell_LowRankApprox(ptr_test_low_rank func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int m = 3; m <= 10; m++)
			for (int n = 1; n <= 10; n++)
			{
				try
				{
					numb++;
					func(m, n, eps, method);
				}
				catch (runtime_error& e)
				{
					++fail_count;
					cerr << test_name << " fail: " << e.what() << endl;
				}
				catch (...) {
					++fail_count;
					cerr << "Unknown exception caught" << endl;
				}
			}
}

void Shell_SymRecCompress(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}
}

void Shell_DiagMult(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}

}

void Shell_RecMultL(ptr_test_mult_diag func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int k = 1; k <= 10; k++)
			{
				try
				{
					numb++;
					func(n, k, eps, method, smallsize);
				}
				catch (runtime_error& e)
				{
					++fail_count;
					cerr << test_name << " fail: " << e.what() << endl;
				}
				catch (...) {
					++fail_count;
					cerr << "Unknown exception caught" << endl;
				}
			}

}


void Shell_Add(ptr_test_add func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-4; eps > 1e-9; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int alpha = -10; alpha < 10; alpha += 2)
				for (int beta = -10; beta < 10; beta += 2)
				{
					if (alpha != 0 && beta != 0)
					{
						try
						{
							numb++;
							func(n, alpha, beta, eps, method, smallsize);
						}
						catch (runtime_error& e)
						{
							++fail_count;
							cerr << test_name << " fail: " << e.what() << endl;
						}
						catch (...) {
							++fail_count;
							cerr << "Unknown exception caught" << endl;
						}
					}
				}
}

void Shell_SymCompUpdate2(ptr_test_update func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-3; eps > 1e-9; eps /= 10)
		for (double alpha = -10; alpha < 10; alpha += 2)
			for (int n = 3; n <= 10; n++)
				for (int k = 1; k <= 10; k++)
				{
					try
					{
						numb++;
						func(n, k, alpha, eps, method, smallsize);
					}
					catch (runtime_error& e)
					{
						++fail_count;
						cerr << test_name << " fail: " << e.what() << endl;
					}
					catch (...) {
						++fail_count;
						cerr << "Unknown exception caught" << endl;
					}
				}
}

void Shell_SymCompRecInv(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}

}

void Shell_LowRankApproxTranspose(ptr_test_mult_diag func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int m = 2; m <= 10; m++)
			{
				try
				{
					numb++;
					func(m, n, eps, method, smallsize);
				}
				catch (runtime_error& e)
				{
					++fail_count;
					cerr << test_name << " fail: " << e.what() << endl;
				}
				catch (...) {
					++fail_count;
					cerr << "Unknown exception caught" << endl;
				}
			}
}

void Test_DirFactFastDiagStructOnline(size_m x, size_m y, size_m z, mnode** Gstr, double *B, double eps, int smallsize)
{
	printf("Testing factorization...\n");
	int n = x.n * y.n;
	int size = n * z.n;
	char bench[255] = "No";
	double *DD = alloc_arr(n * n); int lddd = n;
	double *DR = alloc_arr(n * n); int lddr = n;
	double norm = 0;

	double timer = 0;
	timer = omp_get_wtime();

	GenerateDiagonal2DBlock(0, x, y, z, DD, lddd);

	mnode *DCstr;
	SymCompRecInvStruct(n, Gstr[0], DCstr, smallsize, eps, "SVD");
	SymResRestoreStruct(n, DCstr, DR, lddr, smallsize);

	printf("Block %d. ", 0);
	norm = rel_error(n, n, DR, DD, lddd, eps);

	if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);

	free_arr(&DR);
	FreeNodes(n, DCstr, smallsize);

	for (int k = 1; k < z.n; k++)
	{
		double *DR = alloc_arr(n * n); int lddr = n;
		double *HR = alloc_arr(n * n); int ldhr = n;
		mnode *DCstr, *Hstr;

		printf("Block %d. ", k);

		SymCompRecInvStruct(n, Gstr[k], DCstr, smallsize, eps, "SVD");
		SymResRestoreStruct(n, DCstr, DR, lddr, smallsize);

		CopyStruct(n, Gstr[k - 1], Hstr, smallsize);
		DiagMultStruct(n, Hstr, &B[ind(k - 1, n)], smallsize);
		SymResRestoreStruct(n, Hstr, HR, ldhr, smallsize);

		// Norm of residual
#pragma omp parallel for schedule(static)
		for (int j = 0; j < n; j++)
#pragma omp simd
			for (int i = 0; i < n; i++)
				HR[i + ldhr * j] = HR[i + ldhr * j] + DR[i + lddr * j];

		norm = rel_error(n, n, HR, DD, lddd, eps);

		if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
		else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);

		FreeNodes(n, DCstr, smallsize);
		FreeNodes(n, Hstr, smallsize);
		free_arr(&DR);
		free_arr(&HR);
	}
	timer = omp_get_wtime() - timer;
	printf("Time: %lf\n", timer);

	free_arr(&DD);

}

void Test_DirSolveFactDiagStructConvergence(size_m x, size_m y, size_m z, mnode** Gstr, double thresh, int smallsize)
{
	printf("------------Test convergence-----------\n");
	int n = x.n * y.n;
	double norm = 0;
	for (int i = 1; i < z.n; i++)
	{
		double *GR1 = alloc_arr(n * n); int ldg1 = n;
		double *GR2 = alloc_arr(n * n); int ldg2 = n;
		double *GRL = alloc_arr(n * n); int ldgl = n;
		printf("For block %2d. \n", i);

		SymResRestoreStruct(n, Gstr[i], GR2, ldg2, smallsize);
		SymResRestoreStruct(n, Gstr[i - 1], GR1, ldg1, smallsize);
		norm = rel_error(n, n, GR2, GR1, ldg1, thresh);
		printf("norm |G[%2d] - G[%2d]|/|G[%2d]| = %12.10lf\n", i - 1, i, i - 1, norm);

		SymResRestoreStruct(n, Gstr[z.n - 1], GRL, ldgl, smallsize);
		norm = rel_error(n, n, GR1, GRL, ldgl, thresh);
		printf("norm |G[%2d] - G[%2d]|/|G[%2d]| = %12.10lf\n\n", i - 1, z.n - 1, z.n - 1, norm);

		free_arr(&GR1);
		free_arr(&GR2);
		free_arr(&GRL);
	}
}

void Test_DirSolveFactDiagStructBlockRanks(size_m x, size_m y, size_m z, mnode** Gstr)
{
	printf("----------Trees information-----------\n");
	int *size = (int*)malloc(x.n * sizeof(int));
	int *depth = (int*)malloc(y.n * sizeof(int));

	double time = omp_get_wtime();
	for (int i = 0; i < z.n; i++)
	{
#pragma omp parallel
		{
#pragma omp single
			{
				size[i] = TreeSize(Gstr[i]);
				depth[i] = MaxDepth(Gstr[i]);
			}
		}
	}
	double result = omp_get_wtime() - time;
	printf("Computational time of TreeSize and MaxDepth for all %d trees: %lf\n", x.n, result);

	for (int i = 0; i < z.n; i++)
	{
		printf("For block %2d. Size: %d, MaxDepth: %d, Ranks: ", i, size[i], depth[i]);
		PrintRanksInWidthList(Gstr[i]);
		printf("\n");
	}

	free(size);
	free(depth);

}

void Shell_CopyStruct(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}
}


void Test_SymRecCompress(int n, double eps, char *method, int smallsize)
{
	//printf("*****Test for SymRecCompress  n = %d eps = %e ******* ", n, eps);
	char frob = 'F';
	double norm = 0;

	double *H = alloc_arr(n * n); // init
	double *H1 = alloc_arr(n * n); // compressed
	double *H2 = alloc_arr(n * n); // recovered init

	int ldh = n;
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			H[i + ldh * j] = 1.0 / (i + j + 1);
			H1[i + ldh * j] = 1.0 / (i + j + 1);
		}
#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
#endif
	SymRecCompress(n, H1, ldh, smallsize, eps, "SVD");
	SymResRestore(n, H1, H2, ldh, smallsize);

#ifdef DEBUG
	print(n, n, H1, ldh, "H1 compressed");
	print(n, n, H2, ldh, "H recovered");
#endif

	// Norm of residual || A - L * U ||
	norm = rel_error(n, n, H2, H, ldh, eps);

#ifdef DEBUG
	print(n, n, H, ldh, "H init");
	print(n, n, H2, ldh, "diff");
#endif

	char str[255];
	sprintf(str, "Simple: n = %d ", n);
	AssertLess(norm, eps, str);

	free_arr(&H);
	free_arr(&H2);
	free_arr(&H1);
}

void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize)
{
	//printf("*****Test for SymRecCompressStruct  n = %d eps = %e ******* ", n, eps);
	char frob = 'F';
	double norm = 0;

	double *H = alloc_arr(n * n); // init
	double *H1 = alloc_arr(n * n); // compressed
	double *H2 = alloc_arr(n * n); // recovered init

	int ldh = n;
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			H[i + ldh * j] = 1.0 / (i + j + 1);
			H1[i + ldh * j] = 1.0 / (i + j + 1);
		}

#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
#endif

	mnode *H1str; // pointer to the tree head
	SymRecCompressStruct(n, H1, ldh, H1str, smallsize, eps, "SVD"); // recursive function means recursive allocation of memory for structure fields
	SymResRestoreStruct(n, H1str, H2, ldh, smallsize);

#ifdef DEBUG
	//print(n, n, H1, ldh, "H1 compressed");
	print(n, n, H2, ldh, "H recovered");
#endif

	// Norm of residual || A - L * U ||
	norm = rel_error(n, n, H2, H, ldh, eps);

#ifdef DEBUG
	print(n, n, H, ldh, "H init");
	print(n, n, H2, ldh, "diff");
#endif

	char str[255];
	sprintf(str, "Struct: n = %d ", n);
	AssertLess(norm, eps, str);

	FreeNodes(n, H1str, smallsize);
	free_arr(&H);
	free_arr(&H2);
	free_arr(&H1);
}

void Test_DiagMult(int n, double eps, char *method, int smallsize)
{
	//printf("*****Test for DiagMult  n = %d \n******* ", n);
	double *Hd = alloc_arr(n * n); // diagonal Hd = D * H * D
	double *H1 = alloc_arr(n * n); // compressed H
	double *H2 = alloc_arr(n * n); // recovered H after D * H1 * D
	double *d = alloc_arr(n);
	char str[255];

	double norm = 0;
	int ldh = n;

	for (int j = 0; j < n; j++)
	{
		d[j] = j + 1;
	}

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			Hd[i + ldh * j] = 1.0 / (i + j + 1);
			Hd[i + ldh * j] *= d[j];
			Hd[i + ldh * j] *= d[i];
			H1[i + ldh * j] = 1.0 / (i + j + 1);
		}
#ifdef DEBUG
	print(n, n, H1, ldh, "Initial H");
#endif
	// Compress H1
	SymRecCompress(n, H1, ldh, smallsize, eps, method);

	// Compressed H1 = D * H * D
	DiagMult(n, H1, ldh, d, smallsize);

	// Recove H1 to uncompressed form
	SymResRestore(n, H1, H2, ldh, smallsize);

#ifdef DEBUG
	print(n, n, Hd, ldh, "Initial Hd = D * H * D");
	print(n, n, H2, ldh, "Recovered H2 = (D * H * D)comp");
#endif

	// Compare Hd and H2
	norm = rel_error(n, n, H2, Hd, ldh, eps);

	sprintf(str, "Simple: n = %d ", n);
	AssertLess(norm, eps, str);

	free_arr(&Hd); // diagonal Hd = D * H * D
	free_arr(&H1); // compressed H
	free_arr(&H2); // recovered H after D * H1 * D
	free_arr(&d);
}

void Test_DiagMultStruct(int n, double eps, char *method, int smallsize)
{
	//printf("*****Test for DiagMultStruct  n = %d ******* ", n);
	double *Hd = alloc_arr(n * n); // diagonal Hd = D * H * D
	double *H1 = alloc_arr(n * n); // compressed H
	double *H2 = alloc_arr(n * n); // recovered H after D * H1 * D
	double *d = alloc_arr(n);
	char str[255];

	double norm = 0;
	int ldh = n;

	for (int j = 0; j < n; j++)
	{
		d[j] = j + 1;
	}

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			Hd[i + ldh * j] = 1.0 / (i + j + 1);
			Hd[i + ldh * j] *= d[j];
			Hd[i + ldh * j] *= d[i];
			H1[i + ldh * j] = 1.0 / (i + j + 1);
		}
#ifdef DEBUG
	print(n, n, H1, ldh, "Initial H");
#endif

	mnode *HCstr;
	// Compress H1 to structured form
	SymRecCompressStruct(n, H1, ldh, HCstr, smallsize, eps, method);

	// Compressed H1 = D * H * D
	DiagMultStruct(n, HCstr, d, smallsize);

	// Recove H1 to uncompressed form
	SymResRestoreStruct(n, HCstr, H2, ldh, smallsize);

#ifdef DEBUG
	print(n, n, Hd, ldh, "Initial Hd = D * H * D");
	print(n, n, H2, ldh, "Recovered H2 = (D * H * D)comp");
#endif

	// Compare Hd and H2
	norm = rel_error(n, n, H2, Hd, ldh, eps);

	sprintf(str, "Struct: n = %d ", n);
	AssertLess(norm, eps, str);

	FreeNodes(n, HCstr, smallsize);
	free_arr(&Hd); // diagonal Hd = D * H * D
	free_arr(&H1); // compressed H
	free_arr(&H2); // recovered H after D * H1 * D
	free_arr(&d);
}

/* Тест на сравнение результатов умножения Y = H * X сжимаемой матрицы H на произвольную X.
Сравниваются результаты со сжатием и без */
void Test_RecMultL(int n, int k, double eps, char *method, int smallsize)
{
	//printf("*****Test for RecMultL  n = %d k = %d ******* ", n, k);
	double *H = alloc_arr(n * n); // init and compressed
	double *X = alloc_arr(n * k);
	double *Y = alloc_arr(n * k); // init Y
	double *Y1 = alloc_arr(n * k); // after multiplication woth compressed
	char str[255];

	double norm = 0;
	double alpha = 1.0;
	double beta = 0.0;

	int ldh = n;
	int ldy = n;
	int ldx = n;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			H[i + ldh * j] = 1.0 / (i + j + 1);
		}
		for (int j = 0; j < k; j++)
		{
			X[i + ldx * j] = 1.0 / (i + j + 1);
		}
	}

	dgemm("No", "No", &n, &k, &n, &alpha, H, &ldh, X, &ldx, &beta, Y, &ldy);

#ifdef DEBUG
	print(n, n, H, ldy, "H init");
	print(n, k, X, ldy, "X init");
	print(n, k, Y, ldy, "Y init");
#endif
	// Compress H
	SymRecCompress(n, H, ldh, smallsize, eps, method);

	// RecMult Y1 = comp(H) * X
	RecMultL(n, k, H, ldh, X, ldx, Y1, ldy, smallsize);

	norm = rel_error(n, k, Y1, Y, ldy, eps);

	sprintf(str, "Simple: n = %d k = %d ", n, k);
	AssertLess(norm, eps, str);

#ifdef DEBUG
	print(n, n, H, ldy, "H comp");
	print(n, k, Y1, ldy, "Y1 rec");
#endif

	free_arr(&H);
	free_arr(&X);
	free_arr(&Y);
	free_arr(&Y1);
}

/* Тест на сравнение результатов умножения Y = H * X сжимаемой матрицы H на произвольную X.
Сравниваются результаты со сжатием и без */
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize)
{
	//printf("*****Test for RecMultLStruct  n = %d k = %d ******* ", n, k);
	double *H = alloc_arr(n * n); // init and compressed
	double *X = alloc_arr(n * k);
	double *Y = alloc_arr(n * k); // init Y
	double *Y1 = alloc_arr(n * k); // after multiplication woth compressed
	char str[255];

	double norm = 0;
	double alpha = 1.0;
	double beta = 0.0;

	int ldh = n;
	int ldy = n;
	int ldx = n;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			H[i + ldh * j] = 1.0 / (i + j + 1);
		}
		for (int j = 0; j < k; j++)
		{
			X[i + ldx * j] = 1.0 / (i + j + 1);
		}
	}

	dgemm("No", "No", &n, &k, &n, &alpha, H, &ldh, X, &ldx, &beta, Y, &ldy);

#ifdef DEBUG
	print(n, n, H, ldy, "H init");
	print(n, k, X, ldy, "X init");
	print(n, k, Y, ldy, "Y init");
#endif

	mnode *Hstr;
	// Compress H
	SymRecCompressStruct(n, H, ldh, Hstr, smallsize, eps, method);

	// RecMult Y1 = comp(H) * X
	RecMultLStruct(n, k, Hstr, X, ldx, Y1, ldy, smallsize);

	norm = rel_error(n, k, Y1, Y, ldy, eps);
	sprintf(str, "Struct: n = %d k = %d ", n, k);
	AssertLess(norm, eps, str);

#ifdef DEBUG
	print(n, n, H, ldy, "H comp");
	print(n, k, Y1, ldy, "Y1 rec");
#endif

	FreeNodes(n, Hstr, smallsize);
	free_arr(&H);
	free_arr(&X);
	free_arr(&Y);
	free_arr(&Y1);
}

void Test_LowRankApprox(int m, int n, double eps, char *method)
{
	//printf("Test for LowRankApproximation m = %d n = %d ", m, n);
	// A - matrix in dense order
	double *A = alloc_arr(m * n);
	double *A_init = alloc_arr(m * n);
	double *A_rec = alloc_arr(m * n);
	double *VT;
	char str[255];

	int mn = min(m, n);
	VT = alloc_arr(mn * n);

	int lda = m;
	int ldvt = mn;

	double norm = 0;
	double alpha = 1.0;
	double beta = 0.0;

	int p = 0; // p <= mn

	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			A[i + lda * j] = 1.0 / (n + i + j + 1);
			A_init[i + lda * j] = 1.0 / (n + i + j + 1);
		}

	LowRankApprox(m, n, A, lda, VT, ldvt, p, eps, "SVD");

	dgemm("no", "no", &m, &n, &p, &alpha, A, &lda, VT, &ldvt, &beta, A_rec, &lda);

	norm = rel_error(m, n, A_rec, A_init, lda, eps);
	sprintf(str, "Simple: n = %d m = %d ", n, m);
	AssertLess(norm, eps, str);

	free_arr(&A);
	free_arr(&A_init);
	free_arr(&A_rec);
	free_arr(&VT);

}

void Test_LowRankApproxStruct(int m, int n, double eps, char *method)
{
	// A - matrix in dense order
	double *A = alloc_arr(m * n);
	double *A_init = alloc_arr(m * n);
	double *A_rec = alloc_arr(m * n);
	char str[255];

	int lda = m;

	double norm = 0;
	double alpha = 1.0;
	double beta = 0.0;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			A[i + lda * j] = 1.0 / (n + i + j + 1);
			A_init[i + lda * j] = 1.0 / (n + i + j + 1);
		}

#if 1
	mnode *Astr = (mnode*)malloc(sizeof(mnode));
	//printf("Test for LowRankApproximationStruct m = %d n = %d ", m, n);
	LowRankApproxStruct(m, n, A, lda, Astr, eps, "SVD"); // memory allocation for Astr inside function
#else
	mnode *Astr;
	//printf("Test for LowRankApproximationStruct2 using return m = %d n = %d ", m, n);
	Astr = LowRankApproxStruct2(m, n, A, lda, eps, "SVD"); // memory allocation for Astr inside function
#endif
//	printf("p = %d ", Astr->p);

	dgemm("no", "no", &m, &n, &Astr->p, &alpha, Astr->U, &m, Astr->VT, &Astr->p, &beta, A_rec, &lda);

	norm = rel_error(m, n, A_rec, A_init, lda, eps);
	sprintf(str, "Struct: n = %d m = %d ", n, m);
	AssertLess(norm, eps, str);
	

	free_arr(&A);
	free_arr(&A_init);
	free_arr(&A_rec);
	free_arr(&Astr->U);
	free_arr(&Astr->VT);
	free(Astr);
}

void Test_Add(int n, double alpha, double beta, double eps, char *method, int smallsize)
{
	//printf("*****Test for Add n = %d ******* ", n);
	double *H1 = alloc_arr(n * n);
	double *H2 = alloc_arr(n * n);
	double *G = alloc_arr(n * n);
	double *H1c = alloc_arr(n * n);
	double *H2c = alloc_arr(n * n);
	double *Gc = alloc_arr(n * n);
	double *GcR = alloc_arr(n * n);
	char str[255];

	int ldh = n;
	int ldg = n;
	double norm = 0;

#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			H1[i + ldh * j] = 1.0 / (i + j + 1);
			H2[i + ldh * j] = 1.0 / (i*i + j*j + 1);
			H1c[i + ldh * j] = 1.0 / (i + j + 1);
			H2c[i + ldh * j] = 1.0 / (i*i + j*j + 1);
		}
#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
	print(n, n, H2, ldh, "H2");
#endif

	SymRecCompress(n, H1c, ldh, smallsize, eps, method);
	SymRecCompress(n, H2c, ldh, smallsize, eps, method);

#ifdef DEBUG
	print(n, n, H1c, ldh, "H1c");
	print(n, n, H2c, ldh, "H2c");
#endif

	Add_dense(n, n, alpha, H1, ldh, beta, H2, ldh, G, ldg);
	Add(n, alpha, H1c, ldh, beta, H2c, ldh, Gc, ldg, smallsize, eps, method);

#ifdef DEBUG
	print(n, n, G, ldg, "res_dense");
	print(n, n, Gc, ldg, "res_comp");
#endif

	SymResRestore(n, Gc, GcR, ldg, smallsize);

#ifdef DEBUG
	print(n, n, GcR, ldg, "res_comp_restore");
#endif
	// |GcR - G| / |G|
	norm = rel_error(n, n, GcR, G, ldg, eps);

	sprintf(str, "Simple: n = %d n = %d alpha = %d beta = %d", n, n, alpha, beta);
	AssertLess(norm, eps, str);

	free_arr(&H1);
	free_arr(&H2);
	free_arr(&G);
	free_arr(&H1c);
	free_arr(&H2c);
	free_arr(&Gc);
	free_arr(&GcR);
}

void Test_AddStruct(int n, double alpha, double beta, double eps, char *method, int smallsize)
{
	//printf("*****Test for Add n = %d ******* ", n);
	double *H1 = alloc_arr(n * n);
	double *H2 = alloc_arr(n * n);
	double *G = alloc_arr(n * n);
	double *H1c = alloc_arr(n * n);
	double *H2c = alloc_arr(n * n);
	double *Gc = alloc_arr(n * n);
	double *GcR = alloc_arr(n * n);
	char str[255];

	int ldh = n;
	int ldg = n;
	double norm = 0;

#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			H1[i + ldh * j] = 1.0 / (i + j + 1);
			H2[i + ldh * j] = 1.0 / (i*i + j*j + 1);
			H1c[i + ldh * j] = 1.0 / (i + j + 1);
			H2c[i + ldh * j] = 1.0 / (i*i + j*j + 1);
		}

#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
	print(n, n, H2, ldh, "H2");
#endif

	mnode *H1str, *H2str;
	SymRecCompressStruct(n, H1c, ldh, H1str, smallsize, eps, method);
	SymRecCompressStruct(n, H2c, ldh, H2str, smallsize, eps, method);

#ifdef DEBUG
	print(n, n, H1c, ldh, "H1c");
	print(n, n, H2c, ldh, "H2c");
#endif

	mnode *Gstr;
	Add_dense(n, n, alpha, H1, ldh, beta, H2, ldh, G, ldg);
	AddStruct(n, alpha, H1str, beta, H2str, Gstr, smallsize, eps, method);

#ifdef DEBUG
	print(n, n, G, ldg, "res_dense");
	print(n, n, Gc, ldg, "res_comp");
#endif

	SymResRestoreStruct(n, Gstr, GcR, ldg, smallsize);

#ifdef DEBUG
	print(n, n, GcR, ldg, "res_comp_restore");
#endif
	// |GcR - G| / |G|
	norm = rel_error(n, n, GcR, G, ldg, eps);
	sprintf(str, "Struct: n = %d n = %d alpha = %lf", n, n, alpha, beta);
	AssertLess(norm, eps, str);

	FreeNodes(n, H1str, smallsize);
	FreeNodes(n, H2str, smallsize);
	FreeNodes(n, Gstr, smallsize);
	free_arr(&H1);
	free_arr(&H2);
	free_arr(&G);
	free_arr(&H1c);
	free_arr(&H2c);
	free_arr(&Gc);
	free_arr(&GcR);
}

// B = H - V * Y * VT
void Test_SymCompUpdate2(int n, int k, double alpha, double eps, char* method, int smallsize)
{
//	printf("*****Test for SymCompUpdate2   n = %d k = %d ***** ", n, k);
	double *B = alloc_arr(n * n); int ldb = n;
	double *B_rec = alloc_arr(n * n);
	double *Y = alloc_arr(k * k); int ldy = k;
	double *V = alloc_arr(n * k); int ldv = n; int ldvtr = k;
	double *HC = alloc_arr(n * n); int ldh = n;
	double *H = alloc_arr(n * n);
	double *C = alloc_arr(n * k); int ldc = n;
	char str[255];

	double alpha_one = 1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;
	double norm = 0;

	Hilbert(n, HC, ldh);
	Hilbert(n, H, ldh);

#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < k; i++)
		Y[i + ldy * i] = i + 1;

#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < k; j++)
		for (int i = 0; i < n; i++)
			V[i + ldv * j] = (i + j + 1);

	// C = V * Y
	dsymm("Right", "Up", &n, &k, &alpha_one, Y, &ldy, V, &ldv, &beta_zero, C, &ldc);

	// H = H + alpha * C * VT
	dgemm("No", "Trans", &n, &n, &k, &alpha, C, &ldc, V, &ldv, &beta_one, H, &ldh);

	// Compressed update
	SymRecCompress(n, HC, ldh, smallsize, eps, method);

	SymCompUpdate2(n, k, HC, ldh, alpha, Y, ldy, V, ldv, B, ldb, smallsize, eps, method);
	SymResRestore(n, B, B_rec, ldh, smallsize);

#ifdef DEBUG
	print(n, n, B_rec, ldb, "B_rec");
	print(n, n, H, ldh, "H");
#endif
	// || B_rec - H || / || H ||
	norm = rel_error(n, n, B_rec, H, ldh, eps);
	sprintf(str, "Simple: n = %d k = %d alpha = %lf", n, k, alpha);
	AssertLess(norm, eps, str);

	free_arr(&B);
	free_arr(&B_rec);
	free_arr(&H);
	free_arr(&HC);
	free_arr(&Y);
	free_arr(&C);
	free_arr(&V);
}

// B = H - V * Y * VT
void Test_SymCompUpdate2Struct(int n, int k, double alpha, double eps, char* method, int smallsize)
{
	//printf("*****Test for SymCompUpdate2Struct  n = %d k = %d ***** ", n, k);
	double *B = alloc_arr(n * n); int ldb = n;
	double *B_rec = alloc_arr(n * n);
	double *Y = alloc_arr(k * k); int ldy = k;
	double *V = alloc_arr(n * k); int ldv = n; int ldvtr = k;
	double *HC = alloc_arr(n * n); int ldh = n;
	double *H = alloc_arr(n * n);
	double *C = alloc_arr(n * k); int ldc = n;
	char str[255];

	double alpha_one = 1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;
	double norm = 0;


	Hilbert(n, HC, ldh);
	Hilbert(n, H, ldh);

#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < k; i++)
		Y[i + ldy * i] = i + 1;

#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < k; j++)
		for (int i = 0; i < n; i++)
			V[i + ldv * j] = (i + j + 1);

	// C = V * Y
	dsymm("Right", "Up", &n, &k, &alpha_one, Y, &ldy, V, &ldv, &beta_zero, C, &ldc);

	// H = H + alpha * C * VT
	dgemm("No", "Trans", &n, &n, &k, &alpha, C, &ldc, V, &ldv, &beta_one, H, &ldh);

	mnode *HCstr;
	// Compressed update
	SymRecCompressStruct(n, HC, ldh, HCstr, smallsize, eps, method);

	mnode *Bstr;
	SymCompUpdate2Struct(n, k, HCstr, alpha, Y, ldy, V, ldv, Bstr, smallsize, eps, method);
	SymResRestoreStruct(n, Bstr, B_rec, ldh, smallsize);

#ifdef DEBUG
	print(n, n, B_rec, ldb, "B_rec");
	print(n, n, H, ldh, "H");
#endif

	// || B_rec - H || / || H ||
	norm = rel_error(n, n, B_rec, H, ldh, eps);
	sprintf(str, "Struct: n = %d k = %d alpha = %lf", n, k, alpha);
	AssertLess(norm, eps, str);

	FreeNodes(n, Bstr, smallsize);
	FreeNodes(n, HCstr, smallsize);
	free_arr(&B);
	free_arr(&B_rec);
	free_arr(&H);
	free_arr(&HC);
	free_arr(&Y);
	free_arr(&C);
	free_arr(&V);
}

void Test_SymCompRecInv(int n, double eps, char *method, int smallsize)
{
	//printf("***** Test_SymCompRecInv n = %d eps = %lf **** ", n, eps);
	double *H = alloc_arr(n * n);
	double *Hc = alloc_arr(n * n);
	double *Bc = alloc_arr(n * n);
	double *Brec = alloc_arr(n * n);
	double *Y = alloc_arr(n * n);
	char str[255];

	int ldh = n;
	int ldb = n;
	int ldy = n;

	double alpha_mone = -1.0;
	double beta_one = 1.0;
	double norm = 0;

	Hilbert(n, H, ldh);
	Hilbert(n, Hc, ldh);

	// for stability
	for (int i = 0; i < n; i++)
	{
		H[i + ldh * i] += 1.0;
		Hc[i + ldh * i] += 1.0;
	}

	SymRecCompress(n, Hc, ldh, smallsize, eps, method);
	SymCompRecInv(n, Hc, ldh, Bc, ldb, smallsize, eps, method);
	SymResRestore(n, Bc, Brec, ldb, smallsize);

	Eye(n, Y, ldy);

	// Y = Y - H * Brec
	dgemm("No", "No", &n, &n, &n, &alpha_mone, H, &ldh, Brec, &ldb, &beta_one, Y, &ldy);

	norm = dlange("Frob", &n, &n, Y, &ldy, NULL);
	sprintf(str, "Simple: n = %d", n);
	AssertLess(norm, eps, str);

	//if (norm < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	//else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);
}

void Test_SymCompRecInvStruct(int n, double eps, char *method, int smallsize)
{
	//printf("***** Test_SymCompRecInvStruct n = %d eps = %lf ****", n, eps);
	double *H = alloc_arr(n * n);
	double *Hc = alloc_arr(n * n);
	double *Bc = alloc_arr(n * n);
	double *Brec = alloc_arr(n * n);
	double *Y = alloc_arr(n * n);
	char str[255];

	int ldh = n;
	int ldb = n;
	int ldy = n;

	double alpha_mone = -1.0;
	double beta_one = 1.0;
	double norm = 0;

	Hilbert(n, H, ldh);
	Hilbert(n, Hc, ldh);

	// for stability
	for (int i = 0; i < n; i++)
	{
		H[i + ldh * i] += 1.0;
		Hc[i + ldh * i] += 1.0;
	}

	mnode *HCstr, *BCstr;
	SymRecCompressStruct(n, Hc, ldh, HCstr, smallsize, eps, method);
	SymCompRecInvStruct(n, HCstr, BCstr, smallsize, eps, method);
	SymResRestoreStruct(n, BCstr, Brec, ldb, smallsize);

	Eye(n, Y, ldy);

	// Y = Y - H * Brec
	dgemm("No", "No", &n, &n, &n, &alpha_mone, H, &ldh, Brec, &ldb, &beta_one, Y, &ldy);

	norm = dlange("Frob", &n, &n, Y, &ldy, NULL);
	sprintf(str, "Struct: n = %d", n);
	AssertLess(norm, eps, str);

	//if (norm < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	//else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);

	FreeNodes(n, HCstr, smallsize);
	FreeNodes(n, BCstr, smallsize);
	free_arr(&H);
	free_arr(&Hc);
	free_arr(&Bc);
	free_arr(&Brec);
	free_arr(&Y);
}

void Test_LowRankApproxTranspose(int m, int n, double eps, char *method, int smallsize)
{
	double *H = alloc_arr(m * n); int ldh = m;
	double *Hinit = alloc_arr(m * n);
	double *Hrec = alloc_arr(m * n);
	double *Htr = alloc_arr(n * m); int ldhtr = n;

	int mn = min(m, n);
	double *V = alloc_arr(mn * n); int ldv = mn;
	double *Vtr = alloc_arr(mn * m); int ldvtr = mn;
	char str[255];

	double alpha = 1.0;
	double beta = 0.0;
	double norm = 0;

	int p1 = 0, p2 = 0;
	int k = 0;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			k++;
			H[i + ldh * j] = 1.0 / (i + j + 1 + k);
			Hinit[i + ldh * j] = 1.0 / (i + j + 1 + k);
		}

	Mat_Trans(m, n, H, ldh, Htr, ldhtr);

	LowRankApprox(n, m, Htr, ldhtr, Vtr, ldvtr, p1, eps, "SVD");
	LowRankApprox(m, n, H, ldh, V, ldv, p2, eps, "SVD");

	dgemm("Trans", "Trans", &m, &n, &mn, &alpha, Vtr, &ldvtr, Htr, &ldhtr, &beta, Hrec, &ldh);

	// || H _rec  - H || / || H ||
	norm = rel_error(m, n, Hrec, Hinit, ldh, eps);

	sprintf(str, "Simple: m = %d n = %d", m, n);
	AssertLess(norm, eps, str);
	AssertEqual(p1, p2, str);

	free_arr(&H);
	free_arr(&Hinit);
	free_arr(&Htr);
	free_arr(&Hrec);
	free_arr(&V);
	free_arr(&Vtr);
}

void Test_CopyStruct(int n, double eps, char *method, int smallsize)
{
	double *H = alloc_arr(n * n);
	double *H1 = alloc_arr(n * n);
	double *H2 = alloc_arr(n * n);
	char str[255];

	double norm = 0;
	int ldh = n;

	//printf("***Test CopyStruct n = %d ", n);

	Hilbert(n, H, ldh);
	Hilbert(n, H1, ldh);

	mnode* Hstr, *Hcopy_str;
	SymRecCompressStruct(n, H, ldh, Hstr, smallsize, eps, method);
	CopyStruct(n, Hstr, Hcopy_str, smallsize);
	SymResRestoreStruct(n, Hcopy_str, H2, ldh, smallsize);

	norm = rel_error(n, n, H2, H1, ldh, eps);
	sprintf(str, "Struct: n = %d", n);
	AssertLess(norm, eps, str);

	FreeNodes(n, Hstr, smallsize);
	FreeNodes(n, Hcopy_str, smallsize);
	free_arr(&H2);
	free_arr(&H1);
	free_arr(&H);
}

void Test_RankEqual(mnode *Astr, mnode *Bstr)
{
	try
	{
		char str[255] = "Rank(A01) = Rank(B01) ";
		AssertEqual(Astr->p, Bstr->p, str);
	}
	catch (exception &ex)
	{
		cout << ex.what();
	}

	if (Astr->left != NULL || Bstr->left != NULL)
	{
		Test_RankEqual(Astr->left, Bstr->left);
	}

	if (Astr->right != NULL || Bstr->right != NULL)
	{
		Test_RankEqual(Astr->right, Bstr->right);
	}
}

void Test_RankAdd(mnode *Astr, mnode *Bstr, mnode *Cstr)
{
	try
	{
		char str[255] = "Rank(C01) <= Rank(A01) + Rank(B01) ";
		AssertLess(Astr->p + Bstr->p, Cstr->p, str);
	}
	catch (exception &ex)
	{
		cout << ex.what();
	}

}

void Test_Dense_to_CSR(size_m x, size_m y, size_m z, int non_zeros_in_3diag, double *D, int ldd)
{
	int n = x.n * y.n;
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

	/*
	print(n, n, &B_mat[0], ldb, "B[0]");
	print(size, n, D, ldd, "D");

	printf("all non_zero elements: %d\n", non_zeros_in_block3diag);
	for (int i = 0; i < size + 1; i ++)
	printf("%d: ia = %d value(ia) = %lf diff = %d\n", i, Dcsr->ia[i], Dcsr->values[Dcsr->ia[i] - 1], Dcsr->ia[i+1]- Dcsr->ia[i]);

	print_vec(non_zeros_in_block3diag, Dcsr->ja, Dcsr->values, "ja and values");
	print_map(CSR);*/

}

void Test_QueueList(int n, double eps, char* method, int smallsize)
{
	printf("----Test for queue implementation with LIST---\n n = %d \n", n);
	// A - matrix in dense order
	double *A = alloc_arr(n * n);
	double *A_init = alloc_arr(n * n);
	int lda = n;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			A[i + lda * j] = 1.0 / (i + j + 1);
			A_init[i + lda * j] = 1.0 / (i + j + 1);
		}

	mnode *ACstr;
	SymRecCompressStruct(n, A, lda, ACstr, smallsize, eps, method);

	//printf("first level. p = %d\n", ACstr->p);
	//printf("second level. left: %d right: %d\n", Astr->left->p, Astr->right->p);

	printf("Size: %d, MaxDepth: %d, Ranks: ", TreeSize(ACstr), MaxDepth(ACstr));
	//	PrintRanksInWidth(ACstr);
	printf("List:\n");
	PrintRanksInWidthList(ACstr);

	FreeNodes(n, ACstr, smallsize);
	free_arr(&A);
	free_arr(&A_init);
}

void Test_TransferBlock3Diag_to_CSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_orig, double *f, double eps)
{
	int n = n1 * n2;
	int size = n * n3;
	double RelRes = 0;
	double *g = alloc_arr(size);
	ResidCSR(n1, n2, n3, Dcsr, x_orig, f, g, RelRes);

	if (RelRes < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", RelRes, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", RelRes, eps);

	free_arr(&g);
}

void Test_CompareColumnsOfMatrix(int n1, int n2, int n3, double* D, int ldd, double* B, dcsr* Dcsr, double thresh)
{
	int n = n1 * n2;
	int size = n * n3;
	double RelRes = 0;
	double *g = alloc_arr(size);
	double *f1 = alloc_arr(size);
	double *f2 = alloc_arr(size);
	double Res1 = 0;
	double Res2 = 0;

	for (int i = 0; i < size; i++)
	{
		double *x_test = alloc_arr(size);
		x_test[i] = 1;
		Mult_Au(n1, n2, n3, D, ldd, B, x_test, f1);
		mkl_dcsrgemv("No", &size, Dcsr->values, Dcsr->ia, Dcsr->ja, x_test, f2);
		//	print_vec(size, f1, f2, "f1 and f2");
		rel_error(size, 1, f1, f2, size, thresh);
		free(x_test);
	}

	free_arr(&f1);
	free_arr(&f2);
	free_arr(&g);
}



