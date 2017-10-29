#include "Header.h"


void print(int m, int n, double *u, int ldu, char *mess)
{
	printf("%s\n", mess);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%6.3lf ", u[i + ldu*j]);
			//if (j % N == N - 1) printf("|");
		}
		printf("\n");
		//if (i % N == N - 1) printf("\n");
	}

	printf("\n");
	
}

// Функция выделения памяти под массив
double* alloc_arr(int n)
{
	double *f = (double*)malloc(n * sizeof(double));

#pragma omp parallel for
	for (int i = 0; i < n; i++)
		f[i] = 0.0;

	return f;
}

// Функция освобождения памяти 
void free_arr(double* * arr)
{
	free(*arr);
}

// Инициализация значений массивов. Заданий правой части и краевых условий
void Init_matrix(int nbl, int n2, int n3, double *A, int ldA)
{
	// n2 = nbl * nbl - size of small matrix
	// [n*n] x [3*n] of big matrix A

	int size = n2 * n3;
	double h = 1.0 / (N);
	int lda = nbl;

	A[0:size] = 0;

	for (int j = 0; j < nbl; j++)
	{
		if (j == 0)
		{
			fill_middle(nbl, j, h, &A[j * nbl + ldA * 0], ldA);
			fill_right(nbl, j, h, &A[j * nbl + ldA * nbl], ldA);
		}
		else if (j == nbl - 1)
		{
			fill_left(nbl, j, h, &A[j * nbl + ldA * nbl], ldA);
			fill_middle(nbl, j, h, &A[j * nbl + ldA * 2 * nbl], ldA);
		}
		else
		{
			fill_left(nbl, j, h, &A[j * nbl + ldA * 0], ldA);
			fill_middle(nbl, j, h, &A[j * nbl + ldA * nbl], ldA);
			fill_right(nbl, j, h, &A[j * nbl + ldA * 2 * nbl], ldA);
		}
	}

}

void fill_left(int nbl, int j, double h, double *A, int ldA)
{
	for (int i = 0; i < nbl; i++)
		A[i + ldA*i] = 1.0 / h;
}

void fill_right(int nbl, int j, double h, double *A, int ldA)
{
	for (int i = 0; i < nbl; i++)
		A[i + ldA*i] = 1.0 / h;
}

void fill_middle(int nbl, int j, double h, double *A, int ldA)
{

	A[0 + ldA * 0] = -1.0 / h;
	A[0 + ldA * 1] = 1.0 / h;

	for (int i = 1; i < nbl - 1; i++)
	{
		A[i + ldA*(i - 1)] = 1.0 / h;
		A[i + ldA*(i)] = -4.0 / h;
		A[i + ldA*(i + 1)] = 1.0 / h;
	}

	A[nbl - 1 + ldA * (nbl - 2)] = 1.0 / h;
	A[nbl - 1 + ldA * (nbl - 1)] = -1.0 / h;

}

void sparse_lu(int nbl, int n2, int n3, double *A, int ldA, int **ipiv_mat)
{
	int *ipiv = new int[nbl];
	int info;
	char low = 'L';
	char up = 'U';
	char unit = 'U';
	char no = 'N';
	char all = 'A';
	char left = 'L';
	char right = 'R';
	int ldL = nbl;
	int ldU = nbl;
	int start, start_l, start_r;
	double zero = 0.0;
	double one = 1.0;
	double mone = -1.0;
	int ione = 1;

	double *U = alloc_arr(nbl * nbl);
	double *L = alloc_arr(nbl * nbl);

	U[0:nbl*nbl] = 0;
	L[0:nbl*nbl] = 0;
	int mione = -1;

	for (int j = 0; j < nbl; j++)
	{
		if (j == 0) start = 0 + ldA * 0;
		else if (j == nbl - 1) start = n2 - nbl + ldA * 2 * nbl;
		else start = j * nbl + ldA * nbl;

		// LU of the A11 and copy factors

		dgetrf(&nbl, &nbl, &A[start], &ldA, ipiv, &info);

		// Copy factors
		dlacpy(&low, &nbl, &nbl, &A[start], &ldA, L, &ldL);
		dlacpy(&up, &nbl, &nbl, &A[start], &ldA, U, &ldU);
	
		for (int i = 0; i < nbl; i++)
			if (ipiv[i] != i + 1) printf("Iter: %d, row %d interchanged with row %d\n", j + 1, ipiv[i], i + 1);
	

		for (int i = 0; i < nbl; i++)
		{
			ipiv_mat[j][i] = ipiv[i];
		}
		
		if (j == nbl - 1) break;
		else
		{
			// Apply A21 * U^(-1)

			if (j == nbl - 2) start_l = n2 - nbl + ldA * nbl;
			else start_l = (j + 1) * nbl + ldA * 0;

			// Inversion U

			dtrtri(&up, &no, &nbl, U, &ldU, &info);

			// Apply A21 * U^(-1)

			dtrmm(&right, &up, &no, &no, &nbl, &nbl, &one, U, &ldU, &A[start_l], &ldA);

			// Apply L^(-1) * P ^ (-1) * A12

			if (j == 0) start_r = 0 + ldA * nbl;
			else start_r = (j)* nbl + ldA * 2 * nbl;

			// Inversion of L
			
			dtrtri(&low, &unit, &nbl, L, &ldL, &info);

			// Swap rows of A12 due to P^(-1) * A12

			dlaswp(&nbl, &A[start_r], &ldA, &ione, &nbl, ipiv, &ione);

			// Apply L^(-1) to P^(-1) * A12

			dtrmm(&left, &low, &no, &unit, &nbl, &nbl, &one, L, &ldL, &A[start_r], &ldA);

			// Compute Schur component

			dgemm(&no, &no, &nbl, &nbl, &nbl, &mone, &A[start_l], &ldA, &A[start_r], &ldA, &one, &A[start_l + ldA * nbl], &ldA);
		}

	}

	free_arr(&U);
	free_arr(&L);
}

void change(char right, int nbl, double *L, int ldL, int *ipiv)
{
	if (right == 'R' || right == 'r') // switch columns
	{
		for (int i = 0; i < nbl; i++)
			if (ipiv[i] != (i + 1)) swap_col(nbl, i, ipiv[i] - 1, L, ldL, ipiv);
	}

}

void swap_col(int nbl, int k, int m, double *L, int ldL, int *ipiv)
{
	double *c = alloc_arr(nbl);
	for (int i = 0; i < nbl; i++)
	{
		c[i] = L[i + ldL * m];
		L[i + ldL * m] = L[i + ldL * k];
		L[i + ldL * k] = c[i];
	}
	free_arr(&c);
}

void test(int nbl, int n2, int n3, double *A_f, int ldaf, double *U_f, int lduf, double *L_f, int ldlf, int **ipiv_mat)
{
	char low = 'L';
	char up = 'U';
	char unit = 'U';
	char no = 'N';
	char all = 'A';
	char left = 'L';
	char right = 'R';
	int start, start_l, start_r;
	double zero = 0.0;
	double one = 1.0;
	double mone = -1.0;
	int ione = 1;
	int mione = -1;
	double norm = 0;
	char frob = 'F';

	double *A_res = alloc_arr(n2*n2);
	A_res[0:n2*n2] = 0;

	// Switch rows L = P * L
	for (int j = 0; j < nbl; j++)
	{
		dlaswp(&nbl, &L_f[j * nbl + ldlf * (j) * nbl], &ldlf, &ione, &nbl, ipiv_mat[j], &mione);
	}

	// L = L * U , L is overwritten
	//dtrmm(&right, &up, &no, &no, &n2, &n2, &one, U_f, &lduf, L_f, &ldlf);

	// A = L * U , A - new matrix
	dgemm(&no, &no, &n2, &n2, &n2, &one, L_f, &ldlf, U_f, &lduf, &zero, A_res, &ldaf);

#ifdef DEBUG
	print(n2, n2, A_f, ldaf,"matrix A result: A = L * U" );
#endif

	// Norm of residual || A - L * U ||
	for (int j = 0; j < n2; j++)
		for (int i = 0; i < n2; i++)
			A_res[i + ldaf * j] = A_res[i + ldaf * j] - A_f[i + ldaf * j];

	norm = dlange(&frob, &n2, &n2, A_res, &ldaf, NULL);
	if (norm < EPS) printf("Norm %12.10lf : PASSED\n", norm);
	else printf("Norm %12.10lf : FAILED\n", norm);

	free_arr(&A_res);
}

void construct_A(int nbl, int n2, int n3, double *A, int ldA, double *A_f, int ldaf)
{
	char all = 'A';
	char low = 'L';
	char up = 'U';
	int nbl2 = 2 * nbl;
	int nbl3 = 3 * nbl;

	// Copy block rows 1 and 2
	dlacpy(&all, &nbl2, &nbl3, &A[0 + ldA * 0], &ldA, &A_f[0 + ldaf * 0], &ldaf);

	for (int j = 2; j < nbl - 1; j++)
	{
		// Copy row j
		dlacpy(&all, &nbl, &nbl3, &A[j * nbl + ldA * 0], &ldA, &A_f[j * nbl + ldaf * (j - 1) * nbl], &ldaf);
	}

	// Copy left
	dlacpy(&all, &nbl, &nbl2, &A[(nbl - 1) * nbl + ldA * nbl], &ldA, &A_f[(nbl - 1) * nbl + ldaf * (nbl - 2) * nbl], &ldaf);

#ifdef DEBUG
	print(n2, n2, A_f, ldaf, "A_init");
#endif
}

void construct_L(int nbl, int n2, int n3, double *A, int ldA, double *L_f, int ldlf)
{
	char all = 'A';
	char low = 'L';
	char up = 'U';
	int nbl2 = 2 * nbl;
	int nbl3 = 3 * nbl;

	// Copy j = 0 and j = 1 block rows
	dlacpy(&low, &nbl2, &nbl3, &A[0 + ldA * 0], &ldA, &L_f[0 + ldlf * 0], &ldlf);

	for (int j = 2; j < nbl - 1; j++)
	{
		// Copy left
		dlacpy(&all, &nbl, &nbl, &A[j * nbl + ldA * 0], &ldA, &L_f[j * nbl + ldlf * (j - 1) * nbl], &ldlf);
		// Copy middle
		dlacpy(&low, &nbl, &nbl, &A[j * nbl + ldA * nbl], &ldA, &L_f[j * nbl + ldlf * j * nbl], &ldlf);
	}

	// Copy left
	dlacpy(&all, &nbl, &nbl, &A[(nbl - 1) * nbl + ldA * nbl], &ldA, &L_f[(nbl - 1) * nbl + ldlf * (nbl - 2) * nbl], &ldlf);
	// Copy middle
	dlacpy(&low, &nbl, &nbl, &A[(nbl - 1) * nbl + ldA * 2 * nbl], &ldA, &L_f[(nbl - 1) * nbl + ldlf * (nbl - 1) * nbl], &ldlf);

	for (int i = 0; i < n2; i++)
		L_f[i + ldlf * i] = 1.0;

#ifdef DEBUG
	print(n2, n2, L_f, ldlf,"L_full");
#endif
}

void construct_U(int nbl, int n2, int n3, double *A, int ldA, double *U_f, int lduf)
{
	char all = 'A';
	char low = 'L';
	char up = 'U';
	int nbl2 = 2 * nbl;
	int nbl3 = 3 * nbl;

	// Copy j = 0 and j = 1 block rows
	dlacpy(&up, &nbl2, &nbl3, &A[0 + ldA * 0], &ldA, &U_f[0 + lduf * 0], &lduf);

	for (int j = 2; j < nbl - 1; j++)
	{
		// Copy right
		dlacpy(&all, &nbl, &nbl, &A[j * nbl + ldA * 2 * nbl], &ldA, &U_f[j * nbl + lduf * (j + 1) * nbl], &lduf);
		// Copy middle
		dlacpy(&up, &nbl, &nbl, &A[j * nbl + ldA * nbl], &ldA, &U_f[j * nbl + lduf * j * nbl], &lduf);
	}

	// Copy middle
	dlacpy(&up, &nbl, &nbl, &A[(nbl - 1) * nbl + ldA * 2 * nbl], &ldA, &U_f[(nbl - 1) * nbl + lduf * (nbl - 1) * nbl], &lduf);

#ifdef DEBUG
	print(n2, n2, U_f, lduf, "U_full");
#endif
}

// ---------- Compressed matrices --------------

void SymRecCompress(int n /* order of A */, double *A /* init matrix */, const int lda,  
	const int small_size, int eps, 
	char *method /* SVD or other */)
{
	int ldu, ldv;

	if (n <= small_size)
	{
		return;
	}
	else 
	{
		int n1, n2; // error 3  - неправильное выделение подматриц - похоже на проблему 2
		n2 = (int)ceil(n / 2.0); // округление в большую сторону
		n1 = n - n2; // n2 > n1
		int p = 0; // число значимых сингулярных чисел == rank

#ifdef DEBUG
		printf("SymRecCompress: n = %d n1 = %d n2 = %d\n", n, n1, n2);
#endif

		// LowRank A21
		LowRankApprox(method, n2, n1, &A[n1 + lda * 0], lda, &A[0 + lda * n1], lda, p);

		SymRecCompress(n1, &A[0 + lda * 0], lda, small_size, EPS, method);
		SymRecCompress(n2, &A[n1 + lda * n1], lda, small_size, EPS, method);
	}
}

// Low Rank approximation
void LowRankApprox(char *method, int n2, int n1 /* size of A21 = A */,
					double *A /* A is overwritten by U */, int lda, double *V /* V is stored in A12 */, int ldv, int &p)
{
	char over = 'O';
	char all = 'A';
	char sing = 'S';
	int mn = min(n1, n2);
	int info = 0;
	int lwork = -1;
	p = 0;

	double wkopt;
	double *work;
	double *S = new double[mn];

	if (compare_str(3, method, "SVD"))
	{
		// query 
		dgesvd(&over, &sing, &n2, &n1, A, &lda, S, V, &ldv, V, &ldv, &wkopt, &lwork, &info); // first V - not referenced
		lwork = (int)wkopt;
		work = (double*)malloc(lwork * sizeof(double));
		// error 1

		// A = U1 * S * V1
		dgesvd(&over, &sing, &n2, &n1, A, &lda, S, V, &ldv, V, &ldv, work, &lwork, &info); // first V - not referenced
		// error 2 (как mkl складывает вектора columnwise)

		for (int j = 0; j < mn; j++)
		{
			double s1 = S[j] / S[0];
			if (s1 < EPS)
			{
				break;
			}
			p = j + 1;
			for (int i = 0; i < n2; i++)
				A[i + lda * j] *= S[j];
		}

#ifdef DEBUG
		printf("LowRank after SVD: n2 = %d, n1 = %d, p = %d\n", n2, n1, p);
#endif
							// n1
		for (int j = p; j < mn; j++)   // original part: [n2 x n1], but overwritten part [n2 x min(n2,n1)]
			for (int i = 0; i < n2; i++)
				A[i + lda * j] = 0;

		for (int j = 0; j < n1; j++)   // transposed part: [min(n2,n1) x n1] 
			for (int i = p; i < mn; i++)
				V[i + ldv * j] = 0;

		free(work);
	}
	else
	{
		return;
	}

	free(S);
}

void Test_LowRankApprox_InitA(int m, int n, double eps, char *method)
{
	printf("Test for transposition of LowRank ");
	// A - matrix in dense order
	double *A = alloc_arr(m * n);
	double *A_init = alloc_arr(m * n);
	double *A_rec = alloc_arr(m * n);
	int lda = m;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			A[i + lda * j] = 1.0 / (n + i + j + 1);
			A_init[i + lda * j] = 1.0 / (n + i + j + 1);
		}

	int mn = min(m, n);
	double* V = alloc_arr(mn * n);
	int ldv = mn;
	int p = 0; // p <= mn

	double alpha = 1.0;
	double beta = 0.0;

	LowRankApprox("SVD", m, n, A, lda, V, ldv, p);

	dgemm("no", "no", &m, &n, &mn, &alpha, A, &lda, V, &ldv, &beta, A_rec, &lda);

	rel_error(m, n, A_rec, A_init, lda, eps);

}


/*
Test_LowRankApprox(int m, int n, double *A, bool fill_mat, double eps, char *method)
{
	// A - matrix in dense order

}*/
int compare_str(int n, char *s1, char *s2)
{
	for (int i = 0; i < n; i++)
	{
		if (s1[i] != s2[i]) return 0;
	}
	return 1;
}

void SymResRestore(int n, double *H1 /* compressed */, double *H2 /* recovered */, int ldh, int small_size)
{
	int n1, n2;
	double alpha = 1.0;
	double beta = 0.0;
	char notrans = 'N';
	char trans = 'T';
	char left = 'L';
	char right = 'R';

	n2 = (int)ceil(n / 2.0); // округление в большую сторону
	n1 = n - n2;           

	if (n <= small_size)     // error 4 - не копировалась матрица в этом случае
	{
		dlacpy("All", &n, &n, H1, &ldh, H2, &ldh);
		return;
	}
	else
	{
		// A21 = A21 * A12
		dgemm(&notrans, &notrans, &n2, &n1, &n1, &alpha, &H1[n1 + ldh * 0], &ldh, &H1[0 + ldh * n1], &ldh, &beta, &H2[n1 + ldh * 0], &ldh);

		// A12 = A21*T = A12*T * A21*T
		dgemm(&trans, &trans, &n1, &n2, &n1, &alpha, &H1[0 + ldh * n1], &ldh, &H1[n1 + ldh * 0], &ldh, &beta, &H2[0 + ldh * n1], &ldh);
	

		SymResRestore(n1, &H1[0 + ldh * 0], &H2[0 + ldh * 0], ldh, small_size);
		SymResRestore(n2, &H1[n1 + ldh * n1], &H2[n1 + ldh * n1], ldh, small_size);
	}
}

void Test_SymRecCompress(int n, double eps, char *method, int smallsize)
{
	printf("*****Test for SymRecCompress  n = %d ******* ", n);
	int small_size = 3;
	char frob = 'F';
	double norm = 0;

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
			H[i + ldh * j] = 1.0 / (i * i + j * j + 1);
			H1[i + ldh * j] = 1.0 / (i * i + j * j + 1);
		}
#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
#endif
	SymRecCompress(n, H1, ldh, small_size, EPS, "SVD");
	SymResRestore(n, H1, H2, ldh, small_size);

#ifdef DEBUG
	print(n, n, H1, ldh, "H1 compressed");
	print(n, n, H2, ldh, "H recovered");
#endif

	// Norm of residual || A - L * U ||
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			H2[i + ldh * j] = H2[i + ldh * j] - H[i + ldh * j];

		}
	}

#ifdef DEBUG
	print(n, n, H, ldh,"H init");
	print(n, n, H2, ldh, "diff");
#endif

	norm = dlange(&frob, &n, &n, H2, &ldh, NULL);
	norm = norm / dlange(&frob, &n, &n, H, &ldh, NULL);
	if (norm < EPS) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, EPS);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, EPS);

	free_arr(&H);
	free_arr(&H2);
	free_arr(&H1);
}


// Test for the whole solver

void GenMatrixandRHSandSolution(const int n1, const int n2, const int n3, double *D, int ldd, double *B, int ldb, double *x1, double *f)
{
	// 1. Аппроксимация двумерной плоскости

	int n = n1 * n2; // size of blocks
	int nbr = n3; // number of blocks in one row
	int NBF = nbr * nbr; // full number of blocks in matrix
	double *DD = alloc_arr(n * n); // память под двумерный диагональный блок
	int lddd = n;

	// переделать все это размеры

	// size DD
	int m0 = n;
	int n0 = n;

	// diagonal blocks in dense format
	for (int i = 0; i < n; i++)
		DD[i + lddd * i] = 6.0;

	for (int j = 1; j < n; j++)  // count - until the end
	{
		DD[j + lddd * (j - 1)] = -1.0;
		DD[j - 1 + lddd * j] = -1.0;
	}

	for (int j = n1; j < n; j++) // count - until the end
	{
		DD[j + lddd * (j - n1)] = -1.0;
		DD[j - n1 + lddd * j] = -1.0;
	}

	//print(n, n, DD, lddd);

	// packing into sparse format
	// 5 diagonal matrix with size n2 * nbr + 2 * (n2 * nbr - 1) + 2 * (n2 * nbr - n1)

	int sparse_size = n + 2 * (n - 1) + 2 * (n - n1);
	double *d = (double*)malloc(sparse_size * sizeof(double));
	int *i_ind = (int*)malloc(sparse_size * sizeof(int));
	int *j_ind = (int*)malloc(sparse_size * sizeof(int));

	printf("sparse_size = %d\n", sparse_size);
	map<vector<int>, double> SD;
	SD = dense_to_sparse(n, n, DD, lddd, i_ind, j_ind, d);
	print_map(SD);
	free(DD);

	// Using only sparse matrix D - for 3D. Transfer each 2D block SD to 3D matrix D

	srand((unsigned int)time(0));
	for (int j = 0; j < nbr; j++)
		for (int i = 0; i < n; i++)
	{
			x1[i + n * j] = random(0.0, 1.0);
			if (j < nbr - 1) B[i + n * j] = -1.0;
	}

	// f[i + n * 0]
	// попробуем использовать один и тот же 2D вектор SD для всех блоков NB
	for (int i = 0; i < n; i++)
	{
		f[i + n * 0] += B[i + n * 0] * x1[i + n * 1]; // f[1]
		f[i + n * (nbr - 1)] = B[i + n * (nbr - 2)] * x1[i + n * (nbr - 2)]; // f[NB]
		for (int j = 0; j < n; j++)
		{
			vector<int> vect = { i,  j };
			if (SD.count(vect))
			{
				f[i + n * 0] += SD[vect] * x1[j + n * 0];
				f[i + n * (nbr - 1)] += SD[vect] * x1[j + n * (nbr - 1)];
			}
		}
	}

	for (int blk = 1; blk < nbr - 1; blk++)
		for (int i = 0; i < n; i++)
		{
			f[i + n * blk] += B[i + n * (blk - 1)] * x1[i + n * (blk - 1)] + B[i + n * blk] * x1[i + n * (blk + 1)]; ;
			for (int j = 0; j < n; j++)
			{
				vector<int> vect = { i,  j };
				if (SD.count(vect))
				{
					f[i + n * blk] += SD[vect] * x1[j + n * blk];
				}
			}
	}

}
#if 0
void Block3DSPDSolveFast(int n1, int n2, int n3, double *D, int ldd, double *B, int ldb, double *f, double thresh, int smallsize, int ItRef, char *bench,
			/* output */ double *G, int ldg, double *x, int &success, double &RelRes, int &itcount)
{
	int size = n1 * n2 * n3;
	double tt;
	double tt1;
	double eps = 5e-2;
	tt = omp_get_wtime();
	DirFactFastDiag(D, ldd, B, ldb, eps, smallsize, bench, G);
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display"))
	{
		printf("Total factorization time: %lf\n", tt);
	}

	tt = omp_get_wtime();
	DirSolveFastDiag(G, ldg, B, ldb, f, x);
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display"))
	{
		printf("Solving time: %lf\n", tt);
	}

	double *g = alloc_arr(size);
	double *x1 = alloc_arr(size);
	Resid(D, ldd, B, ldb, x, f, g, RelRes);
	RelRes = 1;
	if (RelRes < thresh)
	{
		success = 1;
		itcount = 0;
	}
	else {
		int success = 0;
		if (ItRef > 0) {
			if (compare_str(7, bench, "display")) printf("Iterative refinement started\n");
			tt1 = omp_get_wtime();
			itcount = 0;
			while ((RelRes > thresh) && (itcount < ItRef))
			{
				tt = omp_get_wtime();
				DirSolveFastDiag(n1, n2, n3, G, ldg, B, ldb, g, x1);
				for (int i = 0; i < size; i++)
					x[i] = x[i] + x1[i];

				Resid(D, ldd, B, ldb, x, f, g, RelRes); // начальное решение f сравниваем с решением A_x0 + A_x1 + A_x2, где
				itcount = itcount + 1;
				tt = omp_get_wtime() - tt;
				if (compare_str(7, bench, "display")) printf("itcount=%d, RelRes=%lf, Time=%lf\n", itcount, RelRes, tt);
			}
			if ((RelRes < thresh) && (itcount < ItRef)) success = 1; // b

			tt1 = omp_get_wtime() - tt1;
			if (compare_str(7, bench, "display")) printf("Iterative refinement total time: %lf\n", tt1);
		}
	}

}

/* Функция вычисления разложения симметричной блочно-диагональной матрицы с использование сжатого формата. 
   Внедиагональные блоки предполагаются диагональными матрицами */
void DirFactFastDiag(int n1, int n2, int n3, double *D, int ldd, double *B, int ldb, double eps, int smallsize, char *bench,
					 double *G /*factorized matrix*/, int ldg)
{
	int n = n1 * n2;
	int NB = n3; // size of D = NB blocks by n elements
	double *DC = alloc_arr(n);  int lddc = n2;
	double *TD1 = alloc_arr(n); int ldtd = n2;
	double *TD = alloc_arr(n);
	

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("Timing DirFactFastDiag\n");
		printf("****************************\n");
	}

	double tt = omp_get_wtime();
	SymRecCompress(n, &D[0 * NB + ldd * 0], ldd, smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display")) printf("Compressing D(0) time: %lf", tt);

	tt = omp_get_wtime();
	SymCompRecInv(DC, lddc, G, ldg, eps, "SVD");
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display")) printf("Computing G(1) time: %lf", tt);


	for (int k = 1; k < NB; k++)
	{
		tt = omp_get_wtime();
		SymRecCompress(n, &D[k * n + ldd * 0], ldd, smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Compressing D(%d) time: %lf", k, tt);
		tt = omp_get_wtime();
		DiagMult(TD1, ldtd, &B[(k - 1) * n + 0 * ldb]);
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Mult D(%d) time: %lf", k, tt);
		tt = omp_get_wtime();
		Add(1, DC, lddc, -1, TD1, ldtd, TD, ldtd, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Add %d time: %lf", k, tt);
		tt = omp_get_wtime();
		SymCompRecInv(TD, ldtd, &G[k * n + ldg * 0], ldg, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Computing G(%d) time: %lf", k, tt);
	}

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("End of DirFactFastDiag\n");
		printf("****************************\n");
	}
}

void DirSolveFastDiag(int n1, int n2, int n3, double *G, int ldg, double *B, int ldb, double *f, double *x)
{
	int n = n1 * n2;
	int NB = n3;
	double *tb = alloc_arr(n * NB);
	double *y = alloc_arr(n);

	for (int i = 0; i < n; i++)
		tb[i] = f[i];

	for (int k = 1; k < NB; k++)
	{
		RecMult(&G[(k - 1) * n + ldg * 0], ldg, &tb[(k - 1) * n], y);
		for (int i = 0; i < n; i++)
			tb[i + k * n] = f[i + k * n] - B[i + (k - 1) * n] * y[i];
	}

	RecMult(&G[(NB - 1) * n + ldg * 0], ldg, &tb[(NB - 1) * n], &x[(NB - 1) * n]);

	for (int k = NB - 2; k > 0; k--)
	{
		for (int i = 0; i < n; i++)
			y[i] = tb[i + k * n] - B[i + k * n] * x[i + (k + 1) * n];
		RecMult(&G[k * n + ldg * 0], ldg, y, &x[k * n]);
	}
}

void SymCompRecInv(double *DC, int lddc, double *G, int ldg, double eps, char *str)
{

}

void RecMult(double *G, int ldg, double *tb /* vector */, double *y)
{

}
#endif
// Рекурсивная функция вычисления DAD, где D - диагональная матрица, а A - сжатая
void DiagMult(int n, double *A, int lda, double *d, int small_size)
{

	if (n <= small_size)     // error 4 - не копировалась матрица в этом случае
	{
		for (int j = 0; j < n; j++)
			for (int i = 0; i < n; i++)
			{
				A[i + j * lda] *= d[j]; // справа D - каждый j - ый столбец A умножается на d[j]
				A[i + j * lda] *= d[i]; // слева D - каждая строка A умножается на d[j]
			}
		return;
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // округление в большую сторону
		int n1 = n - n2;

		DiagMult(n1, &A[0 + lda * 0], lda, &d[0], small_size);
		DiagMult(n2, &A[n1 + lda * n1], lda, &d[n1], small_size);

		// D * U - каждая i-ая строка U умножается на элемент вектора d[i]
		for (int j = 0; j < n1; j++)
			for (int i = 0; i < n2; i++)
				A[i + n1 + lda * (0 + j)] *= d[n1 + i]; // вторая часть массива D

		// VT * D - каждый j-ый столбец умножается на элемент вектора d[j]
		for (int j = 0; j < n2; j++)
			for (int i = 0; i < n1; i++)
				A[i + 0 + lda * (n1 + j)] *= d[j]; 
		// так так вектора матрицы V из разложения A = U * V лежат в транспонированном порядке,
		// то матрицу D стоит умножать на VT слева
	}
}

void Resid(double *D, int ldd, double *B, int ldb, double *x, double *f, double *g, double &RelRes)
{

}

void Test_DiagMult(int n, double eps, char *method, int smallsize)
{
	printf("*****Test for DiagMult  n = %d ******* ", n);
	double *Hd = new double[n*n]; // diagonal Hd = D * H * D
	double *H1 = new double[n*n]; // compressed H
	double *H2 = new double[n*n]; // recovered H after D * H1 * D
	double *d = new double[n];
	double norm = 0;

	Hd[0:n*n] = 0;
	H1[0:n*n] = 0;
	H2[0:n*n] = 0;

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
	rel_error(n, n, H2, Hd, ldh, eps);
	
}

void rel_error(int n, int k, double *Hrec, double *Hinit, int ldh, double eps)
{
	double norm;
	// Norm of residual
	for (int j = 0; j < k; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Hrec[i + ldh * j] = Hrec[i + ldh * j] - Hinit[i + ldh * j];

		}
	}

	norm = dlange("Frob", &n, &k, Hrec, &ldh, NULL);
	norm = norm / dlange("Frob", &n, &k, Hinit, &ldh, NULL);
	if (norm < EPS) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);
}

/* Y = A * X, where A - compressed n * n, X - dense n * m, Y - dense n * m */
void RecMultL(int n, int m, double *A, int lda, double *X, int ldx, double *Y, int ldy, int smallsize)
{
	double alpha = 1.0;
	double beta = 0.0;
	if (n <= smallsize)
	{
		dgemm("No", "No", &n, &m, &n, &alpha, A, &lda, X, &ldx, &beta, Y, &ldy);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // округление в большую сторону
		int n1 = n - n2;
		double *Y12 = alloc_arr(n1 * m);
		double *Y21 = alloc_arr(n2 * m);
		double *Y11 = alloc_arr(n1 * m);
		double *Y22 = alloc_arr(n2 * m);
		double *inter1 = alloc_arr(n2 * n1); // column major - lda = column
		double *inter2 = alloc_arr(n1 * n2);
	
		// A21 = A21 * A12 (the result of multiplication is A21 matrix with size n2 x n1)
		dgemm("No", "No", &n2, &n1, &n1, &alpha, &A[n1 + lda * 0], &lda, &A[0 + lda * n1], &lda, &beta, inter1, &n2);
		
		// Y21 = inter1 (n2 x n1) * X(1...n1, :) (n1 x n)
		dgemm("No", "No", &n2, &m, &n1, &alpha, inter1, &n2, &X[0 + 0 * ldx], &ldx, &beta, Y21, &n2);
	
		// A12 = A21*T = A12*T * A21*T (the result of multiplication is A21 matrix with size n1 x n2)
		dgemm("Trans", "Trans", &n1, &n2, &n1, &alpha, &A[0 + lda * n1], &lda, &A[n1 + lda * 0], &lda, &beta, inter2, &n1);
	
		// Y12 = inter2 (n1 x n2) * X(n1...m, :) (n2 x n)
		dgemm("No", "No", &n1, &m, &n2, &alpha, inter2, &n1, &X[n1 + 0 * ldx], &ldx, &beta, Y12, &n1); // уже транспонировали матрицу в предыдущем dgemm

		RecMultL(n1, m, &A[0 + lda * 0], lda, &X[0 + ldx * 0], ldx, Y11, n1, smallsize);
		RecMultL(n2, m, &A[n1 + lda * n1], lda, &X[n1 + ldx * 0], ldx, Y22, n2, smallsize);

		// first part of Y = Y11 + Y12
		op_mat(n1, m, Y11, Y12, n1, '+');
		dlacpy("All", &n1, &m, Y11, &n1, &Y[0 + ldy * 0], &ldy);

		// second part of Y = Y21 + Y22
		op_mat(n2, m, Y21, Y22, n2, '+');
		dlacpy("All", &n2, &m, Y21, &n2, &Y[n1 + ldy * 0], &ldy);

		free_arr(&Y11);
		free_arr(&Y12);
		free_arr(&Y21);
		free_arr(&Y22);
		free_arr(&inter1);
		free_arr(&inter2);
	
	}
}

/* Тест на сравнение результатов умножения Y = H * X сжимаемой матрицы H на произвольную X.
Сравниваются результаты со сжатием и без */
void Test_RecMultL(int n, int k, double eps, char *method, int smallsize)
{
	printf("*****Test for RecMultL  n = %d k = %d ******* ", n, k);
	double *H = new double[n*n]; // init and compressed
	double *X = new double[n*k];
	double *Y = new double[n*k]; // init Y
	double *Y1 = new double[n*k]; // after multiplication woth compressed
	double norm = 0;
	double alpha = 1.0;
	double beta = 0.0;

	int ldh = n;
	int ldy = n;
	int ldx = n;

	H[0:n*n] = 0;
	X[0:n*k] = 0;
	Y[0:n*k] = 0;
	Y1[0:n*k] = 0;

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

	rel_error(n, k, Y1, Y, ldy, eps);
	
#ifdef DEBUG
	print(n, n, H, ldy, "H comp");
	print(n, k, Y1, ldy, "Y1 rec");
#endif
}

// Функция вычисления линейной комбинации двух сжатых матриц
void Add(int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc, int smallsize, double eps, char *method)
{
	double alpha_loc = 1.0;
	double beta_loc = 0.0;
#ifdef DEBUG
	printf("******Function: Add*******\n");
#endif
	// n - order of A, B and C
	if (n <= smallsize)
	{
#ifdef DEBUG
		printf("Smallsize - doing dense addition for n = %d and smallsize = %d\n", n, smallsize);
#endif
		Add_dense(n, n, alpha, A, lda, beta, B, ldb, C, ldc);
	}
	else
	{
		int p = 0;
		int n2 = (int)ceil(n / 2.0); // округление в большую сторону
		int n1 = n - n2;

		int n1_dbl = 2 * n1;
		int mn = min(n2, n1_dbl);
		
		double *Y21 = alloc_arr(n2 * n1_dbl); int ldy21 = n2;
		double *Y12 = alloc_arr(n1_dbl * n2); int ldy12 = n1_dbl;

		// if m_approx < n_approx => нужно больше памяти
		double *V21 = alloc_arr(n2 * n1_dbl); // p <= n2
		int ldv21 = n2;

		double *V12 = alloc_arr(n1_dbl * n2); // p <= n2
		int ldv12 = n1_dbl;

		// Y = V21'*V12;
		double *Y = alloc_arr(mn * mn); int ldy = mn;

		Add_dense(n2, n1, alpha, &A[n1 + lda * 0], lda, 0.0, B, ldb, &A[n1 + lda * 0], lda);
		Add_dense(n2, n1, beta, &B[n1 + ldb * 0], ldb, 0.0, B, ldb, &B[n1 + ldb * 0], ldb);

		// Y21 = [alpha*A{2,1} beta*B{2,1}];
		dlacpy("All", &n2, &n1, &A[n1 + lda * 0], &lda, &Y21[0 + ldy21 * 0], &ldy21);
		dlacpy("All", &n2, &n1, &B[n1 + ldb * 0], &ldb, &Y21[0 + ldy21 * n1], &ldy21);

		// Y12 = [A{1,2}; B{1,2}];
		dlacpy("All", &n1, &n2, &A[0 + lda * n1], &lda, &Y12[0 + ldy12 * 0], &ldy12);
		dlacpy("All", &n1, &n2, &B[0 + ldb * n1], &ldb, &Y12[n1 + ldy12 * 0], &ldy12);

		// произведение Y21 и Y12 - это матрица n2 x n1
		LowRankApprox("SVD", n2, n1_dbl, Y21, ldy21, V21, ldv21, p); // перезапись Y21
		LowRankApprox("SVD", n1_dbl, n1, Y12, ldy12, V12, ldv12, p);  // перезапись Y12


		// Y = V21'*V12;
		dgemm("No", "No", &mn, &mn, &n1_dbl, &alpha_loc, V21, &ldv21, Y12, &ldy12, &beta_loc, Y, &ldy);

		// C{2,1} = U21*Y;   
		dgemm("No", "No", &n2, &n1, &mn, &alpha_loc, Y21, &ldy21, Y, &ldy, &beta_loc, &C[n1 + ldc * 0], &ldc);

		// C{1,2} = U12';
		dlacpy("All", &n1, &n2, V12, &ldv12, &C[0 + ldc * n1], &ldc);


		Add(n1, alpha, &A[0 + lda * 0], lda, beta, &B[0 + ldb * 0], ldb, &C[0 + ldc * 0], ldc, smallsize, eps, method);
		Add(n2, alpha, &A[n1 + lda * n1], lda, beta, &B[n1 + ldb * n1], ldb, &C[n1 + ldc * n1], ldc, smallsize, eps, method);

		free_arr(&Y21);
		free_arr(&Y12);
		free_arr(&V21);
		free_arr(&V12);
		free_arr(&Y);
	}

}

void Test_add(int n, double alpha, double beta, double smallsize, double eps, char *method)
{
	printf("*****Test for Add n = %d ******* ", n);
	double *H1 = alloc_arr(n*n);
	double *H2 = alloc_arr(n*n);
	double *G = alloc_arr(n*n);
	double *H1c = alloc_arr(n*n);
	double *H2c = alloc_arr(n*n);
	double *Gc = alloc_arr(n*n);
	double *GcR = alloc_arr(n*n);
	int ldh = n;
	int ldg = n;

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
	rel_error(n, n, GcR, G, ldg, eps);

	free_arr(&H1);
	free_arr(&H2);
	free_arr(&G);
	free_arr(&H1c);
	free_arr(&H2c);
	free_arr(&Gc);
	free_arr(&GcR);
}

void Test_transpose(int m, int n, int smallsize, double eps, char *method)
{
	double *H = alloc_arr(m * n); int ldh = m;
	double *H_rec = alloc_arr(m * n);
	int mn = min(m, n);
	double *V = alloc_arr(mn * n); int ldv = mn;
	double *Htr = alloc_arr(n * m); int ldhtr = n;
	double *Hcomp_tr = alloc_arr(m * n);
	double *Vcomp_tr = alloc_arr(n * m);
	double *Vtr = alloc_arr(mn * m); int ldvtr = mn;
	double alpha = 1.0;
	double beta = 0.0;

	int p = 0;
	int k = 0;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			k++;
			Htr[j + ldhtr * i] = 1.0 / (i + j + 1 + k);
			H[i + ldh * j] = 1.0 / (i + j + 1 + k);
		}


	LowRankApprox("SVD", n, m, Htr, ldhtr, Vtr, ldvtr, p);

	dgemm("Trans", "Trans", &m, &n, &mn, &alpha, Vtr, &ldvtr, Htr, &ldhtr, &beta, H_rec, &ldh);
	
	// || H _rec  - H || / || H ||
	rel_error(m, n, H_rec, H, ldh, eps);

}

/* Функция вычисления симметричного малорангового дополнения A:= A + alpha * V * Y * V'
A - симметрическая сжатая (n x n)
Y - плотная симметричная размера k x k, k << n , V - плотная прямоугольная n x k
(n x n) = (n x n) + (n x k) * (k x k) * (k * n) */
void SymCompUpdate2(int n, int k, double *A, int lda, double alpha, double *Y, int ldy, double *V, int ldv, double *B, int ldb, int smallsize, double eps, char* method)
{
	double alpha_one = 1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;

	int p = 0;
	if (n <= smallsize)
	{
		// X = X + alpha * V * Y * VT

		// C = V * Y
		double *C = alloc_arr(n * k); int ldc = n;
		dsymm("Right", "Up", &n, &k, &alpha_one, Y, &ldy, V, &ldv, &beta_zero, C, &ldc);

		// X = X + alpha * C * Vt
		dgemm("No", "Trans", &n, &n, &k, &alpha, C, &ldc, V, &ldv, &beta_one, A, &lda);

		// B = A
		dlacpy("All", &n, &n, A, &lda, B, &ldb);

		free_arr(&C);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // n2 > n1
		int n1 = n - n2;

		int nk = n1 + k;
		// for this division n2 > n1 we can store a low memory
		double *Y12 = alloc_arr(nk * n1); int ldy12 = nk;
		double *Y21 = alloc_arr(n2 * nk); int ldy21 = n2;

		double *V_uptr = alloc_arr(k * n1); int ldvuptr = k;
		double *VY = alloc_arr(n2 * k); int ldvy = n2;

		int mn = min(n2, nk);
		double *V21 = alloc_arr(mn * nk); int ldv21 = mn;
		double *V12 = alloc_arr(n1 * n1); int ldv12 = n1;

		double *VV = alloc_arr(mn * n1); int ldvv = mn;

		dgemm("No", "No", &n2, &k, &k, &alpha, &V[n1 + ldv * 0], &ldv, Y, &ldy, &beta_zero, VY, &ldvy);

		// Y21 = [A{2,1} alpha*V(m:n,:)*Y];
		dlacpy("All", &n2, &n1, &A[n1 + lda * 0], &lda, &Y21[0 + ldy21 * 0], &ldy21);
		dlacpy("All", &n2, &k, VY, &ldvy, &Y21[0 + ldy21 * n1], &ldy21);

		Mat_Trans(n1, k, &V[0 + ldv * 0], ldv, V_uptr, ldvuptr);
	
		// Y12 = [A{1,2} V(1:n1,:)];
		dlacpy("All", &n1, &n1, &A[0 + lda * n1], &lda, &Y12[0 + ldy12 * 0], &ldy12);
		dlacpy("All", &k, &n1, V_uptr, &ldvuptr, &Y12[n1 + ldy21 * 0], &ldy12);

		// [U21,V21] = LowRankApprox (Y21, eps, method);
		LowRankApprox("SVD", n2, nk, Y21, ldy21, V21, ldv21, p);

		// [U12, V12] = LowRankApprox(Y12, eps, method);
		LowRankApprox("SVD", nk, n1, Y12, ldy12, V12, ldv12, p);


		// B{2,1} = U21*(V21'*V12);

		// V21 * Y12
		dgemm("No", "No", &mn, &n1, &nk, &alpha_one, V21, &ldv21, Y12, &ldy12, &beta_zero, VV, &ldvv);
		dgemm("No", "No", &n2, &n1, &mn, &alpha_one, Y21, &ldy21, VV, &ldvv, &beta_zero, &B[n1 + ldb * 0], &ldb);
		
		// B{1,2} = U12;
		dlacpy("All", &n1, &n1, V12, &ldv12, &B[0 + ldb * n1], &ldb);

		// B{1,1} = SymCompUpdate2 (A{1,1}, Y, V(1:n1,:), alpha, eps, method);
		SymCompUpdate2(n1, k, &A[0 + lda * 0], lda, alpha, Y, ldy, &V[0 + ldv * 0], ldv, &B[0 + ldb * 0], ldb, smallsize, eps, method);
		
		// B{2,2} = SymCompUpdate2 (A{2,2}, Y, V(m:n ,:), alpha, eps, method);
		SymCompUpdate2(n2, k, &A[n1 + lda * n1], lda, alpha, Y, ldy, &V[n1 + ldv * 0], ldv, &B[n1 + ldb * n1], ldb, smallsize, eps, method);

		free_arr(&Y21);
		free_arr(&Y12);
		free_arr(&V21);
		free_arr(&V12);
		free_arr(&VY);
		free_arr(&V_uptr);
		free_arr(&VV);
	}
}

// || B - B_rec || / || B ||
void Test_SymCompUpdate2(int n, int k, double alpha, int smallsize, double eps, char* method)
{
	printf("*****Test for SymCompUpdate2   n = %d k = %d ***** ", n, k);
	double alpha_one = 1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;

	// B = H - V * Y * VT
	double *B = alloc_arr(n * n); int ldb = n;
	double *B_rec = alloc_arr(n * n);
	double *Y = alloc_arr(k * k); int ldy = k;
	double *V = alloc_arr(n * k); int ldv = n; int ldvtr = k;
	double *HC = alloc_arr(n * n); int ldh = n;
	double *H = alloc_arr(n * n); 
	double *C = alloc_arr(n * k); int ldc = n;

	Hilbert(n, HC, ldh);
	Hilbert(n, H, ldh);

	for (int i = 0; i < k; i++)
		Y[i + ldy * i] = i + 1;

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
	rel_error(n, n, B_rec, H, ldh, eps);
	
	free_arr(&B);
	free_arr(&B_rec);
	free_arr(&H);
	free_arr(&HC);
	free_arr(&Y);
	free_arr(&C);
	free_arr(&V);
}

// Рекурсивное обращение сжатой матрицы
void SymCompRecInv(int n, double *A, int lda, double *B, int ldb, int smallsize, double eps, char *method)
{
	double alpha_one = 1.0;
	double alpha_mone = -1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;
	int info = 0;
	double wquery = 0;
	int lwork = -1;

	if (n <= smallsize)
	{
		int *ipiv = (int*)malloc(n * sizeof(int));

		// LU factorization of A
		dgetrf(&n, &n, A, &lda, ipiv, &info);

		// space query
		dgetri(&n, A, &lda, ipiv, &wquery, &lwork, &info);

		lwork = (int)wquery;
		double *work = alloc_arr(lwork);

		// inversion of A
		dgetri(&n, A, &lda, ipiv, work, &lwork, &info);

		// dlacpy 
		dlacpy("All", &n, &n, A, &lda, B, &ldb);

		free(ipiv);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // n2 > n1
		int n1 = n - n2;

		double *X11 = alloc_arr(n1 * n1); int ldx11 = n1;
		double *X22 = alloc_arr(n2 * n2); int ldx22 = n2;
		double *Y = alloc_arr(n1 * n1); int ldy = n1;
		double *V = alloc_arr(n1 * n1); int ldv = n1;
		double *B12 = alloc_arr(n1 * n1); int ldb12 = n1;

		// Inversion of A22 to X22
		SymCompRecInv(n2, &A[n1 + lda * n1], lda, X22, ldx22, smallsize, eps, method);

		// Save X22 * U to B{2,1}
		RecMultL(n2, n1, X22, ldx22, &A[n1 + lda * 0], lda, &B[n1 + ldb * 0], ldb, smallsize);

		// Compute Y = UT * X22 * U = | A[2,1]T * B{2,1} | = | (n1 x n2) x (n2 x n1) |
		dgemm("Trans", "No", &n1, &n1, &n2, &alpha_one, &A[n1 + lda * 0], &lda, &B[n1 + ldb * 0], &ldb, &beta_zero, Y, &ldy);

		// Update X11 = A11 - V * UT * X22 * U * VT = | A11 - V * Y * VT | = | (n1 x n1) - (n1 x n1) * (n1 x n1) * (n1 x n1) |
		Mat_Trans(n1, n1, &A[0 + lda * n1], lda, V, ldv);
		SymCompUpdate2(n1, n1, &A[0 + lda * 0], lda, alpha_mone, Y, ldy, V, ldv, X11, ldx11, smallsize, eps, method);

		// Inversion of X11 to B11
		SymCompRecInv(n1, X11, ldx11, &B[0 + ldb * 0], ldb, smallsize, eps, method);

		// Fill B{1,2} as B12 = -B{1,1} * A{1,2} = -X11 * V = (n1 x n1) * (n1 x n1) in real (x2 x x1)
		RecMultL(n1, n1, &B[0 + ldb * 0], ldb, V, ldv, B12, ldb12, smallsize);
		Add_dense(n1, n1, alpha_mone, B12, ldb12, beta_zero, B, ldb, B12, ldb12);

		// B{1,2} = transpose(B12)
		Mat_Trans(n1, n1, B12, ldb12, &B[0 + ldb * n1], ldb);

		// Y = -(A{1,2})' * B{1,2} = -VT * (-X11 * V) = - VT * B12; [n1 x n1] * [n1 x n1]
		dgemm("No", "No", &n1, &n1, &n1, &alpha_mone, &A[0 + lda * n1], &lda, B12, &ldb12, &beta_zero, Y, &ldy);

		// Update X22 + (X22*U) * VT * X11 * V (UT * X22) = X22 + B21 * Y * B21T = (n2 x n2) + (n2 x n1) * (n1 x n1) * (n1 x n2)
		SymCompUpdate2(n2, n1, X22, ldx22, alpha_one, Y, ldy, &B[n1 + ldb * 0], ldb, &B[n1 + ldb * n1], ldb, smallsize, eps, method);

		free_arr(&X11);
		free_arr(&X22);
		free_arr(&Y);
		free_arr(&V);
		free_arr(&B12);
	}
}

void Test_SymCompRecInv(int n, int smallsize, double eps, char *method)
{
	printf("***** Test_SymCompRecInv n = %d **** \n", n);
	double *H = alloc_arr(n * n);
	double *Hc = alloc_arr(n * n);
	double *Bc = alloc_arr(n * n);
	double *Brec = alloc_arr(n * n);
	double *Y = alloc_arr(n * n);

	int ldh = n;
	int ldb = n;
	int ldy = n;

	double alpha_mone = -1.0;
	double beta_one = 1.0;

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
	
	double norm = dlange("Frob", &n, &n, Y, &ldy, NULL);

	if (norm < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);

}

void Eye(int n, double *H, int ldh)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			if (j == i) H[i + ldh * j] = 1.0;
			else H[i + ldh * j] = 0.0;
}

void Hilbert(int n, double *H, int ldh)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			H[i + ldh * j] = 1.0 / (i + j + 1);
}

void Mat_Trans(int m, int n, double *H, int ldh, double *Hcomp_tr, int ldhtr)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
			Hcomp_tr[j + ldhtr * i] = H[i + ldh * j];
}

void Add_dense(int m, int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc)
{

#pragma omp parallel for
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
			C[i + ldc * j] = alpha * A[i + lda * j] + beta * B[i + ldb * j];
	}
}

void op_mat(int m, int n, double *Y11, double *Y12, int ldy, char sign)
{
	if (sign == '+')
	{
		for (int j = 0; j < n; j++)
			for (int i = 0; i < m; i++)
				Y11[i + ldy * j] += Y12[i + ldy *j];
	}
	else if (sign == '-')
	{
		for (int j = 0; j < n; j++)
			for (int i = 0; i < m; i++)
				Y11[i + ldy * j] -= Y12[i + ldy *j];
	}
	else
	{
		printf("Incorrect sign\n");
	}
}


map<vector<int>,double> dense_to_sparse(int m, int n, double *DD, int ldd, int *i_ind, int *j_ind, double *d)
{
	map<vector<int>, double> SD;
	vector<int> v(2);
	double thresh = 1e-8;
	int k = 0;
	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
			if (fabs(DD[i + ldd * j]) != 0)
			{
				d[k] = DD[i + ldd * j];
				i_ind[k] = i;
				j_ind[k] = j;

				v[0] = i;
				v[1] = j;
				SD[v] = DD[i + ldd * j];

				k++;
			}

	return SD;
}

void print_map(const map<vector<int>, double>& SD)
{
	cout << SD.size() << endl;
	for (const auto& item : SD)
		cout << "m = " << item.first[0] << " n = " << item.first[1] << " value = " << item.second << endl;
}

double random(double min, double max)
{
	return (double)(rand()) / RAND_MAX * (max - min) + min;
}

void print_vec(int size, double *vec1, double *vec2, char *name)
{
	printf("%s\n", name);
	for (int i = 0; i < size; i++)
		printf("%d   %lf   %lf\n", i, vec1[i], vec2[i]);
}
