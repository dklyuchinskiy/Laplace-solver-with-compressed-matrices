#include "definitions.h"
#include "templates.h"

using namespace std;

// ---------- Compressed matrices --------------

// Low Rank approximation
void LowRankApprox(int n2, int n1 /* size of A21 = A */,
					double *A /* A is overwritten by U */, int lda, double *V /* V is stored in A12 */, int ldv, int &p, double eps, char *method)
{
	int mn = min(n1, n2);
	int info = 0;
	int lwork = -1;
	p = 0;

	double wkopt;
	double *work;
	double *S;

	if (compare_str(3, method, "SVD"))
	{
		S = alloc_arr(mn);

		// query 
		dgesvd("Over", "Sing", &n2, &n1, A, &lda, S, V, &ldv, V, &ldv, &wkopt, &lwork, &info); // first V - not referenced
		lwork = (int)wkopt;
		work = alloc_arr(lwork);

		// A = U1 * S * V1
		dgesvd("Over", "Sing", &n2, &n1, A, &lda, S, V, &ldv, V, &ldv, work, &lwork, &info); // first V - not referenced
		// error 2 (как mkl складывает вектора columnwise)

		for (int j = 0; j < mn; j++)
		{
			double s1 = S[j] / S[0];
			if (s1 < eps)
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
		
		free_arr(&S);
		free_arr(&work);
	}
	else
	{
		return;
	}
}


void SymRecCompress(int n /* order of A */, double *A /* init matrix */, const int lda,
	const int small_size, double eps,
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
		LowRankApprox(n2, n1, &A[n1 + lda * 0], lda, &A[0 + lda * n1], lda, p, eps, method);

		SymRecCompress(n1, &A[0 + lda * 0], lda, small_size, eps, method);
		SymRecCompress(n2, &A[n1 + lda * n1], lda, small_size, eps, method);
	}
}

void SymResRestore(int n, double *H1 /* compressed */, double *H2 /* recovered */, int ldh, int small_size)
{
	int n1, n2;
	double alpha = 1.0;
	double beta = 0.0;

	if (n <= small_size)     // error 4 - не копировалась матрица в этом случае
	{
		dlacpy("All", &n, &n, H1, &ldh, H2, &ldh);
		return;
	}
	else
	{
		n2 = (int)ceil(n / 2.0); // округление в большую сторону
		n1 = n - n2;

		// A21 = A21 * A12
		dgemm("Notrans", "Notrans", &n2, &n1, &n1, &alpha, &H1[n1 + ldh * 0], &ldh, &H1[0 + ldh * n1], &ldh, &beta, &H2[n1 + ldh * 0], &ldh);

		// A12 = A21*T = A12*T * A21*T
		dgemm("Trans", "Trans", &n1, &n2, &n1, &alpha, &H1[0 + ldh * n1], &ldh, &H1[n1 + ldh * 0], &ldh, &beta, &H2[0 + ldh * n1], &ldh);
	

		SymResRestore(n1, &H1[0 + ldh * 0], &H2[0 + ldh * 0], ldh, small_size);
		SymResRestore(n2, &H1[n1 + ldh * n1], &H2[n1 + ldh * n1], ldh, small_size);
	}
}

// Test for the whole solver

inline int ind(int j, int n)
{
	return n * j;
}

void GenMatrixandRHSandSolution(const int n1, const int n2, const int n3, double *D, int ldd, double *B, double *x1, double *f)
{
	// 1. Аппроксимация двумерной плоскости

	int n = n1 * n2; // size of blocks
	int nbr = n3; // number of blocks in one row
	int NBF = nbr * nbr; // full number of blocks in matrix
	double *DD = alloc_arr(n * n); // память под двумерный диагональный блок
	int lddd = n;
	int size = n * nbr;

	double *f_help = alloc_arr(n);

	double done = 1.0;
	double dzero = 0.0;
	int ione = 1;

	double time1, time2;

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
	//print_map(SD);


	// Using only sparse matrix D - for 3D. Transfer each 2D block SD to 3D matrix D

	GenSolVector(size, x1);

#pragma omp parallel for
	for (int j = 0; j < nbr; j++)
	{
		dlacpy("All", &n, &n, DD, &lddd, &D[ind(j, n) + ldd * 0], &ldd);
		for (int i = 0; i < n; i++)
		{
			if (j < nbr - 1) B[ind(j, n) + i] = -1.0;
		}
	}

#if 0
	time1 = omp_get_wtime();
	// f[i + n * 0]
	// попробуем использовать один и тот же 2D вектор SD для всех блоков NB
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
	{
		f[i + n * 0] = B[i + n * 0] * x1[i + n * 1]; // f[1]
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

#pragma omp parallel for schedule(guided)
	for (int blk = 1; blk < nbr - 1; blk++)
		for (int i = 0; i < n; i++)
		{
			f[i + n * blk] = B[i + n * (blk - 1)] * x1[i + n * (blk - 1)] + B[i + n * blk] * x1[i + n * (blk + 1)];
			for (int j = 0; j < n; j++)
			{
				vector<int> vect = { i,  j };
				if (SD.count(vect))
				{
					f[i + n * blk] += SD[vect] * x1[j + n * blk];
				}
			}
	}
	time1 = omp_get_wtime() - time1;
#endif

	// через mkl

	time2 = omp_get_wtime();
	// f[1] = D[1] * x[1] + diag{B[1]} * x[2]
	DenseDiagMult(n, &B[ind(0, n)], &x1[ind(1, n)], &f[ind(0, n)]);
	dsymv("Up", &n, &done, &D[ind(0, n)], &ldd, &x1[ind(0, n)], &ione, &done, &f[ind(0, n)], &ione);

	// f[N] = diag(B{N-1}) * x{N-1} + D{N} * x{N};
	DenseDiagMult(n, &B[ind(nbr - 2, n)], &x1[ind(nbr - 2, n)], &f[ind(nbr - 1, n)]);
	dsymv("Up", &n, &done, &D[ind(nbr - 1, n)], &ldd, &x1[ind(nbr - 1, n)], &ione, &done, &f[ind(nbr - 1, n)], &ione);


	// f{ i } = diag(B{ i - 1 }) * x{ i - 1 } + D{ i } * x{ i } + diag(B{ i }) * x{ i + 1 };
	for (int blk = 1; blk < nbr - 1; blk++)
	{
		// f{ i } = diag(B{ i - 1 }) * x { i - 1 } + diag(B{ i }) * x { i + 1 };
		DenseDiagMult(n, &B[ind(blk - 1, n)], &x1[ind(blk - 1, n)], &f[ind(blk, n)]);
		DenseDiagMult(n, &B[ind(blk, n)], &x1[ind(blk + 1, n)], f_help);
		daxpby(&n, &done, f_help, &ione, &done, &f[ind(blk, n)], &ione);
		//Add_dense_vect(n, done, &f[ind(blk, n)], done, f_help, &f[ind(blk, n)]);

		// f{i} = f{i} + D{ i } * x{ i } 
		dsymv("Up", &n, &done, &D[ind(blk, n)], &ldd, &x1[ind(blk, n)], &ione, &done, &f[ind(blk, n)], &ione);
	}

	time2 = omp_get_wtime() - time2;

	printf("time_mkl = %lf\n", time2);

	free_arr(&DD);
	free_arr(&d);
	free(i_ind);
	free(j_ind);
	free_arr(&f_help);
}

double F_ex(double x, double y, double z)
{
//	return -12.0 * PI * PI * sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * z);
	return 0;
}

double u_ex(double x, double y, double z)
{
//	return 2.0 + sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * z);
	return x * x + y * y - 2.0 * z * z;
}

void GenerateDiagonal2DBlock(int part_of_field, size_m x, size_m y, size_m z, double *DD, int lddd)
{
	int n = x.n * y.n;
	int size = n * z.n;
	
	// diagonal blocks in dense format
	for (int i = 0; i < n; i++)
	{
		DD[i + lddd * i] = -2.0 * (1.0 / (x.h * x.h) + 1.0 / (y.h * y.h) + 1.0 / (z.h * z.h));
		if (i > 0) DD[i + lddd * (i - 1)] = 1.0 / (x.h * x.h);
		if (i < n - 1) DD[i + lddd * (i + 1)] = 1.0 / (x.h * x.h);
		if (i >= x.n) DD[i + lddd * (i - x.n)] = 1.0 / (y.h * y.h);
		if (i <= n - x.n - 1)  DD[i + lddd * (i + x.n)] = 1.0 / (y.h * y.h);

		if (i % x.n == 0 && i > 0)
		{
			DD[i - 1 + lddd * i] = 0;
			DD[i + lddd * (i - 1)] = 0;
		}
	}

}

void GenRHSandSolution(size_m x, size_m y, size_m z, /* output */ double *u, double *f)
{
	int n = x.n * y.n;

	// approximation of exact right hand side (inner grid points)
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < z.n; k++)
		for (int j = 0; j < y.n; j++)
			for (int i = 0; i < x.n; i++)
				f[k * n + j * x.n + i] = F_ex((i + 1) * x.h, (j + 1) * y.h, (k + 1) * z.h);

	// for boundaries z = 0 and z = Lz we distract blocks B0 and Bm from the RHS
#pragma omp parallel for schedule(dynamic,16)
	for (int j = 0; j < y.n; j++)
		for (int i = 0; i < x.n; i++)
		{
			f[ind(0, n) + ind(j, x.n) + i] -= u_ex((i + 1) * x.h, (j + 1) * y.h, 0) / (z.h * z.h); // u|z = 0
			f[ind(z.n - 1, n) + ind(j, x.n) + i] -= u_ex((i + 1)  * x.h, (j + 1) * y.h, z.l) / (z.h * z.h); // u|z = h
		}


	// for each boundary 0 <= z <= Lz
	// we distract 4 known boundaries f0, fl, g0, gL from right hand side
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < z.n; k++)
	{
		for (int i = 0; i < x.n; i++)
		{
			f[k * n + 0 * x.n + i] -= u_ex((i + 1) * x.h, 0, (k + 1) * z.h) / (y.h * y.h);
			f[k * n + (y.n - 1) * x.n + i] -= u_ex((i + 1) * x.h, y.l, (k + 1) * z.h) / (y.h * y.h);
		}
		for (int j = 0; j < y.n; j++)
		{
			f[k * n + j * x.n + 0] -= u_ex(0, (j + 1) * y.h, (k + 1) * z.h) / (x.h * x.h);
			f[k * n + j * x.n + x.n - 1] -= u_ex(x.l, (j + 1) * y.h, (k + 1) * z.h) / (x.h * x.h);
		}
	}

	// approximation of inner points values
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < z.n; k++)
		for (int j = 0; j < y.n; j++)
			for (int i = 0; i < x.n; i++)
				u[ind(k, n) + ind(j, x.n) + i] = u_ex((i + 1) * x.h, (j + 1) * y.h, (k + 1) * z.h);

	printf("RHS and solution are constructed\n");
}

void GenMatrixandRHSandSolution2(size_m x, size_m y, size_m z,
	/* output */ double *D, int ldd, double *B, double *u, double *f, double thresh)
{
	int n = x.n * y.n; // size of blocks
	int nbr = z.n; // number of blocks in one row
	int lddd = n;
	int size = n * nbr;
	double done = 1.0;
	double dzero = 0.0;
	int ione = 1;

	double *DD = alloc_arr(n * n); // 2D diagonal template block
	double *Au = alloc_arr(size); // mult of generated A and exact solution u

	// n - number of unknowns
	// n * n * n - all unknowns

	// f_rhs = f_inner + f_bound

	// approximation of exact right hand side (inner grid points)
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < nbr; k++)
		for (int j = 0; j < y.n; j++)
			for (int i = 0; i < x.n; i++)
				f[k * n + j * x.n + i] = F_ex((i + 1) * x.h, (j + 1) * y.h, (k + 1) * z.h);

	// for boundaries z = 0 and z = Lz we distract blocks B0 and Bm from the RHS
#pragma omp parallel for schedule(dynamic,16)
	for (int j = 0; j < y.n; j++)
		for (int i = 0; i < x.n; i++)
		{
			f[ind(0, n) + ind(j, x.n) + i] -= u_ex((i + 1) * x.h, (j + 1) * y.h, 0) / (z.h * z.h); // u|z = 0
			f[ind(nbr - 1, n) + ind(j, x.n) + i] -= u_ex((i + 1)  * x.h, (j + 1) * y.h, z.l) / (z.h * z.h); // u|z = h
		}


	// for each boundary 0 <= z <= Lz
	// we distract 4 known boundaries f0, fl, g0, gL from right hand side
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < nbr; k++)
	{
			for (int i = 0; i < x.n; i++)
			{
				f[k * n + 0 * x.n + i] -= u_ex((i + 1) * x.h, 0, (k + 1) * z.h) / (y.h * y.h);
				f[k * n + (y.n - 1) * x.n + i] -= u_ex((i + 1) * x.h, y.l, (k + 1) * z.h) / (y.h * y.h);
			}
			for (int j = 0; j < y.n; j++)
			{
				f[k * n + j * x.n + 0] -= u_ex(0, (j + 1) * y.h, (k  + 1) * z.h) / (x.h * x.h);
				f[k * n + j * x.n + x.n - 1] -= u_ex(x.l, (j + 1) * y.h, (k + 1) * z.h) / (x.h * x.h);
			}
	}
//	if (i % x.n == 0 || (i + 1) % x.n == 0) DD[i + lddd * i] = 1.0;

	// Set vector B
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < z.n - 1; j++)
		for (int i = 0; i < n; i++)
			B[ind(j, n) + i] = 1.0 / (z.h * z.h);

	for (int j = 0; j < nbr; j++)
	{
		GenerateDiagonal2DBlock(j, x, y, z, DD, lddd);
		dlacpy("All", &n, &n, DD, &lddd, &D[ind(j, n) + ldd * 0], &ldd);
	}
	
	// approximation of inner points values
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < nbr; k++)
		for (int j = 0; j < y.n; j++)
			for (int i = 0; i < x.n; i++)
				u[ind(k, n) + ind(j, x.n) + i] = u_ex((i + 1) * x.h, (j + 1) * y.h, (k + 1) * z.h);
	
	Mult_Au(x.n, y.n, z.n, D, ldd, B, u, Au);

#ifdef DEBUG
	print_vec(size - n, B, B, "B_vector");
	print_vec_mat(size, n, D, ldd, u, "D and u");
	print_vec(size, Au, f, "Au_ex vs F_ex");
	system("pause");
#endif

	// check error between Au and F
	rel_error(size, 1.0, Au, f, size, thresh);

	free(Au);
	free(DD);

}


void Mult_Au(int n1, int n2, int n3, double *D, int ldd, double *B, double *u, double *Au /*output*/)
{
	int n = n1 * n2;
	int nbr = n3;
	int size = n * nbr;
	double done = 1.0;
	double dzero = 0.0;
	int ione = 1;
	double *f_help = alloc_arr(n);

	// f[1] = D{1} * x{1} + diag(B{1}) * x{2};
	DenseDiagMult(n, &B[ind(0, n)], &u[ind(1, n)], &Au[ind(0, n)]);
	dgemv("No", &n, &n, &done, &D[ind(0, n)], &ldd, &u[ind(0, n)], &ione, &done, &Au[ind(0, n)], &ione);

	// f[N] = diag(B{N-1}) * x{N-1} + D{N} * x{N};
	DenseDiagMult(n, &B[ind(nbr - 2, n)], &u[ind(nbr - 2, n)], &Au[ind(nbr - 1, n)]);
	dgemv("No", &n, &n, &done, &D[ind(nbr - 1, n)], &ldd, &u[ind(nbr - 1, n)], &ione, &done, &Au[ind(nbr - 1, n)], &ione);

	// f{ i } = diag(B{ i - 1 }) * x{ i - 1 } + D{ i } * x{ i } + diag(B{ i }) * x{ i + 1 };
	for (int blk = 1; blk < nbr - 1; blk++)
	{
		// f{ i } = diag(B{ i - 1 }) * x { i - 1 } + diag(B{ i }) * x { i + 1 };
		DenseDiagMult(n, &B[ind(blk - 1, n)], &u[ind(blk - 1, n)], &Au[ind(blk, n)]);
		DenseDiagMult(n, &B[ind(blk, n)], &u[ind(blk + 1, n)], f_help);
		daxpby(&n, &done, f_help, &ione, &done, &Au[ind(blk, n)], &ione);

		// f{i} = f{i} + D{ i } * x{ i }  matrix D - non symmetric
		dgemv("No", &n, &n, &done, &D[ind(blk, n)], &ldd, &u[ind(blk, n)], &ione, &done, &Au[ind(blk, n)], &ione);
	}

	free_arr(&f_help);
}

inline void Add_dense_vect(int n, double alpha, double *a, double beta, double *b, double *c)
{
#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < n; i++)
		c[i] = alpha * a[i] + beta * b[i];
}

// v[i] = D[i] * v[i]
inline void DenseDiagMult(int n, double *diag, double *v, double *f)
{
#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < n; i++)
		f[i] = diag[i] * v[i];
}

void GenSolVector(int size, double *vector)
{
	srand((unsigned int)time(0));
	for (int i = 0; i < size; i++)
		vector[i] = random(0.0, 1.0);
}

//void Resid_CSR(n1, n2, n3, DI, lddi, B, x_sol, f, g, RelRes);

void Block3DSPDSolveFast(int n1, int n2, int n3, double *D, int ldd, double *B, double *f, double thresh, int smallsize, int ItRef, char *bench,
			/* output */ double *G, int ldg, double *x_sol, int &success, double &RelRes, int &itcount)
{
	int size = n1 * n2 * n3;
	int n = n1 * n2;
	double tt;
	double tt1;
	double *DI = alloc_arr(size * n); int lddi = size;
	dlacpy("All", &size, &n, D, &ldd, DI, &lddi);

	tt = omp_get_wtime();
	DirFactFastDiag(n1, n2, n3, D, ldd, B, G, ldg, thresh, smallsize, bench);
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display"))
	{
		printf("Total factorization time: %lf\n", tt);
	}

	tt = omp_get_wtime();
	DirSolveFastDiag(n1, n2, n3, G, ldg, B, f, x_sol, thresh, smallsize);
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display"))
	{
		printf("Solving time: %lf\n", tt);
	}

	double *g = alloc_arr(size);
	double *x1 = alloc_arr(size);
	RelRes = 1;
	Resid(n1, n2, n3, DI, lddi, B, x_sol, f, g, RelRes);

	printf("RelRes = %lf\n", RelRes);
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
				DirSolveFastDiag(n1, n2, n3, G, ldg, B, g, x1, thresh, smallsize);

#pragma omp parallel for simd schedule(simd:static)
				for (int i = 0; i < size; i++)
					x_sol[i] = x_sol[i] + x1[i];

				Resid(n1, n2, n3, DI, lddi, B, x_sol, f, g, RelRes); // начальное решение f сравниваем с решением A_x0 + A_x1 + A_x2, где
				itcount = itcount + 1;
				tt = omp_get_wtime() - tt;
				if (compare_str(7, bench, "display")) printf("itcount=%d, RelRes=%lf, Time=%lf\n", itcount, RelRes, tt);
			}
			if ((RelRes < thresh) && (itcount < ItRef)) success = 1; // b

			tt1 = omp_get_wtime() - tt1;
			if (compare_str(7, bench, "display")) printf("Iterative refinement total time: %lf\n", tt1);
		}
	}

	free_arr(&DI);
	free_arr(&g);
	free_arr(&x1);
}

// невязка g = Ax - f

// RelRes - относительная невязка = ||g|| / ||f||
void Resid(int n1, int n2, int n3, double *D, int ldd, double *B, double *x_sol, double *f, double *g, double &RelRes)
{
	int n = n1 * n2;
	int size = n * n3;
	double *f1 = alloc_arr(size);
	double done = 1.0;
	int ione = 1;

	Mult_Au(n1, n2, n3, D, ldd, B, x_sol, f1);

#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < size; i++)
		g[i] = f[i] - f1[i];

#ifdef DEBUG
	print_vec(size, f, g, "f and g");
#endif

	RelRes = dlange("Frob", &size, &ione, g, &size, NULL);
	RelRes = RelRes / dlange("Frob", &size, &ione, f, &size, NULL);

	free_arr(&f1);

}

/* Функция вычисления разложения симметричной блочно-диагональной матрицы с использование сжатого формата. 
   Внедиагональные блоки предполагаются диагональными матрицами */
void DirFactFastDiag(int n1, int n2, int n3, double *D, int ldd, double *B, double *G /*factorized matrix*/, int ldg, 
									 double eps, int smallsize, char *bench)
{
	int n = n1 * n2;
	int nbr = n3; // size of D is equal to nbr blocks by n elements
	int size = n * nbr;
	double *TD1 = alloc_arr(n * n); int ldtd = n;
	double *TD = alloc_arr(n * n);

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("Timing DirFactFastDiag\n");
		printf("****************************\n");
	}

	double tt = omp_get_wtime();
	SymRecCompress(n, &D[ind(0, n)], ldd, smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;

	if (compare_str(7, bench, "display")) printf("Compressing D(0) time: %lf\n", tt);

	tt = omp_get_wtime();
	SymCompRecInv(n, &D[ind(0, n)], ldd, &G[ind(0, n)], ldg, smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display")) printf("Computing G(1) time: %lf\n", tt);


	for (int k = 1; k < nbr; k++)
	{
		tt = omp_get_wtime();
		SymRecCompress(n, &D[ind(k, n)], ldd, smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Compressing D(%d) time: %lf\n", k, tt);

		tt = omp_get_wtime();
		dlacpy("All", &n, &n, &G[ind(k - 1, n)], &ldg, TD1, &ldtd);
		DiagMult(n, TD1, ldtd, &B[ind(k - 1, n)], smallsize);
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Mult D(%d) time: %lf\n", k, tt);

		tt = omp_get_wtime();
		Add(n, 1.0, &D[ind(k, n)], ldd, -1.0, TD1, ldtd, TD, ldtd, smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Add %d time: %lf\n", k, tt);

		tt = omp_get_wtime();
		SymCompRecInv(n, TD, ldtd, &G[ind(k,n) + ldg * 0], ldg, smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Computing G(%d) time: %lf\n", k, tt);
		if (compare_str(7, bench, "display")) printf("\n");
	}

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("End of DirFactFastDiag\n");
		printf("****************************\n");
	}

	free_arr(&TD);
	free_arr(&TD1);
}

void DirSolveFastDiag(int n1, int n2, int n3, double *G, int ldg, double *B, double *f, double *x, double eps, int smallsize)
{
	int n = n1 * n2;
	int nbr = n3;
	int size = n * nbr;
	double *tb = alloc_arr(size);
	double *y = alloc_arr(n);

#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < n; i++)
		tb[i] = f[i];

	for (int k = 1; k < nbr; k++)
	{
		RecMultL(n, 1, &G[ind(k - 1, n) + ldg * 0], ldg, &tb[ind(k - 1, n)], size, y, n, smallsize);	
		DenseDiagMult(n, &B[ind(k - 1, n)], y, y);

#pragma omp parallel for simd schedule(simd:static)
		for (int i = 0; i < n; i++)
			tb[ind(k, n) + i] = f[ind(k, n) + i] - y[i];

	}

	RecMultL(n, 1, &G[ind(nbr - 1, n) + ldg * 0], ldg, &tb[ind(nbr - 1, n)], size, &x[ind(nbr - 1, n)], size, smallsize);

	for (int k = nbr - 2; k >= 0; k--)
	{
		DenseDiagMult(n, &B[ind(k, n)], &x[ind(k + 1, n)], y);

#pragma omp parallel for simd schedule(simd:static)
		for (int i = 0; i < n; i++)
			y[i] = tb[ind(k, n) + i] - y[i];

		RecMultL(n, 1, &G[ind(k, n) + ldg * 0], ldg, y, n, &x[ind(k, n)], size, smallsize);
	}

	free_arr(&tb);
	free_arr(&y);
}


// Рекурсивная функция вычисления DAD, где D - диагональная матрица, а A - сжатая
void DiagMult(int n, double *A, int lda, double *d, int small_size)
{

	if (n <= small_size)     // error 4 - не копировалась матрица в этом случае
	{
#pragma omp parallel for simd schedule(simd:static)
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
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n1; j++)
			for (int i = 0; i < n2; i++)
				A[i + n1 + lda * (0 + j)] *= d[n1 + i]; // вторая часть массива D

		// VT * D - каждый j-ый столбец умножается на элемент вектора d[j]
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n2; j++)
			for (int i = 0; i < n1; i++)
				A[i + 0 + lda * (n1 + j)] *= d[j]; 
		// так так вектора матрицы V из разложения A = U * V лежат в транспонированном порядке,
		// то матрицу D стоит умножать на VT слева
	}
}

double rel_error(int n, int k, double *Hrec, double *Hinit, int ldh, double eps)
{
	double norm = 0;

	// Norm of residual
#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < k; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Hrec[i + ldh * j] = Hrec[i + ldh * j] - Hinit[i + ldh * j];
		}
	}

	norm = dlange("Frob", &n, &k, Hrec, &ldh, NULL);
	norm = norm / dlange("Frob", &n, &k, Hinit, &ldh, NULL);

	return norm;
	
	//if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
	//else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);
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
		double *Y12 = alloc_arr(n1 * m); int ldy12 = n1;
		double *Y21 = alloc_arr(n2 * m); int ldy21 = n2;
		double *Y11 = alloc_arr(n1 * m); int ldy11 = n1;
		double *Y22 = alloc_arr(n2 * m); int ldy22 = n2;
		double *inter1 = alloc_arr(n2 * n1); // column major - lda = column
		double *inter2 = alloc_arr(n1 * n2);
	
		// A21 = A21 * A12 (the result of multiplication is A21 matrix with size n2 x n1)
		dgemm("No", "No", &n2, &n1, &n1, &alpha, &A[n1 + lda * 0], &lda, &A[0 + lda * n1], &lda, &beta, inter1, &n2);
		
		// Y21 = inter1 (n2 x n1) * X(1...n1, :) (n1 x n)
		dgemm("No", "No", &n2, &m, &n1, &alpha, inter1, &n2, &X[0 + 0 * ldx], &ldx, &beta, Y21, &ldy21);
	
		// A12 = A21*T = A12*T * A21*T (the result of multiplication is A21 matrix with size n1 x n2)
		dgemm("Trans", "Trans", &n1, &n2, &n1, &alpha, &A[0 + lda * n1], &lda, &A[n1 + lda * 0], &lda, &beta, inter2, &n1);
	
		// Y12 = inter2 (n1 x n2) * X(n1...m, :) (n2 x n)
		dgemm("No", "No", &n1, &m, &n2, &alpha, inter2, &n1, &X[n1 + 0 * ldx], &ldx, &beta, Y12, &ldy12); // уже транспонировали матрицу в предыдущем dgemm

		RecMultL(n1, m, &A[0 + lda * 0], lda, &X[0 + ldx * 0], ldx, Y11, ldy11, smallsize);
		RecMultL(n2, m, &A[n1 + lda * n1], lda, &X[n1 + ldx * 0], ldx, Y22, ldy22, smallsize);

		// first part of Y = Y11 + Y12
		mkl_domatadd('C', 'N', 'N', n1, m, 1.0, Y11, ldy11, 1.0, Y12, ldy12, &Y[0 + ldy * 0], ldy);
		// op_mat(n1, m, Y11, Y12, n1, '+');
		// dlacpy("All", &n1, &m, Y11, &n1, &Y[0 + ldy * 0], &ldy);

		// second part of Y = Y21 + Y22
		mkl_domatadd('C', 'N', 'N', n2, m, 1.0, Y21, ldy21, 1.0, Y22, ldy22, &Y[n1 + ldy * 0], ldy);
		// op_mat(n2, m, Y21, Y22, n2, '+');
		// dlacpy("All", &n2, &m, Y21, &n2, &Y[n1 + ldy * 0], &ldy);

		free_arr(&Y11);
		free_arr(&Y12);
		free_arr(&Y21);
		free_arr(&Y22);
		free_arr(&inter1);
		free_arr(&inter2);
	
	}
}

// Функция вычисления линейной комбинации двух сжатых матриц
void Add(int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc, int smallsize, double eps, char *method)
{
	double alpha_loc = 1.0;
	double beta_loc = 0.0;

	if (fabs(alpha) < eps)
	{
		dlacpy("All", &n, &n, B, &ldb, C, &ldc);
		return;
	}
	else if (fabs(beta) < eps)
	{
		dlacpy("All", &n, &n, A, &lda, C, &ldc);
		return;
	}
	else if (fabs(alpha) < eps && fabs(beta) < eps)
	{
		return;
	}

#ifdef DEBUG
	printf("******Function: Add*******\n");
#endif
	// n - order of A, B and C
	if (n <= smallsize)
	{
#ifdef DEBUG
		printf("Smallsize - doing dense addition for n = %d and smallsize = %d\n", n, smallsize);
#endif
		mkl_domatadd('C', 'N', 'N', n, n, alpha, A, lda, beta, B, ldb, C, ldc);
		//Add_dense(n, n, alpha, A, lda, beta, B, ldb, C, ldc);
	}
	else
	{	
		int p1 = 0, p2 = 0;
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

		mkl_dimatcopy('C', 'N', n2, n1, alpha, &A[n1 + lda * 0], lda, lda);
		mkl_dimatcopy('C', 'N', n2, n1, beta, &B[n1 + ldb * 0], ldb, ldb);
		//Add_dense(n2, n1, alpha, &A[n1 + lda * 0], lda, 0.0, B, ldb, &A[n1 + lda * 0], lda);
		//Add_dense(n2, n1, beta, &B[n1 + ldb * 0], ldb, 0.0, B, ldb, &B[n1 + ldb * 0], ldb);

		// Y21 = [alpha*A{2,1} beta*B{2,1}];
		dlacpy("All", &n2, &n1, &A[n1 + lda * 0], &lda, &Y21[0 + ldy21 * 0], &ldy21);
		dlacpy("All", &n2, &n1, &B[n1 + ldb * 0], &ldb, &Y21[0 + ldy21 * n1], &ldy21);

		// Y12 = [A{1,2}; B{1,2}];
		dlacpy("All", &n1, &n2, &A[0 + lda * n1], &lda, &Y12[0 + ldy12 * 0], &ldy12);
		dlacpy("All", &n1, &n2, &B[0 + ldb * n1], &ldb, &Y12[n1 + ldy12 * 0], &ldy12);

		// произведение Y21 и Y12 - это матрица n2 x n1
		LowRankApprox(n2, n1_dbl, Y21, ldy21, V21, ldv21, p1, eps, "SVD"); // перезапись Y21
		LowRankApprox(n1_dbl, n1, Y12, ldy12, V12, ldv12, p2, eps, "SVD");  // перезапись Y12

		// Y = V21'*V12;
		dgemm("No", "No", &p1, &p2, &n1_dbl, &alpha_loc, V21, &ldv21, Y12, &ldy12, &beta_loc, Y, &ldy); // mn, mn

		// C{2,1} = U21*Y;   
		dgemm("No", "No", &n2, &n1, &p1, &alpha_loc, Y21, &ldy21, Y, &ldy, &beta_loc, &C[n1 + ldc * 0], &ldc); // mn

		// C{1,2} = U12';
		dlacpy("All", &p2, &n1, V12, &ldv12, &C[0 + ldc * n1], &ldc); // n1, n2

		Add(n1, alpha, &A[0 + lda * 0], lda, beta, &B[0 + ldb * 0], ldb, &C[0 + ldc * 0], ldc, smallsize, eps, method);
		Add(n2, alpha, &A[n1 + lda * n1], lda, beta, &B[n1 + ldb * n1], ldb, &C[n1 + ldc * n1], ldc, smallsize, eps, method);

		free_arr(&Y21);
		free_arr(&Y12);
		free_arr(&V21);
		free_arr(&V12);
		free_arr(&Y);
	}

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

	if (fabs(alpha) <= eps)
	{
		dlacpy("All", &n, &n, A, &lda, B, &ldb);
		return;
	}

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

		//mkl_domatcopy('C', 'T', 1.0, n1, k, &V[0 + ldv * 0], ldv, V_uptr, ldvuptr);
		Mat_Trans(n1, k, &V[0 + ldv * 0], ldv, V_uptr, ldvuptr);
	
		// Y12 = [A{1,2} V(1:n1,:)];
		dlacpy("All", &n1, &n1, &A[0 + lda * n1], &lda, &Y12[0 + ldy12 * 0], &ldy12);
		dlacpy("All", &k, &n1, V_uptr, &ldvuptr, &Y12[n1 + ldy21 * 0], &ldy12);

		// [U21,V21] = LowRankApprox (Y21, eps, method);
		LowRankApprox(n2, nk, Y21, ldy21, V21, ldv21, p, eps, "SVD");

		// [U12, V12] = LowRankApprox(Y12, eps, method);
		LowRankApprox(nk, n1, Y12, ldy12, V12, ldv12, p, eps, "SVD");


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


// Recursuve inversion of compressed matrix A to the compressed matrix B
void SymCompRecInv(int n, double *A, int lda, double *B, int ldb, int smallsize, double eps, char *method)
{
	double alpha_one = 1.0;
	double alpha_mone = -1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;
	double wquery = 0;
	int info = 0;
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

		free_arr(&work);
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
		mkl_dimatcopy('C', 'N', n1, n1, -1.0, B12, ldb12, ldb12);
		//Add_dense(n1, n1, alpha_mone, B12, ldb12, beta_zero, B, ldb, B12, ldb12);

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


void Eye(int n, double *H, int ldh)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			if (j == i) H[i + ldh * j] = 1.0;
			else H[i + ldh * j] = 0.0;
}

void Diag(int n, double *H, int ldh, double value)
{
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			if (j == i) H[i + ldh * j] = value;
			else H[i + ldh * j] = 0.0;
}


void Hilbert(int n, double *H, int ldh)
{
#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			H[i + ldh * j] = 1.0 / (i + j + 1);
}

void Mat_Trans(int m, int n, double *H, int ldh, double *Hcomp_tr, int ldhtr)
{
#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
			Hcomp_tr[j + ldhtr * i] = H[i + ldh * j];
}

void Add_dense(int m, int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc)
{
	double dzero = 0.0;

	if (beta == dzero)
	{
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < m; i++)
				C[i + ldc * j] = alpha * A[i + lda * j];
		}
	}
	else if (alpha == dzero)
	{
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < m; i++)
				C[i + ldc * j] = beta * B[i + ldb * j];
		}
	}
	else
	{
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < m; i++)
				C[i + ldc * j] = alpha * A[i + lda * j] + beta * B[i + ldb * j];
		}
	}
}

void op_mat(int m, int n, double *Y11, double *Y12, int ldy, char sign)
{
	if (sign == '+')
	{
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n; j++)
			for (int i = 0; i < m; i++)
				Y11[i + ldy * j] += Y12[i + ldy *j];
	}
	else if (sign == '-')
	{
#pragma omp parallel for simd schedule(simd:static)
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

map<vector<int>, double> dense_to_CSR(int m, int n, double *A, int lda, int *ia, int *ja, double *values)
{
	map<vector<int>, double> CSR;
	vector<int> v(2, 0);
	int k = 0;
	int ik = 0;
	int first_elem_in_row = 0;
	for (int i = 0; i < m; i++)
	{
		first_elem_in_row = 0;
		for (int j = 0; j < n; j++)
		{
			if (fabs(A[i + lda * j]) != 0)
			{
				values[k] = A[i + lda * j];
				if (first_elem_in_row == 0)
				{
					ia[ik] = k + 1;
					ik++;
					first_elem_in_row = 1;
				}
				ja[k] = j + 1;

				v[0] = ia[ik - 1];
				v[1] = ja[k];
				CSR[v] = values[k];

				k++;
			}
		}
	}

	return CSR;
}

void count_dense_elements(int m, int n, double *A, int lda, int& non_zeros)
{
	int k = 0;
#pragma omp parallel for schedule(guided) reduction(+:k)
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (fabs(A[i + lda * j]) != 0)
			{
				k++;
			}
		}
	}
	non_zeros = k;
}

map<vector<int>, double> concat_maps(const map<vector<int>, double>& map1, const map<vector<int>, double>& map2)
{
	map<vector<int>, double> map_res;
	for (const auto& item : map1)
	{
		map_res.insert(item);
	}
	for (const auto& item : map2)
	{
		map_res.insert(item);
	}
	return map_res;
}

void shift_values(map<vector<int>, double>& CSR, int rows, int *ia, int shift_non_zeros, int non_zeros, int *ja, int shift_columns)
{
#pragma omp parallel for schedule(simd:static)
	for (int i = 0; i < rows; i++)
		ia[i] += shift_non_zeros;

#pragma omp parallel for schedule(simd:static)
	for (int i = 0; i < non_zeros; i++)
		ja[i] += shift_columns;
	
}

void construct_block_row(int m, int n, double* BL, int ldbl, double *A, int lda, double *BR, int ldbr, double* Arow, int ldar)
{
	if (BL == NULL)
	{
		dlacpy("All", &m, &n, A, &lda, &Arow[0 + ldar * 0], &ldar);
		dlacpy("All", &m, &n, BR, &ldbr, &Arow[0 + ldar * n], &ldar);
	}
	else if (BR == NULL)
	{
		dlacpy("All", &m, &n, BL, &ldbl, &Arow[0 + ldar * 0], &ldar);
		dlacpy("All", &m, &n, A, &lda, &Arow[0 + ldar * n], &ldar);
	}
	else
	{
		dlacpy("All", &m, &n, BL, &ldbl, &Arow[0 + ldar * 0], &ldar);
		dlacpy("All", &m, &n, A, &lda, &Arow[0 + ldar * n], &ldar);
		dlacpy("All", &m, &n, BR, &ldbr, &Arow[0 + ldar * 2 * n], &ldar);
	}
}

void GenSparseMatrix(size_m x, size_m y, size_m z, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr)
{
	int n = x.n * y.n;
	int size = n * z.n;

	for (int blk = 0; blk < z.n - 1; blk++)
	{
		Diag(n, &BL[ind(blk, n)], ldbl, 1.0 / (z.h * z.h));
		Diag(n, &BR[ind(blk, n)], ldbr, 1.0 / (z.h * z.h));
	}

	map<vector<int>, double> CSR;
	CSR = block3diag_to_CSR(x.n, y.n, z.n, BL, ldbl, A, lda, BR, ldbr, Acsr);
}

void GenSparseMatrix2(size_m x, size_m y, size_m z, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr)
{
	int n = x.n * y.n;
	int size = n * z.n;
	int non_zeros_on_prev_level = 0;
	map<vector<int>, double> CSR;

	Diag(n, BL, ldbl, 1.0 / (z.h * z.h));
	Diag(n, BR, ldbr, 1.0 / (z.h * z.h));

	for (int blk = 0; blk < z.n; blk++)
	{
		GenerateDiagonal2DBlock(blk, x, y, z, A, lda);
		CSR = BlockRowMat_to_CSR(blk, x.n, y.n, z.n, BL, ldbl, A, lda, BR, ldbr, Acsr, non_zeros_on_prev_level); // ВL, ВR and A - is 2D dimensional matrices (n x n)
	}
}

map<vector<int>, double> BlockRowMat_to_CSR(int blk, int n1, int n2, int n3, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr, int& non_zeros_on_prev_level)
{
	map<vector<int>, double> CSR_A;
	map<vector<int>, double> CSR;
	vector<int> v(2, 0);
	int n = n1 * n2;
	int k = 0;
	double *Arow = alloc_arr(n * 3 * n); int ldar = n;

	if (blk == 0)
	{
		construct_block_row(n, n, NULL, ldbl, A, lda, BR, ldbr, Arow, ldar);
		CSR_A = dense_to_CSR(n, 2 * n, Arow, ldar, &Acsr->ia[0], &Acsr->ja[0], &Acsr->values[0]);
		non_zeros_on_prev_level = CSR_A.size();
	}
	else if (blk == n3 - 1)
	{
		construct_block_row(n, n, BL, ldbl, A, lda, NULL, ldbr, Arow, ldar);
		CSR_A = dense_to_CSR(n, 2 * n, Arow, ldar, &Acsr->ia[ind(blk, n)], &Acsr->ja[non_zeros_on_prev_level], &Acsr->values[non_zeros_on_prev_level]);
		shift_values(CSR_A, n, &Acsr->ia[ind(blk, n)], non_zeros_on_prev_level, CSR_A.size(), &Acsr->ja[non_zeros_on_prev_level], n * (blk - 1));
	}
	else
	{
		construct_block_row(n, n, BL, ldbl, A, lda, BR, ldbr, Arow, ldar);
		CSR_A = dense_to_CSR(n, 3 * n, Arow, ldar, &Acsr->ia[ind(blk, n)], &Acsr->ja[non_zeros_on_prev_level], &Acsr->values[non_zeros_on_prev_level]);

		// shift values of arrays according to previous level
		shift_values(CSR_A, n, &Acsr->ia[ind(blk, n)], non_zeros_on_prev_level, CSR_A.size(), &Acsr->ja[non_zeros_on_prev_level], n * (blk - 1));
		non_zeros_on_prev_level += CSR_A.size();
	}

	free(Arow);
	return CSR;
}

map<vector<int>, double> block3diag_to_CSR(int n1, int n2, int blocks, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr)
{
	map<vector<int>, double> CSR_A;
	map<vector<int>, double> CSR;
	vector<int> v(2, 0);
	int n = n1 * n2;
	int k = 0;
	double *AR = alloc_arr(n * 3 * n); int ldar = n;
	int non_zeros_on_prev_level = 0;

	for (int blk = 0; blk < blocks; blk++)
	{
		if (blk == 0)
		{
			construct_block_row(n, n, NULL, ldbl, &A[0], lda, &BR[0], ldbr, AR, ldar);
		//	print(n, n, &AR[0 + ldar * n], ldar, "AR");
			CSR_A = dense_to_CSR(n, 2 * n, AR, ldar, &Acsr->ia[0], &Acsr->ja[0], &Acsr->values[0]);
			non_zeros_on_prev_level = CSR_A.size();
		}
		else if (blk == blocks - 1)
		{
			construct_block_row(n, n, &BL[ind(blk - 1, n)], ldbl, &A[ind(blk, n)], lda, NULL, ldbr, AR, ldar);
			//print(n, 2 * n, AR, ldar, "ldar");
			CSR_A = dense_to_CSR(n, 2 * n, AR, ldar, &Acsr->ia[ind(blk, n)], &Acsr->ja[non_zeros_on_prev_level], &Acsr->values[non_zeros_on_prev_level]);
			shift_values(CSR_A, n, &Acsr->ia[ind(blk, n)], non_zeros_on_prev_level, CSR_A.size(), &Acsr->ja[non_zeros_on_prev_level], n * (blk - 1));
		}
		else
		{
			construct_block_row(n, n, &BL[ind(blk - 1, n)], ldbl, &A[ind(blk, n)], lda, &BR[ind(blk, n)], ldbr, AR, ldar);
			CSR_A = dense_to_CSR(n, 3 * n, AR, ldar, &Acsr->ia[ind(blk, n)], &Acsr->ja[non_zeros_on_prev_level], &Acsr->values[non_zeros_on_prev_level]);

			// shift values of arrays according to previous level
			shift_values(CSR_A, n, &Acsr->ia[ind(blk, n)], non_zeros_on_prev_level, CSR_A.size(), &Acsr->ja[non_zeros_on_prev_level], n * (blk - 1));
			non_zeros_on_prev_level += CSR_A.size();
		}
	}

	free(AR);
	return CSR;
}

void print_map(const map<vector<int>, double>& SD)
{
	cout << "SD size = " << SD.size() << endl;
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

void print_vec(int size, int *vec1, double *vec2, char *name)
{
	printf("%s\n", name);
	for (int i = 0; i < size; i++)
		printf("%d   %d   %lf\n", i, vec1[i], vec2[i]);
}

int compare_str(int n, char *s1, char *s2)
{
	for (int i = 0; i < n; i++)
	{
		if (s1[i] != s2[i]) return 0;
	}
	return 1;
}

void print(int m, int n, double *u, int ldu, char *mess)
{
	printf("%s\n", mess);
	for (int i = 0; i < m; i++)
	{
		printf("%d ", i);
		for (int j = 0; j < n; j++)
		{
			printf("%5.2lf ", u[i + ldu*j]);
		}
		printf("\n");
	}

	printf("\n");

}

void print_vec_mat(int m, int n, double *u, int ldu, double *vec, char *mess)
{
	printf("%s\n", mess);
	for (int i = 0; i < m; i++)
	{
		printf("%d ", i);
		for (int j = 0; j < n; j++)
		{
			printf("%5.2lf ", u[i + ldu*j]);
		}
		printf("  %lf\n", vec[i]);
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