#include "Header.h"
#include "templates.h"

bool lookup(node *node, int value)
{
	if (node == NULL)
	{
		return false;
	}
	else
	{
		if (node->val == value) return true;
		else if (node->val < value) return lookup(node->left, value);
		else return lookup(node->right, value);
	}
}

node* AllocNewNode1(int value)
{
	node* Node = (node*)malloc(sizeof(node));
	Node->val = value;
	Node->left = NULL;
	Node->right = NULL;

	return Node; // return a pointer to the allocated memory inside this function
}

void AllocNewNode2(node* Node, int value)
{
	Node->val = value;
	Node->left = NULL;
	Node->right = NULL;
}


/*
Give a binary search tree and a number, inserts a new node
with the given number in the correct place in the tree.
Returns the new root pointer which the caller should
then use (the standard trick to avoid using reference
parameters).
*/

node* insert1(node* root, int value)
{

	if (root == NULL)
	{
		// этот указатель на память существует внутри функции локально, поэтому его нужно вернуть
		root = AllocNewNode1(value); // положили память в root
		//cout << "root val1 inside = " << root->val << endl;
	}
	else
	{
		if (value <= root->val) root->left = insert1(root->left, value);
		else root->right = insert1(root->right, value);

	}

	return root; // возвращаем память на root (перемещаем память, аллоцированную внутри функции за пределы функции)
}

void insert2(node* *root, int value)
{
	if (*root == NULL)
	{
		// этот указатель на root передан нам по указателю, можем изменять его внутри функции
		*root = AllocNewNode1(value); // положили память в root
		//cout << "root val2 inside = " << (*root)->val << endl;
	}
	else
	{
		if (value <= (*root)->val) insert2(&(*root)->left, value);
		else insert2(&(*root)->right, value);

	}
}

// This problem demonstrates simple binary tree traversal.Given a binary tree, count the number of nodes in the tree.
// recursive function
int TreeSize(node* root)
{
	int size = 0;

	if (root == NULL)
	{
		return 0;
	}
	else
	{
		return TreeSize(root->left) + 1 + TreeSize(root->right);
	}

}

int MaxDepth(node* root)
{
	if (root == NULL)
	{
		return 0;
	}
	else
	{
		int ld = MaxDepth(root->left);
		int rd = MaxDepth(root->right);

		if (ld > rd) return (ld + 1); // + root node
		else return (rd + 1);
	}

}

/* Given a non - empty binary search tree(an ordered binary tree), return the minimum data value found in that tree.
It is not necessary to search the entire tree. A maxValue() function is structurally very similar to this
function. */

int MinValue(node* root)
{
	// we need to find left most leaf
#if 0
	if (root->left == NULL)
	{
		return root->val;
	}
	else
	{
		return MinValue(root->left);
	}
#else
	node* cur = root; // put a pointer to the head of the tree
	while (cur->left != NULL) /* goes throught the tree*/
	{
		cur = cur->left;
	}
	return cur->val;
#endif
}

int MaxValue(node* root)
{
	// we need to find left most leaf
#if 0
	if (root->right == NULL)
	{
		return root->val;
	}
	else
	{
		return MaxValue(root->right);
	}
#else
	node* cur = root; // put a pointer to the head of the tree
	while (cur->right != NULL) /* goes throught the tree*/
	{
		cur = cur->right;
	}
	return cur->val;
#endif
}


/* in*/
void PrintInorder(node* root)
{
	if (root == NULL)
	{
		return;
	}
	else
	{
		PrintInorder(root->left);
		printf("%d ", root->val);
		PrintInorder(root->right);
	}

}

void PrintPostorder(node *root)
{
	if (root == NULL)
	{
		return;
	}
	else
	{
		PrintPostorder(root->left);
		PrintPostorder(root->right);
		printf("%d ", root->val);
	}
}

bool IsBST(node *root)
{
	if (root == NULL)
	{
		return true;
	}

	if (root->left != NULL && MinValue(root->left) > root->val)
		return false;

	if (root->right != NULL && MaxValue(root->right) <= root->val)
		return false;

	if (!IsBST(root->left) || !IsBST(root->right))
		return false;

	return true;

}
#if 0
// ---------- Compressed matrices --------------
mnode* AllocNewMatrixLeaf(int n1, int n2)
{
	mnode* Node = (mnode*)malloc(sizeof(mnode)); // allocated memory for pointers only
	Node->A11 = (double*)malloc(sizeof(n1 * n1));
	Node->A22 = (double*)malloc(sizeof(n2 * n2));
	Node->U = NULL;
	Node->VT = NULL;
	Node->left = NULL;
	Node->right = NULL;

	return Node;
}

mnode* AllocNewMatrixNode(int n1, int n2)
{
	mnode* Node = (mnode*)malloc(sizeof(mnode)); // allocated memory for pointers only
	Node->U = (double*)malloc(sizeof(n2 * n1));
	Node->VT = (double*)malloc(sizeof(n1 * n2));
	Node->A11 = NULL;
	Node->A22 = NULL;
	Node->left = NULL;
	Node->right = NULL;

	return Node;
}
#endif

// ---------------HSS technology----------

void LowRankApproxStruct(int n2, int n1 /* size of A21 = A */,
	double *A /* A is overwritten by U */, int lda, mnode* &Astr, double eps, char *method)
{
	int mn = min(n1, n2);
	int info = 0;
	int lwork = -1;
	Astr->p = 0;

	double wkopt;

	if (compare_str(3, method, "SVD"))
	{
		// mem_alloc
		double *VT = alloc_arr(n1 * n1); int ldvt = n1;
		double *S = alloc_arr(mn);
		dgesvd("Over", "Sing", &n2, &n1, A, &lda, S, VT, &ldvt, VT, &ldvt, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		double *work = alloc_arr(lwork);

		// A = U1 * S * V1
		dgesvd("Over", "Sing", &n2, &n1, A, &lda, S, VT, &ldvt, VT, &ldvt, work, &lwork, &info);
																					

		for (int j = 0; j < mn; j++)
		{
			double s1 = S[j] / S[0];
			if (s1 < eps)
			{
				break;
			}
			Astr->p = j + 1;
			for (int i = 0; i < n2; i++)
				A[i + lda * j] *= S[j];
		}

		// Alloc new node

		Astr->U = alloc_arr(n2 * Astr->p);
		Astr->VT = alloc_arr(Astr->p * n1);
		dlacpy("All", &n2, &Astr->p, A, &lda, Astr->U, &n2);
		dlacpy("All", &Astr->p, &n1, VT, &ldvt, Astr->VT, &Astr->p);
	

#ifdef DEBUG
		printf("LowRankStructure function after SVD: n2 = %d, n1 = %d, p = %d\n", n2, n1, Astr->p);
		print(n2, Astr->p, Astr->U, n2, "U");
		print(Astr->p, n1, Astr->VT, Astr->p, "VT");

#endif
		free_arr(&VT);
		free_arr(&work);
		free_arr(&S);
	}
	return;
}

mnode* LowRankApproxStruct2(int n2, int n1 /* size of A21 = A */,
	double *A /* A is overwritten by U */, int lda, double eps, char *method)
{
	int mn = min(n1, n2);
	int info = 0;
	int lwork = -1;
	mnode* Astr = (mnode*)malloc(sizeof(mnode));
	Astr->p = 0;

	double wkopt;

	if (compare_str(3, method, "SVD"))
	{
		// mem_alloc
		double *VT = alloc_arr(n1 * n1); int ldvt = n1;
		double *S = alloc_arr(mn);
		dgesvd("Over", "Sing", &n2, &n1, A, &lda, S, VT, &ldvt, VT, &ldvt, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		double *work = alloc_arr(lwork);

		// A = U1 * S * V1
		dgesvd("Over", "Sing", &n2, &n1, A, &lda, S, VT, &ldvt, VT, &ldvt, work, &lwork, &info);


		for (int j = 0; j < mn; j++)
		{
			double s1 = S[j] / S[0];
			if (s1 < eps)
			{
				break;
			}
			Astr->p = j + 1;
			for (int i = 0; i < n2; i++)
				A[i + lda * j] *= S[j];
		}

		// Alloc new node

		Astr->U = alloc_arr(n2 * Astr->p);
		Astr->VT = alloc_arr(Astr->p * n1);
		dlacpy("All", &n2, &Astr->p, A, &lda, Astr->U, &n2);
		dlacpy("All", &Astr->p, &n1, VT, &ldvt, Astr->VT, &Astr->p);

		free_arr(&VT);
		free_arr(&work);
		free_arr(&S);
	}
	return Astr;
}

void Test_LowRankApproxStruct(int m, int n, double eps, char *method)
{
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


	double alpha = 1.0;
	double beta = 0.0;


#if 0
	mnode *Astr = (mnode*)malloc(sizeof(mnode));
	printf("Test for LowRankApproximationStruct m = %d n = %d ", m, n);
	LowRankApproxStruct(m, n, A, lda, Astr, eps, "SVD"); // memory allocation for Astr inside function
#else
	mnode *Astr;
	printf("Test for LowRankApproximationStruct2 using return m = %d n = %d ", m, n);
	Astr = LowRankApproxStruct2(m, n, A, lda, eps, "SVD"); // memory allocation for Astr inside function
#endif
	printf("p = %d ", Astr->p);

	dgemm("no", "no", &m, &n, &Astr->p, &alpha, Astr->U, &m, Astr->VT, &Astr->p, &beta, A_rec, &lda);


	rel_error(m, n, A_rec, A_init, lda, eps);
	free_arr(&A);
	free_arr(&A_init);
	free_arr(&A_rec);
	free_arr(&Astr->U);
	free_arr(&Astr->VT);
	free(Astr);
}

void SymRecCompressStruct(int n /* order of A */, double *A /* init matrix */, const int lda, 
	/*output*/ mnode* &ACstr, 
	const int small_size, double eps,
	char *method /* SVD or other */)
{
	ACstr = (mnode*)malloc(sizeof(mnode)); // на каждом шаге мы должны создавать новую структуру, нельзя выносить за функцию этот вызов

	if (n <= small_size)
	{
		alloc_dense_node(n, ACstr);
		dlacpy("All", &n, &n, A, &lda, ACstr->A, &n);
	}
	else
	{
		int n1, n2; // error 3  - неправильное выделение подматриц - похоже на проблему 2
		n2 = (int)ceil(n / 2.0); // округление в большую сторону
		n1 = n - n2; // n2 > n1

		// LowRank A21
		LowRankApproxStruct(n2, n1, &A[n1 + lda * 0], lda, ACstr, eps, method);

#ifdef DEBUG
		printf("SymRecCompressStruct: n = %d n1 = %d n2 = %d p = %d\n", n, n1, n2, ACstr->p);
		print(n1, n1, &A[0 + lda * 0], lda, "Astr");
		print(n2, n2, &A[n1 + lda * n1], lda, "Astr");
#endif

		SymRecCompressStruct(n1, &A[0 + lda * 0], lda, ACstr->left, small_size, eps, method);
		SymRecCompressStruct(n2, &A[n1 + lda * n1], lda, ACstr->right, small_size, eps, method);
	}

}

void SymResRestoreStruct(int n, mnode* H1str, double *H2 /* recovered */, int ldh, int small_size)
{
	double alpha = 1.0;
	double beta = 0.0;

	if (n <= small_size)     // error 4 - не копировалась матрица в этом случае
	{
		dlacpy("All", &n, &n, H1str->A, &n, H2, &ldh);
	}
	else
	{
		int n1, n2;
		n2 = (int)ceil(n / 2.0); // округление в большую сторону
		n1 = n - n2;

		// A21 = A21 * A12
		dgemm("Notrans", "Notrans", &n2, &n1, &H1str->p, &alpha, H1str->U, &n2, H1str->VT, &H1str->p, &beta, &H2[n1 + ldh * 0], &ldh);

		// A12 = A21*T = A12*T * A21*T
		dgemm("Trans", "Trans", &n1, &n2, &H1str->p, &alpha, H1str->VT, &H1str->p, H1str->U, &n2, &beta, &H2[0 + ldh * n1], &ldh);


		SymResRestoreStruct(n1, H1str->left, &H2[0 + ldh * 0], ldh, small_size);
		SymResRestoreStruct(n2, H1str->right, &H2[n1 + ldh * n1], ldh, small_size);
	}
}

void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize)
{
	printf("*****Test for SymRecCompressStruct  n = %d eps = %e ******* ", n, eps);
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
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			H2[i + ldh * j] = H2[i + ldh * j] - H[i + ldh * j];

		}
	}

#ifdef DEBUG
	print(n, n, H, ldh, "H init");
	print(n, n, H2, ldh, "diff");
#endif

	norm = dlange(&frob, &n, &n, H2, &ldh, NULL);
	norm = norm / dlange(&frob, &n, &n, H, &ldh, NULL);
	if (norm < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);

	free_arr(&H);
	free_arr(&H2);
	free_arr(&H1);
}

/* Рекурсивная функция вычисления DAD, где D - диагональная матрица, а Astr - сжатая в структуре */
void DiagMultStruct(int n, mnode* Astr, double *d, int small_size)
{

	if (n <= small_size)     // error 4 - не копировалась матрица в этом случае
	{
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n; j++)
			for (int i = 0; i < n; i++)
			{
				Astr->A[i + j * n] *= d[j]; // справа D - каждый j - ый столбец A умножается на d[j]
				Astr->A[i + j * n] *= d[i]; // слева D - каждая строка A умножается на d[j]
			}
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // округление в большую сторону
		int n1 = n - n2;

		DiagMultStruct(n1, Astr->left, &d[0], small_size);
		DiagMultStruct(n2, Astr->right, &d[n1], small_size);

		// D * U - каждая i-ая строка U умножается на элемент вектора d[i]
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < Astr->p; j++)
			for (int i = 0; i < n2; i++)
				Astr->U[i + n2 * j] *= d[n1 + i]; // вторая часть массива D

														// VT * D - каждый j-ый столбец умножается на элемент вектора d[j]
#pragma omp parallel for simd schedule(simd:static)
		for (int j = 0; j < n1; j++)
			for (int i = 0; i < Astr->p; i++)
				Astr->VT[i + Astr->p * j] *= d[j];
		// так так вектора матрицы V из разложения A = U * V лежат в транспонированном порядке,
		// то матрицу D стоит умножать на VT слева
	}
}

void Test_DiagMultStruct(int n, double eps, char *method, int smallsize)
{
	printf("*****Test for DiagMultStruct  n = %d ******* ", n);
	double *Hd = alloc_arr(n*n); // diagonal Hd = D * H * D
	double *H1 = alloc_arr(n*n); // compressed H
	double *H2 = alloc_arr(n*n); // recovered H after D * H1 * D
	double *d = alloc_arr(n);
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

	mnode *HCstr = (mnode*)malloc(sizeof(mnode));
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
	rel_error(n, n, H2, Hd, ldh, eps);

	free_arr(&Hd); // diagonal Hd = D * H * D
	free_arr(&H1); // compressed H
	free_arr(&H2); // recovered H after D * H1 * D
	free_arr(&d);

}

/* Y = A * X, where A - compressed n * n, X - dense n * m, Y - dense n * m */
void RecMultLStruct(int n, int m, mnode* Astr, double *X, int ldx, double *Y, int ldy, int smallsize)
{
	double alpha = 1.0;
	double beta = 0.0;
	if (n <= smallsize)
	{
		dgemm("No", "No", &n, &m, &n, &alpha, Astr->A, &n, X, &ldx, &beta, Y, &ldy);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // rounding up
		int n1 = n - n2;
		double *Y12 = alloc_arr(n1 * m); int ldy12 = n1;
		double *Y21 = alloc_arr(n2 * m); int ldy21 = n2;
		double *Y11 = alloc_arr(n1 * m); int ldy11 = n1;
		double *Y22 = alloc_arr(n2 * m); int ldy22 = n2;
		double *inter1 = alloc_arr(n2 * n1); // column major - lda = column
		double *inter2 = alloc_arr(n1 * n2);

		// A21 = A21 * A12 (the result of multiplication is A21 matrix with size n2 x n1)
		dgemm("No", "No", &n2, &n1, &Astr->p, &alpha, Astr->U, &n2, Astr->VT, &Astr->p, &beta, inter1, &n2);

		// Y21 = inter1 (n2 x n1) * X(1...n1, :) (n1 x n)
		dgemm("No", "No", &n2, &m, &n1, &alpha, inter1, &n2, &X[0 + 0 * ldx], &ldx, &beta, Y21, &ldy21);

		// A12 = A21*T = A12*T * A21*T (the result of multiplication is A21 matrix with size n1 x n2)
		dgemm("Trans", "Trans", &n1, &n2, &Astr->p, &alpha, Astr->VT, &Astr->p, Astr->U, &n2, &beta, inter2, &n1);

		// Y12 = inter2 (n1 x n2) * X(n1...m, :) (n2 x n)
		dgemm("No", "No", &n1, &m, &n2, &alpha, inter2, &n1, &X[n1 + 0 * ldx], &ldx, &beta, Y12, &ldy12); // we have already transposed this matrix in previous dgemm

		RecMultLStruct(n1, m, Astr->left, &X[0 + ldx * 0], ldx, Y11, ldy11, smallsize);
		RecMultLStruct(n2, m, Astr->right, &X[n1 + ldx * 0], ldx, Y22, ldy22, smallsize);

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

/* Тест на сравнение результатов умножения Y = H * X сжимаемой матрицы H на произвольную X.
Сравниваются результаты со сжатием и без */
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize)
{
	printf("*****Test for RecMultLStruct  n = %d k = %d ******* ", n, k);
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
	mnode *Hstr;
	// Compress H
	SymRecCompressStruct(n, H, ldh, Hstr, smallsize, eps, method);

	// RecMult Y1 = comp(H) * X
	RecMultLStruct(n, k, Hstr, X, ldx, Y1, ldy, smallsize);

	rel_error(n, k, Y1, Y, ldy, eps);

#ifdef DEBUG
	print(n, n, H, ldy, "H comp");
	print(n, k, Y1, ldy, "Y1 rec");
#endif
}

// Функция вычисления линейной комбинации двух сжатых матриц
void AddStruct(int n, double alpha, mnode* Astr, double beta, mnode* Bstr, mnode* &Cstr, int smallsize, double eps, char *method)
{
	double alpha_loc = 1.0;
	double beta_loc = 0.0;
	
	Cstr = (mnode*)malloc(sizeof(mnode));

	// n - order of A, B and C
	if (n <= smallsize)
	{
		alloc_dense_node(n, Cstr);
		mkl_domatadd('C', 'N', 'N', n, n, alpha, Astr->A, n, beta, Bstr->A, n, Cstr->A, n);
		//Add_dense(n, n, alpha, A, lda, beta, B, ldb, C, ldc);
	}
	else
	{
		int p1 = 0, p2 = 0;
		int n2 = (int)ceil(n / 2.0); // округление в большую сторону
		int n1 = n - n2;

		int n1_dbl = Astr->p + Bstr->p;

		double *Y21 = alloc_arr(n2 * n1_dbl); int ldy21 = n2;
		double *Y12 = alloc_arr(n1_dbl * n1); int ldy12 = n1_dbl;

		double *V21 = alloc_arr(n2 * n1_dbl);
		int ldv21 = n2;

		double *V12 = alloc_arr(n1_dbl * n1);
		int ldv12 = n1_dbl;

		mkl_dimatcopy('C', 'N', n2, Astr->p, alpha, Astr->U, n2, n2);
		mkl_dimatcopy('C', 'N', n2, Bstr->p, beta, Bstr->U, n2, n2);
		//Add_dense(n2, n1, alpha, &A[n1 + lda * 0], lda, 0.0, B, ldb, &A[n1 + lda * 0], lda);
		//Add_dense(n2, n1, beta, &B[n1 + ldb * 0], ldb, 0.0, B, ldb, &B[n1 + ldb * 0], ldb);

		// Y21 = [alpha*A{2,1} beta*B{2,1}];
		dlacpy("All", &n2, &Astr->p, Astr->U, &n2, &Y21[0 + ldy21 * 0], &ldy21);
		dlacpy("All", &n2, &Bstr->p, Bstr->U, &n2, &Y21[0 + ldy21 * Astr->p], &ldy21);

		// Y12 = [A{1,2}; B{1,2}];
		dlacpy("All", &Astr->p, &n1, Astr->VT, &Astr->p, &Y12[0 + ldy12 * 0], &ldy12);
		dlacpy("All", &Bstr->p, &n1, Bstr->VT, &Bstr->p, &Y12[Astr->p + ldy12 * 0], &ldy12);

		// произведение Y21 и Y12 - это матрица n2 x n1
		LowRankApprox(n2, n1_dbl, Y21, ldy21, V21, ldv21, p1, eps, "SVD"); // перезапись Y21
		LowRankApprox(n1_dbl, n1, Y12, ldy12, V12, ldv12, p2, eps, "SVD");  // перезапись Y12


		// Y = V21'*V12;
		double *Y = alloc_arr(p1 * p2);
		dgemm("No", "No", &p1, &p2, &n1_dbl, &alpha_loc, V21, &ldv21, Y12, &ldy12, &beta_loc, Y, &p1); // mn, mn

		// C{2,1} = U21*Y;   
		Cstr->U = (double*)malloc(n2 * p2 * sizeof(double));
		dgemm("No", "No", &n2, &p2, &p1, &alpha_loc, Y21, &ldy21, Y, &p1, &beta_loc, Cstr->U, &n2); // mn

		// C{1,2} = U12';
		Cstr->VT = (double*)malloc(p2 * n1 * sizeof(double));
		dlacpy("All", &p2, &n1, V12, &ldv12, Cstr->VT, &p2); // n1, n2
		Cstr->p = p2;

		AddStruct(n1, alpha, Astr->left, beta, Bstr->left, Cstr->left, smallsize, eps, method);
		AddStruct(n2, alpha, Astr->right, beta, Bstr->right, Cstr->right, smallsize, eps, method);

		free_arr(&Y21);
		free_arr(&Y12);
		free_arr(&V21);
		free_arr(&V12);
		free_arr(&Y);
	}

}

void Test_AddStruct(int n, double alpha, double beta, int smallsize, double eps, char *method)
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
	rel_error(n, n, GcR, G, ldg, eps);

	free_arr(&H1);
	free_arr(&H2);
	free_arr(&G);
	free_arr(&H1c);
	free_arr(&H2c);
	free_arr(&Gc);
	free_arr(&GcR);
}

void alloc_dense_node(int n, mnode* &Cstr)
{
	Cstr->A = alloc_arr(n * n);
	Cstr->p = -1;
	Cstr->U = NULL;
	Cstr->VT = NULL;
	Cstr->left = NULL;
	Cstr->right = NULL;
}


/* Функция вычисления симметричного малорангового дополнения A:= A + alpha * V * Y * V'
A - симметрическая сжатая (n x n)
Y - плотная симметричная размера k x k, k << n , V - плотная прямоугольная n x k
(n x n) = (n x n) + (n x k) * (k x k) * (k * n) */
void SymCompUpdate2Struct(int n, int k, mnode* Astr, double alpha, double *Y, int ldy, double *V, int ldv, mnode* &Bstr, int smallsize, double eps, char* method)
{
	double alpha_one = 1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;

	Bstr = (mnode*)malloc(sizeof(mnode));

	int p1 = 0, p2 = 0;
	if (n <= smallsize)
	{
		// X = X + alpha * V * Y * VT

		// C = V * Y
		double *C = alloc_arr(n * k); int ldc = n;
		dsymm("Right", "Up", &n, &k, &alpha_one, Y, &ldy, V, &ldv, &beta_zero, C, &ldc);

		// X = X + alpha * C * Vt
		dgemm("No", "Trans", &n, &n, &k, &alpha, C, &ldc, V, &ldv, &beta_one, Astr->A, &n);

		// B = A
		alloc_dense_node(n, Bstr);
		dlacpy("All", &n, &n, Astr->A, &n, Bstr->A, &n);

		free_arr(&C);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // n2 > n1
		int n1 = n - n2;

		int nk = Astr->p + k;
		// for this division n2 > n1 we can store a low memory
		double *Y12 = alloc_arr(nk * n1); int ldy12 = nk;
		double *Y21 = alloc_arr(n2 * nk); int ldy21 = n2;

		double *V_uptr = alloc_arr(k * n1); int ldvuptr = k;
		double *VY = alloc_arr(n2 * k); int ldvy = n2;

		double *V12 = alloc_arr(nk * n1); int ldv12 = nk;
		double *V21 = alloc_arr(n2 * nk); int ldv21 = n2;

		dgemm("No", "No", &n2, &k, &k, &alpha, &V[n1 + ldv * 0], &ldv, Y, &ldy, &beta_zero, VY, &ldvy);

		// Y21 = [A{2,1} alpha*V(m:n,:)*Y];
		dlacpy("All", &n2, &Astr->p, Astr->U, &n2, &Y21[0 + ldy21 * 0], &ldy21);
		dlacpy("All", &n2, &k, VY, &ldvy, &Y21[0 + ldy21 * Astr->p], &ldy21);

		//mkl_domatcopy('C', 'T', 1.0, n1, k, &V[0 + ldv * 0], ldv, V_uptr, ldvuptr);
		Mat_Trans(n1, k, &V[0 + ldv * 0], ldv, V_uptr, ldvuptr);

		// Y12 = [A{1,2} V(1:n1,:)];
		dlacpy("All", &Astr->p, &n1, Astr->VT, &Astr->p, &Y12[0 + ldy12 * 0], &ldy12);
		dlacpy("All", &k, &n1, V_uptr, &ldvuptr, &Y12[Astr->p + ldy21 * 0], &ldy12);

		// [U21,V21] = LowRankApprox (Y21, eps, method);
		LowRankApprox(n2, nk, Y21, ldy21, V21, ldv21, p1, eps, "SVD");

		// [U12, V12] = LowRankApprox(Y12, eps, method);
		LowRankApprox(nk, n1, Y12, ldy12, V12, ldv12, p2, eps, "SVD");
		Bstr->p = p2;

		// V21 * Y12
		double *VV = alloc_arr(p1 * p2); int ldvv = p1;
		dgemm("No", "No", &p1, &p2, &nk, &alpha_one, V21, &ldv21, Y12, &ldy12, &beta_zero, VV, &ldvv);

		// B{2,1} = U21*(V21'*V12);
		Bstr->U = alloc_arr(n2 * p2);
		dgemm("No", "No", &n2, &p2, &p1, &alpha_one, Y21, &ldy21, VV, &ldvv, &beta_zero, Bstr->U, &n2);

		// B{1,2} = U12;
		Bstr->VT = alloc_arr(p2 * n1);
		dlacpy("All", &p2, &n1, V12, &ldv12, Bstr->VT, &p2);

		// B{1,1} = SymCompUpdate2 (A{1,1}, Y, V(1:n1,:), alpha, eps, method);
		SymCompUpdate2Struct(n1, k, Astr->left, alpha, Y, ldy, &V[0 + ldv * 0], ldv, Bstr->left, smallsize, eps, method);

		// B{2,2} = SymCompUpdate2 (A{2,2}, Y, V(m:n ,:), alpha, eps, method);
		SymCompUpdate2Struct(n2, k, Astr->right, alpha, Y, ldy, &V[n1 + ldv * 0], ldv, Bstr->right, smallsize, eps, method);

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
void Test_SymCompUpdate2Struct(int n, int k, double alpha, int smallsize, double eps, char* method)
{
	printf("*****Test for SymCompUpdate2Struct  n = %d k = %d ***** ", n, k);
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
	rel_error(n, n, B_rec, H, ldh, eps);

	free_arr(&B);
	free_arr(&B_rec);
	free_arr(&H);
	free_arr(&HC);
	free_arr(&Y);
	free_arr(&C);
	free_arr(&V);
}

void SymCompRecInvStruct(int n, mnode* Astr, mnode* &Bstr, int smallsize, double eps, char *method)
{
	double alpha_one = 1.0;
	double alpha_mone = -1.0;
	double beta_zero = 0.0;
	double beta_one = 1.0;
	int info = 0;
	double wquery = 0;
	int lwork = -1;

	Bstr = (mnode*)malloc(sizeof(mnode));

	if (n <= smallsize)
	{
		int *ipiv = (int*)malloc(n * sizeof(int));

		// LU factorization of A
		dgetrf(&n, &n, Astr->A, &n, ipiv, &info);

		// space query
		dgetri(&n, Astr->A, &n, ipiv, &wquery, &lwork, &info);

		lwork = (int)wquery;
		double *work = alloc_arr(lwork);

		// inversion of A
		dgetri(&n, Astr->A, &n, ipiv, work, &lwork, &info);

		// dlacpy
		alloc_dense_node(n, Bstr);
		dlacpy("All", &n, &n, Astr->A, &n, Bstr->A, &n);

		free_arr(&work);
		free(ipiv);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // n2 > n1
		int n1 = n - n2;

		Bstr->p = Astr->p;
		double *X11 = alloc_arr(n1 * n1); int ldx11 = n1;
		double *X22 = alloc_arr(n2 * n2); int ldx22 = n2;
		double *V = alloc_arr(n1 * Astr->p); int ldv = n1;
		double *B12 = alloc_arr(n1 * Bstr->p); int ldb12 = n1;
		double *Y = alloc_arr(Astr->p * Bstr->p); int ldy = Astr->p;
		mnode *X11str, *X22str;

		// Inversion of A22 to X22
		SymCompRecInvStruct(n2, Astr->right, X22str, smallsize, eps, method);

		// Save X22 * U to B{2,1}
		Bstr->U = alloc_arr(n2 * Bstr->p);
		RecMultLStruct(n2, Bstr->p, X22str, Astr->U, n2, Bstr->U, n2, smallsize);

		// Compute Y = UT * X22 * U = | A[2,1]T * B{2,1} | = | (p x n2) x (n2 x p)  = (p x p) |
		dgemm("Trans", "No", &Astr->p, &Bstr->p, &n2, &alpha_one, Astr->U, &n2, Bstr->U, &n2, &beta_zero, Y, &ldy);

		// Update X11 = A11 - V * UT * X22 * U * VT = | A11 - V * Y * VT | = | (n1 x n1) - (n1 x p) * (p x p) * (p x n1) |
		Mat_Trans(Astr->p, n1, Astr->VT, Astr->p, V, ldv);
		SymCompUpdate2Struct(n1, Astr->p, Astr->left, alpha_mone, Y, ldy, V, ldv, X11str, smallsize, eps, method);

		// Inversion of X11 to B11
		SymCompRecInvStruct(n1, X11str, Bstr->left, smallsize, eps, method);

		// Fill B{1,2} as B12 = -B{1,1} * A{1,2} = -X11 * V = (n1 x n1) * (n1 x p) = (n1 x p)
		RecMultLStruct(n1, Bstr->p, Bstr->left, V, ldv, B12, ldb12, smallsize);
		mkl_dimatcopy('C', 'N', n1, Bstr->p, -1.0, B12, ldb12, ldb12);

		// B{1,2} = transpose(B12)
		Bstr->VT = alloc_arr(Bstr->p *n1);
		Mat_Trans(n1, Bstr->p, B12, ldb12, Bstr->VT, Bstr->p);

		// Y = -(A{1,2})' * B{1,2} = -VT * (-X11 * V) = - VT * B12 = (p x n1) * (n1 x p)
		dgemm("No", "No", &Astr->p, &Bstr->p, &n1, &alpha_mone, Astr->VT, &Astr->p, B12, &ldb12, &beta_zero, Y, &ldy);

		// Update X22 + (X22*U) * VT * X11 * V (UT * X22) = X22 + B21 * Y * B21T = (n2 x n2) + (n2 x p) * (p x p) * (p x n2)
		SymCompUpdate2Struct(n2, Bstr->p, X22str, alpha_one, Y, ldy, Bstr->U, n2, Bstr->right, smallsize, eps, method);

		free_arr(&X11);
		free_arr(&X22);
		free_arr(&Y);
		free_arr(&V);
		free_arr(&B12);
	}
}

void Test_SymCompRecInvStruct(int n, int smallsize, double eps, char *method)
{
	printf("***** Test_SymCompRecInvStruct n = %d eps = %lf **** ", n, eps);
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

	mnode *HCstr, *BCstr;
	SymRecCompressStruct(n, Hc, ldh, HCstr, smallsize, eps, method);
	SymCompRecInvStruct(n, HCstr, BCstr, smallsize, eps, method);
	SymResRestoreStruct(n, BCstr, Brec, ldb, smallsize);

	Eye(n, Y, ldy);

	// Y = Y - H * Brec
	dgemm("No", "No", &n, &n, &n, &alpha_mone, H, &ldh, Brec, &ldb, &beta_one, Y, &ldy);

	double norm = dlange("Frob", &n, &n, Y, &ldy, NULL);

	if (norm < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);

}

#if 1
/* Функция вычисления разложения симметричной блочно-диагональной матрицы с использование сжатого формата.
Внедиагональные блоки предполагаются диагональными матрицами */
void DirFactFastDiagStruct(int n1, int n2, int n3, double *D, int ldd, double *B, mnode** &Gstr,
	double eps, int smallsize, char *bench)
{
	int n = n1 * n2;
	int nbr = n3; // size of D is equal to nbr blocks by n elements
	int size = n * nbr;

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("Timing DirFactFastDiag\n");
		printf("****************************\n");
	}

	double tt = omp_get_wtime();
	mnode* *DCstr = (mnode**)malloc(n3 * sizeof(mnode*));
	SymRecCompressStruct(n, &D[ind(0, n)], ldd, DCstr[0], smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;

	if (compare_str(7, bench, "display")) printf("Compressing D(0) time: %lf\n", tt);

	tt = omp_get_wtime();

	Gstr = (mnode**)malloc(n3 * sizeof(mnode*));
	SymCompRecInvStruct(n, DCstr[0], Gstr[0], smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display")) printf("Computing G(1) time: %lf\n", tt);


	for (int k = 1; k < nbr; k++)
	{
		mnode *TDstr, *TD1str;
		tt = omp_get_wtime();
		SymRecCompressStruct(n, &D[ind(k, n)], ldd, DCstr[k], smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Compressing D(%d) time: %lf\n", k, tt);
	
		tt = omp_get_wtime();
		CopyStruct(n, Gstr[k - 1], TD1str, smallsize);
	
		DiagMultStruct(n, TD1str, &B[ind(k - 1, n)], smallsize);
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Mult D(%d) time: %lf\n", k, tt);

		tt = omp_get_wtime();
		AddStruct(n, 1.0, DCstr[k], -1.0, TD1str, TDstr, smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Add %d time: %lf\n", k, tt);

		tt = omp_get_wtime();
		SymCompRecInvStruct(n, TDstr, Gstr[k], smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Computing G(%d) time: %lf\n", k, tt);
		if (compare_str(7, bench, "display")) printf("\n");

		FreeNodes(n, TDstr, smallsize);
		FreeNodes(n, TD1str, smallsize);
	}

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("End of DirFactFastDiag\n");
		printf("****************************\n");
	}

	
	for (int i = n3 - 1; i >= 0; i--)
	{
		FreeNodes(n, DCstr[i], smallsize);
	}

	free(DCstr);

}
#if 1
/* Функция вычисления разложения симметричной блочно-диагональной матрицы с использование сжатого формата.
Внедиагональные блоки предполагаются диагональными матрицами */
void DirFactFastDiagStructOnline(size_m x, size_m y, size_m z, mnode** &Gstr, double *B,
	double eps, int smallsize, char *bench)
{
	int n = x. n * y.n;
	int nbr = z.n; // size of D is equal to nbr blocks by n elements
	int size = n * nbr;
	double *DD = alloc_arr(n * n); int lddd = n;

	// Set vector B
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < z.n - 1; j++)
		for (int i = 0; i < n; i++)
			B[ind(j, n) + i] = 1.0 / (z.h * z.h);

	

	if (compare_str(7, bench, "display"))
	{
		printf("**********************************\n");
		printf("Timing DirFactFastDiagStructOnline\n");
		printf("**********************************\n");
	}

	double tt = omp_get_wtime();
	mnode* *DCstr = (mnode**)malloc(z.n * sizeof(mnode*));
	GenerateDiagonal2DBlock(0, x, y, z, DD, lddd);
	printf("constructrd\n");
	SymRecCompressStruct(n, DD, lddd, DCstr[0], smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;

	if (compare_str(7, bench, "display")) printf("Compressing D(0) time: %lf\n", tt);

	tt = omp_get_wtime();

	Gstr = (mnode**)malloc(z.n * sizeof(mnode*));
	SymCompRecInvStruct(n, DCstr[0], Gstr[0], smallsize, eps, "SVD");
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display")) printf("Computing G(1) time: %lf\n", tt);


	for (int k = 1; k < nbr; k++)
	{
		double *DD = alloc_arr(n * n); int lddd = n;
		mnode *TDstr, *TD1str;
		tt = omp_get_wtime();
		GenerateDiagonal2DBlock(0, x, y, z, DD, lddd);
		SymRecCompressStruct(n, DD, lddd, DCstr[k], smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Compressing D(%d) time: %lf\n", k, tt);

		tt = omp_get_wtime();
		CopyStruct(n, Gstr[k - 1], TD1str, smallsize);

		DiagMultStruct(n, TD1str, &B[ind(k - 1, n)], smallsize);
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Mult D(%d) time: %lf\n", k, tt);

		tt = omp_get_wtime();
		AddStruct(n, 1.0, DCstr[k], -1.0, TD1str, TDstr, smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Add %d time: %lf\n", k, tt);

		tt = omp_get_wtime();
		SymCompRecInvStruct(n, TDstr, Gstr[k], smallsize, eps, "SVD");
		tt = omp_get_wtime() - tt;
		if (compare_str(7, bench, "display")) printf("Computing G(%d) time: %lf\n", k, tt);
		if (compare_str(7, bench, "display")) printf("\n");

		FreeNodes(n, TDstr, smallsize);
		FreeNodes(n, TD1str, smallsize);
		free(DD);
	}

	if (compare_str(7, bench, "display"))
	{
		printf("****************************\n");
		printf("End of DirFactFastDiag\n");
		printf("****************************\n");
	}


	for (int i = z.n - 1; i >= 0; i--)
	{
		FreeNodes(n, DCstr[i], smallsize);
	}

	free(DCstr);
	free(DD);
	free(B);

}
#endif
void DirSolveFastDiagStruct(int n1, int n2, int n3, mnode* *Gstr, double *B, double *f, double *x, double eps, int smallsize)
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
		RecMultLStruct(n, 1, Gstr[k - 1], &tb[ind(k - 1, n)], size, y, n, smallsize);
		DenseDiagMult(n, &B[ind(k - 1, n)], y, y);

#pragma omp parallel for simd schedule(simd:static)
		for (int i = 0; i < n; i++)
			tb[ind(k, n) + i] = f[ind(k, n) + i] - y[i];

	}

	RecMultLStruct(n, 1, Gstr[nbr - 1], &tb[ind(nbr - 1, n)], size, &x[ind(nbr - 1, n)], size, smallsize);

	for (int k = nbr - 2; k >= 0; k--)
	{
		DenseDiagMult(n, &B[ind(k, n)], &x[ind(k + 1, n)], y);

#pragma omp parallel for simd schedule(simd:static)
		for (int i = 0; i < n; i++)
			y[i] = tb[ind(k, n) + i] - y[i];

		RecMultLStruct(n, 1, Gstr[k], y, n, &x[ind(k, n)], size, smallsize);
	}

	free_arr(&tb);
	free_arr(&y);
}


void CopyStruct(int n, mnode* Gstr, mnode* &TD1str, int smallsize)
{
	TD1str = (mnode*)malloc(sizeof(mnode));
	if (n <= smallsize)
	{
		alloc_dense_node(n, TD1str);
		dlacpy("All", &n, &n, Gstr->A, &n, TD1str->A, &n);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // n2 > n1
		int n1 = n - n2;

		TD1str->p = Gstr->p;
		TD1str->U = alloc_arr(n2 * TD1str->p);
		TD1str->VT = alloc_arr(TD1str->p * n1);
		dlacpy("All", &n2, &TD1str->p, Gstr->U, &n2, TD1str->U, &n2);
		dlacpy("All", &TD1str->p, &n1, Gstr->VT, &TD1str->p, TD1str->VT, &TD1str->p);

		CopyStruct(n1, Gstr->left, TD1str->left, smallsize);
		CopyStruct(n2, Gstr->right, TD1str->right, smallsize);
	}
}

void FreeNodes(int n, mnode* &Astr, int smallsize)
{
	if (n <= smallsize)
	{
		free(Astr->A);
	}
	else
	{
		int n2 = (int)ceil(n / 2.0); // n2 > n1
		int n1 = n - n2;

		free(Astr->U);
		free(Astr->VT);

		FreeNodes(n1, Astr->left, smallsize);
		FreeNodes(n2, Astr->right, smallsize);
	}

	free(Astr);
}

void ResidCSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_sol, double *f, double* g, double &RelRes)
{
	int n = n1 * n2;
	int size = n * n3;
	double *f1 = alloc_arr(size);
	int ione = 1;

	
	mkl_dcsrgemv("No", &size, Dcsr->values, Dcsr->ia, Dcsr->ja, x_sol, f1);

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

void Test_TransferBlock3Diag_to_CSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_orig, double *f, double eps)
{
	int n = n1 * n2;
	int size = n * n3;
	double RelRes = 0;
	double *g = alloc_arr(size);
	ResidCSR(n1, n2, n3, Dcsr, x_orig, f, g, RelRes);

	if (RelRes < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", RelRes, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", RelRes, eps);

	free(g);
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

}

void Block3DSPDSolveFastStruct(size_m x, size_m y, size_m z, double *D, int ldd, double *B, double *f, dcsr* Dcsr, double thresh, int smallsize, int ItRef, char *bench,
	/* output */ 	mnode** &Gstr, double *x_sol, int &success, double &RelRes, int &itcount)
{
	int size = x.n * y.n * z.n;
	int n = x.n * y.n;
	double tt;
	double tt1;
#ifndef CSR_FORMAT
	double *DI = alloc_arr(size * n); int lddi = size;
	double *BI = alloc_arr(size - n); // vector of diagonal elementes
	double *x_orig = alloc_arr(size);
	GenMatrixandRHSandSolution2(x, y, z, DI, ldd, BI, x_orig, f, thresh);
//	dlacpy("All", &size, &n, D, &ldd, DI, &lddi);
#endif

	tt = omp_get_wtime();
#ifndef ONLINE
	DirFactFastDiagStruct(x.n, y.n, z.n, D, ldd, B, Gstr, thresh, smallsize, bench);
#else
	DirFactFastDiagStructOnline(x, y, z, Gstr, B, thresh, smallsize, bench);
#endif
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display"))
	{
		printf("Total factorization time: %lf\n", tt);
	}

	tt = omp_get_wtime();
	DirSolveFastDiagStruct(x.n, y.n, z.n, Gstr, BI, f, x_sol, thresh, smallsize);
	tt = omp_get_wtime() - tt;
	if (compare_str(7, bench, "display"))
	{
		printf("Solving time: %lf\n", tt);
	}

	double *g = alloc_arr(size);
	double *x1 = alloc_arr(size);
	RelRes = 1;
#ifndef CSR_FORMAT
	Resid(x.n, y.n, z.n, DI, lddi, BI, x_sol, f, g, RelRes);
#else
	ResidCSR(x.n, y.n, z.n, Dcsr, x_sol, f, g, RelRes);
#endif

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
				DirSolveFastDiagStruct(x.n, y.n, z.n, Gstr, BI, g, x1, thresh, smallsize);

#pragma omp parallel for simd schedule(simd:static)
				for (int i = 0; i < size; i++)
					x_sol[i] = x_sol[i] + x1[i];

				// начальное решение f сравниваем с решением A_x0 + A_x1 + A_x2
#ifndef CSR_FORMAT
				Resid(x.n, y.n, z.n, DI, lddi, BI, x_sol, f, g, RelRes);
#else
				ResidCSR(x.n, y.n, z.n, Dcsr, x_sol, f, g, RelRes);
#endif
				itcount = itcount + 1;
				tt = omp_get_wtime() - tt;
				if (compare_str(7, bench, "display")) printf("itcount=%d, RelRes=%lf, Time=%lf\n", itcount, RelRes, tt);
			}
			if ((RelRes < thresh) && (itcount < ItRef)) success = 1; // b

			tt1 = omp_get_wtime() - tt1;
			if (compare_str(7, bench, "display")) printf("Iterative refinement total time: %lf\n", tt1);
		}
	}

	for (int i = z.n - 1; i >= 0; i--)
	{
		FreeNodes(n, Gstr[i], smallsize);
	}

	free(Gstr);
#ifndef CSR_FORMAT
	free_arr(&DI);
#endif
	free_arr(&g);
	free_arr(&x1);
}

void Test_CopyStruct(int n, double eps, char *method, int smallsize)
{
	double *H = alloc_arr(n * n); int ldh = n;
	double *H1 = alloc_arr(n * n);
	double *H2 = alloc_arr(n * n);

	printf("***Test CopyStruct n = %d ", n);

	mnode* Hstr, *Hcopy_str;

	Hilbert(n, H, ldh);
	Hilbert(n, H1, ldh);

	SymRecCompressStruct(n, H, ldh, Hstr, smallsize, eps, method);
	CopyStruct(n, Hstr, Hcopy_str, smallsize);
	SymResRestoreStruct(n, Hcopy_str, H2, ldh, smallsize);

	rel_error(n, n, H2, H1, ldh, eps);
}



#endif