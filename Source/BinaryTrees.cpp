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
// Low Rank approximation
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
		printf("LowRank after SVD: n2 = %d, n1 = %d, p = %d\n", n2, n1, Astr->p);
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
		ACstr->p = -1;
		ACstr->A = alloc_arr(n * n);
		ACstr->U = NULL;
		ACstr->VT = NULL;
		ACstr->left = NULL;
		ACstr->right = NULL;

		dlacpy("All", &n, &n, A, &lda, ACstr->A, &n);
	}
	else
	{
		int n1, n2; // error 3  - неправильное выделение подматриц - похоже на проблему 2
		n2 = (int)ceil(n / 2.0); // округление в большую сторону
		n1 = n - n2; // n2 > n1

#ifdef DEBUG
		printf("SymRecCompress: n = %d n1 = %d n2 = %d p = %d\n", n, n1, n2, ACstr->p);
		print(n1, n1, &A[0 + lda * 0], lda, "Astr");
		print(n2, n2, &A[n1 + lda * n1], lda, "Astr");
#endif
		// LowRank A21
		LowRankApproxStruct(n2, n1, &A[n1 + lda * 0], lda, ACstr, eps, method);

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


