
void Init_matrix(int nbl, int n2, int n3, double *A, int ldA);
void sparse_lu(int nbl, int n2, int n3, double *A, int ldA, int **ipiv_mat);
void print(int m, int n, double *u, int ldu);
void fill_left(int nbl, int j, double h, double *A, int ldA);
void fill_right(int nbl, int j, double h, double *A, int ldA);
void fill_middle(int nbl, int j, double h, double *A, int ldA);
void change(char right, int nbl, double *L, int ldL, int *ipiv);
void swap_col(int nbl, int k, int m, double *L, int ldL, int *ipiv);
void test(int nbl, int n2, int n3, double *A_f, int ldaf, double *U_f, int lduf, double *L_f, int ldlf, int **ipiv_mat);
void construct_A(int nbl, int n2, int n3, double *A, int ldA, double *A_f, int ldaf);
void construct_L(int nbl, int n2, int n3, double *A, int ldA, double *L_f, int ldlf);
void construct_U(int nbl, int n2, int n3, double *A, int ldA, double *U_f, int lduf);

double* alloc_arr(int n);
void free_arr(double* *arr);

// HSS

void SymRecCompress(int n /* order of A */, double *A /* init matrix */, const int lda,
	const int small_size,
	char *method /* SVD or other */);

void LowRankApprox(char *method, int small_size, int n2, int n1 /* size of A21 = A */,
	double *A /* A is overwritten by U */, int lda, double *V /* V is stored in A12 */, int ldv);

void Test_SymRecCompress(int n, double *H, double *H1, double *H2, int ldh);

int compare_str(int n, char *s1, char *s2);

void SymResRestore(int n, double *H1, double *H2, int ldh, int p);