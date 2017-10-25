
void Init_matrix(int nbl, int n2, int n3, double *A, int ldA);
void sparse_lu(int nbl, int n2, int n3, double *A, int ldA, int **ipiv_mat);
void print(int m, int n, double *u, int ldu, char *mess);
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
	const int small_size, int eps,
	char *method /* SVD or other */);

void LowRankApprox(char *method, int n2, int n1 /* size of A21 = A */,
	double *A /* A is overwritten by U */, int lda, double *V /* V is stored in A12 */, int ldv, int&p);

void Test_SymRecCompress(int n, double *H, double *H1, double *H2, int ldh);

int compare_str(int n, char *s1, char *s2);

void SymResRestore(int n, double *H1, double *H2, int ldh, int p);

// Solver

void GenMatrixandRHSandSolution(const int n1, const int n2, const int n3,
	/* output */ double *D, int ldd, double *B, int ldb, double *x1, double *f);
void Block3DSPDSolveFast(int n1, int n2, int n3, double *D, int ldd, double *B, int ldb, double *f, double thresh, int smallsize, int ItRef, char *bench,
	/* output */ double *G, int ldg, double *x, int &success, double &RelRes, int &itcount);

map<vector<int>,double> dense_to_sparse(int m, int n, double *DD, int ldd, int *i_ind, int *j_ind, double *d);

void print_map(const map<vector<int>,double>& SD);

double random(double min, double max);

void print_vec(int size, double *vec1, double *vec2, char *name);

void DirFactFastDiag(int n1, int n2, int n3, double *D, int ldd, double *B, int ldb, double eps, int smallsize, char *bench, double *G /*factorized matrix*/, int ldg);
void DirSolveFastDiag(int n1, int n2, int n3, double *G, int ldg, double *B, int ldb, double *f, double *x);
void Resid(double *D, int ldd, double *B, int ldb, double *x, double *f, double *g, double &RelRes);
void RecMult(double *G, int ldg, double *tb, double *y);
void SymCompRecInv(double *DC, int lddc, double *G, int ldg, double eps, char *str);

void Test_DiagMult(int n, double eps, char *method, int smallsize);
void rel_error(int n, int k, double *Hrec, double *Hinit, int ldh, double eps);

void RecMultL(int n, int m, double *A, int lda, double *X, int ldx, double *Y, int ldy, int smallsize);
void Test_RecMultL(int n, int k, double eps, char *method, int smallsize);
void op_mat(int n1, int n, double *Y11, double *Y12, int ldy, char sign);

void Add_dense(int m, int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc);
void Test_add(int n, double alpha, double beta, double smallsize, double eps, char *method);
void Test_LowRankApprox_InitA(int m, int n, double eps, char *method);

void Test_transpose(int m, int n, int smallsize, double eps, char *method);
void Mat_Trans(int m, int n, double *H, int ldh, double *Hcomp_tr, int ldhtr);