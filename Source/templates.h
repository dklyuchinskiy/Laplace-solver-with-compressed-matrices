
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