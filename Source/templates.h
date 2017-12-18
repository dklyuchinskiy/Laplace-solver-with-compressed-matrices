// HSS

void SymRecCompress(int n /* order of A */, double *A /* init matrix */, const int lda,
	const int small_size, double eps,
	char *method /* SVD or other */);

void LowRankApprox(int n2, int n1 /* size of A21 = A */, double *A /* A is overwritten by U */, int lda,
							double *V /* V is stored in A12 */, int ldv, int&p, double eps, char *method);

void SymResRestore(int n, double *H1, double *H2, int ldh, int p);
void DiagMult(int n, double *A, int lda, double *d, int small_size);
void RecMultL(int n, int m, double *A, int lda, double *X, int ldx, double *Y, int ldy, int smallsize);
void Add(int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc, int smallsize, double eps, char *method);
void SymCompUpdate2(int n, int k, double *A, int lda, double alpha, double *Y, int ldy, double *V, int ldv, double *B, int ldb, int smallsize, double eps, char* method);
void SymCompRecInv(int n, double *A, int lda, double *B, int ldb, int smallsize, double eps, char *method);

// Solver

void GenMatrixandRHSandSolution(const int n1, const int n2, const int n3,
	/* output */ double *D, int ldd, double *B, double *x1, double *f);
void Block3DSPDSolveFast(int n1, int n2, int n3, double *D, int ldd, double *B, double *f, double thresh, int smallsize, int ItRef, char *bench,
	/* output */ double *G, int ldg, double *x, int &success, double &RelRes, int &itcount);

void DirFactFastDiag(int n1, int n2, int n3, double *D, int ldd, double *B, double *G /*factorized matrix*/, int ldg, double eps, int smallsize, char *bench);
void DirSolveFastDiag(int n1, int n2, int n3, double *G, int ldg, double *B, double *f, double *x, double eps, int smallsize);

void GenMatrixandRHSandSolution2(size_m x, size_m y, size_m z,
	/* output */ double *D, int ldd, double *B, double *x1, double *f, double thresh);



// Tests

void Test_SymRecCompress(int n, double eps, char *method, int smallsize);
void Test_DiagMult(int n, double eps, char *method, int smallsize);
void Test_RecMultL(int n, int k, double eps, char *method, int smallsize);
void Test_Add(int n, double alpha, double beta, int smallsize, double eps, char *method);
void Test_LowRankApprox(int m, int n, double eps, char *method);
void Test_SymCompUpdate2(int n, int k, double alpha, int smallsize, double eps, char* method);
void Test_SymCompRecInv(int n, int smallsize, double eps, char *method);
void Test_Transpose(int m, int n, int smallsize, double eps, char *method);

// Tests - BinaryTrees
void Test_LowRankApproxStruct(int m, int n, double eps, char *method);
void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize);
void Test_DiagMultStruct(int n, double eps, char *method, int smallsize);
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize);
void Test_CopyStruct(int n, double eps, char *method, int smallsize);
void Test_AddStruct(int n, double alpha, double beta, int smallsize, double eps, char *method);
void Test_SymCompRecInvStruct(int n, int smallsize, double eps, char *method);
void Test_SymCompUpdate2Struct(int n, int k, double alpha, int smallsize, double eps, char* method);

// Support

double* alloc_arr(int n);
double random(double min, double max);
double F_ex(double x, double y, double z);
double u_ex(double x, double y, double z);

void free_arr(double* *arr);
void Mat_Trans(int m, int n, double *H, int ldh, double *Hcomp_tr, int ldhtr);
void Hilbert(int n, double *H, int ldh);
void op_mat(int n1, int n, double *Y11, double *Y12, int ldy, char sign);
void Add_dense(int m, int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc);
void rel_error(int n, int k, double *Hrec, double *Hinit, int ldh, double eps);
void Resid(int n1, int n2, int n3, double *D, int ldd, double *B, double *x, double *f, double *g, double &RelRes);
void print_map(const map<vector<int>, double>& SD);
void Eye(int n, double *H, int ldh);
void Diag(int n, double *H, int ldh, double value);
void Add_dense_vect(int n, double alpha, double *a, double beta, double *b, double *c);
void GenSolVector(int size, double *x1);
void DenseDiagMult(int n, double *diag, double *v, double *f);
void Mult_Au(int n1, int n2, int n3, double *D, int ldd, double *B, double *u, double *Au /*output*/);
void print(int m, int n, double *u, int ldu, char *mess);
void print_vec_mat(int m, int n, double *u, int ldu, double *vec, char *mess);
void print_vec(int size, double *vec1, double *vec2, char *name);
void print_vec(int size, int *vec1, double *vec2, char *name);

int compare_str(int n, char *s1, char *s2);
int ind(int j, int n);
map<vector<int>, double> dense_to_sparse(int m, int n, double *DD, int ldd, int *i_ind, int *j_ind, double *d);
map<vector<int>, double> dense_to_CSR(int m, int n, double *A, int lda, int *ia, int *ja, double *values);
map<vector<int>, double> block3diag_to_CSR(int n1, int n2, int blocks, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr);
map<vector<int>, double> concat_maps(const map<vector<int>, double>& map1, const map<vector<int>, double>& map2);
void construct_block_row(int m, int n, double* BL, int ldbl, double *A, int lda, double *BR, int ldbr, double* AR, int ldar);

// BinaryTrees.cpp

node* AllocNewNode1(int value);
node* insert1(node* root, int value);
int TreeSize(node* Node);
int MaxDepth(node* Node);
int MinValue(node* root);
int MaxValue(node* root);
void AllocNewNode2(node* Node, int value);
void insert2(node* *root, int value); // передаем дерево по указателю
void PrintInorder(node* root);
void PrintPostorder(node *root);
void FreeNodes(int n, mnode* &Astr, int smallsize);
void alloc_dense_node(int n, mnode* &Cstr);
void CopyStruct(int n, mnode* Gstr, mnode* &TD1str, int smallsize);
bool IsBST(node *root);
bool lookup(node *node, int value);

 // BinaryTrees - Solver

void LowRankApproxStruct(int n2, int n1 /* size of A21 = A */,
	double *A /* A is overwritten by U */, int lda, mnode* &Astr, double eps, char *method);
mnode* LowRankApproxStruct2(int n2, int n1 /* size of A21 = A */,
	double *A /* A is overwritten by U */, int lda, double eps, char *method);
void SymRecCompressStruct(int n, double *A, const int lda, mnode* &ACstr, const int small_size, double eps, char *method);
void SymResRestoreStruct(int n, mnode* H1str, double *H2 /* recovered */, int ldh, int small_size);
void DiagMultStruct(int n, mnode* Astr, double *d, int small_size);
void RecMultLStruct(int n, int m, mnode* Astr, double *X, int ldx, double *Y, int ldy, int smallsize);
void AddStruct(int n, double alpha, mnode* Astr, double beta, mnode* Bstr, mnode* &Cstr, int smallsize, double eps, char *method);
void SymCompUpdate2Struct(int n, int k, mnode* Astr, double alpha, double *Y, int ldy, double *V, int ldv, mnode* &Bstr, int smallsize, double eps, char* method);
void SymCompRecInvStruct(int n, mnode* Astr, mnode* &Bstr, int smallsize, double eps, char *method);
void DirSolveFastDiagStruct(int n1, int n2, int n3, mnode* *Gstr, double *B, double *f, double *x, double eps, int smallsize);
void DirFactFastDiagStruct(int n1, int n2, int n3, double *D, int ldd, double *B, mnode** &Gstr,
	double eps, int smallsize, char *bench);
void DirFactFastDiagStructOnline(size_m x, size_m y, size_m z, mnode** &Gstr, double *B, double thresh, int smallsize, char *bench);
void Block3DSPDSolveFastStruct(size_m x, size_m y, size_m z, double *D, int ldd, double *B, double *f, dcsr* Dcsr, double thresh, int smallsize, int ItRef, char *bench,
	/* output */ 	mnode** &Gstr, double *x_sol, int &success, double &RelRes, int &itcount);
void ResidCSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_sol, double *f, double* g, double &RelRes);
void Test_TransferBlock3Diag_to_CSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_orig, double *f, double eps);
void Test_CompareColumnsOfMatrix(int n1, int n2, int n3, double* D, int ldd, double* B, dcsr* Dcsr, double thresh);
void GenSparseMatrix(size_m x, size_m y, size_m z, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr);
void GenerateDiagonal2DBlock(int part_of_field, size_m x, size_m y, size_m z, double *DD, int lddd);
void GenRHSandSolution(size_m x, size_m y, size_m z, double *u, double *f);



