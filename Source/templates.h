#pragma once

/****************************
Prototypes for all functions.
****************************/

using namespace std;

// functions.cpp
// HSS

void SymRecCompress(int n, double *A, const int lda, const int small_size, double eps, char *method);
void LowRankApprox(int n2, int n1, double *A , int lda, double *V , int ldv, int&p, double eps, char *method);
void SymResRestore(int n, double *H1, double *H2, int ldh, int p);
void DiagMult(int n, double *A, int lda, double *d, int small_size);
void RecMultL(int n, int m, double *A, int lda, double *X, int ldx, double *Y, int ldy, int smallsize);
void Add(int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc, int smallsize, double eps, char *method);
void SymCompUpdate2(int n, int k, double *A, int lda, double alpha, double *Y, int ldy, double *V, int ldv, double *B, int ldb, int smallsize, double eps, char* method);
void SymCompRecInv(int n, double *A, int lda, double *B, int ldb, int smallsize, double eps, char *method);

// Solver

void GenMatrixandRHSandSolution(const int n1, const int n2, const int n3, double *D, int ldd, double *B, double *x1, double *f);
void Block3DSPDSolveFast(int n1, int n2, int n3, double *D, int ldd, double *B, double *f, double thresh, int smallsize, int ItRef, char *bench,
	 double *G, int ldg, double *x, int &success, double &RelRes, int &itcount);
void DirFactFastDiag(int n1, int n2, int n3, double *D, int ldd, double *B, double *G, int ldg, double eps, int smallsize, char *bench);
void DirSolveFastDiag(int n1, int n2, int n3, double *G, int ldg, double *B, double *f, double *x, double eps, int smallsize);
void GenMatrixandRHSandSolution2(size_m x, size_m y, size_m z, double *D, int ldd, double *B, double *x1, double *f, double thresh);


// Support

double* alloc_arr(int n);
double random(double min, double max);
double F_ex(double x, double y, double z);
double u_ex(double x, double y, double z);
double rel_error(int n, int k, double *Hrec, double *Hinit, int ldh, double eps);

void free_arr(double* *arr);
void Mat_Trans(int m, int n, double *H, int ldh, double *Hcomp_tr, int ldhtr);
void Hilbert(int n, double *H, int ldh);
void op_mat(int n1, int n, double *Y11, double *Y12, int ldy, char sign);
void Add_dense(int m, int n, double alpha, double *A, int lda, double beta, double *B, int ldb, double *C, int ldc);
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
map<vector<int>, double> BlockRowMat_to_CSR(int blk, int n1, int n2, int n3, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr, int& non_zeros_on_prev_level);
map<vector<int>, double> concat_maps(const map<vector<int>, double>& map1, const map<vector<int>, double>& map2);
void construct_block_row(int m, int n, double* BL, int ldbl, double *A, int lda, double *BR, int ldbr, double* AR, int ldar);

// BinaryTrees.cpp

int TreeSize(mnode* root);
int MaxDepth(mnode* Node);
void PrintRanks(mnode* root);
void PrintRanksInWidth(mnode *root);
void CopyStruct(int n, mnode* Gstr, mnode* &TD1str, int smallsize);
void FreeNodes(int n, mnode* &Astr, int smallsize);
void alloc_dense_node(int n, mnode* &Cstr);
void PrintStruct(int n, mnode *root);

// BinaryTrees.cpp

// Solver
void LowRankApproxStruct(int n2, int n1, double *A, int lda, mnode* &Astr, double eps, char *method);
mnode* LowRankApproxStruct2(int n2, int n1, double *A, int lda, double eps, char *method);
void SymRecCompressStruct(int n, double *A, const int lda, mnode* &ACstr, const int smallsize, double eps, char *method);
void SymResRestoreStruct(int n, mnode* H1str, double *H2, int ldh, int smallsize);
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
	mnode** &Gstr, double *x_sol, int &success, double &RelRes, int &itcount);
void ResidCSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_sol, double *f, double* g, double &RelRes);
void GenSparseMatrix(size_m x, size_m y, size_m z, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr);
void GenerateDiagonal2DBlock(int part_of_field, size_m x, size_m y, size_m z, double *DD, int lddd);
void GenRHSandSolution(size_m x, size_m y, size_m z, double* B, double *u, double *f);
void GenSparseMatrixOnline(size_m x, size_m y, size_m z, double *BL, int ldbl, double *A, int lda, double *BR, int ldbr, dcsr* Acsr);
void count_dense_elements(int m, int n, double *A, int lda, int& non_zeros);

// Queue
void init(struct my_queue* &q);
bool my_empty(struct my_queue* q);
void push(struct my_queue* &q, mnode* node);
void pop(struct my_queue* &q);
mnode* front(struct my_queue* q);
void PrintRanksInWidthList(mnode *root);
void print_queue(struct my_queue* q);






