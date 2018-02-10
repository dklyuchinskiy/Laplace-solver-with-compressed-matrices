#pragma once

typedef void(*ptr_test_low_rank)(int, int, double, char *);
typedef void(*ptr_test_sym_rec_compress)(int, double, char *, int);
typedef void(*ptr_test_mult_diag)(int, int, double, char *, int);
typedef void(*ptr_test_add)(int, double, double, double, char *, int);
typedef void(*ptr_test_update)(int, int, double, double, char*, int);
typedef void(*ptr_test)(int, double, double, int, double, char *);
typedef void(*ptr_test_shell)(ptr_test, const string&, int &, int &);

// Tests
void Test_LowRankApprox(int m, int n, double eps, char *method);
void Test_SymRecCompress(int n, double eps, char *method, int smallsize);
void Test_DiagMult(int n, double eps, char *method, int smallsize);
void Test_RecMultL(int n, int k, double eps, char *method, int smallsize);
void Test_Add(int n, double alpha, double beta, double eps, char *method, int smallsize);
void Test_SymCompUpdate2(int n, int k, double alpha, double eps, char* method, int smallsize);
void Test_SymCompRecInv(int n, double eps, char *method, int smallsize);
void Test_LowRankApproxTranspose(int m, int n, double eps, char *method, int smallsize);

// Tests - BinaryTrees
void Test_LowRankApproxStruct(int m, int n, double eps, char *method);
void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize);
void Test_DiagMultStruct(int n, double eps, char *method, int smallsize);
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize);
void Test_AddStruct(int n, double alpha, double beta, double eps, char *method, int smallsize);
void Test_SymCompUpdate2Struct(int n, int k, double alpha, double eps, char* method, int smallsize);
void Test_SymCompRecInvStruct(int n, double eps, char *method, int smallsize);
void Test_CopyStruct(int n, double eps, char *method, int smallsize);
void Test_QueueList(int n, double eps, char* method, int smallsize);
void Test_TransferBlock3Diag_to_CSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_orig, double *f, double eps);
void Test_CompareColumnsOfMatrix(int n1, int n2, int n3, double* D, int ldd, double* B, dcsr* Dcsr, double thresh);
void Test_DirFactFastDiagStructOnline(size_m x, size_m y, size_m z, mnode** Gstr, double *B, double thresh, int smallsize);
void Test_DirSolveFactDiagStructConvergence(size_m x, size_m y, size_m z, mnode** Gstr, double thresh, int smallsize);
void Test_DirSolveFactDiagStructBlockRanks(size_m x, size_m y, size_m z, mnode** Gstr);

// Tests Shells
void TestAll();
void Shell_LowRankApprox(ptr_test_low_rank func, const string& test_name, int &numb, int &fail_count);
void Shell_SymRecCompress(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);
void Shell_DiagMult(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);
void Shell_RecMultL(ptr_test_mult_diag func, const string& test_name, int &numb, int &fail_count);
void Shell_Add(ptr_test_add func, const string& test_name, int &numb, int &fail_count);
void Shell_SymCompUpdate2(ptr_test_update func, const string& test_name, int &numb, int &fail_count);
void Shell_SymCompRecInv(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);
void Shell_LowRankApproxTranspose(ptr_test_mult_diag func, const string& test_name, int &numb, int &fail_count);
void Shell_CopyStruct(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);

// TestOnline

void Test_RankEqual(mnode *Astr, mnode *AIstr);
void Test_RankAdd(mnode *Astr, mnode *Bstr, mnode* Cstr);