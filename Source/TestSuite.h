#pragma once

// Tests
void Test_SymRecCompress(int n, double eps, char *method, int smallsize);
void Test_DiagMult(int n, double eps, char *method, int smallsize);
void Test_RecMultL(int n, int k, double eps, char *method, int smallsize);
void Test_Add(int n, double alpha, double beta, int smallsize, double eps, char *method);
void Test_LowRankApprox(int m, int n, double eps, char *method);
void Test_SymCompUpdate2(int n, int k, double alpha, int smallsize, double eps, char* method);
void Test_SymCompRecInv(int n, int smallsize, double eps, char *method);
void Test_LowRankApproxTranspose(int m, int n, int smallsize, double eps, char *method);

// Tests - BinaryTrees
void Test_LowRankApproxStruct(int m, int n, double eps, char *method);
void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize);
void Test_DiagMultStruct(int n, double eps, char *method, int smallsize);
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize);
void Test_CopyStruct(int n, double eps, char *method, int smallsize);
void Test_AddStruct(int n, double alpha, double beta, int smallsize, double eps, char *method);
void Test_SymCompRecInvStruct(int n, int smallsize, double eps, char *method);
void Test_SymCompUpdate2Struct(int n, int k, double alpha, int smallsize, double eps, char* method);
void Test_QueueList(int n, double eps, char* method, int smallsize);
void Test_TransferBlock3Diag_to_CSR(int n1, int n2, int n3, dcsr* Dcsr, double* x_orig, double *f, double eps);
void Test_CompareColumnsOfMatrix(int n1, int n2, int n3, double* D, int ldd, double* B, dcsr* Dcsr, double thresh);

// Tests Shells
void TestAll();
void Shell_SymRecCompress(int &numb);
void Shell_DiagMult(int &numb);
void Shell_RecMultL(int &numb);
void Shell_LowRankApprox(int &numb);
void Shell_Add(int &numb);
void Shell_SymCompUpdate2(int &numb);
void Shell_SymCompRecInv(int &numb);
void Shell_LowRankApproxTranspose(int &numb);
void Shell_CopyStruct(int &numb);