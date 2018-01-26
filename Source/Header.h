#pragma once

// C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.4.210\windows\mkl\include\mkl.h"

// C++
#include <iostream>
#include <map>
#include <vector>
#include <queue>

using namespace std;

#define PROBLEM 1

#if (PROBLEM == 1)
#define TEST 1
#endif


#define N 10

//#define DEBUG

#if (PROBLEM == 6)
//#define DEBUG
#endif

#define TEST_IT 10

#define EPS 0.00000001

#define min(a,b) ((a) <= (b)) ? (a) : (b)

struct size_m {
	int l;
	int n;
	double h;
};

struct BinaryTreeNode {
	int val;
	struct BinaryTreeNode *left;
	struct BinaryTreeNode *right;
};

typedef struct BinaryTreeNode node;

struct BinaryMatrixTreeNode {

	int p = 0;
	double *U = NULL;
	double *VT = NULL;
	double *A = NULL;
	struct BinaryMatrixTreeNode *left;
	struct BinaryMatrixTreeNode *right;
};

typedef struct BinaryMatrixTreeNode mnode;

struct MatrixCSR {

	int *ia = NULL;
	int *ja = NULL;
	double *values = NULL;
};

typedef struct MatrixCSR dcsr;

struct list {
	mnode* node; // указатель на структуру
	struct list* next;
};

struct my_queue {
	struct list *first, *last;
};

typedef struct list qlist;

#define CSR_FORMAT
#define STRUCT
#define ONLINE
//#define LARGE_SUITE

#define PI 3.141592653589793238462643






