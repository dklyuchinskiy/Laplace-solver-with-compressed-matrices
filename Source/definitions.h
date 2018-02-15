#pragma once

// C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.1.156\windows\mkl\include\mkl.h"

// C++
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


#define TEST 1

#define N 10

//#define DEBUG

#define TEST_IT 10

#define EPS 0.00000001

#define min(a,b) ((a) <= (b)) ? (a) : (b)

struct size_m {
	int l;
	int n;
	double h;
};

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
#define FULL_SVD

//#define COL_UPDATE
//#define COL_ADD


#define PI 3.141592653589793238462643






