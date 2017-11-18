#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <map>
#include <vector>
#include "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.4.210\windows\mkl\include\mkl.h"
using namespace std;

#include <time.h>

#define PROBLEM 2

#if (PROBLEM == 2)
#define TEST 1
#endif


#define N 10

#if (PROBLEM == 5)
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

#include "templates.h"

#define PI 3.141592653589793238462643






