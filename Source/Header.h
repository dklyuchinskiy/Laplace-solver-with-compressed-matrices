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
#include "templates.h"


#define PROBLEM 8


#define N 10

#if (PROBLEM == 6)
//#define DEBUG
#endif

#define TEST_IT 10

#define EPS 0.00000001

#define min(a,b) ((a) <= (b)) ? (a) : (b)




