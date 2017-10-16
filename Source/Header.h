#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <map>
#include <vector>
#include "mkl.h"
using namespace std;

#include <time.h>
#include "templates.h"


#define PROBLEM 2


#define N 10

#if (N <= 5)
#define DEBUG
#endif

#define eps 0.00000001

#define min(a,b) ((a) <= (b)) ? (a) : (b)




