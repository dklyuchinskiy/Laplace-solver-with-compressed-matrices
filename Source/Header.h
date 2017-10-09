#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mkl.h"

#include "templates.h"


#define N 50

#if (N <= 5)
#define DEBUG
#endif

#define eps 0.00000001

