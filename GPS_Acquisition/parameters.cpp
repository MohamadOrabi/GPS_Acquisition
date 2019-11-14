#include "stdafx.h"
#include <fftw3.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <complex>
#include <Eigen/Dense>
#include <cmath>
#include "concurrentqueue.h"
#include "TypeDefs.h"

const fp_precision fss = 5; //in MHz
const int nsampless = 5000;

const fp_precision duration = 2;
const int ntrack = round(duration * (fss * 1E6) / nsampless);