#pragma once

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

extern const fp_precision fs;	//in MHz
extern const int nsamples;

extern const fp_precision duration;
extern const int ntrack;