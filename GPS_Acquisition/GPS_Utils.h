#pragma once
#include "stdafx.h"
#include "TypeDefs.h"
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

using namespace std;

void conduct_FFT(VEC& in, IQData* out);

void conduct_IFFT(VEC& in, IQData* out);

VEC Corr(VEC a_fft, VEC b_fft); //Takes in the fft of the vectors

void findMaxIndexAndMean(const VEC& data, int* maxIndex, fp_precision* maxValue, fp_precision* meanValue);

RealArray creatRealArray(double start, double step, double end);
vector<int> creatIntVector(double start, double step, double end);

IQArray mixSignal(IQArray signal, const double fs, const double frequency, double* initial_theta);

struct PLL {
	Eigen::Matrix2d Ad = Eigen::Matrix2d::Zero();
	Eigen::Vector2d Bd = Eigen::Vector2d::Zero();
	Eigen::RowVector2d Cd = Eigen::RowVector2d::Zero();

	double Dd = 0.0;
	Eigen::Vector2d xk_PLL;	//!< State vector for PLL.
};

PLL configureLoopFilter(const int loop_order_CFO, const int BL_Target_CFO, const double Tsub, const double fd);

double Update_PLL(PLL* loop_filter, const double ek);

void exportVEC(const string fileName, VEC vec_to_export);
void exportvec(const string fileName, vector<fp_precision> vec_to_export);
void exportArray(const string fileName, Eigen::Array<fp_precision, Eigen::Dynamic, Eigen::Dynamic> vec_to_export);
void exportBits(const string fileName, vector<bool> vec_to_export);


void exportIQ(const string fileName, VEC vec_to_export);

VEC shiftedByRows(const VEC& in, int down);

RealVEC shiftedByRowsReal(const RealVEC& in, int down);

VEC resamplePRN(RealArray orginial_PRN, const int nsamples, const double chip_rate, const double fs, const int code_length, const int starting_index);

vector<bool> subvector(vector<bool>, int starting_index, int length);
vector<bool> getBits(vector<bool> main_vector, int starting_index, int length);
void printBits(vector<bool> bits);
double twosCompliment2Dec(vector<bool> twos_compliment);

double bin2dec(vector<bool> binary);

void ecef2lla(const fp_precision x, const fp_precision y, const fp_precision z, fp_precision* latitude,
	fp_precision* longitude, fp_precision* altitude);

extern fp_precision N_Power;