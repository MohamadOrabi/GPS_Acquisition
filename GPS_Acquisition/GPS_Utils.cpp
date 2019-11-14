#define _USE_MATH_DEFINES

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
#include <iomanip>
#include "GPS_Utils.h"

using namespace std;

void conduct_FFT(VEC& in, IQData* out) {
	fftw_plan plan;	//This is the "plan" this contains all the data the FFTW needs to compute FFT
	auto inn = fftw_alloc_complex(in.size());
	auto outt = fftw_alloc_complex(in.size());
	plan = fftw_plan_dft_1d(in.size(), inn, outt, FFTW_FORWARD, FFTW_ESTIMATE);

	if (in.data() == out) {
		// If in-place transform, copy input into temporary array to prevent FFTW from messing up the transform.
		auto input = fftw_alloc_complex(in.size());
		memcpy(input, in.data(), in.size() * sizeof(RealArray));
		fftw_execute_dft(plan, input, reinterpret_cast<fftComplex*>(out));
		fftw_free(input);
	}
	else {
		fftw_execute_dft(plan, reinterpret_cast<fftComplex*>(in.data()), reinterpret_cast<fftComplex*>(out));
	}
	fftw_free(inn); fftw_free(outt);
}

void conduct_IFFT(VEC& in, IQData* out) {
	fftw_plan plan;	//This is the "plan" this contains all the data the FFTW needs to compute FFT
	auto inn = fftw_alloc_complex(in.size());
	auto outt = fftw_alloc_complex(in.size());
	plan = fftw_plan_dft_1d(in.size(), inn, outt, FFTW_BACKWARD, FFTW_ESTIMATE);

	if (in.data() == out) {
		// If in-place transform, copy input into temporary array to prevent FFTW from messing up the transform.
		auto input = fftw_alloc_complex(in.size());
		memcpy(input, in.data(), in.size() * sizeof(RealArray));
		fftw_execute_dft(plan, input, reinterpret_cast<fftComplex*>(out));
		fftw_free(input);
	}
	else {
		fftw_execute_dft(plan, reinterpret_cast<fftComplex*>(in.data()), reinterpret_cast<fftComplex*>(out));
	}
	fftw_free(inn); fftw_free(outt);
}

VEC Corr(VEC a_fft, VEC b_fft) {	//a is data b is code
	//Conducting cross correlation
	VEC corr = a_fft.array() * b_fft.conjugate().array();
	conduct_IFFT(corr, corr.data());	//Inverse FFT
	//corr /= corr.size();	//Normalizing .. this is an issue with FFTW
	corr /= corr.size();
	return corr;
}

void findMaxIndexAndMean(const VEC& data, int* maxIndex, fp_precision* maxValue, fp_precision* meanValue)
{
	Eigen::Array<fp_precision, Eigen::Dynamic, 1> t = data.abs();
	*maxValue = t.maxCoeff(maxIndex);
	*meanValue = t.mean();
}

RealArray creatRealArray(double start, double step, double end) {
	RealArray out;

	const double size = (end - start) / step + 1;
	out.setLinSpaced(size, start, start + step * (size - 1));
	return out;
}

vector<int> creatIntVector(double start, double step, double end) {

	const double size = (end - start) / step + 1;
	vector<int> out(size);

	for (int i = 0; i < size; i++) {
		out[i] = start + i * step;
	}
	return out;
}

IQArray mixSignal(IQArray signal, const double fs, const double frequency, double* initial_theta) {
	IQArray mixedSignal(signal.size());
	//if (time_array.size()-1 == signal.size()) {
		//double theta = *initial_theta;
		for (int i = 0; i < signal.size(); i++) {
			mixedSignal[i] = signal[i] * IQData(cos(*initial_theta), -sin(*initial_theta));
			*initial_theta += 2 * M_PI * frequency / fs / 1E6;
		}
		*initial_theta = fmod(*initial_theta,2 * M_PI);
		return mixedSignal;
	//}
	//else {
	//	cout << "Error in mixSignal signal and time_array don'th have same size: " << endl << signal.size() << endl << time_array.size();
	//}
}

//Configure Loop filter: PLL Parameters and Initial state
PLL configureLoopFilter(const int loop_order_CFO, const int BL_Target_CFO, const double Tsub, const double fd)
{
	double K;
	double A;
	double B;

	Eigen::Vector2d xk_PLL;	

	Eigen::Matrix2d Ad;
	Eigen::Vector2d Bd;
	Eigen::RowVector2d Cd;
	double Dd;

	switch (loop_order_CFO) {
	case 1:
		K = 4 * BL_Target_CFO;

		Ad << 0, 0,
			0, 0;
		Bd << 0,
			0;
		Cd << 0, 0;
		Dd = K;
		break;
	case 2:
		K = 8.0 / 3.0 * BL_Target_CFO;
		A = K / 2;

		Ad << 1, 0,
			0, 0;
		Bd << Tsub,
			0;
		Cd << K * A, 0;
		Dd = K;
		xk_PLL(0) = (2 * M_PI * fd) / Cd(0);
		break;
	case 3:
		A = 1.2 * BL_Target_CFO;
		B = pow(A, 2) / 2;
		K = 2 * A;

		Ad << 1, Tsub,
			0, 1;
		Bd << pow(Tsub, 2) / 2,
			Tsub;
		Cd << K * B, K* A;
		Dd = K;
		break;
	default:
		throw std::runtime_error("Error: configureLoopFilter: loop order must be one of 1, 2, or 3.");
	}

	return PLL{Ad, Bd, Cd, Dd, xk_PLL};
}

//Update PLL state using discriminator output ek and give vk as output
double Update_PLL(PLL* loop_filter, const double ek)
{	
	const double vk = (*loop_filter).Cd * (*loop_filter).xk_PLL + (*loop_filter).Dd * ek;
	(*loop_filter).xk_PLL = (*loop_filter).Ad * (*loop_filter).xk_PLL + (*loop_filter).Bd * ek;
	
	return vk;
}

void exportVEC(const string fileName, VEC vec_to_export) {

	ofstream outputFile;
	
	outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	outputFile.seekp(outputFile.beg);
	//string legendTitle = fileName;
	//legendTitle.erase(legendTitle.end()-4,legendTitle.end());
	//outputFile << legendTitle << endl;
	
	for (int i = 0; i < vec_to_export.size(); i++) {
		outputFile << std::setprecision(12) << real(vec_to_export[i]) << endl;
	}
	
	outputFile.close();
}

void exportvec(const string fileName, vector<fp_precision> vec_to_export) {

	ofstream outputFile;

	outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	outputFile.seekp(outputFile.beg);
	//string legendTitle = fileName;
	//legendTitle.erase(legendTitle.end()-4,legendTitle.end());
	//outputFile << legendTitle << endl;

	for (int i = 0; i < vec_to_export.size(); i++) {
		outputFile << std::setprecision(12) << vec_to_export[i] << endl;
	}

	outputFile.close();
}

void exportArray(const string fileName, Eigen::Array<fp_precision, Eigen::Dynamic, Eigen::Dynamic> vec_to_export)
{
	ofstream outputFile;
	
	outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	outputFile.seekp(outputFile.beg);
	for (int j = 0; j < vec_to_export.rows(); j++) {
		for (int i = 0; i < vec_to_export.cols(); i++) {
			outputFile << std::setprecision(12) << real(vec_to_export(j,i)) << ",";
		}
		outputFile << endl;
	}

	outputFile.close();
}

void exportBits(const string fileName, vector<bool> vec_to_export) {
	ofstream outputFile;

	outputFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	outputFile.seekp(outputFile.beg);
	//string legendTitle = fileName;
	//legendTitle.erase(legendTitle.end()-4,legendTitle.end());
	//outputFile << legendTitle << endl;

	for (int i = 0; i < vec_to_export.size(); i++) {
		outputFile << vec_to_export[i] << endl;
	}

	outputFile.close();
}

void exportIQ(const string fileName, VEC vec_to_export) {

	ofstream outputFile;

	outputFile.open(fileName, std::ofstream::out | std::ofstream::app);
	outputFile.seekp(outputFile.beg);
	string legendTitle = fileName;
	legendTitle.erase(legendTitle.end() - 4, legendTitle.end());
	outputFile << legendTitle << endl;

	for (int i = 0; i < vec_to_export.size(); i++) {
		outputFile << real(vec_to_export[i]) << "," << imag(vec_to_export[i]) << endl;
	}

	outputFile.close();
}

VEC shiftedByRows(const VEC& in, int down)
{	
	//if (down >= in.size()) { cout << "shiftedByRows: shift is >= " << in.size() << endl; }

	if (down == 0) return in;
	VEC out(in.rows(), in.cols());
	if (down > 0) down = down % in.rows();
	else down = in.rows() - (-down % in.rows());
	// We avoid the implementation-defined sign of modulus with negative arg. 
	int rest = in.rows() - down;
	out.topRows(down) = in.bottomRows(down);
	out.bottomRows(rest) = in.topRows(rest);
	return out;
}

RealVEC shiftedByRowsReal(const RealVEC& in, int down)
{
	if (down == 0) return in;
	RealVEC out(in.size());
	if (down >= 0) down = down % in.size();
	else down = in.size() - (-down % in.size());
	// We avoid the implementation-defined sign of modulus with negative arg. 
	int rest = in.size() - down;
	out.topRows(down) = in.bottomRows(down);
	out.bottomRows(rest) = in.topRows(rest);
	return out;
}

VEC resamplePRN(RealArray orginial_PRN, const int nsamples, const double chip_rate, const double fs, const int code_length, const int starting_index) {
	VEC resampled_PRN(nsamples);
	for (int i = 0; i < nsamples; i++) {
		int index = round(chip_rate / fs * (i + starting_index));
		if (index >= 0) index = index % code_length;
		else index = code_length - (-index % code_length);
		//const int index = ((int)(chip_rate / fs * (i + starting_index))) % code_length; // not + 1 here cz array starts at 0
		resampled_PRN[i] = orginial_PRN[index];
	}
	return resampled_PRN;
}

vector<bool> subvector(vector<bool> main_vector, int starting_index, int length) {//Starting index is referred to by the end of the vector
	vector<bool> sub_vector(main_vector.end() - starting_index, main_vector.end() + length - starting_index);
	return sub_vector;
}

vector<bool> getBits(vector<bool> main_vector, int starting_index, int length) {	//First bit is numbered 1
	vector<bool> bits(main_vector.begin() + starting_index - 1, main_vector.begin() + starting_index - 1 + length);
	return bits;
}

void printBits(vector<bool> bits) {
	for (int i = 0; i < bits.size(); i++) {
		cout << bits[i];
	}
	cout << endl;
}

double twosCompliment2Dec(vector<bool> twos_compliment)
{
	double dec_number = 0;

	int sign;

	if (twos_compliment[0]) { 
		sign = -1;
		bool found_first_1 = false;
		for (int i = twos_compliment.size() - 1; i >= 0; i--) {
			if (!found_first_1) {
				if (twos_compliment[i])	found_first_1 = true;
			}
			else {
				twos_compliment[i] = !twos_compliment[i];
			}
		}

		dec_number = bin2dec(twos_compliment);

		return dec_number * sign;
	}
	else { 
		sign = 1;
		return bin2dec(twos_compliment); 
	}
	

}

double bin2dec(vector<bool> binary) {
	double number = 0;
	for (int i = binary.size(); i > 0; i--) {
		number += (int)binary[binary.size()-i] * pow(2, i-1);
	}
	return number;
}

void ecef2lla(const fp_precision x, const fp_precision y, const fp_precision z, fp_precision* latitude,
	fp_precision* longitude, fp_precision* altitude) {

	fp_precision a = 6378137.0;
	fp_precision e = 0.08181919092890624;
	fp_precision p = sqrt(pow(x, 2) + pow(y, 2));


	fp_precision k_new = 1 / (1 - pow(e, 2));
	fp_precision err_threshold = 0.0001;
	fp_precision err = 1;

	while (err > err_threshold) {
		fp_precision k_old = k_new;
		//ci = (p ^ 2 + (1 - e ^ 2) * z ^ 2 * k_old ^ 2) ^ (3 / 2) / (a * e ^ 2);
		fp_precision ci = pow((pow(p, 2) + (1 - pow(e, 2)) * pow(z, 2) * pow(k_old, 2)), 3.0 / 2.0) / (a * pow(e, 2));
		//k_new = 1+(p^2+(1-e^2)*z^2*k_old^3)/(ci-p^2);

		k_new = 1 + (pow(p, 2) + (1 - pow(e, 2)) * pow(z, 2) * pow(k_old, 3)) / (ci - pow(p, 2));
		err = abs(k_new - k_old);
	}

	/*lon(ii) = atan2(y, x);
	lat(ii) = atan2(z * k, p);
	R_N = a / (sqrt(1 - e ^ 2 * sin(lat(ii)) ^ 2));
	if lon(ii) > pi
		lon(ii) = lon(ii) - 2 * pi;
	end
		alt(ii) = p / cos(lat(ii)) - R_N;*/

	*longitude = atan2(y, x);
	*latitude = atan2(z * k_new, p);
	fp_precision R_N = a / sqrt(1 - pow(e, 2) * pow(sin(*latitude), 2));
	if (*longitude > M_PI) { *longitude -= 2 * M_PI; }

	*longitude *= 180.0 / M_PI;
	*latitude *= 180.0 / M_PI;
	*altitude = p / cos(*latitude) - R_N;
}	

fp_precision N_Power = 0;
