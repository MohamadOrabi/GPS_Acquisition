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

#include "GPS_Utils.h"
#include "Channel.h"
#include "FileReader.h"


Channel::Channel(RealArray PRN, const int sat_ID, const fp_precision fs, const fp_precision chip_rate, const int nsamples, const int loop_order_PLL, const fp_precision BPLL, const fp_precision BDLL, fp_precision fd, fp_precision code_shift, const int pp_moving_average, int ntrack)
	: sat_ID(sat_ID), fs(fs), chip_rate(chip_rate), nsamples(nsamples), loop_order_PLL(loop_order_PLL), BPLL(BPLL), BDLL(BDLL), pp_moving_average(pp_moving_average), Tsub (nsamples / (fs * 1E6))
{
	this->fd = fd;
	this->track_acc_PLL = 1;
	this->track_acc_DLL = 2;
	this->loop_filter_PLL = configureLoopFilter(loop_order_PLL, BPLL, track_acc_PLL* Tsub, fd);
	this->code_shift = code_shift;
	this->delta_shift = round(fs / chip_rate / 2);		//This sets the early and late
	this->PRN = PRN;
	this->dec_acc = 20;
	this->track_window = 1;
	this->ek_length = 1;
	this->tracking_threshold = 3.14;
	p = VEC(ntrack);
	e = VEC(ntrack);
	l = VEC(ntrack);
	pp = VEC(ntrack);
	ee = VEC(ntrack);
	ll = VEC(ntrack);
	fd_vec = RealArray(ntrack);
	IQ_Data = VEC(nsamples);

	decoded_bits.clear();

	file_reader_object.addChannelToList(this);
}



void Channel::track() {
	//SampleFrame_Ptr ptr;
	//bool keep_looping = true;
	int k_dec = 0;
	int ek_count = 0;
	IQData pk = 0, ek = 0, lk = 0, pk_dec = 0, pk_acc_PLL = 0, pk_acc_DLL = 0, ek_acc_DLL = 0, lk_acc_DLL = 0;
	RealVEC ek_PLL_vec(ek_length);
	while(k < p.size()) {
		if (workQ.size() > 0) {
			IQ_Data = workQ.front();

			//exportVEC("Data_" + to_string(k) + ".csv",IQ_Data);
			//k = ptr->k;
			//Resampling the PRN and mixing the signal must happen here because we need to keep track of last_t, last_theta, and current_index
			const int current_index = k * nsamples;
			VEC resampled_PRN = resamplePRN(PRN, nsamples, chip_rate, fs, code_length, current_index);
			RealArray t_array = creatRealArray(last_t, 1, last_t + nsamples - 1) / (fs * 1E6);
			this->last_t = t_array[nsamples - 1];	//Save the last value of t_array
			VEC mixed_signal = mixSignal(IQ_Data, t_array, fd, &last_theta);
			VEC shifted_PRN(nsamples);

			// Getting prompt
			double code_shift_samples = code_shift * fs * 1E6;
			shifted_PRN = shiftedByRows(resampled_PRN, round(code_shift_samples));
			IQData pk_now = (mixed_signal * shifted_PRN).sum() / (double)nsamples;
			/*if (is_tracking){
				if (is_decoding) {
					if (k_dec == dec_acc) {
						decoded_bits.push_back(signbit(real(pk_dec)));
						k_dec = 0;
						pk_dec = pk_now;
					}
					else {
						pk_dec += pk_now;
						k_dec++;
					}
				}
				else
				{
					if (real(pk_now) * real(pk) < 0) {
						is_decoding = true;
						pk_dec = pk_now;
					}
				}
			}
			else {
				if (k > ek_PLL_vec.size()) {
					fp_precision var_ek = ek_PLL_vec.norm();
					if (var_ek < tracking_threshold) {
						ek_count++;
					}
					else {
						ek_count = 0;
					}
					if (ek_count == track_window) {
						is_tracking = true;
					}
				}
			}*/
			pk = pk_now;
			p[k] = pk;
			pp[k] = pow(abs(p[k]), 2);

			//cout << "prompt:" << endl;
			//cout << p[k] << endl;
			//cout << pp[k] << endl;

			//Getting early
			shifted_PRN = shiftedByRows(resampled_PRN, round(code_shift_samples - delta_shift));
			ek = (mixed_signal * shifted_PRN).sum() / (double)nsamples;
			e[k] = ek;
			ee[k] = pow(abs(e[k]), 2);

			//cout << "early:" << endl;
			//cout << e[k] << endl;
			//cout << ee[k] << endl;

			//Getting late
			shifted_PRN = shiftedByRows(resampled_PRN, round(code_shift_samples + delta_shift));
			lk = (mixed_signal * shifted_PRN).sum() / (double)nsamples;
			l[k] = lk;
			ll[k] = pow(abs(l[k]), 2);

			//cout << "late:" << endl;
			//cout << l[k] << endl;
			//cout << ll[k] << endl;

			//Getting signal power from moving average
			double C_power;
			if (k < pp_moving_average) { C_power = real(pp.head(k + 1).sum()) / (k + 1); }
			else { C_power = real(pp.segment(k - pp_moving_average, pp_moving_average).sum()) / pp_moving_average; }
			//cout << C_power << endl;

			pk_acc_PLL += pk;
			pk_acc_DLL += pk;
			ek_acc_DLL += ek;
			lk_acc_DLL += lk;
			//if (real(pk) > 0) {
			//	pk_acc_PLL += pk;
			//	pk_acc_DLL += pk;
			//	ek_acc_DLL += ek;
			//	lk_acc_DLL += lk;
			//}
			//else {
			//	pk_acc_PLL += -pk;
			//	pk_acc_DLL += -pk;
			//	ek_acc_DLL += -ek;
			//	lk_acc_DLL += -lk;
			//}
			//Now time to update PLL :D
			if (k % track_acc_PLL == 0) {
				fp_precision ek_PLL = atan(imag(pk_acc_PLL) / real(pk_acc_PLL));	//Output of PLL discriminator
				if (k < ek_PLL_vec.size()) {
					ek_PLL_vec[k] = ek_PLL;
				}
				else
				{
					shiftedByRows(ek_PLL_vec, 1);
					ek_PLL_vec[0] = ek_PLL;
				}
				double vk_PLL = Update_PLL(&loop_filter_PLL, ek_PLL);
				fd = vk_PLL / (2 * M_PI);
				pk_acc_PLL = 0;
			}
			fd_vec[k] = fd;
			//cout << loop_filter_PLL.xk_PLL(0) << endl;

			//DLL
			if (k % track_acc_DLL == 0) {
				fp_precision eek = pow(abs(ek_acc_DLL), 2);
				fp_precision llk = pow(abs(lk_acc_DLL), 2);
				fp_precision ek_DLL = eek - llk;
				//cout << ek_DLL<< endl;
				fp_precision vk_DLL = 4.0 * BDLL * ek_DLL / (2 * sqrt(C_power) * chip_rate * 1E6);
				//cout << "vk_DLL: " << vk_DLL << endl;
				code_shift -= 1E-3 * track_acc_DLL * (vk_DLL + fd / 1575.42 / 1E6);
				pk_acc_DLL = 0;
				ek_acc_DLL = 0;
				lk_acc_DLL = 0;
			}

			//cout << code_shift << endl;
			//code_shift -= vk_DLL/(fs*1E6);
			//1E-3 * fd / 1575.42E6

			//cout << "ek_PLL: " << ek_PLL << endl;
			//cout << "fd: " << fd << endl << endl;

			workQ.pop();
			k++;
		}
	}
}

void Channel::addToWorkQ(VEC vec_to_add) {
	workQ.push(vec_to_add);
}

int Channel::getworkQsize()
{
	return workQ.size();
}


VEC Channel::getp() { return p; }
VEC Channel::gete() { return e; }
VEC Channel::getl() { return l; }
VEC Channel::getpp() { return pp; }
VEC Channel::getee() { return ee; }
VEC Channel::getll() { return ll; }
VEC Channel::getfd_vec() { return fd_vec; }