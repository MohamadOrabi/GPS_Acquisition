#define _USE_MATH_DEFINES

#include "stdafx.h"
#include "TypeDefs.h"
#include <fftw3.h>		//Only included to test the NNDLL
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

#include "GoogleEarthPath.hpp"
#include "GPS_Utils.h"
#include "Channel.h"
#include "FileReader.h"
#include "NavSol.h"



Channel::Channel(RealArray PRN, const int sat_ID, const fp_precision fs, const fp_precision chip_rate, const int nsamples, const int loop_order_PLL, const fp_precision BPLL, const fp_precision BDLL, fp_precision fd, fp_precision code_shift, const int pp_moving_average, int ntrack)
	: sat_ID(sat_ID), fs(fs), chip_rate(chip_rate), nsamples(nsamples), loop_order_PLL(loop_order_PLL), BPLL(BPLL), BDLL(BDLL), pp_moving_average(pp_moving_average), Tsub(nsamples / (fs * 1E6))
{
	this->fd = fd;

	this->n_acc_PLL = 1;
	this->n_acc_DLL = 1;

	this->loop_filter_PLL = configureLoopFilter(loop_order_PLL, BPLL, n_acc_PLL * Tsub, fd);
	if (code_shift <= 0.5E-3) {
		this->code_shift = code_shift;// +1E-3;
	}
	else {
		this->code_shift = code_shift -1E-3;
	}

	this->delta_shift = round(fs / chip_rate);		//This sets the early and late	//TODO: I removed the 2 here so check if it is correct
	//this->delta_shift = 2;
	this->PRN = PRN;
	this->dec_acc = 20;
	this->track_window = 1;
	this->ek_length = 1000;
	this->tracking_threshold = 0.35;

	this->transmission_time = 0;
	this->k_start_transmission_time = 0;
	this->channel_time = 0;

	this->t_OC = NAN; 
	this->a_f2 = NAN;
	this->a_f1 = NAN;
	this->a_f0 = NAN;
	this->deltat_sv = NAN;

	this->x_k = 0;
	this->y_k = 0;
	this->z_k = 0;

	//this->Q_mutex;
	this->done_flag = false;

	//Setting Ephemeris Parameters to NAN
	IODE = NAN;
	C_rs = NAN;
	deltan = NAN;
	M0 = NAN;
	C_uc = NAN;
	eccentricity = NAN;
	C_us = NAN;
	rad_A = NAN;
	t_oe = NAN;
	C_ic = NAN;
	omega_0 = NAN;
	C_is = NAN;
	i_0 = NAN;
	C_rc = NAN;
	Omega = NAN;
	omegadot = NAN;
	IDOT = NAN;


	p = VEC(ntrack);
	e = VEC(ntrack);
	l = VEC(ntrack);
	pp = VEC(ntrack);
	ee = VEC(ntrack);
	ll = VEC(ntrack);
	fd_vec = RealArray(ntrack);
	phase_error = VEC(ntrack);
	sf_decoded = vector<bool>(5,false);

	IQ_Data = VEC(nsamples);

	frames = vector<Subframe>(n_preamble_check, Subframe{vector<bool>(300),vector<bool>(3),vector<bool>(8),vector<bool>(6) });

	current_frame.words.resize(10);

	decoded_bits.clear();

	file_reader_object.addChannelToList(this);
}



void Channel::track() {
	GoogleEarthPath KMlPath("KML_" + to_string(sat_ID)+ ".kml", "Sat " + to_string(sat_ID));
	//SampleFrame_Ptr ptr;
	//bool keep_looping = true;
	int k_dec = 0;
	int ek_count = 0;
	IQData pk = 0, ek = 0, lk = 0, pk_dec = 0, pk_acc_PLL = 0, pk_acc_DLL = 0, ek_acc_DLL = 0, lk_acc_DLL = 0;
	RealVEC ek_PLL_vec(ek_length);	//Used to get the variance of ek_PLL
	//RealArray t_array = last_t + creatRealArray(last_t, 1, last_t + nsamples - 1) / (fs * 1E6);		//Put it here bcz we don't need to track last_t
	RealArray t_array = creatRealArray(0, 1, nsamples) / (fs * 1E6);

	while (k < p.size()) {		
		if (workQ.size() > 0) {
			Q_mutex.lock();
			SampleFrame_Ptr ptr = workQ.front();
			IQ_Data = ptr->sampleArray;
			workQ.pop();
			Q_mutex.unlock();
		
			//Resampling the PRN and mixing the signal must happen here because we need to keep track of last_t, last_theta, and current_index
			const int current_index = k* nsamples;	//Current index could be set to 0, did not make any difference
			// Getting prompt 
			double code_shift_samples = fmod(code_shift * fs * 1E6, nsamples);
			VEC resampled_PRN = resamplePRN(PRN, nsamples, chip_rate, fs, code_length, current_index);	//TODO: Move this outside of this function No need to do it multiple times
			//VEC resampled_PRN = resamplePRN(PRN, nsamples, chip_rate, fs, code_length, current_index - round(code_shift_samples));
			//this->last_t = t_array[nsamples - 1];	//Save the last value of t_array
			//last_t = 0;
			VEC mixed_signal = mixSignal(IQ_Data, fs, fd, &last_theta);
			//cout << last_theta << endl;

			VEC shifted_PRN(nsamples);
			shifted_PRN = shiftedByRows(resampled_PRN, round(code_shift_samples));	//TODO: floor vs round
			IQData pk_now = (mixed_signal * shifted_PRN).sum();// / (double)nsamples;

			//is_tracking and is_decoding logic 
			if (is_tracking){
				if (is_decoding) {
					if (k_dec == dec_acc-1) {
						decoded_bits.push_back(signbit(real(pk_dec)));
						k_dec = 0;
						pk_dec = pk_now;
						
						//if the preambles were detected
						if (preambles_detected) {
							n_bits_transmitted++;	//Getting number of bits transmitted since the beginning of the subframe
						}

						decode();
					}
					else {
						//n_CA_reapeats = (k-k_start_decoding) % dec_acc;	//Getting number of code repeats since detection of current bit
						pk_dec += pk_now;
						k_dec++;
					}
				}
				else
				{
					if (real(pk_now) * real(pk) < 0) {
						is_decoding = true;
						k_start_decoding = k;
						pk_dec = pk_now;
						cout << "PRN " << sat_ID  << " Decoding started at k: " << k << endl;
					}
				}
			}
			else {
				if (k > ek_PLL_vec.size()) {
					fp_precision var_ek = ek_PLL_vec.norm()/sqrt(ek_PLL_vec.size());
					if (var_ek < tracking_threshold) {
						ek_count++;
					}
					else {
						ek_count = 0;
					}
					if (ek_count == track_window) {
						is_tracking = true;
						cout << "PRN " << sat_ID << " Tracking started at k: " << k << endl;
					}
				}
			}

			pk = pk_now;
			p[k] = pk;
			//pp[k] = pow(abs(p[k]), 2);

			//Getting early
			shifted_PRN = shiftedByRows(resampled_PRN, round(code_shift_samples - delta_shift));
			ek = (mixed_signal * shifted_PRN).sum();// / (double)nsamples;
			e[k] = ek;
			ee[k] = pow(abs(e[k]), 2);

			//Getting late
			shifted_PRN = shiftedByRows(resampled_PRN, round(code_shift_samples + delta_shift));
			lk = (mixed_signal * shifted_PRN).sum();// / (double)nsamples;
			l[k] = lk;
			ll[k] = pow(abs(l[k]), 2);


			//Accumulating pk, ek, and lk for PLL and DLL
			if (real(pk) > 0) {
				pk_acc_PLL += pk / (double)n_acc_PLL;
				pk_acc_DLL += pk / (double)n_acc_DLL;
				ek_acc_DLL += ek / (double)n_acc_DLL;
				lk_acc_DLL += lk / (double)n_acc_DLL;
			}
			else {
				pk_acc_PLL += -pk / (double)n_acc_PLL;
				pk_acc_DLL += -pk / (double)n_acc_DLL;
				ek_acc_DLL += -ek / (double)n_acc_DLL;
				lk_acc_DLL += -lk / (double)n_acc_DLL;
			}

			//Getting signal power and variance from moving average
			double C_power;
			if (k % n_acc_DLL == 0) {
				pp[k] = pow(abs(pk_acc_DLL),2);
				if (k < pp_moving_average * n_acc_DLL) {
					C_power = real(pp.head(k*n_acc_DLL + 1).sum()) / (k + 1); 
					//C_power = pow(abs(p.head(k + 1).sum()), 2) / (k + 1);
					C_N = (C_power - N_Power) / (N_Power * n_acc_DLL * nsamples / (fs * 1E6));
					variance = pow(299792458, 2) * (0.005 * 0.8 * pow(1 / 1.023E6, 2) / (2 * C_N) * (1 + 2 / (1E-3 * C_N)));
					variance = 2 + 1000 * variance;
				}
				else {
					C_power = real(pp.segment(k - pp_moving_average * n_acc_DLL, pp_moving_average * n_acc_DLL + 1).sum()) / pp_moving_average;
					//C_power = pow(abs(p.segment(k - pp_moving_average * n_acc_DLL, pp_moving_average * n_acc_DLL).sum()), 2) / (pp_moving_average * n_acc_DLL);
					C_N = (C_power - N_Power) / (N_Power * n_acc_DLL * nsamples / (fs * 1E6));
					variance = pow(299792458, 2) * (0.005 * 0.8 * pow(1 / 1.023E6, 2) / (2 * C_N) * (1 + 2 / (1E-3 * C_N)));
					variance = 2 + 1000 * variance;
				}
			}


			//Now time to update PLL :D
			if (k % n_acc_PLL == 0) {	//if it is time to update the PLL
				fp_precision ek_PLL = atan(imag(pk_acc_PLL) / real(pk_acc_PLL));	//Output of PLL discriminator
				phase_error[k] = ek_PLL;

				//create moving average for ek
				if (k < ek_PLL_vec.size()) {
					ek_PLL_vec[k] = ek_PLL;
				} 
				else {
					ek_PLL_vec = shiftedByRowsReal(ek_PLL_vec, -1);
					ek_PLL_vec[ek_PLL_vec.size()-1] = ek_PLL;
				}

				double vk_PLL = Update_PLL(&loop_filter_PLL, ek_PLL);
				fd = vk_PLL / (2 * M_PI);
				pk_acc_PLL = 0;
			}
			fd_vec[k] = fd;

			//DLL
			//if (k % n_acc_DLL == 0) { //if time to update DLL
			//	fp_precision eek = pow(abs(ek_acc_DLL), 2);
			//	fp_precision llk = pow(abs(lk_acc_DLL), 2);
			//	fp_precision ek_DLL = eek - llk;
			//	//fp_precision ek_DLL = (real(ek_acc_DLL) - real(lk_acc_DLL)) * real(pk_acc_DLL) + (imag(ek_acc_DLL) - imag(lk_acc_DLL)) * imag(pk_acc_DLL);
			//	//phase_error[k] = ek_DLL;
			//	fp_precision vk_DLL = 4.0 * BDLL * ek_DLL / (2 * C_power * chip_rate * 1E6);
			//	code_shift -= 1E-3 * n_acc_DLL * (vk_DLL + 1 * (fd / 1575.42 / 1E6));
			//	//if (code_shift >= 1E-3) { code_shift = fmod(code_shift, 1E-3); }
			//	pk_acc_DLL = 0;
			//	ek_acc_DLL = 0;
			//	lk_acc_DLL = 0;
			//}

			////NNDLL
			VEC IQ_Data_fft(nsamples);		//Data fft
			VEC resampled_PRN_fft(nsamples);
			VEC corr(nsamples);
			RealArray corr_power(nsamples);			//Abs2/nsamples of correlation
			RealArray samplesData(10);
			//Add mutex maybe
			fft_mutex.lock();
			conduct_FFT(mixed_signal, IQ_Data_fft.data());
			conduct_FFT(resampled_PRN, resampled_PRN_fft.data());
			corr = Corr(IQ_Data_fft, resampled_PRN_fft);
			fft_mutex.unlock();
			corr_power = corr.cwiseAbs2();
			for (int i = 0; i < 10; i++) {
				int index = i + code_shift_samples - 5 + 1;
				if (index >= 2500) { index -= 2500; 
				} else if(index < 0) { index += 2500;
				}
				samplesData[i] = corr_power[index];
			}
			fft_mutex.lock();
			fp_precision vk_DLL = -round(predict_shift(samplesData))/(fs*1e6);
			fft_mutex.unlock();
			code_shift -= 1E-3 * n_acc_DLL * (vk_DLL + 1 * (fd / 1575.42 / 1E6));

			//Calculating the transmission time
			if (decoded_ephemeris && decoded_clock_correction) {
				if (TOW % 6 != 0) { cout << "TOW is not a multiple of 6!!!" << endl; }
				transmission_time = 1 * ((1 * TOW) + 1 * n_bits_transmitted * 20E-3 + k_dec * 1E-3 - 1 * code_shift); //fmod(code_shift_mod, 1E-3));
				
				//Check for clock correction parameters & correct clock
				//if (sf_decoded[0]) {
				if (decoded_clock_correction) {
					fp_precision deltat_r = F * eccentricity * rad_A * sin(E_k);
					deltat_sv = 1 * (a_f0 + a_f1 * (transmission_time - t_OC) + a_f2 * pow(transmission_time - t_OC, 2) + deltat_r);
					transmission_time = transmission_time - deltat_sv;
				}

				//Calculate Receiver_time_shift which is the start time of the receiver
				//receiver_time_shift_mutex.lock();	//TODO: Check if this mutex is needed!
				if (receiver_time_shift == 0) {
					//receiver_time_shift = TOW - ptr->k / 1.0E3 - 1 * code_shift - deltat_sv;// +1.0 / 15.0;		//0.0867 is the time needed for light to cross 26000Km (Added for simulation)
					receiver_time_shift = transmission_time - ptr->k/1.0E3;
					cout << "Receiver time shift TOW: " << receiver_time_shift << "s for channel: " << sat_ID;
				}
				//receiver_time_shift_mutex.unlock();

				//fp_precision channel_time = receiver_time_shift + ptr->k / 1.0E3;
				this->channel_time = receiver_time_shift + ptr->k / 1.0E3;

				vec_mutex.lock();	//lock mutex to prevent data race while accessing vectors
				transmission_time_vec.push_back(1 * channel_time-transmission_time);
				variance_vec.push_back(variance);
				//variance_vec.push_back(C_N);
				vec_mutex.unlock();
				//debug_vec.push_back(1 * channel_time-transmission_time);
				getSVCoordinates();
				

				//Adding coordinates to kml file
				if (x_kvec.size() % 1000 == 0) {
					fp_precision longitude = 0; fp_precision latitude = 0; fp_precision altitude = 0;
					ecef2lla(x_kvec[x_kvec.size() - 1], y_kvec[y_kvec.size() - 1], z_kvec[z_kvec.size() - 1], &latitude, &longitude, &altitude);
					KMlPath.addPoint(longitude, latitude);
					//cout << "altitude: " << altitude << endl;
				}

				//If this is the first transmission time value
				if (k_start_transmission_time == 0) {
					k_start_transmission_time = ptr->k;
					navsol_object.AddToChannels(this, k_start_transmission_time);
					n_channels_decoded++;
					if (n_channels_decoded == 4) {
						NavSolThread = thread(&NavSol::Calculate_Position, &navsol_object);
					}
				}

			}
			k++;
		}
	}

	n_channels_decoded--;

	exportArray("Channel_Data_" + to_string(sat_ID) + ".csv", this->exportChannelData());
	exportVEC("fd_" + to_string(sat_ID) + ".csv", this->getfd_vec());
	
	exportVEC("pp_" + to_string(sat_ID) + ".csv", this->getpp());
	exportVEC("ee_" + to_string(sat_ID) + ".csv", this->getee());
	exportVEC("ll_" + to_string(sat_ID) + ".csv", this->getll());
	
	exportVEC("phase_error_" + to_string(sat_ID) + ".csv", this->getphase_error());
	exportBits("decoded_bits_" + to_string(sat_ID) + ".csv", this->getdecodedBits());
	exportvec("debug_vec_" + to_string(sat_ID) + ".csv", this->getdebug_vec());
	exportvec("transmission_time_" + to_string(sat_ID) + ".csv", this->gettransmission_time_vec());

	exportVEC("I_" + to_string(sat_ID) + ".csv", real(this->getp()));
	exportVEC("Q_" + to_string(sat_ID) + ".csv", imag(this->getp()));

	cout << "Data Was exported" << endl;
	done_flag = true;

	Channel::navsolJoinMutex.lock();
	if (NavSolThread.joinable()) {
		NavSolThread.join();
		cout << "NavSolThread was joined" << endl;
	}
	Channel::navsolJoinMutex.unlock();

}

void Channel::addToWorkQ(SampleFrame_Ptr vec_to_add) {
	workQ.push(vec_to_add);
}

int Channel::getworkQsize()
{
	return workQ.size();
}

bool Channel::getNeedDataFlag() {
	if (workQ.size() <= 20) {
		return true;
	}
	return false;
}

bool Channel::check_nav_index(int index)
{
	//cout << "current_index: " << index << endl;
	//cout << "k_start_tras_time: " << k_start_transmission_time << endl;
	index = index - k_start_transmission_time;
	if (index < 0) {
		//cout << "index is less than 0: " << index << endl;
		return true;
	}
	else {
		//cout << "index: " << index << endl;
		return index < z_kvec.size();
	}
}

void Channel::decode() {
	//if preambles detected flag is not set, see if it should be
	if (!preambles_detected) { 
		this->preambles_detected = check_consecutive_preamble();
		//if it didn't work check again with bits flipped
		//if (!preambles_detected) { 
		//	flip = !flip;
		//	preambles_detected = check_consecutive_preamble();
		//}
		//if preambles detected flag is set, then get the past n_preamble_check frames
		if (preambles_detected) {
			this->getFrames();		
			index_of_first_preamble = decoded_bits.size() - 8; //Index of the last preamble in preambles + 300 here is because
			//past_frame = frames[n_preamble_check - 1];

			past_frame = frames[0];
			//Process and decode past_frames
			for (int i = 1; i < n_preamble_check; i++) {
				//current_frame = frames[i];
				
				process_frame(&frames[i]);
				n_bits_transmitted = 300;
				deocode_frame_parameters(frames[i]);
				n_bits_transmitted = 0;
				past_frame = frames[i];
			}
			//TOW = 0;
			n_bits_transmitted = 8; //Because this happened after reading 8 bits of the preamble
			//TOW = 0;
			//decoded_ephemeris = false;
			past_frame = frames[n_preamble_check - 1];

			//cout << "Getting past frame ID: " << past_frame.frameID_int << " at: " << past_frame.index << endl;
			//printBits(past_frame.preamble);

			//cout << "Flip was: " << flip << endl;
			//return;
		}
	} else {
		//vector<bool> error_flags(10,false);	//Error flags for every word in the current frame
		//bool dropped = false;	//Error Flag for the entire frame
		//If consecutive preambles have already been detected
		if ((decoded_bits.size() - index_of_first_preamble) % 300 == 0) {

			getCurrentFrame();	//Gets the current_frame (300 bits)

			process_frame(&current_frame);	//Extracts Data bits and gets received frame ID

			if (current_frame.frameID_int == past_frame.frameID_int % 5 + 1) {
			//if (!error_flags[1] && current_frame.frameID_int < 6 && current_frame.frameID_int > 0) {

				deocode_frame_parameters(current_frame);
			}				
			else {
				//This means that the curent frame ID was bad
				//cout << "Error in 2nd word, frame ID could be wrong!" << endl;
				//current_frame.frameID_int = past_frame.frameID_int % 5 + 1;
				//preambles_detected = false;		//TODO: remove this when flip = !flip works
			}
			
			past_frame = current_frame;		//TODO: put this back inside if consecutive preambles have already been detected
			
		}
	
	}
}

void Channel::process_frame(Subframe *frame_to_process) {
	vector<bool> last_2_bits = getBits(past_frame.bits, 299, 2);	//Stores last 2 bits of previous frame

	//Checking for parity && extracting data bits
	for (int n_word = 0; n_word < 10; n_word++) {

		vector<bool> parity_bits = getBits((*frame_to_process).words[n_word], 25, 6);	//Getting unflipped parity bits

		//Extract Raw data (d) from the frame
		if (last_2_bits[1]) {
			(*frame_to_process).words[n_word].flip();
		}

		vector<bool> calculated_pb = (*frame_to_process).CalculateParityBits(n_word, last_2_bits);

		(*frame_to_process).error_flags[n_word] = parity_bits != calculated_pb;
		if ((*frame_to_process).error_flags[n_word]) { /*dropped = true;*/  cout << "Sat " << sat_ID << " Frame " << past_frame.frameID_int % 5 + 2 << "* was dropped" << endl; }

		//Print parity bits
		cout << "Read pb: \tword#" << n_word << "\t" << "Error: " << (*frame_to_process).error_flags[n_word] << "\t"; printBits(parity_bits);
		cout << "Calculated: \tword#" << n_word << "\t" << "Error: " << (*frame_to_process).error_flags[n_word] << "\t"; printBits(calculated_pb);

		last_2_bits = getBits((*frame_to_process).bits, 30 * n_word + 29, 2);	//Keeping track of last 2 bits in previous words
	}

	//Getting correct frame ID after extracting data bits
	(*frame_to_process).frameID = getBits((*frame_to_process).words[1], 20, 3);
	(*frame_to_process).frameID_int = bin2dec((*frame_to_process).frameID);
	cout << "Current frame ID: " << (*frame_to_process).frameID_int << endl;
}

void Channel::deocode_frame_parameters(Subframe frame_to_decode) {

	//if (!sf_decoded[current_frame.frameID_int - 1] && !dropped) { //&& !dropped) {
	//	sf_decoded[current_frame.frameID_int - 1] = true;	//Set the subframe decoded flag to 1
	//	cout << "Decoded the first frame with ID " << current_frame.frameID_int << endl;
	//	if (current_frame.frameID_int == 1) { cout << "Clock correction starts at k: " << k << endl; }
	//}

	//Get TOW
	vector<bool> error_flags = frame_to_decode.error_flags;
	vector<bool> TOW_bits = getBits(frame_to_decode.words[1], 1, 17);
	TOW_bits.push_back(false); TOW_bits.push_back(false);	//Adding the last 2 LSBs to TOW
	fp_precision new_TOW = bin2dec(TOW_bits) * 1.5;

	//Correcting TOW if necessary
	if (new_TOW == TOW + n_bits_transmitted / 50.0 || TOW == 0 || TOW % 6 != 0) {

		TOW = new_TOW;
		if (TOW % 6 != 0) {
			cout << "TOW is not a multiple of 6" << endl;
		}
		n_bits_transmitted = 0;
		cout << "TOW is: " << TOW << endl;
	}
	else {
		TOW += n_bits_transmitted / 50.0;
		n_bits_transmitted = 0;
		cout << "Error with TOW! New TOW is: " << new_TOW;
	}
	//TOW = new_TOW;
	if (frame_to_decode.frameID_int == past_frame.frameID_int % 5 + 1) {


		//Add Flags for detected errors and drop word if error is detected
		if (frame_to_decode.frameID_int == 1) {	//SUBFRAME1

			if (!error_flags[2]) {
				vector<bool> WeekNumber_bits = getBits(frame_to_decode.words[2], 1, 10);
				WeekNumber = bin2dec(WeekNumber_bits);
				cout << "Week Number: " << WeekNumber << endl;
			}
			if (!error_flags[7]) {
				vector<bool> t_OC_bits = getBits(frame_to_decode.words[7], 9, 16);
				t_OC = bin2dec(t_OC_bits) * pow(2, 4);
			}
			if (!error_flags[8]) {
				vector<bool> a_f2_bits = getBits(frame_to_decode.words[8], 1, 8);
				a_f2 = twosCompliment2Dec(a_f2_bits) * pow(2, -55);
				vector<bool> a_f1_bits = getBits(frame_to_decode.words[8], 9, 16);
				a_f1 = twosCompliment2Dec(a_f1_bits) * pow(2, -43);
			}
			if (!error_flags[9]) {
				vector<bool> a_f0_bits = getBits(frame_to_decode.words[9], 1, 22);
				a_f0 = twosCompliment2Dec(a_f0_bits) * pow(2, -31);
			}

			//a_f2 = 0;
			//a_f1 = -1.932676241267E-12;
			//a_f0 = 2.618110738695E-04;	//Remove thissss!!!!!!!!!

		}
		else if (frame_to_decode.frameID_int == 2) {	//SUBFRAME1

			if (!error_flags[7] && !error_flags[8]) {
				vector<bool> rad_A_bits = getBits(frame_to_decode.words[7], 17, 8);
				vector<bool> rad_A_bits1 = getBits(frame_to_decode.words[8], 1, 24);
				rad_A_bits.insert(rad_A_bits.end(), rad_A_bits1.begin(), rad_A_bits1.end());
				rad_A = bin2dec(rad_A_bits) * pow(2, -19);
				//rad_A = 5.153766977310E03;
				cout << "Semi-major axis" << pow(rad_A, 2) / 1000 << " Km" << endl;
			}

			if (!error_flags[2]) {
				vector<bool> IODE_bits = getBits(frame_to_decode.words[2], 1, 8);
				IODE = bin2dec(IODE_bits);

				vector<bool> C_rs_bits = getBits(frame_to_decode.words[2], 9, 16);
				C_rs = twosCompliment2Dec(C_rs_bits) * pow(2, -5);
				//C_rs = -1.009375000000E01;
			}
			if (!error_flags[3]) {
				vector<bool> deltan_bits = getBits(frame_to_decode.words[3], 1, 16);
				deltan = twosCompliment2Dec(deltan_bits) * pow(2, -43) * M_PI;
				//deltan = 4.650550856967E-09;
			}
			if (!error_flags[3] && !error_flags[4]) {
				vector<bool> M0_bits = getBits(frame_to_decode.words[3], 17, 8);
				vector<bool> M0_bits1 = getBits(frame_to_decode.words[4], 1, 24);
				M0_bits.insert(M0_bits.end(), M0_bits1.begin(), M0_bits1.end());
				M0 = twosCompliment2Dec(M0_bits) * pow(2, -31) * M_PI;	//Pi for normalization
				//M0 = 1.217764630156;
			}
			if (!error_flags[5]) {
				vector<bool> C_uc_bits = getBits(frame_to_decode.words[5], 1, 16);
				C_uc = twosCompliment2Dec(C_uc_bits) * pow(2, -29);
				//C_uc = -4.433095455170E-07;
			}
			if (!error_flags[5] & !error_flags[6]) {
				vector<bool> eccentricity_bits = getBits(frame_to_decode.words[5], 17, 8);
				vector<bool> eccentricity_bits1 = getBits(frame_to_decode.words[6], 1, 24);
				eccentricity_bits.insert(eccentricity_bits.end(), eccentricity_bits1.begin(), eccentricity_bits1.end());
				eccentricity = bin2dec(eccentricity_bits) * pow(2, -33);
				//eccentricity = 8.155977586284E-03;
			}
			if (!error_flags[7]) {
				vector<bool> C_us_bits = getBits(frame_to_decode.words[7], 1, 16);
				C_us = twosCompliment2Dec(C_us_bits) * pow(2, -29);
				//C_us = 5.740672349930E-06;
			}
			if (!error_flags[9]) {
				vector<bool> t_oe_bits = getBits(frame_to_decode.words[9], 1, 16);
				t_oe = bin2dec(t_oe_bits) * pow(2, 4);
				//t_oe = 1.224000000000E05;
				cout << "t_oe: " << t_oe << endl;

			}

		}
		else if (frame_to_decode.frameID_int == 3) {

			if (!error_flags[2]) {
				vector<bool> C_ic_bits = getBits(frame_to_decode.words[2], 1, 16);
				C_ic = twosCompliment2Dec(C_ic_bits) * pow(2, -29);
				//C_ic = 1.322478055954E-07;


			}
			if (!error_flags[2] && !error_flags[3]) {
				vector<bool> omega_0_bits = getBits(frame_to_decode.words[2], 17, 8);
				vector<bool> omega_0_bits1 = getBits(frame_to_decode.words[3], 1, 24);
				omega_0_bits.insert(omega_0_bits.end(), omega_0_bits1.begin(), omega_0_bits1.end());
				omega_0 = twosCompliment2Dec(omega_0_bits) * pow(2, -31) * M_PI;	//Pi for normalization
				//omega_0 = 3.835774861134E-01;
			}
			if (!error_flags[4]) {
				vector<bool> C_is_bits = getBits(frame_to_decode.words[4], 1, 16);
				C_is = twosCompliment2Dec(C_is_bits) * pow(2, -29);
				//C_is = 6.705522537231E-08;
			}
			if (!error_flags[4] && !error_flags[5]) {
				vector<bool> i_0_bits = getBits(frame_to_decode.words[4], 17, 8);
				vector<bool> i_0_bits1 = getBits(frame_to_decode.words[5], 1, 24);
				i_0_bits.insert(i_0_bits.end(), i_0_bits1.begin(), i_0_bits1.end());
				i_0 = twosCompliment2Dec(i_0_bits) * pow(2, -31) * M_PI;
				//i_0 = 9.711274397902E-01;
			}
			if (!error_flags[6]) {
				vector<bool> C_rc_bits = getBits(frame_to_decode.words[6], 1, 16);
				C_rc = twosCompliment2Dec(C_rc_bits) * pow(2, -5);
				//C_rc = 2.689062500000E02;
			}

			if (!error_flags[6] && !error_flags[7]) {
				vector<bool> Omega_bits = getBits(frame_to_decode.words[6], 17, 8);
				vector<bool> Omega_bits1 = getBits(frame_to_decode.words[7], 1, 24);
				Omega_bits.insert(Omega_bits.end(), Omega_bits1.begin(), Omega_bits1.end());
				Omega = twosCompliment2Dec(Omega_bits) * pow(2, -31) * M_PI;
				//Omega = -3.951734701155E-01;
			}
			if (!error_flags[8]) {
				vector<bool> omegadot_bits = getBits(frame_to_decode.words[8], 1, 24);
				omegadot = twosCompliment2Dec(omegadot_bits) * pow(2, -43) * M_PI;
				//omegadot = -8.213913571042E-09;
				cout << "omega dot: " << omegadot << endl;
			}
			if (!error_flags[9]) {
				vector<bool> IODE_bits = getBits(frame_to_decode.words[9], 1, 8);
				IODE = bin2dec(IODE_bits);
			}
			if (!error_flags[10]) {
				vector<bool> IDOT_bits = getBits(frame_to_decode.words[9], 9, 14);
				IDOT = twosCompliment2Dec(IDOT_bits) * pow(2, -43) * M_PI;
				//IDOT = -4.732339978098E-10;
			}
			//getSVCoordinates();
		}
		else if (frame_to_decode.frameID_int == 4) {

		}
		else if (frame_to_decode.frameID_int == 5) {

		}
		else {
			cout << "Frame ID is invalid! (not 1->5): " << frame_to_decode.frameID_int << endl;
		}

		if (!decoded_ephemeris) { check_ephemeris(); }
		if (!decoded_clock_correction) { check_clock_correction(); }
	}
}

bool Channel::check_consecutive_preamble() {
	bool is_preamble = true; //300 bits (6 seconds) is the period of the TLM message also subframe size

	if (decoded_bits.size() > (n_preamble_check) * 300) {	//Decoded bits must be greater than the period to check
		for (int i = 1; i <= n_preamble_check; i++) {	//TODO: See where i should sart and end
			vector<bool> tailbits = subvector(decoded_bits, 8 + 300 * i, 8);

			if ((tailbits == preamble1 || tailbits == preamble2) && is_preamble) {
				cout << "preamble at bit: " << decoded_bits.size() - 8 - 300 * i << endl;
			}
			else { is_preamble = false; }
		}

		if (is_preamble) {
			cout << "Consecutive preambles found at: " << decoded_bits.size() - 8 - 300*(n_preamble_check) << endl;
		}


	} else { is_preamble = false; }
	return is_preamble;
}

void Channel::getCurrentFrame() {
	current_frame.bits = subvector(decoded_bits, 300, 300);
	current_frame.index = decoded_bits.size() - 300;
	//if (flip) { current_frame.bits.flip(); }
	current_frame.preamble = getBits(current_frame.bits, 1, 8);
	if (current_frame.preamble != preamble1) {
		//printBits(current_frame.preamble);
		//current_frame.flipped = true;
		//current_frame.bits.flip();
		//current_frame.preamble.flip();
		if (current_frame.preamble != preamble1) {
			//cout << "This frame is bad" << endl;
			//printBits(current_frame.preamble);
		}
	}
	current_frame.error_flags = vector<bool>(10, false);
	current_frame.setup();
}

void Channel::getSVCoordinates() {
	fp_precision t_k = transmission_time - t_oe - t_flight;// -(channel_time - transmission_time);
	if (t_k > 302400) { t_k -= 604800; }
	else if (t_k < -302400) { t_k += 604800; }

	fp_precision mu = 3.986005E14;
	fp_precision n = sqrt(mu) / pow(rad_A, 3) + deltan;
	fp_precision M_k = M0 + n * t_k;

	//Newton iteration for E_k .. M_k = E_k - e*sin(E_k)
	int max_iterations = 5;
	fp_precision oldx = 1;
	fp_precision error;
	for (int i = 0; i < max_iterations; i++) {
		fp_precision newx = oldx - (oldx - eccentricity * sin(oldx) - M_k) / (1 - eccentricity * cos(oldx));
		error = newx - eccentricity * sin(newx) - M_k;
		E_k = newx;
		oldx = newx;
	}
	if (error >= 1E-6) {
		cout << "Error for E_k is large: " << error << endl;
	}

	fp_precision v_k = atan2(sqrt(1 - pow(eccentricity, 2)) * sin(E_k) , (cos(E_k) - eccentricity));
	fp_precision arg_of_lat = v_k + Omega;
	fp_precision du_k = C_us * sin(2 * arg_of_lat) + C_uc * cos(2 * arg_of_lat);
	fp_precision dr_k = C_rs * sin(2 * arg_of_lat) + C_rc * cos(2 * arg_of_lat);
	fp_precision di_k = C_is * sin(2 * arg_of_lat) + C_ic * cos(2 * arg_of_lat);

	fp_precision corrected_arg_of_lat = arg_of_lat + du_k;
	fp_precision corrected_r_k = pow(rad_A, 2) * (1 - eccentricity * cos(E_k)) + dr_k;
	fp_precision corrected_i_k = i_0 + IDOT * t_k + di_k;


	fp_precision x_kp = corrected_r_k * cos(corrected_arg_of_lat);	//p stands for prime
	fp_precision y_kp = corrected_r_k * sin(corrected_arg_of_lat);

	fp_precision omegadot_e = 7.2921151467E-5;
	fp_precision omega_k = omega_0 + (omegadot - omegadot_e) * t_k - omegadot_e * t_oe;

	x_k = x_kp * cos(omega_k) - y_kp * cos(corrected_i_k) * sin(omega_k);
	y_k = x_kp * sin(omega_k) + y_kp * cos(corrected_i_k) * cos(omega_k);
	//y_k = -y_k;
	z_k = y_kp * sin(corrected_i_k);
	//cout << "Norm of position vector: " << (sqrt(pow(x_k, 2) + pow(y_k, 2) + pow(z_k, 2))) / 1000.0 << "Km" << endl;
	
	vec_mutex.lock();
	x_kvec.push_back(x_k);
	y_kvec.push_back(y_k);
	z_kvec.push_back(z_k);
	vec_mutex.unlock();
}

void Channel::check_ephemeris()
{
	/*bool decoded_ephemeris = (t_OC != 0 &&
	M0 != 0 && deltan != 0 && eccentricity != 0 && rad_A != 0 && omega_0 != 0 && i_0 != 0 && Omega != 0 && omegadot != 0 && IDOT != 0 && C_uc != 0 && 
		C_us != 0 && C_rc != 0 && C_rs != 0 && C_ic != 0 && C_is != 0 && t_oe != 0 && IODE != 0);*/
	//IODE = NAN;
	//C_rs = NAN;
	//deltan = NAN;
	//M0 = NAN;
	//C_uc = NAN;
	//eccentricity = NAN;
	//C_us = NAN;
	//rad_A = NAN;
	//t_oe = NAN;
	//C_ic = NAN;
	//omega_0 = NAN;
	//C_is = NAN;
	//i_0 = NAN;
	//C_rc = NAN;
	//Omega = NAN;
	//omegadot = NAN;
	//IDOT = NAN;
	
	this->decoded_ephemeris = !(isnan(IODE) || isnan(C_rs) || isnan(deltan) || isnan(M0) || isnan(C_uc) || isnan(eccentricity) || isnan(C_us) || isnan(rad_A) || isnan(t_oe)
		|| isnan(C_ic) || isnan(omega_0) || isnan(C_is) || isnan(i_0) || isnan(C_rc) || isnan(Omega) || isnan(omegadot) || isnan(IDOT)); // || not && because there is a not a the beginning
	if (decoded_ephemeris) { cout << "decoded_ephemeris = true at k: " << k << endl; }
}

void Channel::check_clock_correction() {
	/*this->t_OC = NAN;
	this->a_f2 = NAN;
	this->a_f1 = NAN;
	this->a_f0 = NAN;*/

	this->decoded_clock_correction = !(isnan(t_OC) || isnan(a_f2) || isnan(a_f1) || isnan(a_f0));
	if (decoded_clock_correction) { cout << "decoded_clock_correction = true at k: " << k << endl;
	debug_vec.clear();
	debug_vec.push_back(a_f0); debug_vec.push_back(a_f1); debug_vec.push_back(a_f2);
	}
}

void Channel::getFrames() {
	for (int i = 0; i < n_preamble_check; i++) {
		frames[i].bits = subvector(decoded_bits, 8 + 300 * (n_preamble_check - i), 300);
		frames[i].index = decoded_bits.size() - 8 - 300 * (n_preamble_check - i);
		//if (flip) { frames[i].bits.flip(); }
		frames[i].preamble = getBits(frames[i].bits, 1, 8);
		//if (frames[i].preamble != preamble1) { 
		//	frames[i].flipped = true; 
		//	frames[i].bits.flip();
		//	//frames[i].preamble.flip();
		//	if (frames[i].preamble != preamble1) { 
		//		cout << "This frame is bad" << endl; 
		//	}
		//}
		frames[i].error_flags = vector<bool>(10, false);
		frames[i].setup();
		cout << "subframe ID: " << frames[i].frameID[0] << frames[i].frameID[1] << frames[i].frameID[2] << " Decimal: "<< frames[i].frameID_int << endl;
	}
}

Eigen::Matrix<fp_precision, Eigen::Dynamic, 5 > Channel::exportChannelData()
{	
	if (is_tracking) {
		if (is_decoding) {
			if (k_start_transmission_time != 0) {
				int size = p.size() - k_start_transmission_time;
				Eigen::Matrix<fp_precision, Eigen::Dynamic, 5 > ChannelData;
				ChannelData.resize(size, 5);

				ChannelData.col(0) = creatRealArray(k_start_transmission_time + 1, 1, p.size());
				ChannelData.col(1) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(transmission_time_vec.data(), size);
				ChannelData.col(2) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(x_kvec.data(), size);
				ChannelData.col(3) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(y_kvec.data(), size);
				ChannelData.col(4) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(z_kvec.data(), size);


				//for (int i = 0; i < size; i++) {
				//	ChannelData.col(1)[i] = transmission_time_vec[i];
				//}
				//cout << std::setprecision(12) << ChannelData.col(0)[2] << endl;
				//ChannelData.col(2) = 
				return ChannelData;
			} else {
				cout << "PRN# " << sat_ID << " output vecs are bad .. transmission_time = 0" << endl;
			}
		} else {
			cout << "PRN# " << sat_ID << " output vecs are bad .. did not decode" << endl;
		}
	} else {
		cout << "PRN# " << sat_ID << " output vecs are bad .. did not track" << endl;
	}
	return Eigen::Matrix<fp_precision, Eigen::Dynamic, 5 >();
}

VEC Channel::getp() { return p; }
VEC Channel::gete() { return e; }
VEC Channel::getl() { return l; }
VEC Channel::getpp() { return pp; }
VEC Channel::getee() { return ee; }
VEC Channel::getll() { return ll; }

VEC Channel::getfd_vec() { return fd_vec; }
VEC Channel::getphase_error() { return phase_error; }
vector<fp_precision> Channel::getdebug_vec() { return debug_vec; }

vector<fp_precision> Channel::gettransmission_time_vec() { return transmission_time_vec; }
vector<fp_precision> Channel::getx_kvec() { return x_kvec; }
vector<fp_precision> Channel::gety_kvec() { return y_kvec; }
vector<fp_precision> Channel::getz_kvec() { return z_kvec; }

fp_precision Channel::getx_kat(int index) { 
	index = index - k_start_transmission_time;
	if (index < this->x_kvec.size() && index >= 0) {
		return this->x_kvec[index];
	}
	else {
		//navsol_object.set_calc_flag(false);
		//cout << "x index  " << index << "\tx size: " << x_kvec.size() << endl;
		return 0;
	}
}
fp_precision Channel::gety_kat(int index) {
	index = index - k_start_transmission_time;
	if (index < this->y_kvec.size()) {
		return this->y_kvec[index];
	}
	else {
		//navsol_object.set_calc_flag(false);
		cout << "y index  " << index << "\ty size: " << y_kvec.size() << endl;
		return 0;
	}
}
fp_precision Channel::getz_kat(int index) {
	index = index - k_start_transmission_time;
	if (index < this->z_kvec.size()) {
		//cout << "Success for index: " << index << endl;
		return this->z_kvec[index];
	}
	else {
		//navsol_object.set_calc_flag(false);
		cout << "z index  " << index << "\tz size: " << z_kvec.size() << endl;
		return 0;
	}
}
fp_precision Channel::getVarianceat(int index)
{
	index = index - k_start_transmission_time;
	if (index < this->variance_vec.size()) {
		//cout << "Success for index: " << index << endl;
		return this->variance_vec[index];
	}
	else {
		//navsol_object.set_calc_flag(false);
		cout << "C_N index  " << index << "\tz size: " << variance_vec.size() << endl;
		return 0;
	}
}
fp_precision Channel::getrange(int index) {
	index = index - k_start_transmission_time;
	if (index < this->transmission_time_vec.size()) {
		return this->transmission_time_vec[index] * 299792458;
	}
	else {
		//navsol_object.set_calc_flag(false);
		cout << "range index  " << index << "\ttransmision time size: " << transmission_time_vec.size() << endl;
		return 0;
	}
}

void Channel::sett_flight(fp_precision norm) {
	if (this->t_flight == 0) {
		this->t_flight = norm / 299792458.0 - code_shift;
	}
}

vector<bool> Channel::getdecodedBits() { return decoded_bits; }

void Subframe::setup() {
	this->preamble = getBits(this->bits, 1, 8);
	this->frameID = getBits(this->bits, 50, 3);
	this->frameID_int = bin2dec(this->frameID);
	if (this->words.size() == 0) { this->words.resize(10); }

	for (int i = 0; i < 10; i++) {	//Every frame has 10 words
		vector<bool> word(30);
		word = getBits(this->bits, 30 * i + 1, 30);
		this->words[i] = word;
	}
}

vector<bool> Subframe::getParityBitsforWord(int n_word) { //n_word starts at 0
	vector<bool> parity_bits = subvector(this->bits, (n_word) * 30 + 25, 6);
	return parity_bits;
}

vector<bool> Subframe::CalculateParityBits(int n_word, vector<bool> last_2_bits) {

	vector<bool> calculated_parity_bits(6);
	calculated_parity_bits[0] = (last_2_bits[0] + words[n_word][0] + words[n_word][1] + words[n_word][2] + words[n_word][4] + words[n_word][5] + words[n_word][9] + words[n_word][10] + words[n_word][11] + words[n_word][12] + words[n_word][13] + words[n_word][16] + words[n_word][17] + words[n_word][19] + words[n_word][22])%2;
	calculated_parity_bits[1] = (last_2_bits[1] + words[n_word][1] + words[n_word][2] + words[n_word][3] + words[n_word][5] + words[n_word][6] + words[n_word][10] + words[n_word][11] + words[n_word][12] + words[n_word][13] + words[n_word][14] + words[n_word][17] + words[n_word][18] + words[n_word][20] + words[n_word][23])%2;
	calculated_parity_bits[2] = (last_2_bits[0] + words[n_word][0] + words[n_word][2] + words[n_word][3] + words[n_word][4] + words[n_word][6] + words[n_word][7] + words[n_word][11] + words[n_word][12] + words[n_word][13] + words[n_word][14] + words[n_word][15] + words[n_word][18] + words[n_word][19] + words[n_word][21])%2;
	calculated_parity_bits[3] = (last_2_bits[1] + words[n_word][1] + words[n_word][3] + words[n_word][4] + words[n_word][5] + words[n_word][7] + words[n_word][8] + words[n_word][12] + words[n_word][13] + words[n_word][14] + words[n_word][15] + words[n_word][16] + words[n_word][19] + words[n_word][20] + words[n_word][22])%2;
	calculated_parity_bits[4] = (last_2_bits[1] + words[n_word][0] + words[n_word][2] + words[n_word][4] + words[n_word][5] + words[n_word][6] + words[n_word][8] + words[n_word][9] + words[n_word][13] + words[n_word][14] + words[n_word][15] + words[n_word][16] + words[n_word][17] + words[n_word][20] + words[n_word][21] + words[n_word][23])%2;
	calculated_parity_bits[5] = (last_2_bits[0] + words[n_word][2] + words[n_word][4] + words[n_word][5] + words[n_word][7] + words[n_word][8] + words[n_word][9] + words[n_word][10] + words[n_word][12] + words[n_word][14] + words[n_word][18] + words[n_word][21] + words[n_word][22] + words[n_word][23]) %2;


	return calculated_parity_bits;
}

double Channel::receiver_time_shift = 0;
int Channel::n_channels_decoded = 0;
thread Channel::NavSolThread;
mutex Channel::navsolJoinMutex;
mutex Channel::receiver_time_shift_mutex;
mutex Channel::fft_mutex;

double predict_shift(RealArray Inputs)
{
	double variable_1 = Inputs[0];

	double variable_2 = Inputs[1];

	double variable_3 = Inputs[2];

	double variable_4 = Inputs[3];

	double variable_5 = Inputs[4];

	double variable_6 = Inputs[5];

	double variable_7 = Inputs[6];

	double variable_8 = Inputs[7];

	double variable_9 = Inputs[8];

	double variable_10 = Inputs[9];



	double scaled_variable_1 = (variable_1 + 0.0336529) / 0.0250872;

	double scaled_variable_2 = (variable_2 + 0.0245872) / 0.0217108;

	double scaled_variable_3 = (variable_3 - 0.0143649) / 0.0568902;

	double scaled_variable_4 = (variable_4 - 0.202) / 0.135375;

	double scaled_variable_5 = (variable_5 - 0.590248) / 0.176932;

	double scaled_variable_6 = (variable_6 - 0.924538) / 0.16001;

	double scaled_variable_7 = (variable_7 - 0.590505) / 0.177078;

	double scaled_variable_8 = (variable_8 - 0.202189) / 0.136049;

	double scaled_variable_9 = (variable_9 - 0.0142492) / 0.0569561;

	double scaled_variable_10 = (variable_10 + 0.0246912) / 0.021866;

	double y_1_1 = tanh(0.0103027 + (scaled_variable_1 * 0.034145) + (scaled_variable_2 * 0.0610118) + (scaled_variable_3 * -0.036235) + (scaled_variable_4 * -0.096343) + (scaled_variable_5 * -0.11895) + (scaled_variable_6 * 0.0171221) + (scaled_variable_7 * 0.137965) + (scaled_variable_8 * 0.108012) + (scaled_variable_9 * 0.0432789) + (scaled_variable_10 * 0.0220822));

	double y_1_2 = tanh(-0.277061 + (scaled_variable_1 * -0.0197795) + (scaled_variable_2 * -0.0189173) + (scaled_variable_3 * 0.0513916) + (scaled_variable_4 * 0.102931) + (scaled_variable_5 * 0.122833) + (scaled_variable_6 * -0.15521) + (scaled_variable_7 * -0.114664) + (scaled_variable_8 * -0.0308528) + (scaled_variable_9 * 0.0454386) + (scaled_variable_10 * -0.0193236));

	double y_1_3 = tanh(-0.34393 + (scaled_variable_1 * 0.00123173) + (scaled_variable_2 * 0.00410929) + (scaled_variable_3 * 0.0333932) + (scaled_variable_4 * -0.0256943) + (scaled_variable_5 * -0.122572) + (scaled_variable_6 * -0.155014) + (scaled_variable_7 * 0.158431) + (scaled_variable_8 * 0.127239) + (scaled_variable_9 * 0.0262074) + (scaled_variable_10 * 0.0162523));

	double y_1_4 = tanh(0.325237 + (scaled_variable_1 * 0.0133848) + (scaled_variable_2 * -0.032028) + (scaled_variable_3 * -0.0985864) + (scaled_variable_4 * -0.172562) + (scaled_variable_5 * -0.194299) + (scaled_variable_6 * 0.180506) + (scaled_variable_7 * 0.163639) + (scaled_variable_8 * 0.0610796) + (scaled_variable_9 * -0.00680306) + (scaled_variable_10 * 0.00310898));

	double y_1_5 = tanh(0.0689462 + (scaled_variable_1 * 0.0014953) + (scaled_variable_2 * 0.0697997) + (scaled_variable_3 * 0.0632471) + (scaled_variable_4 * 0.0927943) + (scaled_variable_5 * 0.114189) + (scaled_variable_6 * 0.0517562) + (scaled_variable_7 * -0.0863052) + (scaled_variable_8 * -0.0993178) + (scaled_variable_9 * -0.0659923) + (scaled_variable_10 * 0.0159323));

	double y_1_6 = tanh(0.329877 + (scaled_variable_1 * -0.0105319) + (scaled_variable_2 * -0.0079829) + (scaled_variable_3 * -0.0487476) + (scaled_variable_4 * 0.0409321) + (scaled_variable_5 * 0.143477) + (scaled_variable_6 * 0.172765) + (scaled_variable_7 * -0.140972) + (scaled_variable_8 * -0.127302) + (scaled_variable_9 * -0.0717496) + (scaled_variable_10 * 0.0472764));

	double y_1_7 = tanh(0.259699 + (scaled_variable_1 * -0.0425024) + (scaled_variable_2 * 0.0351423) + (scaled_variable_3 * -0.0942927) + (scaled_variable_4 * -0.109853) + (scaled_variable_5 * -0.0344637) + (scaled_variable_6 * 0.339436) + (scaled_variable_7 * -0.0281772) + (scaled_variable_8 * -0.107878) + (scaled_variable_9 * -0.114199) + (scaled_variable_10 * 0.0150968));

	double y_1_8 = tanh(-0.0844164 + (scaled_variable_1 * 0.00183961) + (scaled_variable_2 * -0.0517248) + (scaled_variable_3 * -0.0301948) + (scaled_variable_4 * -0.0297692) + (scaled_variable_5 * -0.0376777) + (scaled_variable_6 * -0.0894228) + (scaled_variable_7 * 0.0305809) + (scaled_variable_8 * 0.0602627) + (scaled_variable_9 * 0.0542816) + (scaled_variable_10 * -0.0250551));

	double y_1_9 = tanh(-0.355683 + (scaled_variable_1 * 0.0162187) + (scaled_variable_2 * -0.00916782) + (scaled_variable_3 * 0.0923017) + (scaled_variable_4 * 0.160392) + (scaled_variable_5 * 0.159833) + (scaled_variable_6 * -0.180647) + (scaled_variable_7 * -0.131761) + (scaled_variable_8 * -0.0683242) + (scaled_variable_9 * 0.00705421) + (scaled_variable_10 * 0.00289191));

	double y_1_10 = tanh(0.371377 + (scaled_variable_1 * -0.00587178) + (scaled_variable_2 * 0.00547511) + (scaled_variable_3 * -0.0266422) + (scaled_variable_4 * 0.0448954) + (scaled_variable_5 * 0.144715) + (scaled_variable_6 * 0.173567) + (scaled_variable_7 * -0.160896) + (scaled_variable_8 * -0.140279) + (scaled_variable_9 * -0.0341507) + (scaled_variable_10 * 0.00118841));

	double y_2_1 = tanh(-0.515753 + (y_1_1 * -0.17442) + (y_1_2 * 0.442767) + (y_1_3 * 0.0696405) + (y_1_4 * -0.513157) + (y_1_5 * 0.0903226) + (y_1_6 * -0.0899238) + (y_1_7 * -0.379632) + (y_1_8 * 0.0393267) + (y_1_9 * 0.581704) + (y_1_10 * -0.103306));

	double y_2_2 = tanh(-0.133852 + (y_1_1 * 0.0909715) + (y_1_2 * 0.00557959) + (y_1_3 * 0.0468805) + (y_1_4 * 0.00879384) + (y_1_5 * -0.0575361) + (y_1_6 * -0.0717132) + (y_1_7 * 0.0366168) + (y_1_8 * 0.0160666) + (y_1_9 * -0.027252) + (y_1_10 * -0.0330968));

	double y_2_3 = tanh(0.112987 + (y_1_1 * -0.0365617) + (y_1_2 * 7.13841e-05) + (y_1_3 * -0.000541131) + (y_1_4 * 0.0164929) + (y_1_5 * 0.00161918) + (y_1_6 * 0.0431421) + (y_1_7 * -0.0450221) + (y_1_8 * 0.0105385) + (y_1_9 * 0.0181873) + (y_1_10 * 0.00453016));

	double y_2_4 = tanh(0.0215439 + (y_1_1 * 0.0829699) + (y_1_2 * -0.0511604) + (y_1_3 * 0.0346931) + (y_1_4 * 0.0360999) + (y_1_5 * -0.0565902) + (y_1_6 * -0.0401599) + (y_1_7 * -0.014192) + (y_1_8 * 0.0143465) + (y_1_9 * -0.0771263) + (y_1_10 * -0.0148117));

	double y_2_5 = tanh(-0.106444 + (y_1_1 * -0.116992) + (y_1_2 * 0.0917893) + (y_1_3 * -0.0991975) + (y_1_4 * -0.0977002) + (y_1_5 * 0.101384) + (y_1_6 * 0.146897) + (y_1_7 * 0.095398) + (y_1_8 * -0.0253341) + (y_1_9 * 0.124232) + (y_1_10 * 0.0890511));

	double y_2_6 = tanh(-0.112415 + (y_1_1 * -0.0740131) + (y_1_2 * 0.0357311) + (y_1_3 * 0.0379817) + (y_1_4 * -0.0927012) + (y_1_5 * 0.0435112) + (y_1_6 * -0.062825) + (y_1_7 * -0.0468188) + (y_1_8 * -0.0308812) + (y_1_9 * 0.166953) + (y_1_10 * -0.054236));

	double y_2_7 = tanh(0.0172275 + (y_1_1 * 0.0217104) + (y_1_2 * -0.0140042) + (y_1_3 * 0.0215908) + (y_1_4 * 0.0309656) + (y_1_5 * -0.0303642) + (y_1_6 * -0.00977788) + (y_1_7 * -0.0108599) + (y_1_8 * 0.0100145) + (y_1_9 * -0.0271988) + (y_1_10 * -0.0198237));

	double y_2_8 = tanh(0.0317485 + (y_1_1 * 0.054853) + (y_1_2 * -0.0419356) + (y_1_3 * 0.0417897) + (y_1_4 * 0.0542852) + (y_1_5 * -0.0576472) + (y_1_6 * -0.0220825) + (y_1_7 * -0.0352424) + (y_1_8 * 0.0260448) + (y_1_9 * -0.0643854) + (y_1_10 * -0.0316045));

	double y_2_9 = tanh(0.0410604 + (y_1_1 * -0.00848072) + (y_1_2 * -0.00833515) + (y_1_3 * -0.00624634) + (y_1_4 * 0.0119827) + (y_1_5 * 0.00577651) + (y_1_6 * 0.0134792) + (y_1_7 * -0.00941924) + (y_1_8 * 0.0005248) + (y_1_9 * -0.0055533) + (y_1_10 * 0.00448629));

	double y_2_10 = tanh(-0.472503 + (y_1_1 * 0.165291) + (y_1_2 * 0.0288089) + (y_1_3 * 0.50445) + (y_1_4 * -0.0988261) + (y_1_5 * -0.209268) + (y_1_6 * -0.466544) + (y_1_7 * -0.348859) + (y_1_8 * 0.129309) + (y_1_9 * 0.0731362) + (y_1_10 * -0.522845));

	double scaled_variable_11 = (0.0106483 + (y_2_1 * -0.793459) + (y_2_2 * -0.193532) + (y_2_3 * -0.00236684) + (y_2_4 * -0.172372) + (y_2_5 * 0.194053) + (y_2_6 * 0.171904) + (y_2_7 * -0.0934002) + (y_2_8 * -0.209814) + (y_2_9 * 0.000357772) + (y_2_10 * 0.765557));

	double variable_11 = (0.5 * (scaled_variable_11 + 1.0) * (1 + 1) - 1);

	return variable_11;

}
