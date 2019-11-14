#define _USE_MATH_DEFINES

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
#include <thread>
#include <chrono>
#include <iomanip>


#include "TypeDefs.h"
#include "GPS_Utils.h"
#include "Channel.h"
#include "FileReader.h"
#include "parameters.h"


//#include "GPS_functions.h"

using namespace std;

const int code_length = 1023;	// 1023 is a constant for GPS
vector<char> PRNS;		//PRNS is vector that holds all the PRNS

RealArray getPRN(const int sat_id) {
	RealArray PRN_ID(code_length);

	int index = (sat_id - 1) * code_length;

	for (int i = 0; i < code_length; i++) {
		PRN_ID[i] = PRNS[i + index];
	}

	return PRN_ID;
}

//This function prints arrays for debugging purposes
template <typename T>
void printArray(T Array_to_print) {
	for (int i = 0; i < Array_to_print.size(); i++) {
		cout << i << ": " << Array_to_print[i] << "\t";
	}
	cout << endl << endl;
}

void SampleIQDataToIQArray(vector<SampleIQData>& A1, IQArray& A2) {
	for (int i = 0; i < A1.size(); i++) {
		A2[i] = A1[i];
	}
}

template <typename T>
void loadBinaryToVector(vector<T>& Data_vector, const string FileName, int nread, const int current_index) {	//if nread = 0, then vector is resized to length of file .. nread IS IN BYTES
	ifstream file(FileName, ios::in | ios::binary);

	//If file was found and opened
	if (file.is_open()) {
		//cout << FileName << " File Opened!" << endl;

		//Getting file length
		file.seekg(0, file.end);
		int file_length = file.tellg();
		file.seekg(current_index, file.beg);
		//cout << "Length of file is: " << file_length << endl << endl << endl;	//Print file length
		if (nread == 0) {
			nread = file_length;
			Data_vector.resize(nread);			//It only resizes the PRN vector .. when nread = 0
		}

		file.read((char*)& Data_vector[0], nread);


		/*
		cout << FileName << endl;
		cout << "current index was: " << current_index << endl;
		cout << "tellg: " << file.tellg() << endl;
		cout << "Data type size: " << sizeof(T) << endl;
		*/

	//If file didn't open (ex. was not found)
	}
	else {
		cout << FileName << " File didn't open!" << endl;
	}
}

void loadData(const string FileName, IQArray& IQ_Data, const int nsamples, const int current_index) {
	vector<SampleIQData> Data(nsamples);
	loadBinaryToVector<SampleIQData>(Data, FileName, 2 * 2 * Data.size(), current_index); // *2*2 is because 4 bytes for complex
	
	SampleIQDataToIQArray(Data, IQ_Data);
}

void printPRN(const int satID) {
	RealArray PRN_to_print = getPRN(satID);	//Use getPRN function to get specific PRN

	cout << "PRN#" << satID << ": ";
	for (int i = 0; i < code_length; i++) {
		cout << (int)PRN_to_print[i];
	}
	cout << endl;
}


//**************************************************
//This is where main starts
//**************************************************

int main() {
	/*vector<bool> batata{0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0};
	fp_precision batata_int = twosCompliment2Dec(batata);
	cout << "batata: " << batata_int << endl;*/

	const string PRNFileName = "PRNSFile.bin";	//Location of file to get PRNS from
	//const string DataFileName = "Test3_2p5M.bin"; //Location of file to get Data from
	//const string DataFileName = "Test1_5MHz.bin";
	const string DataFileName = fileName;
	const double duration = global_duration;	//Duration to run tracking for

	const fp_precision fs = global_fs;		//Sampling frequency IN MHZ
	const fp_precision chip_rate = 1.023; //Chip_rate of sattelite IN MHZ
	const int naccumulations = 5;	//Number of accumulations for acquisition (Maybe tracking too)

	const int nsamples =(int) (fs * code_length / chip_rate);	//Number of samples in one code length
	const fp_precision Tsub = nsamples / (fs * 1E6);	//Tsub of the loop filter 

	const fp_precision f_lo = -6E3;	//In Hz
	const fp_precision f_hi = 6E3;
	const fp_precision f_step = 25;
	RealArray f_array = creatRealArray(f_lo, f_step, f_hi);
	RealArray t_array = creatRealArray(0, 1, nsamples - 1) / (fs*1E6);
	
	//Initialization of PRNS vector
	loadBinaryToVector<char>(PRNS, PRNFileName, 0, 0);	//Loads PRNS from file

	//******************************************************************************************************************************************************************
	//1,15,18,21,22,24,26,27,29
	//vector<int> sats{1,14,22,31};
	//vector<int> sats{ 1,11,14,26};	//This gives correct navigation solution
	//vector<int> sats = {1,3,10,11,14,22,31}; // 3
	vector<int> sats = { 1,3,11,14,22,31 }; // 3

	//vector<int> sats = {1, 11, 14, 26};
	//vector<int> sats{26}; //3
	//vector<int> sats = creatIntVector(1, 1, 32);

	//vector<int> detectedsats;

	vector<RealArray> PRNforSat;
	vector<fp_precision> fdforSat;
	vector<fp_precision> code_shiftforSat;
	vector<bool> detectedforSat;
	int ndetected = 0;
	fp_precision last_t = 0;

	const float threshold = 0;	//Threshold for detection

	auto startAcquisition = std::chrono::high_resolution_clock::now();

	const int noise_PRN = 37;	//This PRN is used to get an estimate of the noise power
	bool Noise_run = true;
	fp_precision N_Powerdb = 0;

	for (int i = 0; i < sats.size(); i++) {
		if (Noise_run) {
			auto it = sats.begin();
			sats.insert(it, noise_PRN);
		}
		PRNforSat.push_back(2 * getPRN(sats[i]) - 1); // Converting 0,1 to -1,1

		VEC IQ_Data(nsamples);		//Complex vector storing Data in I +jQ
		VEC IQ_Data_fft(nsamples);		//Data fft
		VEC resampled_PRN(nsamples);	//Resampled PRN
		VEC resampled_PRN_fft(nsamples);	//Resampled PRN fft
		VEC corr(nsamples);					//Correlation Result
		RealArray corr_power(nsamples);			//Abs2/nsamples of correlation

		//Start Acquisition 
		Eigen::Array<fp_precision, Eigen::Dynamic, Eigen::Dynamic> corrs;
		corrs.setZero(nsamples, f_array.size());
		int corr_index = 0;
		int corr_doppler = 0;
		fp_precision corr_max = 0, corr_mean = 0;

		//cout << "Satellite #" << sats[i];
		for (int k = 0; k < naccumulations; k++) {
			cout << "\rSatellite #" << sats[i]<< "\tAccumulation #" << k + 1;
			//cout << "\tAccumulation #" << k + 1 << endl;

			//Load new data and resample PRNS
			int current_index = k * nsamples;
			loadData(DataFileName, IQ_Data, nsamples, current_index * sizeof(SampleIQData)); //Load IQData from file
			resampled_PRN = resamplePRN(PRNforSat[i], nsamples, chip_rate, fs, code_length, current_index);
			conduct_FFT(resampled_PRN, resampled_PRN_fft.data());

			//Search for frequency
			for (int findex = 0; findex < f_array.size(); findex++) {
				fp_precision df = f_array[findex];	//Getting doppler for current iteration

				double initialphase = 0.00;
				VEC mixed_Data(nsamples);
				RealArray t_array = creatRealArray(last_t, 1, last_t + nsamples - 1) / (fs * 1E6);
				last_t = t_array[nsamples - 1];	//Save the last value of t_array
				mixed_Data = mixSignal(IQ_Data, fs, df, &initialphase); //IQ_Data .* exp(2*pi*df*t + initial_theta)

				//Conducting cross correlation
				conduct_FFT(mixed_Data, IQ_Data_fft.data());
				corr = Corr(IQ_Data_fft, resampled_PRN_fft);
				corr_power = corr.cwiseAbs2();
				corrs.col(findex) += corr_power / naccumulations;	//This is for accumulation
			}
		}
		cout << endl;


		corr_max = corrs.maxCoeff(&corr_index, &corr_doppler);

		//corr_doppler -= 1;
		corr_mean = corrs.sum() / (corrs.cols() * corrs.rows());

		fp_precision corr_maxdb = 10 * log10(corr_max);
		fp_precision corr_meandb = 10 * log10(corr_mean);

		if (Noise_run) {
			N_Power = corr_mean;
			N_Powerdb = corr_meandb;
			cout << "\nNoise Power = " << N_Powerdb << endl;
			Noise_run = false;
			auto it = sats.begin();
			sats.erase(it);
			PRNforSat.clear();
			i -= 1;
		}
		else {
			//exportArray("corrs_" + to_string(sats[i]) + ".csv", corrs);	//For Debugging, takes a lot of time!

			fp_precision CNR =10*log10((corr_max - N_Power) / (N_Power * nsamples / (fs * 1E6)));

			if ((CNR) > threshold) { detectedforSat.push_back(true); ndetected += 1; }
			else { detectedforSat.push_back(false); }
			cout << "Max power(dB): \t\t" << corr_maxdb << "\nMean power(dB): \t" << corr_meandb << "\ncode shift: \t\t" << corr_index << endl;
			//fp_precision CNR = (std::pow(corr_maxdb, 2) - 2 * std::pow(corr_meandb, 2)) * 15e3 / (2 * std::pow(corr_meandb, 2)); //CNR in Matrix??
			//CNR = 10 * log10(CNR);
			cout << "CNR(dB): \t\t" << CNR << endl << "Max(dB) - Mean(dB): \t" << corr_maxdb - corr_meandb << endl;

			//OUTPUTS:
			fdforSat.push_back(f_array[corr_doppler]);
			code_shiftforSat.push_back((corr_index) / (fs * 1E6));
			cout << "fd: " << fdforSat[i] << endl << endl;
			cout << "Detected: " << detectedforSat[i] << endl << endl;

		}

		

	}
	auto finishAcquisition = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedAcquisition = finishAcquisition - startAcquisition;
	std::cout << "Acqusition time: " << elapsedAcquisition.count() << " s\n" << endl;


	//**********************************************************************************************************
	//Tracking																								   /
	//**********************************************************************************************************

	//Channel Paramteres
	const int loop_order_PLL = 2;
	const fp_precision BPLL = 40;
	//const fp_precision BPLL = 50;
	const fp_precision BDLL = 0.005;

	int ntrack = duration*1000;	//codelengths to track for (duration * 1000)
	int pp_moving_average = 1000;	//Moving average for carrier power

	//Initializing channels and threads vectors
	vector<Channel*> channels;
	channels.clear();
	vector<thread> threads;
	threads.clear();

	//Creating channels and pushing them to channels vector
	cout << "Creating channels for detected satellites: ";
	for (int i = 0; i < sats.size(); i++) {
		if (detectedforSat[i] == true) {
			cout << sats[i] << " ";
			Channel* chan = new Channel(PRNforSat[i], sats[i], fs, chip_rate, nsamples, loop_order_PLL,
				BPLL, BDLL, fdforSat[i], code_shiftforSat[i], pp_moving_average, ntrack);
			//cout << "fdforSat: " << fdforSat[i] << endl;
			//cout << "sats: " << sats[i] << endl;
			//cout << "code_shiftforSat: " << code_shiftforSat[i]*fs*1E6 << endl;

			channels.push_back(chan);
		}
	}
	cout << endl;
	//Start the file reading thread
	//std::this_thread::sleep_for(6s);
	thread FileReaderThread(&FileReader::gatherData, &file_reader_object);
	cout << "File reader thread done" << endl;

	//std::this_thread::sleep_for(1s);
	auto startTracking = std::chrono::high_resolution_clock::now();

	//Start tracking threads for all channels
	for (int i = 0; i < ndetected; i++) {
		threads.push_back(thread(&Channel::track, (channels[i])));
	}

	//Wait for all threads to finish before exiting
	FileReaderThread.join();	//This should be removed so that reading and tracking happen at the same time
	for (auto& th : threads) {
		th.join();
	}

	auto finishTracking = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedTracking = finishTracking - startTracking;
	std::cout << "Tracking time: " << elapsedTracking.count() << " s\n" << endl;

	//Save vectors for all channels for debugginf purposes
	int index = 0;

	//for (int i = 0; i < sats.size(); i++) {
	//	if (detectedforSat[i] == true) {
	//		exportArray("Channel_Data_" + to_string(sats[i]) + ".csv", (*channels[index]).exportChannelData());

	//		exportVEC("pp_" + to_string(sats[i]) + ".csv", (*channels[index]).getpp());
	//		exportVEC("ee_" + to_string(sats[i]) + ".csv", (*channels[index]).getee());
	//		exportVEC("ll_" + to_string(sats[i]) + ".csv", (*channels[index]).getll());
	//		exportVEC("fd_" + to_string(sats[i]) + ".csv", (*channels[index]).getfd_vec());
	//		exportVEC("phase_error_" + to_string(sats[i]) + ".csv", (*channels[index]).getphase_error());
	//		exportBits("decoded_bits_" + to_string(sats[i]) + ".csv", (*channels[index]).getdecodedBits());

	//		//exportvec("debug_vec_" + to_string(sats[i]) + ".csv", (*channels[index]).getdebug_vec());
	//		exportvec("transmission_time_" + to_string(sats[i]) + ".csv", (*channels[index]).gettransmission_time_vec());
	//		exportvec("x_kvec_" + to_string(sats[i]) + ".csv", (*channels[index]).getx_kvec());
	//		exportvec("y_kvec_" + to_string(sats[i]) + ".csv", (*channels[index]).gety_kvec());
	//		exportvec("z_kvec_" + to_string(sats[i]) + ".csv", (*channels[index]).getz_kvec());

	//		exportVEC("I_" + to_string(sats[i]) + ".csv", real((*channels[index]).getp()));
	//		exportVEC("Q_" + to_string(sats[i]) + ".csv", imag((*channels[index]).getp()));

	//		index += 1;
	//	}
	//}

	cout << "\n\n" << "Press Enter to Exit" << endl;
	getchar();
	return 0;
}


