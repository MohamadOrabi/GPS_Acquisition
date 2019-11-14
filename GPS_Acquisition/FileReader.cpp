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
#include "parameters.h"

fp_precision global_duration = 292;//292;// 100;
fp_precision global_fs = 2.5;
string fileName = "Test3_2p5M.bin";// "Test1_5MHz.bin"	"Test3_2p5M.bin"	"Kimia_GNSS_GPS__test1.bin"		GPS-22-7-2019-20Msps.bin
FileReader file_reader_object(fileName, global_fs * 1E3, global_fs, global_duration);
//fp_precision skip_time = 40; //in seconds

FileReader::FileReader(string fileName, int nsamples, double fs, double duration) {

	this->Channel_list.clear();
	this->fileName = fileName;
	this->file = ifstream(fileName, ios::in | ios::binary);
	this->nsamples = nsamples;
	this->sampledData_vector = vector<SampleIQData>(nsamples);
	this->Data_vector = VEC(nsamples);
	this->fs = fs;
	this->duration = duration;
	this->shouldReadData = false;

	this->duration_in_codelengths = round(duration * (fs * 1E6) / nsamples);
	
	if (file.is_open()) {
		cout << fileName << " was opened successfully" << endl;
		cout << fileName << " length is: " << this->getFileLength() << endl;
		//file.seekg(0, (skip_time * 1e3) * 2 * 2 * sampledData_vector.size());
	}
	else cout << fileName << " failed to open" << endl;
}

VEC FileReader::gatherData() {
	if (file.is_open()) {
		int k = 0;
		while (k <= duration_in_codelengths) {
			if (!file.good()) { cout << "Reached end of file at: "<< k/1E3 << " (probably)" << endl; }
			if (Channel_list.size() > 0 && shouldReadData_check()) {
				file.read((char*)& sampledData_vector[0], 2 * 2 * sampledData_vector.size()); //2*2 is because complex double
				//this->Data_vector = VEC(nsamples);
				SampleIQDataToIQArray(sampledData_vector, Data_vector);
				Data_vectors.push_back(Data_vector);
				//auto shared_ptr = make_shared<VEC>(&Data_vectors[k]);
				SampleFrame_Ptr p = std::make_shared<SampleFrame>();
				p->sampleArray = std::move(Data_vectors[k]);
				p->current_index = k * nsamples;
				p->k = k;

				for (int i = 0; i < nchannels; i++) {
					(*Channel_list[i]).Q_mutex.lock();
					(*Channel_list[i]).addToWorkQ(p);
					(*Channel_list[i]).Q_mutex.unlock();
				}
				k++;
			}
		}
	}
	return Data_vector;
}

int FileReader::getFileLength() {
	int tempg = file.tellg();	//Saves current g pointer

	file.seekg(0, file.end);
	int file_length = file.tellg();
	file.seekg(tempg, file.beg);	//Restors g pointer

	return file_length;
}

void FileReader::addChannelToList(Channel* channel_to_add) {
	Channel_list.push_back(channel_to_add);
	nchannels++;
}

void FileReader::SampleIQDataToIQArray(vector<SampleIQData>& A1, VEC& A2) {
	for (int i = 0; i < A1.size(); i++) {
		A2[i] = A1[i];
	}
}

bool FileReader::shouldReadData_check() {	
	for (int i = 0; i < Channel_list.size(); i++) {
		if ((*Channel_list[i]).getNeedDataFlag() == true) {
			return true;
		}
	}
	return false;
}
