#pragma once
#include "stdafx.h"
#include "TypeDefs.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <complex>
#include <Eigen/Dense>
#include <cmath>

#include <queue>
#include <deque>

#include "GPS_Utils.h"
#include "Channel.h"

using namespace std;

class FileReader {

private:

	queue<VEC> reading_queue;
	string fileName;
	ifstream file;

	int nsamples;		//nsamples is the number of samples to read at a time 
	vector<SampleIQData> sampledData_vector;	//vector storing the last read Data NOT VEC
	VEC Data_vector;	//vector storing the last read Data VEC
	vector<VEC> Data_vectors;

	int nchannels = 0;	//Number of channels the FileReader is adding to their work queues
	vector<Channel*> Channel_list; //List of Channels the FileReader is adding to their work queues
	double fs;
	double duration;
	int duration_in_codelengths;

	bool shouldReadData;

	void SampleIQDataToIQArray(vector<SampleIQData>& A1, VEC& A2); //Method used internally to convert from SampleIQData vector to IQData vector

	bool shouldReadData_check();

public:
	//bool should_gather_data = true;
	FileReader(string fileName, int nsamples, double fs, double duration);

	VEC gatherData();

	int getFileLength();

	void addChannelToList(Channel* channel_to_add);
};

extern FileReader file_reader_object;
extern fp_precision global_duration;
extern fp_precision global_fs;
//extern fp_precision global_skip_time;
extern string fileName;
