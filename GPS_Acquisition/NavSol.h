#pragma once
#include "stdafx.h"
#include "TypeDefs.h"
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <cmath>
#include <queue>



#include "GPS_Utils.h"
#include "Channel.h"
#include "FileReader.h"

class NavSol {
private:
	vector<Channel*> Channels;
	bool calc_flag;
	bool initialization;
	int current_index;
	int time_step;

	fp_precision threshold;
	fp_precision ex;
	fp_precision ey;
	fp_precision ez;
	fp_precision b;

	vector<fp_precision> ex_vec;
	vector<fp_precision> ey_vec;
	vector<fp_precision> ez_vec;
	vector<fp_precision> b_vec;
	vector<fp_precision> t_flight_codelengths;
	Eigen::Matrix<fp_precision, 4, 1> sol;

	mutex AddToChannelsMutex;
	mutex ChangingChannelsMutex;

public:
	
	NavSol();

	void Calculate_Position();
	void AddToChannels(Channel* channel_to_add, int index);
	void update_calc_flag();
	void set_calc_flag(bool flag);

	void exportNavSol();
	void setInitialization(bool boolean);

};

extern NavSol navsol_object;