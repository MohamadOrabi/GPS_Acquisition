#include "stdafx.h"
#include "TypeDefs.h"
//#include <fftw3.h>
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

NavSol navsol_object;

NavSol::NavSol() {
	this->current_index = 0;
	this->time_step = 100; //Current index step (in codelengths)
	this->threshold = 1E-3;
	this->calc_flag = true;
	this->initialization = true;
}

void NavSol::Calculate_Position() {
	GoogleEarthPath receiverPath("KML_navsol.kml", "Receiver");
	while ((*Channels[0]).done_flag == false) {
		update_calc_flag();
		if (calc_flag && Channels.size() >= 4) {		

			Eigen::Matrix<fp_precision, Eigen::Dynamic, 4> H_k;
			Eigen::Matrix<fp_precision, Eigen::Dynamic, Eigen::Dynamic> variance;
			Eigen::Matrix<fp_precision, 4, 1> deltax;
			Eigen::Matrix<fp_precision, Eigen::Dynamic, 1> h_k;
			Eigen::Matrix<fp_precision, Eigen::Dynamic, 1> range;

			fp_precision error = 100;
			int count = 1;

			vector<fp_precision> xs;
			vector<fp_precision> ys;
			vector<fp_precision> zs;
			vector<fp_precision>range_vec;
			vector<fp_precision> variance_vec;

			//Getting Satellite positions
			for (int i = 0; i < Channels.size(); i++) {
				//cout << "Channel:\t" << i << "for current_index: " << current_index << endl;
				(*Channels[i]).vec_mutex.lock();	//This line could easily be moved to channel.cpp
				if ((*Channels[i]).getx_kat(current_index) != 0) {

					xs.push_back((*Channels[i]).getx_kat(current_index));
					ys.push_back((*Channels[i]).gety_kat(current_index));
					zs.push_back((*Channels[i]).getz_kat(current_index));

					//range.conservativeResize(xs.size(), 1);
					//cout << "range at i: " << xs.size() << "\t" << (*Channels[i]).getrange(current_index) << endl;;
					//range(i) = (*Channels[i]).getrange(current_index);

					/*variance.conservativeResize(xs.size(), xs.size());
					variance(i,i) = (*Channels[i]).getC_Nat(current_index);*/
					variance_vec.push_back((*Channels[i]).getVarianceat(current_index));
					range_vec.push_back((*Channels[i]).getrange(current_index));// +t_flight_codelengths[i]);

				
					
					//range_vec.push_back((*Channels[i]).getrange(current_index));
					//TODO: Some error was happening when range was here and there was not time_step (time_step = 1)
				}
				(*Channels[i]).vec_mutex.unlock();

				
			}

			if (xs.size() >= 4) {
				H_k.resize(xs.size(), 4);
				h_k.resize(xs.size(), 1);
				variance.resize(xs.size(), xs.size());
				range.resize(xs.size(), 1);
				variance.setZero();
				//range.resize(xs.size(), 1);

				//Initializing estimated position
				if (initialization) {
					ex = 0; ey = 0; ez = 0; b = 0;
					fp_precision size = xs.size();

					for (int i = 0; i < xs.size(); i++) {
						ex += xs[i] / size;
						ey += ys[i] / size;
						ez += zs[i] / size;
					}

					fp_precision norm = sqrt(pow(ex, 2) + pow(ey, 2) + pow(ez, 2));
					ex *= 6371E3 / norm;
					ey *= 6371E3 / norm;
					ez *= 6371E3 / norm;

					sol(0) = ex;
					sol(1) = ey;
					sol(2) = ez;
					//sol(3) = b;

					//TEST
					/*for (int i = 0; i < xs.size(); i++) {
						fp_precision norm = sqrt(pow(ex - xs[i], 2) + pow(ey - ys[i], 2) + pow(ez - zs[i], 2));
						fp_precision flight_time = norm / 299792458;
						flight_time -= fmod(flight_time, 1E-3);
						range_vec[i] += flight_time * 299792458;
						(*Channels[i]).sett_flight(norm);
					}*/

					for (int i = 0; i < xs.size(); i++) {
						fp_precision calc_range = sqrt(pow(ex - xs[i], 2) + pow(ey - ys[i], 2) + pow(ez - zs[i], 2));
						b += (range_vec[i] - calc_range) / size;
					}

					sol(3) = b;

					initialization = false;
				}
				
				bool first_run = true;	//flag used to populate the variance and range matrices

				//Estimating position from pseudoranges using Newton Iteration
				while (abs(error) > threshold && count < 20) {

					for (int i = 0; i < xs.size(); i++) {
						fp_precision norm = sqrt(pow(ex - xs[i], 2) + pow(ey - ys[i], 2) + pow(ez - zs[i], 2));
						H_k(i, 0) = (ex - xs[i]) / norm;
						H_k(i, 1) = (ey - ys[i]) / norm;
						H_k(i, 2) = (ez - zs[i]) / norm;
						H_k(i, 3) = 1.0;
						h_k(i, 0) = norm + b;
						if (first_run) {
							//variance(i, i) = 2 + 1000*pow(299792458,2)*(0.005*0.8*pow(1/1.023E6,2)/(2*variance_vec[i]) * (1 + 2/(4E-3*variance_vec[i])));
							variance(i, i) = variance_vec[i];
							range(i) = range_vec[i];
						}
					}
					//cout << variance << endl;

					//cout << variance << endl;

					//Least Sqaures Estimator
					//deltax = (H_k.transpose() * H_k).inverse() * H_k.transpose() * (range - h_k);	//Check if inverse() does its job here
					//Maximum Likelihood Estimator
					deltax = (H_k.transpose()* variance.inverse() * H_k).inverse() * H_k.transpose()* variance.inverse() * (range - h_k);	//Check if inverse() does its job here

					sol += deltax;
					error = deltax.norm();

					ex = sol(0);
					ey = sol(1);
					ez = sol(2);
					b = sol(3);

					count++;
					first_run = false;
					//cout << "count:\t"<< count << endl;
				}

				//Checking for error after Newton Iteration (Helps with debugging)
				if (abs(error) > threshold) {
					cout << "Navsol Error is large: " << error << endl;
				}
				if (isnan(error)) {
					//initialization = true;
					cout << "Navsol Error is NAN!" << endl;
				}

				//Stroing estimates in vectors to be exported
				ex_vec.push_back(ex);
				ey_vec.push_back(ey);
				ez_vec.push_back(ez);
				b_vec.push_back(b);

				for (int i = 0; i < Channels.size(); i++) {
					fp_precision norm = sqrt(pow(ex - xs[i], 2) + pow(ey - ys[i], 2) + pow(ez - zs[i], 2));
					if (t_flight_codelengths[i] == 0) {
						t_flight_codelengths[i] = round(norm / 299792458/1E-3);
						t_flight_codelengths[i] *= 299792458;
						//initialization = true;

					}
					(*Channels[i]).sett_flight(norm);
				}




				//cout << variance << endl;

				fp_precision lon = 0; fp_precision lat = 0; fp_precision alt = 0;
				ecef2lla(ex, ey, ez, &lat, &lon, &alt);
				receiverPath.addPoint(lon, lat);
			
			}
			current_index += time_step;
		}
		else {
			//update_calc_flag();
		}
	}
	//this->exportNavSol();
}

void NavSol::AddToChannels(Channel* channel_to_add, int index) {
	AddToChannelsMutex.lock();	//Prevent Data Race
	//TODO: Maybe add a channelsToBeAdded vector!
	this->Channels.push_back(channel_to_add);	//TODO: Maybe add a mutex here to prevent writing to channels while trying to get data from it!
	this->t_flight_codelengths.push_back(0);

	if (Channels.size() <= 4) {
		if (index > current_index) { current_index = index; }	//TODO: this will cause problems when adding channels and channels.size() > 4
																//I think this was fixed by if statement at line 60
	}
	AddToChannelsMutex.unlock();
}

void NavSol::update_calc_flag() {
	calc_flag = true;
	for (int i = 0; i < Channels.size() && calc_flag; i++) {
		bool available = (*Channels[i]).check_nav_index(current_index);
		if (!available) {
			calc_flag = false;
			//cout << "Navsol is waiting for data!" << endl;
		}
	}
}

void NavSol::set_calc_flag(bool flag) {
	calc_flag = flag;
}

void NavSol::exportNavSol() {

	Eigen::Matrix<fp_precision, Eigen::Dynamic, 4> NavSolData;
	NavSolData.resize(ex_vec.size(), 4);

	NavSolData.col(1) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(ex_vec.data(), ex_vec.size());
	NavSolData.col(2) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(ey_vec.data(), ey_vec.size());
	NavSolData.col(3) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(ez_vec.data(), ez_vec.size());
	NavSolData.col(4) = Eigen::Map<Eigen::Matrix<fp_precision, Eigen::Dynamic, 1>>(b_vec.data(), b_vec.size());

	exportArray("NavSol.csv", NavSolData);
	cout << "exported" << endl;
}

void NavSol::setInitialization(bool boolean) {
	this->initialization = true;
}
