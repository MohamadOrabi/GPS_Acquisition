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

//Code copied from matrix project so we follow the same convention



using fp_precision = double;

using IQData = std::complex<fp_precision>;
using SampleIQData = std::complex<short>;

//using SampleArray = Eigen::Array<SampleIQData, Eigen::Dynamic, 1>;
using IQArray = Eigen::Array<IQData, Eigen::Dynamic, 1>;
using RealArray = Eigen::Array<fp_precision, Eigen::Dynamic, 1>;

using VEC = IQArray;
using MAT = Eigen::Matrix<IQData, Eigen::Dynamic, Eigen::Dynamic>;
using RealMAT = Eigen::Matrix<fp_precision, Eigen::Dynamic, Eigen::Dynamic>;
using RealVEC = Eigen::VectorXd;
using fftPlan = fftw_plan;
using fftComplex = fftw_complex;

struct  SampleFrame {
	VEC sampleArray; //!< The actual I samples
	int current_index;
	int k;
	SampleFrame() = default;

	/**
	* \brief Copy Constructor Disabled
	*
	* Copies over the sample array and the sequence number from the other sample frame.
	* \param other The target sample frame from which to copy from.
	*/
	SampleFrame(SampleFrame& other) = delete;

	// Assignment is disabled.
	SampleFrame& operator=(const SampleFrame& sf) = delete;
};

// This is the data structure that we pass around.
using SampleFrame_Ptr = std::shared_ptr<SampleFrame>;

// This is a work queue, from here data is dequeued.
using WorkQueue = moodycamel::ConcurrentQueue<SampleFrame_Ptr>;

