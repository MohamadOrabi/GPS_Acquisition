#pragma once
#define _USE_MATH_DEFINES

#include "stdafx.h"
#include "TypeDefs.h"
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <cmath>
#include <queue>
#include <mutex>
#include <optional>

#include "GPS_Utils.h"
using namespace std;

struct Subframe {
	vector<bool> bits;	//This contains the recieved bits 
	vector<bool> frameID;
	vector<bool> preamble;
	vector<bool> parity_bits;

	vector<vector<bool>> words;	//This contains the actual data bits
	vector<bool> error_flags;
	int frameID_int;
	int index;
	//bool flipped = false;

	void setup();
	vector<bool> getParityBitsforWord(int n_word);
	vector<bool> CalculateParityBits(int n_word, vector<bool> last_2_bits);
};

class Channel {
private:
	int k = 0; //Tracking iteration
	int dec_acc, n_acc_PLL, n_acc_DLL, track_window, ek_length;
	fp_precision tracking_threshold;
	fp_precision channel_time;
	//queue<VEC> workQ;
	queue<SampleFrame_Ptr> workQ;

	//Flags
	vector<bool> sf_decoded;
	const int sat_ID;	//ID of sattelit to track 
	const fp_precision fs;		//Sampling frequency IN MHZ

	const fp_precision chip_rate; //Chip_rate of sattelite IN MHZ
	const int code_length = 1023;

	const int nsamples;	//Number of samples in one code length
	const fp_precision Tsub;	//Tsub of the loop filter 

	fp_precision C_N;
	fp_precision variance;

	//Tracking
	const int loop_order_PLL;
	const fp_precision BPLL, BDLL;
	fp_precision fd, code_shift;
	RealArray PRN, fd_vec;
	vector<bool> decoded_bits;

	//Configure Loop Filter
	PLL loop_filter_PLL = configureLoopFilter(loop_order_PLL, BPLL, Tsub, fd);

	//PLL Discriminator 
	fp_precision last_t = 0;	//Stores the last value in t_array
	fp_precision last_theta = 0;
	int delta_shift;		//This sets the early and late

	const int pp_moving_average;

	VEC IQ_Data;

	
	//These are to be removed (Just for debugging), they don't need to be lists
	VEC p;
	VEC e;
	VEC l;
	VEC pp;
	VEC ee;
	VEC ll;
	VEC phase_error;
	IQData pk;
	vector<fp_precision> debug_vec;
	vector<fp_precision> x_kvec, y_kvec, z_kvec, variance_vec, transmission_time_vec;

	bool is_decoding = false, is_tracking = false;
	//Satellite coordinates
	fp_precision x_k, y_k, z_k;

	//For Decoding
	const vector<bool> preamble1{ 1,0,0,0,1,0,1,1 };	//First 8 bits of TLM subframe
	const vector<bool> preamble2{ 0,1,1,1,0,1,0,0 };	//First 8 bits of TLM subframe
	const int n_preamble_check = 4;
	//bool flip = false;
	bool preambles_detected = false;	//Consecutive preambles detected flag
	bool decoded_ephemeris = false;
	bool decoded_clock_correction = false;
	int index_of_first_preamble = 0;

	int n_bits_transmitted = 0;
	//int n_CA_reapeats = 0;
	int n_CA_chips = 0;
	int n_chip_fractions = 0;
	//double time_since_detection = 0;
	int k_start_decoding;
	int k_start_transmission_time = 0;

	vector<Subframe> frames;	//Holds n_preamble_check frames that were detected by check_consecutive_preamble
	Subframe current_frame;		//Holds last frame the was read from the binary file
	Subframe past_frame;		//Hold the past frame that was processed and decoded from the binary file

	const fp_precision F = -4.442807633E-10;
	fp_precision transmission_time, deltat_sv;
	fp_precision t_OC, a_f2, a_f1, a_f0, E_k;	//Parameters for clock correction
	fp_precision M0, deltan, eccentricity, rad_A, omega_0, i_0, Omega, omegadot, IDOT, C_uc, C_us, C_rc, C_rs, C_ic, C_is, t_oe, IODE;	//Parameters of subframe 2 and 3
	int TOW = 0;
	int WeekNumber = 0;

	//Private methods
	void decode();
	void process_frame(Subframe *frame_to_process);	//Process a frame: Calculating PB bits and extracting data bits from transmitted bits
	void deocode_frame_parameters(Subframe frame_to_decode);	//Extract parameters from an already processed frame
	bool check_consecutive_preamble();	//Checking for consecutive preambles to detect the start of a frame
	void getFrames();	//Get n_preamble_check frames for the consecutive preambles detected
	void getCurrentFrame();	//Get the last frame that was sampled
	void getSVCoordinates();	//Calculates sattelite coordinates using parameters from subframes 2 and 3
	void check_ephemeris();		//Checks if ephemeris paramters were decoded
	void check_clock_correction();	//Checks if clock correction paramters were decoded

	static mutex navsolJoinMutex;	//Mutex for joining the navsol thread
	static mutex receiver_time_shift_mutex;	//Mutex for setting the receiver starting time
	static mutex fft_mutex;


public:

	//Constructor
	Channel(RealArray PRN, const int sat_ID, const fp_precision fs, const fp_precision chip_rate, const int nsamples, const int loop_order_PLL,
		const fp_precision BPLL, const fp_precision BDLL, fp_precision fd, fp_precision code_shift, int pp_moving_average, int ntrack);
	
	void track(); //mixed_Data should be size of nsamples!

	void addToWorkQ(SampleFrame_Ptr vec_to_add);
	int getworkQsize();
	bool getNeedDataFlag();
	bool check_nav_index(int index);

	mutex Q_mutex;
	mutex vec_mutex;

	fp_precision t_flight = 0;


	Eigen::Matrix<fp_precision, Eigen::Dynamic, 5 > exportChannelData();

	VEC getp();
	VEC gete();
	VEC getl();

	VEC getpp();
	VEC getee();
	VEC getll();

	VEC getfd_vec();
	VEC getphase_error();
	vector<fp_precision> getdebug_vec();
	vector<fp_precision> gettransmission_time_vec();
	vector<fp_precision> getx_kvec();
	vector<fp_precision> gety_kvec();
	vector<fp_precision> getz_kvec();
	vector<bool> getdecodedBits();

	//functions needed by NavSol
	fp_precision getx_kat(int index);
	fp_precision gety_kat(int index);
	fp_precision getz_kat(int index);
	fp_precision getVarianceat(int index);
	fp_precision getrange(int index);

	void sett_flight(fp_precision norm);	//TODO: Change this function so that it takes the x,y,z of receiver and call it whenever channel is added to NavSol object

	static int n_channels_decoded;	//Keeps track of the number of channels that were decoded
	static double receiver_time_shift;	//Receiver starting time common to all of the channels
	bool done_flag;
	static thread NavSolThread;	//NavSol thread common to all the channels
};

double predict_shift(RealArray Input);