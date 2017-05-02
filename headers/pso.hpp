#include <vector>
#include <deque>
#include <algorithm>
#include <ctime>
#include <random>

#include <iostream>

#include <mutex>
#include <condition_variable>

#include "../headers/additional.hpp"

#define HISTORY_NUM 5

std::vector<double> first_derivative(std::vector<double> &data, double step);

class PSO {
public:
	size_t id;
	std::pair<bool, size_t> helper;
	bool initialized = false;

	static std::vector<bool> freq_set;
	bool get_freq_set() { return freq_set[id]; }
	static std::vector<bool> freq_ready;

	size_t numofparticles;
	size_t numofdims;
	size_t numofiterations;
	size_t numofsamples;

	size_t numofsamples_2;

	std::vector< std::vector<double> > V;
	std::vector< std::vector<double> > X;
	std::vector<double> Xmax;
	std::vector<double> Xmin;
	std::vector<double> Vmax;
	std::vector<double> Vmin;
    std::vector< std::vector<double> > pbests;
    std::vector<double> pbestfits;
    std::vector<double> fitnesses;
    double gbestfit;

	std::vector<double> gbest;
	std::vector<double> worsts;
	std::vector<double> meanfits;
	std::vector<double> bests;
	double c1, c2;

    double w;
    double minfit;
    int   minfitidx;
    int   minfithistoryidx;

    std::vector<double> *realdata;
    std::vector<double> *time;
    double fs;	// sampling frequency
    std::vector<double> A;
    std::vector<double> A_gauss;
    std::vector<double> P;

    double max_A;

//    std::vector<double> *found_freqs;
    std::vector<std::vector<double>> *founds;
    std::vector<std::vector<size_t>> to_skip;
    std::vector<std::mutex> *m;
    std::vector<std::condition_variable> *cv;

    double init_param(size_t j);

    void init_max_velocities();

    void initpopulation();

    void addparticle();

    double fitnessfunc_singleparticle(size_t p);

    double fitnessfunc_singleparticle(std::vector<double> &p);

    double calc_response(double amp, double omega, double phase, double bamp, double t);

    void calc_response(std::vector<std::vector<double>> results, size_t numofsamples, double ts, std::vector<double> &response);

    void fitnessfunc();

    void fitnessfunc_thread(size_t start, size_t end);

    void fitnessfunc_multi();

    void calcgbest(bool first = false);

    void update();

    double update_velocity(double w, double X, double V, double V_min, double V_max,
    		double gbest, double pbests, double c1, double c2);

    double update_position(double X, double V, double X_min, double X_max);

    void do_nothing();

public:
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;
	size_t skip_low;
	size_t skip_high;

	PSO(	size_t numofparticles,
			size_t numofdims,
			std::vector<double> Xmin,
			std::vector<double> Xmax,
			std::vector<double> *time,
			std::vector<double> *realdata,
			size_t numofiterations,
			size_t id,
//			std::vector<double> *found_freqs,
			std::vector<std::vector<double>> *founds,
			std::vector<std::mutex> *m,
			std::vector<std::condition_variable> *cv,
			size_t skip_low = 0,
			size_t skip_high = 0,
			std::pair<bool, size_t> helper = std::pair<bool, size_t> (false, -1),
			double c1 = 2, double c2 = 2);

	void run();

	std::vector<double> getgbest();

	double getgbestfit();

	size_t getId() { return id; }

	bool should_skip(size_t freq) {
		double freq_by_idx = freq*fs/(numofsamples_2-1)/2;
		if(freq_by_idx > skip_low && freq_by_idx < skip_high)
			return true;
		else
			return false;
	}
//	bool should_skip(size_t f){
//		for (size_t i = 0; i < to_skip.size(); ++i)
//			if (f >= to_skip[i].first && f <= to_skip[i].second)
//				return true;
//		return false;
//	}

	bool should_skip_2(size_t freq) {
		double freq_by_idx = freq*fs/(numofsamples_2-1)/2;
		if(freq_by_idx > 0.7*skip_high && freq_by_idx < 1.3*skip_high)
			return true;
		else
			return false;
	}

	std::vector<double> getXat(size_t i) { return X.at(i); }

	void setLowSkip(size_t skip) { skip_low = skip; }

	void setHighSkip(size_t skip) { skip_high = skip; }
};
