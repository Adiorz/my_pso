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
protected:
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

	std::vector< std::vector<float> > V;
	std::vector< std::vector<float> > X;
	std::vector<float> Xmax;
	std::vector<float> Xmin;
	std::vector<float> Vmax;
	std::vector<float> Vmin;
    std::vector< std::vector<float> > pbests;
    std::vector<float> pbestfits;
    std::vector<double> fitnesses;
    float gbestfit;

	std::vector<float> gbest;
	std::vector<float> worsts;
	std::vector<float> meanfits;
	std::vector<float> bests;
	float c1, c2;

    float w;
    float minfit;
    int   minfitidx;
    int   minfithistoryidx;

    std::vector<double> *realdata;
    std::vector<double> *time;
    float fs;	// sampling frequency
    std::vector<float> A;
    std::vector<float> A_gauss;
    std::vector<float> P;

    float max_A;

    std::vector<float> *found_freqs;
    std::vector<std::vector<size_t>> to_skip;
    std::vector<std::mutex> *m;
    std::vector<std::condition_variable> *cv;

    float init_param(size_t j);

    void init_max_velocities();

    void initpopulation();

    void addparticle();

    float fitnessfunc_singleparticle(size_t p);

    float fitnessfunc_singleparticle(std::vector<float> &p);

    float calc_response(float amp, float omega, float phase, float bump, float t);

    void fitnessfunc();

    void fitnessfunc_thread(size_t start, size_t end);

    void fitnessfunc_multi();

    void calcgbest(bool first = false);

    void update();

    float update_velocity(float w, float X, float V, float V_min, float V_max,
    		float gbest, float pbests, float c1, float c2);

    float update_position(float X, float V, float X_min, float X_max);

    void do_nothing();

public:
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution;
	size_t skip_low;
	size_t skip_high;

	PSO(	size_t numofparticles,
			size_t numofdims,
			std::vector<float> Xmin,
			std::vector<float> Xmax,
			std::vector<double> *time,
			std::vector<double> *realdata,
			size_t numofiterations,
			size_t id,
			std::vector<float> *found_freqs,
			std::vector<std::mutex> *m,
			std::vector<std::condition_variable> *cv,
			size_t skip_low = 0,
			size_t skip_high = 0,
			std::pair<bool, size_t> helper = std::pair<bool, size_t> (false, -1),
			float c1 = 2, float c2 = 2);

	void run();

	std::vector<float> getgbest();

	float getgbestfit();

	size_t getId() { return id; }

	bool should_skip(size_t freq) {
		float freq_by_idx = freq*fs/(numofsamples_2-1)/2;
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
		float freq_by_idx = freq*fs/(numofsamples_2-1)/2;
		if(freq_by_idx > 0.7*skip_high && freq_by_idx < 1.3*skip_high)
			return true;
		else
			return false;
	}

	std::vector<float> getXat(size_t i) { return X.at(i); }

	void setLowSkip(size_t skip) { skip_low = skip; }

	void setHighSkip(size_t skip) { skip_high = skip; }
};
