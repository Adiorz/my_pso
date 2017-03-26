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
	bool initialized = false;

	static bool freq_set;
	static bool freq_ready;

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
//	std::vector< std::deque< std::vector<float> > > pbest_history;
//	std::vector< std::deque<float> > pbest_fits_history;
    std::vector<float> pbestfits;
    std::vector<float> fitnesses;
    float gbestfit;

	std::vector<float> gbest;
//	std::deque< std::vector<float> > gbest_history;
//	std::deque<float> gbest_fits_history;
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
    std::vector<float> P;

    std::vector<float> *found_freqs;
	std::mutex *m;
	std::condition_variable *cv;

    float init_param(size_t j);

    void init_max_velocities();

    void initpopulation();

    void addparticle();

    float fitnessfunc_singleparticle(size_t p);

    float calc_response(float amp, float omega, float phase, float bump, float t);

    void fitnessfunc();

    void fitnessfunc_thread(size_t start, size_t end);

    void fitnessfunc_multi();

    void calcgbest(bool first = false);

    void update();

//    void update_gbest_history(std::vector<float> gbests, float fit) {
//    	if (gbest_history.size() >= HISTORY_NUM) {
//    		gbest_history.pop_front();
//    		gbest_fits_history.pop_front();
//    	}
//    	gbest_history.push_back(gbests);
//		gbest_fits_history.push_back(fit);
//    }
//
//    void update_pbest_history(size_t p, std::vector<float> pbests, float fit) {
//    	if (pbest_history.at(p).size() >= HISTORY_NUM) {
//    		pbest_history.at(p).pop_front();
//    		pbest_fits_history.at(p).pop_front();
//    	}
//    	pbest_history.at(p).push_back(pbests);
//		pbest_fits_history.at(p).push_back(fit);
//    }

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
			std::vector<float> &Xmin,
			std::vector<float> &Xmax,
			std::vector<double> *time,
			std::vector<double> *realdata,
			size_t numofiterations,
			size_t id,
			std::vector<float> *found_freqs,
			std::mutex *m,
			std::condition_variable *cv,
			size_t skip_low = 0,
			size_t skip_high = 0,
			float c1 = 2, float c2 = 2);

	void run();

	std::vector<float> getgbest();

//    std::vector<float> get_gbest_from_history() {
//    	std::deque<float>::iterator minfit_it = std::min_element(std::begin(gbest_fits_history), std::end(gbest_fits_history));
//    	minfithistoryidx = std::distance(std::begin(gbest_fits_history), minfit_it);
//    			return gbest_history.at(minfithistoryidx);
//    }
//
//    std::vector< std::vector<float> > get_pbest_from_history() {
//    	size_t idx = 0;
//    	std::vector< std::vector<float> > pbests;
//    	for (size_t p = 0; p < numofparticles; ++p) {
//    		std::deque<float>::iterator minfit_it = std::min_element(std::begin(pbest_fits_history.at(p)), std::end(pbest_fits_history.at(p)));
//    		idx = std::distance(std::begin(pbest_fits_history.at(p)), minfit_it);
//    		pbests.push_back(pbest_history.at(p).at(idx));
//    	}
//    	return pbests;
//    }

	float getgbestfit();

	size_t getId() { return id; }

	bool should_skip(size_t freq) {
		float freq_by_idx = freq*fs/(numofsamples_2-1)/2;
		if(freq_by_idx > skip_low && freq_by_idx < skip_high)
			return true;
		else
			return false;
	}

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

class PSO_improve: public PSO {
	float calc_response(float amp, float omega, float phase, float bump, float t);
	float fitnessfunc_singleparticle(size_t p);
};
