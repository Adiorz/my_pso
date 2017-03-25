#include "../headers/pso.hpp"
#include "../headers/fft.hpp"
#include <algorithm>
#include <random>
#include <chrono>
#include <complex>
#include <map>
#include <iostream>

#include <thread>

#include <mutex>

#define rand_01 ((float)rand() / (float)RAND_MAX)

std::mutex m;

std::vector<double> first_derivative(std::vector<double> &data, double step) {
	std::vector<double> derivative(data.size());
	for (size_t i = 0; i < data.size(); ++i) {
		if ((i == 0) || (i == data.size() - 1))
			derivative[i] = data[i];
		else {
			double tmp = data[i-1] + data[i+1];
			tmp /= (2*step);
			derivative[i] = tmp;
		}
	}
	return derivative;
}

float PSO::init_param(size_t j) {
    	return PSO::distribution(generator)*(Xmax.at(j)-Xmin.at(j)) + Xmin.at(j);
    }

void PSO::init_max_velocities() {
	for(size_t i = 0; i < numofdims; i++)
	{
		Vmax.at(i) = 0.2 * (Xmax.at(i) - Xmin.at(i));
		Vmin.at(i) = -Vmax.at(i);
	}
}

void PSO::initpopulation() {
	for(size_t i = 0; i < numofparticles; i++) {
		X.at(i) = std::vector<float>(numofdims);
		V.at(i) = std::vector<float>(numofdims);
		pbests.at(i) = std::vector<float>(numofdims);
		for(size_t j = 0; j < numofdims; j++) {
			X.at(i).at(j) = init_param(j);
			// damping factor cannot be bigger than the frequency (dimensionless damping factor cannot be > 1)
			while ((j == 3) && (X.at(i).at(j) > X.at(i).at(1))) {
				X.at(i).at(j) = init_param(j);
			}
		}
	}
}

void PSO::addparticle() {
	numofparticles++;
	numofsamples_2 = size_t(numofsamples/2)+1;
	std::vector<float> new_x(numofdims);
	std::vector<float> new_v(numofdims, 0);
	for (size_t i = 0; i < numofdims; ++i) {
		new_x[i] = init_param(i);
	}
	X.push_back(new_x);
	V.push_back(new_v);
	pbests.push_back(new_x);
	float fitness = fitnessfunc_singleparticle(numofparticles-1);
	fitnesses.push_back(fitness);
	pbestfits.push_back(fitness);
}

float PSO::fitnessfunc_singleparticle(size_t p) {
	float fitness = 0.f;
	while (X[p][3] > X[p][1]) {
		float high = Xmax[3];
		float low = Xmin[3];
		float rand = distribution(generator);
		float init_position = rand * (high-low) + low;
		X[p][3] = init_position;
	}
	float amp = X[p][0];
	float omega = X[p][1];
	float phase = X[p][2];
	float bump = X[p][3];

	std::vector<double> response(numofsamples, 0);
	for(size_t j = 0; j < numofsamples; ++j) {
		float t = time->at(j);
		response[j] += calc_response(amp, omega, phase, bump, t);
	}

	std::vector<float> A;
	std::vector<float> P;
	fft(response, A, P);

	float freq_found_by_first;
	//get freq found by the first thread
//	if (id == 1) {
//		std::unique_lock<std::mutex> lk(*m);
//		cv->wait(lk, []{return freq_ready;});
//		freq_found_by_first = found_freqs->at(0);
//		freq_set = true;
//		lk.unlock();
//		cv->notify_one();
//	}


	size_t f;
	size_t f_07;
	size_t f_13;
	if (id == 1) {
		f = skip_high*fs/(numofsamples_2-1)/2;
		f_07 = 0.7*skip_high*fs/(numofsamples_2-1)/2;
		f_13 = 1.3*skip_high*fs/(numofsamples_2-1)/2;
		findMinima(this->A, f, f_07, f_13);
	}

	for(size_t j = 0; j < numofsamples_2; ++j) {
		float residue = 0.;
		residue = this->A[j] - A[j];
		fitness += residue*residue*residue*residue;

		//////////// additional penalty for the second thread at freq found by the first one
		if (id == 1) {
			if (should_skip_2(j)) {
				float penalty = 0.;
				if (j < f) {
					penalty = (j-f_07)/(float)(f-f_07);
				}
				else {
					penalty = (f_13-j)/(float)(f_13-f);
				}
				fitness *= abs(A[j])*(1. + penalty);
			}
		}
	}

	return fitness;
}

float PSO::calc_response(float amp, float omega, float phase, float bump, float t) {
	return amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*t+phase)*exp(-bump*t);
}

void PSO::fitnessfunc() {
	for(size_t p = 0; p < numofparticles; p++)
	{
		fitnesses[p] = fitnessfunc_singleparticle(p);
	}
}

void PSO::fitnessfunc_thread(size_t start, size_t end) {
	for(size_t p = start; p < end; p++) {
		fitnesses[p] = fitnessfunc_singleparticle(p);
	}
}

void PSO::fitnessfunc_multi() {
	size_t t_num = std::thread::hardware_concurrency();
	t_num --;
	int* part_size = new int(t_num);
	for (size_t i = 0; i < t_num - 1; ++i)
		part_size[i] = ceil(numofparticles / t_num);

	part_size[t_num - 1] = numofparticles - (t_num - 1) * ceil(numofparticles / t_num);

	std::vector<std::thread> threads(t_num);

	for (size_t i = 0; i < t_num; ++i) {
		size_t start = i * ceil(numofparticles / t_num);
		size_t end = i * ceil(numofparticles / t_num) + part_size[i];
		threads[i] = std::move(std::thread(&PSO::fitnessfunc_thread, this, start, end));
	}

	for (size_t i = 0; i < t_num; ++i)
		threads[i].join();

	delete part_size;
}

void PSO::calcgbest(bool first) {
	std::vector<float>::iterator minfit_it = std::min_element(std::begin(fitnesses), std::end(fitnesses));
	minfit = *minfit_it;
	minfitidx = std::distance(std::begin(fitnesses), minfit_it);
	if(first || minfit < gbestfit)
	{
		gbestfit = minfit;
		// change for fast vector copying
		for (size_t i = 0; i < numofdims; ++i)
			gbest.at(i) = X.at(minfitidx).at(i);
	}
}

void PSO::update() {
	for(size_t i = 0; i < numofparticles; i++) {
		for(size_t j = 0; j < numofdims; j++) {
			// update velocity
			V.at(i).at(j) = update_velocity(w, X.at(i).at(j), V.at(i).at(j), Vmin[j], Vmax[j], gbest.at(j), pbests.at(i).at(j), c1, c2);
			// update position
			X.at(i).at(j) = update_position(X.at(i).at(j), V.at(i).at(j), Xmin[j], Xmax[j]);
		}
	}
}

float PSO::update_velocity(float w, float X, float V, float V_min, float V_max,
		float gbest, float pbests, float c1, float c2) {
	return std::min(std::max((w * V + rand_01 * c1 * (pbests - X)
			+ rand_01 * c2 * (gbest - X)), V_min), V_max);
}

float PSO::update_position(float X, float V, float X_min, float X_max) {
	return std::min(std::max((X + V), X_min), X_max);
}

PSO::PSO(	size_t numofparticles,
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
			size_t skip_low,
			size_t skip_high,
			float c1, float c2) {

	this->id = id;
	this->numofparticles = numofparticles;
	this->numofdims = numofdims;
	this->numofiterations = numofiterations;
	this->numofsamples = realdata->size();
	this->numofsamples_2 = size_t(numofsamples/2)+1;
	this->skip_low = skip_low;
	this->skip_high = skip_high;
	this->c1 = c1;
	this->c2 = c2;

	V = std::vector< std::vector<float> >(numofparticles);
	X = std::vector< std::vector<float> >(numofparticles);
	this->Xmax = Xmax;
	this->Xmin = Xmin;
	Vmax = std::vector<float>(numofdims);
	Vmin = std::vector<float>(numofdims);
	pbests = std::vector< std::vector<float> >(numofparticles);
	pbestfits = std::vector<float>(numofparticles);
	fitnesses = std::vector<float>(numofparticles);

	gbest = std::vector<float>(numofdims);
	worsts = std::vector<float>(numofiterations);
	meanfits = std::vector<float>(numofiterations);
	bests = std::vector<float>(numofiterations);

	this->time = time;
	this->realdata = realdata;

	this->found_freqs = found_freqs;
	this->m = m;
	this->cv = cv;

	float ts = time->at(1) - this->time->at(0);
	fs = 1/ts;

	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();

	//generator = std::default_random_engine(seed);
	generator = std::default_random_engine();
	distribution = std::uniform_real_distribution<float>(0.0, 1.0);

	fft(*this->realdata, A, P);

	init_max_velocities();
	initpopulation();
}

void PSO::run() {
	if (!initialized) {
		fitnessfunc();
		calcgbest(true);
		initialized = true;
	}
	std::map<size_t, bool> prog;
	for(size_t t = 0; t < numofiterations; t++)
	{
		size_t progress = t*100 / numofiterations;
		if ((progress % 10 == 0) && (!prog[progress])) {
			std::cout << "id: " << id << ": " << progress << ": " << gbestfit << std::endl;
			prog[progress] = true;
			addparticle();
			addparticle();
		}
		if (numofiterations)
			w = 0.9 - 0.7 * t / numofiterations;

		for(size_t i = 0; i < numofparticles; i++)
		{
			if(fitnesses[i] < pbestfits[i])
			{
				pbestfits[i] = fitnesses[i];
				for(size_t j = 0; j < numofdims; j++)
					pbests[i][j] = X[i][j];
			}
		}
		update();

		fitnessfunc_multi();
		calcgbest();
		std::vector<float>::iterator maxfit_it = std::max_element(std::begin(fitnesses), std::end(fitnesses));
		worsts[t] = *maxfit_it;
		bests[t] = gbestfit;

		//make sure first freq is ready before 2nd thread tries to use it
//		if (id == 0) {
//		    std::lock_guard<std::mutex> lk(*m);
//			found_freqs->at(0) = gbest.at(0);
//			freq_ready = true;
//		    cv->notify_one();
//		}
//		//wait for the freq to be read
//		if (id == 0) {
//			std::unique_lock<std::mutex> lk(*m);
//			cv->wait(lk, []{return freq_set;});
//		}

		meanfits[t] = std::accumulate(fitnesses.begin(), fitnesses.end(), 0)/numofparticles;
	}
}

std::vector<float> PSO::getgbest() { return gbest; }

float PSO::getgbestfit() { return gbestfit; }

void PSO::do_nothing() {
	std::cout << "nothing" << std::endl;
}

bool PSO::freq_set = false;
bool PSO::freq_ready = false;
