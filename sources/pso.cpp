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
    	// rand * (high-low) + low
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
	//response = first_derivative(response, time->at(1) - time->at(0));
	fft(response, A, P);

	for(size_t j = 0; j < numofsamples_2; ++j) {
		float this_a = this->A[j];
		float a = A[j];
		float residue = this->A[j] - A[j];
		fitness += residue*residue*residue*residue;
	}

	return fitness;
}

float PSO::calc_response(float amp, float omega, float phase, float bump, float t) {
	return amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*t+phase)*exp(-bump*t);
}

void PSO::fitnessfunc() {
//	float amp_first = 72.3051;
//	float omega_first = 124*2*M_PI;
//	float phase_first = 0;
//	float bump_first = 0.00332;

	//fitnesses = std::vector<float>(numofparticles, 0.f);

	for(size_t p = 0; p < numofparticles; p++)
	{
		fitnesses[p] = fitnessfunc_singleparticle(p);
	}
}

//    void fitnessfunc_multi() {
//		int* part_size = new int(t_num);
//		for (size_t i = 0; i < t_num - 1; ++i)
//			part_size[i] = ceil(numofparticles / t_num);
//
//		part_size[t_num - 1] = numofparticles - (t_num - 1) * ceil(numofparticles / t_num);
//
//		threads = std::vector<std::thread> (t_num);
//
//		for (size_t i = 0; i < t_num; ++i) {
//			size_t start = i * ceil(numofparticles / t_num);
//			size_t end = i * ceil(numofparticles / t_num) + part_size[i];
//			threads[i] = std::move(std::thread(fitnessfunc_thread, &X, std::ref(fitnesses), realdata,
//					start, end));
//		}
//
//		for (size_t i = 0; i < t_num; ++i)
//			threads[i].join();
//
//		delete part_size;
//    }

void PSO::fitnessfunc_thread(size_t start, size_t end) {
	for(size_t p = start; p < end; p++) {
//		while (X[p][3] > X[p][1]) {
//			float high = Xmax[3];
//			float low = Xmin[3];
//			m.lock();
//			float rand = distribution(generator);
//			m.unlock();
//			float init_position = rand * (high-low) + low;
//			X[p][3] = init_position;
//		}
//		float amp = X[p][0];
//		float omega = X[p][1];
//		float phase = X[p][2];
//		float bump = X[p][3];
//		std::vector<double> response(numofsamples, 0);
//		for(size_t j = 0; j < numofsamples; ++j) {
//			float t = time->at(j);
//			response[j] += calc_response(amp, omega, phase, bump, t);
//		}
//		std::vector<float> A;
//		std::vector<float> P;
//		response = first_derivative(response, time->at(1) - time->at(0));
//		fft(response, A, P);
//		for(size_t j = 0; j < numofsamples_2; ++j) {
//			float residue = this->A[j] - A[j];
//			fitnesses[p] += residue*residue*residue*residue;
		fitnesses[p] = fitnessfunc_singleparticle(p);
//		}
	}
}

void PSO::fitnessfunc_multi() {
	size_t t_num = std::thread::hardware_concurrency();
	//t_num --;
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
			float c1, float c2) {

	this->id = id;
	this->numofparticles = numofparticles;
	this->numofdims = numofdims;
	this->numofiterations = numofiterations;
	this->numofsamples = realdata->size();
	this->numofsamples_2 = size_t(numofsamples/2)+1;
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

	float ts = time->at(1) - this->time->at(0);
	fs = 1/ts;

	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();

	//generator = std::default_random_engine(seed);
	generator = std::default_random_engine();
	distribution = std::uniform_real_distribution<float>(0.0, 1.0);
	this->skip_low  = 0;
	this->skip_high = 0;

	std::cout << "this: " << this << ", id: " << id << ", A: " << &A << std::endl;

	fft(*this->realdata, A, P);

	std::cout << "this: " << this << ", id: " << id << ", A: " << &A << std::endl;

	std::cout << skip_low << std::endl;
	std::cout << skip_high << std::endl;
	std::cout << "-----------------------------------------" << std::endl;
	for (size_t i = 0; i < numofsamples_2; ++i) {
		float freq_by_idx = i*fs/(numofsamples_2-1)/2;
		if(freq_by_idx > skip_low && freq_by_idx < skip_high) {
			A[i] = 0;
			//std::cout << freq_by_idx << " skipped" <<std::endl;
		}
	}
	std::cout << "-----------------------------------------" << std::endl;

	init_max_velocities();
	initpopulation();

	fitnessfunc();
//	fitnessfunc_multi();
	calcgbest(true);
}

//PSO::PSO(const PSO& c) {
//	std::cout << "Copy constructor" << std::endl;
//
//	id = c.id;
//	numofparticles = c.numofparticles;
//	numofdims = c.numofdims;
//	numofiterations = c.numofiterations;
//	numofsamples = c.numofsamples;;
//
//	numofsamples_2 = c.numofsamples_2;
//
//	V = c.V;
//	X = c.X;
//	Xmax = c.Xmax;
//	Xmin = c.Xmin;
//	Vmax = c.Vmax;
//	Vmin = c.Vmin;
//	pbests = c.pbests;
//	pbestfits = c.pbestfits;
//	fitnesses = c.fitnesses;
//	gbestfit = c.gbestfit;
//
//	gbest = c.gbest;
//	worsts = c.worsts;
//	meanfits = c.meanfits;
//	bests = c.bests;
//	c1 = c.c1;
//	c2 = c.c2;
//
//    w = c.w;
//    minfit = c.minfit;
//    minfitidx = c.minfitidx;
//
//    realdata = c.realdata;
//    time = c.time;
//    fs = c.fs;	// sampling frequency
//    A = c.A;
//    P = c.P;
//}

void PSO::run() {
	std::map<size_t, bool> prog;
	for(size_t t = 0; t < numofiterations; t++)
	{
		size_t progress = t*100 / numofiterations;
		if ((progress % 10 == 0) && (!prog[progress])) {
			std::cout << progress << ": " << gbestfit << std::endl;
			prog[progress] = true;
//			addparticle();
//			addparticle();
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

//		fitnessfunc();
		fitnessfunc_multi();
		calcgbest();
		std::vector<float>::iterator maxfit_it = std::max_element(std::begin(fitnesses), std::end(fitnesses));
		worsts[t] = *maxfit_it;
		bests[t] = gbestfit;
		meanfits[t] = std::accumulate(fitnesses.begin(), fitnesses.end(), 0)/numofparticles;
	}
//	std::thread t = std::move(std::thread(&PSO::do_nothing, this));
//	t.join();
}

std::vector<float> PSO::getgbest() { return gbest; }

float PSO::getgbestfit() { return gbestfit; }

void PSO::do_nothing() {
	std::cout << "nothing" << std::endl;
}
