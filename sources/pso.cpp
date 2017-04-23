#include "../headers/pso.hpp"
#include "../headers/fft.hpp"
#include <algorithm>
#include <random>
#include <chrono>
#include <complex>
#include <map>
#include <iostream>

#include <fstream>

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
	double fitness = fitnessfunc_singleparticle(X[numofparticles-1]);
	fitnesses.push_back(fitness);
	pbestfits.push_back(fitness);
}

//float PSO::fitnessfunc_singleparticle(size_t p) {
//	float fitness = 0.f;
//	while (X[p][3] > X[p][1]) {
//		float high = Xmax[3];
//		float low = Xmin[3];
//		float rand = distribution(generator);
//		float init_position = rand * (high-low) + low;
//		X[p][3] = init_position;
//	}
//	float amp = X[p][0];
//	float omega = X[p][1];
//	float phase = X[p][2];
//	float bump = X[p][3];
//
//	std::vector<double> response(numofsamples, 0);
//	for(size_t j = 0; j < numofsamples; ++j) {
//		float t = time->at(j);
//		response[j] += calc_response(amp, omega, phase, bump, t);
//	}
//
//	std::vector<float> A;
//	std::vector<float> P;
//	fft(response, A, P);
//
//	float freq_found;
//	// get freq found by the first thread
//	if (id > 0 && !helper.first) {
//		std::unique_lock<std::mutex> lk(m->at(id-1));
//		cv->at(id-1).wait(lk, [&]{
//				return freq_ready[id-1];});
//		freq_found = found_freqs->at(id-1);
//		freq_set[id-1] = true;
//		lk.unlock();
//		cv->at(id-1).notify_one();
//	}
//
//	size_t f;
//	size_t f_l;
//	size_t f_r;
//	if (id > 0 && !helper.first) {
//		f = freq_found*numofsamples/fs;
//		std::vector<size_t> idxL;
//		std::vector<size_t> idxR;
//		findMinimas(A_gauss, 0, f, idxL);
//		findMinimas(A_gauss, f, numofsamples-1, idxR);
//		f_l = idxL[idxL.size() - 1];
//		f_r = idxR[0];
//	}
//
//	if (helper.first) {
//		for (size_t j = 0; j < numofsamples_2; ++j) {
//			float residue = 0.;
//			//focus only on the specified range
//			if (j >= skip_low && j <= skip_high) {
//				residue = this->A[j] - A[j];
//				fitness += residue*residue*residue*residue;
//			}
//		}
//	}
//	else
//		for(size_t j = 0; j < numofsamples_2; ++j) {
//			float residue = 0.;
//			residue = this->A[j] - A[j];
//			fitness += residue*residue*residue*residue;
//
//			//////////// additional penalty for the second thread at freq found by the first one
//			if (id > 0 && !helper.first) {
//				if (j >= f_l && j < f_r ) {
//					if (id == 2)
//						std::cout << j << std::endl;
//					float penalty = 0.;
//					if (j < f) {
//						penalty = (j-f_l)/(float)(f-f_l);
//					}
//					else {
//						penalty = (f_r-j)/(float)(f_r-f);
//					}
//					fitness *= abs(A[j])*(1. + penalty);
//				}
//			}
////			if (id > 0 && !helper.first) {
////				for (size_t i = 0; i < to_skip.size(); ++i) {
////					if (j >= to_skip[i][0] && j <= to_skip[i][2]) {
////						float penalty = 0.;
////						if (j < to_skip[i][1]) {
////							penalty = (j-to_skip[i][0])/(float)(to_skip[i][1]-to_skip[i][0]);
////						}
////						else {
////							penalty = (to_skip[i][2]-j)/(float)(to_skip[i][2]-to_skip[i][1]);
////						}
////						fitness *= abs(A[j])*(1. + penalty);
////					}
////				}
////			}
//		}
//
//	return fitness;
//}

float PSO::fitnessfunc_singleparticle(std::vector<float> &p) {
//	std::cout << "reference" << std::endl;
	while (p[3] > p[1]) {
		float high = Xmax[3];
		float low = Xmin[3];
		float rand = distribution(generator);
		float init_position = rand * (high-low) + low;
		p[3] = init_position;
	}

	float amp = p[0];
	float omega = p[1];
	float phase = p[2];
	float bump = p[3];

	std::vector<double> response(numofsamples, 0);
	for(size_t j = 0; j < numofsamples; ++j) {
		float t = time->at(j);
		response[j] += calc_response(amp, omega, phase, bump, t);
//		if (id == 2 && initialized)
//			std::cout << response[j] << std::endl;
	}

	std::vector<float> A;
	std::vector<float> P;
	fft(response, A, P);

//	float freq_found;
//	// get freq found by the first thread
//	if (id > 0 && !helper.first) {
//		std::unique_lock<std::mutex> lk(m->at(id-1));
//		cv->at(id-1).wait(lk, [&]{
//				return freq_ready[id-1];});
//		freq_found = found_freqs->at(id-1);
//		freq_set[id-1] = true;
//		lk.unlock();
//		cv->at(id-1).notify_one();
//	}
//
//	size_t f;
//	size_t f_l;
//	size_t f_r;
//	if (id > 0 && !helper.first) {
//		f = freq_found*numofsamples/fs;
//		std::vector<size_t> idxL;
//		std::vector<size_t> idxR;
//		findMinimas(A_gauss, 0, f, idxL);
//		findMinimas(A_gauss, f, numofsamples-1, idxR);
//		f_l = idxL[idxL.size() - 1];
//		f_r = idxR[0];
//	}
	float fitness = 0.f;
	if (helper.first) {
		for (size_t j = 0; j < numofsamples_2; ++j) {
			//focus only on the specified range
			if (j >= skip_low && j <= skip_high) {
				float residue = this->A[j] - A[j];
				//fitness += residue*residue*residue*residue;
				float temp_fit = residue*residue*residue*residue/max_A;
				fitness += temp_fit;
			}
		}
	}
	else
		for(size_t j = 0; j < numofsamples_2; ++j) {
			float residue = this->A[j] - A[j];
			float temp_fit = residue*residue*residue*residue/max_A;
			if (id > 0 && !helper.first) {
				size_t l_skip;
				size_t m_skip;
				size_t h_skip;
				bool b_skip = false;
				for (size_t i = 0; i < to_skip.size(); ++i) {
					if (j >= to_skip[i][0] && j <= to_skip[i][2]) {
						b_skip = true;
						l_skip = to_skip[i][0];
						m_skip = to_skip[i][1];
						h_skip = to_skip[i][2];
						break;
					}
				}
					if (b_skip) {
						float penalty = 0.;
						if (j < m_skip) {
							penalty = (j-l_skip)/(float)(m_skip-l_skip);
						}
						else {
							penalty = (h_skip-j)/(float)(h_skip-m_skip);
						}
						temp_fit *= abs(A[j])*(1.f + penalty);
				}
			}
		fitness += temp_fit;
	}
	return fitness;
}

float PSO::calc_response(float amp, float omega, float phase, float bump, float t) {
	float ret = amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*t+phase)*
			exp(-bump*t);
	return ret;
}

void PSO::fitnessfunc() {
	for(size_t p = 0; p < numofparticles; p++)
	{
		fitnesses[p] = fitnessfunc_singleparticle(X[p]);
//		fitnesses[p] = fitnessfunc_singleparticle(p);
	}
}

void PSO::fitnessfunc_thread(size_t start, size_t end) {
	std::cout << "multi" << std::endl;
	for(size_t p = start; p < end; p++) {
//		fitnesses[p] = fitnessfunc_singleparticle(X[p]);
//		fitnesses[p] = fitnessfunc_singleparticle(p);
		;
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
	std::vector<double>::iterator minfit_it = std::min_element(std::begin(fitnesses), std::end(fitnesses));
	minfit = *minfit_it;
	minfitidx = std::distance(std::begin(fitnesses), minfit_it);
	if (id > 0 && !helper.first) {
		gbestfit = fitnessfunc_singleparticle(gbest);
	}
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
			std::vector<float> Xmin,
			std::vector<float> Xmax,
			std::vector<double> *time,
			std::vector<double> *realdata,
			size_t numofiterations,
			size_t id,
			std::vector<float> *found_freqs,
			std::vector<std::mutex> *m,
			std::vector<std::condition_variable> *cv,
			size_t skip_low,
			size_t skip_high,
			std::pair<bool, size_t> helper,
			float c1, float c2) {

	this->id = id;
	this->helper = helper;
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
	if (helper.first){
		std::cout << "l: " << skip_low << std::endl;
		std::cout << "h: " << skip_high << std::endl;
		Xmin[1] = skip_low*2*M_PI*fs/numofsamples;
		Xmax[1] = skip_high*2*M_PI*fs/numofsamples;
	}
	Vmax = std::vector<float>(numofdims);
	Vmin = std::vector<float>(numofdims);
	pbests = std::vector< std::vector<float> >(numofparticles);
	pbestfits = std::vector<float>(numofparticles);
	fitnesses = std::vector<double>(numofparticles);

	gbest = std::vector<float>(numofdims);
	worsts = std::vector<float>(numofiterations);
	meanfits = std::vector<float>(numofiterations);
	bests = std::vector<float>(numofiterations);

	this->time = time;
	this->realdata = realdata;

	this->found_freqs = found_freqs;
	to_skip = std::vector<std::vector<size_t>>();
	this->m = m;
	this->cv = cv;

	float ts = time->at(1) - this->time->at(0);
	fs = 1/ts;

	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();

	generator = std::default_random_engine(seed);
	//generator = std::default_random_engine();
	distribution = std::uniform_real_distribution<float>(0.0, 1.0);

	fft(*this->realdata, A, P);
	size_t max_A_idx = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
	max_A = A[max_A_idx];
	A_gauss = gaussian_filter(A, 6);

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

		to_skip.clear();
		if (id > 0 && !helper.first) {
			std::unique_lock<std::mutex> lk(m->at(id-1));
			cv->at(id-1).wait(lk, [&]{
					return freq_ready.at(id-1);});
			//float freq_found = found_freqs->at(id-1);

			for (size_t i = 0; i < id; ++i) {
				size_t f = found_freqs->at(i)*numofsamples/fs;
				std::vector<size_t> idxL;
				std::vector<size_t> idxR;
				findMinimas(A_gauss, 0, f, idxL);
				findMinimas(A_gauss, f, numofsamples-1, idxR);
				to_skip.push_back(std::vector<size_t>{idxL[idxL.size() - 1], f, idxR[0]});
			}
			freq_set[id-1] = true;
			lk.unlock();
			cv->at(id-1).notify_one();
		}

		//fitnessfunc_multi();
		fitnessfunc();
//		if (id == 2 && t == numofiterations - 10) {
//			std::cout << "entered" << std::endl;
//			std::vector<float> temp_particle(this->gbest);
//			std::vector<float> freqs;
//			get_frequencies(freqs, numofsamples, fs);
//			std::cout << "got frequencies" << std::endl;
//
//			std::ofstream myfile;
//			for (size_t f = 0; f < freqs.size(); ++f) {
//				float freq = freqs[f];
//				temp_particle[1] = freq;
//				double penalty = fitnessfunc_singleparticle(temp_particle);
//
//				std::cout << f << "/" << freqs.size() << ", " << freq << ": " << penalty << std::endl;
//				myfile.open ("fitnesses.log", std::ios::out | std::ios::app);
//				myfile  << freq << ";" << penalty << std::endl;
//				myfile.close();
//			}
//		}
		calcgbest();
		std::vector<double>::iterator maxfit_it = std::max_element(std::begin(fitnesses), std::end(fitnesses));
		worsts[t] = *maxfit_it;
		bests[t] = gbestfit;

		//make sure first freq is ready before 2nd thread tries to use it
		if (id < found_freqs->size() && !helper.first) {
		    std::lock_guard<std::mutex> lk(m->at(id));
			found_freqs->at(id) = gbest.at(1)/2/M_PI;
			freq_ready[id] = true;
		    cv->at(id).notify_one();
		}
		//wait for the freq to be read
		if (id < found_freqs->size() && !helper.first) {
			std::unique_lock<std::mutex> lk(m->at(id));
			cv->at(id).wait(lk, [&]{
					return freq_set[id];
			});
		}

		meanfits[t] = std::accumulate(fitnesses.begin(), fitnesses.end(), 0)/numofparticles;
	}
}

std::vector<float> PSO::getgbest() { return gbest; }

float PSO::getgbestfit() { return gbestfit; }
