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
#include <iomanip>

#define rand_01 ((double)rand() / (double)RAND_MAX)

std::mutex m;

double PSO::init_param(size_t j) {
		if (j == 1)
	    	return PSO::distribution(generator)*PSO::distribution(generator)*(Xmax.at(j)-Xmin.at(j)) + Xmin.at(j);
		if (j == 3)
	    	return PSO::distribution(generator)*PSO::distribution(generator)*PSO::distribution(generator)*(Xmax.at(j)-Xmin.at(j)) + Xmin.at(j);
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
		X.at(i) = std::vector<double>(numofdims);
		V.at(i) = std::vector<double>(numofdims);
		pbests.at(i) = std::vector<double>(numofdims);
		for(size_t j = 0; j < numofdims; j++) {
			X.at(i).at(j) = init_param(j);
			// damping factor cannot be bigger than the frequency (dimensionless damping factor cannot be > 1)
			while ((j == 3) && (X.at(i).at(j) > X.at(i).at(1))) {
				X.at(i).at(j) = init_param(j);
			}
		}
		if (helper.first && i == 0) {
			std::cout << founds->size() << std::endl;
			std::cout << founds->at(id)[0] << std::endl;
			std::cout << founds->at(id)[1] << "(" << founds->at(id)[1]/2/M_PI << ")" << std::endl;
			std::cout << founds->at(id)[2] << std::endl;
			std::cout << founds->at(id)[3] << "(" << founds->at(id)[3]/founds->at(id)[1] << ")" << std::endl;
			X.at(i) = std::vector<double> (founds->at(id));
		}
	}
}

void PSO::addparticle() {
	numofparticles++;
	numofsamples_2 = size_t(numofsamples/2)+1;
	std::vector<double> new_x(numofdims);
	std::vector<double> new_v(numofdims, 0);
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

double PSO::fitnessfunc_singleparticle(std::vector<double> &p) {
	while (p[3] > p[1]) {
		double high = Xmax[3];
		double low = Xmin[3];
		double rand = distribution(generator);
		double init_position = rand * (high-low) + low;
		p[3] = init_position;
	}

	double amp = p[0];
	double omega = p[1];
	double phase = p[2];
	double damp = p[3];

	std::vector<double> response(numofsamples, 0);
	if (helper.first) {
		std::vector<std::vector<double>> parameters;
		for (size_t i = 0; i < id; ++i)
			parameters.push_back(local_copy_of_founds.at(i));
		parameters.push_back(p);
		calc_response(parameters, numofsamples, 1/fs, response);
	}
	else
		for(size_t j = 0; j < numofsamples; ++j) {
			double t = time->at(j);
			response[j] += calc_response(amp, omega, phase, damp, t);
		}

	std::vector<double> A;
	std::vector<double> P;
	fft(response, A, P);
//    std::vector<double>::iterator max_abs_it;
//    max_abs_it = std::max_element(realdata->begin(), realdata->end(), abs_compare);
//	size_t max_abs_idx = std::distance(realdata->begin(), max_abs_it);
//	double max_abs = abs(realdata->at(max_abs_idx));

	double fitness = 0.f;
	if (helper.first) {
		for (size_t j = 0; j < numofsamples_2; ++j) {
			//focus only on the specified range
			if (j >= skip_low && j <= skip_high) {
				double residue = this->A[j] - A[j];
				//fitness += residue*residue*residue*residue;
				double temp_fit = residue*residue*residue*residue/max_A;
				fitness += temp_fit;
			}
		}
	}
	else
		for(size_t j = 0; j < numofsamples_2; ++j) {
			double residue = this->A[j] - A[j];
			double temp_fit = residue*residue*residue*residue/max_A;
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
						double penalty = 0.;
						if (j < m_skip) {
							penalty = (j-l_skip)/(double)(m_skip-l_skip);
						}
						else {
							penalty = (h_skip-j)/(double)(h_skip-m_skip);
						}
//						temp_fit *= abs(local_copy_of_founds[0][0]/max_abs*A[j])*(1.f + penalty);
						temp_fit *= abs(local_copy_of_founds[0][0]/max_A*A[j])*(1.f + penalty);
//						temp_fit *= abs(A[j])*(1.f + penalty);
				}
			}
		fitness += temp_fit;
	}
	return fitness;
}

double PSO::calc_response(double amp, double omega, double phase, double damp, double t) {
//	double ret = amp*sin(omega*sqrt(1-(damp/omega)*(damp/omega))*t+phase);
	double ret = amp*sin(omega*sqrt(1-(damp/omega)*(damp/omega))*t);
	ret *= exp(-damp*t);
	return ret;
}

void PSO::calc_response(std::vector<std::vector<double>> results, size_t numofsamples, double ts, std::vector<double> &response) {
	response = std::vector<double> (numofsamples, 0.);
	for (size_t i = 0; i < results.size(); ++i) {
		double amp = results[i][0];
		double omega = results[i][1];
		double phase = results[i][2];
		double damp = results[i][3];
		for (size_t j = 0; j < numofsamples; ++j) {
			double t = j*ts;
			response[j] += amp*sin(omega*sqrt(1-(damp/omega)*(damp/omega))*t)*
					exp(-damp*t);
		}
	}
}

void PSO::fitnessfunc() {
	for(size_t p = 0; p < numofparticles; p++)
	{
		fitnesses[p] = fitnessfunc_singleparticle(X[p]);
	}
}

void PSO::fitnessfunc_thread(size_t start, size_t end) {
	std::cout << "multi" << std::endl;
	for(size_t p = start; p < end; p++) {
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

double PSO::update_velocity(double w, double X, double V, double V_min, double V_max,
		double gbest, double pbests, double c1, double c2) {
	return std::min(std::max((w * V + rand_01 * c1 * (pbests - X)
			+ rand_01 * c2 * (gbest - X)), V_min), V_max);
}

double PSO::update_position(double X, double V, double X_min, double X_max) {
	return std::min(std::max((X + V), X_min), X_max);
}

PSO::PSO(	size_t numofparticles,
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
			size_t skip_low,
			size_t skip_high,
			std::pair<bool, size_t> helper,
			double c1, double c2) {

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

	V = std::vector< std::vector<double> >(numofparticles);
	X = std::vector< std::vector<double> >(numofparticles);
	Vmax = std::vector<double>(numofdims);
	Vmin = std::vector<double>(numofdims);
	pbests = std::vector< std::vector<double> >(numofparticles);
	pbestfits = std::vector<double>(numofparticles);
	fitnesses = std::vector<double>(numofparticles);

	gbest = std::vector<double>(numofdims);
	worsts = std::vector<double>(numofiterations);
	meanfits = std::vector<double>(numofiterations);
	bests = std::vector<double>(numofiterations);

	this->time = time;
	this->realdata = realdata;

//	this->found_freqs = found_freqs;
	this->founds = founds;
	this->local_copy_of_founds = *founds;
	to_skip = std::vector<std::vector<size_t>>();
	this->m = m;
	this->cv = cv;

	double ts = time->at(1) - this->time->at(0);
	fs = 1/ts;

	this->Xmax = Xmax;
	this->Xmin = Xmin;
	if (helper.first){
		std::cout << "l: " << skip_low << std::endl;
		std::cout << "h: " << skip_high << std::endl;
		float new_min = skip_low*2*M_PI*fs/numofsamples;
		float new_max = skip_high*2*M_PI*fs/numofsamples;
		if (new_min > this->Xmin[1])
			this->Xmin[1] = new_min;
		if (new_max < this->Xmax[1])
			this->Xmax[1] = new_max;
	}

	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();

	generator = std::default_random_engine(seed);
	//generator = std::default_random_engine();
	distribution = std::uniform_real_distribution<double>(0.0, 1.0);

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
		if (progress % 2 == 0) {
			addparticle();
			addparticle();
		}
		if ((progress % 10 == 0) && (!prog[progress])) {
			std::cout << "id: " << id << ": " << progress << ": " << std::fixed << std::setprecision(17) << gbestfit << std::setprecision(5) << std::endl;//<< ", the worst: " << fitnesses[the_worst_fit_idx] << std::endl;
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

		to_skip.clear();
		if (id > 0) {
			std::unique_lock<std::mutex> lk(m->at(id-1));
			cv->at(id-1).wait(lk, [&]{
					return freq_ready.at(id-1);});
			this->local_copy_of_founds = *founds;
			freq_set[id-1] = true;
			lk.unlock();
			cv->at(id-1).notify_one();

			for (size_t i = 0; i < id; ++i) {
				size_t f = local_copy_of_founds.at(i)[1]/2/M_PI*numofsamples/fs;
				std::vector<size_t> idxL;
				std::vector<size_t> idxR;
				findMinimas(A_gauss, 0, f, idxL);
				findMinimas(A_gauss, f, numofsamples-1, idxR);
				if (idxR.size() == 0)
					idxR.push_back(numofsamples-1);
				to_skip.push_back(std::vector<size_t>{idxL[idxL.size() - 1], f, idxR[0]});
			}
		}

		fitnessfunc();
		calcgbest();
		std::vector<double>::iterator maxfit_it = std::max_element(std::begin(fitnesses), std::end(fitnesses));
		worsts[t] = *maxfit_it;
		bests[t] = gbestfit;

		//make sure first freq is ready before 2nd thread tries to use it
		if (id < founds->size()) {
		    std::lock_guard<std::mutex> lk(m->at(id));
		    founds->at(id) = std::vector<double>(gbest);
		    // TODO: do not divide
			//founds->at(id)[1] /= 2*M_PI;
			freq_ready[id] = true;
		    cv->at(id).notify_one();
		}
		//wait for the freq to be read
		if (id < founds->size() - 1) {
			std::unique_lock<std::mutex> lk(m->at(id));
			cv->at(id).wait(lk, [&]{
					return freq_set[id];
			});
		}

		meanfits[t] = std::accumulate(fitnesses.begin(), fitnesses.end(), 0)/numofparticles;
	}
}

std::vector<double> PSO::getgbest() { return gbest; }

double PSO::getgbestfit() { return gbestfit; }
