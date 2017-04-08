//============================================================================
// Name        : my_pso.cpp
// Author      : Adrian Orzechowski
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;

#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>

#include <fstream>
#include <sstream>

#include <thread>

#include <mutex>
#include <condition_variable>

#define PLOT
#ifdef PLOT
#include "../headers/scatter.h"
#endif

#include "../headers/pso.hpp"
#include "../headers/fft.hpp"
#include "../headers/additional.hpp"

const size_t numofdims = 4;
const float samplingrate = 50; // samples per second

using namespace std;


int main(int argc, const char ** argv)
{

    size_t numofparticles = atoi(argv[1]);
    size_t numofiterations = atoi(argv[2]);
    std::string file_name(argv[3]);
    std::string deriv(argv[4]);

//	DataStream time("data/input.lvm", 0, 5000);
//	DataStream channel0("data/input.lvm", 1, 5000);
//	DataStream time("data/bp2_mod_45_01.lvm", 0, 5000);
//	DataStream channel0("data/bp2_mod_45_01.lvm", 1, 5000);
//	DataStream time("data/bp2_mod_90_01.lvm", 0, 5000);
//	DataStream channel0("data/bp2_mod_90_01.lvm", 1, 5000);
	DataStream time(file_name, 0, 5000);
	DataStream channel0(file_name, 1, 5000);
	size_t numofsamples = time.size();
	std::vector<double> channel0_real	(channel0.get());
	std::vector<double> time_real		(time.get());
	float ts = time_real[1] - time_real[0];
	float fs = 1/ts;
	std::vector<double> data_derivative(first_derivative(channel0_real, ts));
	std::vector<float> freq;
	get_frequencies(freq, numofsamples, fs);
	std::vector<float> A;
	std::vector<float> P;
	std::vector<double> *used_data = &channel0_real;
//	std::vector<double> *used_data = &data_derivative;
	fft(*used_data, A, P);

	if (deriv == "deriv") {
		std::cout << "use signal derivative as input" << std::endl;
		used_data = &data_derivative;
	}
	size_t max_A_idx = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
	std::cout << "Original dominating amp: " << A[max_A_idx] << std::endl;
	std::cout << "Original dominating freq: " << freq[max_A_idx] << std::endl;
    std::cout << "Number of particles: " << numofparticles << std::endl;
    std::cout << "Number of iterations: " << numofiterations << std::endl;

	size_t t_num = std::thread::hardware_concurrency();
	t_num = 2;
	//std::cout << "Number of detected threads: " << t_num << std::endl;

	std::vector< std::vector<float> > init(numofparticles);
	std::vector<float> xmin, xmax;
	init_minmax(xmin, xmax, numofdims, channel0_real);

    std::vector<std::thread> threads(t_num);

    std::vector<float> *found_freqs = new std::vector<float>(2);
	std::mutex m;
	std::condition_variable cv;

    std::vector<PSO*> *psos = new std::vector<PSO*>();
	psos->push_back(new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 0, found_freqs, &m, &cv));
    for (size_t p = 1; p < t_num; ++p) {
    	PSO *pso = new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, p, found_freqs, &m, &cv);
    	psos->push_back(pso);
    	threads[p] = std::thread(&PSO::run, std::ref(*pso));
    }

    psos->at(0)->run();
    for (size_t p = 1; p < psos->size(); ++p) {
    	threads[p].join();
    }

    used_data = &channel0_real;
	fft(*used_data, A, P);
	std::vector<float> A_gauss(gaussian_filter(A, 6));
	size_t f = psos->at(1)->getgbest().at(1)/2/M_PI*numofsamples/fs;
	std::cout << "fixed: " << f << std::endl;
	std::vector<size_t> idxL;
	std::vector<size_t> idxR;
	findMinimas(A_gauss, 0, f, idxL);
	findMinimas(A_gauss, 0, f, idxR);
	transform(idxR.begin(), idxR.end(), idxR.begin(),
	          bind2nd(std::plus<size_t>(), f));
	size_t f_l;
	size_t f_r;
	f_l = idxL[idxL.size() - 1];
	f_r = idxR[0];
	std::cout << f_l << std::endl;
	std::cout << f_r << std::endl;
	PSO third = PSO(1.5*numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 2, found_freqs, &m, &cv, f_l, f_r);
	third.run();
	f = psos->at(0)->getgbest().at(1)/2/M_PI*numofsamples/fs;
	std::cout << "fixed: " << f << std::endl;
	findMinimas(A_gauss, 0, f, idxL);
	findMinimas(A_gauss, 0, f, idxR);
	transform(idxR.begin(), idxR.end(), idxR.begin(),
	          bind2nd(std::plus<size_t>(), f));
	f_l = idxL[idxL.size() - 1];
	f_r = idxR[0];
	std::cout << f_l << std::endl;
	std::cout << f_r << std::endl;
	PSO forth = PSO(1.5*numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 3, found_freqs, &m, &cv, f_l, f_r);
	forth.run();

    std::cout << "first:" << std::endl;
    for (size_t p = 0; p < t_num; ++p) {
    	for (size_t d = 0; d < numofdims; ++d)
    		std::cout << psos->at(p)->getgbest().at(d) << std::endl;
    }
    for (size_t p = 0; p < t_num; ++p) {
    	std::cout << "freq: " << psos->at(p)->getgbest().at(1)/2/M_PI << std::endl;
    	std::cout << "dumping: " << psos->at(p)->getgbest().at(3)/psos->at(p)->getgbest().at(1) << std::endl;
        std::cout << "fitness: " << psos->at(p)->getgbestfit() <<endl;
    }
	std::cout << "amp: " << third.getgbest().at(0) << std::endl;
	std::cout << "freq: " << third.getgbest().at(1)/2/M_PI << std::endl;
	std::cout << "dumping: " << third.getgbest().at(3)/third.getgbest().at(1) << std::endl;
    std::cout << "fitness: " << third.getgbestfit() <<endl;
	std::cout << "amp: " << forth.getgbest().at(0) << std::endl;
	std::cout << "freq: " << forth.getgbest().at(1)/2/M_PI << std::endl;
	std::cout << "dumping: " << forth.getgbest().at(3)/forth.getgbest().at(1) << std::endl;
    std::cout << "fitness: " << forth.getgbestfit() <<endl;

#ifdef PLOT
    t_num = 4;
    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
    double **A_approx = new double*[t_num+1];
    double mins[t_num];
    t_num = 2;
    for (size_t i = 0; i < t_num; ++i) {
    	A_approx[i] = new double[numofsamples_2];
		approximate_amp(
				std::vector<std::vector<float>> {psos->at(i)->getgbest()},
				time_real,
				numofsamples,
				A_approx[i]);
		mins[i] = min(A_approx[i], numofsamples_2);
    }
    t_num = 2;
	A_approx[2] = new double[numofsamples_2];
	approximate_amp(
			std::vector<std::vector<float>> {third.getgbest()},
			time_real,
			numofsamples,
			A_approx[2]);
	mins[2] = min(A_approx[2], numofsamples_2);
	A_approx[3] = new double[numofsamples_2];
	approximate_amp(
			std::vector<std::vector<float>> {forth.getgbest()},
			time_real,
			numofsamples,
			A_approx[3]);
	mins[3] = min(A_approx[3], numofsamples_2);

	A_approx[4] = new double[numofsamples_2];
	std::vector<std::vector<float>> factors;
	factors.push_back(third.getgbest());
	factors.push_back(forth.getgbest());
    approximate_amp(
    		factors,
			time_real,
			numofsamples,
			A_approx[4]);
	std:: cout << A_approx[4][123] << std::endl;

    double min_val = min(mins, t_num);
    double g[numofsamples_2];
    double y[numofsamples_2];
    double x[numofsamples_2];
	for(size_t i = 0; i < numofsamples_2; ++i) {
		g[i] = 20*log(A_gauss[i]);
		x[i] = freq[i];
		y[i] = 20*log(A[i]);
	}
    ScatterPlot *w = new ScatterPlot(x, y, (size_t)numofsamples_2);
	//w->addData(x, g);

	//w->addData(x, A_approx[1]);
    w->setXArray(freq[0], freq[numofsamples_2-1]);
    //std::cout << "min: " << min_val << std::endl;
    //w->setYArray(min_val, 1.1*20*log(A[max_A_idx]));
    w->setYArray(-100, 100);
//    for (size_t i = 0; i < t_num; ++i)
//    	w->addData(x, A_approx[i]);
    w->addData(x, A_approx[2]);
    w->addData(x, A_approx[3]);
    w->addData(x, A_approx[4]);
    w->display(argc-3, argv+4);
    delete w;
    delete found_freqs;
    delete[] A_approx;
#endif

    for (size_t i = 0; i < t_num; ++i)
    	delete psos->at(i);
    ofstream myfile;
    myfile.open (file_name + ".log", ios::in | ios::app);
    myfile  << third.getgbest().at(0) << ";"
			<< third.getgbest().at(1)/2/M_PI << ";"
			<< third.getgbest().at(3)/third.getgbest().at(1) << ";"
			<< forth.getgbest().at(0) << ";"
			<< forth.getgbest().at(1)/2/M_PI << ";"
			<< forth.getgbest().at(3)/forth.getgbest().at(1) << ";" << std::endl;

    myfile.close();
    return 0;
}
