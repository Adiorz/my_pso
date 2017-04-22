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
try {
    size_t numofparticles = atoi(argv[1]);
    size_t numofiterations = atoi(argv[2]);
    std::string file_name(argv[3]);
    std::string deriv(argv[4]);
//    file_name = "data/s10_mod_04.lvm";
//    file_name = "data/bp2_mod_45_01.lvm";
//    file_name = "data/bp2_mod_90_01.lvm";

//	DataStream time("data/s10_mod_04.lvm", 0, 5000);
//	DataStream channel0("data/s10_mod_04.lvm", 1, 5000);
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
	t_num = 3;
	//std::cout << "Number of detected threads: " << t_num << std::endl;

	std::vector< std::vector<float> > init(numofparticles);
	std::vector<float> xmin, xmax;
	init_minmax(xmin, xmax, numofdims, channel0_real);

    std::vector<std::thread> threads(t_num);
    std::vector<std::thread> helper_threads(t_num);

    std::vector<float> *found_freqs = new std::vector<float>(t_num-1);
	std::vector<std::mutex> mutexes(t_num);
	std::vector<std::condition_variable> cvs(t_num);

    std::vector<PSO*> *psos = new std::vector<PSO*>();
	psos->push_back(new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 0, found_freqs, &mutexes, &cvs));
    for (size_t p = 1; p < t_num; ++p) {
    	PSO *pso = new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, p, found_freqs, &mutexes, &cvs);
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
	std::vector<PSO*> *helper_psos = new std::vector<PSO*>();
	size_t f;
	size_t f_l;
	size_t f_r;
	std::vector<size_t> idxL;
	std::vector<size_t> idxR;
	for (size_t p = 0; p < t_num; ++p) {
		f = psos->at(p)->getgbest().at(1)/2/M_PI*numofsamples/fs;
		std::cout << "fixed: " << f << std::endl;
		findMinimas(A_gauss, 0, f, idxL);
		findMinimas(A_gauss, f, numofsamples-1, idxR);
		f_l = idxL[idxL.size() - 1];
		f_r = idxR[0];
		std::cout << "L: " << f_l << std::endl;
		std::cout << "R: " << f_r << std::endl;
		PSO *helper_pso = new PSO(1.5*numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, t_num + p, found_freqs, &mutexes, &cvs, f_l, f_r, std::pair<bool, size_t> (true, p));
		helper_psos->push_back(helper_pso);
		helper_threads[p] = std::thread(&PSO::run, std::ref(*helper_pso));
//		std::cout << "p: " << p << "here" << std::endl;
	}
//    for (size_t p = 0; p < helper_psos->size(); ++p) {
//		std::cout << "p: " << p << "here" << std::endl;
//		helper_threads[p] = std::thread(&PSO::run, std::ref(*helper_psos->at(p)));
//		std::cout << "p: " << p << "here" << std::endl;
//    }
    for (size_t p = 0; p < helper_psos->size(); ++p) {
		helper_threads[p].join();
    }

    std::cout << "first:" << std::endl;
    for (size_t p = 0; p < t_num; ++p) {
    	for (size_t d = 0; d < numofdims; ++d)
    		std::cout << psos->at(p)->getgbest().at(d) << std::endl;
    }
    for (size_t p = 0; p < t_num; ++p) {
    	std::cout << "amp: " << psos->at(p)->getgbest().at(0) << std::endl;
    	std::cout << "freq: " << psos->at(p)->getgbest().at(1)/2/M_PI << std::endl;
    	std::cout << "dumping: " << psos->at(p)->getgbest().at(3)/psos->at(p)->getgbest().at(1) << std::endl;
        std::cout << "fitness: " << psos->at(p)->getgbestfit() <<endl;
    }
    std::cout << "helpers:" << std::endl;
    for (size_t p = 0; p < t_num; ++p) {
    	std::cout << "amp: " << helper_psos->at(p)->getgbest().at(0) << std::endl;
    	std::cout << "freq: " << helper_psos->at(p)->getgbest().at(1)/2/M_PI << std::endl;
    	std::cout << "phase: " << helper_psos->at(p)->getgbest().at(2) << std::endl;
    	std::cout << "dumping: " << helper_psos->at(p)->getgbest().at(3)/helper_psos->at(p)->getgbest().at(1) << std::endl;
        std::cout << "fitness: " << helper_psos->at(p)->getgbestfit() <<endl;
    }


    ofstream myfile;
    myfile.open (file_name + ".log", ios::out | ios::app);
    for (size_t i = 0; i < t_num; ++i) {
    myfile  << helper_psos->at(i)->getgbest().at(0) << ";"
			<< helper_psos->at(i)->getgbest().at(1)/2/M_PI << ";"
			<< helper_psos->at(i)->getgbest().at(2) << ";"
			<< helper_psos->at(i)->getgbest().at(3)/helper_psos->at(i)->getgbest().at(1) << ";";
    }

    myfile << std::endl;
    myfile.close();

    std::vector<std::vector<float>> results;
    std::vector<float> response;
    for (size_t i = 0; i < t_num; ++i)
    	results.push_back(helper_psos->at(i)->getgbest());
    calc_response(results, numofsamples, ts, response);

    myfile.open (file_name + "_response_in_time.log", ios::out);
    for (size_t i = 0; i < numofsamples; ++i) {
    	float real = used_data->at(i);
    	float resp = response[i];
    	myfile  << real << ";" << resp;
    	myfile  << std::endl;
    }
    myfile.close();

#ifdef PLOT
    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
    std::vector<double *> A_approx(t_num+1);
    double mins[t_num];

    for (size_t i = 0; i < t_num; ++i) {
		A_approx[i] = new double[numofsamples_2];
		approximate_amp(
				std::vector<std::vector<float>> {helper_psos->at(i)->getgbest()},
				time_real,
				numofsamples,
				A_approx[i]);
		mins[i] = min(A_approx[i], numofsamples_2);
    }

	A_approx[t_num] = new double[numofsamples_2];
	std::vector<std::vector<float>> factors;
	for (size_t i = 0; i < t_num; ++i)
		factors.push_back(helper_psos->at(i)->getgbest());
    approximate_amp(
    		factors,
			time_real,
			numofsamples,
			A_approx[t_num]);

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

    w->setXArray(freq[0], freq[numofsamples_2-1]);
    //std::cout << "min: " << min_val << std::endl;
    //w->setYArray(min_val, 1.1*20*log(A[max_A_idx]));
    w->setYArray(-100, 100);
    for (size_t i = 0; i < t_num; ++i) {
    	w->addData(x, A_approx[i]);
    }
    w->addData(x, A_approx[t_num]);
    w->display(argc-3, argv+4);
    delete w;
    for (size_t i = 0; i < A_approx.size(); ++i)
    	delete A_approx[i];
#endif

    for (size_t i = 0; i < psos->size(); ++i)
    	delete psos->at(i);
    for (size_t i = 0; i < helper_psos->size(); ++i)
    	delete helper_psos->at(i);
    delete found_freqs;

	}
    catch (std::exception& e)
    {
      cout << e.what() << '\n';
    }
    return 0;
}
