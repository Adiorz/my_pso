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

//	DataStream time("data/input.lvm", 0, 5000);
//	DataStream channel0("data/input.lvm", 1, 5000);
//	DataStream time("data/bp2_mod_45_01.lvm", 0, 5000);
//	DataStream channel0("data/bp2_mod_45_01.lvm", 1, 5000);
	DataStream time("data/bp2_mod_90_01.lvm", 0, 5000);
	DataStream channel0("data/bp2_mod_90_01.lvm", 1, 5000);
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
	//std::vector<double> *used_data = &data_derivative;
	fft(*used_data, A, P);
	//fft(data_derivative, A, P);

	size_t max_A_idx = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
	std::cout << "Original dominating amp: " << A[max_A_idx] << std::endl;
	std::cout << "Original dominating freq: " << freq[max_A_idx] << std::endl;
    std::cout << "Number of particles: " << numofparticles << std::endl;
    std::cout << "Number of iterations: " << numofiterations << std::endl;

	size_t t_num = std::thread::hardware_concurrency();
	t_num = 1;
	std::cout << "Number of detected threads: " << t_num << std::endl;

	std::vector< std::vector<float> > init(numofparticles);
	std::vector<float> xmin, xmax;
	init_minmax(xmin, xmax, numofdims, channel0_real);

    std::vector<std::thread> threads(t_num);
    std::vector<PSO*> *psos = new std::vector<PSO*>();
	psos->push_back(new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 0));
	psos->at(0)->setPSOsVector(psos);
    for (size_t p = 1; p < t_num; ++p) {
    	PSO *pso = new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, p);
    	pso->setPSOsVector(psos);
    	psos->push_back(pso);
    	threads[p] = std::thread(&PSO::run, std::ref(*pso));
    }

    psos->at(0)->run();
    for (size_t p = 1; p < psos->size(); ++p) {
    	threads[p].join();
    }

//	float freq_by_idx = (psos->at(0)->getgbest().at(1)/2/M_PI)*fs/(((size_t)(numofsamples/2)+1)-1)/2;
//	size_t peak = findPeakUtil(*used_data, 0, freq_by_idx);
//    std::cout << "peak: " << peak << std::endl;
    //std::cout << "peak: " << peak << std::endl;

    //size_t low = 0.7*psos->at(0)->getgbest().at(1)/2/M_PI;
    //size_t high = 1.3*psos->at(0)->getgbest().at(1)/2/M_PI;
    size_t high = psos->at(0)->getgbest().at(1)/2/M_PI;

    PSO second = PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 1, 0, high);
    second.run();

//    for (size_t p = 0; p < t_num; ++p) {
//    	for (size_t d = 0; d < numofdims; ++d)
//    		std::cout << psos[p].getgbest().at(d) << std::endl;;
//    }
    for (size_t p = 0; p < t_num; ++p) {
    	std::cout << "freq: " << psos->at(p)->getgbest().at(1)/2/M_PI << std::endl;
        std::cout << "fitness: " << psos->at(p)->getgbestfit() <<endl;
    }
	std::cout << "freq: " << second.getgbest().at(1)/2/M_PI << std::endl;
    std::cout << "fitness: " << second.getgbestfit() <<endl;

#ifdef PLOT
    t_num = 2;
    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
    double **A_approx = new double*[t_num];
    double mins[t_num];
    t_num = 1;
    for (size_t i = 0; i < t_num; ++i) {
    	A_approx[i] = new double[numofsamples_2];
		approximate_amp(
				psos->at(i)->getgbest().at(0),
				psos->at(i)->getgbest().at(1),
				psos->at(i)->getgbest().at(2),
				psos->at(i)->getgbest().at(3),
				time_real,
				numofsamples,
				A_approx[i]);
		mins[i] = min(A_approx[i], numofsamples_2);
	    //std::cout << "min: " << mins[i] << std::endl;
    }
//    std::cout << second.getgbest().at(0) << std::endl;
//    std::cout << second.getgbest().at(1) << std::endl;
//    std::cout << second.getgbest().at(2) << std::endl;
//    std::cout << second.getgbest().at(3) << std::endl;

	A_approx[1] = new double[numofsamples_2];
    approximate_amp(
			second.getgbest().at(0),
			second.getgbest().at(1),
			second.getgbest().at(2),
			second.getgbest().at(3),
			time_real,
			numofsamples,
			A_approx[1]);
	mins[1] = min(A_approx[1], numofsamples_2);

    double min_val = min(mins, t_num);
    double y[numofsamples_2];
    double x[numofsamples_2];
	for(size_t i = 0; i < numofsamples_2; ++i) {
		x[i] = freq[i];
		y[i] = 20*log(A[i]);
	}
    ScatterPlot *w = new ScatterPlot(x, y, (size_t)numofsamples_2);

    for (size_t i = 0; i < t_num; ++i)
    	w->addData(x, A_approx[i]);
	w->addData(x, A_approx[1]);
    w->setXArray(freq[0], freq[numofsamples_2-1]);
    //std::cout << "min: " << min_val << std::endl;
    w->setYArray(min_val, 1.1*20*log(A[max_A_idx]));
    w->display(argc-1, argv+2);
    delete w;
    delete[] A_approx;
#endif
    for (size_t i = 0; i < t_num; ++i)
    	delete psos->at(i);;
    return 0;
}
