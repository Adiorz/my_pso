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
	DataStream time("data/bp2_mod_45_01.lvm", 0, 5000);
	DataStream channel0("data/bp2_mod_45_01.lvm", 1, 5000);
	size_t numofsamples = time.size();
	std::vector<double> channel0_real	(channel0.get());
	std::vector<double> time_real		(time.get());
	float ts = time_real[1] - time_real[0];
	float fs = 1/ts;
	//std::vector<double> data_derivative(first_derivative(channel0_real, ts));
	std::vector<float> freq;
	get_frequencies(freq, numofsamples, fs);
	std::vector<float> A;
	std::vector<float> P;
	fft(channel0_real, A, P);
	//fft(data_derivative, A, P);

	size_t max_A_idx = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
	std::cout << "Original dominating amp: " << A[max_A_idx] << std::endl;
	std::cout << "Original dominating freq: " << freq[max_A_idx] << std::endl;
    std::cout << "Number of particles: " << numofparticles << std::endl;
    std::cout << "Number of iterations: " << numofiterations << std::endl;

	size_t t_num = std::thread::hardware_concurrency();
	t_num = 3;
	std::cout << "Number of detected threads: " << t_num << std::endl;

	std::vector< std::vector<float> > init(numofparticles);
	std::vector<float> xmin, xmax;
	init_minmax(xmin, xmax, numofdims, channel0_real);

    std::vector<std::thread> threads(t_num);
    std::vector<PSO*> psos;
	psos.push_back(new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &channel0_real, numofiterations, 0));
    for (size_t p = 1; p < t_num; ++p) {
    	PSO *pso = new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &channel0_real, numofiterations, p);
    	psos.push_back(pso);
    	//threads[p] = std::move(std::thread(&PSO::run, std::ref(psos[p])));
    	//threads[p] = std::move(std::thread(&PSO::run, std::ref(psos[p])));
    	std::cout << "p: " << p << ", addr: " << pso << std::endl;
    	std::cout << "p: " << p << ", addr: " << psos[p] << std::endl;
    	threads[p] = std::thread(&PSO::run, std::ref(*pso));
    }

    psos[0]->run();
    for (size_t p = 1; p < psos.size(); ++p) {
    	std::cout << "p: " << p << std::endl;
    	threads[p].join();
    }


//    for (size_t p = 0; p < t_num; ++p) {
//    	for (size_t d = 0; d < numofdims; ++d)
//    		std::cout << psos[p].getgbest().at(d) << std::endl;;
//    }
    for (size_t p = 0; p < t_num; ++p) {
    	std::cout << "freq: " << psos[p]->getgbest().at(1)/2/M_PI << std::endl;
        std::cout << "fitness: " << psos[p]->getgbestfit() <<endl;
    }

#ifdef PLOT

    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
    double **A_approx = new double*[t_num];
    double mins[t_num];
    for (size_t i = 0; i < t_num; ++i) {
    	A_approx[i] = new double[numofsamples_2];
		approximate_amp(
				psos[0]->getgbest().at(0),
				psos[0]->getgbest().at(1),
				psos[0]->getgbest().at(2),
				psos[0]->getgbest().at(3),
				time_real,
				numofsamples,
				A_approx[i]);
		mins[i] = min(A_approx[i], numofsamples_2);
	    std::cout << "min: " << mins[i] << std::endl;
    }

    double min_val = min(mins, t_num);
    double y[numofsamples_2];
    double x[numofsamples_2];
	for(size_t i = 0; i < numofsamples_2; ++i) {
		x[i] = freq[i];
		y[i] = 20*log(A[i]);
	}
    ScatterPlot *w = new ScatterPlot(x, y, (size_t)numofsamples_2);

    w->addData(x, A_approx[0]);
    w->addData(x, A_approx[1]);
    w->setXArray(freq[0], freq[numofsamples_2-1]);
    std::cout << "min: " << min_val << std::endl;
    w->setYArray(min_val, 1.1*20*log(A[max_A_idx]));
    w->display(argc-1, argv+2);
    delete w;
    delete[] A_approx;
#endif
    return 0;
}
