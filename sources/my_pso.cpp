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
const double samplingrate = 50; // samples per second

using namespace std;

std::vector<bool> PSO::freq_set = std::vector<bool> ();
std::vector<bool> PSO::freq_ready = std::vector<bool> ();

int main(int argc, const char ** argv)
{
try {
    size_t numofparticles = atoi(argv[1]);
    numofparticles = 20;
    size_t numofiterations = atoi(argv[2]);
    std::string file_name(argv[4]);
    std::cout << "Filename: " << file_name << std::endl;
    std::string window_name(argv[5]);
    size_t modes_num = atoi(argv[3]);
//    file_name = "data/s10_mod_04.lvm";
//    file_name = "data/bp2_mod_45_01.lvm";
//    file_name = "data/bp2_mod_90_01.lvm";

	DataStream time(file_name, 0, 5000, 0);
//	DataStream channel0(file_name, 1, 5000, 166);
	DataStream channel0(file_name, 1, 5000, 0);
	size_t numofsamples = channel0.size();
	std::vector<double> channel0_real	(channel0.get());
	std::vector<double> time_real		(time.get());
	double ts = time_real[1] - time_real[0];
	double fs = 1/ts;
	std::vector<double> freq;
	get_frequencies(freq, numofsamples, fs);
	std::vector<double> A;
	std::vector<double> P;

//    std::vector<std::vector<double>> sim {{10., 200, 0., 0.005},{3, 350, 0., 0.01},{5, 750, 0., 0.02}};
    std::vector<std::vector<double>> sim {{10., 200, 0., 0.005},{3, 350, 0., 0.01},{5, 750, 0., 0.02},{1, 700, 0., 0.01}};
    for (size_t i = 0; i < sim.size(); ++i) {
        sim[i][1] *= 2*M_PI;
        sim[i][3] *= sim[i][1];
    }
    for (size_t i = 0; i < sim.size(); ++i)
    	for (size_t j = 0; j < numofdims; ++j)
    		std::cout << std::fixed << sim[i][j] << std::endl;
    std::vector<double> simulated;
    calc_response(sim, numofsamples, ts, simulated);

	std::vector<double> *used_data = &channel0_real;
//	std::vector<double> *used_data = &simulated;
	fft(*used_data, A, P);

	size_t max_A_idx = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
	size_t max_D_idx = std::distance(used_data->begin(), std::max_element(used_data->begin(), used_data->end()));
	std::cout << "Original dominating amp: " << A[max_A_idx] << std::endl;
	std::cout << "Original dominating amp: " << used_data->at(max_D_idx) << std::endl;
	std::cout << "Original dominating freq: " << freq[max_A_idx] << std::endl;
    std::cout << "Number of particles: " << numofparticles << std::endl;
    std::cout << "Number of iterations: " << numofiterations << std::endl;

	size_t t_num = std::thread::hardware_concurrency();
//	t_num = 4;
	t_num = modes_num;
	PSO::freq_set = std::vector<bool> (t_num-1, false);
	PSO::freq_ready = std::vector<bool> (t_num-1, false);
	//std::cout << "Number of detected threads: " << t_num << std::endl;

	std::vector< std::vector<double> > init(numofparticles);
	std::vector<double> xmin, xmax;
	init_minmax(xmin, xmax, numofdims, channel0_real, fs);//freq[freq.size()-1]);

    std::vector<std::thread> threads(t_num);
    std::vector<std::thread> helper_threads(t_num);

    std::vector<std::vector<double>> *founds = new std::vector<std::vector<double>>(t_num);
    for (size_t i = 0; i < t_num; ++i)
    	founds->at(i) = std::vector<double>(numofdims);
	std::vector<std::mutex> mutexes(t_num);
	std::vector<std::condition_variable> cvs(t_num);

    std::vector<PSO*> *psos = new std::vector<PSO*>();
	psos->push_back(new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, 0, founds, &mutexes, &cvs));
    for (size_t p = 1; p < t_num; ++p) {
    	PSO *pso = new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, p, founds, &mutexes, &cvs);
    	psos->push_back(pso);
    	threads[p] = std::thread(&PSO::run, std::ref(*pso));
    }

    psos->at(0)->run();
    for (size_t p = 1; p < psos->size(); ++p) {
    	threads[p].join();
    }

	fft(*used_data, A, P);
	std::vector<double> A_gauss(gaussian_filter(A, 5));
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
//		PSO *helper_pso = new PSO(numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), numofiterations, p, founds, &mutexes, &cvs, f_l, f_r, std::pair<bool, size_t> (true, p));
		PSO *helper_pso = new PSO(1.5*numofparticles, numofdims, xmin, xmax, &time_real, &(*used_data), 1.5*numofiterations, p, founds, &mutexes, &cvs, f_l, f_r, std::pair<bool, size_t> (true, p));
		helper_psos->push_back(helper_pso);
		helper_threads[p] = std::thread(&PSO::run, std::ref(*helper_pso));
	}
    for (size_t p = 0; p < helper_psos->size(); ++p) {
		helper_threads[p].join();
    }

    std::cout << "-------" << std::endl;
    std::cout << "main swarms:" << std::endl;
    for (size_t p = 0; p < t_num; ++p) {
        std::cout << "-------" << std::endl;
        std::cout << "swarm: " << p << std::endl;
    	std::cout << "amp: " << psos->at(p)->getgbest().at(0) << std::endl;
    	std::cout << "freq: " << psos->at(p)->getgbest().at(1)/2/M_PI << std::endl;
    	std::cout << "phase: " << psos->at(p)->getgbest().at(2) << std::endl;
    	std::cout << std::fixed;
    	std::cout << "damping: " << psos->at(p)->getgbest().at(3)/psos->at(p)->getgbest().at(1) << std::endl;
        std::cout << "fitness: " << psos->at(p)->getgbestfit() <<endl;
    }
    std::cout << "-------" << std::endl;
    std::cout << "helpers:" << std::endl;
    for (size_t p = 0; p < helper_psos->size(); ++p) {
        std::cout << "-------" << std::endl;
    	std::cout << "amp: " << helper_psos->at(p)->getgbest().at(0) << std::endl;
    	std::cout << "freq: " << helper_psos->at(p)->getgbest().at(1)/2/M_PI << std::endl;
    	std::cout << "phase: " << helper_psos->at(p)->getgbest().at(2) << std::endl;
    	std::cout << std::fixed;
    	std::cout << "damping: " << helper_psos->at(p)->getgbest().at(3)/helper_psos->at(p)->getgbest().at(1) << std::endl;
        std::cout << "fitness: " << helper_psos->at(p)->getgbestfit() <<endl;
    }


    ofstream myfile;
    myfile.open (file_name + ".log", ios::out | ios::app);
    for (size_t i = 0; i < helper_psos->size(); ++i) {
    myfile  << helper_psos->at(i)->getgbest().at(0) << ";"
			<< helper_psos->at(i)->getgbest().at(1)/2/M_PI << ";"
			<< helper_psos->at(i)->getgbest().at(2) << ";"
			<< helper_psos->at(i)->getgbest().at(3)/helper_psos->at(i)->getgbest().at(1) << ";";
    }

    myfile << std::endl;
    myfile.close();

    std::vector<std::vector<double>> results;
    std::vector<double> response;
    for (size_t i = 0; i < helper_psos->size(); ++i)
    	results.push_back(helper_psos->at(i)->getgbest());
    calc_response(results, numofsamples, ts, response);

    myfile.open (file_name + "_response_in_time.log", ios::out);
    for (size_t i = 0; i < numofsamples; ++i) {
    	double real = used_data->at(i);
    	double resp = response[i];
    	myfile  << real << ";" << resp;
    	myfile  << std::endl;
    }
    myfile.close();

#ifdef PLOT
    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
    std::vector<double *> A_approx(helper_psos->size()+1);

    for (size_t i = 0; i < helper_psos->size(); ++i) {
		A_approx[i] = new double[numofsamples_2];
		approximate_amp(
				std::vector<std::vector<double>> {helper_psos->at(i)->getgbest()},
				time_real,
				numofsamples,
				A_approx[i]);
    }

	A_approx[helper_psos->size()] = new double[numofsamples_2];
	std::vector<std::vector<double>> factors;
	for (size_t i = 0; i < helper_psos->size(); ++i)
		factors.push_back(helper_psos->at(i)->getgbest());
//		factors.push_back(psos->at(i)->getgbest());
    approximate_amp(
    		factors,
			time_real,
			numofsamples,
			A_approx[helper_psos->size()]);

    double y[numofsamples_2];
    double x[numofsamples_2];
	for(size_t i = 0; i < numofsamples_2; ++i) {
		x[i] = freq[i];
		y[i] = 20*log(A[i]);
	}
    ScatterPlot *w = new ScatterPlot(x, y, (size_t)numofsamples_2);

    w->setXArray(freq[0], freq[numofsamples_2-1]);
    //std::cout << "min: " << min_val << std::endl;
    //w->setYArray(min_val, 1.1*20*log(A[max_A_idx]));
    w->setYArray(-300, 100);
    for (size_t i = 0; i < helper_psos->size(); ++i) {
    	w->addData(x, A_approx[i]);
    }
    w->addData(x, A_approx[helper_psos->size()]);

    for (size_t i = 0; i < psos->size(); ++i) {
    	delete psos->at(i);
    }
    for (size_t i = 0; i < helper_psos->size(); ++i) {
    	delete helper_psos->at(i);
    }
    delete founds;
#endif

#ifdef PLOT
    //    double y_[0];
    //    double x_[0];
    //    ScatterPlot *w1 = new ScatterPlot(x_, y_, 0);
    ////    w1->display(argc-4, argv+5);
    //    w1->display(4, argv+5);
    //    delete w1;
    w->display(4, argv+5);
    delete w;
    for (size_t i = 0; i < A_approx.size(); ++i) {
    	delete A_approx[i];
    }
#endif
	}
    catch (std::exception& e)
    {
      cout << e.what() << '\n';
    }
    return 0;
}
