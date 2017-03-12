#include "../headers/additional.hpp"
#include "../headers/fft.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

//void calculaterealdata(std::vector<float> &target, size_t numofsamples) {
//	for(size_t i = 0; i < numofsamples; ++i) {
//		float t = i/samplingrate;
//		target.at(i) = sin(t);
//	}
//}

void donothing(std::vector<float> *realdata, std::vector<float> &fitnesses, std::vector< std::vector<float> > &X) {
//		,
//		size_t start, size_t end) {
	;
}

DataStream::DataStream(std::string file_name, size_t column, size_t number) {
	std::ifstream data_stream(file_name);
	if (data_stream)
		parseFile(data_stream, column, number);
	else
		std::cout << "Cannot open the " << file_name << " file" << std::endl;
}
size_t DataStream::size() {
	return data.size();
}
const std::vector<double> &DataStream::get() {
	return data;
}
const double DataStream::operator[] (size_t n) const {
	return data[n];
}
double & DataStream::operator[](size_t n) {
	return data[n];
}

void DataStream::parseFile(std::ifstream& file, size_t column, size_t number)
{
	size_t read_number = 0;
	for (std::string line; std::getline(file, line); )
	{
		if (read_number >= number)
			break;
		std::istringstream in(line);
		double piece; // of data
		size_t counter = 0;
		while (in >> piece) {
			if (counter == (column)) {
				data.push_back(piece);
				read_number++;
			}
			counter++;
		}
	}
}

static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}

//std::vector<double> first_derivative(std::vector<double> &data, double step) {
//	std::vector<double> derivative(data.size());
//	for (size_t i = 0; i < data.size(); ++i) {
//		if ((i == 0) || (i == data.size() - 1))
//			derivative[i] = data[i];
//		else {
//			double tmp = data[i-1] + data[i+1];
//			tmp /= (2*step);
//			derivative[i] = tmp;
//		}
//	}
//	return derivative;
//}

void init_minmax(std::vector<float> &xmin, std::vector<float> &xmax, size_t numofdims, std::vector<double> &data) {
    xmin = std::vector<float>(numofdims);
    xmax = std::vector<float>(numofdims);
    std::vector<double>::iterator max_abs_it;
    max_abs_it = std::max_element(data.begin(), data.end(), abs_compare);
	size_t max_abs_idx = std::distance(data.begin(), max_abs_it);
	float max_abs = std::abs(data[max_abs_idx]);
    xmin[0] = 0; //amp
    xmax[0] = 2*max_abs;
    xmin[1] = 2*20*M_PI; //omega
    xmax[1] = 2*1500*M_PI;//359.99f;
    xmin[2] = 0; //phase
    xmax[2] = 0.2*M_PI;//1000;
    xmin[3] = 0.0001*xmin[1]; //bumping
    xmax[3] = 0.95*xmax[1];//1000;
}

void get_frequencies(std::vector<float> &freq, size_t numofsamples, float fs) {
	size_t numofsamples_2 = (size_t)(numofsamples/2)+1;

	freq = std::vector<float>(numofsamples_2);

	for (size_t i = 0; i < freq.size(); ++i)
			freq[i] = i*fs/(numofsamples_2-1)/2;
}

void approximate_amp(
		float amp, float omega, float phase, float bump,
		std::vector<double> &time_real, size_t numofsamples, double *out) {
	std::vector<double> appr(numofsamples);
	for(size_t i = 0; i < numofsamples; ++i) {
		appr[i] = amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*time_real[i]+phase);
		appr[i] *= exp(-bump*time_real[i]);
    }

    std::vector<float> A_appr;
    std::vector<float> P_appr;
    fft(appr, A_appr, P_appr);
    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
	for(size_t i = 0; i < numofsamples_2; ++i)
		out[i] = 20*log(A_appr[i]);
}

double min(double *values, size_t size) {
	size_t counter = 0;
	double min = values[0];
	while (isinf(min)) {
		min = values[counter];
		++ counter;
	}
	for (size_t i = counter; i < size; ++i) {
	    if (values[i] < min and !isinf(values[i]))
	        min = values[i];
	}
	return min;
}
