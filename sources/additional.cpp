#include "../headers/additional.hpp"
#include "../headers/fft.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

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

//void findMaximas(std::vector<double> &v, std::vector<size_t> &idx) {
//	for (unsigned int i = 1; i < v.size(); ++i) {
//		if (v[i] > v[i-1]) {
//		   unsigned int j = i;
//		   while (v[j] == v[j+1])
//			   ++j;
//		   ++j;
//		   if (v[j] < v[i]) {
//			   idx.push_back(i);
//			   i = j;
//		   }
//		}
//	}
//}

void findMinimas(std::vector<float> &v, size_t start, size_t end, std::vector<size_t> &idx) {
	for (unsigned int i = start+1; i < end; ++i) {
	   if (v[i] < v[i-1]) {
		   unsigned int j = i;
		   while (v[j] == v[j+1]) {
			   ++j;
		   }
		   ++j;
		   if (v[j] > v[i]) {
			   idx.push_back(i);
			   i = j;
		   }
	   }
	}
}

void findMinima(std::vector<float> &v, size_t maxIDx, size_t &idxL, size_t &idxR) {
   if (maxIDx <= 3)
       idxL = 0;
   else {
       std::vector<size_t> minsIDXLeft;
       std::vector<size_t> minsIDXLeft2;
       std::vector<float> minsLeft;
       findMinimas(v, 0, maxIDx, minsIDXLeft);
       if (minsIDXLeft.size() <= 3)
           idxL = minsIDXLeft[minsIDXLeft.size()-1];
       else {
           for (size_t i = 0; i < minsIDXLeft.size(); ++i) {
               minsLeft.push_back(v[minsIDXLeft[i]]);
           }
           findMinimas(minsLeft, 0, minsLeft.size(), minsIDXLeft2);
           idxL = minsIDXLeft[minsIDXLeft2[minsIDXLeft2.size() - 1]];
       }
   }

   if (v.size() - maxIDx <= 3)
       idxR = 0;
   else {
       std::vector<size_t> minsIDXRight;
       std::vector<size_t> minsIDXRight2;
       std::vector<float> minsRight;
       findMinimas(v, maxIDx, v.size(), minsIDXRight);
       //std::cout << minsIDXRight.size() << std::endl;
       if (minsIDXRight.size() <= 3) {
    	   if (minsIDXRight.size() == 0)
    		   idxR = 0;
    	   else
    		   idxR = minsIDXRight[minsIDXRight.size()-1];
       }
       else {
           for (size_t i = 0; i < minsIDXRight.size(); ++i) {
               minsRight.push_back(v[maxIDx+minsIDXRight[i]]);
           }
           findMinimas(minsRight, 0, minsRight.size(), minsIDXRight2);
           idxR = minsIDXRight[minsIDXRight2[minsIDXRight2.size() - 1]];
       }
   }
}

std::vector<float>  gaussian_filter(std::vector<float> &input, size_t sigma) {
	size_t n = 6*sigma;
	std::vector<float> kernel(n);
	double half_n = (n - 1) / 2.0;
	double sum = 0.0;
	for (size_t i = 0;  i < n;  ++i)
	{
	    double x = i - half_n;
	    kernel[i] = 1.0/sqrtf(2*M_PI*sigma*sigma)*exp(-x*x/(2*sigma*sigma));// use formula with x
	    //std::cout << kernel[i] << std::endl;
	    sum += kernel[i];
	}
    std::cout << sum << std::endl;
	//std::vector<float> kernel = {0.006, 0.061, 0.242, 0.383, 0.242, 0.061, 0.006};// find center position of kernel (half of kernel size)
	std::vector<float> output = std::vector<float>(input.size());
	size_t cols = input.size();
	size_t kCols = kernel.size();
	size_t kCenterX = kernel.size() / 2;
	for(size_t j=0; j < cols; ++j)          // columns
	{
		for(size_t n=0; n < kCols; ++n) // kernel columns
		{
			size_t nn = kCols - 1 - n;  // column index of flipped kernel

			// index of input signal, used for checking boundary
			size_t jj = j + (n - kCenterX);

			// ignore input samples which are out of bound
			if(jj >= 0 && jj < cols )
				output[j] += input[jj] * kernel[nn];
		}
	}
	return output;
}
