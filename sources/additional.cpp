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

void init_minmax(std::vector<float> &xmin, std::vector<float> &xmax, size_t numofdims, std::vector<double> &data, size_t max_freq) {
    xmin = std::vector<float>(numofdims);
    xmax = std::vector<float>(numofdims);
    std::vector<double>::iterator max_abs_it;
    max_abs_it = std::max_element(data.begin(), data.end(), abs_compare);
	size_t max_abs_idx = std::distance(data.begin(), max_abs_it);
	float max_abs = std::abs(data[max_abs_idx]);
    xmin[0] = 0; //amp
    xmax[0] = 2*max_abs;
    xmin[1] = 2*20*M_PI; //omega
    if (max_freq == -1)
    	xmax[1] = 2*3000*M_PI;//359.99f;
    else
    	xmax[1] = 2*M_PI*max_freq;
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

void calc_response(std::vector<std::vector<float>> results, size_t numofsamples, float ts, std::vector<float> &response) {
	response = std::vector<float> (numofsamples, 0.f);
	for (size_t i = 0; i < results.size(); ++i) {
		float amp = results[i][0];
		float omega = results[i][1];
		float phase = results[i][2];
		float bump = results[i][3];
		for (size_t j = 0; j < numofsamples; ++j) {
			double t = j*ts;
			response[j] += amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*t+phase)*
					exp(-bump*t);
		}
	}
}

void calc_response(float amp, float omega, float phase, float bump, size_t numofsamples, float ts, std::vector<double> &response) {
	response = std::vector<double> (numofsamples, 0.f);
	for (size_t j = 0; j < numofsamples; ++j) {
		double t = j*ts;
		response[j] += amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*t+phase)*
				exp(-bump*t);
	}
}

float fitnessfunc(std::vector<float> &p, size_t numofsamples, float ts, std::vector<float> A_target, float max_A) {
//	std::cout << "entered" << std::endl;
//	while (p[3] > p[1]) {
//		float high = Xmax[3];
//		float low = Xmin[3];
//		float rand = distribution(generator);
//		float init_position = rand * (high-low) + low;
//		p[3] = init_position;
//	}

	float amp = p[0];
	float omega = p[1];
	float phase = p[2];
	float bump = p[3];

	std::vector<double> response(numofsamples, 0);
	calc_response(amp, omega, phase, bump, numofsamples, ts, response);

	std::vector<float> A;
	std::vector<float> P;
	fft(response, A, P);

	size_t numofsamples_2 = size_t(numofsamples/2)+1;

	float fitness = 0.f;
	std::vector<std::vector<size_t>> to_skip;
	to_skip.push_back(std::vector<size_t>{52, 216, 336});
	to_skip.push_back(std::vector<size_t>{367, 413, 764});
	for(size_t j = 0; j < numofsamples_2; ++j) {
		float residue = A_target[j] - A[j];
		float temp_fit = residue*residue*residue*residue/max_A;
		size_t l_skip;
		size_t m_skip;
		size_t h_skip;
		bool b_skip = false;
//		for (size_t i = 0; i < to_skip.size(); ++i) {
//			if (j >= to_skip[i][0] && j <= to_skip[i][2]) {
//				b_skip = true;
//				l_skip = to_skip[i][0];
//				m_skip = to_skip[i][1];
//				h_skip = to_skip[i][2];
//				break;
//			}
//		}
//		if (b_skip) {
//			float penalty = 0.;
//			if (j < m_skip) {
//				penalty = (j-l_skip)/(float)(m_skip-l_skip);
//			}
//			else {
//				penalty = (h_skip-j)/(float)(h_skip-m_skip);
//			}
//			temp_fit *= abs(A[j])*(1.f + penalty);
//		}
//		else
//			std::cout << "no add penalty: " << j << std::endl;
		fitness += temp_fit;
	}
	return fitness;
}

void approximate_amp(
		float amp, float omega, float phase, float bump,
		std::vector<double> &time_real, size_t numofsamples, double *out) {
	std::vector<double> appr(numofsamples, 0.0);
	for(size_t i = 0; i < numofsamples; ++i) {
		appr[i] += amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*time_real[i]+phase) * exp(-bump*time_real[i]);
    }

    std::vector<float> A_appr;
    std::vector<float> P_appr;
    fft(appr, A_appr, P_appr);
    size_t numofsamples_2 = (size_t)(numofsamples/2)+1;
	for(size_t i = 0; i < numofsamples_2; ++i)
		out[i] = 20*log(A_appr[i]);
}

void approximate_amp(
		std::vector<std::vector<float>> factors,
		std::vector<double> &time_real, size_t numofsamples, double *out) {
	std::vector<double> appr(numofsamples, 0.0);
	for (size_t i = 0; i < factors.size(); ++i) {
		float amp = factors[i][0];
		float omega = factors[i][1];
		float phase = factors[i][2];
		float bump = factors[i][3];
		for(size_t j = 0; j < numofsamples; ++j) {
			appr[j] += amp*sin(omega*sqrt(1-(bump/omega)*(bump/omega))*time_real[j]+phase) * exp(-bump*time_real[j]);
	    }
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

void findMinimas(std::vector<float> &v, size_t start, size_t end, std::vector<size_t> &idx) {
	idx = std::vector<size_t> ();
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
	    kernel[i] = 1.0/sqrtf(2*M_PI*sigma*sigma)*exp(-x*x/(2*sigma*sigma));
	    sum += kernel[i];
	}
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

bool should_skip(size_t f, std::vector<std::pair<size_t, size_t>> &to_skip){
	for (size_t i = 0; i < to_skip.size(); ++i)
		if (f >= to_skip[i].first && f <= to_skip[i].second)
			return true;
	return false;
}
