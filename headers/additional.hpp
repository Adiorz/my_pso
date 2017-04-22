#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>

#ifndef ADDITIONAL_H_
#define ADDITIONAL_H_

void donothing(std::vector<float> *realdata, std::vector<float> &fitnesses, std::vector< std::vector<float> > &X);

static bool abs_compare(int a, int b);

void init_minmax(std::vector<float> &xmin, std::vector<float> &xmax, size_t numofdims, std::vector<double> &data);
void get_frequencies(std::vector<float> &freq, size_t numofsamples, float fs);

void calc_response(std::vector<std::vector<float>> results, size_t numofsamples, float ts, std::vector<float> &response);
void approximate_amp(std::vector<std::vector<float>> factors, std::vector<double> &time_real, size_t numofsamples, double *out);

class DataStream {
public:
	DataStream(std::string file_name, size_t column, size_t number);
	size_t size();
	const std::vector<double> &get();

	const double operator[] (size_t n) const;
	double & operator[](size_t n);
private:
	std::vector<double> data;
	void parseFile(std::ifstream& file, size_t column, size_t number);
};

double min(double *values, size_t size);

void findMinimas(std::vector<float> &v, size_t start, size_t end, std::vector<size_t> &idx);

void findMinima(std::vector<float> &v, size_t maxIDx, size_t &idxL, size_t &idxR);

std::vector<float> gaussian_filter(std::vector<float> &input, size_t sigma);

bool should_skip(size_t f, std::vector<std::pair<size_t, size_t>> &to_skip);

#endif /* ADDITIONAL_H_ */
