#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>

#ifndef ADDITIONAL_H_
#define ADDITIONAL_H_

void donothing(std::vector<double> *realdata, std::vector<double> &fitnesses, std::vector< std::vector<double> > &X);

void init_minmax(std::vector<double> &xmin, std::vector<double> &xmax, size_t numofdims, std::vector<double> &data, size_t max_freq = -1);
void get_frequencies(std::vector<double> &freq, size_t numofsamples, double fs);

void calc_response(std::vector<std::vector<double>> results, size_t numofsamples, double ts, std::vector<double> &response);
void approximate_amp(std::vector<std::vector<double>> factors, std::vector<double> &time_real, size_t numofsamples, double *out);

class DataStream {
public:
	DataStream(std::string file_name, size_t column, size_t number, size_t to_skip);
	size_t size();
	const std::vector<double> &get();

	const double operator[] (size_t n) const;
	double & operator[](size_t n);
private:
	std::vector<double> data;
	void parseFile(std::ifstream& file, size_t column, size_t number, size_t to_skip);
};

double min(double *values, size_t size);

void findMinimas(std::vector<double> &v, size_t start, size_t end, std::vector<size_t> &idx);

void findMinima(std::vector<double> &v, size_t maxIDx, size_t &idxL, size_t &idxR);

std::vector<double> gaussian_filter(std::vector<double> &input, size_t sigma);

bool should_skip(size_t f, std::vector<std::pair<size_t, size_t>> &to_skip);

#endif /* ADDITIONAL_H_ */
