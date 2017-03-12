#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>

#ifndef ADDITIONAL_H_
#define ADDITIONAL_H_

//void calculaterealdata(std::vector<float> &target, size_t numofsamples);
void donothing(std::vector<float> *realdata, std::vector<float> &fitnesses, std::vector< std::vector<float> > &X);

static bool abs_compare(int a, int b);

//std::vector<double> first_derivative(std::vector<double> &data, double step);

void init_minmax(std::vector<float> &xmin, std::vector<float> &xmax, size_t numofdims, std::vector<double> &data);
void get_frequencies(std::vector<float> &freq, size_t numofsamples, float fs);

void approximate_amp(float amp, float omega, float phase, float bump, std::vector<double> &time_real, size_t numofsamples, double *out);

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

#endif /* ADDITIONAL_H_ */
