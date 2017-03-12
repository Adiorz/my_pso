#include <vector>
#include <ctime>
#include <random>

std::vector<double> first_derivative(std::vector<double> &data, double step);

class PSO {
private:
	size_t id;

	size_t numofparticles;
	size_t numofdims;
	size_t numofiterations;
	size_t numofsamples;

	size_t numofsamples_2;

	std::vector< std::vector<float> > V;
	std::vector< std::vector<float> > X;
	std::vector<float> Xmax;
	std::vector<float> Xmin;
	std::vector<float> Vmax;
	std::vector<float> Vmin;
    std::vector< std::vector<float> > pbests;
    std::vector<float> pbestfits;
    std::vector<float> fitnesses;
    float gbestfit;

	std::vector<float> gbest;
	std::vector<float> worsts;
	std::vector<float> meanfits;
	std::vector<float> bests;
	float c1, c2;

    float w;
    float minfit;
    int   minfitidx;

    std::vector<double> *realdata;
    std::vector<double> *time;
    float fs;	// sampling frequency
    std::vector<float> A;
    std::vector<float> P;

    float init_param(size_t j);

    void init_max_velocities();

    void initpopulation();

    void addparticle();

    float fitnessfunc_singleparticle(size_t p);

    float calc_response(float amp, float omega, float phase, float bump, float t);

    void fitnessfunc();

    void fitnessfunc_thread(size_t start, size_t end);

    void fitnessfunc_multi();

    void calcgbest(bool first = false);

    void update();

    float update_velocity(float w, float X, float V, float V_min, float V_max,
    		float gbest, float pbests, float c1, float c2);

    float update_position(float X, float V, float X_min, float X_max);

    void do_nothing();

public:
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution;
	size_t skip_low;
	size_t skip_high;

	PSO(	size_t numofparticles,
			size_t numofdims,
			std::vector<float> &Xmin,
			std::vector<float> &Xmax,
			std::vector<double> *time,
			std::vector<double> *realdata,
			size_t numofiterations,
			size_t id,
			float c1 = 2, float c2 = 2);

//	PSO(const PSO& c);

//	PSO & operator=(const PSO &) {
//		std::cout << "Assignment operator" << std::endl;
//	}

	void run();

	std::vector<float> getgbest();

	float getgbestfit();
};