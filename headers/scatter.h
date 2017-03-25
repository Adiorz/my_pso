/*
 * scatter.h
 *
 *  Created on: Jan 21, 2017
 *      Author: adiorz
 */

#ifndef HEADERS_SCATTER_H_
#define HEADERS_SCATTER_H_

#include "plc++demos.h"
#include "plplot.h"
#include <vector>

class ScatterPlot {
public:
	ScatterPlot(PLFLT *x, PLFLT *y, size_t n);
	void display(int, const char **);
    size_t getN() { return n; }
    void addData(PLFLT *x, PLFLT *y);
    void addPoint(PLFLT x, PLFLT y);
    void setXArray(PLFLT xmin = 0., PLFLT xmax = 10.0);
    void setYArray(PLFLT ymin = -10., PLFLT ymax = 10.0);
private:
    // Class data
    plstream         *pls;

    size_t n;
    std::vector<PLFLT *> x;
    std::vector<PLFLT *> y;
    PLFLT xmin;
    PLFLT xmax;
    PLFLT ymin;
    PLFLT ymax;
};

void ScatterPlot::setXArray(PLFLT xmin, PLFLT xmax) {
	this->xmin = xmin;
	this->xmax = xmax;
}
void ScatterPlot::setYArray(PLFLT ymin, PLFLT ymax) {
	this->ymin = ymin;
	this->ymax = ymax;
}
void ScatterPlot::addData(PLFLT *x, PLFLT *y) {
	this->x.push_back(x);
	this->y.push_back(y);
}

void ScatterPlot::addPoint(PLFLT x, PLFLT y) {
	double X[1] = {x};
	double Y[1] = {y};
	pls->poin(6, X, Y, 9);
	plcol0(5);
}

void ScatterPlot::display(int argc, const char **argv) {
    pls = new plstream();

    // Parse and process command line arguments
    pls->parseopts( &argc, argv, PL_PARSE_FULL );

    // Initialize plplot
    // change background colour to grey
    plscolbg(128,128,128);
    pls->init();
    // change axes colour to black
    //plscol0(0, 0, 0, 0);
    plcol0(15);

    // Create a labelled box to hold the plot.
    pls->env(xmin, xmax, ymin, ymax, 0, 0);
    //pls->lab("x", "y=100 x#u2#d", "Simple PLplot demo of a 2D line plot");
    pls->lab("freq", "amp", "Spectrum of found mode");

    // Plot the data that was prepared above.
    for (size_t i = 0; i < x.size(); ++i) {
    	//plscol0(0, 0, 255, 0);
		plcol0(i+1);
    	pls->line(n, x.at(i), y.at(i));
    }

    // In C++ we don't call plend() to close PLplot library
    // this is handled by the destructor
    delete pls;

}

ScatterPlot::ScatterPlot(PLFLT *x, PLFLT *y, size_t n)
{
	this->n = n;
	this->x = std::vector<PLFLT *>();
	this->y = std::vector<PLFLT *>();
	addData(x, y);
}

#endif /* HEADERS_SCATTER_H_ */
