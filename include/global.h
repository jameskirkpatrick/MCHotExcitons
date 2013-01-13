#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

const double hbar = 6.57E-16;
const double Pi=3.141592653589793;

int factorial(int);
double factln(const int &);
double binomial(const int &, const int &);
int KrDe(const int &, const int &);
void normalise(vector<double> *);
int min(const int &, const int &);
int max(const int &, const int &);
int invGeom(const double & res, const double & a, const double & r);
double Geom(const int & n, const double & a, const double & r);
double Gaussian(const double &mu=0.0, const double &sigma=1.0);

#endif
