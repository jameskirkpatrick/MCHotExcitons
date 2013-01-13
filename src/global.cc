#include "global.h"

int min(const int &a, const int & b){
	return a <= b ? a : b;
}

int max(const int &a, const int & b){
	return a >= b ? a : b;
}
int factorial(int n){
    static int ntop=4;
    static double a[33] = {1.,1.,2.,6.,24.};
    int j;
    if (n<0) cerr << "negative factorial???" <<endl;
    if (n>32) cerr << "Overflow factorial" <<endl;
    while (ntop <n){
    	j=ntop++;
	a[ntop]=a[j]*ntop;
    }
    return a[n];
}

double gammaln (const double & xx){
	int j;
	double x,y,tmp,ser;
	static const double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 
	0.1208650973866179E-02, -0.5395239384953E-5};
	y=x=xx;
	tmp = x+5.5;
	tmp -= (x+0.5) *log(tmp);
	ser = 1.000000000190015;
	for (j=0; j<6;j++) ser +=cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double factln(const int & n){
    static double a[101];

    if (n<0) cerr << "factorial of a negative numero" <<endl;
    if (n<=1) return 0.0;
    if (n<=100) return (a[n] != 0.0 ? a[n] : (a[n] = gammaln(n+1.)));
    else return gammaln(n+1.);
}

double binomial(const int &n, const int &k){
	return floor(0.5+exp(factln(n)- factln(k) - factln(n-k)));
}

int KrDe(const int &a , const int &b){
	if(a==b) return 1;
	else return 0;
}

void normalise (vector <double> * a){
	vector <double>::iterator ita=a->begin();
	double max = 0.;
	for ( ; ita!= a->end(); ita++) max += *ita;
//	cout << " max " << max <<endl;
        for (int i=0 ; i <  a->size(); i++) a->at(i) = a->at(i)/max;

}


int invGeom(const double & res, const double & a, const double & r){
	double rn = 1 - res /(a*(1-r)) 	;
	return floor(log(rn) / log(r));
}

double Geom(const int & n, const double & a, const double & r){
	return a * (1-pow(r,n)) * (1-r);
}

/******************************************************************************/
/* randn()
 *
 * Normally (Gaussian) distributed random numbers, using the Box-Muller
 * transformation.  This transformation takes two uniformly distributed deviates
 * within the unit circle, and transforms them into two independently
 * distributed normal deviates.  Utilizes the internal rand() function; this can
 * easily be changed to use a better and faster RNG.
 *
 * The parameters passed to the function are the mean and standard deviation of
 * the desired distribution.  The default values used, when no arguments are
 * passed, are 0 and 1 - the standard normal distribution.
 *
 *
 * Two functions are provided:
 *
 * The first uses the so-called polar version of the B-M transformation, using
 * multiple calls to a uniform RNG to ensure the initial deviates are within the
 * unit circle.  This avoids making any costly trigonometric function calls.
 *
 * The second makes only a single set of calls to the RNG, and calculates a
 * position within the unit circle with two trigonometric function calls.
 *
 * The polar version is generally superior in terms of speed; however, on some
 * systems, the optimization of the math libraries may result in better
 * performance of the second.  Try it out on the target system to see which
 * works best for you.  On my test machine (Athlon 3800+), the non-trig version
 * runs at about 3x10^6 calls/s; while the trig version runs at about
 * 1.8x10^6 calls/s (-O2 optimization).
 *
 *
 * Example calls:
 * randn_notrig();	//returns normal deviate with mean=0.0, std. deviation=1.0
 * randn_notrig(5.2,3.0);	//returns deviate with mean=5.2, std. deviation=3.0
 *
 *
 * Dependencies - requires <cmath> for the sqrt(), sin(), and cos() calls, and a
 * #defined value for PI.
 */

/******************************************************************************/


/******************************************************************************/
//	Standard version with trigonometric calls
#define PI 3.14159265358979323846

double Gaussian(const double &mu, const double &sigma) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;//	deviate from previous calculation
	double dist, angle;

	//	If no deviate has been stored, the standard Box-Muller transformation is
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {

		//	choose a pair of uniformly distributed deviates, one for the
		//	distance and one for the angle, and perform transformations
		dist=sqrt( -2.0 * log(double(drand48()) ) );
		angle=2.0 * PI * (double(drand48()) );

		//	calculate and store first deviate and set flag
		storedDeviate=dist*cos(angle);
		deviateAvailable=true;

		//	calcaulate return second deviate
		return dist * sin(angle) * sigma + mu;
	}

	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}