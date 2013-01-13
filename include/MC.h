
#ifndef _MC_H
#define _MC_H

#include <stdlib.h>
#include <ctime>

using namespace std;

#include "State.h"
#include "Parameters.h" 
#include "global.h"
#include "Neighbours.h"

class MC{
    private:
	map < int,  vector <int>* > _NList;
	double _time;
        bool _run;
	int _tbin;
	Neighbours * _Rates;
	double _ratevib;
	vector <  vector<double> *  >  _elstate;
	vector <  vector<double> *  > _vibnrg;
	double _tmax;
	double _dt;
	/// the ratio of successive time increments
	double _tau;
	vector <double> _timebins;
        vector <double> _CollectionTimes;
        vector <double> _RecombTimes;
	void IncrementTime(const double &, State *);
        vector <string> _prop;
        bool _recordHistograms;
        bool _randomizeSeed;
        int _totRun;
        

    public:
	MC(){}
	MC(Parameters * para, Neighbours * neigh);
	~MC();
	const double & Time() {
	    return _time;
	}

	void Init(Parameters *);
	void NextState(State *);
	void Run(const State &,const int & n, string);
	void Run(State *);
	friend ostream& operator<< (ostream &, MC &);

        double PCollected(){
            return (double) _CollectionTimes.size()/_totRun;
        }
        double PRecomb(){
            return (double) _RecombTimes.size()/_totRun;
        }

};


#endif
