#ifndef PARAMETERS_H
#define PARAMETERS_H


#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>


#include "State.h"
#include "EnergyFunction.h"
#include "NonIsoPol.h"
#include "libxml/parser.h"

class Parameters{
    private:


        struct _point{
            const _point & operator = (_point  &r);
            bool operator == ( const _point & r){
            return this->_x1 == r._x1 && this->_y1 == r._y1 &&this->_z1 == r._z1
                 &&this->_x2 == r._x2 && this->_y2 == r._y2 &&this->_z2 == r._z2;
            }
            
            int _x1,_y1,_z1,_x2,_y2,_z2;
        };
       
        

        vector <_point> neigh(_point &);
        vector <_point> neigh1Crg(_point &); // no pbc in xy
        vector <_point> neigh1Crg(_point &, const int & a); // with pbc in xyz

	void Read(const char *, const char *);


        /// still too lazy to write the accessors :D
    
	// the list of energies in state [i]  (in eV)
	vector <double> _energies;
	//The list of electronic states connected  to each other
	vector < pair <int, int> > _top;
	vector < double > _transfers;
        // a list of properties for each state
        vector <string> _properties;
	//the number of nuclear modes per site 
	int _nmodemax;
        /// the disorder in log(J)
        double _sigmaJ;
	//the Huang Rhys factor
	double _S;
	//the energy of a phonon
	double _hbarOmega;
	//the outer sphere reorganization energy
	double _lambdaO;
	// temperature
	double _kT;
        ///the stopping time for simulation
        double _tmax;
        ///the deltat for MC
        double _dt;
	///the ratio of successive intervals
	double _tau;
        /// the rate for vibrational relaxation
        double _ratevib;
        /// the number of phonons to start with
        string _nstartphonon;
        /// A label for the run
        int _label;

        friend class MC;
        friend class Neighbours;

        void Read(const char *, const char *,const char *);
	//create the state which is labelled "Generator"
        
        void Init2Charges(const int & n, const double &d,
        const double &eps, const double &F,const double &J, const double &Jrec, const double & gsnrg  );
        void Init1D(const int & n, const double &d,
        const double &eps, const double &F,
        const double &J, const double &Jrec,
        const double & gsnrg  );

        void Init1DS(const int & n, const double &d,
        const double &eps, const double &kT,
        const double &J, const double &Jrec,
        const double & gsnrg  );
        void Init3D(const int & n, const double & a, EnergyFunction & EF, const double & J,
        const double & Jrec, const double & gsnrg);

public:
        ///Change this to allow for changing values
	Parameters(){
	    _S =1 ;
	    _kT=0.025;
	    _lambdaO=0.1;
	    _hbarOmega=0.2;
            _tmax=1E-10;
            _dt =1E-16;
	    _tau =2;
            _ratevib=1E14;
        }
	

	~Parameters(){
	//    cerr << "Clearing parameters" <<endl;
	    _energies.clear();
	    _top.clear();	
	}

        void Init(string);
        State NN(const int &n);
        

        string GetNStartPhonon() {
            return _nstartphonon;
        }

        void Print(string &);

};

#endif 


