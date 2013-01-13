#ifndef _NEIGHBOURS_H
#define _NEIGHBOURS_H

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>



using namespace std;

#include "State.h"
#include "Parameters.h" 
#include "global.h"

//because the Franck Condon factors are hard coded, might as well use static arrays
const int  NMODE=20;

class Neighbours{
    private:

	   struct _NeighType{
	   	int _ElFrom;
	        int _ElTo;
		pair <int, int> _NuFrom;
	   };
	   struct _Neigh{
		vector< pair <double, pair<int, int> > > _Rates;
		double _totRate;
           };


	   //these two containers will contain all my precomputed data
	   static vector <_NeighType * > _AllNeighTypes;
           // the first pair gives the electronic states (from and to) and the second pair gives the nuclear states (to and from)
	   static map < pair <pair <int, int>, pair<int, int> > , _Neigh*  > _AllRates;	
	   static double ** _FC;
	   static map < pair< int, int> ,  vector<double> > _AllProbs;
	   static map < pair<int, int> , double **> _PreComRates; 

	   bool init;
	   bool initFC;
	   bool initPreComRates;
		
	   void ClearFC();
	   void ClearPreComRates();

	   // compute all the hashed rates etc
	   // the silly container filler
	   void ComputeAllRates(Parameters *);
	   void InitFC(Parameters *, const int &n);
	   void InitPreComRates( Parameters *);
	   void InitAllProbs();

	   // the actual rate computing algs
	   double ComputeRate(Parameters *, int fromS, int toS );
	   double Prefactor(Parameters *, double);

	   double Exponential(Parameters * , const int &, const int & , const int &, const int & );
//	   vector <double> Prob(const int &, const int &);
	   double Prob(const int &, const int &, const int &);

	   //sort and prune the neighbour list
	   void SortAndPrune(map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itRat);
	   void TotRate(map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itRat);

    public:
	   vector <double> Prob(const int &, const int &);
	   Neighbours(){
	   	init = false;
		initFC = false;
		initPreComRates=false;
	   }
	   Neighbours( Parameters *para ){
	   	init = false;
		initFC = false;
		initPreComRates=false;
		Init( para);
	   }

	   ~Neighbours();
				
	   void  Init( Parameters* );

	   // The important accesor routine is to get the totalrate for a give elfrom, elto, nufrom and nuto. 
	   const double & TotRate(const int & elfrom, const int & elto, const int & nufrom, const int & nuto);
	   const pair<int, int> &  RandomDestination(const int & elfrom, const int & elto, const int & nufrom, const int & nuto);
           
           //either returns the rate if precompute or computes it
           map < pair <pair <int, int>, pair<int, int> > , Neighbours::_Neigh*  >::iterator  GetItRate(const int & elfrom, const int & elto, const int & nufrom, const int & nuto);


	   // Debug routines
	   void PrintAllNeighTypes(ostream &);
	   void PrintAllRates(ostream &);
	   void PrintProbs(ostream &);
	   void PrintCalcRates(ostream & );


};

#endif
