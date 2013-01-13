#ifndef _STATE_H
#define _STATE_H

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

class State{
	private:
	    int _dim; // the number of oscillators
	    int * _data;
	public:
	    State(const int & dim){
		Init(dim);
	    }
	    void  Init( int n){
		_data = new int[n + 1];
		_dim = n;
	    
	    }
	    State(){
		_dim=0;
	    }
	    ~State(){
	//	cerr << "Clearing state info" <<endl;
		if (_dim > 0 ) delete []_data;
	    }
	    State (const State &a){
		_dim = a._dim;
		_data = new int[_dim+1];
		for (int i=0;i<_dim+1;++i) _data[i] = a._data[i];
	    }
	
	    int  ElState() {
		return *_data;
	    }
	    int * NuclearVib() {
		return _data+1;
	    }
	    const int & NLevels(){ 
		return _dim;
	    }

	    void SetEl(const int & el){
		_data[0] = el;
	    }
	    void SetNucl (int * nu ){
	    	for (int i=0; i<_dim; i++) _data[i+1] = nu[i];
	    }
	    const int & GetNucl (const int & i ){
	    	return _data[i+1];
	    }

	    void ChangePairNus(const int &  elfrom, const int & elto, const int & nufrom, const int & nuto){
		ChangeNu(elfrom, nufrom);
		ChangeNu(elto, nuto);
	    }

	    void ChangeNu(const int &elfrom, const int & nufrom){
		_data[elfrom+1] = nufrom;
	    }

	    int NuExcited(){
		int res=0;
	    	for (int i =1; i<_dim+1 ;++i){
			if (_data[i] != 0) res++;
		}
		return res;
	    }
	    int TotNu(){
		int res=0;
	    	for (int i =1; i<_dim+1 ;++i){
		    res += _data[i];
		}
		return res;
	    }
	    void DeExciteRandom(){
		vector <int> exc;
		for (int i =1; i<_dim+1 ;++i){
			if (_data[i] != 0) exc.push_back(i);
		}
		int choice = floor( drand48() * exc.size() );
		if (exc.size() != 0 ){
		    _data[exc[choice]]--;
		}
	    }


};

#endif
