#include "MC.h"

MC::MC(Parameters * para, Neighbours * neigh){
        Init(para);
        _recordHistograms=false;
        _time=0;
        _Rates  = neigh;
        _ratevib=para->_ratevib;
        _tmax = para->_tmax;
        _dt = para->_dt;
	_tau = para->_tau;
        int n = invGeom(_tmax, _dt, _tau);
	//check
	cout << "Initialise " << n << " time bins with ratio " << _tau << " and max time " << Geom(n, _dt, _tau) <<endl;
        for (int i =1; i<=n;i++){
		_timebins.push_back(Geom(i,_dt, _tau));
                vector <double> * pv = new vector < double >;
                vector <double> * pi = new vector <double>;

                _elstate.push_back(pv );
                _vibnrg.push_back( pi );

                //init the times to 0 for starters
                for (int j=0;j< para->_energies.size() ;j++) pv->push_back(0.);
                for (int j=0;j< 2*NMODE ;j++) pi->push_back(0.);
        }
}

MC::~MC(){
//        cerr << "DESTROYING MC" <<endl;
        if (_NList.size() > 0 ){
                map <int, vector <int> * >::iterator itmap = _NList.begin();
                for ( ; itmap!= _NList.end(); itmap++){
                        itmap->second->clear();
                        delete itmap->second;
                }
                _NList.clear();
        }
        vector < vector<double> *  >::iterator ittime=_elstate.begin();
        for ( ; ittime != _elstate.end() ; ittime++){
                (*ittime)->clear();
                delete *ittime;
        }
        _elstate.clear();
        for (ittime = _vibnrg.begin() ; ittime != _vibnrg.end() ; ittime++){
                (*ittime)->clear();
                delete *ittime;
        }
        _vibnrg.clear();
        _CollectionTimes.clear();
        _RecombTimes.clear();

        _randomizeSeed = true;
        if (_randomizeSeed){
            srand48( unsigned(time(NULL)));
        }
}
void MC::Init(Parameters * para){
	vector <pair <int, int> >::iterator itTop = para->_top.begin();
	map < int, vector <int>* >::iterator itMap;
	for ( ; itTop != para->_top.end() ; ++itTop ){
	    int one = itTop->first;
	    itMap=_NList.find(one);
	    if (itMap == _NList.end() ){
		vector <int> * pV = new vector <int>;
		_NList.insert(make_pair(one, pV));
		pV->push_back(itTop->second);
	    }
	    else{
	    	itMap->second->push_back(itTop->second);
	    }
	}
        for ( vector<string>::iterator its = para->_properties.begin();
                its!=para->_properties.end(); ++its){
            _prop.push_back(*its);
        }
}

void MC::NextState(State *a){
    	int elfrom = a->ElState();
	int nufrom = a->GetNucl(elfrom);
	vector <int> elneigh = *(_NList[elfrom]); 
	double totRate = 0;
	vector <double> rates;
	vector <int>::iterator itn=elneigh.begin();
	for ( ; itn!=elneigh.end(); ++itn){
		int telto = *itn;
		int tnuto = a->GetNucl(telto);
//		if (elfrom == 0 && telto == 9 && nufrom == 8 && tnuto == 33 ){
//			cerr <<"DBG?" <<endl;
//		}
		double rate = _Rates->TotRate(elfrom, telto, nufrom, tnuto);
		rates.push_back(rate);
		totRate += rate;
	}

	//increment the total rate by the amount of vibrational decay we expect
	int totNNu= a->NuExcited();
	//int totNu = a->TotNu();
	double ratevib = _ratevib * double(totNNu);
	totRate += ratevib;
	/// increment time and bins
	IncrementTime(totRate, a);
	///mod state
	if ( drand48() >= ratevib/totRate ) {
	    //pick which neighbour it goes to
	    double ran = drand48()*(totRate-ratevib);
	    int elto;
	    vector <double>::iterator itR=rates.begin();
	    itn=elneigh.begin();
	    for ( ; itR != rates.end() ;++itR, itn++){
		if (*itR >= ran){
		    elto = *itn;
		    break;
		}
		else {
		    ran -= *itR;
		}
	    }

	    // pick which nuclear vibrational state it goes to
	    int nuto = a->GetNucl(elto);
	    pair <int, int> NuTo = _Rates->RandomDestination(elfrom, elto, nufrom, nuto);
	    //modify State
	    a->ChangePairNus(elfrom, elto, NuTo.first, NuTo.second);
	    a->SetEl(elto);
	}
	else {
	    a->DeExciteRandom();
	}
}

void MC::IncrementTime(const double & totRate, State * a){
	
        if (_prop[a->ElState()] == string("Collector") ){
            _CollectionTimes.push_back(_time);
            _run = false;
        }
        if (_prop[a->ElState()] == string("GroundState") ){
            _RecombTimes.push_back(_time);
            _run = false;
        }
        double dt= -log(drand48()) / totRate;

        double newt = _time+dt;
        if (_recordHistograms){
            int totvib= a->TotNu();
            int elfrom = a->ElState();

            //find the new tbin we are sitting in:
            int newtbin ;
            for (newtbin=_tbin; newtbin<_timebins.size() ;++newtbin){
                    if (_timebins[newtbin] > _time+dt)
                            break;
            }
            for (int i =_tbin; i<=newtbin ; ++i){
                    if (i > _timebins.size()-1) break;
                    double dtbin = min(_timebins[i], _time+dt) - _time;
                    _time+=dtbin;
                    (_elstate[i])->at(elfrom) += dtbin;
                    (_vibnrg[i])->at(totvib) += dtbin;
            }
            _tbin=newtbin;
        }
	_time=newt;
}

void MC::Run(State * a){
    _time=0;
    _run = true;
    _tbin=0;
    while ( _time < _tmax && _run == true) {
	NextState(a);
//	cout << *this;
    }
}

void MC::Run(const State & a, const int & n, string fout){
    cout <<"About to start" ;
    double avnufin=0.;
    _totRun = 0;
    const double toll = 1E-4;
    
    int imin=n;
    bool conv = false;
    double oldcollected = -1000;
    while (_totRun<n && conv== false){
	if (_totRun%1000==0) {
	    cout << '\r' << "Start run " << _totRun << " " << 100.0*double(_totRun)/double(n) <<"% done" ;
	    cout.flush();	
	}

	State init(a);
	Run(&init );
        _totRun++;
        avnufin += init.NuExcited();
        
        double newcollected = double (_CollectionTimes.size()) / _totRun ;
        double diff  = (oldcollected - newcollected) / newcollected;
        if (diff  < toll && _totRun > imin){
            conv=true;
        }
        oldcollected = newcollected;
    }
    ofstream out (fout.c_str() );
    out << *this;
    cout << '\n' << "So many charges collected " << PCollected() << endl;
    cout << '\n' << "So many charges recombined " << PRecomb() << endl;
    cout << '\n' << "So many phonons left: " << avnufin/double(_totRun) <<endl;
    
}
ostream& operator<< (ostream &out, MC &mc){
/*
    out << "#--Electronic State Occupation--" <<endl;
    out << "#timebin state prob" <<endl;
    vector <double>::iterator itt=mc._timebins.begin();
    vector <   vector<double> * >::iterator itel= mc._elstate.begin();
    for (; itt!= mc._timebins.end(); itel++, itt++){
	vector<double>::iterator itp = (*itel)->begin();
	double tottinbin=0.;
        for ( ; itp != (*itel)->end(); ++itp)tottinbin+=*itp;
	for (itp = (*itel)->begin() ; itp < (*itel)->end() ; ++itp){
        	out << *itt <<  '\t' << itp-(*itel)->begin() << '\t' << *itp/tottinbin<<'\n' ;
        } 	 
    }
    out << '\n' << '\n' ;
*/
    vector <double>::iterator ittim;
    out << "#Arrival Times" << endl;
    for ( ittim = mc._CollectionTimes.begin(); ittim!=mc._CollectionTimes.end(); ++ittim ){
        out << *ittim<<'\n';
    }
    out << endl;
/*
    out << "#Recombination times" <<endl;
    for ( ittim = mc._RecombTimes.begin(); ittim!=mc._RecombTimes.end(); ++ittim ){
        out << *ittim<<'\n';
    }
    out << endl;
*/

  /*  out << "--Vibrational State Occupation--" <<endl;
    out << "#timebin state prob" <<endl;
    itt=mc._timebins.begin();
    vector <   vector<double> * >::iterator itvib= mc._vibnrg.begin();
    for (; itt!= mc._timebins.end(); itvib++, itt++){
	vector<double>::iterator itp = (*itvib)->begin();
	double tottinbin=0.;
        for ( ; itp != (*itvib)->end(); ++itp)tottinbin+=*itp;
	for (itp = (*itvib)->begin() ; itp < (*itvib)->end() ; ++itp){
        	out << *itt <<  '\t' << itp-(*itvib)->begin() << '\t' << *itp/tottinbin<<'\n' ;
        } 	 
    }
    out << '\n' << '\n' ;
*/
}


