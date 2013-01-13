#include "Neighbours.h"

vector <Neighbours::_NeighType * > Neighbours::_AllNeighTypes;
map < pair <pair <int, int>, pair<int, int> >  , Neighbours::_Neigh*  > Neighbours::_AllRates;
double ** Neighbours::_FC;
map < pair<int, int> , double **> Neighbours::_PreComRates; 
map <  pair< int, int>  , vector<double> >  Neighbours::_AllProbs;


Neighbours::~Neighbours(){
//     cerr << "Clearing rates etc " <<endl;
    if (init == true) {
        _AllNeighTypes.clear();
        _AllRates.clear();
    }
    ClearFC();
    ClearPreComRates();
    init = false;
}

void Neighbours::InitAllProbs(){
    cout << "Start HASHing the probabilities " <<endl;
    pair <int, int> from;
    int to;
    for (from.first=0; from.first<2*NMODE+2; from.first++){
    	for (from.second=0; from.second<2*NMODE+2; from.second++){
	    _AllProbs.insert(make_pair(from, Prob(from.first, from.second) ) ) ;
	}
    }

}


void Neighbours::ClearFC(){
    if (initFC == true){
            delete []_FC[0];
            delete []_FC;
    }
    initFC=false;
}
void Neighbours::ClearPreComRates(){
    if (initPreComRates == true){
            map <pair<int, int> , double ** > ::iterator itMap = _PreComRates.begin();
           for ( ; itMap != _PreComRates.end(); ++itMap){
            double ** tmpAr = itMap-> second;
            delete []tmpAr[0];
            delete []tmpAr;
           }
    }
    initPreComRates=false;
}
void Neighbours::Init ( Parameters * para){
	vector <pair <int, int> >::iterator itTop = para->_top.begin();
	for ( ; itTop != para->_top.end() ; ++itTop ){
	   pair < int , int> nufrom;
	   // now a quad loop to the para._nmodemax for from and to
	   for ( nufrom.first=0; nufrom.first< NMODE ; nufrom.first++){
	       for ( nufrom.second=0; nufrom.second< NMODE; nufrom.second++){

		   _NeighType *From;
		   From = new _NeighType;

		   From->_ElFrom = itTop -> first;
		   From->_ElTo   = itTop -> second;
		   From->_NuFrom.first = nufrom.first;
		   From->_NuFrom.second = nufrom.second;
				
		   _AllNeighTypes.push_back(From);
	       }
	   } 
	}
	cout << "precompute the base rates" <<endl;
	
	//initialise the FC factors
	cout << "Compute FC factors" <<endl;
	ClearFC();
	InitFC(para,NMODE+1);
	//initialise the rates for symmetric modes
	cout << "Init rates" <<endl;
	ClearPreComRates();
	InitPreComRates(para);
	//init probs
	InitAllProbs();

        init=true;
} 

bool cmp (pair<double, pair <int,int> > A, pair <double, pair <int, int> > B ){
	return A.first > B.first;
}

void Neighbours::SortAndPrune(map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itRat){
    
    
    sort(itRat->second->_Rates.begin(), itRat->second->_Rates.end(), cmp);
    
    //truncate them for a fixed value
    const double tol = 1E-20;
    vector < pair <double, pair <int, int > > >::iterator itn = itRat->second->_Rates.begin();
    vector < pair <double, pair <int, int > > >::iterator itv = itRat->second->_Rates.end();
    for ( ; itn != itRat->second->_Rates.end() ; ++itn){
        if (itn->first < tol) break;
    }
    itRat->second->_Rates.erase(itn, itv);
    TotRate(itRat);
}

void Neighbours::TotRate(map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itRat) {
    // Compute the total Rate
    _Neigh * Dest = itRat->second;
    double totRate=0.;
    vector <pair <double, pair <int, int> > >::iterator itR = Dest->_Rates.begin();
    for ( ; itR!= Dest->_Rates.end() ;++itR ) {
        totRate += itR->first;
    }
    Dest -> _totRate = totRate;
    for ( itR=Dest->_Rates.begin() ; itR!= Dest->_Rates.end() ;++itR ) {
        itR->first /= totRate;
    }
}

// this probability gives you the chance of starting 
// in modes i and j and ending up in symmetric modes [l]
double Neighbours::Prob(const int & i, const int & j, const int & m){
	int max = i+j;
	double val=0;
	for(int k=0;k<=j;k++){
		for(int l=0;l <=i;l++){
                    double corr = KrDe(m,k+l)* double(binomial(i,l) * binomial(j,k) ) *
                            exp( (factln(l+k) + factln(max-l-k) - factln(i) - factln(j))/2. );
		    if (k%2==0) val += corr;
		    else val -= corr;
		}
	}
	return (val*val / pow(2.0, max) );

}

vector <double> Neighbours::Prob(const int & i, const int & j){
	int max = i+j;
	vector <double> res;
	for (int m=0; m<=min(max,NMODE);m++){
		res.push_back(Prob(i,j,m) );
	}
	//normalise(&res);
	return res;

}

double Neighbours::Prefactor(Parameters * para, double J){
	return J * J / hbar * sqrt(Pi / (para->_lambdaO * para->_kT) );
}

double Neighbours::Exponential(Parameters * para, 
		const int & NuFrom, const int & NuTo,
		const int & ElFrom, const int & ElTo){
	double DE =  ( para-> _energies)[ElTo]  -( para-> _energies)[ElFrom] ;
	return _FC[NuFrom][NuTo]*_FC[NuFrom][NuTo]*exp(- 
		       (DE + para->_lambdaO+(NuTo-NuFrom)*para->_hbarOmega)*
		       (DE + para->_lambdaO+(NuTo-NuFrom)*para->_hbarOmega)
		       /(4 * para->_lambdaO * para->_kT));	
}

void Neighbours::InitPreComRates( Parameters * para){
	if (initPreComRates == false){
		

		vector < pair <int, int> >::iterator itTop = para->_top.begin();
		vector <double>::iterator itTrans = para->_transfers.begin();

		for ( ; itTop != para->_top.end();itTop++, itTrans++){
			pair <int, int> El = (*itTop);
//			if (El.first == 2 ){
//				cerr << "WTF???" <<endl;
//			}
			double pre  = Prefactor(para, *itTrans);
			//assign memory to this biatch
			double ** vals;
			vals = new double * [NMODE+1];
			vals[0] = new double [(NMODE+1)*(NMODE+1)];
			for (int i =1;i <NMODE+1 ; i++) vals[i] = vals[i-1] +NMODE+1;	
			for (int i=0; i<NMODE+1;i++){
				for(int j=0;j<NMODE+1;j++){
                                  /*  if (El.second == para->_energies.size() -1  ){ // HACK!
                                        vals[i][j] =1E13/NMODE;
                                    }
                                    else{*/
					vals[i][j] = pre * Exponential(para, i,j,El.first, El.second );
                                    //}
				}
			}
			_PreComRates.insert(make_pair(El, vals));
		}
		initPreComRates=true;
	}
}

/*
void Neighbours::ComputeAllRates(Parameters * para ) {
	
	cout << "HASH it" <<endl;
	vector <_NeighType*>::iterator itNeigh = _AllNeighTypes.begin();
	int ntyp = _AllNeighTypes.size();
	int c=0;
	for ( ; itNeigh < _AllNeighTypes.end(); ++itNeigh,++c ){
	    	if (c%100==0){
			cout << '\r' << " HASH " << int(100.* double(c)/double(ntyp)) 
			    << "% complete               " ;
			cout.flush();
		}
		int elfrom=(*itNeigh)->_ElFrom;
//		if (elfrom == 2 ){
//			cout << "WTF???" <<endl;
//		}
		int elto=(*itNeigh)->_ElTo;
		int nufrom = (*itNeigh)->_NuFrom.first;
		int nuto = (*itNeigh)->_NuFrom.second;
		_Neigh * NDest;
		NDest = new _Neigh;


		vector <double> prob = _AllProbs[(*itNeigh)->_NuFrom];

//		if (nufrom == 8 && nuto== 0){
//		    cerr <<"DBG"<<endl;
//		}
		// now we need to double loop over the nuclear modes the sys is going to
		for (int i=0;i<NMODE;i++){
			for (int j=0;j<NMODE;j++){
//				if (i==3 && j==10 && nufrom ==0 && nuto==0 &&elfrom==1 &&elto==2){
//					cout <<"BEG ONLY"<<endl;
//				}
				pair <int , int>  To;
				To.first=i;
				To.second=j;
				double Rate=0.;

				//cycle over start symmetric mode
				for ( int fromSymm=0; fromSymm <  prob.size();fromSymm++){
					//end symm mode must by
				        // toSymm + fromAsymm = i+j
				    	// toSymm + nufrom +nuto - fromSymm = i+j
					int toSymm = j+i+fromSymm-nufrom-nuto;
			
					// the is an overflow here when fromSymm or toSymm>Nmode
					if (toSymm < NMODE && toSymm >= 0){
					    	double prob1 = prob[fromSymm];
						double barerate = _PreComRates[make_pair(elfrom,elto)][fromSymm][toSymm];
						double prob2 = _AllProbs[make_pair(nufrom+nuto-fromSymm, toSymm)][i];
						Rate += prob1 * prob2 *barerate;
				//		if (toSymm == 15 && fromSymm == 3 && nuto == 0 &&  nufrom == 3 ){
			//			    cerr << "ETF> ???" <<endl;
		//				}
					}
							
				}

				//create the bleeding Neigh type and load it onto the static array
				NDest->_Rates.push_back(make_pair(Rate, To));
//				cerr << i << " " << j << " " << nufrom << " " << nuto << " " << endl;
//				if ( i== 10 && j == 7 && nufrom == 1 && nuto == 4) {
//				    cerr << "FUCK!!!" <<endl;
//				}
			}
		}

		pair< pair <int, int>, pair<int, int> > pairrate = make_pair ( make_pair(elfrom , elto ), make_pair(nufrom, nuto));

		//Add the map entry
		_AllRates.insert(make_pair(pairrate, NDest));
	}
	cout <<endl;
	cout <<"Start Sorting the HASHed neighbour list" <<endl;
	SortAndPrune();
}
*/
//debug routines 
void Neighbours::PrintAllRates( ostream & out ){
	map < pair<pair<int, int> , pair<int, int> >, _Neigh*>::iterator itR = _AllRates.begin();
	for (; itR!=_AllRates.end();++itR){
		pair < pair<int, int>, pair<int,int> > NT = (itR->first);
		out << "El state from: " << NT.first.first << '\t' <<
		"El state to: " << NT.first.second << '\t' <<
		"Nuclear modes from " <<NT.second.first <<" " << NT.second.second <<'\n';	

		_Neigh N = *(itR->second);
		out << "Total Rate : " << N._totRate << endl;
		vector < pair < double, pair <int, int> > > rates = (N._Rates);
		for (int i=0;i < rates.size() ; i++){
			out << "	Rate to:" << rates[i].first << " state to " << rates[i].second.first 
				<< ", " << rates[i].second.second <<endl;
		}
	}
}

void Neighbours::PrintAllNeighTypes( ostream & out ){
	vector < _NeighType*>::iterator itNT = _AllNeighTypes.begin();
	for (; itNT!=_AllNeighTypes.end();++itNT){
		out << (*itNT) << '\t' << (*itNT)->_ElFrom	
			    << '\t' << (*itNT)->_ElTo
			    << '\t' << (*itNT)->_NuFrom.first
			    << '\t' << (*itNT)->_NuFrom.second << '\n';
	}
}
void Neighbours::PrintProbs( ostream & out ){
        static map < pair< int, int> ,  vector<double> >::iterator itP= _AllProbs.begin();
        for ( ; itP != _AllProbs.end() ; ++itP){
            vector <double>::iterator itprob = itP->second.begin();
            out << "From nucl " << itP->first.first << " and " << itP->first.second << " ";
            for ( ; itprob != itP->second.end(); ++itprob){
                out << *itprob << " ";
            }
            out << '\n';
        }
        out << endl;
}

void Neighbours::PrintCalcRates(ostream & out) {
	vector < _NeighType*>::iterator itNT = _AllNeighTypes.begin();
        for ( ; itNT != _AllNeighTypes.end() ; itNT++ ){
            int elfrom = (*itNT)->_ElFrom;
            int elto = (*itNT) -> _ElTo;
            int fromSymm = (*itNT)->_NuFrom.first;
            int toSymm   = (*itNT)->_NuFrom.second;
            double barerate = _PreComRates[make_pair(elfrom,elto)][fromSymm][toSymm];
            out << "From elstate " << elfrom << " to elstate " << elto << " with symmetric from: "
                    << fromSymm <<  " to symmetric mode: " << toSymm << " rate: " << barerate << '\n';
         
        }
        out << endl;
        
}



// The important accesor routine is to get the totalrate for a give elfrom, elto, nufrom and nuto.
const double & Neighbours::TotRate(const int & elfrom, const int & elto, const int & nufrom, const int & nuto){
    map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itmap=GetItRate(elfrom, elto, nufrom, nuto);
    return (*itmap).second->_totRate;

}
const pair<int, int> &  Neighbours::RandomDestination(const int & elfrom, const int & elto, const int & nufrom, const int & nuto){
    map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itmap=GetItRate(elfrom, elto, nufrom, nuto);

    double ran = drand48();
    vector < pair < double, pair <int, int > > >::iterator itR=itmap->second->_Rates.begin();
    for ( ; itR != itmap->second->_Rates.end(); itR++){
        if (itR->first >= ran) return itR->second;
        else ran -= itR->first;
    }
    cerr << "Error in finding a random destination" <<endl;

}

map < pair <pair <int, int>, pair<int, int> > , Neighbours::_Neigh*  >::iterator Neighbours::GetItRate(const int & elfrom, const int & elto, const int & nufrom, const int & nuto){
    pair< pair <int, int>, pair<int, int> > pairrate = make_pair ( make_pair(elfrom , elto ), make_pair(nufrom, nuto));

    map < pair <pair <int, int>, pair<int, int> > , _Neigh*  >::iterator itmap=_AllRates.find(pairrate) ;
    if (itmap == _AllRates.end() ) {        
        _Neigh * NDest;
        NDest = new _Neigh;

        vector <double> prob = _AllProbs[make_pair(nufrom, nuto)];
        for (int i=0;i<NMODE;i++){
            for (int j=0;j<NMODE;j++){
                pair <int , int>  To;
                To.first=i;
                To.second=j;
                double Rate=0.;

                //cycle over start symmetric mode
                for ( int fromSymm=0; fromSymm <  prob.size();fromSymm++){
                        //end symm mode must by
                        //toSymm + fromAsymm = i+j
                        //toSymm  = i+j - fromAsymm
                        //start modes must obey:
                        // fromSymm + fromAsymm = nufrom +nuto
                        // fromAsymm = nufrom +nuto - fromSymm
                        // so:
                        // toSymm  = i+j - (nufrom+nuto-fromsymm)
                        int toSymm = j+i+fromSymm-nufrom-nuto;

                        // the is an overflow here when fromSymm or toSymm>Nmode
                        if (toSymm < NMODE && toSymm >= 0){
                            double prob1 = prob[fromSymm];
                            double barerate = _PreComRates[make_pair(elfrom,elto)][fromSymm][toSymm];
                            double prob2 = _AllProbs[make_pair(nufrom+nuto-fromSymm, toSymm)][i];
                            Rate += prob1 * prob2 *barerate;
                        }
                }

                //create the bleeding Neigh type and load it onto the static array
                NDest->_Rates.push_back(make_pair(Rate, To));
            }
        }

       
        //Add the map entry
        _AllRates.insert(make_pair(pairrate, NDest));
        //do total rate and pruning
        itmap=_AllRates.find(pairrate);
        SortAndPrune(itmap);
    }
    return itmap;
}

void Neighbours::InitFC( Parameters * para, const int & n ){
    const double e = 2.718281828459045;
	if (initFC == false ){
		initFC=true;
		double S = para->_S;
		_FC = new double * [n];
		_FC[0] = new double [n*n];
		for (int i=1;i<n;i++) _FC[i] = _FC[i-1]+n;

		//double TMPFC[13][13] = {{pow(2.718281828459045,-0.5*S),(-1.*sqrt(S))/pow(2.718281828459045,0.5*S),(0.7071067811865476*S)/pow(2.718281828459045,0.5*S),(-0.408248290463863*pow(S,1.5))/pow(2.718281828459045,0.5*S),(0.2041241452319315*pow(S,2))/pow(2.718281828459045,0.5*S),(-0.09128709291752768*pow(S,2.5))/pow(2.718281828459045,0.5*S),(0.03726779962499649*pow(S,3))/pow(2.718281828459045,0.5*S),(-0.014085904245475275*pow(S,3.5))/pow(2.718281828459045,0.5*S),(0.004980119205559973*pow(S,4))/pow(2.718281828459045,0.5*S),(-0.0016600397351866577*pow(S,4.5))/pow(2.718281828459045,0.5*S),(0.0005249506569572601*pow(S,5))/pow(2.718281828459045,0.5*S),(-0.00015827857841616382*pow(S,5.5))/pow(2.718281828459045,0.5*S),(0.00004569108992776174*pow(S,6))/pow(2.718281828459045,0.5*S)},{(-1.*sqrt(S))/pow(2.718281828459045,0.5*S),(-1.*(-1. + S))/pow(2.718281828459045,0.5*S),(0.7071067811865476*(-2. + S)*sqrt(S))/pow(2.718281828459045,0.5*S),(-0.408248290463863*(-3. + S)*S)/pow(2.718281828459045,0.5*S),(0.2041241452319315*(-4. + S)*pow(S,1.5))/pow(2.718281828459045,0.5*S),(-0.09128709291752768*(-5. + S)*pow(S,2))/pow(2.718281828459045,0.5*S),(0.03726779962499649*(-6. + S)*pow(S,2.5))/pow(2.718281828459045,0.5*S),(-0.014085904245475275*(-7. + S)*pow(S,3))/pow(2.718281828459045,0.5*S),(0.004980119205559973*(-8. + S)*pow(S,3.5))/pow(2.718281828459045,0.5*S),(-0.0016600397351866577*(-9. + S)*pow(S,4))/pow(2.718281828459045,0.5*S),(0.0005249506569572601*(-10. + S)*pow(S,4.5))/pow(2.718281828459045,0.5*S),(-0.00015827857841616382*(-11. + S)*pow(S,5))/pow(2.718281828459045,0.5*S),(0.00004569108992776174*(-12. + S)*pow(S,5.5))/pow(2.718281828459045,0.5*S)},{(0.7071067811865476*S)/pow(2.718281828459045,0.5*S),(0.7071067811865476*(-2. + S)*sqrt(S))/pow(2.718281828459045,0.5*S),(0.5*(2. - 4.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.28867513459481287*sqrt(S)*(6. - 6.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.14433756729740643*S*(12. - 8.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.06454972243679029*pow(S,1.5)*(20. - 10.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.026352313834736494*pow(S,2)*(30. - 12.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.009960238411119947*pow(S,2.5)*(42. - 14.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.003521476061368819*pow(S,3)*(56. - 16.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.0011738253537896062*pow(S,3.5)*(72. - 18.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.0003711961693228117*pow(S,4)*(90. - 20.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.00011191985611463616*pow(S,4.5)*(110. - 22.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.00003230847952772469*pow(S,5)*(132. - 24.*S + pow(S,2)))/pow(2.718281828459045,0.5*S)},{(-0.408248290463863*pow(S,1.5))/pow(2.718281828459045,0.5*S),(-0.408248290463863*(-3. + S)*S)/pow(2.718281828459045,0.5*S),(-0.28867513459481287*sqrt(S)*(6. - 6.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.16666666666666666*(6. - 18.*S + 9.*pow(S,2) - 1.*pow(S,3)))/pow(2.718281828459045,0.5*S),(0.08333333333333333*sqrt(S)*(-24. + pow(-6. + S,2)*S))/pow(2.718281828459045,0.5*S),(-0.03726779962499649*S*(-60. + 60.*S - 15.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.015214515486254613*pow(S,1.5)*(-120. + 90.*S - 18.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.005750546327852952*pow(S,2)*(-210. + 126.*S - 21.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.0020331251519761107*pow(S,2.5)*(-336. + 168.*S - 24.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.0006777083839920369*pow(S,3)*(-504. + 216.*S - 27.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.00021431020828068323*pow(S,3.5)*(-720. + 270.*S - 30.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.00006461695905544937*pow(S,4)*(-990. + 330.*S - 33.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.00001865330935243936*pow(S,4.5)*(-1320. + 396.*S - 36.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S)},{(0.2041241452319315*pow(S,2))/pow(2.718281828459045,0.5*S),(0.2041241452319315*(-4. + S)*pow(S,1.5))/pow(2.718281828459045,0.5*S),(0.14433756729740643*S*(12. - 8.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.08333333333333333*sqrt(S)*(-24. + pow(-6. + S,2)*S))/pow(2.718281828459045,0.5*S),(0.041666666666666664*(24. - 96.*S + 72.*pow(S,2) - 16.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(-0.018633899812498245*sqrt(S)*(120. - 240.*S + 120.*pow(S,2) - 20.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(0.007607257743127306*S*(360. - 480.*S + 180.*pow(S,2) - 24.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(-0.002875273163926476*pow(S,1.5)*(840. - 840.*S + 252.*pow(S,2) - 28.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(0.0010165625759880554*pow(S,2)*(1680. - 1344.*S + 336.*pow(S,2) - 32.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(-0.00033885419199601845*(-6. + S)*pow(S,2.5)*(-504. + 252.*S - 30.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.00010715510414034161*pow(S,3)*(5040. - 2880.*S + 540.*pow(S,2) - 40.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(-0.00003230847952772469*pow(S,3.5)*(7920. - 3960.*S + 660.*pow(S,2) - 44.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(9.32665467621968e-6*pow(S,4)*(11880. - 5280.*S + 792.*pow(S,2) - 48.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S)},{(-0.09128709291752768*pow(S,2.5))/pow(2.718281828459045,0.5*S),(-0.09128709291752768*(-5. + S)*pow(S,2))/pow(2.718281828459045,0.5*S),(-0.06454972243679029*pow(S,1.5)*(20. - 10.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.03726779962499649*S*(-60. + 60.*S - 15.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.018633899812498245*sqrt(S)*(120. - 240.*S + 120.*pow(S,2) - 20.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(0.008333333333333333*(120. - 1.*S*(30. - 15.*S + pow(S,2))*(20. - 10.*S + pow(S,2))))/pow(2.718281828459045,0.5*S),(0.0034020690871988586*sqrt(S)*(-720. + 1800.*S - 1200.*pow(S,2) + 300.*pow(S,3) - 30.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(-0.0012858612496840993*S*(-2520. + 4200.*S - 2100.*pow(S,2) + 420.*pow(S,3) - 35.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(0.0004546206046583175*pow(S,1.5)*(-6720. + 8400.*S - 3360.*pow(S,2) + 560.*pow(S,3) - 40.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(-0.0001515402015527725*pow(S,2)*(-15120. + 15120.*S - 5040.*pow(S,2) + 720.*pow(S,3) - 45.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(0.000047921219398774604*pow(S,2.5)*(-30240. + 25200.*S - 7200.*pow(S,2) + 900.*pow(S,3) - 50.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(-0.000014448791294730538*pow(S,3)*(-55440. + 39600.*S - 9900.*pow(S,2) + 1100.*pow(S,3) - 55.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(4.1710067717387e-6*pow(S,3.5)*(-95040. + 59400.*S - 13200.*pow(S,2) + 1320.*pow(S,3) - 60.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S)},{(0.03726779962499649*pow(S,3))/pow(2.718281828459045,0.5*S),(0.03726779962499649*(-6. + S)*pow(S,2.5))/pow(2.718281828459045,0.5*S),(0.026352313834736494*pow(S,2)*(30. - 12.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.015214515486254613*pow(S,1.5)*(-120. + 90.*S - 18.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.007607257743127306*S*(360. - 480.*S + 180.*pow(S,2) - 24.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(0.0034020690871988586*sqrt(S)*(-720. + 1800.*S - 1200.*pow(S,2) + 300.*pow(S,3) - 30.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(0.001388888888888889*(720. + (-6. + S)*S*(720. - 780.*S + 270.*pow(S,2) - 30.*pow(S,3) + pow(S,4))))/pow(2.718281828459045,0.5*S),(-0.0005249506569572601*sqrt(S)*(5040. - 15120.*S + 12600.*pow(S,2) - 4200.*pow(S,3) + 630.*pow(S,4) - 42.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(0.00018559808466140585*S*(20160. - 40320.*S + 25200.*pow(S,2) - 6720.*pow(S,3) + 840.*pow(S,4) - 48.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(-0.00006186602822046861*pow(S,1.5)*(60480. - 90720.*S + 45360.*pow(S,2) - 10080.*pow(S,3) + 1080.*pow(S,4) - 54.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(0.000019563755896493438*pow(S,2)*(151200. - 181440.*S + 75600.*pow(S,2) - 14400.*pow(S,3) + 1350.*pow(S,4) - 60.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(-5.898694345342888e-6*pow(S,2.5)*(332640. - 332640.*S + 118800.*pow(S,2) - 19800.*pow(S,3) + 1650.*pow(S,4) - 66.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(1.70280638407552e-6*pow(S,3)*(665280. - 570240.*S + 178200.*pow(S,2) - 26400.*pow(S,3) + 1980.*pow(S,4) - 72.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S)},{(-0.014085904245475275*pow(S,3.5))/pow(2.718281828459045,0.5*S),(-0.014085904245475275*(-7. + S)*pow(S,3))/pow(2.718281828459045,0.5*S),(-0.009960238411119947*pow(S,2.5)*(42. - 14.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.005750546327852952*pow(S,2)*(-210. + 126.*S - 21.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.002875273163926476*pow(S,1.5)*(840. - 840.*S + 252.*pow(S,2) - 28.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(-0.0012858612496840993*S*(-2520. + 4200.*S - 2100.*pow(S,2) + 420.*pow(S,3) - 35.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(-0.0005249506569572601*sqrt(S)*(5040. - 15120.*S + 12600.*pow(S,2) - 4200.*pow(S,3) + 630.*pow(S,4) - 42.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(0.0001984126984126984*(5040. - 35280.*S + 52920.*pow(S,2) - 29400.*pow(S,3) + 7350.*pow(S,4) - 882.*pow(S,5) + 49.*pow(S,6) - 1.*pow(S,7)))/pow(2.718281828459045,0.5*S),(0.00007014948226057019*sqrt(S)*(-40320. + 141120.*S - 141120.*pow(S,2) + 58800.*pow(S,3) - 11760.*pow(S,4) + 1176.*pow(S,5) - 56.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(-0.0000233831607535234*S*(-181440. + 423360.*S - 317520.*pow(S,2) + 105840.*pow(S,3) - 17640.*pow(S,4) + 1512.*pow(S,5) - 63.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(7.394404687499305e-6*pow(S,1.5)*(-604800. + 1.0584e6*S - 635040.*pow(S,2) + 176400.*pow(S,3) - 25200.*pow(S,4) + 1890.*pow(S,5) - 70.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(-2.2294968996800335e-6*pow(S,2)*(-1.6632e6 + 2.32848e6*S - 1.16424e6*pow(S,2) + 277200.*pow(S,3) - 34650.*pow(S,4) + 2310.*pow(S,5) - 77.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(6.436003175938517e-7*pow(S,2.5)*(-3.99168e6 + 4.65696e6*S - 1.99584e6*pow(S,2) + 415800.*pow(S,3) - 46200.*pow(S,4) + 2772.*pow(S,5) - 84.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S)},{(0.004980119205559973*pow(S,4))/pow(2.718281828459045,0.5*S),(0.004980119205559973*(-8. + S)*pow(S,3.5))/pow(2.718281828459045,0.5*S),(0.003521476061368819*pow(S,3)*(56. - 16.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.0020331251519761107*pow(S,2.5)*(-336. + 168.*S - 24.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.0010165625759880554*pow(S,2)*(1680. - 1344.*S + 336.*pow(S,2) - 32.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(0.0004546206046583175*pow(S,1.5)*(-6720. + 8400.*S - 3360.*pow(S,2) + 560.*pow(S,3) - 40.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(0.00018559808466140585*S*(20160. - 40320.*S + 25200.*pow(S,2) - 6720.*pow(S,3) + 840.*pow(S,4) - 48.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(0.00007014948226057019*sqrt(S)*(-40320. + 141120.*S - 141120.*pow(S,2) + 58800.*pow(S,3) - 11760.*pow(S,4) + 1176.*pow(S,5) - 56.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(0.0000248015873015873*(40320. - 322560.*S + 564480.*pow(S,2) - 376320.*pow(S,3) + 117600.*pow(S,4) - 18816.*pow(S,5) + 1568.*pow(S,6) - 64.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(-8.267195767195768e-6*sqrt(S)*(362880. - 1.45152e6*S + 1.69344e6*pow(S,2) - 846720.*pow(S,3) + 211680.*pow(S,4) - 28224.*pow(S,5) + 2016.*pow(S,6) - 72.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(2.6143168486841763e-6*S*(1.8144e6 - 4.8384e6*S + 4.2336e6*pow(S,2) - 1.69344e6*pow(S,3) + 352800.*pow(S,4) - 40320.*pow(S,5) + 2520.*pow(S,6) - 80.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(-7.882461881990679e-7*pow(S,1.5)*(6.6528e6 - 1.33056e7*S + 9.31392e6*pow(S,2) - 3.10464e6*pow(S,3) + 554400.*pow(S,4) - 55440.*pow(S,5) + 3080.*pow(S,6) - 88.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(2.275470744722141e-7*pow(S,2)*(1.99584e7 - 3.193344e7*S + 1.862784e7*pow(S,2) - 5.32224e6*pow(S,3) + 831600.*pow(S,4) - 73920.*pow(S,5) + 3696.*pow(S,6) - 96.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S)},{(-0.0016600397351866577*pow(S,4.5))/pow(2.718281828459045,0.5*S),(-0.0016600397351866577*(-9. + S)*pow(S,4))/pow(2.718281828459045,0.5*S),(-0.0011738253537896062*pow(S,3.5)*(72. - 18.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.0006777083839920369*pow(S,3)*(-504. + 216.*S - 27.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.00033885419199601845*(-6. + S)*pow(S,2.5)*(-504. + 252.*S - 30.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.0001515402015527725*pow(S,2)*(-15120. + 15120.*S - 5040.*pow(S,2) + 720.*pow(S,3) - 45.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(-0.00006186602822046861*pow(S,1.5)*(60480. - 90720.*S + 45360.*pow(S,2) - 10080.*pow(S,3) + 1080.*pow(S,4) - 54.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(-0.0000233831607535234*S*(-181440. + 423360.*S - 317520.*pow(S,2) + 105840.*pow(S,3) - 17640.*pow(S,4) + 1512.*pow(S,5) - 63.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(-8.267195767195768e-6*sqrt(S)*(362880. - 1.45152e6*S + 1.69344e6*pow(S,2) - 846720.*pow(S,3) + 211680.*pow(S,4) - 28224.*pow(S,5) + 2016.*pow(S,6) - 72.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(2.7557319223985893e-6*(362880. - 1.*(-6. + S)*S*(-544320. + 997920.*S - 680400.*pow(S,2) + 204120.*pow(S,3) - 29484.*pow(S,4) + 2142.*pow(S,5) - 75.*pow(S,6) + pow(S,7))))/pow(2.718281828459045,0.5*S),(8.714389495613921e-7*sqrt(S)*(-3.6288e6 + 1.63296e7*S - 2.17728e7*pow(S,2) + 1.27008e7*pow(S,3) - 3.81024e6*pow(S,4) + 635040.*pow(S,5) - 60480.*pow(S,6) + 3240.*pow(S,7) - 90.*pow(S,8) + pow(S,9)))/pow(2.718281828459045,0.5*S),(-2.6274872939968927e-7*S*(-1.99584e7 + 5.98752e7*S - 5.98752e7*pow(S,2) + 2.794176e7*pow(S,3) - 6.98544e6*pow(S,4) + 997920.*pow(S,5) - 83160.*pow(S,6) + 3960.*pow(S,7) - 99.*pow(S,8) + pow(S,9)))/pow(2.718281828459045,0.5*S),(7.584902482407137e-8*pow(S,1.5)*(-7.98336e7 + 1.796256e8*S - 1.4370048e8*pow(S,2) + 5.588352e7*pow(S,3) - 1.197504e7*pow(S,4) + 1.49688e6*pow(S,5) - 110880.*pow(S,6) + 4752.*pow(S,7) - 108.*pow(S,8) + pow(S,9)))/pow(2.718281828459045,0.5*S)},{(0.0005249506569572601*pow(S,5))/pow(2.718281828459045,0.5*S),(0.0005249506569572601*(-10. + S)*pow(S,4.5))/pow(2.718281828459045,0.5*S),(0.0003711961693228117*pow(S,4)*(90. - 20.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.00021431020828068323*pow(S,3.5)*(-720. + 270.*S - 30.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(0.00010715510414034161*pow(S,3)*(5040. - 2880.*S + 540.*pow(S,2) - 40.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(0.000047921219398774604*pow(S,2.5)*(-30240. + 25200.*S - 7200.*pow(S,2) + 900.*pow(S,3) - 50.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(0.000019563755896493438*pow(S,2)*(151200. - 181440.*S + 75600.*pow(S,2) - 14400.*pow(S,3) + 1350.*pow(S,4) - 60.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(7.394404687499305e-6*pow(S,1.5)*(-604800. + 1.0584e6*S - 635040.*pow(S,2) + 176400.*pow(S,3) - 25200.*pow(S,4) + 1890.*pow(S,5) - 70.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(2.6143168486841763e-6*S*(1.8144e6 - 4.8384e6*S + 4.2336e6*pow(S,2) - 1.69344e6*pow(S,3) + 352800.*pow(S,4) - 40320.*pow(S,5) + 2520.*pow(S,6) - 80.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(8.714389495613921e-7*sqrt(S)*(-3.6288e6 + 1.63296e7*S - 2.17728e7*pow(S,2) + 1.27008e7*pow(S,3) - 3.81024e6*pow(S,4) + 635040.*pow(S,5) - 60480.*pow(S,6) + 3240.*pow(S,7) - 90.*pow(S,8) + pow(S,9)))/pow(2.718281828459045,0.5*S),(2.755731922398589e-7*(3.6288e6 - 3.6288e7*S + 8.1648e7*pow(S,2) - 7.2576e7*pow(S,3) + 3.1752e7*pow(S,4) - 7.62048e6*pow(S,5) + 1.0584e6*pow(S,6) - 86400.*pow(S,7) + 4050.*pow(S,8) - 100.*pow(S,9) + pow(S,10)))/pow(2.718281828459045,0.5*S),(-8.308844372182639e-8*sqrt(S)*(3.99168e7 - 1.99584e8*S + 2.99376e8*pow(S,2) - 1.99584e8*pow(S,3) + 6.98544e7*pow(S,4) - 1.397088e7*pow(S,5) + 1.6632e6*pow(S,6) - 118800.*pow(S,7) + 4950.*pow(S,8) - 110.*pow(S,9) + pow(S,10)))/pow(2.718281828459045,0.5*S),(2.398556767467177e-8*S*(2.395008e8 - 7.98336e8*S + 8.98128e8*pow(S,2) - 4.790016e8*pow(S,3) + 1.397088e8*pow(S,4) - 2.395008e7*pow(S,5) + 2.4948e6*pow(S,6) - 158400.*pow(S,7) + 5940.*pow(S,8) - 120.*pow(S,9) + pow(S,10)))/pow(2.718281828459045,0.5*S)},{(-0.00015827857841616382*pow(S,5.5))/pow(2.718281828459045,0.5*S),(-0.00015827857841616382*(-11. + S)*pow(S,5))/pow(2.718281828459045,0.5*S),(-0.00011191985611463616*pow(S,4.5)*(110. - 22.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(-0.00006461695905544937*pow(S,4)*(-990. + 330.*S - 33.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(-0.00003230847952772469*pow(S,3.5)*(7920. - 3960.*S + 660.*pow(S,2) - 44.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(-0.000014448791294730538*pow(S,3)*(-55440. + 39600.*S - 9900.*pow(S,2) + 1100.*pow(S,3) - 55.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(-5.898694345342888e-6*pow(S,2.5)*(332640. - 332640.*S + 118800.*pow(S,2) - 19800.*pow(S,3) + 1650.*pow(S,4) - 66.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(-2.2294968996800335e-6*pow(S,2)*(-1.6632e6 + 2.32848e6*S - 1.16424e6*pow(S,2) + 277200.*pow(S,3) - 34650.*pow(S,4) + 2310.*pow(S,5) - 77.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(-7.882461881990679e-7*pow(S,1.5)*(6.6528e6 - 1.33056e7*S + 9.31392e6*pow(S,2) - 3.10464e6*pow(S,3) + 554400.*pow(S,4) - 55440.*pow(S,5) + 3080.*pow(S,6) - 88.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(-2.6274872939968927e-7*S*(-1.99584e7 + 5.98752e7*S - 5.98752e7*pow(S,2) + 2.794176e7*pow(S,3) - 6.98544e6*pow(S,4) + 997920.*pow(S,5) - 83160.*pow(S,6) + 3960.*pow(S,7) - 99.*pow(S,8) + pow(S,9)))/pow(2.718281828459045,0.5*S),(-8.308844372182639e-8*sqrt(S)*(3.99168e7 - 1.99584e8*S + 2.99376e8*pow(S,2) - 1.99584e8*pow(S,3) + 6.98544e7*pow(S,4) - 1.397088e7*pow(S,5) + 1.6632e6*pow(S,6) - 118800.*pow(S,7) + 4950.*pow(S,8) - 110.*pow(S,9) + pow(S,10)))/pow(2.718281828459045,0.5*S),(2.505210838544172e-8*(3.99168e7 - 4.390848e8*S + 1.097712e9*pow(S,2) - 1.097712e9*pow(S,3) + 5.48856e8*pow(S,4) - 1.5367968e8*pow(S,5) + 2.561328e7*pow(S,6) - 2.6136e6*pow(S,7) + 163350.*pow(S,8) - 6050.*pow(S,9) + 121.*pow(S,10) - 1.*pow(S,11)))/pow(2.718281828459045,0.5*S),(7.231920760051229e-9*sqrt(S)*(-4.790016e8 + 2.6345088e9*S - 4.390848e9*pow(S,2) + 3.293136e9*pow(S,3) - 1.3172544e9*pow(S,4) + 3.0735936e8*pow(S,5) - 4.390848e7*pow(S,6) + 3.9204e6*pow(S,7) - 217800.*pow(S,8) + 7260.*pow(S,9) - 132.*pow(S,10) + pow(S,11)))/pow(2.718281828459045,0.5*S)},{(0.00004569108992776174*pow(S,6))/pow(2.718281828459045,0.5*S),(0.00004569108992776174*(-12. + S)*pow(S,5.5))/pow(2.718281828459045,0.5*S),(0.00003230847952772469*pow(S,5)*(132. - 24.*S + pow(S,2)))/pow(2.718281828459045,0.5*S),(0.00001865330935243936*pow(S,4.5)*(-1320. + 396.*S - 36.*pow(S,2) + pow(S,3)))/pow(2.718281828459045,0.5*S),(9.32665467621968e-6*pow(S,4)*(11880. - 5280.*S + 792.*pow(S,2) - 48.*pow(S,3) + pow(S,4)))/pow(2.718281828459045,0.5*S),(4.1710067717387e-6*pow(S,3.5)*(-95040. + 59400.*S - 13200.*pow(S,2) + 1320.*pow(S,3) - 60.*pow(S,4) + pow(S,5)))/pow(2.718281828459045,0.5*S),(1.70280638407552e-6*pow(S,3)*(665280. - 570240.*S + 178200.*pow(S,2) - 26400.*pow(S,3) + 1980.*pow(S,4) - 72.*pow(S,5) + pow(S,6)))/pow(2.718281828459045,0.5*S),(6.436003175938517e-7*pow(S,2.5)*(-3.99168e6 + 4.65696e6*S - 1.99584e6*pow(S,2) + 415800.*pow(S,3) - 46200.*pow(S,4) + 2772.*pow(S,5) - 84.*pow(S,6) + pow(S,7)))/pow(2.718281828459045,0.5*S),(2.275470744722141e-7*pow(S,2)*(1.99584e7 - 3.193344e7*S + 1.862784e7*pow(S,2) - 5.32224e6*pow(S,3) + 831600.*pow(S,4) - 73920.*pow(S,5) + 3696.*pow(S,6) - 96.*pow(S,7) + pow(S,8)))/pow(2.718281828459045,0.5*S),(7.584902482407137e-8*pow(S,1.5)*(-7.98336e7 + 1.796256e8*S - 1.4370048e8*pow(S,2) + 5.588352e7*pow(S,3) - 1.197504e7*pow(S,4) + 1.49688e6*pow(S,5) - 110880.*pow(S,6) + 4752.*pow(S,7) - 108.*pow(S,8) + pow(S,9)))/pow(2.718281828459045,0.5*S),(2.398556767467177e-8*S*(2.395008e8 - 7.98336e8*S + 8.98128e8*pow(S,2) - 4.790016e8*pow(S,3) + 1.397088e8*pow(S,4) - 2.395008e7*pow(S,5) + 2.4948e6*pow(S,6) - 158400.*pow(S,7) + 5940.*pow(S,8) - 120.*pow(S,9) + pow(S,10)))/pow(2.718281828459045,0.5*S),(7.231920760051229e-9*sqrt(S)*(-4.790016e8 + 2.6345088e9*S - 4.390848e9*pow(S,2) + 3.293136e9*pow(S,3) - 1.3172544e9*pow(S,4) + 3.0735936e8*pow(S,5) - 4.390848e7*pow(S,6) + 3.9204e6*pow(S,7) - 217800.*pow(S,8) + 7260.*pow(S,9) - 132.*pow(S,10) + pow(S,11)))/pow(2.718281828459045,0.5*S),(2.08767569878681e-9*(4.790016e8 - 5.7480192e9*S + 1.58070528e10*pow(S,2) - 1.7563392e10*pow(S,3) + 9.879408e9*pow(S,4) - 3.16141056e9*pow(S,5) + 6.1471872e8*pow(S,6) - 7.527168e7*pow(S,7) + 5.8806e6*pow(S,8) - 290400.*pow(S,9) + 8712.*pow(S,10) - 144.*pow(S,11) + pow(S,12)))/pow(2.718281828459045,0.5*S)}};
// S == lambda^2
                double pref = pow(e, -0.5*S);
                for (int i=0; i<n;i++){
                    for (int j=0;j<n;j++){

                        double tot = 0;
                        for (int k =0; k <= j; k++){
                            if (i-j+k >= 0 ){
                                double t1 = pow(-sqrt(S),k) *pow(sqrt(S),i-j+k) ;
                                //double t2 = sqrt(factorial(i)*factorial(j)) /(factorial(j-k) *factorial(k) *factorial(i-j+k) );
                                double tt = (factln (i) + factln(j))/2 - factln (j-k) - factln(k) -factln(i-j+k);
                                double t2 = ( exp(tt));
                                double t3 = sqrt(factorial(i)*factorial(j)) /(factorial(j-k) *factorial(k) *factorial(i-j+k) );
                                tot += (t1*t2);
                            }
                            else {
                                tot+=0;
                            }
                        }
                        _FC[i][j] = pref *tot;
                    }
                }

	}
}
