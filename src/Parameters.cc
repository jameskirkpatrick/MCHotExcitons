#include "Parameters.h"
#include <algorithm>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <libxml2/libxml/tree.h>

using namespace boost;

void Parameters::Read(const char * vertex, const char * edge){
	ifstream ver(vertex);
	ifstream edg(edge);
        string flushline;
	//read in the edges
        if ( !ver || !edg){
            throw runtime_error(string("could not find the vertex file: ") +string(vertex)
                    +string(" or edge file") +string(edge) );
        }
	while (ver){
		double nrg;
		ver >> nrg;
		_energies.push_back(nrg);
                getline(ver,flushline);
	}
	while (edg){
		int from, to;
		double J;
		edg >> from >> to >> J;
		_top.push_back(make_pair(from, to));
		_transfers.push_back(J);
                getline(edg,flushline);
	}
}
void Parameters::Read(const char * vertex, const char * edge, const char *proper){
    Read(vertex, edge);
    ifstream prop(proper);
    string flushline;
    while (prop){
        string property ;
        prop >> property;
        _properties.push_back(property);
        getline(prop, flushline);
    }
}
State Parameters::NN(const int &n){
    int nmol = _properties.size();
    State res(nmol);
    int *data;
    data = new int [nmol];
    vector <string>::iterator itmin = find(_properties.begin() , _properties.end(), string("Generator"));
    int imin = itmin - _properties.begin();
    for (int i=0;i<nmol;i++){
        if (i==imin) data[i] =n;
        else data[i]=0;
    }
    res.SetEl(imin);
    res.SetNucl(data);
    delete [] data;
    return res;
}

void Parameters::Print(string & name){
    ofstream out ((name+ string("nrg.dat")).c_str());
    vector < double>::iterator itd;
    for (itd = _energies.begin(); itd != _energies.end(); ++itd){
        out << *itd <<'\n';
    }
    out.close();

    out.open((name+string("prop.dat")).c_str());
    vector <string>::iterator its;
    for (its = _properties.begin(); its!= _properties.end();++its){
        out << *its <<'\n';
    }
    out.close();

    out.open((name+string("edges.dat")).c_str());
    vector<pair<int, int> >::iterator itp;
    for (itp = _top.begin(),itd=_transfers.begin(); itp != _top.end() ; ++itp, ++itd){
        out << itp-> first << '\t' << itp->second << '\t' << *itd<<'\n';
    }
    out.close();
}


void Parameters::Init(string namexml){
    xmlDocPtr doc;
    xmlNodePtr node;
    xmlChar *key;

    /* this var can take values: 0,1,2
     * 0 == read
     * 1 == init2crg
     * 2 == init1d
     * 3 == init3d
     */
    int internalinit=0;
    doc = xmlParseFile(namexml.c_str());
    if (doc == NULL)
        throw runtime_error(string("Error on opening file") + namexml);
    node = xmlDocGetRootElement(doc);

    if(node == NULL) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, empty xml document: " + namexml);
    }

    if(xmlStrcmp(node->name, (const xmlChar *)"MCParameters")) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, xml file not labeled crgunit_type: " + namexml);
    }

    string nameen, nameed, namepro;
    double J,F, Jrec, eps, eps1, eps2, gsnrg,d,sigma;
    
    _sigmaJ=0.;
    int n;
    // parse xml tree
    for(node = node->xmlChildrenNode; node != NULL; node = node->next) {
            key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);

            if (!xmlStrcmp( node->name, (const xmlChar *) "S")){
                _S = lexical_cast<double> (key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "kT")){
                _kT = lexical_cast<double> (key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "sigmaJ")){
                _sigmaJ = lexical_cast<double> (key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "hbaromega")){
                _hbarOmega = lexical_cast<double> (key);
            }

            else if (!xmlStrcmp( node->name, (const xmlChar *) "lambdaO")){
                _lambdaO = lexical_cast<double>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "ratevib")){
                _ratevib = lexical_cast<double>(key);
            }

            else if (!xmlStrcmp( node->name, (const xmlChar *) "tmax")){
                _tmax = lexical_cast<double>(key);
            }

            else if (!xmlStrcmp( node->name, (const xmlChar *) "tau")){
                _tau = lexical_cast<double>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "dt")){
                _dt = lexical_cast<double>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "nameenergy")){
                nameen = lexical_cast<string>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "nameedges")){
                nameed = lexical_cast<string>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "nameproperties")){
                namepro = lexical_cast<string>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "nstartphonon")){
                _nstartphonon = lexical_cast<string>(key);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "seed")){
                int seed = lexical_cast<int>(key);
                srand48(seed);
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "generate2Crg")){
                xmlNodePtr noder2 = node->xmlChildrenNode;

                internalinit=1;
                for (; noder2 != NULL ; noder2= noder2->next){
                    key = xmlNodeListGetString(doc, noder2->xmlChildrenNode, 1);
                    if (!xmlStrcmp(noder2->name , (const xmlChar*) "J")){
                       J = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "n")){
                       n = lexical_cast<int> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "d")){
                       d = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "F")){
                       F = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "Jrec")){
                       Jrec = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "eps")){
                       eps = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "gsnrg")){
                       gsnrg = lexical_cast<double> (key);
                    }

                }
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "generate1D")){
                xmlNodePtr noder2 = node->xmlChildrenNode;

                internalinit=2;
                for (; noder2 != NULL ; noder2= noder2->next){
                    key = xmlNodeListGetString(doc, noder2->xmlChildrenNode, 1);
                    if (!xmlStrcmp(noder2->name , (const xmlChar*) "J")){
                       J = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "n")){
                       n = lexical_cast<int> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "d")){
                       d = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "F")){
                       F = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "Jrec")){
                       Jrec = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "eps")){
                       eps = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "gsnrg")){
                       gsnrg = lexical_cast<double> (key);
                    }

                }
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "generate1D")){
                xmlNodePtr noder2 = node->xmlChildrenNode;

                internalinit=2;
                for (; noder2 != NULL ; noder2= noder2->next){
                    key = xmlNodeListGetString(doc, noder2->xmlChildrenNode, 1);
                    if (!xmlStrcmp(noder2->name , (const xmlChar*) "J")){
                       J = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "n")){
                       n = lexical_cast<int> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "d")){
                       d = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "F")){
                       F = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "Jrec")){
                       Jrec = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "eps")){
                       eps = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "gsnrg")){
                       gsnrg = lexical_cast<double> (key);
                    }

                }
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "generateS")){
                xmlNodePtr noder2 = node->xmlChildrenNode;

                internalinit=4;
                for (; noder2 != NULL ; noder2= noder2->next){
                    key = xmlNodeListGetString(doc, noder2->xmlChildrenNode, 1);
                    if (!xmlStrcmp(noder2->name , (const xmlChar*) "J")){
                       J = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "n")){
                       n = lexical_cast<int> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "d")){
                       d = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "Jrec")){
                       Jrec = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "eps")){
                       eps = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "gsnrg")){
                       gsnrg = lexical_cast<double> (key);
                    }

                }
            }
            else if (!xmlStrcmp( node->name, (const xmlChar *) "generate3D")){
                xmlNodePtr noder2 = node->xmlChildrenNode;

                internalinit=3;
                for (; noder2 != NULL ; noder2= noder2->next){
                    key = xmlNodeListGetString(doc, noder2->xmlChildrenNode, 1);
                    if (!xmlStrcmp(noder2->name , (const xmlChar*) "J")){
                       J = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "n")){
                       n = lexical_cast<int> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "d")){
                       d = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "Jrec")){
                       Jrec = lexical_cast<double> (key);
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "eps")){
                       eps = lexical_cast<double> (key);
                       eps1 = eps;
                       eps2 = eps;
                    }
                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "gsnrg")){
                       gsnrg = lexical_cast<double> (key);
                    }

                    else if (!xmlStrcmp(noder2->name , (const xmlChar*) "sigma")){
                       sigma = lexical_cast<double> (key);
                    }

                }
            }
    }
    if( internalinit==0)    Read(nameen.c_str(), nameed.c_str(), namepro.c_str());
    else if ( internalinit==1){
        Init2Charges(n,d,eps,F, J,Jrec, gsnrg);
        
        cout << "Saving 2 Charges data" <<endl;
        string out ("tmp");
        Print(out);
        cout << "done " <<endl;
    }
    else if ( internalinit==2){
        Init1D(n,d,eps,F, J,Jrec, gsnrg);
     
        cout << "Saving 1D Data" <<endl;
        string out ("save");
        Print(out);
        cout << "done " <<endl;
    }
    else if (internalinit==3){
        eps1 /= (27.21 * 0.529);
        eps2 /= (27.21 * 0.529);
        NonIsoPol fun(F, eps1, eps2, sigma);
        Init3D(n,d, fun, J, Jrec,gsnrg );
        cout << "Saving 3D Data" <<endl;
        string out ("save");
        Print(out);
        cout << "done " <<endl;
    }
    else if (internalinit==4){
        
        Init1DS(n,d,eps,_kT, J,Jrec, gsnrg);
        cout << "Saving S Data" <<endl;
        string out ("save");
        Print(out);
        cout << "done " <<endl;
    }
}
vector<Parameters::_point> Parameters::neigh1Crg(_point & a, const int & n){
    vector<_point> res = neigh1Crg(a);
    vector<_point>::iterator itn= res.begin();
    for ( ; itn != res.end(); itn++){
        if( (*itn)._x1  < 0 ) (*itn)._x1 +=n;
        if( (*itn)._y1  < 0 ) (*itn)._y1 +=n;
        if( (*itn)._x1  >= n ) (*itn)._y1 -=n;
        if( (*itn)._y1  >= n ) (*itn)._y1 -=n;
    }
    return res;
}

vector<Parameters::_point> Parameters::neigh1Crg(_point & a){
    vector <_point> res;
    _point n1, n2, n3,n4,n5,n6;

    n1 = a;
    n1._x1 -= 1;
    n2 = a;
    n2._x1 += 1;

    n3 = a;
    n3._y1 -= 1;
    n4 = a;
    n4._y1 += 1;

    n5 =a;
    n5._z1 -=1;
    n6 = a;
    n6._z1 +=1;

    res.push_back(n1);
    res.push_back(n2);
    res.push_back(n3);
    res.push_back(n4);
    res.push_back(n5);
    res.push_back(n6);

    return res;
}

const Parameters::_point & Parameters::_point::operator =( Parameters::_point &r){
    this->_x1 = r._x1;
    this->_y1 = r._y1;
    this->_z1 = r._z1;
    this->_x2 = r._x2;
    this->_y2 = r._y2;
    this->_z2 = r._z2;
    return *this;
}

vector<Parameters::_point> Parameters::neigh(_point & a){
    vector <_point> res;
    _point n1, n2, n3,n4,n5,n6,n7,n8,n9,n10,n11,n12;

    n1 = a;
    n1._x1 -= 1;
    n2 = a;
    n2._x1 += 1;

    n3 = a;
    n3._x2 -=1;
    n4 = a;
    n4._x2 +=1;

    n5 = a;
    n5._y1 -= 1;
    n6 = a;
    n6._y1 += 1;

    n7 = a;
    n7._y2 -=1;
    n8 = a;
    n8._y2 +=1;

    n9 =a;
    n9._z1 -=1;
    n10 = a;
    n10._z1 +=1;

    n11 = a;
    n11._z2 -=1;
    n12 = a;
    n12._z2 +=1;

    res.push_back(n1);
    res.push_back(n2);
    res.push_back(n3);
    res.push_back(n4);
    res.push_back(n5);
    res.push_back(n6);
    res.push_back(n7);
    res.push_back(n8);
    res.push_back(n9);
    res.push_back(n10);
    res.push_back(n11);
    res.push_back(n12);

    return res;
}

void Parameters::Init2Charges(const int & n, const double &d,
        const double &eps, const double &F,const double &J, const double &Jrec, const double & gsnrg  ){
    
    vector < _point> allpoints;
    vector < _point> recomb;
    _point gen;
    gen._x1 = gen._x2 = gen._y1 = gen._y2= (n-1)/2;
    gen._z1 = gen._z2 = 0;

    /// generate all the points 
    _point tp;
    cout << "Generate all points" <<endl;
    for (tp._x1 =0 ; tp._x1 < n; tp._x1++){
        for (tp._y1 =0 ; tp._y1 < n; tp._y1++){
            for (tp._z1 =0 ; tp._z1 < n; tp._z1++){
                for (tp._x2 =0 ; tp._x2 < n; tp._x2++){
                    for (tp._y2 =0 ; tp._y2 < n; tp._y2++){
                        for (tp._z2 =0 ; tp._z2 < n; tp._z2++){
                            allpoints.push_back(tp);
                        }
                    }
                }
            }
        }
    }

    cout << "Generate all energies and properties" <<endl;
    vector < _point>::iterator itp;
    for (itp = allpoints.begin(); itp!= allpoints.end() ; itp++ ){
        //energy
        double dist = d * sqrt(double((itp->_x1-itp->_x2)*(itp->_x1-itp->_x2)+
        (itp->_y1-itp->_y2)*(itp->_y1-itp->_y2)+
        (itp->_z1+itp->_z2+1)*(itp->_z1+itp->_z2+1)) );
        double energy = -27.21 * 0.529 / ( eps * dist) - F * dist;
        _energies.push_back(energy);

        //properties
        if (itp->_z1 == n-1 || itp->_z2 == n-1){
            _properties.push_back(string("Collector"));
        }
        else if (*itp == gen){
            _properties.push_back(string("Generator"));
        }
        else {
            _properties.push_back(string ("Nothing"));
        }

        //a potential recombination site?
        if (itp->_z1 == itp->_z2 && itp->_x1 == itp->_x2 && itp->_y1 == itp->_y2){
            recomb.push_back(*itp);
        }

    }

    cout << "now make all neighbours" <<endl;
    for (itp = allpoints.begin(); itp!= allpoints.end() ; itp++ ){

        //neighbours
        vector <_point> neighbour = neigh(*itp);
        vector <_point>::iterator itnp =neighbour.begin();
        for (; itnp != neighbour.end(); itnp++){
            if (itnp->_x1 < n &&itnp->_x2 < n &&itnp->_y1 < n && itnp->_y2 < n && itnp->_z1 < n && itnp->_z2 < n &&
                itnp->_x1 >= 0 &&itnp->_x2 >= 0 &&itnp->_y1 >= 0 && itnp->_y2 >= 0 && itnp->_z1 >= 0 && itnp->_z2 >= 0){
                vector<_point>::iterator ittmp = find(allpoints.begin(), allpoints.end(), *itnp);
                if (ittmp== allpoints.end()){
                    cerr << "What the fuck I cannot fintd the following point:" <<
                            itnp->_x1 << '\t' << itnp->_y1 << '\t' <<  itnp->_z1 << '\t' <<
                            itnp->_x2 << '\t' << itnp->_y2 << '\t' <<  itnp->_z2 <<endl;
                }
                int from = itp - allpoints.begin();
                int to  =ittmp - allpoints.begin();
                if (to > from){
                    double thisJ=J*exp(Gaussian(0, _sigmaJ));
                    _top.push_back(make_pair(from, to));
                    _transfers.push_back(thisJ);
                    _top.push_back(make_pair(to,from));
                    _transfers.push_back(thisJ);
                }
            }
        }
        neighbour.clear();
    }
    _energies.push_back(gsnrg);
    _properties.push_back("GroundState");

    vector <_point>::iterator itr = recomb.begin();
    for ( ; itr != recomb.end(); ++itr){

        int igen = find(allpoints.begin(), allpoints.end() , *itr) - allpoints.begin();
        _top.push_back(make_pair(igen, _energies.size() - 1));
        _top.push_back(make_pair(_energies.size() - 1, igen));
        _transfers.push_back(Jrec);
        _transfers.push_back(Jrec);
    }
}

void Parameters::Init1DS(const int & n, const double &d,
        const double &eps, const double &kT,
        const double &J, const double &Jrec,
        const double & gsnrg  ){

    cout << "Generate all energies and properties" <<endl;


    for (int i =0; i< n;i++){
        //energy
        double dist = d *(i+1);

        double energy = -27.21 * 0.529 / ( eps * dist) - kT * log(2*Pi*double((i+1)*(i+1))) ;
        _energies.push_back(energy);

        //properties
        if (i==0){
            _properties.push_back(string("Generator"));
            _top.push_back(make_pair(i, i+1));
            _transfers.push_back(J);
        }
        else if (i==n-1){
            _properties.push_back(string("Collector"));
            _top.push_back(make_pair(i, i-1));
            _transfers.push_back(J);
        }
        else {
            _properties.push_back(string ("Nothing"));
            _top.push_back(make_pair(i, i+1));
            _transfers.push_back(J);
            _top.push_back(make_pair(i, i-1));
            _transfers.push_back(J);
        }
    }

    _energies.push_back(gsnrg);
    _properties.push_back("GroundState");

    _top.push_back(make_pair(n,0));
    _transfers.push_back(Jrec);
    _top.push_back(make_pair(0,n));
    _transfers.push_back(Jrec);

}

void Parameters::Init1D(const int & n, const double &d,
        const double &eps, const double &F,
        const double &J, const double &Jrec,
        const double & gsnrg  ){

    cout << "Generate all energies and properties" <<endl;


    for (int i =0; i< n;i++){
        //energy
        double dist = d *(i+1);
        
        double energy = -27.21 * 0.529 / ( eps * dist) - F * dist;
        _energies.push_back(energy);

        //properties
        if (i==0){
            _properties.push_back(string("Generator"));
            _top.push_back(make_pair(i, i+1));
            _transfers.push_back(J);
        }
        else if (i==n-1){
            _properties.push_back(string("Collector"));
            _top.push_back(make_pair(i, i-1));
            _transfers.push_back(J);
        }
        else {
            _properties.push_back(string ("Nothing"));
            _top.push_back(make_pair(i, i+1));
            _transfers.push_back(J);
            _top.push_back(make_pair(i, i-1));
            _transfers.push_back(J);
        }
    }

    _energies.push_back(gsnrg);
    _properties.push_back("GroundState");
    
    _top.push_back(make_pair(n,0));
    _transfers.push_back(Jrec);
    _top.push_back(make_pair(0,n));
    _transfers.push_back(Jrec);
    
}

 void Parameters::Init3D(const int & n, const double & a, EnergyFunction & EF, const double & J,
        const double & Jrec, const double & gsnrg){
    // generate all points
    vector < _point> allpoints;
    vector < _point> recomb;
    _point gen;
    gen._x1 = gen._x2 = gen._y1 = gen._y2= (n-1)/2;
    gen._z1 = gen._z2 = 0;

    /// generate all the points
    _point tp;
    tp._x2 = (n-1)/2;
    tp._y2 = (n-1)/2;
    tp._z2 = 0;
    cout << "Generate all points" <<endl;
    for (tp._x1 =0 ; tp._x1 < n; tp._x1++){
        for (tp._y1 =0 ; tp._y1 < n; tp._y1++){
            for (tp._z1 =0 ; tp._z1 < n/2; tp._z1++){
                allpoints.push_back(tp);
            }
        }
    }
    /// energies
    vector < _point>::iterator itp= allpoints.begin();
    for ( ;itp!= allpoints.end(); ++itp){
        double z1 = a *double(itp->_z1) + a/2.;
        double z2 = a/2.;
        double dr = a * sqrt(double (itp->_x1 -itp->_x2 )*(itp->_x1 -itp->_x2 )+
                             double (itp->_y1 -itp->_y2 )*(itp->_y1 -itp->_y2 ));
        _energies.push_back(EF(dr,z1,z2));

         //properties
        if (itp->_z1 == n/2-1 ||itp->_x1 == 0||itp->_x1 == n-1
                ||itp->_y1 == 0||itp->_y1 == n-1){
            _properties.push_back(string("Collector"));
        }
        else if (*itp == gen){
            _properties.push_back(string("Generator"));
        }
        else {
            _properties.push_back(string ("Nothing"));
        }

        //a potential recombination site?
        if (itp->_z1 == itp->_z2 && itp->_x1 == itp->_x2 && itp->_y1 == itp->_y2){
            recomb.push_back(*itp);
        }
    }

    /// neighbours
    for (itp = allpoints.begin() ;itp!= allpoints.end(); ++itp){
         //neighbours
        vector <_point> neighbour = neigh1Crg(*itp,n);
        vector <_point>::iterator itnp =neighbour.begin();
        for (; itnp != neighbour.end(); itnp++){
            if (itnp->_x1 < n  &&itnp->_y1 < n && itnp->_z1 < n/2  &&
                itnp->_x1 >= 0 &&itnp->_y1 >= 0&& itnp->_z1 >= 0 ){
                vector<_point>::iterator ittmp = find(allpoints.begin(), allpoints.end(), *itnp);
                if (ittmp== allpoints.end()){
                    cerr << "What the fuck I cannot fintd the following point:" <<
                            itnp->_x1 << '\t' << itnp->_y1 << '\t' <<  itnp->_z1 << '\t' <<
                            itnp->_x2 << '\t' << itnp->_y2 << '\t' <<  itnp->_z2 <<endl;
                }
                int from = itp - allpoints.begin();
                int to  =ittmp - allpoints.begin();
                if (to > from){
                    double thisJ=J*exp(Gaussian(0, _sigmaJ));
                    _top.push_back(make_pair(from, to));
                    _transfers.push_back(thisJ);
                    _top.push_back(make_pair(to,from));
                    _transfers.push_back(thisJ);
                }
            }
        }
        neighbour.clear();
    }

    _energies.push_back(gsnrg);
    _properties.push_back("GroundState");

    vector <_point>::iterator itr = recomb.begin();
    for ( ; itr != recomb.end(); ++itr){

        int igen = find(allpoints.begin(), allpoints.end() , *itr) - allpoints.begin();
        _top.push_back(make_pair(igen, _energies.size() - 1));
        _top.push_back(make_pair(_energies.size() - 1, igen));
        _transfers.push_back(Jrec);
        _transfers.push_back(Jrec);
    }
 }
