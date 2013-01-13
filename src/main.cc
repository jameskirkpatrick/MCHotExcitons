#include <iostream>
#include <fstream>

#include "State.h"
#include "Neighbours.h"
#include "MC.h"

#include <boost/lexical_cast.hpp>

using namespace boost;

using namespace std;

int main(int argc , char * argv[]){
	
    vector <double> collected;
    vector <double> recomb;
    Parameters par;
    const char * xmlf = argv[1];
    par.Init(xmlf);

    Neighbours neigh(&par);
    cout << "Finishe computing fixed rates" <<endl;

    string s_phonons=par.GetNStartPhonon();
    stringstream snphot;
    snphot << s_phonons;
    int n_phonon;
    int nrun =lexical_cast<int> (argv[2]);
    while (snphot>> n_phonon){
        cout << endl << "Start with so many phonons " << n_phonon << endl;
        MC simexc(&par, &neigh);
        State a= par.NN(n_phonon);
        string nameout = lexical_cast<string>(n_phonon)+string(argv[1]) + string(".out");
        simexc.Run(a, nrun, nameout);

        collected.push_back(simexc.PCollected());
        recomb.push_back(simexc.PRecomb());

        if (argc==4 ){
            string namerates=lexical_cast<string>(n_phonon)+string(argv[1]) + string("checkrates.tmp");
            ofstream out(namerates.c_str());
            neigh.PrintCalcRates(out);
            out.close();
        }
    }

        /*cout << '\n' << "Save rates" <<endl;
        ofstream out("checkrates.tmp");
	neigh.PrintAllRates(out);
        out.close();
        out.open("probabilites.tmp");
	neigh.PrintProbs(out);
        out.close();
        out.open("barerates.tmp");
        neigh.PrintCalcRates(out );
        out.close();*/

}
