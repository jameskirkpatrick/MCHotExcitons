/* 
 * File:   EnergyFunction.h
 * Author: james
 *
 * Created on March 23, 2010, 4:10 PM
 */

#ifndef _ENERGYFUNCTION_H
#define	_ENERGYFUNCTION_H

class EnergyFunction{
public:
    EnergyFunction(){}
    virtual ~EnergyFunction(){}
    double virtual operator ()(const double & dx, const double & dy, const double & dz){
        cerr << "Should always call the derived operator" <<endl;
    }

};


#endif	/* _ENERGYFUNCTION_H */

