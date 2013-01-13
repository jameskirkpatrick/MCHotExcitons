/* 
 * File:   NonIsoPol.h
 * Author: james
 *
 * Created on March 23, 2010, 4:11 PM
 */

#ifndef _NONISOPOL_H
#define	_NONISOPOL_H
#include "EnergyFunction.h"
#include "global.h"
#include <cmath>
#include <stdlib.h>

class NonIsoPol: public EnergyFunction{
private:
    double _F;
    double _epsR;
    double _epsZ;
    double _sigma;
    
public:
    NonIsoPol(){}
    NonIsoPol(const double &F, const double &epsR, const double &epsZ, const double&sigma):
    _F(F),_epsR(epsR), _epsZ(epsZ), _sigma(sigma){}
    ~NonIsoPol(){}

    double operator()(const double &dx, const double &dy, const double &dz);
};

inline double NonIsoPol::operator()(const double &dr, const double &z1, const double &z2){
    double gamma  = sqrt(_epsR/_epsZ);
    double aveps = 1 / (sqrt(_epsR*_epsZ) + _epsR);
    double reeps = (_epsR - sqrt(_epsR*_epsZ) )/(sqrt(_epsR*_epsZ) + _epsR);
    double res = reeps /(_epsR * 2 * z1)
            - aveps * 1/sqrt(dr*dr+ (z1 *gamma+z2)*(z1 *gamma+z2))
            - aveps * 1/sqrt(dr*dr+(z1+z2)*(z1+z2)*gamma*gamma)
            - _F * (z1+z2);
    double random = Gaussian(0, _sigma);
    return res + random;


}

#endif	/* _NONISOPOL_H */

