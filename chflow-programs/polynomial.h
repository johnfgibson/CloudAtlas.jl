#ifndef CHANNELFLOW_POLYNOMIAL_H
#define CHANNELFLOW_POLYNOMIAL_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"

namespace chflow {

class Polynomial {
public:
    Polynomial();
    Polynomial(int order);
    Polynomial(const Vector& coeffs);
    Polynomial(Real c0, Real c1=0.0, Real c2=0.0, Real c3=0.0, Real c4=0.0, Real c5=0.0);

    Polynomial(const Polynomial& p);
    Polynomial& operator=(const Polynomial& p);
    
    Real operator()(Real x) const;   // evaluate polynomial at x

    Polynomial& operator/=(Real a);  // divide by const
    Polynomial& operator*=(Real a);  // mult by const

    int order() const;   // polynomial order       (zero or greater)     
    int length() const;  // length of coeff vector (can be zero for null poly P(x) == 0)

    // could be private
    Vector c;        // polynomial coeffs, P(x) = c0 + c1 x + c2 x^2 + ...
};

Polynomial operator*(Real a, const Polynomial& p); 
Polynomial operator+(const Polynomial& p, const Polynomial& q); 
Polynomial operator-(const Polynomial& p, const Polynomial& q); 
Polynomial operator*(const Polynomial& p, const Polynomial& q); 

Polynomial diff(Polynomial& p);                     // return dp/dx
Polynomial integrate(const Polynomial& p, int n=1); // return in definite integral of p(x), setting const to zero

Real L2innerproduct(const Polynomial& p, const Polynomial& q);
Real L2Norm2(const Polynomial& p, const Polynomial& q);
Real L2Dist2(const Polynomial& p, const Polynomial& q);
Real L2Norm2(const Polynomial& p, const Polynomial& q);
Real L2Dist2(const Polynomial& p, const Polynomial& q);



std::ostream& operator<<(std::ostream& os, const Polynomial& p);

}
#endif 
