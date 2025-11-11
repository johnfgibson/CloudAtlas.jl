#include "polynomial.h"

namespace chflow {

Polynomial::Polynomial()
    : c(0)
{}

Polynomial::Polynomial(int order)
    : c(order+1)
{}

Polynomial::Polynomial(const Vector& coeffs)
    : c(coeffs)
{}

Polynomial::Polynomial(const Polynomial& p)
    : c(p.c)
{}

// array literals in C++ are so awkward, hence this clunkiness
Polynomial::Polynomial(Real c0, Real c1, Real c2, Real c3, Real c4, Real c5)
    :
    c(0)
{

    if (c5 != 0.0) {
	c.resize(6);
	c[0]=c0; c[1]=c1; c[2]=c2; c[3]=c3; c[4]=c4; c[5]=c5;
    }
    else if (c4 != 0.0) {
	c.resize(5);
	c[0]=c0; c[1]=c1; c[2]=c2; c[3]=c3; c[4]=c4;
    }
    else if (c3 != 0.0) {
	c.resize(4);
	c[0]=c0; c[1]=c1; c[2]=c2; c[3]=c3;
    }
    else if (c2 != 0.0) {
	c.resize(3);
	c[0]=c0; c[1]=c1; c[2]=c2;
    }
    else if (c1 != 0.0) {
	c.resize(2);
	c[0]=c0; c[1]=c1;
    }
    else if (c0 != 0.0) {
	c.resize(1);
	c[0]=c0;
    }
    else
	c.resize(0);
}

Polynomial& Polynomial::operator=(const Polynomial& p) {
    c = p.c;
}


int Polynomial::order() const {
    return c.length() == 0 ? 0 : c.length()-1;
}


int Polynomial::length() const {
    return c.length();
}


// P(x) = 0                             N==0
// P(x) = c0                            N==1
// P(x) = c0 + x*(c1)                   N==2
// P(x) = c0 + x*(c1 + x*(c2)           N==3
// P(x) = c0 + x*(c1 + x*(c2 + x*(c3))  N==4

// for N==4 case
// n == N:    px == 0
// N == N-1:  
Real Polynomial::operator()(Real x) const {
    int N = c.length();
    Real px = 0.0;
    for (int n=N-1; n>=0; --n) 
	px = c[n] + x*px;
    
    return px;
}

Polynomial& Polynomial::operator/=(Real a) {
    c /= a;
    return *this;
}
Polynomial& Polynomial::operator*=(Real a) {
    c *= a;
    return *this;
}
Polynomial operator*(Real a, const Polynomial& p) {
    Polynomial ap(p);
    ap *= a;
    return ap;
}

Polynomial operator+(const Polynomial& p_, const Polynomial& q_) {
    // relabel to make p the higher-order polynomial
    const Polynomial& p = p_.length() > q_.length() ? p_ : q_;
    const Polynomial& q = p_.length() > q_.length() ? q_ : p_;
    
    Polynomial s(p); // copy p into s, then add q
    for (int n=0; n<q.length(); ++n)
	s.c[n] += q.c[n];
    return s;
}

Polynomial operator-(const Polynomial& p_, const Polynomial& q_) {
    // relabel to make p the higher-order polynomial
    const Polynomial& p = p_.length() > q_.length() ? p_ : q_;
    const Polynomial& q = p_.length() > q_.length() ? q_ : p_;
    
    Polynomial s(p); // copy p into s, then subtract q
    for (int n=0; n<q.length(); ++n)
	s.c[n] -= q.c[n];
    return s;
}


//  p(x) = a0 + a1 x + a2 x^2
//  q(x) = b0 + b1 x + b2 x^2
// pq(x) = a0*b0 + (a0 b1 + a1 b0) x + (a0 b2 + a1 b1 + a2 b0) x^2 + (a1 b2 + a2 b1 + a2 b0) x^3 + a2 b2 x^4

Polynomial operator*(const Polynomial& p, const Polynomial& q) {
    int corder = p.order() + q.order();
    int clength = (p.length() > 0 && q.length() > 0) ? corder+1 : 0;

    Vector c(clength);
    for (int n=0; n<p.length(); ++n)
	for (int m=0; m<q.length(); ++m)
	    c[m+n] += p.c[n] * q.c[m];
    return Polynomial(c);
}

// p(x) =     p0 +    p1  x +    p2  x^2  +  p3 x^3     plength == 4
// dpdx = (1) p1 + (2 p2) x + (3 p3) x^2                clength == 3
// dpdx =     c0 +    c1  x +    c2  x^2         
Polynomial diff(Polynomial& p) {
    int plength = p.length();
    int clength = plength > 0 ? plength-1 : 0;

    Vector c(clength);
    for (int n=0; n<clength; ++n)
	c[n] = (n+1)*p.c[n+1];
    return Polynomial(c);
}

// p(x)    = p0 + p1 x +     p2 x^2 +     p3 x^3                 plength == 4
// intp(x) = 0  + p0 x + 1/2 p1 x^2 + 1/3 p2 x^3 + (1/4 p3) x^4  clength == 5
// intp(x) = c0 + c1 x +     c2 x^2 +     c3 x^3 +      c4  x^5
Polynomial integrate(Polynomial& p) {
    int plength = p.length();
    int clength = plength > 0 ? plength+1 : 0;

    Vector c(clength);
    for (int n=1; n<clength; ++n)
	c[n] = p.c[n-1]/n;
    return Polynomial(c);
}

Real L2innerproduct(const Polynomial& p, const Polynomial& q) {
    Polynomial pq = p*q;
    Polynomial integral_pq = integrate(pq);
    return 0.5*(integral_pq(1.0) - integral_pq(-1.0)); // HARDCODED for domain [-1,1]
}

Real L2Norm2(const Polynomial& p) {
    return L2innerproduct(p,p);
}
Real L2Norm(const Polynomial& p) {return sqrt(L2Norm2(p));}

Real L2Dist2(const Polynomial& p, const Polynomial& q) {
    Polynomial d = p-q;
    return L2Norm2(d);
}

Real L2Dist(const Polynomial& p, const Polynomial& q) {return sqrt(L2Dist2(p,q));}


std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
    if (p.length() == 0)
	os << "0";
    else {
	os << p.c[0];
	
	for (int n=1; n<p.length(); ++n) {
	    if (p.c[n] != 0.0) {
		os << ((p.c[n] < 0) ? " - " : " + ");
		if (abs(p.c[n]) != 1.0) 
		    cout << abs(p.c[n]);
		if (n==1)
		    cout << " y";
		else if (n>=2)
		    cout << " y^" << n;
	    }
	}
    }
    return os;
}
		 
    
	
}
