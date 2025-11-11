#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

#include <stdlib.h>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "cfbasics/mathdefs.h"
#include "cfbasics/cfvector.h"

#include "polynomial.h"

using namespace std;
using namespace chflow;

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;

cfarray<Polynomial> P;
cfarray<Polynomial> S;
cfarray<Polynomial> Sprime;

Real E(int j, Real ax);
    
// return mth component of Psi_{ijkl} at (x,y,z), i.e. (Psi{ijkl}_m)(x,y,z)
Real Psi(int i, int j, int k, int l, int m, Real x, Real y, Real z, Real alpha, Real gamma);

MatrixXi readijkl(const string& filename);

#endif // BASISFUNCTIONS_H
