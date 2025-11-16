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
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Vector<double, Eigen::Dynamic> VectorXd;

cfarray<Polynomial> P;
cfarray<Polynomial> S;
cfarray<Polynomial> Sprime;

Real E(int j, Real ax);
    
// return mth component of Psi_{ijkl} at (x,y,z), i.e. (Psi{ijkl}_m)(x,y,z)
Real Psi(int i, int j, int k, int l, int m, Real x, Real y, Real z, Real alpha, Real gamma);

// Some new save/load functions to ease pain of comment character diffs btwn chflow and Julia
MatrixXi readijkl(const string& filename, char comment_char);
void save(const MatrixXd& A, const string& filename, char comment_char);
void save(const VectorXd& x, const string& filename, char comment_char);
void load(MatrixXd& A, const string& filename, char comment_char);
void load(VectorXd& x, const string& filename, char comment_char);


#endif // BASISFUNCTIONS_H
