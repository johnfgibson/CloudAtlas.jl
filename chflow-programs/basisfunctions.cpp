#include "basisfunctions.h"

MatrixXi readijkl(const string& filebase) {
    ifstream is;
    string filename = ifstreamOpen(is, filebase, ".asc");
    if (!is) {
        cerr << "readijkl(filebase) : can't open file " << filebase << " or " << (filebase + ".asc")
             << endl;
        exit(1);
    }

    // Read in header. Form is "% N"
    string comment;
    int N = 0;
    is >> comment >> N;

    cout << "reading N == " << N << " ijkl indices" << endl;
    MatrixXi ijkl(N,4);
    int i,j,k,l;
    for (int n=0; n<N; ++n) {
	is >> ijkl(n,0) >> ijkl(n,1) >> ijkl(n,2) >> ijkl(n,3);
	//cout << ijkl(n,0) << ' ' << ijkl(n,1) << ' ' << ijkl(n,2) << ' ' << ijkl(n,3) << endl;
	
	if (!is.good()) {
	    cerr << "warning: bad istream in readijkl from file " << filebase << endl;
	    is.close();
	}
										 
    }
    //cout << "ikl == " << endl << ijkl << endl;
    return ijkl;
}

Real E(int j, Real ax) {
    Real Ejax = 1.0;
    if (j<0)
	Ejax = cos(j*ax);
    else if (j>0)
	Ejax = sin(j*ax);
    return Ejax;
}

// return mth component of Psi_{ijkl} at (x,y,z), i.e. (Psi{ijkl}_m)(x,y,z)
// m uses zero-based indices, following channelflow vectors
// i uses one-based indices, following Julia code used to develop Psi
Real Psi(int i, int j, int k, int l, int m, Real x, Real y, Real z, Real alpha, Real gamma) { 
    Real rtn = 0.0;

    if (i==1 && m == 0)
	rtn = E(k, gamma*z) * Sprime[l](y);

    else if (i==2 && m == 1)
	rtn = gamma*k* E(k, gamma*z) * S[l](y);
    else if (i==2 && m == 2)
	rtn = E(-k, gamma*z) * Sprime[l](y);
    
    else if (i==3 && m == 2)
	rtn = E(j, alpha*x) * Sprime[l](y);
    
    else if (i==4 && m == 0)
	rtn =  E(-j, alpha*x) * Sprime[l](y);
    else if (i==4 && m == 1)
	rtn =  alpha*j * E(j, alpha*x) * S[l](y);
    
    else if (i==5 && m == 0)
	rtn =  gamma*k * E(-j, alpha*x) * E(k, gamma*z) * Sprime[l](y);
    else if (i==5 && m == 2)
	rtn = -alpha*j * E(j, alpha*x) * E(-k, gamma*z) * Sprime[l](y);

    else if (i==6 && m == 0)
	rtn = gamma*k * E(-j, alpha*x) * E(k, gamma*z) * Sprime[l](y);
    else if (i==6 && m == 1)
	rtn = 2*alpha*gamma*j*k * E(j, alpha*x) * E(k, gamma*z) * S[l](y);
    else if (i==6 && m == 2)
	rtn = alpha*j * E(j, alpha*x) * E(-k, gamma*z) * Sprime[l](y);

    return rtn;
}


MatrixXi readijkl(const string& filebase, char comment_char) {
    ifstream is;
    string filename = ifstreamOpen(is, filebase, ".asc");
    if (!is) {
        cerr << "readijkl(filebase) : can't open file " << filebase << " or " << (filebase + ".asc")
             << endl;
        exit(1);
    }

    // Read in header. Form is "% N"
    string comment;
    int N = 0;
    is >> comment >> N;

    cout << "reading N == " << N << " ijkl indices" << endl;
    MatrixXi ijkl(N,4);
    int i,j,k,l;
    for (int n=0; n<N; ++n) {
	is >> ijkl(n,0) >> ijkl(n,1) >> ijkl(n,2) >> ijkl(n,3);
	//cout << ijkl(n,0) << ' ' << ijkl(n,1) << ' ' << ijkl(n,2) << ' ' << ijkl(n,3) << endl;
	
	if (!is.good()) {
	    cerr << "warning: bad istream in readijkl from file " << filebase << endl;
	    is.close();
	}
										 
    }
    //cout << "ikl == " << endl << ijkl << endl;
    return ijkl;
}

void save(const MatrixXd& A, const string& filebase, char comment_char) {
    int taskid = mpirank();

    if (taskid != 0) {
        return;
    }
    std::string filename = appendSuffix(filebase, ".asc");
    std::ofstream os(filename.c_str());
    os << std::setprecision(17);
    os << comment_char << ' ' << A.rows() << ' ' << A.cols() << '\n';
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j)
            os << A(i, j) << ' ';
        os << '\n';
    }
}


void save(const VectorXd& x, const string& filebase, char comment_char) {
    int taskid = mpirank();
    if (taskid != 0) {
        return;
    }
    std::string filename = appendSuffix(filebase, ".asc");
    std::ofstream os(filename.c_str());
    os << std::setprecision(17);
    os << comment_char << ' ' << x.size() << '\n';
    for (int i = 0; i < x.size(); ++i)
        os << x(i) << '\n';
}

void load(MatrixXd& A, const string& filebase, char comment_char) {
    std::string filename = appendSuffix(filebase, ".asc");
    std::ifstream is(filename.c_str());
    int M, N;
    char c;
    is >> c;
    if (c != comment_char) {
        std::string message("load(Matrix&, filebase): bad header in file, wrong comment character");
        message += filename;
        cferror(message);
    }
    is >> M >> N;
    A.resize(M, N);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            is >> A(i, j);
}

void load(VectorXd& x, const string& filebase, char comment_char) {
    std::string filename = appendSuffix(filebase, ".asc");
    std::ifstream is(filename.c_str());
    int N;
    char c;
    is >> c;
    if (c != comment_char) {
        std::string message("load(VectorXd&, filebase): bad header in file, wrong comment character");
        message += filename;
        cferror(message);
    }
    is >> N;
    x.resize(N);
    for (int i = 0; i < x.size(); ++i)
        is >> x(i);
}
