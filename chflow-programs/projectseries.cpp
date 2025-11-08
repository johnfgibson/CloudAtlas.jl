// License declaration at end of file

#include <stdlib.h>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/utilfuncs.h"
#include "cfbasics/mathdefs.h"
#include "polynomial.h"
#include "basisfunctions.h"

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) { cfMPI_Init(&argc, &argv); {

  string purpose("produce I,D,E, time series and means from u(x,t) time series");
  ArgList args(argc, argv, purpose);

  const Real T0        = args.getreal("-T0", "--T0", 0,   "start time");
  const Real T1        = args.getreal("-T1", "--T1", 500, "end time");
  /***/ Real dT        = args.getreal("-dT", "--dT", 1.0, "save interval");
   const bool adjustdT = args.getflag("-adT","--adjustdT", "adjust dT so that it evenly divides (T1-T0)");
  const string udir    = args.getpath("-d",  "--datadir", "data/", "flowfield series directory");
  const string ulabel  = args.getstr ("-ul", "--ulabel",  "u",  "flowfield filename label");
  //const bool turbstats = args.getflag("-t",  "--turbstats",  "compute reynolds stresses");

  //const string Uname   = args.getstr("-U", "--Ubase", "", "base flow, one of [zero|linear|parabolic|<filename>]");
  //const bool stddev    = args.getflag("-sd",  "--stddev",  "compute std dev <|u-<u>|^2>^(1/2)");
  /***/ int digits     = args.getint ("-dg", "--digits",  8,  "# digits in output");
  const bool normalize  = args.getflag("-nrm", "--normalize", "normalize basis elements");
  const string ijklname = args.getstr(2, "<ijklfile>", "filename for ijkl indices of basis functions");
  const string outlabel = args.getstr(1, "<outlabel>",  "filename stub for outputs e.g. foo, produces xfoo.nc and tfoo.asc");
  
  args.check();
  args.save("./");

  const char s = ' ';

  if (adjustdT) {
    // Use timestep dt simply to adjust dT to evenly fit (T1-T0)
    TimeStep dt(0.01, 0.001, 0.1, dT, 0, 1, false);
    dt.adjust_for_T(T1-T0, true);
    cout << setprecision(16);
    cout << "Adjusting time steps to evenly divide T1-T0:\n";
    cout << "dt == " << dt.dt() << endl;
    cout << "dT == " << dt.dT() << endl;
    cout << setprecision(6);    
    dT = dt.dT();
  }
  const bool inttime = (abs(dT - int(dT)) < 1e-12) ? true : false;
  digits = lesser(abs(digits), 17);
  const int width = digits + 6;

  FlowField u0(udir + ulabel + t2s(T0, inttime));
  const int Nx = u0.Nx();
  const int Ny = u0.Ny();
  const int Nz = u0.Nz();
  const Real Lx = u0.Lx();
  const Real Lz = u0.Lz();
  const Real ua = u0.a();           // should check these are -1,1 
  const Real ub = u0.b();
  const Real alpha = 2*pi/Lx;
  const Real gamma = 2*pi/Lz;

  cout << setprecision(16);
  cout << "alpha, gamma == " << alpha << ", " << gamma << endl;
  cout << "Nx, Ny, Nz == " << Nx << ", " << Ny << ", "<< Nz << endl;

  cout << "Reading ijkl indices of basis set from file" << endl;
  MatrixXi ijkl = readijkl(ijklname);
  cout << "ijkl[0] == " << ijkl(0,0) << ' ' << ijkl(0,1) << ' ' << ijkl(0,2) << ' ' << ijkl(0,3) << endl;

  // Figure out the largest value of L in index set
  int L=0;
  for (int n=0; n<ijkl.rows(); ++n)
      L = ijkl(n,3) > L ? ijkl(n,3) : L;
	
  cout << "L == max l == " << L << endl;
	
  cout << "Constructing Legendre polynomials" << endl;
  P.resize(L+1);
  P[0] = Polynomial(1.0);
  if (L >= 1)
      P[1] = Polynomial(0.0, 1.0);
  for (int l = 2; l<=L; ++l) 
      P[l] = ((2*l-1.0)/l)*P[1]*P[l-1] - ((l-1.0)/l)*P[l-2];

  //cout << "Print Legendre polynomials ..." << endl;
  //for (int l = 0; l<=L; ++l)
  //cout << "P[" << l << "] == " << P[l] << endl;

  cout << "Constructing S, and Sprime polynomials" << endl;

  Polynomial smasher(1.0, 0.0, -1.0);
  Polynomial smasher2 = smasher*smasher;
  Polynomial zeropoly = 0.0*smasher;

  S.resize(L+1);
  Sprime.resize(L+1);

  S[0] = Polynomial(0.0, 1.0, 0.0, -1.0/3.0);  // S0(y) = y - y^3/3
  //cout << "2S[0] == " << 2*S[0] << endl;
  for (int l=1; l<=L; ++l) {
      S[l] = P[l-1]*smasher2;
      //cout << "2S[" << l << "] == " << 2*S[l] << endl;
  }
  for (int l=0; l<=L; ++l) {
      Sprime[l] = diff(S[l]);
      //cout << "2Sprime[" << l << "] == " << 2*Sprime[l] << endl;
  }

  Vector xgrid = periodicpoints(Nx, Lx);
  Vector ygrid = chebypoints(Ny, ua, ub);	
  Vector zgrid = periodicpoints(Nz, Lz);
	
  // build H-symmetric fundamental basis functions Psi(i,j,k,l) using i,j,k,l, indices
  // (generated in Julia) from a file. 

  // A wasteful way to generate / evaluate basis set as FlowFields. Could
  // evaluate on low-order grid then expand with FlowField::resize. Will see
  // how long this takes for 50 basis functions

  int Npsi = ijkl.rows();
  cfarray<FlowField> psi(Npsi);

  // Evaluate Psi(i,j,k,l,m, x,y,z) over grid and assign to flowfield psi[n]
  // where ijkl[n,:] = [i,j,k,l]
	
  cout << "Generating psi[i]..." << flush;
  for (int n=0; n<Npsi; ++n) {
      psi[n] = FlowField(Nx,Ny,Nz,3,Lx,Lz, ua, ub);
      psi[n].setState(Physical, Physical);
	    
      int i = ijkl(n,0);
      int j = ijkl(n,1);
      int k = ijkl(n,2);
      int l = ijkl(n,3);
      //cout << "psi(" << i << ", " << j << ", " << k << ", " << l << ")" << endl;
      // m == 0,1,2 is the index of the components of psi vector 
      for (int m=0; m<3; ++m) {
	  for (int nx=0; nx<Nx; ++nx) {
	      Real x = xgrid[nx];
	      for (int nz=0; nz<Nz; ++nz) {
		  Real z = zgrid[nz];		
		  for (int ny=0; ny<Ny; ++ny) {
		      Real y = ygrid[ny];
		      psi[n](nx,ny,nz,m) = Psi(i,j,k,l, m, x,y,z, alpha, gamma);
		  }
	      }
	  }
      }

      psi[n].makeSpectral();
      if (normalize)
	  psi[n] *= 1.0/L2Norm(psi[n]);

  }

  // INEFFICIENCY: IP matrix should be sparse
  cout << "Computing inner product matrix..." << endl;
  Eigen::MatrixXd IP(Npsi, Npsi);
  for (int m=0; m<Npsi; ++m) {
      cout << m << ' ' << flush;
      int mj = ijkl(m,1); // fourier j,k indices of Psi[m]
      int mk = ijkl(m,2);		    
      for (int n=0; n<m; ++n) {
	  int nj = ijkl(n,1); // fourier j,k indices of Psi[n]
	  int nk = ijkl(n,2);		    

	  // nonmatching Fourier modes are orthogonal. Short-circuit eval those IPs to zero.
	  Real ip = (mj == nj && mk == nk) ? L2IP(psi[m], psi[n]) : 0.0;
	  ip = (fabs(ip) > 1e-15) ? ip : 0.0;
	  IP(m,n) = ip;
	  IP(n,m) = ip;
      }
      IP(m,m) = L2IP(psi[m], psi[m]);
  }
  cout << endl;

  cout << "Computing LU decomp of inner product matrix..." << endl;
  Eigen::PartialPivLU<Eigen::MatrixXd> IPlu =  IP.partialPivLu();
  cout << endl;

  
  ofstream xos(("x" + outlabel + ".asc").c_str());
  ofstream tos(("t" + outlabel + ".asc").c_str());

  xos << "# ";
  for (int n=0; n<argc; ++n)
      xos << argv[n] << s;
  xos << endl;
  xos << setprecision(digits);

  tos << "# ";
  for (int n=0; n<argc; ++n)
      tos << argv[n] << s;
  tos << endl;
  tos << setprecision(digits);
  
  for (Real t=T0; t<=T1+dT/2; t += dT) {
      cout << t << ' ' << flush;
      string uname = ulabel + t2s(t, inttime);
      //cout << uname << s << flush;
      FlowField u(udir + uname);

      Eigen::VectorXd ip_psi_u(Npsi);
      for (int i=0; i<Npsi; ++i)
	  ip_psi_u[i] = L2IP(psi[i], u);

      Eigen::VectorXd x = IPlu.solve(ip_psi_u);

      tos << t << '\n';
      for (int i=0; i<Npsi; ++i)
	  xos << setw(width) << x[i] << (i<Npsi-1 ? '\t' : '\n');
  }
  cout << endl;
}}


#include "polynomial.cpp"
#include "basisfunctions.cpp"


/* project.cpp: project time series of fields onto a given basis
 * channelflow-1.1 PCF-utils
 *
 * Copyright (C) 2001-2007  John F. Gibson
 *
 * gibson@cns.physics.gatech.edu  jfg@member.fsf.org
 *
 * Center for Nonlinear Science
 * School of Physics
 * Georgia Institute of Technology
 * Atlanta, GA 30332-0430
 * 404 385 2509
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation version 2
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */
