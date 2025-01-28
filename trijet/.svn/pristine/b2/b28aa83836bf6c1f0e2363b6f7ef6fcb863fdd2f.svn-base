#include <iostream>
#include <iomanip>
#include <complex>
#include <stdlib.h>
#include <nlo-innerprod.h>
#include <nlo-event.h>
#include <hep-lorentzvector.h>

#include <ampg5.h>
#include <ampq2g3.h>
#include <ampq4g1.h>

#include <ampg6.h>
#include <ampq2g4.h>
#include <ampq4g2.h>
#include <ampq6.h>

#include "random.h"
#include <nlo-color.h>

using namespace std;
using namespace nlo;

//extern const unsigned Nf = 5U;

extern "C" {




// Wrapper to the 3-jet production Born contribution of NLOJET++
// INPUT
//     subproc: the specific subprocess
//     nleg   : the number of external legs (i.e. nlegborn)
//     pmom   : the array containing the physical momenta
// OUTPUT
//     out    : the Born contribution
  void nlojet3jetborn_(const int &subproc, const int &nleg, double *pmom, 
                       double &out) {

// mandatory declarations:
// lorentz vector with double kind:
    typedef lorentzvector<double> _Lv;
// event datatype to hold momenta:
    typedef hadronic_event_traits<0U,2U,0U> _EventTraits;
    hadronic_event<_Lv,_EventTraits> p(3);
// Innerproducts (A,B and S):
    innerprod<_Lv> _M_ip;
// Random number generator, not needed really:
    rng_mt19937 rng((unsigned long int) rand());
//  pi:
    double pi = 3.14159265358979323846;
// Nf:
//    unsigned int Nf = 5U;
//    unsigned int Nf = nf;
// flavor:
//    double nf = Nf, nf1 = Nf-1;

    double sub1[3];

// We fill up the event:
    for (int ileg = 0; ileg < nleg; ileg++) {
//      cout << "ileg: " << ileg << endl;
      p[ileg-1] = _Lv(pmom[1 + ileg*5],pmom[2 + ileg*5],pmom[3 + ileg*5],pmom[0 + ileg*5]);
//      for (int mu = 0; mu < 5; mu++) {
//	cout << "p" << mu << ": " << pmom[mu + ileg*5] << endl;
//      }
    }
// Calculate inner products:
    _M_ip.calculate(p);
// We select the subprocess:
// gg -> ggg:
   if (subproc == 0) {
//     cout << "Calculating gg -> ggg" << endl;
     ampg5 *_M_g5;
     _M_g5 = new ampg5(_M_ip,rng);
     out = (_M_g5 -> su3_tree(-1,0,1,2,3))/(6.0*Na2);
     out /= 4.;
     delete _M_g5;
     return;
   }
// qqb -> ggg:
   else if (subproc == 1) {
//     cout << "Calculating qqb -> ggg" << endl;
     ampq2g3 *_M_q2g3;
     _M_q2g3 = new ampq2g3(_M_ip,rng);
     out = (_M_q2g3 -> su3_tree(0,-1,1,2,3))/(6.0*Nc2);
     out /= 4.;
     delete _M_q2g3;
     return;
   }
// qqb -> qqbg:
   else if (subproc == 2) {
//     cout << "Calculating qqb -> qqbg" << endl;
     ampq4g1 *_M_q4g1;
     _M_q4g1 = new ampq4g1(_M_ip,rng);
     _M_q4g1 -> su3_tree(0,-1,1,2,3,sub1);
     out = sub1[1]/Nc2;
     out /= 4.;
     delete _M_q4g1;
     return;
   }
// qqb -> rrbg:
   else if (subproc == 3) {
//     cout << "Calculating qqb -> rrbg" << endl;
     ampq4g1 *_M_q4g1;
     _M_q4g1 = new ampq4g1(_M_ip,rng);
     _M_q4g1 -> su3_tree(0,-1,1,2,3,sub1);
     out = sub1[0]/Nc2;
     out /= 4.;
     delete _M_q4g1;
     return;
   }
  }






// Wrapper to the 3-jet production virtual contribution of NLOJET++
// INPUT
//     subproc: the specific subprocess
//     nleg   : the number of external legs (i.e. nlegborn)
//     pmom   : the array containing the physical momenta
//     mur    : the renormalization scale
//     nf     : the number of light flavors
// OUTPUT
//     out    : the virtual contribution
  void nlojet3jetvirt_(const int &subproc, const int &nleg, double *pmom, 
                       const double &mur, const int &nf,
                       double &out) {

// mandatory declarations:
// lorentz vector with double kind:
    typedef lorentzvector<double> _Lv;
// event datatype to hold momenta:
    typedef hadronic_event_traits<0U,2U,0U> _EventTraits;
    hadronic_event<_Lv,_EventTraits> p(3);
// Innerproducts (A,B and S):
    innerprod<_Lv> _M_ip;
// Random number generator, not needed really:
    rng_mt19937 rng((unsigned long int) rand());
//  pi:
    double pi = 3.14159265358979323846;
// Nf:
//    unsigned int Nf = 5U;
    unsigned int Nf = nf;
// flavor:
//    double nf = Nf, nf1 = Nf-1;

    double sub1[3];

// We fill up the event:
    for (int ileg = 0; ileg < nleg; ileg++) {
//      cout << "ileg: " << ileg << endl;
      p[ileg-1] = _Lv(pmom[1 + ileg*5],pmom[2 + ileg*5],pmom[3 + ileg*5],pmom[0 + ileg*5]);
//      for (int mu = 0; mu < 5; mu++) {
//	cout << "p" << mu << ": " << pmom[mu + ileg*5] << endl;
//      }
    }
// Calculate inner products:
    _M_ip.calculate(p);
// We select the subprocess:
// gg -> ggg:
   if (subproc == 0) {
//     cout << "Calculating gg -> ggg" << endl;
     ampg5 *_M_g5;
     _M_g5 = new ampg5(_M_ip,rng);
// We set up the scale:
     _M_g5->set_scale(mur);
// The full finite term of the virtual part:
     out = (_M_g5 -> su3_1loop(Nf, -1,0,1,2,3))/(6.0*Na2);
     //     outtree = (_M_g5 -> su3_tree(-1,0,1,2,3))/(6.0*Na2);
     out /= 32.*pi*pi;
     //     outtree /= 4.;
     delete _M_g5;
     return;
   }
// qqb -> ggg:
   else if (subproc == 1) {
//     cout << "Calculating qqb -> ggg" << endl;
     ampq2g3 *_M_q2g3;
     _M_q2g3 = new ampq2g3(_M_ip,rng);
// We set up the scale:
     _M_q2g3->set_scale(mur);
// The full finite term of the virtual part:
     out = (_M_q2g3 -> su3_1loop(Nf,0,-1,1,2,3))/(6.0*Nc2);
     //     outtree = (_M_q2g3 -> su3_tree(0,-1,1,2,3))/(6.0*Nc2);
     out /= 32.*pi*pi;
     //     outtree /= 4.;
     delete _M_q2g3;
     return;
   }
// qqb -> qqbg:
   else if (subproc == 2) {
//     cout << "Calculating qqb -> qqbg" << endl;
     ampq4g1 *_M_q4g1;
     _M_q4g1 = new ampq4g1(_M_ip,rng);
// We set up the scale:
     _M_q4g1->set_scale(mur);
// The full finite term of the virtual part:
     _M_q4g1 -> su3_1loop(Nf,0,-1,1,2,3,sub1);
     out = sub1[1]/Nc2;
     //     _M_q4g1 -> su3_tree(0,-1,1,2,3,sub1);
     //     outtree = sub1[1]/Nc2;
     out /= 32.*pi*pi;
     //     outtree /= 4.;
     delete _M_q4g1;
     return;
   }
// qqb -> rrbg:
   else if (subproc == 3) {
//     cout << "Calculating qqb -> rrbg" << endl;
     ampq4g1 *_M_q4g1;
     _M_q4g1 = new ampq4g1(_M_ip,rng);
// We set up the scale:
     _M_q4g1->set_scale(mur);
// The full finite term of the virtual part:
     _M_q4g1 -> su3_1loop(Nf,0,-1,1,2,3,sub1);
     out = sub1[0]/Nc2;
     //     _M_q4g1 -> su3_tree(0,-1,1,2,3,sub1);
     //     outtree = sub1[0]/Nc2;
     out /= 32.*pi*pi;
     //     outtree /= 4.;
     delete _M_q4g1;
     return;
   }
  }


// Wrapper to the 4-jet production real contribution of NLOJET++
// INPUT
//     subproc: the specific subprocess
//     nleg   : the number of external legs (i.e. nlegreal)
//     pmom   : the array containing the physical momenta
// OUTPUT
//     out    : the real contribution
  void nlojet3jetreal_(const int &subproc, const int &nleg, double *pmom, 
                       double &out) {

// mandatory declarations:
// lorentz vector with double kind:
    typedef lorentzvector<double> _Lv;
// event datatype to hold momenta:
    typedef hadronic_event_traits<0U,2U,0U> _EventTraits;
    hadronic_event<_Lv,_EventTraits> p(4);
// Innerproducts (A,B and S):
    innerprod<_Lv> _M_ip;
// Random number generator, not needed really:
    rng_mt19937 rng((unsigned long int) rand());
//  pi:
    double pi = 3.14159265358979323846;
// Nf:
//    unsigned Nf = 5U;
// flavor:
//    double nf = Nf, nf1 = Nf-1;

    double sub1[3], sub2[5];

// We fill up the event:
    for (int ileg = 0; ileg < nleg; ileg++) {
//      cout << "ileg: " << ileg << endl;
      p[ileg-1] = _Lv(pmom[1 + ileg*5],pmom[2 + ileg*5],pmom[3 + ileg*5],pmom[0 + ileg*5]);
//      for (int mu = 0; mu < 5; mu++) {
//	cout << "p" << mu << ": " << pmom[mu + ileg*5] << endl;
//      }
    }
// Calculate inner products:
    _M_ip.calculate(p);
// We select the subprocess:
// gg -> gggg:
   if (subproc == 0) {
     ampg6 *_M_g6;
     _M_g6 = new ampg6(_M_ip,rng);
     out = (_M_g6 -> su3_tree(-1,0,1,2,3,4))/(24.0*Na2);
     out /= 4.;
     delete _M_g6;
     return;
   }
// qqb -> gggg:
   else if (subproc == 1) {
     ampq2g4 *_M_q2g4;
     _M_q2g4 = new ampq2g4(_M_ip,rng);
     out = (_M_q2g4 -> su3_tree(0,-1,1,2,3,4))/(24.0*Nc2);
     out /= 4.;
     delete _M_q2g4;
     return;
   }
// qqb -> qqbgg:
   else if (subproc == 2) {
     ampq4g2 *_M_q4g2;
     _M_q4g2 = new ampq4g2(_M_ip,rng);
     _M_q4g2 -> su3_tree(0,-1,1,2,3,4,sub1);
     out = sub1[1]/(2.0*Nc2);
     out /= 4.;
     delete _M_q4g2;
     return;
   }
// qqb -> rrbgg:
   else if (subproc == 3) {
     ampq4g2 *_M_q4g2;
     _M_q4g2 = new ampq4g2(_M_ip,rng);
     _M_q4g2 -> su3_tree(0,-1,1,2,3,4,sub1);
     out = sub1[0]/(2.0*Nc2);
     out /= 4.;
     delete _M_q4g2;
     return;
   }
// qqb -> qqbqqb:
   else if (subproc == 4) {
     ampq6 *_M_q6;
     _M_q6 = new ampq6(_M_ip,rng);
     _M_q6 -> su3_tree(0,-1,1,2,3,4,"11011",sub2);
     out = sub2[4]/(4.0*Nc2);
     out /= 4.;
     delete _M_q6;
     return;
   }
// qqb -> qqbrrb:
   else if (subproc == 5) {
     ampq6 *_M_q6;
     _M_q6 = new ampq6(_M_ip,rng);
     _M_q6 -> su3_tree(0,-1,1,2,3,4,"11011",sub2);
     out = sub2[1]/(Nc2);
     out /= 4.;
     delete _M_q6;
     return;
   }
// qqb -> rrbssb:
   else if (subproc == 6) {
     ampq6 *_M_q6;
     _M_q6 = new ampq6(_M_ip,rng);
     _M_q6 -> su3_tree(0,-1,1,2,3,4,"11011",sub2);
     out = sub2[0]/(Nc2);
     out /= 4.;
     delete _M_q6;
     return;
   }
// qqb -> rrbrrb:
   else if (subproc == 7) {
     ampq6 *_M_q6;
     _M_q6 = new ampq6(_M_ip,rng);
     _M_q6 -> su3_tree(0,-1,1,2,3,4,"11011",sub2);
     out = sub2[3]/(4.0*Nc2);
     out /= 4.;
     delete _M_q6;
     return;
   }
  }
}
