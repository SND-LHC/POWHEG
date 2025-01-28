// Only for debug purposes, delete when not needed...
#include <iostream>
#include <sstream>  // needed for internal io
#include <iomanip>  
// <--

#include "splitandfilt.hh"
#include "fastjet/tools/Filter.hh"

using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);

// Short description of the Split and filter routine:
string splitandfilt::description() const{
  ostringstream oss;

  oss << "This is the Split and Filter routine with mu_cut = " << _mu_cut
      << " and v_cut = " << _v_cut;

  return oss.str();

}

PseudoJet splitandfilt::result(const PseudoJet &jet) const{

  PseudoJet j1,j2;
  PseudoJet filtered_jet;

//  cout << "===================================================" << endl;
//  cout << "This is the main Split and Filter routine..." << endl;
//  cout << "We got a jet: " << jet << endl;

  _recurse_through_jet(jet,j1,j2,jet);

// If there is no parent, the jet does not qualify for one which does
// have substructure:
  if (j1 == 0) return PseudoJet();

// We recluster the constituents of the jet using the CA-algorithm:
// We extract the constituents:
  vector<PseudoJet> particles = jet.constituents();
//  cout << "The number of constituents: " << particles.size() << endl;

// We set up the jet algorithm:
  double Rfilt = min(_dRj1j2,sqrt(j1.squared_distance(j2))/2.0);
//  cout << "Rfilt is: " << Rfilt << endl;
  JetDefinition jet_def(cambridge_algorithm,Rfilt);

// We perform the clustering:
  ClusterSequence cs(particles,jet_def);
// We sort the subjets according to their pt:
  vector<PseudoJet> subjets = sorted_by_pt(cs.inclusive_jets());
  int nsubjets = subjets.size();
//  cout << "Number of subjets is: " << nsubjets << endl;
  for (int i = 0; i < min(nsubjets,3); i++){
//    cout << "subjet no. " << i << " has " << subjets[i].constituents().size() << " constituents" << endl;
    if (i == 0) filtered_jet = subjets[i];
    if (i != 0) filtered_jet = join(filtered_jet,subjets[i]);
  }
//  cout << "The number of constituents for the new jet: " << filtered_jet.constituents().size() << endl;
//  cout << "The filtered jet is: " << filtered_jet << endl;

// Give back the new filtered jet:
  return filtered_jet;
    
}

inline void splitandfilt::_recurse_through_jet(const PseudoJet &jet, 
                                               PseudoJet &j1, PseudoJet &j2, 
					       const PseudoJet &original_jet) const{

  PseudoJet parent1 = PseudoJet();
  PseudoJet parent2 = PseudoJet();

// We undo the last clustering step to obtain the parents of jet,
// if there is no one simply return:
//  cout << "This is the recursion..." << endl;
//  cout << "The jet is:          " << jet << endl;
//  cout << "The original jet is: " << original_jet << endl;
//  cout << "Do we have parents: " << jet.has_parents(parent1, parent2) << endl;
  if (! jet.has_parents(parent1, parent2)) return;
// We need one extra condition, there should be a minimal separation 
// between the parents:
  double dRj1j2 = sqrt(parent1.squared_distance(parent2));
//  cout << "The distance between the parents: " << dRj1j2 << endl;
  if (dRj1j2 < _dRj1j2) return;

// The parents should be ordered in mass:
  if (parent1.m() < parent2.m()) swap(parent1,parent2);

// We calculate mu and v:
  double mu = parent1.m() / jet.m();
  double v  = min(parent1.perp2(),parent2.perp2())/jet.m2() *
              parent1.squared_distance(parent2);

//  cout << "mu = " << mu << " v = " << v << endl;
// if mu > _mu_cut and/or v is < _v_cut we have to run the declustering
// again, but with j1:
  if (mu > _mu_cut || v < _v_cut) _recurse_through_jet(parent1,j1,j2,jet);
// Otherwise the two subjets are parent1 and parent2, put them into j1, j2:
  else {
    j1 = parent1;
    j2 = parent2;
//    cout << "The parents are: " << endl;
//    cout << "j1: " << j1 << endl;
//    cout << "j2: " << j2 << endl;
    return;
  }
}
