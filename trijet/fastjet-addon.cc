// Jet algorithm with special features:
//
#include <iostream>
#include <sstream>  // needed for internal io
#include <iomanip>  
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "splitandfilt.hh"
using namespace fastjet;
using namespace std;

//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
// forward declaration for ktscales:
//----------------------------------------------------------------------
void ktscales(const PseudoJet &,const double &,double &,double &);


//----------------------------------------------------------------------
// forward declaration for nsubjettiness:
//----------------------------------------------------------------------
void nsubjettiness(const PseudoJet &, const double &, const int &, double &);

void PressEnterToContinue()
  {
        std::cout << "Press ENTER to continue... " << flush;
	std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
  }


extern "C" { // !extern

void fastjetalg_(const double *p, const int &npart, const double &R, 
                 const double &jetalg, 
		 const bool &filt, const bool &ktsplit, const bool &subjet,
		 const double &ptmin, const double &ymax, 
		 double *pjets, int &njets, 
		 double *sd12arr, double *sd23arr, 
		 double *tau21arr,double *tau32arr)
{
  signed itrack;
  double sd12,sd23;
// We fill up the tracks:
  vector<PseudoJet> particles;
  for (itrack = 0; itrack < npart; itrack++){
    valarray<double> mom(4);
    for (int mu = 0; mu < 4; mu++){
      mom[mu] = *(p++);
    }
//  cout << " " << mom[0] << " " << mom[1] << " " << mom[2] << " " << mom[3] << endl;
  particles.push_back(PseudoJet(mom));
  }

// We have to select a jet algorithm according to the value of jetalg: 
// jetalg = 1 -> kt
// jetalg = 0 -> CA
// jetalg = -1 -> anti-kt
  JetDefinition jet_def(genkt_algorithm,R,jetalg);

// We perform the clustering:
  ClusterSequence cs(particles,jet_def);
// Sorting the jets according to their pt:
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
// We setup our cuts (pt and rapidity):
  Selector select_rapidity = SelectorAbsRapMax(ymax);
  Selector select_pt       = SelectorPtMin(ptmin);
  Selector select_both     = select_pt && select_rapidity;

// We apply the set of cuts:
  vector<PseudoJet> cut_jets = select_both(jets);
// Counting down the surviving jets:
  njets = cut_jets.size();

// We define the filtering routine:
  splitandfilt filter(0.67,0.09,0.3);
// Temporary jet is declared:
  PseudoJet jet_tmp;
// The filtered jets are stored elsewhere:
  vector<PseudoJet> filtered_jets;
// The number of jets can also change:
  int nfiltjets = 0;

  if (njets == 0) return;
  if (filt){
// If filtering is requested we loop over the jets:
//    cout << "Filtering is on..." << endl;
//    cout << "Number of jets is: " << njets << endl;
    for (int i =0; i < njets; i++){
// We filter each and every jet:
//      cout << "The jet under filtering is: " << cut_jets[i] << endl;
      jet_tmp = filter(cut_jets[i]);
// If the jet had substructure and suitable for filtering it is
// added to filtered jets:
      if (jet_tmp != 0){
	nfiltjets += 1;
	filtered_jets.push_back(jet_tmp);
      }
    }
// The number of jets is made equal to the number of filtered jets:
    njets = nfiltjets;
//    cout << "Number of filtered jets is: " << njets << endl;
  }

// transfer jets -> pjets[4*ijet+0..3]
  for (int i=0; i<njets; i++) {
    for (int j=0;j<=3; j++) {
      if (!filt){
        *pjets = cut_jets[i][j];
      }
      else{
        *pjets = filtered_jets[i][j];
      }
      pjets++;
    }
  }

// We calculate the kt splitting scales if needed:
  if (ktsplit){
    for (int i = 0; i < njets; i++){
// Only when jets are not filtered!
      if (!filt){
	ktscales(cut_jets[i],R,sd12,sd23);
	*sd12arr++ = sd12;
	*sd23arr++ = sd23;
      }
      else{
	ktscales(filtered_jets[i],R,sd12,sd23);
	*sd12arr++ = sd12;
	*sd23arr++ = sd23;
      }
    }
  }

// We calculate the subjettiness if needed:
  if (subjet){
    double tau1,tau2,tau3;
    for (int i = 0; i < njets; i++){
      if (!filt){
	nsubjettiness(cut_jets[i],R,1,tau1);
	nsubjettiness(cut_jets[i],R,2,tau2);
	nsubjettiness(cut_jets[i],R,3,tau3);
	if ((tau1 > 0) && ( tau2 > 0)){
	  *tau21arr++ = tau2/tau1;
	}
	else{
	  *tau21arr++ = 0.0;
        }
	if ((tau2 > 0) && ( tau3 > 0)){
	  *tau32arr++ = tau3/tau2;
	}
	else{
	  *tau32arr++ = 0.0;
	}
      }
      else{
	nsubjettiness(filtered_jets[i],R,1,tau1);
	nsubjettiness(filtered_jets[i],R,2,tau2);
	nsubjettiness(filtered_jets[i],R,3,tau3);
	if ((tau1 > 0) && ( tau2 > 0)){
	  *tau21arr++ = tau2/tau1;
	}
	else{
	  *tau21arr++ = 0.0;
        }
	if ((tau2 > 0) && ( tau3 > 0)){
	  *tau32arr++ = tau3/tau2;
	}
	else{
	  *tau32arr++ = 0.0;
	}
      }
    }
  }

//  PressEnterToContinue();

  return;
}

void fastjetppgenkt_pty_(const double * p, const int & npart,                   
                     const double & R, const double & palg, 
		     const double & ptmin, const double & maxrap,
                     double * f77jets, int & njets, int * f77jetvec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
    
    // prepare jet def and run fastjet
    JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = JetDefinition(kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = JetDefinition(cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = JetDefinition(antikt_algorithm, R);
    } else {
      jet_def = JetDefinition(genkt_algorithm, R, palg);
    }

    
    // perform clustering
    ClusterSequence cs(input_particles, jet_def);
    // We have an additional rapidity cut as well:
    Selector select_rap = SelectorAbsRapMax(maxrap);
    Selector select_pt = SelectorPtMin(ptmin);
    Selector select_both = select_rap && select_pt;
    vector<PseudoJet> jets = select_both(cs.inclusive_jets());
    jets = sorted_by_pt(jets);
    njets = jets.size();

    // find particles inside i-th jet
    vector<PseudoJet> *constit;
    constit=new vector<PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }



    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    // clean up
    delete [] constit;
    
   }

void fastjetppgenkt_pteta_(const double * p, const int & npart,                   
                     const double & R, const double & palg, 
		     const double & ptmin, const double & maxeta,
                     double * f77jets, int & njets, int * f77jetvec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
    
    // prepare jet def and run fastjet
    JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = JetDefinition(kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = JetDefinition(cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = JetDefinition(antikt_algorithm, R);
    } else {
      jet_def = JetDefinition(genkt_algorithm, R, palg);
    }

    
    // perform clustering
    ClusterSequence cs(input_particles, jet_def);
    // We have an additional rapidity cut as well:
    Selector select_eta = SelectorAbsEtaMax(maxeta);
    Selector select_pt = SelectorPtMin(ptmin);
    Selector select_both = select_eta && select_pt;
    vector<PseudoJet> jets = select_both(cs.inclusive_jets());
    jets = sorted_by_pt(jets);
    njets = jets.size();

    // find particles inside i-th jet
    vector<PseudoJet> *constit;
    constit=new vector<PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }



    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    // clean up
    delete [] constit;
    
   }

} // !extern

// This routine calculates the kt splitting scales: \sqrt{d_{12}} and
// \sqrt{d_{23}}:
void ktscales(const PseudoJet &jet,const double &R,double &sd12,double &sd23){

// If the jet does not have any constituent we return zeros:
  if (!jet.has_constituents()){
    sd12 = -1.0;
    sd23 = -1.0;
    return;
  }

// We extract all the constituents:
  vector<PseudoJet> particles = jet.constituents();

// Setting up the jet algorithm:
  JetDefinition jet_def(kt_algorithm, R);

// We perform the clustering:
  ClusterSequence cs(particles,jet_def);

// It is possible that a jet contains less than 3 particles in
// this case sd23 cannot be calculated, when the jet only contains
// 1 particle not even sd12 is calculable:
  if (particles.size() < 2){
    sd12 = -1.0;
    sd23 = -1.0;
    return;
  }

// We extract two jets:
  vector<PseudoJet> subjet12 = cs.exclusive_jets(2);
//  cout << "The content of subjet12: " << endl;
//  cout << "subjet no. 1: " << subjet12[0] << endl;
//  cout << "subjet no. 2: " << subjet12[1] << endl;
  sd12 = min(subjet12[0].perp2(),subjet12[1].perp2()) *
	 subjet12[0].squared_distance(subjet12[1]);
// So far this is the square of the scale:
  sd12 = sqrt(sd12);

// When we have less than 3 particles sd23 cannot be calculated:
  if (particles.size() < 3){
    sd23 = -1.0;
    return;
  }

// we extract three jets:
  vector<PseudoJet> subjet123 = cs.exclusive_jets(3);
//  cout << "The content of subjet123: " << endl;
//  cout << "subjet no. 1: " << subjet123[0] << endl;
//  cout << "subjet no. 2: " << subjet123[1] << endl;
//  cout << "subjet no. 3: " << subjet123[2] << endl;
//  A careful look at the subjets revealed that the pseudojets
//  combined in the penultimate step are situating at the first 
//  two positions:
  sd23 = min(subjet123[0].perp2(),subjet123[1].perp2()) *
	 subjet123[0].squared_distance(subjet123[1]);
// So far this is the square of the scale:
  sd23 = sqrt(sd23);
//  cout << "sd12, sd23: " << sd12 << " " << sd23 << endl;

  return;
}

void nsubjettiness(const PseudoJet &jet, const double &R, const int &n, double &taun){

// We have to extract the constituents of the jet:
  vector<PseudoJet> particles = jet.constituents();

// Setting up the jet algorithm:
  JetDefinition jet_def(kt_algorithm, R);

// We perform the clustering:
  ClusterSequence cs(particles,jet_def);

// It is possible that we ask for more subjets than the number of particles
// building up the original one, if this happens we give back a negative taun:
  if (particles.size() < n){
    taun = -10000.0;
    return;
  }

// From the cluster we extract n subjets:
  vector<PseudoJet> subjets = cs.exclusive_jets(n);

// We calculate d0:
  double d0 = 0;
  for (int k = 0; k < particles.size(); k++){
    d0 += R*particles[k].perp();
  }

// To calculate tauN we have to go through all the constituents:
  taun = 0.0;
  for (int k = 0; k < particles.size(); k++){
// We have to determine the minimal distance of the kth constituent from the N subjets:
    double minR = 1000000.0;
    for (int i = 0; i < subjets.size(); i++){
      minR = min(minR,sqrt(particles[k].squared_distance(subjets[i])));
    }
    taun += 1/d0*particles[k].perp()*minR;
  }

  return;
}

//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  ostr << "pt, y, phi =" 
       << " " << setw(10) << jet.perp() 
       << " " << setw(6) <<  jet.rap()  
       << " " << setw(6) <<  jet.phi()  
       << ", mass = " << setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;
}
