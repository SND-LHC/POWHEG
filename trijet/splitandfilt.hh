#ifndef __SPLITANDFILT_HH__
#define __SPLITANDFILT_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Transformer.hh"
#include "fastjet/WrappedStructure.hh"

using namespace std;
using namespace fastjet;

class splitandfilt;
class splitandfiltStructure;

class splitandfilt : public Transformer {
public:

  splitandfilt(double mu_cut = 0.67, double v_cut = 0.09, 
               double dRj1j2 = 0.3)
    : _mu_cut(mu_cut), _v_cut(v_cut), _dRj1j2(dRj1j2){};

  virtual PseudoJet result(const PseudoJet & jet) const;

// the type of Structure returned
  typedef splitandfiltStructure StructureType;

  virtual std::string description() const;

protected:

  void _recurse_through_jet(const PseudoJet & current_jet,
                            PseudoJet &j1, PseudoJet &j2,
			    const PseudoJet & original_jet) const;

  double _mu_cut;
  double _v_cut;
  double _dRj1j2;
};

class splitandfiltStructure : public WrappedStructure{
public:
  /// default ctor
  ///  \param result_jet   the jet for which we have to keep the
  ///                      structure
  splitandfiltStructure(const PseudoJet & result_jet)
    : WrappedStructure(result_jet.structure_shared_ptr()){}

protected:
  friend class splitandfilt; ///< to allow setting the internal information
};

#endif // __SPLITANDFILT_HH__ 
