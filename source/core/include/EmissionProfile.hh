#ifndef EmissionProfile_H
#define EmissionProfile_H 1

#include <iostream>
#include <string>
#include <vector>
#include "MaterialProperties.hh"
#include "TChain.h"
class TChain;

class EmissionProfile
{
public:
  void setMaterialProperties(double r, double h, double vz, double tau, double kTe);
  void readRoot(const std::vector<std::string>& filelist);
  void calculateWeight();
  void calculateBeamPattern();
  void generateSpectrum(double emin, double emax, double ratio, int bp, const std::string& filename);
  void generateEmissionProfile(double phimin, double phimax, int num_bins, double emin, double emax, int bp, const std::string& filename);
  void generateEmissionProfileWithEnergyRanges(double phimin, double phimax, int num_bins, std::vector<double>& energy_array, int bp, const std::string& filename_base);
  
  void setSeedGamma_(double v) { seedGamma_ = v; }

  double SeedGamma() { return seedGamma_; }
  std::vector<double>& Weight() { return weight_; }
    
private:
  std::unique_ptr<MaterialProperties> mp_;
  TChain* etree_;
  double seedGamma_ = 1.0;
  std::vector<double> weight_;
  std::vector<int> beamPattern_;
  
  float energy_ = 0.0;
  float iniEnergy_ = 0.0;
  float dirx_ = 0.0;
  float diry_ = 0.0;
  float dirz_ = 0.0;
  float iniDirx_ = 0.0;
  float iniDiry_ = 0.0;
  float iniDirz_ = 0.0;
  float posx_ = 0.0;
  float posy_ = 0.0;
  float posz_ = 0.0;
};

#endif
