#ifndef EmissionProfile_H
#define EmissionProfile_H 1

#include <iostream>
#include <string>
#include <vector>
#include "PulseProfile.hh"


class EmissionProfile
{
public:
  void clear();
  void readPulseProfile(const std::string& filename);
  void readEmissionProfile(const std::string& filename);
  void fit(double alpha, double beta, double phi_min, double phi_max, double& score, const std::string& filename);
  void convertPulseProfile(std::vector<double>& model_y, const std::vector<bool>& valid);
  double calculateChiSquare(const std::vector<double>& my, const std::vector<bool>& valid);
    
private:
  std::unique_ptr<PulseProfile> pp_ = nullptr;
  std::vector<double> phaseArray_;
  std::vector<double> epx_;
  std::vector<double> epy_;
  std::vector<double> angleArray_;
  int numSample_ = 1000;

};

#endif
