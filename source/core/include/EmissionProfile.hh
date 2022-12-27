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
  void readEmissionProfileComponents(const std::vector<std::string>& filelist);
  void filterReflection(int index, double ang_min, double ang_max);
  void fit(double alpha, double beta, double phi_min, double phi_max, double& a, double& b, double& score, bool fixed_offset=false);
  void convertPulseProfile(std::vector<double>& model_y, const std::vector<bool>& valid, double& a, double& b, bool fixed_offset=false);
  double calculateChiSquare(const std::vector<double>& my, const std::vector<bool>& valid);
  void writeFile(const std::string& filename);

  std::vector<double>& ModelY() { return modelY_; }
    
private:
  std::unique_ptr<PulseProfile> pp_ = nullptr;
  std::vector<double> phaseArray_;
  std::vector<double> epx_;
  std::vector<double> epy_;
  std::vector<double> angleArray_;
  int numSample_ = 1000;
  std::vector<double> modelY_;
  std::vector<std::vector<double>> epycomp_;
  std::vector<std::vector<double>> modelYcomp_;
};

#endif
