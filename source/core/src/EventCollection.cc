#include <cmath>
#include <fstream>
#include "EventCollection.hh"

namespace
{
void writeFile(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename);
void makeHistogramLog(double xmin, double xmax, double ratio, std::vector<double>& x_array, std::vector<double>& y_array);
void makeHistogramLinear(double xmin, double xmax, int num_bins, std::vector<double>& x_array, std::vector<double>& y_array);
}

void EventCollection::setMaterialProperties(double r, double h, double vz, double tau, double kTe)
{
  mp_ = std::make_unique<MaterialProperties>();
  mp_ -> setMaterialProperties(r, h, vz, tau, kTe);
}

void EventCollection::readRoot(const std::vector<std::string>& filelist)
{
  etree_ = new TChain("etree");
  for (const std::string& filename: filelist) {
    etree_->Add(filename.c_str());
  }
  
  etree_->SetBranchAddress("energy", &energy_);
  etree_->SetBranchAddress("ini_energy", &iniEnergy_);
  etree_->SetBranchAddress("dirx", &dirx_);
  etree_->SetBranchAddress("diry", &diry_);
  etree_->SetBranchAddress("dirz", &dirz_);
  etree_->SetBranchAddress("ini_dirx", &iniDirx_);
  etree_->SetBranchAddress("ini_diry", &iniDiry_);
  etree_->SetBranchAddress("ini_dirz", &iniDirz_);
  etree_->SetBranchAddress("posx", &posx_);
  etree_->SetBranchAddress("posy", &posy_);
  etree_->SetBranchAddress("posz", &posz_);
}

void EventCollection::calculateWeight()
{
  const int n = etree_->GetEntries();
  weight_.resize(n);

  const double tau = mp_->Tau();
  const double kTe = mp_->kTe();

  for (int i=0; i<n; i++) {
    etree_ -> GetEntry(i);
    const double bremss = (tau*tau) * std::sqrt(kTe) * std::exp(-iniEnergy_/kTe) / iniEnergy_;
    const double powerlaw = std::pow(iniEnergy_, -seedGamma_);
    const double sintheta = std::sqrt(1.0-dirz_*dirz_);
    const double w = (bremss/powerlaw)/sintheta;
    weight_[i] = w;
  }
}

void EventCollection::calculateBeamPattern()
{
  const int n = etree_->GetEntries();
  beamPattern_.resize(n);
  
  const double radius = mp_->Radius();
  const double height = mp_->Height();

  for (int i=0; i<n; i++) {
    etree_->GetEntry(i);
    int id = 0;
    if (posz_<0.5*height && posz_>-0.5*height) {
      if (dirz_>0) {
        const double dz = 0.5*height - posz_;
        const double dx = dz * (dirx_/dirz_);
        const double dy = dz * (diry_/dirz_);
        const double x = posx_ + dx;
        const double y = posy_ + dy;
        const double r2 = x*x + y*y;
        if (r2<radius*radius) {
          id = 1<<0; //pencil
        }
        else {
          id = 1<<1; //fan
        }
      }
      else {
        const double dz = -0.5*height-posz_;
        const double dx = dz * (dirx_/dirz_);
        const double dy = dz * (diry_/dirz_);
        const double x = posx_ + dx;
        const double y = posy_ + dy;
        const double r2 = x*x + y*y;
        if (r2>radius*radius) {
          id = 1<<1; //fan
        }
      }
    }
    else {
      if (dirz_>0) {
        id |= 1<<2; //reflection
      }
    }
    beamPattern_[i] = id;
  }
}

void EventCollection::generateSpectrum(double emin, double emax, double ratio, int bp, const std::string& filename)
{
  std::vector<double> energy_array;
  std::vector<double> photons_array;
  makeHistogramLog(emin, emax, ratio, energy_array, photons_array);
  const int n = etree_->GetEntries();
  for (int i=0; i<n; i++) {
    etree_->GetEntry(i);
    if ((bp & beamPattern_[i])==0) {
      continue;
    }
    auto it = std::lower_bound(energy_array.begin(), energy_array.end(), energy_);
    if (it==energy_array.begin() || it==energy_array.end()) {
      continue;
    }
    const int index = it - energy_array.begin() - 1;
    photons_array[index] += weight_[i];
  }
  writeFile(energy_array, photons_array, filename);
}

void EventCollection::generateEmissionProfile(double phimin, double phimax, int num_bins, double emin, double emax, int bp, const std::string& filename)
{
  std::vector<double> angle_array;
  std::vector<double> photons_array;
  makeHistogramLinear(phimin, phimax, num_bins, angle_array, photons_array);
  const int n = etree_->GetEntries();
  for (int i=0; i<n; i++) {
    etree_->GetEntry(i);
    if ((bp & beamPattern_[i])==0) {
      continue;
    }
    if (energy_<emin || energy_>emax) {
      continue;
    }
    double angle = std::acos(dirz_)*180.0/M_PI;
    for (int j=0; j<2; j++) {
      auto it = std::lower_bound(angle_array.begin(), angle_array.end(), angle);
      if (it==angle_array.begin() || it==angle_array.end()) {
        continue;
      }
      const int index = it - angle_array.begin() - 1;
      photons_array[index] += weight_[i];
      angle *= -1.0;
    }
  }
  writeFile(angle_array, photons_array, filename);  
}

void EventCollection::generateEmissionProfileWithEnergyRanges(double phimin, double phimax, int num_bins, std::vector<double>& energy_array, int bp, const std::string& filename_base)
{
  const int num_energy_bins = energy_array.size() - 1;
  std::vector<double> angle_array;
  std::vector<std::vector<double>> photons_array(num_energy_bins);

  for (int i=0; i<num_energy_bins; i++) {
    makeHistogramLinear(phimin, phimax, num_bins, angle_array, photons_array[i]);
  }
  
  const int n = etree_->GetEntries();
  for (int i=0; i<n; i++) {
    etree_->GetEntry(i);
    if ((bp & beamPattern_[i])==0) {
      continue;
    }
    auto it1 = std::lower_bound(energy_array.begin(), energy_array.end(), energy_);
    if (it1==energy_array.begin() || it1==energy_array.end()) {
      continue;
    }
    const int energy_index = it1 - energy_array.begin() - 1;

    double angle = std::acos(dirz_)*180.0/M_PI;
    for (int j=0; j<2; j++) {
      auto it2 = std::lower_bound(angle_array.begin(), angle_array.end(), angle);
      if (it2==angle_array.begin() || it2==angle_array.end()) {
        continue;
      }
      const int angle_index = it2 - angle_array.begin() - 1;
      photons_array[energy_index][angle_index] += weight_[i];
      angle *= -1.0;
    }
  }

  for (int i=0; i<num_energy_bins; i++) {
    const int emin = energy_array[i];
    const int emax = energy_array[i+1];
    std::string filename = filename_base + "_" + std::to_string(emin) + "-" + std::to_string(emax) + "keV.txt";
    writeFile(angle_array, photons_array[i], filename);
  }
}


namespace
{
void makeHistogramLinear(double xmin, double xmax, int num_bins, std::vector<double>& x_array, std::vector<double>& y_array)
{
  x_array.resize(num_bins+1);
  x_array[0] = xmin;
  for (int i=1; i<=num_bins; i++) {
    x_array[i] = xmin + (xmax-xmin)*i/num_bins;
  }
  y_array.resize(num_bins);
}

void makeHistogramLog(double xmin, double xmax, double ratio, std::vector<double>& x_array, std::vector<double>& y_array)
{
  x_array.resize(1);
  x_array[0] = xmin;
  while (x_array.back()<xmax) {
    x_array.push_back(x_array.back()*ratio);
  }
  y_array.resize(x_array.size()-1, 0.0);
}

void writeFile(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename)
{
  std::ofstream fout(filename);
  const int nx = x.size();
  const int ny = y.size();

  if (nx==ny) {
    for (int i=0; i<ny; i++) {
      fout << x[i] << " " << y[i] << std::endl;
    }
    fout.close();
    return;
  }
  
  if (nx==ny+1) {
    std::vector<double> new_x(ny);
    for (int i=0; i<ny; i++) {
      new_x[i] = 0.5 * (x[i]+x[i+1]);
    }
    for (int i=0; i<ny; i++) {
      fout << new_x[i] << " " << y[i] << std::endl;
    }
    fout.close();
    return;
  }
  std::cerr << "ERROR: different size for x and y" << std::endl;
}
}





