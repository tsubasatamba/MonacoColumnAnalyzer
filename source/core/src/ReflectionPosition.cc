#include <cmath>
#include <fstream>
#include "ReflectionPosition.hh"

namespace
{
void writeFile(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename);
void makeHistogramLog(double xmin, double xmax, double ratio, std::vector<double>& x_array, std::vector<double>& y_array);
void makeHistogramLinear(double xmin, double xmax, int num_bins, std::vector<double>& x_array, std::vector<double>& y_array);
}


void ReflectionPosition::calculatePosition(double phimin, double phimax, int num_bins, double emin, double emax, int bp, const std::string& filename)
{
  TChain* etree = Etree();
  std::vector<double> weight = Weight();
  std::vector<int> beam_pattern = BeamPattern();
  
  std::vector<double> angle_array;
  std::vector<double> rsum_array;
  std::vector<double> wsum_array;
  makeHistogramLinear(phimin, phimax, num_bins, angle_array, rsum_array);
  wsum_array.resize(rsum_array.size());
  const int n = etree->GetEntries();
  
  
  for (int i=0; i<n; i++) {
    etree->GetEntry(i);
    double energy = Energy();
    double posx = Posx();
    double posy = Posy();
    double dirz = Dirz();
    if ((bp & beam_pattern[i])==0) {
      continue;
    }
    if (energy<emin || energy>emax) {
      continue;
    }
    double angle = std::acos(dirz)*180.0/M_PI;
    for (int j=0; j<2; j++) {
      auto it = std::lower_bound(angle_array.begin(), angle_array.end(), angle);
      if (it==angle_array.begin() || it==angle_array.end()) {
        continue;
      }
      const int index = it - angle_array.begin() - 1;
      double r = std::sqrt(posx*posx+posy*posy);
      rsum_array[index] += r*weight[i];
      wsum_array[index] += weight[i];
      angle *= -1.0;
    }
  }

  const int m = rsum_array.size();
  for (int i=0; i<m; i++) {
    if (wsum_array[i]==0.0) continue;
    rsum_array[i] /= wsum_array[i];
  }
  
  writeFile(angle_array, rsum_array, filename);  
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





