#include <cmath>
#include <fstream>
#include <algorithm>
#include "EmissionProfile.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "TGraphErrors.h"
#include "TF1.h"

namespace
{
double average(const std::vector<double>& v, const std::vector<bool>& valid);
double stddev(const std::vector<double>& v, const std::vector<bool>& valid);
}

void EmissionProfile::clear()
{
  pp_ = nullptr;
  epx_.clear();
  epy_.clear();
  angleArray_.clear();
  modelY_.clear();
  epycomp_.clear();
  modelYcomp_.clear();
}

void EmissionProfile::readPulseProfile(const std::string& filename)
{
  pp_ = std::make_unique<PulseProfile>();
  int cx = 0;
  int cy = 8;
  int cyerr = 9;
  int header = 3;
  pp_->readQDP(filename, cx, cy, cyerr, header);
  std::vector<double>& x = pp_->Phase();
  const int n = x.size();
  phaseArray_.push_back(x[0]-0.5*(x[1]-x[0]));
  for (int i=0; i<n-1; i++) {
    phaseArray_.push_back(0.5*(x[i]+x[i+1]));
  }
  phaseArray_.push_back(x[n-1]+0.5*(x[n-1]-x[n-2]));
  modelY_.resize(n, 0.0);
}

void EmissionProfile::readEmissionProfile(const std::string& filename)
{
  std::ifstream fin(filename);
  double x, y;
  while (fin >> x >> y) {
    epx_.push_back(x);
    epy_.push_back(y);
  }
  fin.close();
  const int n = epx_.size();
  angleArray_.push_back(epx_[0]-0.5*(epx_[1]-epx_[0]));
  for (int i=0; i<n-1; i++) {
    angleArray_.push_back(0.5*(epx_[i]+epx_[i+1]));
  }
  angleArray_.push_back(epx_[n-1]+0.5*(epx_[n-1]-epx_[n-2]));
}

void EmissionProfile::readEmissionProfileComponents(const std::vector<std::string>& filelist)
{
  for (const std::string& filename: filelist) {
    epycomp_.push_back(std::vector<double>());
    std::ifstream fin(filename);
    double x, y;
    while (fin >> x >> y) {
      if (y>1E18) y = 0.0;
      epycomp_.back().push_back(y);
    }
    fin.close();
  }
  const int n = pp_->Phase().size();
  modelYcomp_.resize(filelist.size(), std::vector<double>(n, 0.0));
}

void EmissionProfile::filterReflection(int index, double ang_min, double ang_max)
{
  const int n = epy_.size();
  double selected = 0.0;
  double all = 0.0;
  for (int i=0; i<n; i++) {
    all += epycomp_[index][i];
    if (epx_[i]>ang_min && epx_[i]<ang_max) {
      selected += epycomp_[index][i];
    }
    else {
      epy_[i] -= epycomp_[index][i];
      epycomp_[index][i] = 0.0;
    }
  }
  const double ratio = all/selected;
  
  for (int i=0; i<n; i++) {
    if (epx_[i]>ang_min && epx_[i]<ang_max) {
      double increase = epycomp_[index][i]*(ratio-1.0);
      epycomp_[index][i] *= ratio;
      epy_[i] += increase;
    }
  }
}

void EmissionProfile::fit(double alpha, double beta, double phi_min, double phi_max, double& a, double& b, double& score, bool fixed_offset/*=false*/)
{
  using std::cos;
  using std::sin;
  using std::sqrt;
  using std::acos;
  const double twopi = 2.0*M_PI;
  
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;

  const G4ThreeVector los(sin(beta), 0.0, cos(beta));

  const int n = pp_->Data().size();
  std::vector<bool> valid(n, true);

  for (int i=0; i<n; i++) {
    const double p = pp_->Phase()[i];
    if (p<phi_min || p>phi_max) {
      valid[i] = false;
    }
  }
  
  for (int i=0; i<n; i++) {
    const double phi0 = (phaseArray_[i]-0.26) * twopi;
    const double phi1 = (phaseArray_[i+1]-0.26) * twopi;
    for (int j=0; j<numSample_; j++) {
      const double r = G4UniformRand();
      const double phi = phi0 + (phi1-phi0)*r;
      const G4ThreeVector col_dir(sin(alpha)*cos(phi), sin(alpha)*sin(phi), cos(alpha));
      const double dot = los.dot(col_dir);
      double ang = acos(dot)*180.0/M_PI;
      auto it = std::lower_bound(angleArray_.begin(), angleArray_.end(), ang);
      if (it==angleArray_.begin() || it==angleArray_.end()) {
        continue;
      }
      const int index = it - angleArray_.begin() - 1;
      modelY_[i] += epy_[index];
      for (int k=0; k<static_cast<int>(epycomp_.size()); k++) {
        modelYcomp_[k][i] += epycomp_[k][index];
      }
    }
  }

  
  for (int i=0; i<static_cast<int>(modelYcomp_.size()); i++) {
    for (int j=0; j<static_cast<int>(modelYcomp_[i].size()); j++) {
      if (modelY_[j]==0.0) {
        modelYcomp_[i][j] = 0.0;
        continue;
      }
      modelYcomp_[i][j] /= modelY_[j];
    }
  }
  
  convertPulseProfile(modelY_, valid, a, b, fixed_offset);
  
  for (int i=0; i<static_cast<int>(modelYcomp_.size()); i++) {
    for (int j=0; j<static_cast<int>(modelYcomp_[i].size()); j++) {
      modelYcomp_[i][j] *= modelY_[j];
    }
  }
  
  score = calculateChiSquare(modelY_, valid);
}

void EmissionProfile::convertPulseProfile(std::vector<double>& model_y, const std::vector<bool>& valid, double& a, double& b, bool fixed_offset/*=false*/)
{
  const int n = model_y.size();
  int num = 0;
  for (int i=0; i<n; i++) {
    if (valid[i]) num++;
  }

  TGraphErrors* gr = new TGraphErrors(num);
  int index = 0;
  for (int i=0; i<n; i++) {
    if (!valid[i]) continue;
    gr->SetPoint(index, model_y[i], pp_->Data()[i]);
    gr->SetPointError(index, 0.0, pp_->Error()[i]);
    index++;
  }
  TF1* f = new TF1("f1", "[0]*x+[1]");
  f->SetParameter(0, 1E-9);
  f->SetParLimits(0, 1E-10, 1E-7);
  f->SetParameter(1, 0.0);
  if (fixed_offset) {
    f->SetParLimits(1, 0.0, 1E-15);
  }
  gr->Fit("f1");
  
  a = f->GetParameter(0);
  b = f->GetParameter(1);
  
  for (int i=0; i<n; i++) {
    model_y[i] = a*model_y[i]+b;
  }
}

double EmissionProfile::calculateChiSquare(const std::vector<double>& my, const std::vector<bool>& valid)
{
  double cs = 0.0;
  const int n = my.size();
  for (int i=0; i<n; i++) {
    if (!valid[i]) continue;
    const double chi = (my[i]-pp_->Data()[i])/pp_->Error()[i];
    const double chi2 = chi*chi;
    cs += chi2;
  }
  return cs;
}

void EmissionProfile::writeFile(const std::string& filename)
{
  std::vector<double>& data_x = pp_->Phase();
  std::vector<double>& data_y = pp_->Data();
  std::vector<double>& data_yerr = pp_->Error();
  std::vector<double>& model_y = ModelY();
  
  std::ofstream fout(filename);
  const int n = data_x.size();

  if (n!=static_cast<int>(data_y.size()) || n!=static_cast<int>(data_yerr.size()) | n!=static_cast<int>(model_y.size())) {
    std::cerr << "ERROR: different size for x and y" << std::endl;
    return;
  }
  for (int i=0; i<n; i++) {
    fout << data_x[i] << " " << data_y[i] << " " << data_yerr[i] << " " << model_y[i];
    for (int j=0; j<static_cast<int>(modelYcomp_.size()); j++) {
      fout << " " << modelYcomp_[j][i];
    }
    fout << "\n";
  }

  fout.close(); 
}

namespace
{
double average(const std::vector<double>& v, const std::vector<bool>& valid)
{
  const int n = v.size();
  double sum = 0.0;
  int num = 0;
  for (int i=0; i<n; i++) {
    if (!valid[i]) continue;
    sum += v[i];
    num++;
  }
  if (num>0) {
    sum /= num;
  }
  return sum;
}

double stddev(const std::vector<double>& v, const std::vector<bool>& valid)
{
  const int n = v.size();
  double sum = 0.0;
  double sum2 = 0.0;
  int num = 0;
  for (int i=0; i<n; i++) {
    if (!valid[i]) continue;
    sum += v[i];
    sum2 += v[i]*v[i];
    num++;
  }
  if (num==0) {
    return 0.0;
  }
  const double ret = std::sqrt(sum2/num - (sum/num)*(sum/num));
  return ret;
}

}





