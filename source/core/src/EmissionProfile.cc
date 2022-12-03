#include <cmath>
#include <fstream>
#include <algorithm>
#include "EmissionProfile.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

namespace
{
void writeFile(const std::vector<double>& data_x, const std::vector<double>& data_y, const std::vector<double>& data_yerr, const std::vector<double>& model_y, const std::vector<bool>& valid, const std::string& filename);
double average(const std::vector<double>& v, const std::vector<bool>& valid);
double stddev(const std::vector<double>& v, const std::vector<bool>& valid);
}

void EmissionProfile::clear()
{
  pp_ = nullptr;
  epx_.clear();
  epy_.clear();
  angleArray_.clear();
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
}

void EmissionProfile::readEmissionProfile(const std::string& filename)
{
  std::ifstream fin(filename);
  double x, y;
  while (fin >> x >> y) {
    epx_.push_back(x);
    epy_.push_back(y);
  }
  const int n = epx_.size();
  angleArray_.push_back(epx_[0]-0.5*(epx_[1]-epx_[0]));
  for (int i=0; i<n-1; i++) {
    angleArray_.push_back(0.5*(epx_[i]+epx_[i+1]));
  }
  angleArray_.push_back(epx_[n-1]+0.5*(epx_[n-1]-epx_[n-2]));
}

void EmissionProfile::fit(double alpha, double beta, double phi_min, double phi_max, double& score, const std::string& filename)
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
  std::vector<double> model_y(n, -1E10);
  std::vector<bool> valid(n, true);
  
  for (int i=0; i<n; i++) {
    const double phi0 = (phaseArray_[i]-0.25) * twopi;
    const double phi1 = (phaseArray_[i+1]-0.25) * twopi;
    const double phi_center = 0.5*(phi0+phi1);
    if (phi_center<phi_min || phi_center>phi_max) {
      valid[i] = false;
      continue;
    }
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
      const int index = it - angleArray_.begin() -1;
      model_y[i] += epy_[index];
    }
  }

  convertPulseProfile(model_y, valid);

  score = calculateChiSquare(model_y, valid);

  writeFile(pp_->Phase(), pp_->Data(), pp_->Error(), model_y, valid, filename);
}

void EmissionProfile::convertPulseProfile(std::vector<double>& model_y, const std::vector<bool>& valid)
{
  const double model_mu = average(model_y, valid);
  const double model_sigma = stddev(model_y, valid);
  const double data_mu = average(pp_->Data(), valid);
  const double data_sigma = stddev(pp_->Data(), valid);
  const double a = data_sigma/model_sigma;
  const double b = data_mu/(model_mu*a);
  
  const int n = model_y.size();
  for (int i=0; i<n; i++) {
    if (!valid[i]) {
      continue;
    }
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

namespace
{
void writeFile(const std::vector<double>& data_x, const std::vector<double>& data_y, const std::vector<double>& data_yerr, const std::vector<double>& model_y, const std::vector<bool>& valid, const std::string& filename)
{
  std::ofstream fout(filename);
  const int n = data_x.size();

  if (n!=static_cast<int>(data_y.size()) || n!=static_cast<int>(data_yerr.size()) | n!=static_cast<int>(model_y.size())) {
    std::cerr << "ERROR: different size for x and y" << std::endl;
    return;
  }

  for (int i=0; i<n; i++) {
    if (!valid[i]) continue;
    fout << data_x[i] << " " << data_y[i] << " " << data_yerr[i] << " " << model_y[i] << "\n";
  }
  fout.close(); 
}

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





