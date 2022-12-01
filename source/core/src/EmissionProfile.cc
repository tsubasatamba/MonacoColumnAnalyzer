#include <cmath>
#include <fstream>
#include "EmissionProfile.hh"

namespace
{
void writeFile(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename);
}




namespace
{
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





