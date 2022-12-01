#include <cmath>
#include <fstream>
#include "MaterialProperties.hh"


void MaterialProperties::setMaterialProperties(double r, double h, double vz, double tau, double kTe)
{
  setRadius(r);
  setHeight(h);
  setBulkVelocityZ(vz);
  setTau(tau);
  setkTe(kTe);
}
