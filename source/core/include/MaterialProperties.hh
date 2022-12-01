#ifndef MaterialProperties_H
#define MaterialProperties_H 1

#include <iostream>

class MaterialProperties
{
public:
  void setMaterialProperties(double r, double h, double vz, double tau, double kTe);
  
  void setRadius(double v) { radius_ = v; }
  void setHeight(double v) { height_ = v; }
  void setBulkVelocityZ(double v) { bulkVelocityZ_ = v; }
  void setTau(double v) { tau_ = v; }
  void setkTe(double v) { kTe_ = v; }

  double Radius() { return radius_; }
  double Height() { return height_; }
  double BulkVelocityZ() { return bulkVelocityZ_; }
  double Tau() { return tau_; }
  double kTe() { return kTe_; }
    
private:
  double radius_ = 0.0;
  double height_ = 0.0;
  double bulkVelocityZ_ = 0.0;
  double tau_ = 0.0;
  double kTe_ = 0.0;
};

#endif
