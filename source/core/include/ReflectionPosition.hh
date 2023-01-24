#ifndef ReflectionPosition_H
#define ReflectionPosition_H 1

#include <iostream>
#include <string>
#include <vector>
#include "MaterialProperties.hh"
#include "EventCollection.hh"
#include "TChain.h"
class TChain;

class ReflectionPosition : public EventCollection
{
public:
  void calculatePosition(double phimin, double phimax, int num_bins, double emin, double emax, int bp, const std::string& filename);
    
private:

};

#endif
