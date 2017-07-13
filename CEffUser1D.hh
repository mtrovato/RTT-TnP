#ifndef CEFFUSER1D_HH
#define CEFFUSER1D_HH

#include <TGraphAsymmErrors.h>
#include <iostream>

class CEffUser1D
{
public:
  CEffUser1D();
  ~CEffUser1D();
  
  void  loadEff(TGraphAsymmErrors* gr);
  float getEff(const double x);
  float getErrLow(const double x);
  float getErrHigh(const double x);    
  void  printEff(std::ostream& os);
  void  printErrLow(std::ostream& os);
  void  printErrHigh(std::ostream& os);

protected:
  Int_t getBin(const double x);
  void  print(const double *yval, std::ostream& os);  
  
  TGraphAsymmErrors *graph;
};

#endif
