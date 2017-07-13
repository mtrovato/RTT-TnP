#include "CEffUser2D.hh"
#include <iomanip>
#include <cassert>

using namespace std;

CEffUser2D::CEffUser2D():hEff(0),hErrl(0),hErrh(0){}
CEffUser2D::~CEffUser2D(){}

//--------------------------------------------------------------------------------------------------
void CEffUser2D::loadEff(TH2D* eff, TH2D* errl, TH2D* errh)
{
  hEff = eff;
  hErrl = errl;
  hErrh = errh;
}

//--------------------------------------------------------------------------------------------------
float CEffUser2D::getEff(const double x, const double y)
{
  if(!hEff) {
    cout << "Efficiency table not loaded! Aborting..." << endl;
    assert(0);
  }
  return getValue(hEff,x,y);
}

//--------------------------------------------------------------------------------------------------
float CEffUser2D::getErrLow(const double x, const double y)
{
  if(!hErrl) {
    cout << "Low errors table not loaded! Aborting..." << endl;
    assert(0);
  }
  return getValue(hErrl,x,y);
}

//--------------------------------------------------------------------------------------------------
float CEffUser2D::getErrHigh(const double x, const double y)
{
  if(!hErrh) {
    cout << "High errors table not loaded! Aborting..." << endl;
    assert(0);
  }
  return getValue(hErrh,x,y);
}

//--------------------------------------------------------------------------------------------------  
void CEffUser2D::printEff(ostream& os)
{
  if(!hEff) {
    cout << "Efficiency table not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "Efficiency Table:" << endl;
  os << "-----------------" << endl;
  printHist2D(hEff,os);   
}

//--------------------------------------------------------------------------------------------------
void CEffUser2D::printErrLow(ostream& os)
{
  if(!hErrl) {
    cout << "Error (low) table not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "Low Errors Table:" << endl;
  os << "-----------------" << endl;
  printHist2D(hErrl,os); 
}

//--------------------------------------------------------------------------------------------------
void CEffUser2D::printErrHigh(ostream& os)
{
  if(!hErrh) {
    cout << "Error (high) table not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "High Errors Table:" << endl;
  os << "------------------" << endl;
  printHist2D(hErrh,os);
}

//--------------------------------------------------------------------------------------------------
void CEffUser2D::printHist2D(const TH2D* h, ostream& os) {
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  
  for(int iy=0; iy<=ny; iy++) {
    for(int ix=0; ix<=nx; ix++) {
      if(ix==0 && iy==0) {
        os << setw(11) << "";
      } else if(ix==0) {
        os << "[" << setw(4) << h->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << h->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      } else if(iy==0) { 
        os << "[" << setw(4) << h->GetXaxis()->GetBinLowEdge(ix) << "," << setw(4) << h->GetXaxis()->GetBinLowEdge(ix+1) << "]";
      } else {
        ios_base::fmtflags flags = os.flags();
	os.precision(7);
	os << " " << setw(9) << fixed << h->GetCellContent(ix,iy) << " ";
	os.flags(flags);
      }
    }
    os << endl;
  }
}

//--------------------------------------------------------------------------------------------------
float CEffUser2D::getValue(const TH2D* h, const double x, const double y)
{
  int ix=0;
  int iy=0;
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  
  for(int i=1; i<=nx; i++) {
    if((x >= h->GetXaxis()->GetBinLowEdge(i)) && (x < h->GetXaxis()->GetBinLowEdge(i+1))) {
      ix=i;
      break;
    }
  }
  
  for(int i=1; i<=ny; i++) {
    if((y >= h->GetYaxis()->GetBinLowEdge(i)) && (y < h->GetYaxis()->GetBinLowEdge(i+1))) {
      iy=i;
      break;
    }
  }
  
  if(ix>0 && iy>0)
    return h->GetCellContent(ix,iy);
  else 
    return -1;
}
