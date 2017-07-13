#include "TROOT.h"
#include "TH1D.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooVoigtianShape.h"
#include "RooKeysPdf.h"

class CSignalModel
{
public:
  CSignalModel():model(0){}
  virtual ~CSignalModel(){ delete model; }
  RooAbsPdf *model;
};

// class CBreitWignerConvCrystalBall : public CSignalModel
// {
// public:
//   CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass);
//   ~CBreitWignerConvCrystalBall();
//   RooRealVar     *mass, *width;
//   RooBreitWigner *bw;
//   RooRealVar     *mean, *sigma, *alpha, *n;
//   RooCBShape     *cb;
// };

class CMCTemplateConvGaussian : public CSignalModel
{
public:
  CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const int pass, RooRealVar *sigma0=0, int intOrder=1);
  ~CMCTemplateConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
};


//--------------------------------------------------------------------------------------------------
CMCTemplateConvGaussian::CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const int pass, RooRealVar *sigma0, int intOrder)
{  
  char name[10];
  if(pass==0) sprintf(name,"%s","Nuis1");
  else if(pass==1)    sprintf(name,"%s","Nuis2");
  else if(pass==2)    sprintf(name,"%s","Nuis3");


  
  char vname[50];  
  
  sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2,2); 
  if(sigma0) { sigma = sigma0; }
  else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,4); } 
  sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  //}


  sprintf(vname,"inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);
  
  sprintf(vname,"dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
  sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);
  sprintf(vname,"signal%s",name);   model    = new RooFFTConvPdf(vname,vname,m,*histPdf,*gaus);


}

CMCTemplateConvGaussian::~CMCTemplateConvGaussian()
{
  if(mean){
    delete mean;     
    mean=0;
  }
  if(sigma){
    delete sigma;    
    sigma=0;
  }
  if(gaus){
    delete gaus;     
    gaus=0;
  }
  if(inHist){    
    delete inHist;
    inHist=0;
  }
  if(dataHist){
    delete dataHist; 
    dataHist=0;
  }
  if(histPdf){
    delete histPdf;  
    histPdf=0;
  }
}

