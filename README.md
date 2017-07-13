i# Overview
This code package performs the extraction of efficiencies on a given ROOT ntuple of probes from a tag-and-probe selection. After execution, the efficiencies are stored as `TH2D` or `TGraphAsymErrors` in a ROOT file. The source code for the executable is `effZFit.cc`. The class that does all the heavy lifting is `CEffZFitter`. Signal and background models are found in `ZSignals.hh` and `ZBackgrounds.hh`. Classes such as `CEffUser1D`, `CEffUser2D`, and `CPlot` provide functions for handling the set of efficiency numbers and making plots.


## Inputs
Two major inputs to the executable is a binning config file and a ROOT ntuple of probes.


### Binning config file
The binning config file is a text file that specifies in which variables to bin for the efficiency extraction and the boundaries defining the bins. The file should have 5 sections:

1. Toggle to extract efficiencies in pT-bins, &eta;-bins, &phi;-bins, nPV-bins, (&eta;,pT)-bins, or (&eta;,&phi;)-bins
2. Binning in pT
3. Binning in &eta; (option for |&eta;|)
4. Binning in &phi; (option for |&phi;|)
5. Binning in number of primary vertices

The upper and lower bounds for all variables affect the probes used for efficiency extraction, even if the extraction is not done with respect to some variables. For example, if the limits of &eta; are between -1.2 and +1.2, then only probles satisfying those limits are used even if efficiencies are computed with respect to pT.

The parsing of the binning config file is found in `CEffZFitter.cc`. The parser is pretty basic, so the format of the config file is fairly strict. The '\#' character at the beginning of a line denotes a comment line. The '%' character is the delimiter between sections of the config file.

Examples are provided in the tutorial.


### Probes ntuple data format
The expected data format is,
```
unsigned int runNum, lumiSec, evtNum;   // event ID
unsigned int npv;                       // number of primary vertices
unsigned int pass;                      // whether probe passes requirements
float        npu;                       // mean number of expected pileup
float        scale1fb;                  // event weight per 1/fb
float        mass;                      // tag-probe mass
int          qtag, qprobe;              // tag, probe charge
TLorentzVector *tag=0, *probe=0;        // tag, probe 4-vector
```
stored in a TTree named "Events",
```
TFile *infile = new TFile(infname.c_str());    assert(infile);
TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
intree->SetBranchAddress("runNum",   &runNum);
intree->SetBranchAddress("lumiSec",  &lumiSec);
intree->SetBranchAddress("evtNum",   &evtNum);
intree->SetBranchAddress("npv",      &npv);
intree->SetBranchAddress("pass",     &pass);
intree->SetBranchAddress("npu",      &npu);
intree->SetBranchAddress("scale1fb", &scale1fb);
intree->SetBranchAddress("mass",     &mass);
intree->SetBranchAddress("qtag",     &qtag);
intree->SetBranchAddress("qprobe",   &qprobe);
intree->SetBranchAddress("tag",      &tag);
intree->SetBranchAddress("probe",    &probe);
```

The data format expected by the code is found in `CEffZFitter.cc`.


## Running the code
If the executable has not been built yet, do so by calling,
```
source build.sh
```

Set the `LD_LIBRARY_PATH` environment variable to include the current directory:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
```
this needs to be done every time you open a new shell.

The executable call is:
```
./effZFit  <conf>  <sig-pass>  <bkg-pass>  <sig-fail>  <bkg-fail>  <in-fname>  <outdir>  <doPU>  <charge>  <template-file>
```

where the arguments are:
- **conf:** name of the binning config file
- **sig-pass:** signal model number for passing probes
- **bkg-pass:** background model number for passing probes
- **sig-fail:** signal model number for failing probes
- **bkg-fail:** background model number for failing probes
- **in-fname:** input probes ntuple filename
- **outdir:** output directory name
- **doPU:** do pileup reweighting? (0 or 1)
- **charge:** charge requirement on probe (0=none, +1, -1)
- **template-file:** probes ntuple filename used to build templates (needed for some signal models)

The model numbers corresponding to various signal and background models can be found in `CEffZFitter.hh`.

By default, the fits are done with Strategy=1. Using Strategy=2 and MINOS for more accurate minimization and uncertainty can be done by uncommenting a couple lines in `CEffZFitter.cc`.

*(RooFit may throw a bunch of messages at the end of execution. It seems to be an issue of end program clean up but it is not yet clear what is the problem specifically.)*

### Re-doing fits
Sometimes a fit does not converge and different initial conditions or limits on the parameters need to be used. The output of the fit to a bin includes a text file printout of the result with name pattern, `fitresults*.txt`. When the code is executed to do fitting, it looks for the fit results printout first, and simply extracts the efficiency from the printout if the file exists, thereby not repeating any unnecessary fitting to save time. If a fit needs to be re-done, delete the corresponding printout file.


## Outputs
All outputs are placed in the specified output directory. The ROOT file of efficiencies is `eff.root`. There are HTML files to allow for quick browsing of plots. The image files and fit result printouts are placed under the `plots/` folder in the output directory. A summary text file listing all the effciencies is found at `plots/summary.txt`.


## Tutorial
*Tested on lxplus with CMSSW_5_3_18. The code does not involve any CMSSW, but this is the easiest way to pick up ROOT and RooFit. Set up a CMSSW release,*
```
cmsrel CMSSW_5_3_18
cd CMSSW_5_3_18/src
cmsenv
```
*or go to your favorite release and do* `cmsenv`.


1. Check out and compile the tag-and-probe code.
    ```
    git clone https://github.com/ksung25/TagAndProbe.git
    cd TagAndProbe
    source build.sh
    ```

2. Copy tutorial files and ntuples. The ntuples correspond to Run-1 data and MC.
  * The run script, binning files, and pileup reweighting file

        ```
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/run_example.sh ./
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/IsoMu24_eta2p1.bins ./
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/muon_selection.bins ./
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/PUWeights_2012.root ./
        ```
  * Probes ntuples for measuring efficiency of `HLT_IsoMu24_eta2p1` trigger

        ```
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/SingleMu_2012-22Jan2013_smubits.root ./
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/Summer12_DYJetsToLL_M-50_TuneZ2Star_smubits.root ./
        ```
  * Probes ntuples for measuring efficiency of tight muon selection (ID and isolation)

        ```
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/SingleMu_2012-22Jan2013_muselbits.root ./
        eos cp /eos/cms/store/user/ksung/TagAndProbeExample/Summer12_DYJetsToLL_M-50_TuneZ2Star_muselbits.root ./
        ```

3. Execute the run script.

    ```
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
    source run_example.sh
    ```

  The example performs counting for extracting efficiencies from MC, as well as for trigger efficiencies in data. Fitting is performed for measuring the muon selection efficiency in data. It takes about 20-30 mins to complete all the fits.

4. Look at the output. There should be 4 directories created with the outputs from the efficiency extraction. Start browsing though the summary plots from `plots.html`.
