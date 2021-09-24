#include <iostream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

using namespace std;

#include "helpers/Common.h"
#include "helpers/Table2.h"


void GenerateTableFromData(bool likeSign = false)
{

  string dataDir = "/data/fmazzasc/PbPb_2body/no_pt_cut";
  string tableDir = "/data/fmazzasc/PbPb_2body/no_pt_cut";

  string lsString = likeSign ? "LS.root" : ".root";

  string inFileNameQ = "HyperTritonTree_18q";
  string inFileArgQ = dataDir + "/" + inFileNameQ + lsString;

  string inFileNameR = "HyperTritonTree_18r";
  string inFileArgR = dataDir + "/" + inFileNameR + lsString;

  TChain inputChain("_custom/fTreeV0");
  // inputChain.AddFile(inFileArgQ.data());
  inputChain.AddFile(inFileArgR.data());

  string outFileName = "DataTable_18r";
  string outFileArg = tableDir + "/" + outFileName + lsString;


  TTreeReader fReader(&inputChain);
  TFile outFile(outFileArg.data(), "RECREATE");

  TTreeReaderArray<RHyperTritonHe3pi> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"}; 
  Table2 tree("DataTable", "Data Table");
  int counter = 0;
  while (fReader.Next() and counter < 10)
  {
    for (auto &RHyper : RHyperVec)
      tree.Fill(RHyper, *RColl);
    counter ++;
  }

  outFile.cd();
  tree.Write();
  outFile.Close();
  std::cout << "\nDerived tree from Data generated!\n" << std::endl;

}
