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

  string dataDir = "/data/fmazzasc/PbPb_2body/trees";
  string tableDir = "/data/fmazzasc/PbPb_2body/old_stuff";

  string lsString = likeSign ? "_pass3_LS.root" : "_pass3.root";

  string inFileNameQ = "HyperTritonTree_18q";
  string inFileArgQ = dataDir + "/" + inFileNameQ + lsString;

  string inFileNameR = "HyperTritonTree_18r";
  string inFileArgR = dataDir + "/" + inFileNameR + lsString;

  TChain inputChain("_custom/fTreeV0");
  inputChain.AddFile(inFileArgQ.data());
  inputChain.AddFile(inFileArgR.data());

  string outFileName = "DataTable_18qr";
  string outFileArg = tableDir + "/" + outFileName + lsString;


  TTreeReader fReader(&inputChain);
  TFile outFile(outFileArg.data(), "RECREATE");

  TTreeReaderArray<RHyperTritonHe3pi> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"}; 
  Table2 tree("DataTable", "Data Table");
  while (fReader.Next())
  {
    for (auto &RHyper : RHyperVec)
      tree.Fill(RHyper, *RColl);
  }

  outFile.cd();
  tree.Write();
  outFile.Close();
  std::cout << "\nDerived tree from Data generated!\n" << std::endl;

}
