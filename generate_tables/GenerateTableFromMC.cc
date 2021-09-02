#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "AliAnalysisTaskHyperTriton2He3piML.h"
#include "AliPID.h"

#include "helpers/Common.h"
#include "helpers//GenTable2.h"
#include "helpers/Table2.h"

void GenerateTableFromMC()
{
  gRandom->SetSeed(1995);

  string inFileName = "HyperTritonTree_20g7.root";
  string inFileArg = "/data/fmazzasc/PbPb_2body/trees/new_trees/" + inFileName;

  string outFileName = "SignalTable_20g7_flat_pt.root";

  string bwFileName;
  int cent_num;
  string outFileArg = "/data/fmazzasc/PbPb_2body/" + outFileName;

  TFile *inFile = new TFile(inFileArg.data(), "READ");
  TTreeReader fReader("_default/fTreeV0", inFile);
  TTreeReaderArray<RHyperTritonHe3piFull> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderArray<SHyperTritonHe3pi> SHyperVec = {fReader, "SHyperTriton"};
  TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"};

  // new flat tree with the features
  TFile outFile(outFileArg.data(), "RECREATE");
  Table2 table("SignalTable", "Signal Table");
  GenTable2 genTable("GenTable", "Generated particle table");

  while (fReader.Next())
  {
    for (auto &SHyper : SHyperVec)
    {
      int ind = SHyper.fRecoIndex;
      if (ind >= 0)
      {
        auto &RHyper = RHyperVec[ind];
        table.Fill(SHyper, RHyper, *RColl);
      }
      else
      {
        table.Fill(SHyper, *RColl);
      }
    }
  }

  outFile.cd();
  table.Write();
  outFile.Close();
}
