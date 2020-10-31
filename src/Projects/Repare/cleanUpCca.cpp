//////////////////////////////////////////////////////////////////////////////
// MolBioLib: A C++11 framework for rapidly developing bioinformatics tasks //
// Copyright (C)  2018  Repare Therapeutics                                 //
//////////////////////////////////////////////////////////////////////////////


/** \file cleanUpCca.cpp
 * Clean up CCA file by removing unused columns.
 */


#include "src/include/MolBioLib.hpp"



int main(int argc, char* argv[])
{

  BeginCommandArguments
  CommandArgumentStringDoc("IN", IN, "CCA");
  CommandArgumentFileOut("OUT", OUT);
  EndCommandArguments

  ReadOnlyTSVFile ifp(IN, true, true);
  ifp.readRow(0);
  size_t test_med = 0, ctrl_med = 0, score = 0, rank = 0,
         z_score = 0, z_rank = 0, jenks = 0, z_jenks = 0,
         uni_gene_start = 0, start = 0;
  if (ifp.tokens[2] == "Positive Control Gene") {
    start = 5;
    if (ifp.tokens[4] != "Number of sgRNAs") {
      cerr << "Error!  The fifth column should be Number of sgRNAs.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
  } else {
    start = 4;
    if (ifp.tokens[3] != "Number of sgRNAs") {
      cerr << "Error!  The fourth column should be Number of sgRNAs.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
  }
  for (size_t i = start; i < ifp.tokens.size(); ++i) {
    if (ifp.tokens[i] == "Test median depletion") test_med = i;
    else if (ifp.tokens[i] == "Ctrl median depletion") ctrl_med = i;
    else if (ifp.tokens[i] == "Score") score = i;
    else if (ifp.tokens[i] == "Rank") rank = i;
    else if (ifp.tokens[i] == "Negative non-parameteric Z-score") z_score = i;
    else if (ifp.tokens[i] == "Non-parametric Z-score rank") z_rank = i;
    else if (ifp.tokens[i] == "Jenks class") jenks = i;
    else if (ifp.tokens[i] == "Jenks class for non-parametric Z-scores") z_jenks = i;
    else if ((ifp.tokens[i].find("Uniprot") != string::npos ||
              ifp.tokens[i].find("Gene_ontology") != string::npos) &&
             uni_gene_start == 0) uni_gene_start = i;
  }

  Ofstream ofp(OUT);
  for (size_t j = 0; j < start; ++j) {
    ofp << ifp.tokens[j] << "\t";
  }
  ofp << ifp.tokens[test_med] << "\t"
      << ifp.tokens[ctrl_med] << "\t"
      << ifp.tokens[score] << "\t"
      << ifp.tokens[rank] << "\t"
      << ifp.tokens[z_score] << "\t"
      << ifp.tokens[z_rank] << "\t"
      << ifp.tokens[jenks] << "\t"
      << ifp.tokens[z_jenks];
  for (size_t j = uni_gene_start; j < ifp.tokens.size(); ++j) {
    ofp << "\t" << ifp.tokens[j];
  }
  ofp << endl;
  for (size_t i = 1; !ifp.fail(); ++i) {
    ifp.readRow(i);
    for (size_t j = 0; j < start; ++j) {
      ofp << ifp.tokens[j] << "\t";
    }
    ofp << ifp.tokens[test_med] << "\t"
        << ifp.tokens[ctrl_med] << "\t"
        << ifp.tokens[score] << "\t"
        << ifp.tokens[rank] << "\t"
        << ifp.tokens[z_score] << "\t"
        << ifp.tokens[z_rank] << "\t"
        << ifp.tokens[jenks] << "\t"
        << ifp.tokens[z_jenks];
    for (size_t j = uni_gene_start; j < ifp.tokens.size(); ++j) {
      ofp << "\t" << ifp.tokens[j];
    }
    ofp << endl;
  }
  ofp.close();
  ifp.close();


  return 0;
}

