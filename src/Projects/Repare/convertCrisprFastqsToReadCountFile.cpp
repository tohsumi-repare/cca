//////////////////////////////////////////////////////////////////////////////
// MolBioLib: A C++11 framework for rapidly developing bioinformatics tasks //
// Copyright (C)  2018  Repare Therapeutics                                 //
//////////////////////////////////////////////////////////////////////////////


/** \file convertCrisprFastqsToReadCountFile.cpp
 * Given fastq files, convert to a read count file.
 */


#include "src/include/MolBioLib.hpp"



int main(int argc, char* argv[])
{

  BeginCommandArguments
  CommandArgumentBoolDefaultDoc("AT_HOME", AT_HOME, true, "If true, then all files are assumed to be in /home, so /home/ is prepended to the input and output.  For Docker version only.");
  CommandArgumentStringDoc("IN", IN, "CSV of FASTQ files ending in .fastq");
  CommandArgumentFileInDoc("REF", REF, "FASTA reference file of sgRNAs.  Header must be of form gene_anything.");
  CommandArgumentSize_tDefaultDoc("GENE_POS_IN_REF", GENE_POS_IN_REF, 1, "0-based index of gene position in reference labels with underscore as delimiter.");
  CommandArgumentFileOut("OUT", OUT);
  CommandArgumentFileOut("OUT_GENES", OUT_GENES);
  // CommandArgumentFileOutDefaultDoc("MULTIPLE_ALIGN", MULTIPLE_ALIGN, "", "By default, bowtie2 doesn't output multiple alignments, only the best one.  Thus, use this only if -k N, where N > 1, is used as a parameter of bowtie2.");
  // CommandArgumentFileOut("NO_ALIGN", NO_ALIGN);
  CommandArgumentFileOut("SUMMARY", SUMMARY);
  EndCommandArguments

  if (AT_HOME) {
    OUT = "/home/" + OUT;
    OUT_GENES = "/home/" + OUT_GENES;
    SUMMARY = "/home/" + SUMMARY;
  }

  vector<string> fastqs, basenames;
  splitString(IN, ",", fastqs);
  size_t n = fastqs.size();
  basenames.resize(fastqs.size());
  for (size_t i = 0; i < n; ++i) {
    if (AT_HOME) fastqs[i] = "/home/" + fastqs[i];
    if (fileSize(fastqs[i]) == static_cast<ifstream::pos_type>(-1)) {
      cerr << "Error!  " << fastqs[i] << " does not exist.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    if (fastqs[i].size() < 7 || fastqs[i].substr(fastqs[i].size() - 6) != ".fastq") {
      cerr << "Error!  File " << fastqs[i] << " does not end with a .fastq extension.  Exiting." << endl;
    }
    basenames[i] = beforeString(fastqs[i], ".fastq");
  }


  set<string> sgrnas, genes;
  unordered_map<string, string> seq2gene, sgrna2seq, sgrna2gene;
  vector< unordered_map<string, size_t> > seq_hits(n), gene_hits(n);
  ReadOnlyStringFile ifpR(REF, true, true);
  for (size_t i = 0; !ifpR.fail(); i += 2) {
    ifpR.readRow(i);
    if (ifpR.fail()) break;
    string curr_sgrna = ifpR.line.substr(1);
    vector<string> tokens;
    splitString(curr_sgrna, "_", tokens);
    string curr_gene = tokens[GENE_POS_IN_REF];
    sgrnas.insert(curr_sgrna);
    genes.insert(curr_gene);
    sgrna2gene[curr_sgrna] = curr_gene;
    ifpR.readRow(i+1);
    for (size_t j = 0; j < n; ++j) {
      seq_hits[j][ifpR.line] = 0;
      gene_hits[j][curr_gene] = 0;
      sgrna2seq[curr_sgrna] = ifpR.line;
      seq2gene[ifpR.line] = curr_gene;
    }
  }
  ifpR.close();


  vector<size_t> num_reads(n), num_no_align(n);
  for (size_t i = 0; i < n; ++i) {
    ReadOnlyStringFile ifp(fastqs[i], true, true);
    cout << "Current working on " << fastqs[i] << endl << flush;
    for (size_t j = 0; !ifp.fail(); j += 4) {
      ifp.readRow(j);
      if (ifp.fail()) break;
      // Actually need sequence info...
      ifp.readRow(j+1);
      if (ifp.fail()) break;
      ++num_reads[i];
      if (seq_hits[i].find(ifp.line) == seq_hits[i].end()) {
        ++num_no_align[i];
      } else {
        ++seq_hits[i][ifp.line];
        ++gene_hits[i][seq2gene[ifp.line]];
      }
    }
    ifp.close();
  }


  Ofstream ofp(OUT);
  ofp << "sgRNA\tGene";
  for (size_t i = 0; i < n; ++i) ofp << "\t" << basenames[i];
  ofp << endl;
  for (auto j = sgrnas.begin(); j != sgrnas.end(); ++j) {
    ofp << (*j) << "\t" << sgrna2gene[(*j)];
    for (size_t i = 0; i < n; ++i) {
      ofp << "\t" << seq_hits[i][sgrna2seq[(*j)]];
    }
    ofp << endl;
  }
  ofp.close();


  // Only output if row is non-zero.   This is because
  // RNA-Seq packages do not handle zero rows well.
  Ofstream ofpG(OUT_GENES);
  ofpG << "Gene";
  for (size_t i = 0; i < n; ++i) ofpG << "\t" << basenames[i];
  ofpG << endl;
  for (auto j = genes.begin(); j != genes.end(); ++j) {
    ofpG << (*j);
    bool has_nonzero = false;
    for (size_t i = 0; i < n; ++i) {
      if (gene_hits[i][(*j)] != 0) {
        has_nonzero = true;
      }
    }
    if (has_nonzero) {
      for (size_t i = 0; i < n; ++i) {
        ofpG << "\t" << gene_hits[i][(*j)];
      }
      ofpG << endl;
    }
  }
  ofpG.close();


  Ofstream ofpS(SUMMARY);
  ofpS << "Base\tNum reads\tNum not aligned\tPercent aligned" << endl;
  for (size_t i = 0; i < n; ++i) {
    ofpS << basenames[i] << "\t" 
         << num_reads[i] << "\t" 
         << num_no_align[i] << "\t"
         << 100.0*(num_reads[i] - num_no_align[i])/num_reads[i]
         << endl;
  }
  ofpS.close();


  return 0;
}

