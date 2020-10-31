//////////////////////////////////////////////////////////////////////////////
// MolBioLib: A C++11 framework for rapidly developing bioinformatics tasks //
// Copyright (C)  2019  Repare Therapeutics                                 //
//////////////////////////////////////////////////////////////////////////////


/** \file ccaGatherQC.cpp
 * Gather CCA QC files and plots and output a HTML file.
 * Assumes a standard CCA output directory with OUT_MAX used.
 */

// May want to either pass the compiler flag -DDEBUG and/or -DPROG_DEBUG or
// uncomment either or both of the below lines
// #define DEBUG
// #define PROG_DEBUG


#include "src/include/MolBioLib.hpp"



int main(int argc, char* argv[])
{
    BeginCommandArguments
    CommandArgumentStringDoc("COUNTS", COUNTS, "Pathway and filename to the input COUNTS file.");
    CommandArgumentStringDoc("CORREL_NORMED_COUNTS", CORREL_NORMED_COUNTS, "Path and filename to the correl.normed_counts.tsv file.  This is used since every run, regardless of resistors, etc. has this as one of the files.");
    CommandArgumentSize_tDefaultDoc("BARPLOT_TOP_N", BARPLOT_TOP_N, 25, "Top number of genes to bar plot.");
    CommandArgumentSize_tDefaultDoc("RESISTOR_TOP_N", RESISTOR_TOP_N, 100, "Top number of genes to output in resistor screens.");
    CommandArgumentSize_tDefaultDoc("GSEA_TOP_N", GSEA_TOP_N, 100, "Top number of pathways to output.");
    CommandArgumentDoubleDefaultDoc("MIN_T0_CORREL", MIN_T0_CORREL, 0.9, "Minimum T0 correlation not to give a warning.");
    CommandArgumentStringDefaultDoc("PREFIX", PREFIX, "", "Prepend prefix (e.g.  if moving files to a new subdirectory) to all output file references.");
    CommandArgumentBoolDefaultDoc("AT_HOME", AT_HOME, false, "If true, add /home to COUNTS to find it.");
    CommandArgumentFileOut("HTML", HTML);
    CommandArgumentFileOut("TSV", TSV);
    EndCommandArguments

    ifstream::pos_type no_file = static_cast<ifstream::pos_type>(-1);

    bool doCounts = true;
    if (COUNTS.find(",") != string::npos) {
      cerr << "Warning.  Skipping plots of counts because there are multiple COUNT files." << endl;
      doCounts = false;
    } else if (fileSize(CORREL_NORMED_COUNTS) == no_file) {
      cerr << "Error!  " << CORREL_NORMED_COUNTS << " does not exist.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    string head = beforeString(CORREL_NORMED_COUNTS, ".correl.normed_counts.tsv");
    string counts_home = AT_HOME ? "/home/" + COUNTS : COUNTS;
    string counts_prefix = "";
    if (COUNTS.find("/") != string::npos) {
      counts_prefix = PREFIX + afterString(COUNTS, "/");
    } else {
      counts_prefix = PREFIX + COUNTS;
    }
    string ofp_head = PREFIX + afterString(head, "/");

    Ofstream ofp(HTML);
    Ofstream tsv(TSV);
    ofp << "<html>" << endl;
    ofp << "<head>" << endl 
        << "<title>" << endl;
    ofp << ofp_head << endl;
    ofp << "</title>" << endl;
    ofp << "<style type=\"text/css\">" << endl;
    ofp << "table.hovertable {" << endl;
    ofp << "    font-family: verdana,arial,sans-serif;" << endl;
    ofp << "    font-size:11px;" << endl;
    ofp << "    color:#333333;" << endl;
    ofp << "    border-width: 1px;" << endl;
    ofp << "    border-color: #999999;" << endl;
    ofp << "    border-collapse: collapse;" << endl;
    ofp << "}" << endl;
    ofp << "table.hovertable th {" << endl;
    ofp << "    background-color:#c3dde0;" << endl;
    ofp << "    border-width: 1px;" << endl;
    ofp << "    padding: 8px;" << endl;
    ofp << "    border-style: solid;" << endl;
    ofp << "    border-color: #a9c6c9;" << endl;
    ofp << "}" << endl;
    ofp << "table.hovertable tr {" << endl;
    ofp << "    background-color:#d4e3e5;" << endl;
    ofp << "}" << endl;
    ofp << "table.hovertable td {" << endl;
    ofp << "    border-width: 1px;" << endl;
    ofp << "    padding: 8px;" << endl;
    ofp << "    border-style: solid;" << endl;
    ofp << "    border-color: #a9c6c9;" << endl;
    ofp << "}" << endl;
    ofp << "</style>" << endl;
    ofp << "</head>" << endl;
    ofp << "<body>" << endl;
    if (doCounts && fileSize(counts_home + ".init_input_qc.Total_reads_per_sample.png") != no_file && fileSize(counts_home + ".init_input_qc.Distribution_of_total_reads_per_sample.png") != no_file) {
      ofp << "<h1>Inital QC of raw read counts.</h1><br>" << endl;
      ofp << "<h2>Total reads per sample on raw data and distribution:</h2>" << endl
          << "<img src=\"" << counts_prefix << ".init_input_qc.Total_reads_per_sample.png\", height = 500> &nbsp;&nbsp;&nbsp; "
          << "<img src=\"" << counts_prefix << ".init_input_qc.Distribution_of_total_reads_per_sample.png\", height = 500><br>" << endl;
      ofp << "<h3>All bars above should be above the red line.</h3><br>" << endl;
      if (doCounts && fileSize(counts_home + ".init_input_qc.Dist_noness_reads_per_sample.png") != no_file && fileSize(counts_home + ".init_input_qc.Dist_ess_reads_per_sample.png") != no_file) {
        ofp << "<h2>Distribution of essential and nonessential sgRNAs:</h2>" << endl
            << "<img src=\"" << counts_prefix << ".init_input_qc.Dist_ess_reads_per_sample.png\", height = 500> &nbsp;&nbsp;&nbsp " 
            << "<img src=\"" << counts_prefix << ".init_input_qc.Dist_noness_reads_per_sample.png\", height = 500><br>" << endl;
        ofp << "<h3>It is all right if the essential genes' median is below the line.  This happens when a lot of the essential genes drop out.</h3><br>" << endl;
      }
      if (doCounts && fileSize(counts_home + ".init_input_qc.Clustering_of_samples.png") != no_file) {
        ofp << "<h2>Clustering of raw read counts in linear space:</h2>" << endl
            << "<img src=\"" << counts_prefix << ".init_input_qc.Clustering_of_samples.png\", height = 1000><br>" << endl;
        ofp << "<h3>All samples should cluster logically - T0s should cluster, MUT and WTs should cluster separately, time points should cluster, and finally replicates should cluster.</h3><br>" << endl;
      }
      ofp << "<br><hr><br>" << endl;
    }


    if (fileSize(head + ".correl.normed_counts.cor.png") != no_file) {
      ofp << "<h1>Analysis of normed read counts and foldchange:</h1><br>" << endl;
      ofp << "<h2>Correlation of normed read counts:</h2><br>" << endl;
#if 0
      ofp << "<img src=\"" << ofp_head << ".correl.normed_counts.cor.png\", height = 1000><br>" << endl;
      ofp << "<h3>The sample names are on the diagonal.  Axes values alternate between left and right and top and bottom.</h3><br>" << endl;
#endif
      ReadOnlyTSVFile ifp(head + ".correl.normed_counts.cor.tsv", true, true);
      ifp.readRow(0);
      vector<string> header = ifp.tokens;
      for (size_t i = 1; !ifp.fail(); ++i) {
        ifp.readRow(i);
        for (size_t j = 1; j < i; ++j) {
          if (header[j].find("T0.of.") != string::npos &&
              ifp.tokens[0].find("T0.of.") != string::npos &&
              convertFromString<double>(ifp.tokens[j]) < MIN_T0_CORREL) {
            ofp << "<h4><font color=\"red\">Warning!  " << header[j] << " vs " << ifp.tokens[0] << " has a correlation of " << ifp.tokens[j] << ".</font><br>" << endl;
            tsv << "T0 correlation greater than " << convertToString<double>(MIN_T0_CORREL) << "\t" << header[j] << "\t" << ifp.tokens[0] << "\t" << ifp.tokens[j] << endl;
          }
        }
      }
      ifp.close();
      ofp << "<br><hr><br>" << endl;
    }

    if (fileSize(head + ".killing.essential.killing_statistics.tsv") != no_file || fileSize(head + ".killing.test.depletion_all_top_hits_as_group_med.png") != no_file) {
      // ofp << "<h1>Depletion and killing distributions:</h1><br>" << endl;
      ofp << "<h1>Depletion distributions:</h1><br>" << endl;
      if ( fileSize(head + ".killing.depletion_noness_med.png") != no_file &&
          fileSize(head + ".killing.depletion_ess_med.png") != no_file ) {
        ofp << "<h2>Comparative plots of control and test samples.  Ordering is: 1. essential genes distribution, 2. nonessential genes distribution, 3. median depletion of control and test samples for the top 300 hits:</h2>" << endl;
        ofp << "<img src=\"" << ofp_head << ".killing.depletion_ess_med.png\", height = 300> &nbsp;&nbsp;&nbsp;&nbsp; <img src=\"" << ofp_head <<".killing.depletion_noness_med.png\", height = 300> &nbsp;&nbsp;&nbsp; <img src=\"" << ofp_head << ".killing.depletion_top_300_hits_med.png\", height = 300> <br>" << endl;
        ofp << "<h3>In the above plots, we expect the test to have greater median depletion than control for isogenic and chemogenomic screens.</h3><br>" << endl;
      }

      if (fileSize(head + ".killing.control.depletion_all_top_hits_as_group_med.png") != no_file) {
        ofp << "<h2>Comparative plots of control median depletions among various groups:</h2>" << endl;
        ofp << "<img src=\"" << ofp_head << ".killing.control.depletion_all_top_hits_as_group_med.png\", height = 750> <br>" << endl;
        ofp << "<h3>For isogenic and chemogenomic screens, the essential genes should be the right-most, followed by other genes and the nonessential genes.  Ideally, the top 300 hits should overlap with the distribution of the other genes and nonessential genes.<h3><br>" << endl;
      }
      if (fileSize(head + ".killing.test.depletion_all_top_hits_as_group_med.png") != no_file) {
        ofp << "<h2>Comparative plots of test median depletions among various groups:</h2>" << endl;
        ofp << "<img src=\"" << ofp_head << ".killing.test.depletion_all_top_hits_as_group_med.png\", height = 750> <br>" << endl;
        ofp << "<h3>For isogenic and chemogenomic screens, the essential genes should be the right-most followed by the top 300 hits, followed by other genes, and with the nonessential genes being the left-most.<h3><br>" << endl;
      }
      ofp << "<br><hr><br>" << endl;
   }


   if (fileSize(head + ".stratas.png") != no_file) {
     ofp << "<h1>Cutoff stratas.</h1>" << endl;
     if (fileSize(head + ".stratas.png") != no_file) {
       if (fileSize(head + ".stratas_z.png") == no_file) {
         ofp << "<h2>Stratification of top CCA score hits plots with the top-most hit on the right:</h2>" << endl;
         ofp << "<img src=\"" << ofp_head << ".stratas.png\", height = 500><br>" << endl;
       } else {
         ofp << "<h2>Stratification of top CCA score and Z-score hits plots with the top-most hit on the right:</h2>" << endl;
         ofp << "<img src=\"" << ofp_head << ".stratas.png\", height = 400> &nbsp;&nbsp&nbsp; <img src =\"" << ofp_head << ".stratas_z.png\", height = 400><br>" << endl;
       }
       ofp << "<h3>Points on a vertical line mean that the point belongs to the set to the right of that line.</h3><br>" << endl;
     }
     ofp << "<br><hr><br>" << endl;
   }

   string cerr_file = "";
   if (fileSize(head + ".cerr") != no_file) cerr_file = head + ".cerr";
   else if (fileSize(head + ".cca.cerr") != no_file) cerr_file = head + ".cca.cerr";
   if (cerr_file != "") {
     vector<string> warning_lines, error_lines;
     bool has_count_median_diff = false;
     ReadOnlyStringFile ifp(cerr_file, true, true);
     for (size_t i = 0; !ifp.fail(); ++i) {
       ifp.readRow(i);
       string upper_line = convertStringToUpperCase(ifp.line);
       if (upper_line.find("ERROR") != string::npos) {
         error_lines.push_back(ifp.line);
       } else if (upper_line.find("WARNING") != string::npos) {
         if (upper_line.find("QUANTILE SIZE IS ZERO") == string::npos &&
             upper_line.find("WARNING MESSAGE:") == string::npos) {
           warning_lines.push_back(ifp.line);
           if (upper_line.find("WARNING: COUNT MEDIAN DIFFERENCE IN") != string::npos) {
             has_count_median_diff = true;
           }
         }
       }
     }
     ifp.close(); 

     if (error_lines.size() > 0 || warning_lines.size() > 0) {
       ofp << "<h1>CCA warnings and errors.</h1>" << endl;
       if (error_lines.size() > 0) {
         ofp << "<h2>Errors:</h2>" << endl;
         for (size_t j = 0; j < error_lines.size(); ++j) {
           ofp << "<h3><font color=\"red\">" << error_lines[j] << "</font></h3>" << endl;
         }
       }
       if (warning_lines.size() > 0) {
         ofp << "<h2>Warnings:</h2>" << endl;
         for (size_t j = 0; j < warning_lines.size(); ++j) {
           ofp << "<h3><font color=\"red\">" << warning_lines[j] << "</font></h3>" << endl;
         }
         if (has_count_median_diff) {
           ofp << "<h3><p>A count median difference warning means that the difference in median counts between samples is too large.  Check the file " << cerr_file << " to find the samples and medians and see in which samples the biggest differences exists.  This problem could be addressed by rerunning CCA with the option MEDIAN_NORMALIZATION=true, but experience shows this could actually make things worse.  This problem could be due to something suboptimal in the biological experiment.  The warning also often occurs when mixing multiple cell lines in a CCA run or when using targeted panels.</p></h3>" << endl;
         }
       }
       ofp << "<br><hr><br>" << endl;
     }
   }

   if (fileSize(head + ".cca.tsv") != no_file) {
     vector<string> top_hits;
     ReadOnlyTSVFile ifp(head + ".cca.tsv", true, true);
     ifp.readRow(0);
     vector<string> header = ifp.tokens;
     bool has_pc = false, has_jenks = false, has_detrend = false,
          has_z = false, has_zjenks = false, has_kill = false,
          has_uniprot = false, has_cortellis = false;
     size_t gene = 0, ess = 0, pc = 0, tsg = 0, num = 0, rank = 0, zrank = 0,
            jenks = 0, zjenks = 0, kill = 0, test_var = 0, ctrl_var = 0,
            test_t0_var = 0, ctrl_t0_var = 0, detrend = 0, 
            uniprot = 0, cortellis = 0,
            paralogs = 0;
     for (size_t i = 0; i < header.size(); ++i) {
       if (header[i] == "Gene") gene = i;
       else if (header[i] == "Core Ess Or Noness") ess = i;
       else if (header[i] == "Positive Control Gene") {
         has_pc = true;
         pc = i;
       } else if (header[i] == "Tumor Suppressor Gene") tsg = i;
       else if (header[i] == "Number of sgRNAs") num = i;
       else if (header[i] == "Rank") rank = i;
       else if (header[i] == "Non-parametric Z-score rank") {
         has_z = true;
         zrank = i;
       } else if (header[i] == "Jenks class") {
         has_jenks = true;
         jenks = i;
       } else if (header[i] == "Jenks class for non-parametric Z-scores") {
         has_zjenks = true;
         zjenks = i;
       } else if (header[i] == "Percent of essential genes with a lower median test depletion") {
	 has_kill = true;
	 kill = i;
       }
       else if (header[i] == "IQR/(Q3 + Q1) of test depletion") test_var = i;
       else if (header[i] == "IQR/(Q3 + Q1) of ctrl depletion") ctrl_var = i;
       else if (header[i] == "IQR/(Q3 + Q1) test T0 count") test_t0_var = i;
       else if (header[i] == "IQR/(Q3 + Q1) ctrl T0 count") ctrl_t0_var = i;
       else if (header[i] == "Detrended Rank") {
         has_detrend = true;
         detrend = i;
       } else if (header[i] == "Ensembl paralogs") paralogs = i;
       else if (header[i] == "Uniprot_function") {
         has_uniprot = true;
         if (header[i+1] != "Gene_ontology_biological_process") {
           cerr << "Warning!  Non-standard uniprot annotation.  Not used." << endl;
           has_uniprot = false;
         }
         uniprot = i;
       } else if (header[i] == "Cortellis drug summary") {
         has_cortellis = true;
         cortellis = i;
       }
     }
     size_t max_rank = 300, max_zrank = 300;
     if (has_jenks || has_zjenks) {
       if (has_jenks) max_rank = 0;
       if (has_zjenks) max_zrank = 0;
       for (size_t i = 1; !ifp.fail(); ++i) {
         ifp.readRow(i);
         if (has_jenks) {
           size_t curr_rank = convertFromString<size_t>(ifp.tokens[rank]);
           // if (has_jenks && ifp.tokens[jenks] == "Beyond p = 0.05 cutoff - 1" &&
           if (has_jenks && ifp.tokens[jenks] == "3" &&
               max_rank < curr_rank) max_rank = curr_rank;
         }
         if (has_zjenks) {
           size_t curr_zrank = convertFromString<size_t>(ifp.tokens[zrank]);
           // if (has_zjenks && ifp.tokens[zjenks] == "Beyond p = 0.05 cutoff - 1" &&
           if (has_zjenks && ifp.tokens[zjenks] == "3" &&
               max_zrank < curr_zrank) max_zrank = curr_zrank;
           }
       }
       ifp.resetToBegin();
     }

     for (size_t i = 1; !ifp.fail(); ++i) {
       ifp.readRow(i);
       if (i <= BARPLOT_TOP_N) {
         top_hits.push_back(ifp.tokens[gene]);
       }
     }
     ifp.resetToBegin();

     ofp << "<h1>Bar plots of top " << top_hits.size() << " hits:</h1><br>" << endl;
     for (size_t i = 0; i < top_hits.size(); ++i) {
       if (fileSize(head + ".foldchange_barplot." + convertToString<size_t>(i+1) + "." + top_hits[i] + ".png") != no_file) {
         ofp << "<h2>Gene " << top_hits[i] << ":</h2>" << endl;
         ofp << "<img src=\"" << ofp_head << ".foldchange_barplot." << convertToString<size_t>(i+1) + "." << top_hits[i] << ".png\", width = 1750><br><br><br><br>" << endl;
       }
     }

     //    rank     line
     Table<size_t, string> good, goodz;
     for (size_t i = 1; !ifp.fail(); ++i) {
       ifp.readRow(i);
       size_t curr_rank = convertFromString<size_t>(ifp.tokens[rank]);
       size_t curr_zrank = has_z ? convertFromString<size_t>(ifp.tokens[zrank]) : 0;
       string line = ifp.tokens[gene] + "\t" + ifp.tokens[ess];
       if (has_pc) line += "\t" + ifp.tokens[pc];
       line += "\t" + ifp.tokens[tsg] + "\t" + ifp.tokens[num] + 
               "\t" + ifp.tokens[rank];
       if (has_jenks) line += "\t" + ifp.tokens[jenks];
       if (has_z) {
         line += "\t" + ifp.tokens[zrank];
         if (has_zjenks) line += "\t" + ifp.tokens[zjenks];
       }
       line += "\t";
       if (has_kill) {
         if (convertFromString<double>(ifp.tokens[kill]) < 20.0) line += "Y";
	 line += "\t";
       } 
       double var = 0.0;
       if (trimSpacesString(ifp.tokens[test_var]) != "") {
         var = convertFromString<double>(ifp.tokens[test_var]);
         if (var > 2.0) line += "Y";
         line += "\t";
         if (var > 10.0) line += "Y";
         line += "\t";
       } else {
         line += "\t\t";
       }
       if (trimSpacesString(ifp.tokens[ctrl_var]) != "") {
         var = convertFromString<double>(ifp.tokens[ctrl_var]);
         if (var > 2.0) line += "Y";
         line += "\t";
         if (var > 10.0) line += "Y";
         line += "\t";
       } else {
         line += "\t\t";
       }
       if (trimSpacesString(ifp.tokens[test_t0_var]) != "") {
         var = convertFromString<double>(ifp.tokens[test_t0_var]);
         if (var > 2.0) line += "Y";
       }
       line += "\t";
       if (trimSpacesString(ifp.tokens[ctrl_t0_var]) != "") {
         var = convertFromString<double>(ifp.tokens[ctrl_t0_var]);
         if (var > 2.0) line += "Y";
       }
       if (has_detrend) {
         line += "\t";
         if (trimSpacesString(ifp.tokens[detrend]) != "") {
           size_t detrend_rank = convertFromString<size_t>(ifp.tokens[detrend]);
           if (has_jenks && detrend_rank > max_rank) {
             line += ifp.tokens[detrend];
           }
         }
       }
       line += "\t" + ifp.tokens[paralogs];
       if (has_uniprot) {
         for (size_t j = uniprot; j < (uniprot+2); ++j) {
           line += "\t" + ifp.tokens[j];
         }
       }
       if (has_cortellis) {
         line += "\t" + ifp.tokens[cortellis];
       }

       if (curr_rank <= max_rank) good.push_back(curr_rank, line);
       else if (has_z && curr_zrank <= max_zrank) goodz.push_back(curr_zrank, line);
     }
     ifp.close();
     if (has_z) sortTable<0>(goodz);   // Only need to sort Z - rest in order.

     ofp << "<h1>CCA top hits QC:</h1>" << endl;
     tsv << "CCA top hits QC:" << endl;
     ofp << "<table class = \"hovertable\">" << endl;
     ofp << "<tr><th>Gene</th><th>Core Ess or Noness</th>";
     tsv << "Gene\tCore Ess or Noness";
     if (has_pc) {
       ofp << "<th>Positive Control</th>";
       tsv << "\tPositive Control";
     }
     ofp << "<th>Tumor Suppressor Gene</th><th>Number of sgRNAs</th>"
         << "<th>Rank</th>";
     tsv << "\tTumor Suppressor Gene\tNumber of sgRNAs\tRank";
     if (has_jenks) {
       ofp << "<th>Jenks class</th>";
       tsv << "\tJenks class";
     }
     if (has_z) {
       ofp << "<th>Non-parametric Z-score rank</th>";
       tsv << "\tNon-parametric Z-score rank";
       if (has_zjenks) {
         ofp << "<th>Jenks class for non-parametric Z-scores</th>";
         tsv << "\tJenks class for non-parametric Z-scores";
       }
     }
     if (has_kill) {
       ofp << "<th>Less than 20 percent of essential genes with a lower killing measure</th>";
     }
     ofp << "<th>Has &gt; 2 variance in TEST</th><th>Has &gt; 10 variance in TEST</th><th>Has &gt; 2 variance in CTRL</th><th>Has &gt; 10 variance in CTRL</th><th>Has &gt; 2 variance in TEST T0</th><th>Has &gt; 2 variance in CTRL T0</th>";
     if (has_kill) {
       tsv << "\tLess than 20 percent of essential genes with a lower killing measure";
     }
     tsv << "\tHas > 2 variance in TEST\tHas > 10 variance in TEST\tHas > 2 variance in CTRL\tHas > 10 variance in CTRL\tHas > 2 variance in TEST T0\tHas > 2 variance in CTRL T0";
     if (has_detrend) {
       ofp << "<th>Detrended rank if detrending pushes genes out of top (p &lt; 0.05) genes</th>";
       tsv << "\tDetrended rank if detrending pushes genes out of top (p < 0.05) genes";
     }
     ofp << "<th>Ensembl paralogs</th>";
     tsv << "\tEnsembl paralogs";
     if (has_uniprot) {
       ofp << "<th>Uniprot function</th><th>Gene ontology biological process</th>";
       tsv << "\tUniprot function<\tGene ontology biological process";
     }
     if (has_cortellis) {
       ofp << "<th>Cortellis drug summary</th>";
       tsv << "\tCortellis drug summary";
     }
     ofp << "</tr>" << endl;
     tsv << endl;
  

     size_t wrap_cols = has_uniprot ? 3 : 1;
     if (has_cortellis) ++wrap_cols;
     for (size_t i = 0; i < good.size(); ++i) {
       size_t rank;
       string line = "";
       good.getRow(i, rank, line);
       ofp << "<tr onmouseover=\"this.style.backgroundColor='#ffff66';\" onmouseout=\"this.style.backgroundColor='#d4e3e5';\">";
       vector<string> tokens;
       splitString(line, "\t", tokens);
       for (size_t j = 0; j < tokens.size() - wrap_cols; ++j) {
         ofp << "<td>" << tokens[j] << "</td>";
       }
       // Force word wrap for paralogs and possibly uniprot
       for (size_t j = tokens.size() - wrap_cols; j < tokens.size(); ++j) {
         ofp << "<td><div style = \"width:200px; word-wrap: break-word\">" << tokens[j] << "</div></td>";
       }
       ofp << "</tr>" << endl;
       tsv << line << endl;
     }
     if (has_z) {
       for (size_t i = 0; i < goodz.size(); ++i) {
         size_t rank;
         string line = "";
         goodz.getRow(i, rank, line);
         ofp << "<tr onmouseover=\"this.style.backgroundColor='#ffff66';\" onmouseout=\"this.style.backgroundColor='#d4e3e5';\">";
         vector<string> tokens;
         splitString(line, "\t", tokens);
         for (size_t j = 0; j < tokens.size() - wrap_cols; ++j) {
           ofp << "<td>" << tokens[j] << "</td>";
         }
         // Force word wrap for paralogs and possibly uniprot
         for (size_t j = tokens.size() - wrap_cols; j < tokens.size(); ++j) {
           ofp << "<td><div style = \"width:200px; word-wrap: break-word\">" << tokens[j] << "</div></td>";
         }
         ofp << "</tr>" << endl;
         tsv << line << endl;
       }
     }
     ofp << "</table>" << endl;
     ofp << "<br><hr><br>" << endl;

   } else if (fileSize(head + ".resistor_cca.tsv") != no_file) {
     bool has_pc = false, has_uniprot = false;
     size_t gene = 0, num = 0, ess = 0, pc = 0, tsg = 0, var = 0, var_t0 = 0,
            uniprot = 0;
     ReadOnlyTSVFile ifp(head + ".resistor_cca.tsv", true, true);
     ifp.readRow(0);
     for (size_t i = 0; i < ifp.tokens.size(); ++i) {
       if (ifp.tokens[i] == "Gene") gene = i;
       else if (ifp.tokens[i] == "Core Ess Or Noness") ess = i;
       else if (ifp.tokens[i] == "Number of sgRNAs per gene") num = i;
       else if (ifp.tokens[i] == "Positive Control Gene") {
         has_pc = true;
         pc = i;
       } else if (ifp.tokens[i] == "Tumor Suppressor Gene") tsg = i;
       else if (ifp.tokens[i] == "IQR/(Q3 + Q1)") var = i;
       else if (ifp.tokens[i] == "IQR/(Q3 + Q1) of T0") var_t0 = i;
       else if (ifp.tokens[i] == "Uniprot_function") {
         has_uniprot = true;
         if (ifp.tokens[i+1] != "Gene_ontology_biological_process") {
           cerr << "Warning!  Non-standard uniprot annotation.  Not used." << endl;
           has_uniprot = false;
         }
         uniprot = i;
       }
     }
     
     ofp << "<h1>CCA top resistors QC:</h1>" << endl;
     ofp << "<table class = \"hovertable\">" << endl;
     ofp << "<tr><th>Gene</th><th>Core Ess or Noness</th>";
     tsv << "Gene\tCore Ess or Noness";
     if (has_pc) {
       ofp << "<th>Positive Control</th>";
       tsv << "\tPositive Control";
     }
     ofp << "<th>Tumor Suppressor Gene</th><th>Number of sgRNAs</th>"
         << "<th>Rank</th>";
     tsv << "\tTumor Suppressor Gene\tNumber of sgRNAs\tRank";
     ofp << "<th>Has &gt; 2 variance in TEST</th><th>Has &gt; 10 variance in TEST</th><th>Has &gt; 2 variance in T0</th>";
     tsv << "\tHas > 2 variance in TEST\tHas > 10 variance in TEST\tHas > 2 variance in T0";
     if (has_uniprot) {
       ofp << "<th>Uniprot function</th><th>Gene ontology biological process</th>";
       tsv << "Uniprot function<\tGene ontology biological process";
     }
     tsv << endl;
     ofp << "</tr>" << endl;

     for (size_t i = 1; i <= RESISTOR_TOP_N && !ifp.fail(); ++i) {
       ifp.readRow(i);
       ofp << "<tr onmouseover=\"this.style.backgroundColor='#ffff66';\" onmouseout=\"this.style.backgroundColor='#d4e3e5';\">";
       ofp << "<td>" << ifp.tokens[gene] << "</td>";
       ofp << "<td>" << ifp.tokens[ess] << "</td>";
       string line = ifp.tokens[gene] + "\t" + ifp.tokens[ess];
       if (has_pc) {
         ofp << "<td>" << ifp.tokens[pc] << "</td>";
         line += "\t" + ifp.tokens[pc];
       }
       ofp << "<td>" << ifp.tokens[tsg] << "</td>";
       ofp << "<td>" << ifp.tokens[num] << "</td>";
       ofp << "<td>" << i << "</td>";
       line += "\t" + ifp.tokens[tsg] + "\t" + ifp.tokens[num] +
               "\t" + convertToString<size_t>(i);
       double val = 0.0;
       if (trimSpacesString(ifp.tokens[var]) != "") {
         val = convertFromString<double>(ifp.tokens[var]);
       }
       line += "\t";
       if (val > 2.0) {
         ofp << "<td>Y</td>";
         line += "Y";
       } else {
         ofp << "<td></td>";
       }
       line += "\t";
       if (val > 10.0) {
         ofp << "<td>Y</td>";
         line += "Y";
       } else {
         ofp << "<td></td>";
       }
       val = 0.0;
       if (trimSpacesString(ifp.tokens[var_t0]) != "") {
         val = convertFromString<double>(ifp.tokens[var_t0]);
       }
       line += "\t";
       if (val > 2.0) {
         ofp << "<td>Y</td>";
         line += "Y";
       } else {
         ofp << "<td></td>";
       }
       if (has_uniprot) {
         ofp << "<td>" << ifp.tokens[uniprot] << "</td><td>" << ifp.tokens[uniprot+1] << "</td>";
         line += ifp.tokens[uniprot] + "\t" + ifp.tokens[uniprot+1];
       }
       ofp << "</tr>" << endl;
       tsv << line << endl;
     }
     ifp.close();
     ofp << "</table>" << endl;
   }

   ofp << "</body>" << endl;
   ofp << "</html>" << endl;
   tsv.close();
   ofp.close();



   return 0;
}

