//////////////////////////////////////////////////////////////////////////////
// MolBioLib: A C++11 framework for rapidly developing bioinformatics tasks //
// Copyright (C)  2018  Repare Therapeutics                                 //
//////////////////////////////////////////////////////////////////////////////


/** \file crisprCountsAnalysis.cpp
 * Given table of counts + which are controls and WT and T0,
 * get gene top hits.
 */


#include "src/include/MolBioLib.hpp"
#include "src/Projects/Repare/repare.hpp"
#include <thread>
#include <numeric>

bool isADouble(string& str)
{
    double temp = 0.0;
    bool result = true;
    try
    {
        temp = stod(str);
        ++temp;   // To avoid non-usage of temp warning.
    }
    catch (invalid_argument& err)
    {
        result = false;
    }
    return result;
}


float diff_func(float curr, float t0) {
  // Using log2f(curr/t0) messed up p-value separations, so not used.
  return (t0 - curr)/t0;
}


void read_counts(string file,
                 unordered_map<string, bool>& all_samples,
                 bool ALLOW_SUPERSET,
                 unordered_map<string, string>& sample2T0,
                 unordered_set<string>& T0s,
                 unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                 unordered_map< string, unordered_map< string, unordered_map<string, float> > >& count,
                 double min_T0, double normalized_count,
                 float max_median_diff) {
  ReadOnlyTSVFile ifp(file, true, true);
  ifp.readRow(0);
  vector<string> header = ifp.tokens;
  // First check that all the samples are in the REPMAP
  for (size_t j = 2; j < header.size(); ++j) {
    if (all_samples.find(header[j]) == all_samples.end()) {
      if (!ALLOW_SUPERSET) {
        cerr << "Error!  In file " << file << ", sample " << header[j] << " does not appear in REPMAP.  If intended, rerun adding ALLOW_SUPERSET=true.  Exiting." << endl;
        ifp.close();
        exit(EXIT_FAILURE);
      } else {
        cerr << "Warning!  In file " << file << ", sample " << header[j] << " does not appear in REPMAP.  Continuing since ALLOW_SUPERSET is true." << endl;
      }
    } else {
      all_samples[header[j]] = true;
    }
  }
  unordered_map<string, size_t> sample2T0index;
  unordered_set<size_t> t0indices;
  for (auto i : sample2T0) {
    size_t target = 0;
    for (size_t j = 2; target == 0 && j < header.size(); ++j) {
      if (header[j] == i.second) target = j;
    }
    sample2T0index[i.first] = target;
  }
  for (size_t j = 2; j < header.size(); ++j) {
    if (T0s.find(header[j]) != T0s.end()) t0indices.insert(j);
  }
  vector<float> total_count(header.size(), 0.0);
  for (size_t i = 1; !ifp.fail(); ++i) {
    ifp.readRow(i);
    for (size_t j = 2; j < ifp.tokens.size(); ++j) {
      if (!isADouble(ifp.tokens[j])) continue;
      float curr_val = convertFromString<float>(ifp.tokens[j]);
      total_count[j] += curr_val;
    }
  }
  ifp.resetToBegin(); 

  // Below is to check if the medians of the samples are similar;
  //   sample    values
  map< string, vector<float> > t0_values, count_values;
  for (size_t i = 1; !ifp.fail(); ++i) {
    ifp.readRow(i);
    string& sgRNA = ifp.tokens[0];
    string& gene = ifp.tokens[1];
    for (size_t j = 2; j < ifp.tokens.size(); ++j) {
      if (t0indices.find(j) != t0indices.end()) continue;
      string& sample = header[j];
      size_t t0index = sample2T0index[sample];
      if (!isADouble(ifp.tokens[t0index])) continue;
      size_t t0count = convertFromString<size_t>(ifp.tokens[t0index]);     
      if (t0count < min_T0) continue;
      if (!isADouble(ifp.tokens[j])) continue;
      float curr_val = convertFromString<float>(ifp.tokens[j]);
      count[gene][sgRNA][sample] = normalized_count*curr_val/total_count[j];
      count_values[sample].push_back(count[gene][sgRNA][sample]);
      if (t0[gene][sgRNA].find(sample) == t0[gene][sgRNA].end()) {
        t0[gene][sgRNA][sample] = normalized_count*t0count/total_count[t0index];
        t0_values[sample].push_back(t0[gene][sgRNA][sample]);
      }
    }
  }
  ifp.close();

  // Compute median values and compare largest with smallest.
  // If the difference is too large, flag.
  Table<string, float> count_medians, t0_medians;
  for (map< string, vector<float> >::iterator i = count_values.begin(); 
       i != count_values.end(); ++i) {
    count_medians.push_back(i->first, vecMedian(i->second));
  }
  sortTable<1>(count_medians);
  if (count_medians.size() > 1) {
    string curr_sample_temp = "";
    float max_med = 0, min_med = 0;
    count_medians.getRow(count_medians.size() - 1, curr_sample_temp, max_med);
    count_medians.getRow(0, curr_sample_temp, min_med);
    float curr_diff = max_med - min_med;
    if (curr_diff > max_median_diff) {
      cerr << "Warning: count median difference in " << file << " is " 
           << curr_diff << " > max = " << max_median_diff << endl;
      cerr << "Samples and their medians after normalization are:" << endl; 
      for (size_t i = 0; i < count_medians.size(); ++i) {
        string curr_sample = "";
        float curr_med = 0;
        count_medians.getRow(i, curr_sample, curr_med);
        cerr << curr_sample << "\t" << curr_med << endl;
      } 
      cerr << endl;
    }
  }
  for (map< string, vector<float> >::iterator i = t0_values.begin(); 
       i != t0_values.end(); ++i) {
    t0_medians.push_back(i->first, vecMedian(i->second));
  }
  sortTable<1>(t0_medians);
  if (t0_medians.size() > 1) {
    string curr_sample_temp = "";
    float max_med = 0, min_med = 0;
    t0_medians.getRow(t0_medians.size() - 1, curr_sample_temp, max_med);
    t0_medians.getRow(0, curr_sample_temp, min_med);
    float curr_diff = max_med - min_med;
    if (curr_diff > max_median_diff) {
      cerr << "Warning: T0 median difference in " << file << " is " 
           << curr_diff << " > max = " << max_median_diff << endl;
      cerr << "Samples and their medians after normalization are:" << endl; 
      for (size_t i = 0; i < t0_medians.size(); ++i) {
        string curr_sample = "";
        float curr_med = 0;
        t0_medians.getRow(i, curr_sample, curr_med);
        cerr << curr_sample << "\t" << curr_med << endl;
      } 
      cerr << endl;
    }
  }
}

void median_normalize_counts(unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                             unordered_map< string, unordered_map< string, unordered_map<string, float> > >& count) {
  // Maps are gene/sgRNA/sample/count.
  // Get all values of t0 and counts by samples, take medians
  unordered_map< string, vector<float> > t0_vals, count_vals;
  for (auto i : t0) {
    for (auto j : i.second) {
      for (auto k : j.second) {
        t0_vals[k.first].push_back(k.second);
      }
    }
  }
  vector<float> medians;
  unordered_map<string, float> sample_count_to_median, sample_t0_to_median;
  for (auto i : t0_vals) {
    sample_t0_to_median[i.first] = vecMedian(i.second);
    medians.push_back(sample_t0_to_median[i.first]);
  }

  for (auto i : count) {
    for (auto j : i.second) {
      for (auto k : j.second) {
        count_vals[k.first].push_back(k.second);
      }
    }
  }
  for (auto i : count_vals) {
    sample_count_to_median[i.first] = vecMedian(i.second);
    medians.push_back(sample_count_to_median[i.first]);
  }

  // Find largest median
  sort(medians.begin(), medians.end());
  float max_median = medians[medians.size() - 1];

  // Add (largest median - sample median) to all values;
  for (unordered_map< string, unordered_map< string, unordered_map<string, float> > >::iterator i = t0.begin(); i != t0.end(); ++i) {
    for (unordered_map< string, unordered_map<string, float> >::iterator j = (i->second).begin(); j != (i->second).end(); ++j) {
      for (unordered_map<string, float>::iterator k = (j->second).begin(); k != (j->second).end(); ++k) {
        string curr_sample = k->first;
        float curr_val = k->second;
        float curr_median = sample_t0_to_median[curr_sample];
        float new_val = curr_val + (max_median - curr_median);
        k->second = new_val;
      }
    }
  }
  for (unordered_map< string, unordered_map< string, unordered_map<string, float> > >::iterator i = count.begin(); i != count.end(); ++i) {
    for (unordered_map< string, unordered_map<string, float> >::iterator j = (i->second).begin(); j != (i->second).end(); ++j) {
      for (unordered_map<string, float>::iterator k = (j->second).begin(); k != (j->second).end(); ++k) {
        string curr_sample = k->first;
        float curr_val = k->second;
        float curr_median = sample_count_to_median[curr_sample];
        float new_val = curr_val + (max_median - curr_median);
        k->second = new_val;
      }
    }
  }
  // Results of normalizing was checked by printing out pre and post 
  // normalization of medians of all samples.
}


void cap_control_ratio(unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                       unordered_map< string, unordered_map< string, unordered_map<string, float> > >& count,
                       unordered_map<string, string>& sample2T0,
                       unordered_set<string>& T0s,
                       unordered_map<string, bool>& is_control,
                       double& max_control_ratio) {
  for (auto i : count) {
    string curr_gene = i.first;
    for (auto j : i.second) {
      string curr_sgrna = j.first;
      for (auto k : j.second) {
        string curr_sample = k.first;
        if (is_control[curr_sample] && T0s.find(curr_sample) == T0s.end()) {
          float curr_count = k.second;
          float curr_t0 = t0[curr_gene][curr_sgrna][curr_sample];
          if (curr_count > (max_control_ratio * curr_t0)) {
            count[curr_gene][curr_sgrna][curr_sample] = max_control_ratio * curr_t0;
          }
        }
      }
    }
  }
}


void counts_corr(unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                 unordered_map< string, unordered_map< string, unordered_map<string, float> > >& count,
                 string correlation_R,
                 string correlation_prefix,
                 double min_output,
                 bool chemo, bool chemo_paired,
                 map<string, string>& test2ctrl) {
  //             gene                   sgRNA                 sample  count 

  set<string> t0s, samples;
  map< string, set<string> > gene2sgRNAs;
  unordered_set<string> all_sgrnas;
  for (auto i : t0) {
    string curr_gene = i.first;
    for (auto j : i.second) {
      string curr_sgrna = j.first;
      all_sgrnas.insert(curr_sgrna);
      gene2sgRNAs[curr_gene].insert(curr_sgrna);
      for (auto k : j.second) {
        t0s.insert(k.first);
      }
    }
  }
  for (auto i : count) {
    string curr_gene = i.first;
    for (auto j : i.second) {
      string curr_sgrna = j.first;
      all_sgrnas.insert(curr_sgrna);
      gene2sgRNAs[curr_gene].insert(curr_sgrna);
      for (auto k : j.second) {
        samples.insert(k.first);
      }
    }
  }
  vector<string> sorted_samples, sorted_t0s, sorted_t0s_out;
  convertSetToVector<string>(t0s, sorted_t0s);
  for (size_t i = 0; i < sorted_t0s.size(); ++i) {
    sorted_t0s_out.push_back("T0 of " + sorted_t0s[i]);
  }
  convertSetToVector<string>(samples, sorted_samples);

  size_t num_sgrnas = all_sgrnas.size(), num_out = 0;
  Ofstream ofp(correlation_prefix + ".normed_counts.tsv");
  Ofstream ofp2(correlation_prefix + ".foldchange.tsv");
  Ofstream ofp2a(correlation_prefix + ".depletion.tsv");
  Ofstream ofp3(correlation_prefix + ".log2foldchange.tsv");
  ofp << "sgRNA\tGene\t" << joinStrings(sorted_t0s_out) 
                         << "\t" << joinStrings(sorted_samples) << endl;
  ofp2 << "sgRNA\tGene\t" << joinStrings(sorted_samples) << endl;
  ofp2a << "sgRNA\tGene\t" << joinStrings(sorted_samples) << endl;
  ofp3 << "sgRNA\tGene\t" << joinStrings(sorted_samples) << endl;
  for (auto i : gene2sgRNAs) {
    string curr_gene = i.first;
    for (auto j : i.second) {
      string curr_sgrna = j;
      vector<string> out_row, out_row2, out_row2a, out_row3;
      bool not_found = false, not_found2 = false;
      for (size_t k = 0; !not_found && k < sorted_t0s.size(); ++k) {
        if (t0[curr_gene][curr_sgrna].find(sorted_t0s[k]) != 
            t0[curr_gene][curr_sgrna].end()) {
          out_row.push_back(convertToString<float>(t0[curr_gene][curr_sgrna][sorted_t0s[k]]));
        } else {
          not_found = true;
        }
      }
      for (size_t k = 0; !not_found && k < sorted_samples.size(); ++k) {
        if (count[curr_gene][curr_sgrna].find(sorted_samples[k]) != 
            count[curr_gene][curr_sgrna].end()) {
          float curr_value = count[curr_gene][curr_sgrna][sorted_samples[k]];
          out_row.push_back(convertToString<float>(curr_value));
          if (t0[curr_gene][curr_sgrna].find(sorted_samples[k]) !=
              t0[curr_gene][curr_sgrna].end()) {
            float curr_t0 = 
              t0[curr_gene][curr_sgrna][sorted_samples[k]];
            if (curr_t0 != 0.0) {
               float curr_frac = curr_value/curr_t0;
               out_row2.push_back(convertToString<float>(curr_frac));
               out_row2a.push_back(convertToString<float>(1.0 - curr_frac));
              if (curr_frac > 0.0) {
                out_row3.push_back(convertToString<float>(log2f(curr_frac)));
              } else {
                float small_num = numeric_limits<float>::min();
                out_row3.push_back(convertToString<float>(log2f(small_num)));
              }
            } else {
              not_found2 = true;
            }
          } else {
            not_found2 = true;
          }
        } else {
          not_found = true;
        }
      }
      if (!not_found) {
        ++num_out;
        ofp << curr_sgrna << "\t" << curr_gene << "\t" 
            << joinStrings(out_row) << endl;
      }
      if (!not_found && !not_found2) {
        ofp2 << curr_sgrna << "\t" << curr_gene << "\t" 
             << joinStrings(out_row2) << endl;
        ofp2a << curr_sgrna << "\t" << curr_gene << "\t" 
              << joinStrings(out_row2a) << endl;
        ofp3 << curr_sgrna << "\t" << curr_gene << "\t" 
             << joinStrings(out_row3) << endl;
      }
    }
  }
  ofp3.close();
  ofp2a.close();
  ofp2.close();
  ofp.close();

  Ofstream ofp6(correlation_prefix + ".normed_counts.all_entries.tsv");
  Ofstream ofp7(correlation_prefix + ".foldchange.all_entries.tsv");
  Ofstream ofp7a(correlation_prefix + ".depletion.all_entries.tsv");
  Ofstream ofp8(correlation_prefix + ".log2foldchange.all_entries.tsv");
  ofp6 << "sgRNA\tGene\t" << joinStrings(sorted_t0s_out) 
                          << "\t" << joinStrings(sorted_samples) << endl;
  ofp7 << "sgRNA\tGene\t" << joinStrings(sorted_samples) << endl;
  ofp7a << "sgRNA\tGene\t" << joinStrings(sorted_samples) << endl;
  ofp8 << "sgRNA\tGene\t" << joinStrings(sorted_samples) << endl;
  for (auto i : gene2sgRNAs) {
    string curr_gene = i.first;
    for (auto j : i.second) {
      string curr_sgrna = j;
      vector<string> out_row, out_row2, out_row2a, out_row3;
      for (size_t k = 0; k < sorted_t0s.size(); ++k) {
        if (t0[curr_gene][curr_sgrna].find(sorted_t0s[k]) != 
            t0[curr_gene][curr_sgrna].end()) {
          out_row.push_back(convertToString<float>(t0[curr_gene][curr_sgrna][sorted_t0s[k]]));
        } else {
          out_row.push_back("");
        }
      }
      for (size_t k = 0; k < sorted_samples.size(); ++k) {
        if (count[curr_gene][curr_sgrna].find(sorted_samples[k]) != 
            count[curr_gene][curr_sgrna].end()) {
          float curr_value = count[curr_gene][curr_sgrna][sorted_samples[k]];
          out_row.push_back(convertToString<float>(curr_value));
          if (t0[curr_gene][curr_sgrna].find(sorted_samples[k]) !=
              t0[curr_gene][curr_sgrna].end()) {
            float curr_t0 = 
              t0[curr_gene][curr_sgrna][sorted_samples[k]];
            if (curr_t0 != 0.0) {
              float curr_frac = curr_value/curr_t0;
              out_row2.push_back(convertToString<float>(curr_frac));
              out_row2a.push_back(convertToString<float>(1.0 - curr_frac));
              if (curr_frac > 0.0) {
                out_row3.push_back(convertToString<float>(log2f(curr_frac)));
              } else {
                float small_num = numeric_limits<float>::min();
                out_row3.push_back(convertToString<float>(log2f(small_num)));
              }
            } else {
              out_row2.push_back("inf");
              out_row2a.push_back("-inf");
              out_row3.push_back("inf");
            }
          } else {
            out_row2.push_back("");
            out_row2a.push_back("");
            out_row3.push_back("");
          }
        } else {
          out_row.push_back("");
          out_row2.push_back("");
          out_row2a.push_back("");
          out_row3.push_back("");
        }
      }
      ofp6 << curr_sgrna << "\t" << curr_gene << "\t" 
           << joinStrings(out_row) << endl;
      ofp7 << curr_sgrna << "\t" << curr_gene << "\t" 
           << joinStrings(out_row2) << endl;
      ofp7a << curr_sgrna << "\t" << curr_gene << "\t" 
            << joinStrings(out_row2a) << endl;
      ofp8 << curr_sgrna << "\t" << curr_gene << "\t" 
           << joinStrings(out_row3) << endl;
    }
  }
  ofp8.close();
  ofp7a.close();
  ofp7.close();
  ofp6.close();

  if (chemo && chemo_paired) {
    ReadOnlyTSVFile ifp2a(correlation_prefix + ".depletion.tsv", true, true);
    ifp2a.readRow(0);
    unordered_map<string, size_t> sample2col;
    for (size_t i = 2; i < ifp2a.tokens.size(); ++i) {
      sample2col[ifp2a.tokens[i]] = i;
    }
    Ofstream ofp4(correlation_prefix + ".depletion_test_minus_ctrl.tsv");
    ofp4 << "sgRNA\tGene";
    for (auto j : test2ctrl) {
      ofp4 << "\t" << j.first << " - " << j.second;
    }
    ofp4 << endl;
    for (size_t i = 1; !ifp2a.fail(); ++i) {
      ifp2a.readRow(i);
      ofp4 << ifp2a.tokens[0] << "\t" << ifp2a.tokens[1];
      for (auto j : test2ctrl) {
        double test_val = convertFromString<double>(ifp2a.tokens[sample2col[j.first]]);
        double ctrl_val = convertFromString<double>(ifp2a.tokens[sample2col[j.second]]);
        ofp4 << "\t" << test_val - ctrl_val;
      }
      ofp4 << endl;
    }
    ofp4.close();
    ifp2a.close();

    ReadOnlyTSVFile ifp3(correlation_prefix + ".log2foldchange.tsv", true, true);
    ifp3.readRow(0);
    sample2col.clear();
    for (size_t i = 2; i < ifp3.tokens.size(); ++i) {
      sample2col[ifp3.tokens[i]] = i;
    }
    Ofstream ofp5(correlation_prefix + ".log2foldchange_test_minus_ctrl.tsv");
    ofp5 << "sgRNA\tGene";
    for (auto j : test2ctrl) {
      ofp5 << "\t" << j.first << " - " << j.second;
    }
    ofp5 << endl;
    for (size_t i = 1; !ifp3.fail(); ++i) {
      ifp3.readRow(i);
      ofp5 << ifp3.tokens[0] << "\t" << ifp3.tokens[1];
      for (auto j : test2ctrl) {
        double test_val = convertFromString<double>(ifp3.tokens[sample2col[j.first]]);
        double ctrl_val = convertFromString<double>(ifp3.tokens[sample2col[j.second]]);
        ofp5 << "\t" << test_val - ctrl_val;
      }
      ofp5 << endl;
    }
    ofp5.close();
    ifp3.close();
  }
  if (static_cast<double>(num_out)/static_cast<double>(num_sgrnas) < min_output) {
    cerr << "Error!  Less than " << 100*min_output << " percent of sgRNAs output.  Did you use this on multiple libraries?  No plot produced." << endl;
  } else {
    system("Rscript --vanilla \"" + correlation_R + "\" \"" + 
            correlation_prefix + ".normed_counts.tsv\" \"" +
            correlation_prefix + ".normed_counts.cor\" 0");
    system("Rscript --vanilla \"" + correlation_R + "\" \"" + 
            correlation_prefix + ".foldchange.tsv\" \"" +
            correlation_prefix + ".foldchange.cor\" 1");
    // Note that the correlation of the depletion
    // is the same as for the foldchange.
  }
}


void detrend_sgrnas(unordered_map<string, string>& sgrna_category,
                    unordered_map<string, size_t>& num_in_category,
                    unordered_map<string, double>& pct_in_category,
                    double& detrend_min_pct,
                    unordered_set<string>& T0s,
                    unordered_map<string, bool>& isControl,
                    unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                    unordered_map< string, unordered_map< string, unordered_map<string, float> > >& detrend_count,
                    unordered_map<string, double>& sgrna_detrend_val_mut,
                    unordered_map<string, double>& sgrna_detrend_val_wt,
                    string detrend_base, string detrend_out) {

  float epsilon = numeric_limits<float>::min();

  // 1. Collect foldchanges per category.
  unordered_map< string, vector<float> > category_fcs_mut, category_fcs_wt;
  for (auto gene : detrend_count) {
    string curr_gene = gene.first;
    for (auto sgRNA : gene.second) {
      string curr_sgRNA = sgRNA.first;
      for (auto sample : sgRNA.second) {
        string curr_sample = sample.first;
        // Do not include T0s in detrending!
        if (T0s.find(curr_sample) != T0s.end()) continue;
        float curr_count =
            static_cast<float>(detrend_count[curr_gene][curr_sgRNA][curr_sample]);
        float curr_t0 =
            static_cast<float>(t0[curr_gene][curr_sgRNA][curr_sample]);
        string curr_cat = sgrna_category[curr_sgRNA];
        if (isControl[curr_sample]) {
          category_fcs_wt[curr_cat].push_back(curr_count/curr_t0);
        } else {
          category_fcs_mut[curr_cat].push_back(curr_count/curr_t0);
        }
      }
    }
  }


  // 2. Compute median foldchange per category.
  unordered_map<string, float> category_median_fcs_mut, category_median_fcs_wt;
  for (auto cat : category_fcs_mut) {
    category_median_fcs_mut[cat.first] = vecMedian(cat.second);
  }
  if (category_median_fcs_mut[detrend_base] == 0.0) {
    cerr << "Error!  Base TEST category's median is 0.  Exiting." << endl;
    exit(EXIT_FAILURE);
  }
  for (auto cat : category_fcs_wt) {
    category_median_fcs_wt[cat.first] = vecMedian(cat.second);
  }
  if (category_median_fcs_wt[detrend_base] == 0.0) {
    cerr << "Error!  Base CTRL category's median is 0.  Exiting." << endl;
    exit(EXIT_FAILURE);
  }


  // 3. Compute detrend
  unordered_map<string, float> category_detrend_mut, category_detrend_wt;
  for (auto cat : category_median_fcs_mut) {
    category_detrend_mut[cat.first] = category_median_fcs_mut[detrend_base]/
                                      category_median_fcs_mut[cat.first];
  }
  for (auto cat : category_median_fcs_wt) {
    category_detrend_wt[cat.first] = category_median_fcs_wt[detrend_base]/
                                     category_median_fcs_wt[cat.first];
  }


  // 4. Output file
  if (detrend_out != "") {
    Ofstream ofp(detrend_out);
    ofp << "Category\tNumber of sgRNAs in category\tPercent of sgRNAs in category\tTEST normalized counts correction multiplicative factor\tTEST log2(correction factor)\tCTRL normalized counts correction multiplicative factor\tCTRL log2(correction factor)\tUsing factor" << endl;
    for (auto cat : category_detrend_mut) {
      string good = (pct_in_category[cat.first] >= detrend_min_pct) ? "Y" : "N";
      ofp << cat.first << "\t"
          << num_in_category[cat.first] << "\t"
          << pct_in_category[cat.first] << "\t"
          << cat.second << "\t"
          << log2(cat.second + epsilon) << "\t"
          << category_detrend_wt[cat.first] << "\t"
          << log2(category_detrend_wt[cat.first] + epsilon) << "\t"
          << good 
          << endl;
    }
    ofp.close();
  }


  // 5. Detrend counts
  for (auto gene : detrend_count) {
    string curr_gene = gene.first;
    for (auto sgRNA : gene.second) {
      string curr_sgRNA = sgRNA.first;
      string curr_category = sgrna_category[curr_sgRNA];
      if (pct_in_category[curr_category] < detrend_min_pct) continue;
      sgrna_detrend_val_mut[curr_sgRNA] = category_detrend_mut[curr_category];
      sgrna_detrend_val_wt[curr_sgRNA] = category_detrend_wt[curr_category];
      for (auto sample : sgRNA.second) {
        string curr_sample = sample.first;
        if (T0s.find(curr_sample) != T0s.end()) continue;
        float curr_count =
            static_cast<float>(detrend_count[curr_gene][curr_sgRNA][curr_sample]);
        float curr_detrend = 1.0;
        if (isControl[curr_sample]) {
          curr_detrend = category_detrend_wt[curr_category];
        } else {
          curr_detrend = category_detrend_mut[curr_category];
        }
        float curr_detrend_count = curr_count * curr_detrend;
        detrend_count[curr_gene][curr_sgRNA][curr_sample] = curr_detrend_count;
      }
    }
  }
}


void thread_output_lines(unordered_map<string, size_t>& gene_thread,
                         size_t thread_id,
                         unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                         unordered_map< string, unordered_map< string, unordered_map<string, float> > >& count,
                         unordered_map<string, bool>& isControl,
                         map<string, string>& test2ctrl,
                         vector<float>& ess_diff, 
                         vector<float>& noness_diff,
                         set<string>& ess,
                         set<string>& ess_ribo,
                         set<string>& noness,
                         set<string>& pr_pos,
                         set<string>& tsg,
                         double LARGE_NEG, size_t PARAMS,
                         bool ADD_RESISTORS, bool ADD_RIBO, bool has_pr_pos,
                         bool ADD_COEFF_VAR,
                         bool CHEMO, bool CHEMO_PAIRED, bool CHEMO_NO_DEPLETION,
                         bool SINGLE_CELL_LINE, bool MIXED_LIBRARIES,
                         bool SGRNA_PAIRED_DIFFERENCES, bool ADD_TOP_TWO,
                         Table<double, string, string, size_t>& output_lines) {
  for (auto gene : count) {
    string curr_gene = gene.first;
    // Only work on the genes for this thread.
    if (gene_thread[curr_gene] != thread_id) continue;
    vector<float> all_test, all_test_t0, test_diff, 
                  all_ctrl, all_ctrl_t0, ctrl_diff;
    //             sgRNA     values
    unordered_map< string, vector<float> > test_vals, ctrl_vals;
    // TO: I know above is redundant to below, but the memory usage and
    //     CPU impact is small and it simplifies the logic further down.  
    //             sgRNA                  sample  value 
    unordered_map< string, unordered_map< string, float> > test_sample_val, ctrl_sample_val;
    set<string> sgRNAs;
    for (auto sgRNA : gene.second) {
      string curr_sgRNA = sgRNA.first;
      sgRNAs.insert(curr_sgRNA);
      for (auto sample : sgRNA.second) {
        string curr_sample = sample.first;
        float curr_count =
            static_cast<float>(count[curr_gene][curr_sgRNA][curr_sample]);
        float curr_t0 =
            static_cast<float>(t0[curr_gene][curr_sgRNA][curr_sample]);
        if (isControl[curr_sample]) {
          ctrl_diff.push_back(diff_func(curr_count, curr_t0));
          ctrl_vals[curr_sgRNA].push_back(diff_func(curr_count, curr_t0));
          ctrl_sample_val[curr_sgRNA][curr_sample] = diff_func(curr_count, curr_t0);
          if (ADD_COEFF_VAR) {
            all_ctrl.push_back(curr_count);
            all_ctrl_t0.push_back(curr_t0);
          }
        } else {
          test_diff.push_back(diff_func(curr_count, curr_t0));
          test_vals[curr_sgRNA].push_back(diff_func(curr_count, curr_t0));
          test_sample_val[curr_sgRNA][curr_sample] = diff_func(curr_count, curr_t0);
          if (ADD_COEFF_VAR) {
            all_test.push_back(curr_count);
            all_test_t0.push_back(curr_t0);
          }
        }
      }
    }
     
    // false = sort from biggest to lowest.
    // SyncSort sorts the keys, e.g. ctrl_diff, upon initialization.
    SyncSort<float> ctrl_sort(ctrl_diff, false);
    double ctrl_q1, ctrl_median, ctrl_q3;
    // Note for the below tableQ123 calls, we reverse the q3 and q1
    // variables, since the input array is sorted in reverse order.
    tableQ123(ctrl_diff, ctrl_median, ctrl_q3, ctrl_q1);
    SyncSort<float> test_sort(test_diff, false);
    double test_q1, test_median, test_q3;
    tableQ123(test_diff, test_median, test_q3, test_q1);
    double both, left, right;
    mannWhitneyU(test_diff, ctrl_diff, both, left, right);

    // Compute difference between test and control below
    double median_diff, q3_diff, q1_diff;
    vector<float> all_diffs, median_diffs, q3_diffs, q1_diffs;
    double fraction_pos_diff = -1.0;

    double med_depl_top_1 = LARGE_NEG, med_depl_top_2 = LARGE_NEG,
           med_diff_top_1 = LARGE_NEG, med_diff_top_2 = LARGE_NEG,
           sum_med_diff_top_1_and_2 = LARGE_NEG;
    if (ADD_TOP_TWO) {
      // Some repetition between this and SGRNA_PAIRED_DIFFERENCES, but
      // clarity over efficiency as the hit on CPU will be low.
      //    t_depl  diff
      Table<double, double> top_data;
      for (auto test_sgrnas : test_vals) {
        string curr_test_sgrna = test_sgrnas.first;
        if (ctrl_vals.find(curr_test_sgrna) != ctrl_vals.end()) {
          vector<float> curr_test_vals = test_sgrnas.second;
          vector<float> curr_ctrl_vals = ctrl_vals[curr_test_sgrna];
          sort(curr_test_vals.begin(), curr_test_vals.end()); 
          sort(curr_ctrl_vals.begin(), curr_ctrl_vals.end()); 
          // Below are median depletion for this sgRNA.
          double temp_q1, temp_q3, test_depl, ctrl_depl;
          tableQ123(curr_test_vals, test_depl, temp_q1, temp_q3);
          tableQ123(curr_ctrl_vals, ctrl_depl, temp_q1, temp_q3);
          double delta_depl = test_depl - ctrl_depl;
          top_data.push_back(test_depl, delta_depl);   
        }
      }
      sortTable<1>(top_data, true);
      if (top_data.size() > 0) {
        top_data.getRow(0, med_depl_top_1, med_diff_top_1);
        sum_med_diff_top_1_and_2 = med_diff_top_1;
      }
      if (top_data.size() > 1) {
        top_data.getRow(1, med_depl_top_2, med_diff_top_2);
        sum_med_diff_top_1_and_2 += med_diff_top_2;
      }
    }
    
    // The below only works when the libraries aren't mixed.
    if (SGRNA_PAIRED_DIFFERENCES && !(CHEMO && CHEMO_PAIRED)) {
      for (auto test_sgrnas : test_vals) {
        string curr_test_sgrna = test_sgrnas.first;
        if (ctrl_vals.find(curr_test_sgrna) != ctrl_vals.end()) {
          vector<float> curr_test_vals = test_sgrnas.second;
          vector<float> curr_ctrl_vals = ctrl_vals[curr_test_sgrna];
          sort(curr_test_vals.begin(), curr_test_vals.end()); 
          double curr_test_q1, curr_test_q2, curr_test_q3;
          tableQ123(curr_test_vals, curr_test_q2, curr_test_q1, curr_test_q3);
          sort(curr_ctrl_vals.begin(), curr_ctrl_vals.end()); 
          double curr_ctrl_q1, curr_ctrl_q2, curr_ctrl_q3;
          tableQ123(curr_ctrl_vals, curr_ctrl_q2, curr_ctrl_q1, curr_ctrl_q3);
          median_diffs.push_back(curr_test_q2 - curr_ctrl_q2);
          q3_diffs.push_back(curr_test_q3 - curr_ctrl_q3);
          q1_diffs.push_back(curr_test_q1 - curr_ctrl_q1);
        }
      } 
      double dummy1, dummy2;   // dummy variables never used.
      sort(median_diffs.begin(), median_diffs.end());
      tableQ123(median_diffs, median_diff, dummy1, dummy2);
      sort(q3_diffs.begin(), q3_diffs.end());
      tableQ123(q3_diffs, dummy1, dummy2, q3_diff);
      sort(q1_diffs.begin(), q1_diffs.end());
      tableQ123(q1_diffs, dummy1, q1_diff, dummy2);
    } else if (!MIXED_LIBRARIES || (CHEMO && CHEMO_PAIRED)) {
      if (!MIXED_LIBRARIES) {
        // Do all pairs difference.
        for (auto test_sgrnas : test_vals) {
          string curr_test_sgrna = test_sgrnas.first;
          vector<float>& curr_test_vals = test_sgrnas.second;
          if (ctrl_vals.find(curr_test_sgrna) != ctrl_vals.end()) {
            vector<float>& curr_ctrl_vals = ctrl_vals[curr_test_sgrna];
            for (size_t i = 0; i < curr_test_vals.size(); ++i) {
              for (size_t j = 0; j < curr_ctrl_vals.size(); ++j) {
                all_diffs.push_back(curr_test_vals[i] - curr_ctrl_vals[j]);
              }
            }
          }
        } 
      } else if (CHEMO && CHEMO_PAIRED) {
        // CHEMO and CHEMO_PAIRED
        // Here, just do the paired samples.
        double num_pos_diff = 0.0;
        for (auto test_sgrnas : test_sample_val) {
          string curr_test_sgrna = test_sgrnas.first;
          for (auto test_sample : test_sgrnas.second) {
            string curr_test_sample = test_sample.first;
            if (test2ctrl.find(curr_test_sample) != test2ctrl.end()) {
              string curr_ctrl_sample = test2ctrl[curr_test_sample];
              if (ctrl_sample_val[curr_test_sgrna].find(curr_ctrl_sample) != ctrl_sample_val[curr_test_sgrna].end()) {
                
                float curr_val = test_sample.second - ctrl_sample_val[curr_test_sgrna][curr_ctrl_sample];
                all_diffs.push_back(curr_val);
                if (curr_val >= 0.0) {
                  num_pos_diff += 1.0;
                }
              }
            }
          } 
        }
        fraction_pos_diff = num_pos_diff/static_cast<double>(all_diffs.size());
      }
      // Because we're sorting forward here, in tableQ123 
      // we have q1_diff then q3_diff, unlike the rest of the code.
      SyncSort<float> all_diffs_sort(all_diffs);
      tableQ123(all_diffs, median_diff, q1_diff, q3_diff);
    } else {
      // MIXED_LIBRARIES (default true) and not CHEMO and CHEMO_PAIRED
      median_diff = test_median - ctrl_median;
      q3_diff = test_q3 - ctrl_q3;
      q1_diff = test_q1 - ctrl_q1;
    }


    // For the gene to be essential in test, test against non-essential 
    // and want nright <= 0.0001;
    double nboth, nleft, nright;
    if (!CHEMO_NO_DEPLETION) {
      mannWhitneyU(test_diff, noness_diff, nboth, nleft, nright);
    } else {
      nboth = -2.0;
      nleft = -2.0;
      nright = -2.0;
    }

    // For the gene to be nonessential in control, test against essential and
    // want eleft <= 0.0001;
    double eboth, eleft, eright;
    if (!CHEMO) {
      mannWhitneyU(ctrl_diff, ess_diff, eboth, eleft, eright);
    } else {
      eboth = -2.0;
      eleft = -2.0;
      eright = -2.0;
    }

    float median_test = 0.0, median_ctrl = 0.0, 
          median_test_t0 = 0.0, median_ctrl_t0 = 0.0; 
    string str_coeff_var_test = "", str_coeff_var_ctrl = "",
           str_coeff_var_diff = "",
           str_var_test = "", str_var_ctrl = "",
           str_var_test_t0 = "", str_var_ctrl_t0 = ""; 
    if (ADD_COEFF_VAR) {
      if (test_q3 + test_q1 != 0.0) { 
        float coeff_var_test = (test_q3 - test_q1)/(test_q3 + test_q1);
        str_coeff_var_test = convertToString<float>(coeff_var_test);
      }
      if (ctrl_q3 + ctrl_q1 != 0.0) { 
        float coeff_var_ctrl = (ctrl_q3 - ctrl_q1)/(ctrl_q3 + ctrl_q1);
        str_coeff_var_ctrl = convertToString<float>(coeff_var_ctrl);
      }
      double coeff_var_diff_median, coeff_var_diff_q1, coeff_var_diff_q3;
      if (MIXED_LIBRARIES && !CHEMO) {
        // Need to populate all_diffs vector anyway to get this statistic.
        // Do all pairs difference.
        for (auto test_sgrnas : test_vals) {
          string curr_test_sgrna = test_sgrnas.first;
          vector<float>& curr_test_vals = test_sgrnas.second;
          if (ctrl_vals.find(curr_test_sgrna) != ctrl_vals.end()) {
            vector<float>& curr_ctrl_vals = ctrl_vals[curr_test_sgrna];
            for (size_t i = 0; i < curr_test_vals.size(); ++i) {
              for (size_t j = 0; j < curr_ctrl_vals.size(); ++j) {
                all_diffs.push_back(curr_test_vals[i] - curr_ctrl_vals[j]);
              }
            }
          }
        } 
        SyncSort<float> all_diffs_sort(all_diffs);
        tableQ123(all_diffs, coeff_var_diff_median,
                  coeff_var_diff_q1, coeff_var_diff_q3);
      } else {
        coeff_var_diff_median = median_diff;
        coeff_var_diff_q1 = q1_diff;
        coeff_var_diff_q3 = q3_diff;
      }
      if (coeff_var_diff_q3 + coeff_var_diff_q1 != 0.0) {
        float coeff_var_diff_val = (coeff_var_diff_q3 - coeff_var_diff_q1)/
                                   (coeff_var_diff_q3 + coeff_var_diff_q1);
        str_coeff_var_diff = convertToString<float>(coeff_var_diff_val);
      }

      double temp_q1, temp_median, temp_q3;
      SyncSort<float> all_test_sort(all_test, false);
      tableQ123(all_test, temp_median, temp_q3, temp_q1);
      median_test = temp_median;
      if (temp_q3 + temp_q1 != 0.0) {
        float var_test = (temp_q3 - temp_q1)/(temp_q3 + temp_q1);
        str_var_test = convertToString<float>(var_test);
      }

      SyncSort<float> all_test_t0_sort(all_test_t0, false);
      tableQ123(all_test_t0, temp_median, temp_q3, temp_q1);
      median_test_t0 = temp_median;
      if (temp_q3 + temp_q1 != 0.0) {
        float var_test_t0 = (temp_q3 - temp_q1)/(temp_q3 + temp_q1);
        str_var_test_t0 = convertToString<float>(var_test_t0);
      }
     
      SyncSort<float> all_ctrl_sort(all_ctrl, false);
      tableQ123(all_ctrl, temp_median, temp_q3, temp_q1);
      median_ctrl = temp_median;
      if (temp_q3 + temp_q1 != 0.0) {
        float var_ctrl = (temp_q3 - temp_q1)/(temp_q3 + temp_q1);
        str_var_ctrl = convertToString<float>(var_ctrl);
      }

      SyncSort<float> all_ctrl_t0_sort(all_ctrl_t0, false);
      tableQ123(all_ctrl_t0, temp_median, temp_q3, temp_q1);
      median_ctrl_t0 = temp_median;
      if (temp_q3 + temp_q1 != 0.0) {
        float var_ctrl_t0 = (temp_q3 - temp_q1)/(temp_q3 + temp_q1);
        str_var_ctrl_t0 = convertToString<float>(var_ctrl_t0);
      }
    }

    double rboth = 0.0, rleft = 0.0, rright = 0.0;
    if (ADD_RESISTORS) {
      // For a *resistor* gene to be nonessential in test, test against 
      // essential and want rleft <= 0.0001;
      mannWhitneyU(test_diff, ess_diff, rboth, rleft, rright);
    }

    string essential = "";
    if (ess.find(curr_gene) != ess.end()) {
      essential = "Essential";
    } else if (noness.find(curr_gene) != noness.end()) {
      essential = "Nonessential";
    }
    string ribosomal = "";
    if (ADD_RIBO && ess_ribo.find(curr_gene) != ess_ribo.end()) {
      ribosomal = "Essential ribosomal";
    }
    string pr_pos_string = "";
    if (has_pr_pos && pr_pos.find(curr_gene) != pr_pos.end()) {
      pr_pos_string = "Y";
    }
    // Only print things with a good p-value.
    if (right == -2) continue;

    double correction = (test_median < 0.0) ? LARGE_NEG : 0.0;

    // Fres = F for resistors.  In V2, F=0, so we need something
    // to hold values of F for resistors.
    double A, B, C, D, E, F, G, H, Fres, K;
    // double A, B, C, D, E, F, G, H, II, J;
    if (PARAMS == 1 && !CHEMO) {
      if (SINGLE_CELL_LINE) {
        // From v27
        A = 8.596739653;
        B = 9.344152324;
        C = 0.584706108;
        D = 2.805848566;
        E = 2.781028203;
        F = 0.206275346;
        Fres = F;
        G = 3.294868408;
        H = 2.88214286;
        K = 0.0;
        // II = 0.0; J = 0.0;
      } else {
        // From v27
        A = 8.596739653;
        B = 9.344152324;
        C = 0.584706108;
        D = 2.805848566;
        E = 2.781028203;
        F = 0.206275346;
        Fres = F;
        G = 3.294868408;
        H = 2.88214286;
        K = 0.0;
        // II = 0.0; J = 0.0;
      }
    } else if (PARAMS == 1 && CHEMO) {
      if (!CHEMO_NO_DEPLETION) {
        if (SINGLE_CELL_LINE) {
          // From v27 
          A = 1.041577465;
          B = 4.111248228;
          C = 2.06410473;
          D = 1.680487337;
          E = 0.679885697;
          F = 0.0;
          Fres = 0.206275346;   // Use F of isogenic single cell line
          G = 1.365148475;
          H = 1.534659974;
          K = 0.0;
          // II = 0.0; J = 0.0;
       } else {
          // From v27 
          A = 1.041577465;
          B = 4.111248228;
          C = 2.06410473;
          D = 1.680487337;
          E = 0.679885697;
          F = 0.0;
          Fres = 0.206275346;   // Use F of isogenic combined cell line
          G = 1.365148475;
          H = 1.534659974;
          K = 0.0;
          // II = 0.0; J = 0.0;
       }
      } else {
        // CHEMO_NO_DEPLETION
        if (SINGLE_CELL_LINE) {
          A = 0.0;
          B = 0.0;
          C = 2.06410473;
          D = 1.680487337;
          E = 0.0;
          F = 0.0;
          Fres = 0.206275346;   // Use F of isogenic single cell line
          G = 1.365148475;
          H = 1.534659974;
          K = 0.0;
          // II = 0.0; J = 0.0;
        } else {
          A = 0.0;
          B = 0.0;
          C = 2.06410473;
          D = 1.680487337;
          E = 0.0;
          F = 0.0;
          Fres = 0.206275346;   // Use F of isogenic combined cell line
          G = 1.365148475;
          H = 1.534659974;
          K = 0.0;
          // II = 0.0; J = 0.0;
        }
      }
    } else if (PARAMS == 2) {
      // All PARAMS == 2 values are from v2 of version2
      // Directory 20190314_optimize_CCA_params_round_2_v2
      if (CHEMO) {
        if (CHEMO_NO_DEPLETION) {
          if (SINGLE_CELL_LINE) {
            A = 0.0;
            B = 0.0;
            C = 3.468453878;
            D = 1.0;
            E = 0.0;
            F = 0.0;
            Fres = 0.354296361;   // Use F of single isogenic for resistors.
            G = 0.010789329;
            H = 0.013745531;
            K = 0.0;
          } else {
            // Combined
            A = 0.0;
            B = 0.0;
            C = 5.004139989;
            D = 1.0;
            E = 0.0;
            F = 0.0;
            Fres = 0.180033491;   // Use F of combined isogenic for resistors.
            G = 0.010700926;
            H = 0.010827516;
            K = 0.0;
          }
        } else {
          // Take into account depletion
          if (SINGLE_CELL_LINE) {
            A = 0.485518506;
            B = 3.96557834;
            C = 0.799534557;
            D = 1.0;
            E = 0.289668341;
            F = 0;
            Fres = 0.354296361;   // Use F of single isogenic for resistors.
            G = 1.34481421;
            H = 0.580980998;
            K = 0.0;   // XXX Probably change this
          } else {
            // Combined cell lines
            A = 5.340589607;
            B = 5.115177494;
            C = 9.704818376;
            D = 1.0;
            E = 0.169660711;
            F = 0.0;
            Fres = 0.180033491;   // Use F of combined isogenic for resistors.
            G = 0.918500145;
            H = 0.058745284;
            K = 0.0;   // XXX Probably change this
          }
        }
      } else {
        // Isogenic
        if (SINGLE_CELL_LINE) {
          A = 2.015163396;
          B = 0.017303765;
          C = 0.020330701;
          D = 1.0;
          E = 8.763376388;
          F = 0.354296361;
          Fres = F;
          G = 7.062717517;
          H = 0.221024716;
          K = 0.0;
        } else {
          // Combined
          A = 5.637716436;
          B = 0.549799684;
          C = 1.198982817;
          D = 1.0;
          E = 7.654839527;
          F = 0.180033491;
          Fres = F;
          G = 7.698926805;
          H = 0.341913052;
          K = 0.0;
        }
      }

    } else if (PARAMS == 3) {
      // All PARAMS == 3 values are from 
      // directory 20190626_optimize_CCA_params_round_4_new_diffs
      // Only for SGRNA_PAIRED_DIFFERENCES = true mode.
      if (!CHEMO) {
        // Isogenic
        if (SINGLE_CELL_LINE) {
          A = 5.267169707;
          B = 8.364181233;
          C = 2.268950324;
          D = 1.0;
          E = 5.30206221;
          F = 0.005299121;
          Fres = F;
          G = 5.06729317;
          H = 0.996932512;
          K = 0.0;
        } else {
          // Combined
          A = 8.936605889;
          B = 2.499972998;
          C = 5.757248283;
          D = 1.0;
          E = 4.826243161;
          F = 0.437423772;
          Fres = F;
          G = 0.128810004;
          H = 0.377570551;
          K = 0.0;
        }
      } else {
        cerr << "Error!  PARAMS=3 but CHEMO=true.  Exiting" << endl;
        exit(EXIT_FAILURE);
      }
    } else if (PARAMS == 4) {
      // All PARAMS == 4 values are from 
      // directory 20190626_optimize_CCA_params_round_4_new_diffs
      // Only for SGRNA_PAIRED_DIFFERENCES = true mode.
      if (!CHEMO) {
        // Isogenic
        A = 1.618701972;
        B = 7.070265348;
        C = 1.441648997;
        D = 1.0;
        E = 5.99515171;
        F = 0.252835296;
        Fres = F;
        G = 0.084206047;
        H = 6.668960692;
        K = 0.0;   // XXX Probably change this
      } else {
        cerr << "Error!  PARAMS=4 but CHEMO=true.  Exiting" << endl;
        exit(EXIT_FAILURE);
      }
    } else {
      // Original parameters here. 
      // Originally used for both isogenic and chemogenomic screens.
      A = 3.0; B = 1.5; C = 0.5; D = 1.0; 
      E = 0.01; F = 0.01; G = 0.25; H = 0.25;
      Fres = F;
      K = 0.0;
      // II = 0.0; J = 0.0;
    }
    double score = A*test_median + B*test_q3 + C*median_diff + D*q3_diff;
    // double score = A*test_median + B*test_q3 + C*median_diff + D*q3_diff +
                     // II*ctrl_median + J*ctrl_q3;
    if (CHEMO && CHEMO_PAIRED) {
      score += K*fraction_pos_diff;
    }
#if 0
    double v_num_sgrnas = min(4.0, static_cast<double>(sgRNAs.size()));
    double rate = 0.1375;  // last 0.1, then 0.2, then 0.15 (bad gene too high), then 0.1475 (bad gene too high), then 0.125 (too far down), then 0.135 (too far down), 0.14 (too high)
    v_num_sgrnas = pow(v_num_sgrnas, rate);
    G = G/(v_num_sgrnas/(pow(4.0, rate)));
    H = H/(v_num_sgrnas/(pow(4.0, rate)));
    E = E/(v_num_sgrnas/(pow(4.0, rate)));
    F = F/(v_num_sgrnas/(pow(4.0, rate)));
#endif
    if (both >= 0.0) score = score*(1 - pow(both, G));
    if (right >= 0.0) score = score*(1 - pow(right, H));
    // Yes, nright is associated with E and eleft with F, for historical reason.
    if (CHEMO && !CHEMO_NO_DEPLETION && nright >= 0.0) score = score*(1 - pow(nright, E));
    if (!CHEMO && eleft >= 0.0) score = score*(1 - pow(eleft, F));
    // Only add below if CHEMO mode is false.
    if (!CHEMO) score += correction;
    double resistor_score = 0.0;
    string tsg_string = "";
    if (ADD_RESISTORS) {
      if (tsg.find(curr_gene) != tsg.end()) {
        tsg_string = "Y";
      }
      // Correction like above likely not hold for genes like CTC1 and TEN1 in
      // BRCA1 -/- lesions.
      resistor_score = (-C*median_diff - D*q1_diff);
      if (both >=- 0.0) resistor_score *= 1.0 - pow(both, G);
      if (left >=- 0.0) resistor_score *= 1.0 - pow(left, H);
      if (rleft >=- 0.0) resistor_score *= 1.0 - pow(rleft, Fres);
    }
    string line = curr_gene + "\t" +
                  essential + "\t";
    if (ADD_RIBO) {
      line += ribosomal + "\t";
    }
    if (has_pr_pos) {
      line += pr_pos_string + "\t";
    }
    if (ADD_RESISTORS) {
      line += tsg_string + "\t";
    }
    line += convertToString<size_t>(sgRNAs.size()) + "\t";
    if (ADD_RESISTORS) {
      line += convertToString<double>(test_q1) + "\t" +
              convertToString<double>(ctrl_q1) + "\t" +
              convertToString<double>(q1_diff) + "\t";
    }
    line +=  convertToString<double>(test_median) + "\t" +
             convertToString<double>(ctrl_median) + "\t" +
             convertToString<double>(median_diff) + "\t" +
             convertToString<double>(test_q3) + "\t" +
             convertToString<double>(ctrl_q3) + "\t" +
             convertToString<double>(q3_diff) + "\t" +
             convertToString<double>(both) + "\t" +
             convertToString<double>(left) + "\t" +
             convertToString<double>(right) + "\t" +
             convertToString<double>(nright) + "\t" +
             convertToString<double>(eleft) + "\t" +
             convertToString<double>(nright + eleft) + "\t";
    if (ADD_TOP_TWO) {
      string mde1 = "", mdi1 = "", mde2 = "", mdi2 = "", smdt12 = "";
      if (med_depl_top_1 != LARGE_NEG) {
        mde1 = convertToString<double>(med_depl_top_1);
        mdi1 = convertToString<double>(med_diff_top_1);
        smdt12 = convertToString<double>(sum_med_diff_top_1_and_2);
      }
      if (med_depl_top_2 != LARGE_NEG) {
        mde2 = convertToString<double>(med_depl_top_2);
        mdi2 = convertToString<double>(med_diff_top_2);
      }
      line += mde1 + "\t" + mdi1 + "\t" + 
              mde2 + "\t" + mdi2 + "\t" + smdt12 + "\t";
    }
    if (CHEMO && CHEMO_PAIRED) {
      line += convertToString<double>(fraction_pos_diff) + "\t";
    }
    if (ADD_COEFF_VAR) {
      line += str_coeff_var_test + "\t" +
              str_coeff_var_ctrl + "\t" +
              str_coeff_var_diff + "\t" +
              convertToString<float>(median_test) + "\t" +
              str_var_test + "\t" +
              convertToString<float>(median_test_t0) + "\t" +
              str_var_test_t0 + "\t" +
              convertToString<float>(median_ctrl) + "\t" +
              str_var_ctrl + "\t" +
              convertToString<float>(median_ctrl_t0) + "\t" +
              str_var_ctrl_t0 + "\t";
    }
    if (ADD_RESISTORS) {
      line += convertToString<double>(rleft) + "\t";
      line += convertToString<double>(resistor_score) + "\t";
    } 
    if (!CHEMO) {
      line += convertToString<double>(score - correction) + "\t0\t";
    }
    line += convertToString<double>(score);
    output_lines.push_back(score, curr_gene, line, 0);
  }
}


void compute_output_lines(set<string>& ess,
                          set<string>& ess_ribo,
                          set<string>& noness,
                          set<string>& pr_pos,
                          set<string>& tsg,
                          unordered_map<string, bool>& isControl,
                          map<string, string>& test2ctrl,
                          unordered_map< string, unordered_map< string, unordered_map<string, float> > >& t0,
                          unordered_map< string, unordered_map< string, unordered_map<string, float> > >& count,
                          Table<double, string, string, size_t>& output_lines,
                          double LARGE_NEG, size_t PARAMS,
                          bool ADD_RESISTORS, bool ADD_RIBO, bool has_pr_pos,
                          bool ADD_COEFF_VAR, 
                          bool CHEMO, bool CHEMO_PAIRED, bool CHEMO_NO_DEPLETION,
                          bool SINGLE_CELL_LINE, bool MIXED_LIBRARIES,
                          bool SGRNA_PAIRED_DIFFERENCES, bool ADD_TOP_TWO,
                          size_t NUM_THREADS) {
  vector<float> ess_diff, noness_diff;
  for (auto gene : ess) {
    string curr_gene = gene;
    auto& forGene = count[curr_gene];
    for (auto sgRNA : forGene) {
      string curr_sgRNA = sgRNA.first;
      for (auto sample : sgRNA.second) {
        string curr_sample = sample.first;
        float curr_count =
            static_cast<float>(count[curr_gene][curr_sgRNA][curr_sample]);
        float curr_t0 =
            static_cast<float>(t0[curr_gene][curr_sgRNA][curr_sample]);
        ess_diff.push_back(diff_func(curr_count, curr_t0));
      }
    }
  }
  sort(ess_diff.begin(), ess_diff.end());
  for (auto gene : noness) {
    string curr_gene = gene;
    auto& forGene = count[curr_gene];
    for (auto sgRNA : forGene) {
      string curr_sgRNA = sgRNA.first;
      for (auto sample : sgRNA.second) {
        string curr_sample = sample.first;
        float curr_count =
            static_cast<float>(count[curr_gene][curr_sgRNA][curr_sample]);
        float curr_t0 =
            static_cast<float>(t0[curr_gene][curr_sgRNA][curr_sample]);
        noness_diff.push_back(diff_func(curr_count, curr_t0));
      }
    }
  }
  sort(noness_diff.begin(), noness_diff.end());

  unordered_map<string, size_t> gene_thread;
  size_t thread_count = 0;
  for (auto gene : count) {
    gene_thread[gene.first] = thread_count;
    ++thread_count;
    if (thread_count >= NUM_THREADS) {
      thread_count = 0;
    }
  }

  vector< Table<double, string, string, size_t> > thread_outputs(NUM_THREADS);

  vector<thread> threads(NUM_THREADS);
  for (size_t i = 0; i < NUM_THREADS; ++i) {
    threads[i] = thread(ref(thread_output_lines),
                        ref(gene_thread), i, 
                        ref(t0), ref(count), 
                        ref(isControl), ref(test2ctrl),
                        ref(ess_diff), ref(noness_diff), 
                        ref(ess), ref(ess_ribo), ref(noness), 
                        ref(pr_pos), ref(tsg),
                        LARGE_NEG, PARAMS,
                        ADD_RESISTORS, ADD_RIBO, has_pr_pos, ADD_COEFF_VAR,
                        CHEMO, CHEMO_PAIRED, CHEMO_NO_DEPLETION,
                        SINGLE_CELL_LINE, MIXED_LIBRARIES,
                        SGRNA_PAIRED_DIFFERENCES, ADD_TOP_TWO,
                        ref(thread_outputs[i]));
  }
  for (size_t i = 0; i < NUM_THREADS; ++i) {
    threads[i].join();
  }
  for (size_t i = 0; i < NUM_THREADS; ++i) {
    Table<double, string, string, size_t>& curr_table = thread_outputs[i];
    tuple<double, string, string, size_t> row;
    for (size_t j = 0; j < curr_table.size(); ++j) {
      curr_table.getRow(j, row);
      output_lines.push_back(row);
    }
  }

  // true below = sort in reverse order.  Largest score first.
  sortTable<0>(output_lines, true);
  for (size_t i = 0; i < output_lines.size(); ++i) {
    // This adds the rank.
    get<3>(output_lines.theData)[i] = i+1;
  }
}



int main(int argc, char* argv[])
{
  BeginCommandArguments
  CommandArgumentBoolDefaultDoc("AT_HOME", AT_HOME, true, "If true, then all files are assumed to be in /home, so /home/ is prepended to the input and output.  For Docker version only.");
  CommandArgumentStringDoc("COUNTS", COUNTS, "CSV of files, each of form sgRNA[tab]Gene[tab]Sample1[tab]Sample2[tab]...");
  CommandArgumentStringDoc("REPMAP", REPMAP, "CSV of files, each of form Replicate[tab]Sample (CTRL or TEST)[tab]T0 sample name - if Replicate is T0, use same name here {for CHEMO_PAIRED [tab]CTRL sample name for TEST}.");
  CommandArgumentBoolDefaultDoc("ALLOW_SUPERSET", ALLOW_SUPERSET, false, "If true and if there are samples in COUNTS that aren't in REPMAP, allow it.");
  CommandArgumentSize_tDefaultDoc("PARAMS", PARAMS, 2, "Param set.");
  CommandArgumentBoolDefaultDoc("SINGLE_CELL_LINE", SINGLE_CELL_LINE, true, "If true, use parameters for single cell line runs, else use parameters for mixed cell lines.");
  CommandArgumentBoolDefaultDoc("CHEMO", CHEMO, false, "If true, use chemogenomic mode.  Currently, this means a different set of parameters, ignoring the check of whether control is essential, ignoring if the median depletion is negative.");
  CommandArgumentBoolDefaultDoc("CHEMO_PAIRED", CHEMO_PAIRED, true, "If true and CHEMO is true, use paired information in REPMAP to get the median and Q3 difference in depletion.");
  CommandArgumentBoolDefaultDoc("CHEMO_NO_DEPLETION", CHEMO_NO_DEPLETION, false, "If true, do not use test depletion or check if test is nonessential in this mode.");
  CommandArgumentBoolDefaultDoc("RESISTOR_SCREEN", RESISTOR_SCREEN, false, "If true, the COUNTS and REPMAP file is for a resistor screen.   Thus, a full CCA analysis is not done, but rater a resistor analysis is done.");
  CommandArgumentDoubleDefaultDoc("MAX_CONTROL_RATIO", MAX_CONTROL_RATIO, 1.0, "If > 0, then cap the input control count such that the non-T0 controls are at most this time as much as the T0 count.");
  CommandArgumentFileOutDefaultDoc("RESISTOR_OUT", RESISTOR_OUT, "", "If RESISTOR_SCREEN is true, then output the resistor analysis to this file.   Not using OUT to avoid confusion.");
  CommandArgumentBoolDefaultDoc("ADD_COEFF_VAR", ADD_COEFF_VAR, true, "If true, then add the coefficient of variation and distance of the mean coefficient of variation.  This is to find genes whose depletions vary a lot, when the distance is positively large.");
  CommandArgumentDoubleDefaultDoc("MAX_MEDIAN_DIFF", MAX_MEDIAN_DIFF, 30.0, "If the median of the normalized counts differs by more than this, a warning message is produced.");
  CommandArgumentBoolDefaultDoc("MEDIAN_NORMALIZATION", MEDIAN_NORMALIZATION, false, "If true, normalize all counts by the maximum median taken over all samples.");
  CommandArgumentFileInDefaultDoc("CORRELATION_R", CORRELATION_R, "/root/src/Projects/Repare/cca_sgrna_gene_correl.R", "R script to compute correlation and plot.");
  CommandArgumentStringDefaultDoc("CORRELATION_PREFIX", CORRELATION_PREFIX, "", "If non-empty, output normalized counts correlation matrix and plot.");
  CommandArgumentBoolDefaultDoc("STOP_AFTER_CORRELATION", STOP_AFTER_CORRELATION, false, "If true, stop after initial QC and correlation plots.  No analysis is done.  For obtaining foldchange matrices, etc.");
  CommandArgumentStringDefaultDoc("BARPLOT_HEAD", BARPLOT_HEAD, "", "If specified and if CORRELATION_PREFIX is specified, then barplot BARPLOT_TOP_N genes by CCA score.");
  CommandArgumentBoolDefaultDoc("BARPLOT_ADD_RANK", BARPLOT_ADD_RANK, true, "If true, add rank to the output file name before the gene name.");
  CommandArgumentSize_tDefaultDoc("BARPLOT_TOP_N", BARPLOT_TOP_N, 25, "Top number of genes to bar plot.");
  CommandArgumentDoubleDefaultDoc("BARPLOT_MIN", BARPLOT_MIN, 0.0, "Start of barplot y-axis.");
  CommandArgumentDoubleDefaultDoc("BARPLOT_MAX", BARPLOT_MAX, 3.0, "End of barplot y-axis.");
  CommandArgumentDoubleDefaultDoc("BARPLOT_BY", BARPLOT_BY, 0.25, "Tick spacing on barplot y-axis.");
  CommandArgumentBoolDefaultDoc("BARPLOT_NEG_LOG2", BARPLOT_NEG_LOG2, false, "If true, plot -log2(foldchange).");
  CommandArgumentDoubleDefaultDoc("BARPLOT_EPSILON", BARPLOT_EPSILON, numeric_limits<double>::min(), "Small value to pass to bar plot when plotting -log2(foldchange).");
  CommandArgumentBoolDefaultDoc("ADD_Z_SCORES", ADD_Z_SCORES, true, "If true, add Z-scores to output.");
  CommandArgumentDoubleDefaultDoc("Z_SCORES_MIN", Z_SCORES_MIN, 1.0, "Add this much to numerator and denominator of log-foldchange ratio between test and control.");
  CommandArgumentSize_tDefaultDoc("Z_SCORES_HALF_WINDOW", Z_SCORES_HALF_WINDOW, 0, "If specified, use as half size of the window used to calculate the deviation.  Otherwise, computed from minimizing the distance of the D'Agostino score from 2.");
  CommandArgumentSize_tDefaultDoc("Z_SCORES_WINDOW_MIN", Z_SCORES_WINDOW_MIN, 100, "If specified, then use this as the minimum window size.");
  CommandArgumentDoubleDefaultDoc("Z_SCORES_L", Z_SCORES_L, 0.0, "Multiply the Z-score by this amount and add to the score to get the final score and rank.");
  CommandArgumentDoubleDefaultDoc("MIN_CORRELATION_OUTPUT", MIN_CORRELATION_OUTPUT, (2.0/3.0), "Ratio of sgRNAs that must appear in the output for a plot to be done.  Default is two-thirds.");
  CommandArgumentFileInDefaultDoc("DETREND", DETREND, "", "If non-empty, the file contains information about the alignment properties of the sgRNAs.   Specify the columns to use in DETREND_COLUMNS and minimum percentages DETREND_MIN_PCT.  If having multiple libraries, this file must be the concatenation of all the libaries DETREND file.");
  CommandArgumentBoolDefaultDoc("DETREND_AUTO", DETREND_AUTO, true, "If true, the sample names are standardized and must be all of the same library.  The correct detrending file will then be picked.");
  CommandArgumentStringDefaultDoc("DETREND_DIR", DETREND_DIR, "/root/data/libraries", "Full path of libraries in Repare's CRISPR screens.");
  CommandArgumentStringDefaultDoc("DETREND_COLUMNS", DETREND_COLUMNS, "4,5,6,7,8,9", "If empty, use all columns except first two (sgRNA & gene).  Otherwise, CSV of 0-based column indices to use.");
  CommandArgumentStringDefaultDoc("DETREND_BASE", DETREND_BASE, "0_0_0_1_0_0", "This is the base category against which other categories are to be detrended.");
  CommandArgumentFileOutDefaultDoc("DETREND_OUT", DETREND_OUT, "", "If specified, output file with detrending values.");
  CommandArgumentDoubleDefaultDoc("DETREND_MIN_PCT", DETREND_MIN_PCT, 0.3, "Minimum percentage of sgRNAs must be in a category for it to be used.");
  CommandArgumentFileInDefaultDoc("ESS", ESS, "/root/data/essential_genes/CEGv2.txt", "Essential gene list");
  CommandArgumentBoolDefaultDoc("ADD_RIBO", ADD_RIBO, false, "Add column whether something is an essential ribosomal gene.");
  CommandArgumentFileInDefaultDoc("ESS_RIBO", ESS_RIBO, "/root/data/essential_genes/RIBO.txt", "Essential ribosomal gene list");
  CommandArgumentFileInDefaultDoc("NONESS", NONESS, "/root/data/essential_genes/NEGv1.txt", "Nonessential gene list");
  CommandArgumentFileInDefaultDoc("INITQC", INITQC, "/root/src/Projects/Repare/ccaInitQc.R", "R script to do simple inital QC.");
  CommandArgumentStringDefaultDoc("INITQC_INPUT_PREFIX", INITQC_INPUT_PREFIX, "", "If specified, run initial QC on input counts.  Output head = input file + INITQC_INPUT_PREFIX, thus the first character should be a delimiter like a period.");
  CommandArgumentStringDefaultDoc("INITQC_PREFIX", INITQC_PREFIX, "", "If specified and CORRELATION_PREFIX is specified, run initial QC on normalized counts.");
  CommandArgumentSize_tDefaultDoc("NUM_THREADS", NUM_THREADS, thread::hardware_concurrency(), "Number of threads to use.");
  CommandArgumentBoolDefaultDoc("MIXED_LIBRARIES", MIXED_LIBRARIES, true, "If true, use med diff depletion = med(test depl) - med(ctrl depl) instead of med(test_i - ctrl_j) for all i, j, such that they are the same sgRNA.  Generally, no reason these days to make this false.");
  CommandArgumentBoolDefaultDoc("SGRNA_PAIRED_DIFFERENCES", SGRNA_PAIRED_DIFFERENCES, false, "If true, enable differences a statistics of differences within sgRNA depletions.  Incompatible with CHEMO, so either this is true or CHEMO is true or neither is true.  This method assumes the libraries are not mixed.");
  CommandArgumentBoolDefaultDoc("ADD_TOP_TWO", ADD_TOP_TWO, true, "If true, add statistics on top two differential sgRNAs for each gene.");
  CommandArgumentFileOutDefaultDoc("KILLING_PLOTS", KILLING_PLOTS, "", "If specified, then add measure of killing metric based on median depletion of essential genes and add plots of goodness of fits of killing are output to this file header.");
  CommandArgumentFileInDefaultDoc("EVAL_KILLING", EVAL_KILLING, "/root/src/Projects/Repare/ccaBetaDistrKill.R", "R script to fit the beta distribution to CCA essential genes output and provide a measure of killing.");
  CommandArgumentFileInDefaultDoc("FC_KILLING", FC_KILLING, "/root/src/Projects/Repare/ccaFcKill.R", "R script to plot distribution of median test depletions.");
  CommandArgumentSize_tDefaultDoc("TOP_FC_KILL_GENES", TOP_FC_KILL_GENES, 300, "Number of top genes to consider in plotting foldchange kill plots.");
  CommandArgumentDoubleDefaultDoc("NORMALIZED_COUNT", NORMALIZED_COUNT, 10000000.0, "Normalize the readcounts to this.");
  CommandArgumentSize_tDefaultDoc("MIN_T0", MIN_T0, 30, "If T0 is not atleast this, ignore the sgRNA.");
  CommandArgumentDoubleDefaultDoc("LARGE_NEG", LARGE_NEG, -10000.0, "Large negative number, an order of magnitude larger than any depletion.");
  CommandArgumentBoolDefaultDoc("ADD_RESISTORS", ADD_RESISTORS, true, "If true, add resistors score columns.   EVAL_BETA and CUTOFFS still done one depletors.");
  CommandArgumentFileInDefaultDoc("TSG_FILE", TSG_FILE, "/root/data/TSG/TSG.tsv", "File with list of tumor suppressor genes.");
  CommandArgumentFileInDefaultDoc("ENSEMBL_PARALOGS", ENSEMBL_PARALOGS, "/root/data/Ensembl/human_paralogs_list.tsv", "A TSV file with header and 2 columns.  1. gene, 2. CSV of paralog genes.");
  CommandArgumentBoolDefaultDoc("ADD_PARALOGS", ADD_PARALOGS, true, "If true, add paralog annotations to the output.   Important since CRISPR may target both the current gene and the paralogs.");
  CommandArgumentFileInDefaultDoc("EVAL_BETA", EVAL_BETA, "/root/src/Projects/Repare/ccaBetaDistrEval.R", "R script to fit the beta distribution to CCA output and provide cutoffs.");
  CommandArgumentFileInDefaultDoc("EVAL_BETA_Z", EVAL_BETA_Z, "/root/src/Projects/Repare/ccaBetaDistrEvalZscores.R", "R script to fit the beta distribution to CCA non-parametric Z-score output and provide cutoffs.");
  CommandArgumentSize_tDefaultDoc("BETA_TOP_N", BETA_TOP_N, 3000, "Number of top genes to use in determining the beta distribution");
  CommandArgumentDoubleDefaultDoc("CUTOFF_P", CUTOFF_P, 0.05, "Cutoff p-value.");
  CommandArgumentFileOutDefaultDoc("CUTOFF", CUTOFF, "", "If specified, output beta distribution parameters and maximum rank and cutoff.");
  CommandArgumentFileOutDefaultDoc("CUTOFF_Z", CUTOFF_Z, "", "If specified, output beta distribution parameters and maximum rank and cutoff.");
  CommandArgumentFileOutDefaultDoc("STRATAS", STRATAS, "", "If specified and CUTOFF is specified, then output Jenks natural cuts for top CUTOFF + STRATA_ADD genes.");
  CommandArgumentFileOutDefaultDoc("STRATAS_Z", STRATAS_Z, "", "If specified and CUTOFF is specified, then output Jenks natural cuts for top CUTOFF + STRATA_ADD genes.");
  CommandArgumentSize_tDefaultDoc("NUM_STRATA", NUM_STRATA, 4, "Group the scores into this many classes.");
  CommandArgumentSize_tDefaultDoc("STRATA_ADD", STRATA_ADD, 20, "Add this many genes to the CUTOFF to get the Jenks cuts.");
  CommandArgumentStringDefaultDoc("PRECISION_RECALL", PRECISION_RECALL, "", "Output PR plots to this prefix.");
  CommandArgumentFileInDefaultDoc("POSITIVE_CONTROLS", POSITIVE_CONTROLS, "", "If specificed, then this is a file of genes expected to come up in the CCA run.  There must be a header line.  The genes are assumed to be in the first column.");
  CommandArgumentBoolDefaultDoc("POSITIVE_AUTO", POSITIVE_AUTO, true, "If true, scan sample names for a single type of modification.  If only a single one is found, then use.");
  CommandArgumentStringDefaultDoc("POSITIVE_DIR", POSITIVE_DIR, "/root/data/positive_controls", "Top-level directory of the positive controls.");
  CommandArgumentFileInDefaultDoc("EVAL_PR", EVAL_PR, "/root/src/Projects/Repare/ccaPrPlots.R", "R script to plot PR curves.");
  CommandArgumentFileOutDefaultDoc("OUT", OUT, "", "If specified, assume all other files paths specified.");
  CommandArgumentStringDefaultDoc("OUT_MAX", OUT_MAX, "", "If specified, use this as head path to output files.   Output as much as possible.  This overrides OUT.");
  CommandArgumentBoolDefaultDoc("RUN_GATHER_QC", RUN_GATHER_QC, false, "If true and CORRELATION_PREFIX is specified, output QC to GATHER_QC_HTML and GATHER_QC_TSV.");
  CommandArgumentFileOutDefaultDoc("GATHER_QC_HTML", GATHER_QC_HTML, "", "Must specify when RUN_GATHER_QC is true.");
  CommandArgumentFileOutDefaultDoc("GATHER_QC_TSV", GATHER_QC_TSV, "", "Must specify when RUN_GATHER_QC is true.");
  CommandArgumentBoolDefaultDoc("ADD_UNIPROT", ADD_UNIPROT, true, "If true, add Uniprot and GO annotations to each gene.   Warning!  This will greatly increase the size of the CCA output table.");
  CommandArgumentFileInDefaultDoc("UNIPROT", UNIPROT, "/root/data/Uniprot/uniprot-human-gene-function-go.txt", "Uniprot file.");
  EndCommandArguments


  if (SGRNA_PAIRED_DIFFERENCES) {
    if (CHEMO) {
      cerr << "Error!  CHEMO and SGRNA_PAIRED_DIFFERENCES are incompatible.  Use one or the other.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    if (PARAMS != 3 && PARAMS != 4) {
      cerr << "PARAMS set to 3 since SGRNA_PAIRED_DIFFERENCES is enabled." 
           << endl;
    }
    PARAMS = 3;
  } else if (PARAMS == 3 || PARAMS == 4) {
     cerr << "Error!  PARAMS=3 and PARAMS=4 are only for SGRNA_PAIRED_DIFFERENCES mode.  Exiting." << endl;
     exit(EXIT_FAILURE);
  }

  if (OUT == "" && OUT_MAX == "") {
    cerr << "Error!  OUT or OUT_MAX must be specified.  Exiting." << endl;
    exit(EXIT_FAILURE);
  }

  if (AT_HOME) {
    if (OUT != "") OUT = "/home/" + OUT;
    if (OUT_MAX != "") OUT_MAX = "/home/" + OUT_MAX;
  }

  if (OUT_MAX != "") {
    INITQC_INPUT_PREFIX = ".init_input_qc";
    CORRELATION_PREFIX = OUT_MAX + ".correl";
    if (BARPLOT_NEG_LOG2) {
      BARPLOT_HEAD = OUT_MAX + ".neglog2foldchange_barplot";
    } else {
      BARPLOT_HEAD = OUT_MAX + ".foldchange_barplot";
    }
    INITQC_PREFIX = OUT_MAX + ".init_qc";
    DETREND_AUTO = false;
    // DETREND_OUT = OUT_MAX + ".detrend.tsv";
    KILLING_PLOTS = OUT_MAX + ".killing";
    CUTOFF = OUT_MAX + ".cutoff.tsv";
    CUTOFF_Z = OUT_MAX + ".cutoff_z.tsv";
    STRATAS = OUT_MAX + ".stratas";
    STRATAS_Z = OUT_MAX + ".stratas_z";
    // PRECISION_RECALL = OUT_MAX + ".precision_recall";
    OUT = OUT_MAX + ".cca.tsv";
    RESISTOR_OUT = OUT_MAX + ".resistor_cca.tsv";
    // if (fileSize(DETREND_OUT) != static_cast<ifstream::pos_type>(-1)) {
      // cerr << "Error!  " << DETREND_OUT << " exists.  Exiting." << endl;
      // exit(EXIT_FAILURE);
    // }
    if (fileSize(CUTOFF) != static_cast<ifstream::pos_type>(-1)) {
      cerr << "Error!  " << CUTOFF << " exists.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    if (fileSize(CUTOFF_Z) != static_cast<ifstream::pos_type>(-1)) {
      cerr << "Error!  " << CUTOFF_Z << " exists.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    if (fileSize(OUT) != static_cast<ifstream::pos_type>(-1)) {
      cerr << "Error!  " << OUT << " exists.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    RUN_GATHER_QC = true;
    GATHER_QC_HTML = OUT_MAX + ".qc.html";
    GATHER_QC_TSV = OUT_MAX + ".qc.tsv";
  }

  if (RESISTOR_SCREEN) {
    DETREND_AUTO = false;
    ADD_RESISTORS = true;
    STOP_AFTER_CORRELATION = true;
    if (RESISTOR_OUT == "") {
      cerr << "Error!  RESISTOR_OUT must be specified when RESISTOR_SCREEN is true.  Exiting."  << endl;
      exit(EXIT_FAILURE);
    }
  }

  if (RESISTOR_SCREEN && CORRELATION_PREFIX == "") {
    cerr << "Error!  A CORRELATION_PREFIX must be specified when doing an analyis on a resistor screen.  Exiting." << endl;
    exit(EXIT_FAILURE);
  }

  set<string> ess, ess_ribo, noness, tsg, pr_pos;
  // false below = there is a header line.
  readTSVSet<string, 0>(ess, ESS, false);
  readTSVSet<string, 0>(ess_ribo, ESS_RIBO, false);
  readTSVSet<string, 0>(noness, NONESS, false);
  if (ADD_RESISTORS) {
    readTSVSet<string, 0>(tsg, TSG_FILE, false);
  }
  
  unordered_map<string, string> sample2T0;
  map<string, string> test2ctrl;
  unordered_set<string> T0s;
  unordered_map<string, bool> all_samples, isControl;
  vector<string> repmaps;
  splitString(REPMAP, ",", repmaps, true);
  for (size_t k = 0; k < repmaps.size(); ++k) {
    if (AT_HOME) repmaps[k] = "/home/" + repmaps[k];
    if (fileSize(repmaps[k]) == static_cast<ifstream::pos_type>(-1)) {
      cerr << "Error!  REPMAP argument " << k << " = " << repmaps[k] << " does not exist.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
  }


  vector<string> count_files;
  splitString(COUNTS, ",", count_files, true);
  for (size_t k = 0; k < count_files.size(); ++k) {
    if (AT_HOME) count_files[k] = "/home/" + count_files[k];
    if (fileSize(count_files[k]) == static_cast<ifstream::pos_type>(-1)) {
      cerr << "Error!  COUNTS argument " << k << " = " << count_files[k] << " does not exist.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
  }
  unordered_set<string> all_count_samples;
  for (size_t k = 0; k < count_files.size(); ++k) {
    ReadOnlyTSVFile ifpTempCount(count_files[k], true, true);
    ifpTempCount.readRow(0);
    for (size_t k2 = 2; k2 < ifpTempCount.tokens.size(); ++k2) {
      all_count_samples.insert(ifpTempCount.tokens[k2]);
    }
    ifpTempCount.close();
  }
  for (size_t k = 0; k < repmaps.size(); ++k) {
    bool found_test = false, found_ctrl = false, found_t0 = false;
    ReadOnlyTSVFile ifpR(repmaps[k], true, true); 
    for (size_t i = 1; !ifpR.fail(); ++i) {
      ifpR.readRow(i);
      if (ifpR.tokens.size() < 3) {
        cerr << "Error!  " << repmaps[k] << " must have at least 3 columns.  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
      string& sample_name = ifpR.tokens[0];
      if (all_count_samples.find(sample_name) == all_count_samples.end()) {
        cerr << "Error!  "  << repmaps[k] << " on line " << i+1 << " does not have a corresponding sample in the COUNTS file(s).  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
      if (all_count_samples.find(ifpR.tokens[2]) == all_count_samples.end()) {
        cerr << "Error!  "  << repmaps[k] << " on line " << i+1 << " does not have a corresponding T0 sample in the COUNTS file(s).  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
      all_samples[sample_name] = false;
      if (ifpR.tokens[2] == sample_name) {
        T0s.insert(sample_name);
        found_t0 = true;
      } else {
        sample2T0[sample_name] = ifpR.tokens[2];
        if (ifpR.tokens[1] == "CTRL") {
          isControl[sample_name] = true;
          found_ctrl = true;
        } else {
          // sample is TEST
          isControl[sample_name] = false;
          found_test = true;
          if (CHEMO && CHEMO_PAIRED) {
            if (ifpR.tokens.size() < 4) {
              cerr << "Error!  In CHEMO_PAIRED mode, " << repmaps[k] << " must have at least 4 columns with the fourth column being the control pair for test.  Exiting." << endl;
              exit(EXIT_FAILURE);
            }
            if (all_count_samples.find(ifpR.tokens[3]) == all_count_samples.end()) {
              cerr << "Error!  "  << repmaps[k] << " on line " << i+1 << " does not have a corresponding paired control sample in the COUNTS file(s).  Exiting." << endl;
              exit(EXIT_FAILURE);
            }
            test2ctrl[sample_name] = ifpR.tokens[3];
          }
        }
      }
    }
    ifpR.close(); 
    if (!found_t0) {
      cerr << "Error!  No T0 found in " << repmaps[k] << ".  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    if (!found_ctrl && !RESISTOR_SCREEN) {
      cerr << "Error!  No CTRL found in " << repmaps[k] << ".  Exiting." << endl;
      exit(EXIT_FAILURE);
    } else if (found_ctrl && RESISTOR_SCREEN) {
      cerr << "Error!  CTRL found in " << repmaps[k] << " in a resistor screen.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    if (!found_test) {
      cerr << "Error!  No TEST (not CTRL) found in " << repmaps[k] << ".  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
  }
  if (CHEMO && CHEMO_PAIRED) {
    for (auto i : test2ctrl) {
      if (all_samples.find(i.second) == all_samples.end()) {
        cerr << "Error!  The pair " << i.second << " is not a known sample.  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
      if (!isControl[i.second]) {
        cerr << "Error!  The pair " << i.second << " is not a control sample.  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  bool has_pr_pos = false;
  if ((POSITIVE_CONTROLS != "") || POSITIVE_AUTO) {
    if (POSITIVE_CONTROLS == "") {
      // POSITIVE_AUTO is true
      unordered_set<string> modifiers;
      string the_mod = "", test_mod = "";
      for (auto i : all_samples) {
        vector<string> sample_tokens;
        splitString(i.first, "_", sample_tokens);
        if (sample_tokens.size() != 7) {
          cerr << "Error!  The input sample " << i.first << " does not conform to the standardized sample format = {screen id}_{library}_{cell line}_{lesion}_{modifier}_{clone}_T{time}{replicate}.  Exiting." << endl;
          exit(EXIT_FAILURE);
        }
        the_mod = sample_tokens[4];
        if (the_mod == "NT.2017.04" || the_mod == "NT.2017.07") {
          the_mod = "WT";   // Special case of WT sample-type naming.
        } else if (the_mod.find("WT.") != string::npos) {
          the_mod = afterString(the_mod, "WT.");
        } else if (the_mod.find(".WT") != string::npos) {
          the_mod = beforeString(the_mod, ".WT");
        } else if (the_mod.find("wt.") != string::npos) {
          the_mod = afterString(the_mod, "wt.");
        } else if (the_mod.find(".wt") != string::npos) {
          the_mod = beforeString(the_mod, ".wt");
        } else if (the_mod.find("dmso.") != string::npos) {
          the_mod = afterString(the_mod, "dmso.");
        } else if (the_mod.find(".dmso") != string::npos) {
          the_mod = beforeString(the_mod, ".dmso");
        } else if (the_mod.find("DMSO.") != string::npos) {
          the_mod = afterString(the_mod, "DMSO.");
        } else if (the_mod.find(".DMSO") != string::npos) {
          the_mod = beforeString(the_mod, ".DMSO");
        } else if (the_mod.find("nt.") != string::npos) {
          the_mod = afterString(the_mod, "nt.");
        } else if (the_mod.find(".nt") != string::npos) {
          the_mod = beforeString(the_mod, ".nt");
        } else if (the_mod.find("NT.") != string::npos) {
          the_mod = afterString(the_mod, "NT.");
        } else if (the_mod.find(".NT") != string::npos) {
          the_mod = beforeString(the_mod, ".NT");
        }
        if (the_mod != "WT" && the_mod != "control" && the_mod != "nt" &&
            the_mod != "nt" && the_mod != "NT" &&
            the_mod != "dmso" && the_mod != "DMSO" &&
            the_mod != "DMSO" && the_mod != "dmso") {
          test_mod = the_mod;
          modifiers.insert(the_mod);
        }
      }
      if (modifiers.size() != 1) {
        cerr << "Warning!  For positive controls, cannot have more than one modification.  One may manually specify POSITIVE_CONTROLS. Continuing without positive controls." << endl;
        POSITIVE_CONTROLS = "";  
      } else { 
        POSITIVE_CONTROLS = POSITIVE_DIR + "/" + test_mod + "/positive_controls.txt";
        if (fileSize(POSITIVE_CONTROLS) != static_cast<ifstream::pos_type>(-1)) {
          // cerr << "Using positive controls " << POSITIVE_CONTROLS << endl;
          ;  // NOP
        } else {
          POSITIVE_CONTROLS = "";
        }
      }
    }
    if (POSITIVE_CONTROLS != "") {
      has_pr_pos = true;
      readTSVSet<string, 0>(pr_pos, POSITIVE_CONTROLS, false);
    }
  }

  unordered_map<string, string> ensembl_paralogs;
  if (ADD_PARALOGS) {
    ReadOnlyTSVFile ifpEP(ENSEMBL_PARALOGS, true, true);
    for (size_t i = 1; !ifpEP.fail(); ++i) {
      ifpEP.readRow(i);
      ensembl_paralogs[ifpEP.tokens[0]] = ifpEP.tokens[1];
    }
    ifpEP.close();
  }

  bool has_detrend = (DETREND != "") || DETREND_AUTO ? true : false;
  unordered_map<string, string> sgrna_category;
  unordered_map< string, set<string> > gene_to_sgrnas;
  unordered_map<string, size_t> num_in_category;
  unordered_map<string, double> pct_in_category;
  unordered_map<string, double> sgrna_detrend_val_mut, sgrna_detrend_val_wt;
  if (has_detrend) {
    if (DETREND == "") {
      // DETREND_AUTO is true.
      unordered_set<string> libraries;
      string the_lib = "";
      for (auto i : all_samples) {
        vector<string> sample_tokens;
        splitString(i.first, "_", sample_tokens);
        if (sample_tokens.size() != 7) {
          cerr << "Error!  The input sample does not conform to the standardized sample format = {screen id}_{library}_{cell line}_{lesion}_{modifier}_{clone}_T{time}{replicate}.  Exiting." << endl;
          exit(EXIT_FAILURE);
        }
        the_lib = sample_tokens[1];
        if (the_lib == "RKOv1.RFP" || the_lib == "RKOv1.Cas9") {
          the_lib = "RKOv1";
        } else if (the_lib == "RFKOv1.RFP" || the_lib == "RFKOv1.Cas9") {
          the_lib = "RFKOv1";
        }
        libraries.insert(the_lib);
      }
      if (libraries.size() != 1) {
        cerr << "Warning!  For detrending, cannot have more than one library.  Continuing without detrending." << endl;
        DETREND = "";
        has_detrend = false;
      } else {
        DETREND = DETREND_DIR + "/" + the_lib + "/" + the_lib + ".sgrnaSummary.tsv";
        if (fileSize(DETREND) != static_cast<ifstream::pos_type>(-1)) {
          // cerr << "Using detrend file " << DETREND << endl;
          ;  // NOP
        } else {
          cerr << "Warning!  There is so such detrend file, " << DETREND << ", for the library " << the_lib << ".  Continuing without detrending." << endl;
          DETREND = "";
          has_detrend = false;
        }
      }
    }
    if (DETREND != "") {
      ReadOnlyTSVFile ifpD(DETREND, true, true);
      ifpD.readRow(0);
      vector<size_t> cols;
      if (DETREND_COLUMNS == "") {
        for (size_t i = 2; i < ifpD.tokens.size(); ++i) {
          cols.push_back(i);
        }
      } else {
        vector<string> tokens;
        splitString(DETREND_COLUMNS, ",", tokens);
        for (size_t i = 0; i < tokens.size(); ++i) {
          cols.push_back(convertFromString<size_t>(tokens[i]));
        }
      }

      size_t tot_sgrnas = 0;
      for (size_t i = 1; !ifpD.fail(); ++i) {
        ifpD.readRow(i);
        ++tot_sgrnas;
        string curr_sgrna = ifpD.tokens[0];
        gene_to_sgrnas[ifpD.tokens[1]].insert(curr_sgrna);
        string curr_cat = ifpD.tokens[cols[0]];
        for (size_t j = 1; j < cols.size(); ++j) {
          curr_cat += "_" + ifpD.tokens[cols[j]];
        }
        sgrna_category[curr_sgrna] = curr_cat; 
        ++num_in_category[curr_cat];
      }
      ifpD.close();

      if (num_in_category.find(DETREND_BASE) == num_in_category.end()) {
        cerr << "Error!  " << DETREND_BASE << " not found as a category in " << DETREND << ".  Exiting." << endl;
        exit(EXIT_FAILURE);
      }

      for (auto i : num_in_category) {
        pct_in_category[i.first] = 100.0*static_cast<double>(i.second)/
                                         static_cast<double>(tot_sgrnas);
      }
    }
  }


  //             gene                   sgRNA                 sample  count 
  unordered_map< string, unordered_map< string, unordered_map<string, float> > > count, t0, detrend_count;
  for (size_t k = 0; k < count_files.size(); ++k) {
    read_counts(count_files[k], all_samples, ALLOW_SUPERSET,
                sample2T0, T0s, t0, count, 
                MIN_T0, NORMALIZED_COUNT, MAX_MEDIAN_DIFF);
    if (INITQC_INPUT_PREFIX != "") {
      system("Rscript --vanilla \"" + INITQC + "\" \"" + 
             ESS + "\" \"" + NONESS + "\" \"" + count_files[k] + 
             "\" \"" + count_files[k] + INITQC_INPUT_PREFIX + "\"");
    }    
  }

  for (auto i : all_samples) {
    if (!i.second) {
      cerr << "Error!  Sample " << i.first << " never appears in any of the input COUNTS file(s).  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
  }

  if (MEDIAN_NORMALIZATION) {
    median_normalize_counts(t0, count);
  }

  if (MAX_CONTROL_RATIO > 0.0) {
    cap_control_ratio(t0, count, sample2T0, T0s, isControl, MAX_CONTROL_RATIO);
  }

  if (CORRELATION_PREFIX != "") {
    counts_corr(t0, count, CORRELATION_R, CORRELATION_PREFIX, MIN_CORRELATION_OUTPUT, CHEMO, CHEMO_PAIRED, test2ctrl);
    if (INITQC_PREFIX != "") {
      system("Rscript --vanilla \"" + INITQC + "\" \"" + 
             ESS + "\" \"" + NONESS + "\" \"" + CORRELATION_PREFIX +
             ".normed_counts.tsv\" \"" + INITQC_PREFIX + "\" norm");
    }
  }


  if (RESISTOR_SCREEN) {
    // Do primary resistor screen analysis here.

    Table<double, string> out_lines;
    for (auto gene : count) {
      string curr_gene = gene.first;
      string essential = "";
      if (ess.find(curr_gene) != ess.end()) {
        essential = "Essential";
      } else if (noness.find(curr_gene) != noness.end()) {
        essential = "Nonessential";
      }
      string pr_pos_string = "";
      if (has_pr_pos && pr_pos.find(curr_gene) != pr_pos.end()) {
        pr_pos_string = "Y";
      }
      string tsg_string = "";
      if (tsg.find(curr_gene) != tsg.end()) {
        tsg_string = "Y";
      }
       
      vector<float> depletions, curr_t0s;
      size_t num_sgrnas = 0;
      for (auto sgRNA : gene.second) {
        string curr_sgRNA = sgRNA.first;
        ++num_sgrnas;
        for (auto sample : sgRNA.second) {
          string curr_sample = sample.first;
          float curr_count =
              static_cast<float>(count[curr_gene][curr_sgRNA][curr_sample]);
          float curr_t0 =
              static_cast<float>(t0[curr_gene][curr_sgRNA][curr_sample]);
          depletions.push_back(diff_func(curr_count, curr_t0));
          curr_t0s.push_back(curr_t0);
        }
      } 

      sort(depletions.begin(), depletions.end());
      double q1_depl = 0.0, med_depl = 0.0, q3_depl = 0.0,
             q1_t0 = 0.0, med_t0 = 0.0, q3_t0 = 0.0; 
      tableQ123(depletions, med_depl, q1_depl, q3_depl);
      double iqr = q3_depl - q1_depl;
      sort(curr_t0s.begin(), curr_t0s.end());
      tableQ123(curr_t0s, med_t0, q1_t0, q3_t0);
      double iqr_t0 = q3_t0 - q1_t0;

      string line = curr_gene + "\t" + 
                    convertToString<size_t>(num_sgrnas) + "\t" + 
                    essential + "\t";
      if (has_pr_pos) {
        line += pr_pos_string + "\t";
      }
      line += tsg_string + "\t" + 
              convertToString<double>(q1_depl) + "\t" + 
              convertToString<double>(med_depl) + "\t" +
              convertToString<double>(q3_depl) + "\t" + 
              convertToString<double>(iqr) + "\t" +
              convertToString<double>(iqr/(q3_depl + q1_depl)) + "\t" +
              convertToString<double>(q1_t0) + "\t" + 
              convertToString<double>(med_t0) + "\t" +
              convertToString<double>(q3_t0) + "\t" + 
              convertToString<double>(iqr_t0) + "\t" +
              convertToString<double>(iqr_t0/(q3_t0 + q1_t0));
      out_lines.push_back(med_depl, line);
    }

    sortTable<0>(out_lines);

    Ofstream ofpR(RESISTOR_OUT);
    ofpR << "Gene\tNumber of sgRNAs per gene\tCore Ess Or Noness";
    if (has_pr_pos) {
      ofpR << "\tPositive Control Gene";
    }
    ofpR << "\tTumor Suppressor Gene\tQ1 depletion\tMedian depletion\tQ3 depletion\tIQR\tIQR/(Q3 + Q1)\tQ1 T0\tMedian T0\tQ3 T0\tIQR T0\tIQR/(Q3 + Q1) of T0" << endl;
    double curr_depl = 0.0;
    string curr_line = "";
    for (size_t i = 0; i < out_lines.size(); ++i) {
      out_lines.getRow(i, curr_depl, curr_line);
      ofpR << curr_line << endl;
    }
    ofpR.close();


    if (ADD_UNIPROT) {
      ReadOnlyTSVFile ifpU(UNIPROT, true, true);
      ifpU.readRow(0);
      bool gene_found = false;
      size_t gene_col = 0;
      for (size_t i = 0; !gene_found && i < ifpU.tokens.size(); ++i) {
        if (ifpU.tokens[i] == "Gene names") {
          gene_col = i;
          gene_found = true;
        }
      }
      if (!gene_found) {
        cerr << "Error!  Gene names column not found in " << UNIPROT << ".  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
      vector<string> tokens, blanks;
      for (size_t j = gene_col + 1; j < ifpU.tokens.size(); ++j) {
        tokens.push_back(ifpU.tokens[j]);
        blanks.push_back("");
      }
      string header = joinStrings(tokens);
      string blank_uni = joinStrings(blanks);
      unordered_map<string, string> gene_uniprot;
      for (size_t i = 1; !ifpU.fail(); ++i) {
        ifpU.readRow(i);
        tokens.clear();
        for (size_t j = gene_col + 1; j < ifpU.tokens.size(); ++j) {
          tokens.push_back(ifpU.tokens[j]);
        }
        string curr_uniprot = joinStrings(tokens);
        tokens.clear();
        splitStringOnSpaces(ifpU.tokens[gene_col], tokens);
        for (size_t j = 0; j < tokens.size(); ++j) {
          if ((j == 0) || 
              (gene_uniprot.find(trimSpacesString(tokens[j])) == 
               gene_uniprot.end())) {
            gene_uniprot[trimSpacesString(tokens[j])] = curr_uniprot;
          }
        }
      }
      ifpU.close();
  
      system("mv \"" + RESISTOR_OUT + "\" \"" + RESISTOR_OUT + ".temp_before_uniprot\"");
      system("sync");
      ReadOnlyTSVFile ifpUT(RESISTOR_OUT + ".temp_before_uniprot", true, true);
      ifpUT.readRow(0);
      if (ifpUT.tokens[0] != "Gene") {
        cerr << "Error!  CCA output file first column is not Gene.  Exiting." << endl;
        exit(EXIT_FAILURE);
      }
      Ofstream ofpUT(RESISTOR_OUT);
      ofpUT << ifpUT.line << "\t" << header << endl;
      for (size_t i = 1; !ifpUT.fail(); ++i) {
        ifpUT.readRow(i); 
        ofpUT << ifpUT.line << "\t";
        if (gene_uniprot.find(trimSpacesString(ifpUT.tokens[0])) != gene_uniprot.end()) {
          ofpUT << gene_uniprot[trimSpacesString(ifpUT.tokens[0])] << endl;
        } else {
          ofpUT << blank_uni << endl;
        }
      }
      ofpUT.close();
      ifpUT.close();
      system("rm -f \"" + RESISTOR_OUT + ".temp_before_uniprot\"");
    }
  }


  if (STOP_AFTER_CORRELATION) {
    // Below is part of the Imagemagick package.
    // mogrify trim will remove all white borders from the *.png files
    string OUT_DIR = "./";
    if (OUT.find("/") != string::npos) {
      // false below finds the last instance of / in the OUT path.
      OUT_DIR = beforeString(OUT, "/", false);
    }
    system("mogrify -quiet -trim " + OUT_DIR + "/*.png");
    if (RUN_GATHER_QC && CORRELATION_PREFIX != "" &&
        GATHER_QC_HTML != "" && GATHER_QC_TSV != "") {
      string at_home = AT_HOME ? "AT_HOME=true" : "";
      system("ccaGatherQC COUNTS=\"" + COUNTS + "\" CORREL_NORMED_COUNTS=\"" + CORRELATION_PREFIX + ".normed_counts.tsv\" TSV=\"" + GATHER_QC_TSV + "\" HTML=\"" + GATHER_QC_HTML + "\" " + at_home);
    }
    exit(EXIT_SUCCESS);
  }

  //    score   gene   output line  rank
  Table<double, string, string, size_t> output_lines, detrend_output_lines;
  map<string, size_t> gene_to_detrend_outline;
  compute_output_lines(ess, ess_ribo, noness, pr_pos, tsg, isControl, test2ctrl, t0, count, output_lines, LARGE_NEG, PARAMS, ADD_RESISTORS, ADD_RIBO, has_pr_pos, ADD_COEFF_VAR, CHEMO, CHEMO_PAIRED, CHEMO_NO_DEPLETION, SINGLE_CELL_LINE, MIXED_LIBRARIES, SGRNA_PAIRED_DIFFERENCES, ADD_TOP_TWO, NUM_THREADS);


  if (has_detrend) {
    detrend_count = count;
    // We never detrend t0.  We detrend only what comes afterwards.
    detrend_sgrnas(sgrna_category, num_in_category, pct_in_category, DETREND_MIN_PCT, T0s, isControl, t0, detrend_count, sgrna_detrend_val_mut, sgrna_detrend_val_wt, DETREND_BASE, DETREND_OUT);
    compute_output_lines(ess, ess_ribo, noness, pr_pos, tsg, isControl, test2ctrl, t0, detrend_count, detrend_output_lines, LARGE_NEG, PARAMS, ADD_RESISTORS, ADD_RIBO, has_pr_pos, ADD_COEFF_VAR, CHEMO, CHEMO_PAIRED, CHEMO_NO_DEPLETION, SINGLE_CELL_LINE, MIXED_LIBRARIES, SGRNA_PAIRED_DIFFERENCES, ADD_TOP_TWO, NUM_THREADS);
    for (size_t i = 0; i < detrend_output_lines.size(); ++i) {
      gene_to_detrend_outline[get<1>(detrend_output_lines.theData)[i]] = i;
    }
  }


  Ofstream ofp(OUT);
  ofp << "Gene\tCore Ess Or Noness";
  if (ADD_RIBO) {
    ofp << "\tRibosomal Essential Gene";
  }
  if (has_pr_pos) {
    ofp << "\tPositive Control Gene";
  }
  if (ADD_RESISTORS) {
    ofp << "\tTumor Suppressor Gene";
  }
  ofp << "\tNumber of sgRNAs";
  if (ADD_RESISTORS) {
    ofp << "\tTest Q1\tCtrl Q1\tTest-Ctrl Q1";
  }
  ofp << "\tTest median depletion\tCtrl median depletion\tTest-Ctrl median\tTest Q3\tCtrl Q3\tTest-Ctrl Q3\tBoth p-value\tLeft p-value\tRight p-value\tP-value of if non-essential in test - smaller is better\tP-value of if essential in ctrl - smaller is better\tSum of ess_noness p-vals";
  if (ADD_TOP_TWO) {
    ofp << "\tMedian test depletion top diff sgRNA\tMedian test-ctrl depletion top diff sgRNA\tMedian test depletion in 2nd top diff sgRNA\tMedian test-ctrl depletion 2nd top diff sgRNA\tSum median test-ctrl depletion of top and 2nd top diff sgRNAs";
  }
  if (CHEMO && CHEMO_PAIRED) {
    ofp << "\tFraction of sgRNA pairs with positive difference in depletion";
  }
  if (ADD_COEFF_VAR) {
    // IQR = interquartile range = Q3 - Q1
    // IQR/(Q3 + Q1) = quartile coefficient of dispersion
    //               ~ similar to stdev.  
    //               More robust than using Q2 = Median, since median may be 0.
    //               Above *100 = quartile variation coefficient.
    ofp << "\tIQR/(Q3 + Q1) of test depletion\tIQR/(Q3 + Q1) of ctrl depletion\tIQR/(Q3 + Q1) of test-ctrl depletion\tMedian test count\tIQR/(Q3 + Q1) test count\tMedian test T0 count\tIQR/(Q3 + Q1) test T0 count\tMedian ctrl count\tIQR/(Q3 + Q1) ctrl count\tMedian ctrl T0 count\tIQR/(Q3 + Q1) ctrl T0 count";
  }
  if (ADD_RESISTORS) {
    ofp << "\tP-value of it essential in test - smaller is better for resistors\tResistor score";
  }
  if (!CHEMO) {
    ofp << "\tScore without proliferation penalty\tRank without proliferation penalty";
  }
  ofp << "\tScore\tRank";
  if (ADD_PARALOGS) {
    ofp << "\tEnsembl paralogs";
  }
  if (has_detrend) {
    ofp << "\tDetrend categories\tDetrended Test amounts\tDetrended Ctrl amounts";
    if (ADD_RESISTORS) {
      ofp << "\tDetrended resistor score";
    }
    ofp << "\tDetrended Score\tDetrended Rank\t100 * (Detrended - original score) / original score\tDetrended - original rank";
  }
  ofp << endl;
  for (size_t i = 0; i < output_lines.size(); ++i) {
    string curr_line = get<2>(output_lines.theData)[i];
    ofp << curr_line << "\t"
        << get<3>(output_lines.theData)[i];

    vector<string> tokens;
    splitString(curr_line, "\t", tokens);
    string curr_gene = tokens[0];

    if (ADD_PARALOGS) {
      string curr_ensembl_paralogs = "";
      if (ensembl_paralogs.find(curr_gene) != ensembl_paralogs.end()) {
        curr_ensembl_paralogs = ensembl_paralogs[curr_gene];
      }
      ofp << "\t" << curr_ensembl_paralogs;
    }

    if (has_detrend) {
      set<string> curr_cats;
      set<string> curr_sgrnas = gene_to_sgrnas[curr_gene];
      string detrend_amt_mut = "", detrend_amt_wt = "";
      set<double> amt_done_mut, amt_done_wt;
      for (auto j : curr_sgrnas) {
        if (sgrna_detrend_val_mut.find(j) != sgrna_detrend_val_mut.end()) {
          double curr_val = sgrna_detrend_val_mut[j];
          if (curr_val != 1.0) {
            if (amt_done_mut.find(curr_val) != amt_done_mut.end()) continue;
            amt_done_mut.insert(curr_val);
            if (sgrna_category.find(j) != sgrna_category.end()) curr_cats.insert(sgrna_category[j]);
            if (detrend_amt_mut == "") {
              detrend_amt_mut = convertToString<double>(curr_val);
            } else {
              detrend_amt_mut += "," + convertToString<double>(curr_val);
            }
          }
        }
        if (sgrna_detrend_val_wt.find(j) != sgrna_detrend_val_wt.end()) {
          double curr_val = sgrna_detrend_val_wt[j];
          if (curr_val != 1.0) {
            if (amt_done_wt.find(curr_val) != amt_done_wt.end()) continue;
            amt_done_wt.insert(curr_val);
            if (sgrna_category.find(j) != sgrna_category.end()) curr_cats.insert(sgrna_category[j]);
            if (detrend_amt_wt == "") {
              detrend_amt_wt = convertToString<double>(curr_val);
            } else {
              detrend_amt_wt += "," + convertToString<double>(curr_val);
            }
          }
        }
      }
      string detrend_cat = "";
      for (auto j : curr_cats) {
        if (j == DETREND_BASE) continue;
        if (detrend_cat == "") {
          detrend_cat = j;
        } else {
          detrend_cat += "," + j;
        }
      }

      string detrend_line = get<2>(detrend_output_lines.theData)[gene_to_detrend_outline[curr_gene]];
      vector<string> detrend_line_tokens;
      splitString(detrend_line, "\t", detrend_line_tokens);
      string detrend_score = detrend_line_tokens[detrend_line_tokens.size() - 1];
      string detrend_resistor_score = "";
      if (ADD_RESISTORS) {
        detrend_resistor_score = detrend_line_tokens[detrend_line_tokens.size() - 2];
      }
      string detrend_score_diff = "";
      if (isADouble(detrend_score)) {
        double detrend_score_d = convertFromString<double>(detrend_score);
        double org_score = get<0>(output_lines.theData)[i];
        double diff_pct = 100.0*(detrend_score_d - org_score)/org_score;
        double epsilon = static_cast<double>(2.0e-8);
        if (fabs(diff_pct) < epsilon) diff_pct = 0.0;
        detrend_score_diff = convertToString<double>(diff_pct);
      }

      size_t detrend_rank = get<3>(detrend_output_lines.theData)[gene_to_detrend_outline[curr_gene]];
      ofp << "\t" << detrend_cat << "\t" 
          << detrend_amt_mut << "\t" << detrend_amt_wt << "\t";
      if (ADD_RESISTORS) {
        ofp << detrend_resistor_score << "\t";
      } 
      ofp << detrend_score << "\t" << detrend_rank << "\t"
          << detrend_score_diff << "\t"
          << static_cast<long>(detrend_rank) -
             static_cast<long>(get<3>(output_lines.theData)[i]);
    }
    ofp << endl;
  }
  ofp.close();


  if (BARPLOT_HEAD != "" && CORRELATION_PREFIX != "") {
    if (REPMAP.find(",") != string::npos) {
      cerr << "Warning!  Barplots are not available when using more than one REPMAP.  Skipping." << endl;
    } else {
      string neglog2 = BARPLOT_NEG_LOG2 ? "NEG_LOG2=true EPSILON=" + convertToString<double>(BARPLOT_EPSILON) : "";
      string barrepmap = AT_HOME ? "/home/" + REPMAP : REPMAP;
      string addrank = BARPLOT_ADD_RANK ? "ADD_RANK=true" : "ADD_RANK=false";
      system("barplotGenesFromCcaFoldchange AT_HOME=false GENES=\"" + OUT + "\" OUT_HEAD=\"" + BARPLOT_HEAD + "\" FOLDCHANGE=\"" + CORRELATION_PREFIX + ".foldchange.all_entries.tsv\" REPMAP=\"" + barrepmap + "\" NUM_GENES=" + convertToString<size_t>(BARPLOT_TOP_N) + " MIN=" + convertToString<double>(BARPLOT_MIN) + " MAX=" + convertToString<double>(BARPLOT_MAX) + " BY=" + convertToString<double>(BARPLOT_BY) + " " + neglog2 + " " + addrank);
    }
  }


  // The logic is easier to process the output and sort and rank the
  // score without a proliferation penalty, than to add it earlier.
  if (!CHEMO) {
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_rank_score_wo_prolif\"");
    system("sync");
    //    wo_prolif  org_rank
    ReadOnlyTSVFile ifpWP(OUT + ".temp_before_rank_score_wo_prolif", true, true);
    ifpWP.readRow(0);
    vector<string> header = ifpWP.tokens;
    size_t rank_col = numeric_limits<size_t>::max(),
           wo_prolif_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < header.size(); ++i) {
      if (header[i] == "Rank") rank_col = i;
      else if (header[i] == "Score without proliferation penalty") wo_prolif_col = i;
    }
    Table<double, size_t, vector<string> > wp_data;
    for (size_t i = 1; !ifpWP.fail(); ++i) {
      ifpWP.readRow(i);
      wp_data.push_back(convertFromString<double>(ifpWP.tokens[wo_prolif_col]),
                        convertFromString<size_t>(ifpWP.tokens[rank_col]),
                        ifpWP.tokens);
    }
    ifpWP.close();

    // True = sort by largest score without proliferation penality first.
    sortTable<0>(wp_data, true);
    for (size_t i = 0; i < wp_data.size(); ++i) {
      double curr_score = 0.0;
      size_t curr_rank = 0;
      vector<string> tokens; 
      wp_data.getRow(i, curr_score, curr_rank, tokens);
      // Rank is right after the score column.
      tokens[wo_prolif_col + 1] = convertToString<size_t>(i+1);
      wp_data.setRow(i, curr_score, curr_rank, tokens);
    }
    // Sort by rank.
    sortTable<1>(wp_data);
    Ofstream ofpWP(OUT);
    ofpWP << joinStrings(header) << endl;
    for (size_t i = 0; i < wp_data.size(); ++i) {
      double curr_score = 0.0;
      size_t curr_rank = 0;
      vector<string> tokens; 
      wp_data.getRow(i, curr_score, curr_rank, tokens);
      ofpWP << joinStrings(tokens) << endl;
    }
    ofpWP.close();
    
    system("rm \"" + OUT + ".temp_before_rank_score_wo_prolif\"");
  }


  if (POSITIVE_CONTROLS != "") {
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_precision_recall\"");
    system("sync");
    ReadOnlyTSVFile ifpPR(OUT + ".temp_before_precision_recall", true, true);
    ifpPR.readRow(0);
    if (ifpPR.tokens[0] != "Gene") {
      cerr << "Error!  The first column is no longer the Gene column.  This cannot happen.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    Ofstream ofpPR(OUT);
    ofpPR << ifpPR.line << "\tCCA score recall\tCCA score precision" << endl;
    double totNumPos = static_cast<double>(pr_pos.size());
    double cumulative_tp = 0.0, cumulative_fp = 0.0, 
           precision = 1.0, recall = 0.0;
    for (size_t i = 1; !ifpPR.fail(); ++i) {
      ifpPR.readRow(i);
      string& curr_gene = ifpPR.tokens[0];
      if (pr_pos.find(curr_gene) != pr_pos.end()) {
        cumulative_tp += 1.0;
      } else if (noness.find(curr_gene) != noness.end()) {
        cumulative_fp += 1.0;
      }
      recall = cumulative_tp / totNumPos;
      if (cumulative_tp > 0 || cumulative_fp > 0) {
        precision = cumulative_tp/(cumulative_tp +  cumulative_fp);
      }
      ofpPR << ifpPR.line << "\t" << recall << "\t" << precision << endl;
    }
    ofpPR.close();
    ifpPR.close();
    system("rm \"" + OUT + ".temp_before_precision_recall\"");

    system("Rscript --vanilla \"" + EVAL_PR + "\" \"" + OUT + "\" \"" +
            PRECISION_RECALL + ".CCA_Score\" 0");
  }
 

  if (ADD_Z_SCORES && CHEMO && CHEMO_PAIRED && CORRELATION_PREFIX != "") {
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_paired_diff_zvals\"");
    system("sync");

#if 0
    // 1. compute Q1, Q2, Q3 of depletion.
    map< size_t, vector<double> > sample2vals;
    ReadOnlyTSVFile ifp4z(CORRELATION_PREFIX + ".depletion_test_minus_ctrl.tsv", true, true);
    ifp4z.readRow(0);
    vector<string> header = ifp4z.tokens;
    for (size_t i = 1; !ifp4z.fail(); ++i) {
      ifp4z.readRow(i);
      for (size_t j = 2; j < ifp4z.tokens.size(); ++j) {
        sample2vals[j].push_back(convertFromString<double>(ifp4z.tokens[j]));
      }
    }
    ifp4z.resetToBegin();

    map<size_t, double> sample2q1, sample2med, sample2q3;
    for (auto i : sample2vals) {
      sort(i.second.begin(), i.second.end());
      tableQ123(i.second, sample2med[i.first], sample2q1[i.first], sample2q3[i.first]);
    }

    // 2. compute sum Z-score and divide by the number of samples in each gene 
    //            gene
    unordered_map<string, vector<double> > depl_zscores;
    for (size_t i = 1; !ifp4z.fail(); ++i) {
      ifp4z.readRow(i);
      string& curr_gene = ifp4z.tokens[0];
      for (size_t j = 2; j < ifp4z.tokens.size(); ++j) {
        double curr_val = convertFromString<double>(ifp4z.tokens[j]);
        // Below is still 1e-32 or so.
        double iqr = 1000000.0*numeric_limits<double>::min();
        if (sample2q3[j] != sample2q1[j]) iqr = sample2q3[j] - sample2q1[j];
        double curr_z = (curr_val - sample2med[j])/iqr;
        depl_zscores[curr_gene].push_back(curr_z);
      } 
    }
    ifp4z.close();

    // 3. compute Q1, Q2, Q3 of logfoldchange.
    sample2vals.clear();
    ReadOnlyTSVFile ifp5z(CORRELATION_PREFIX + ".log2foldchange_test_minus_ctrl.tsv", true, true);
    ifp5z.readRow(0);
    header.clear();
    header = ifp5z.tokens;
    for (size_t i = 1; !ifp5z.fail(); ++i) {
      ifp5z.readRow(i);
      for (size_t j = 2; j < ifp5z.tokens.size(); ++j) {
        sample2vals[j].push_back(convertFromString<double>(ifp5z.tokens[j]));
      }
    }
    ifp5z.resetToBegin();

    sample2q1.clear();  sample2med.clear();  sample2q3.clear();
    map<size_t, double> sample2mean, sample2stdev;
    for (auto i : sample2vals) {
      sort(i.second.begin(), i.second.end());
      tableQ123(i.second, sample2med[i.first], sample2q1[i.first], sample2q3[i.first]);
      sample2mean[i.first] = vecMean(i.second);
      sample2stdev[i.first] = vecStdev(i.second, sample2mean[i.first]);
    }

    // 4. compute sum Z-score and divide by the number of samples in each gene 
    //            gene
    unordered_map<string, vector<double> > fc_zscores;
    for (size_t i = 1; !ifp5z.fail(); ++i) {
      ifp5z.readRow(i);
      string& curr_gene = ifp5z.tokens[0];
      for (size_t j = 2; j < ifp5z.tokens.size(); ++j) {
        double curr_val = convertFromString<double>(ifp5z.tokens[j]);
        // Below is still 1e-32 or so.
        // double iqr = 1000000.0*numeric_limits<double>::min();
        // if (sample2q3[j] != sample2q1[j]) iqr = sample2q3[j] - sample2q1[j];
        // double curr_z = (curr_val - sample2med[j])/iqr;
        // double curr_z = (curr_val - sample2mean[j])/sample2stdev[j];
        // fc_zscores[curr_gene].push_back(curr_z);
        fc_zscores[curr_gene].push_back(curr_val);
      } 
    }
    ifp5z.close();
#endif


    //            test sample      gene   sgRNA  ctl cnt  lg-ratio  deviation
    unordered_map< string, Table<string, string, float, float, float> > zdata;
    ReadOnlyTSVFile ifpz(CORRELATION_PREFIX + ".normed_counts.tsv", true, true);
    ifpz.readRow(0);
    vector<string> header = ifpz.tokens;
    map<size_t, size_t> test2ctrl_col;
    for (size_t j = 2; j < header.size(); ++j) {
      if (test2ctrl.find(header[j]) != test2ctrl.end()) {
        string curr_ctrl = test2ctrl[header[j]];
        bool found = false;
        for (size_t k = 2; !found && k < header.size(); ++k) {
          if (header[k] == curr_ctrl) {
            found = true;
            test2ctrl_col[j] = k;
          }
        }
      }
    }
    for (size_t i = 1; !ifpz.fail(); ++i) {
      ifpz.readRow(i);
      string& curr_sgrna = ifpz.tokens[0];
      string& curr_gene = ifpz.tokens[1];
      for (auto j : test2ctrl_col) {
        string& curr_test_sample = header[j.first];
        float curr_test_val = convertFromString<float>(ifpz.tokens[j.first]);
        float curr_ctrl_val = convertFromString<float>(ifpz.tokens[j.second]);
        // We can do this instead of dividing each by T0 because for paired
        // chemogenomic screens, T0 is always shared between the pairs.
        // float log_ratio = log2f(curr_test_val/curr_ctrl_val);
         float log_ratio = log2f((curr_test_val+Z_SCORES_MIN)/(curr_ctrl_val+Z_SCORES_MIN));
        zdata[curr_test_sample].push_back(curr_gene, curr_sgrna, curr_ctrl_val, log_ratio, 0.0);
      }
    }
    ifpz.close();

    //             gene                   sgRNA                 sample
    unordered_map< string, unordered_map< string, unordered_map<string, float> > > zscores;
    for (auto i : zdata) {
      string curr_sample = i.first;
      Table<string, string, float, float, float>& curr_table = i.second;
      // 1. Get median of log fold-changes.  Sort on fcs and get median.
      sortTable<3>(curr_table);
      float lg_fc_median = vecMedian(get<3>(curr_table.theData));
      // float lg_fc_mean = vecMean(get<3>(curr_table.theData));

      // Sort on decreasing control counts.
      size_t n = curr_table.size();
      sortTable<2>(curr_table, false);

      // 2. Compute window size
      if (Z_SCORES_HALF_WINDOW == 0) {
        double bestDist = numeric_limits<double>::max();
        Z_SCORES_HALF_WINDOW = 1;
        for (size_t j = Z_SCORES_WINDOW_MIN; j < n; ++j) {
          double currDA = tableColumnDAgostinoK2<3>(curr_table, 0, j);
          double currDist = fabs(2.0 - currDA);
          if (currDist < bestDist) {
            bestDist = currDist;
            Z_SCORES_HALF_WINDOW = j/2;
          }
        }
        cerr << "For " << curr_sample << ", best D'Agostino test distance away = " << bestDist << " with window size = " << 2*Z_SCORES_HALF_WINDOW << endl;
      } else {
        Z_SCORES_HALF_WINDOW = max(Z_SCORES_HALF_WINDOW, 
                                   Z_SCORES_WINDOW_MIN/2);
      }

      // float last_deviation = -1.0;
      // 3. Compute deviation on windows.
      for (size_t j = 0; j < n; ++j) {
        size_t start = 0, end = n-1;
        if (j < Z_SCORES_HALF_WINDOW) {
          // If too close to the beginning, compensate by having more above.
          end = Z_SCORES_HALF_WINDOW + (Z_SCORES_HALF_WINDOW - j); 
        } else if (j >= (n-Z_SCORES_HALF_WINDOW-1)) {
          // If too close to the end, compensate by having more below.
          start = j - Z_SCORES_HALF_WINDOW - (j + Z_SCORES_HALF_WINDOW - n);
        } else {
          start = j - Z_SCORES_HALF_WINDOW;
          end = j + Z_SCORES_HALF_WINDOW;
        }
        size_t m = end - start + 1; 
        vector<float> temp_vec(m, 0.0);
        curr_table.getCol<3>(temp_vec, start, end);
        sort(temp_vec.begin(), temp_vec.end());
        double temp_q1, temp_q2, temp_q3;
        tableQ123(temp_vec, temp_q2, temp_q1, temp_q3);
        float deviation = static_cast<float>(temp_q3 - temp_q1); 
        // float mean = vecMean(temp_vec);
        // float deviation = vecStdev(temp_vec, mean);
        // We have to put a lower-limit on the deviation.
        // if (deviation == 0.0) deviation = numeric_limits<float>::min();
        // if (deviation < 1.0e-5) deviation = 1.0e-5;
        // Below always forces the deviation in the denominator to go up.
        // Similar to Traver Hart's drugZ method.
        // if (deviation >= last_deviation) {
          // last_deviation = deviation;
        // } else {
          // deviation = last_deviation;
        // }
        curr_table.setElement<4>(j, deviation);
      }
      // Smooth the deviation
      Table<double> smoothed_deviation;
      movingAverageTableSmoother<4>(curr_table, Z_SCORES_HALF_WINDOW, smoothed_deviation);

      for (size_t j = 0; j < n; ++j) {
        string curr_gene = "", curr_sgrna = "";
        float dummy_count = 0.0, curr_lg_ratio = 0.0,
              dummy_deviation = 0.0;
        curr_table.getRow(j, curr_gene, curr_sgrna, dummy_count, curr_lg_ratio,
                             dummy_deviation);
        float numerator = curr_lg_ratio - lg_fc_median;
        // float numerator = curr_lg_ratio - lg_fc_mean;
        // We have to put a lower-limit on the deviation.
        if (smoothed_deviation[j] < 1.0e-6) smoothed_deviation[j] = 1.0e-6;
        float curr_zscore = numerator/smoothed_deviation[j];

        zscores[curr_gene][curr_sgrna][curr_sample] = curr_zscore;
      }       
    }


    // 5. Read in CCA output and add normalized Z-scores.
    ReadOnlyTSVFile ifpZ(OUT + ".temp_paired_diff_zvals", true, true);
    ifpZ.readRow(0);
    size_t gene_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpZ.tokens.size(); ++i) {
      if (ifpZ.tokens[i] == "Gene") gene_col = i;
    }
    Ofstream ofpZ(OUT + ".temp_paired_before_zscore_sort");
    // ofpZ << ifpZ.line << "\tNumber of Z-score observations\tDepletion normalized non-parameteric Z-score\tLog fold-change normalized parameteric Z-score" << endl;
    ofpZ << ifpZ.line << "\tNumber of Z-score observations\tNegative non-parameteric Z-score" << endl;
    for (size_t i = 1; !ifpZ.fail(); ++i) {
      ifpZ.readRow(i);
      string& curr_gene = ifpZ.tokens[gene_col];
#if 0
      vector<double>& curr_depl_zscores = depl_zscores[curr_gene];
      double num_z = static_cast<double>(curr_depl_zscores.size());
      double sum_depl_zscores = accumulate(curr_depl_zscores.begin(),
                                           curr_depl_zscores.end(), 0.0);
      vector<double>& curr_fc_zscores = fc_zscores[curr_gene];
      double sum_fc_zscores = accumulate(curr_fc_zscores.begin(),
                                         curr_fc_zscores.end(), 0.0);
      ofpZ << ifpZ.line << "\t"
           << num_z << "\t"
           << sum_depl_zscores/num_z << "\t"
           << sum_fc_zscores/sqrt(num_z)
           << endl;
#endif
      // Sum all sgRNA Z scores to get a gene-level Z score.
      float num_z = 0.0;
      float sum_z = 0.0;
      for (auto j : zscores[curr_gene]) {
        for (auto k : j.second) {
          num_z += 1.0;
          sum_z += k.second;
        }
      }
      if (num_z > 0.0) {
        float final_z = sum_z/sqrtf(num_z);
        ofpZ << joinStrings(ifpZ.tokens) << "\t" << num_z << "\t" << -1.0*final_z << endl;
      } else {
        // Put some junk value, LARGE_NEG, that clearly shows this was because
        // num_z is 0.   This can happen because the CTRL sgRNAs disappear due
        // to low T0 counts.
        ofpZ << joinStrings(ifpZ.tokens) << "\t" << num_z << "\t" << LARGE_NEG << endl;
      }
    }
    ofpZ.close();
    ifpZ.close();


    //    Score   Z-score                 Z-score rank
    Table< double, double, vector<string>, size_t > temp_paired;
    ReadOnlyTSVFile ifpZ2(OUT + ".temp_paired_before_zscore_sort", true, true);
    ifpZ2.readRow(0);
    string temp_header = ifpZ2.line;
    size_t score_col = numeric_limits<size_t>::max(),
           rank_col = numeric_limits<size_t>::max(),
           zscore_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpZ2.tokens.size(); ++i) {
      if (ifpZ2.tokens[i] == "Score") score_col = i;
      else if (ifpZ2.tokens[i] == "Rank") rank_col = i;
      else if (ifpZ2.tokens[i] == "Negative non-parameteric Z-score") zscore_col = i;
    }
    for (size_t i = 1; !ifpZ2.fail(); ++i) {
      ifpZ2.readRow(i); 
      double score = convertFromString<double>(ifpZ2.tokens[score_col]);
      double zscore = convertFromString<double>(ifpZ2.tokens[zscore_col]);
      temp_paired.push_back(score, zscore, ifpZ2.tokens, 0);
    }
    ifpZ2.close();

    // Sort table in reverse Z-score order.
    sortTable<1>(temp_paired, true);
    for (size_t i = 0; i < temp_paired.size(); ++i) {
      // Add Z-score rank.
      temp_paired.setElement<3>(i, i+1);
      if (Z_SCORES_L != 0.0) {
        // While we're at it, compute new score, if integrating Z-scores.
        double curr_score = temp_paired.getElement<0>(i);
        double curr_zscore = temp_paired.getElement<1>(i);
        // TO: I did try the idea of multiplying the curr_zscore by
        //     (1-both^G)*(1-right^H), but it didn't work as well.   Keep
        //     the Z-score separate from the probabilities.
        curr_score += Z_SCORES_L*curr_zscore;
        temp_paired.setElement<0>(i, curr_score);
      }
    }

    // Below sort in reverse order (biggest score first).
    sortTable<0>(temp_paired, true);

    Ofstream ofpZ2(OUT);
    ofpZ2 << temp_header << "\tNon-parametric Z-score rank" << endl;
    for (size_t i = 0; i < temp_paired.size(); ++i) {
      double curr_score = 0.0, curr_zscore = 0.0;
      vector<string> curr_line;
      size_t curr_zrank = 0;
      temp_paired.getRow(i, curr_score, curr_zscore, curr_line, curr_zrank);
      if (Z_SCORES_L != 0.0) {
        curr_line[score_col] = convertToString<double>(curr_score);
        curr_line[rank_col] = convertToString<size_t>(i+1);
      }
      ofpZ2 << joinStrings(curr_line) << "\t" << curr_zrank << endl;
    }
    ofpZ2.close();

    system("rm \"" + OUT + ".temp_paired_diff_zvals\" \"" + OUT + ".temp_paired_before_zscore_sort\"");

  } else if (ADD_Z_SCORES && (!CHEMO || !CHEMO_PAIRED) && CORRELATION_PREFIX != "") {

    system("mv \"" + OUT + "\" \"" + OUT + ".temp_diff_zvals\"");
    system("sync");

    // Since we have no paired genes, we'll do this at the gene level using
    // log(median(test_fcs)/median(ctrl_fcs))
    float epsilon = static_cast<float>(1.0/NORMALIZED_COUNT);

    //             gene    
    unordered_map< string, vector<float> > test_fcs, ctrl_fcs;
    unordered_map<string, float> num_fcs;
    ReadOnlyTSVFile ifpz(CORRELATION_PREFIX + ".foldchange.tsv", true, true);
    ifpz.readRow(0);
    vector<string> header = ifpz.tokens;
    unordered_set<size_t> test_cols, ctrl_cols;
    for (size_t j = 2; j < header.size(); ++j) {
      if (isControl.find(header[j]) != isControl.end()) {
        if (isControl[header[j]]) {
          ctrl_cols.insert(j);
        } else {
          test_cols.insert(j);
        }
      }
    }
    for (size_t i = 1; !ifpz.fail(); ++i) {
      ifpz.readRow(i);
      // string& curr_sgrna = ifpz.tokens[0];
      string& curr_gene = ifpz.tokens[1];
      for (auto j : test_cols) {
        test_fcs[curr_gene].push_back(convertFromString<float>(ifpz.tokens[j]));
      }
      for (auto j : ctrl_cols) {
        ctrl_fcs[curr_gene].push_back(convertFromString<float>(ifpz.tokens[j]));
      }
    }
    ifpz.close();
    //    gene  med ctrl lg-ratio  deviation
    Table<string, float, float, float> zdata;
    for (auto i : test_fcs) {
      string curr_gene = i.first;
      float med_test = vecMedian(i.second),
            med_ctrl = vecMedian(ctrl_fcs[curr_gene]);
      float lg_ratio = log2f((med_test+epsilon)/(med_ctrl+epsilon));
      // Median/average number of observations.
      float med_num_obvs = static_cast<float>(i.second.size() + ctrl_fcs[curr_gene].size())/2.0;
      num_fcs[curr_gene] = med_num_obvs;
      zdata.push_back(curr_gene, med_ctrl, lg_ratio, 0.0);
    }

    //             gene
    unordered_map<string, float> zscores;
    // 1. Get median of log fold-changes.  Sort on fcs and get median.
    sortTable<2>(zdata);
    float lg_fc_median = vecMedian(get<2>(zdata.theData));
    // Sort on decreasing control counts.
    size_t n = zdata.size();
    sortTable<1>(zdata, false);

    // 2. Compute window size
    if (Z_SCORES_HALF_WINDOW == 0) {
      double bestDist = numeric_limits<double>::max();
      Z_SCORES_HALF_WINDOW = 1;
      for (size_t j = Z_SCORES_WINDOW_MIN; j < n; ++j) {
        double currDA = tableColumnDAgostinoK2<2>(zdata, 0, j);
        double currDist = fabs(2.0 - currDA);
        if (currDist < bestDist) {
          bestDist = currDist;
          Z_SCORES_HALF_WINDOW = j/2;
        }
      }
      cerr << "Best D'Agostino test distance away = " << bestDist << " with window size = " << 2*Z_SCORES_HALF_WINDOW << endl;
    } else {
      Z_SCORES_HALF_WINDOW = max(Z_SCORES_HALF_WINDOW, 
                                 Z_SCORES_WINDOW_MIN/2);
    }
    // float last_deviation = -1.0;
    // 3. Compute deviation on windows.
    for (size_t j = 0; j < n; ++j) {
      size_t start = 0, end = n-1;
      if (j < Z_SCORES_HALF_WINDOW) {
        // If too close to the beginning, compensate by having more above.
        end = Z_SCORES_HALF_WINDOW + (Z_SCORES_HALF_WINDOW - j); 
      } else if (j >= (n-Z_SCORES_HALF_WINDOW-1)) {
        // If too close to the end, compensate by having more below.
        start = j - Z_SCORES_HALF_WINDOW - (j + Z_SCORES_HALF_WINDOW - n);
      } else {
        start = j - Z_SCORES_HALF_WINDOW;
        end = j + Z_SCORES_HALF_WINDOW;
      }
      size_t m = end - start + 1; 
      vector<float> temp_vec(m, 0.0);
      zdata.getCol<2>(temp_vec, start, end);
      sort(temp_vec.begin(), temp_vec.end());
      double temp_q1, temp_q2, temp_q3;
      tableQ123(temp_vec, temp_q2, temp_q1, temp_q3);
      float deviation = static_cast<float>(temp_q3 - temp_q1); 
      // float mean = vecMean(temp_vec);
      // float deviation = vecStdev(temp_vec, mean);
      // We have to put a lower-limit on the deviation.
      // if (deviation == 0.0) deviation = numeric_limits<float>::min();
      // if (deviation < 1.0e-5) deviation = 1.0e-5;
      // Below always forces the deviation in the denominator to go up.
      // Similar to Traver Hart's drugZ method.
      // if (deviation >= last_deviation) {
        // last_deviation = deviation;
      // } else {
        // deviation = last_deviation;
      // }
      zdata.setElement<3>(j, deviation);
    }

    // Smooth the deviation
    Table<double> smoothed_deviation;
    movingAverageTableSmoother<3>(zdata, Z_SCORES_HALF_WINDOW, smoothed_deviation);

    for (size_t j = 0; j < n; ++j) {
      string curr_gene = "";
      float dummy_count = 0.0, curr_lg_ratio = 0.0,
            dummy_deviation = 0.0;
      zdata.getRow(j, curr_gene, dummy_count, curr_lg_ratio, dummy_deviation);
      float numerator = curr_lg_ratio - lg_fc_median;
      // float numerator = curr_lg_ratio - lg_fc_mean;
      // We have to put a lower-limit on the deviation.
      // Below should be 6 orders of magnitude smaller than paired version
      // since we're doing fc instead of read counts
      if (smoothed_deviation[j] < 1.0e-12) smoothed_deviation[j] = 1.0e-12;
      zscores[curr_gene] = numerator/smoothed_deviation[j];
    }       

    // 5. Read in CCA output and add normalized Z-scores.
    ReadOnlyTSVFile ifpZ(OUT + ".temp_diff_zvals", true, true);
    ifpZ.readRow(0);
    size_t gene_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpZ.tokens.size(); ++i) {
      if (ifpZ.tokens[i] == "Gene") gene_col = i;
    }
    Ofstream ofpZ(OUT + ".temp_before_zscore_sort");
    ofpZ << ifpZ.line << "\tMedian number of Z-score observations of test and control\tNegative non-parameteric Z-score" << endl;
    for (size_t i = 1; !ifpZ.fail(); ++i) {
      ifpZ.readRow(i);
      string& curr_gene = ifpZ.tokens[gene_col];
      if (num_fcs[curr_gene] > 0.0) {
        float final_z = zscores[curr_gene]/sqrtf(num_fcs[curr_gene]);
        ofpZ << joinStrings(ifpZ.tokens) << "\t" << num_fcs[curr_gene] << "\t" << -1.0*final_z << endl;
      } else {
        // Put some junk value, LARGE_NEG, that clearly shows this was because
        // num_z is 0.   This can happen because the CTRL sgRNAs disappear due
        // to low T0 counts.
        ofpZ << joinStrings(ifpZ.tokens) << "\t" << num_fcs[curr_gene] << "\t" << LARGE_NEG << endl;
      }
    }
    ofpZ.close();
    ifpZ.close();


    //    Score   Z-score                 Z-score rank
    Table< double, double, vector<string>, size_t > temp_paired;
    ReadOnlyTSVFile ifpZ2(OUT + ".temp_before_zscore_sort", true, true);
    ifpZ2.readRow(0);
    string temp_header = ifpZ2.line;
    size_t score_col = numeric_limits<size_t>::max(),
           rank_col = numeric_limits<size_t>::max(),
           zscore_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpZ2.tokens.size(); ++i) {
      if (ifpZ2.tokens[i] == "Score") score_col = i;
      else if (ifpZ2.tokens[i] == "Rank") rank_col = i;
      else if (ifpZ2.tokens[i] == "Negative non-parameteric Z-score") zscore_col = i;
    }
    for (size_t i = 1; !ifpZ2.fail(); ++i) {
      ifpZ2.readRow(i); 
      double score = convertFromString<double>(ifpZ2.tokens[score_col]);
      double zscore = convertFromString<double>(ifpZ2.tokens[zscore_col]);
      temp_paired.push_back(score, zscore, ifpZ2.tokens, 0);
    }
    ifpZ2.close();

    // Sort table in reverse Z-score order.
    sortTable<1>(temp_paired, true);
    for (size_t i = 0; i < temp_paired.size(); ++i) {
      // Add Z-score rank.
      temp_paired.setElement<3>(i, i+1);
      if (Z_SCORES_L != 0.0) {
        // While we're at it, compute new score, if integrating Z-scores.
        double curr_score = temp_paired.getElement<0>(i);
        double curr_zscore = temp_paired.getElement<1>(i);
        // TO: I did try the idea of multiplying the curr_zscore by
        //     (1-both^G)*(1-right^H), but it didn't work as well.   Keep
        //     the Z-score separate from the probabilities.
        curr_score += Z_SCORES_L*curr_zscore;
        temp_paired.setElement<0>(i, curr_score);
      }
    }

    // Below sort in reverse order (biggest score first).
    sortTable<0>(temp_paired, true);

    Ofstream ofpZ2(OUT);
    ofpZ2 << temp_header << "\tNon-parametric Z-score rank" << endl;
    for (size_t i = 0; i < temp_paired.size(); ++i) {
      double curr_score = 0.0, curr_zscore = 0.0;
      vector<string> curr_line;
      size_t curr_zrank = 0;
      temp_paired.getRow(i, curr_score, curr_zscore, curr_line, curr_zrank);
      if (Z_SCORES_L != 0.0) {
        curr_line[score_col] = convertToString<double>(curr_score);
        curr_line[rank_col] = convertToString<size_t>(i+1);
      }
      ofpZ2 << joinStrings(curr_line) << "\t" << curr_zrank << endl;
    }
    ofpZ2.close();

    system("rm \"" + OUT + ".temp_diff_zvals\" \"" + OUT + ".temp_before_zscore_sort\"");
  }


  if (ADD_COEFF_VAR) {
    system("mv " + OUT + " " + OUT + ".temp_coeff_var");
    system("sync");
    vector<float> test_cvs, ctrl_cvs, diff_cvs;
    ReadOnlyTSVFile ifpV(OUT + ".temp_coeff_var", true, true);
    ifpV.readRow(0);
    size_t tcol = numeric_limits<size_t>::max(), 
           ccol = numeric_limits<size_t>::max(),
           dcol = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpV.tokens.size(); ++i) {
      if (ifpV.tokens[i] == "IQR/(Q3 + Q1) of test depletion") {
        tcol = i;
      } else if (ifpV.tokens[i] == "IQR/(Q3 + Q1) of ctrl depletion") {
        ccol = i;
      } else if (ifpV.tokens[i] == "IQR/(Q3 + Q1) of test-ctrl depletion") {
        dcol = i;
      }
    }
    for (size_t i = 1; !ifpV.fail(); ++i) {
      ifpV.readRow(i);
      if (isADouble(ifpV.tokens[tcol])) {
        test_cvs.push_back(convertFromString<float>(ifpV.tokens[tcol]));
      }
      if (isADouble(ifpV.tokens[ccol])) {
        ctrl_cvs.push_back(convertFromString<float>(ifpV.tokens[ccol]));
      }
      if (isADouble(ifpV.tokens[dcol])) {
        diff_cvs.push_back(convertFromString<float>(ifpV.tokens[dcol]));
      }
    }
    float test_median_cvs = vecMedian(test_cvs), 
          ctrl_median_cvs = vecMedian(ctrl_cvs),
          diff_median_cvs = vecMedian(diff_cvs);
    ifpV.resetToBegin();
    Ofstream ofpV(OUT);
    ifpV.readRow(0);
    ofpV << ifpV.line 
         << "\tDist test IQR/(Q3 + Q1) from IQR/(Q3 + Q1) median\tDist ctrl IQR/(Q3 + Q1) from IQR/(Q3 + Q1) median\tDist test-ctrl IQR/(Q3 + Q1) from IQR/(Q3 + Q1) median" 
         << endl;
    for (size_t i = 1; !ifpV.fail(); ++i) {
      ifpV.readRow(i);
      string t_diff_str = "", c_diff_str = "", d_diff_str = "";
      if (isADouble(ifpV.tokens[tcol])) {
        float t_val = convertFromString<float>(ifpV.tokens[tcol]);
        t_diff_str = convertToString<float>(t_val - test_median_cvs);
      }
      if (isADouble(ifpV.tokens[ccol]) ) {
        float c_val = convertFromString<float>(ifpV.tokens[ccol]);
        c_diff_str = convertToString<float>(c_val - ctrl_median_cvs);
      }
      if (isADouble(ifpV.tokens[dcol]) ) {
        float d_val = convertFromString<float>(ifpV.tokens[dcol]);
        d_diff_str = convertToString<float>(d_val - diff_median_cvs);
      }
      ofpV << ifpV.line << "\t" << t_diff_str << "\t" << c_diff_str 
                        << "\t" << d_diff_str << endl;
    }
    ofpV.close();
    ifpV.close();
    system("rm \"" + OUT + ".temp_coeff_var\"");
  }


  if (KILLING_PLOTS != "") {
    system("mv " + OUT + " " + OUT + ".temp_before_killing");
    system("sync");
    system("Rscript --vanilla \"" + EVAL_KILLING + "\" \"" + 
           OUT + ".temp_before_killing\" \"" +
           OUT + ".temp_killing" + "\" \"" + KILLING_PLOTS + "\"");
    if (fileSize(OUT + ".temp_killing") == 
        static_cast<ifstream::pos_type>(-1)) {         
      system("mv " + OUT + ".temp_before_killing " + OUT);
      system("sync");
    } else {
      unordered_map<string, double> killing, lower_kill;
      ReadOnlyTSVFile ifpK(OUT + ".temp_killing", true, true);
      for (size_t i = 1; !ifpK.fail(); ++i) {
        ifpK.readRow(i);
        killing[ifpK.tokens[0]] = convertFromString<double>(ifpK.tokens[1]);
        lower_kill[ifpK.tokens[0]] = convertFromString<double>(ifpK.tokens[2]);
      }
      ifpK.close();
  
      ReadOnlyTSVFile ifpK2(OUT + ".temp_before_killing", true, true);
      vector<double> essential_depletion;
      //    row     non/ess  depl   recall  precision  q3   q3 rec   q3 prec
      Table<size_t, string, double, double, double, double, double, double> all_depletion;
      size_t gene_col = numeric_limits<size_t>::max(),
             depl_col = numeric_limits<size_t>::max(),
             q3_col = numeric_limits<size_t>::max(),
             ess_col = numeric_limits<size_t>::max();
      ifpK2.readRow(0);
      for (size_t i = 0; i < ifpK2.tokens.size(); ++i) {
        if (ifpK2.tokens[i] == "Gene") gene_col = i;
        else if (ifpK2.tokens[i] == "Test median depletion") depl_col = i;
        else if (ifpK2.tokens[i] == "Test Q3") q3_col = i;
        else if (ifpK2.tokens[i] == "Core Ess Or Noness") ess_col = i;
      }
      double totNumEssentials = 0.0;
      for (size_t i = 1; !ifpK2.fail(); ++i) {
        ifpK2.readRow(i);
        double depl_val = convertFromString<double>(ifpK2.tokens[depl_col]);
        double q3_val = convertFromString<double>(ifpK2.tokens[q3_col]);
        all_depletion.push_back(i, ifpK2.tokens[ess_col], depl_val, 0.0, 1.0,
                                                          q3_val, 0.0, 1.0);
        if (ess.find(ifpK2.tokens[gene_col]) != ess.end()) {
          totNumEssentials += 1.0;
          essential_depletion.push_back(depl_val);
        }
      }
      sort(essential_depletion.begin(), essential_depletion.end());
      double min_ess = essential_depletion[0],
             max_ess = essential_depletion[essential_depletion.size() - 1],
             size_ess = static_cast<double>(essential_depletion.size());
  
      // true = reverse sort on median test depletion...
      sortTable<2>(all_depletion, true);
      double cumulative_tp = 0.0, cumulative_fp = 0.0;
      for (size_t i = 0; i < all_depletion.size(); ++i) {
        size_t row = 0;
        string curr_ess = "";
        double depl_val = 0.0, curr_recall = 0.0, curr_precision = 1.0,
               q3_val = 0.0, q3_recall = 0.0, q3_precision =  1.0;
        all_depletion.getRow(i, row, curr_ess, depl_val,
                             curr_recall, curr_precision,
                             q3_val, q3_recall, q3_precision);
        if (curr_ess == "Essential") {
          cumulative_tp += 1.0;
        } else if (curr_ess == "Nonessential") {
          cumulative_fp += 1.0;
        }
        curr_recall = cumulative_tp / totNumEssentials;
        if (cumulative_tp > 0 || cumulative_fp > 0) {
          curr_precision = cumulative_tp/(cumulative_tp +  cumulative_fp);
        }
        all_depletion.setRow(i, row, curr_ess, depl_val,
                             curr_recall, curr_precision,
                             q3_val, q3_recall, q3_precision);
      } 
      // true = reverse sort on Q3 test depletion...
      sortTable<5>(all_depletion, true);
      cumulative_tp = 0.0, cumulative_fp = 0.0;
      for (size_t i = 0; i < all_depletion.size(); ++i) {
        size_t row = 0;
        string curr_ess = "";
        double depl_val = 0.0, curr_recall = 0.0, curr_precision = 1.0,
               q3_val = 0.0, q3_recall = 0.0, q3_precision =  1.0;
        all_depletion.getRow(i, row, curr_ess, depl_val,
                             curr_recall, curr_precision,
                             q3_val, q3_recall, q3_precision);
        if (curr_ess == "Essential") {
          cumulative_tp += 1.0;
        } else if (curr_ess == "Nonessential") {
          cumulative_fp += 1.0;
        }
        q3_recall = cumulative_tp / totNumEssentials;
        if (cumulative_tp > 0 || cumulative_fp > 0) {
          q3_precision = cumulative_tp/(cumulative_tp +  cumulative_fp);
        }
        all_depletion.setRow(i, row, curr_ess, depl_val,
                             curr_recall, curr_precision,
                             q3_val, q3_recall, q3_precision);
      } 
      // Sort table back on row number... 
      sortTable<0>(all_depletion);
    
      ifpK2.resetToBegin();
      Ofstream ofpK(OUT);
      ifpK2.readRow(0);
      ofpK << ifpK2.line << "\tKilling measure\tPercent of essential genes with a lower killing measure\tPercent of essential genes with a lower median test depletion";
      if (PRECISION_RECALL != "") {
        ofpK << "\tMedian test depletion recall\tMedian test depletion precision\tQ3 test depletion recall\tQ3 test depletion precision";
      }
      ofpK << endl;
      for (size_t i = 1; !ifpK2.fail(); ++i) {
        ifpK2.readRow(i);
        string kill_measure = "", pct_lower_kill = "";
        if (killing.find(ifpK2.tokens[0]) != killing.end()) {
          kill_measure = convertToString<double>(killing[ifpK2.tokens[0]]);
          pct_lower_kill = convertToString<double>(lower_kill[ifpK2.tokens[0]]);
        }
        double curr_depl = convertFromString<double>(ifpK2.tokens[depl_col]);
        double pct_low_depl = -100.0;
        if (curr_depl >= max_ess) pct_low_depl = 100.0;
        else if (curr_depl < min_ess) pct_low_depl = 0.0;
        else {
          vector<double>::iterator up = upper_bound(essential_depletion.begin(),
                                                    essential_depletion.end(),
                                                    curr_depl);
          double curr_pos = static_cast<double>(up - essential_depletion.begin());
          pct_low_depl = 100.0*(curr_pos)/size_ess;
        }
        ofpK << ifpK2.line << "\t" << kill_measure << "\t" << pct_lower_kill 
             << "\t" << convertToString<double>(pct_low_depl);
        if (PRECISION_RECALL != "") {
          size_t row = 0;
          string curr_ess = "";
          double depl_val = 0.0, curr_recall = 0.0, curr_precision = 1.0,
                 q3_val = 0.0, q3_recall = 0.0, q3_precision = 1.0;
          all_depletion.getRow(i-1, row, curr_ess, depl_val,
                               curr_recall, curr_precision,
                               q3_val, q3_recall, q3_precision);
          ofpK << "\t" << curr_recall << "\t" << curr_precision
               << "\t" << q3_recall << "\t" << q3_precision;
        }
        ofpK << endl; 
      }
      ofpK.close();
      system("rm \"" + OUT + ".temp_before_killing\" \"" + 
                       OUT + ".temp_killing\"");

      system("Rscript --vanilla \"" + FC_KILLING + "\" \"" +
             OUT + "\" \"" + KILLING_PLOTS + "\" " +
             convertToString<size_t>(TOP_FC_KILL_GENES));
  
      if (PRECISION_RECALL != "") {
        system("Rscript --vanilla \"" + EVAL_PR + "\" \"" + OUT + "\" \"" +
                PRECISION_RECALL + ".median_test_depletion\" 2");
        system("Rscript --vanilla \"" + EVAL_PR + "\" \"" + OUT + "\" \"" +
                PRECISION_RECALL + ".q3_test_depletion\" 3");
      }
    }
  }


  if (CUTOFF != "") {
    string strata = "";
    if (STRATAS != "") {
      strata = "\"" + STRATAS + "\" " + 
               convertToString<size_t>(NUM_STRATA) + " " +
               convertToString<size_t>(STRATA_ADD) + " \"" +
               STRATAS + ".temp_rank_classes\"";
    }
    system("Rscript --vanilla \"" + EVAL_BETA + "\" \"" + OUT + "\" \"" +
            CUTOFF + "\" " + convertToString<size_t>(BETA_TOP_N) + " " +
            convertToString<double>(CUTOFF_P) + " TRUE " + strata);
    if (STRATAS != "" &&
        fileSize(CUTOFF) != static_cast<ifstream::pos_type>(-1) &&
        fileSize(STRATAS + ".temp_rank_classes") != 
        static_cast<ifstream::pos_type>(-1)) {
      ReadOnlyTSVFile ifpC(CUTOFF, true, true);
      ifpC.readRow(1);
      size_t max_rank = convertFromString<size_t>(ifpC.tokens[2]);
      ifpC.close();
      unordered_map<size_t, string> rank_to_class;
      ReadOnlyTSVFile ifpT(STRATAS + ".temp_rank_classes", true, true);
      size_t max_i = 0;
      for (size_t i = 1; !ifpT.fail(); ++i) {
        ifpT.readRow(i);
        rank_to_class[i] = ifpT.tokens[2];
        max_i = i;
      }
      ifpT.close();
      system("rm \"" + STRATAS + ".temp_rank_classes\"");
      system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_rank_classes\"");
      system("sync");
      ReadOnlyTSVFile ifpB(OUT + ".temp_before_rank_classes", true, true);
      ifpB.readRow(0);
      Ofstream ofp2(OUT);
      ofp2 << ifpB.line << "\tJenks class" << endl;
      for (size_t i = 1; !ifpB.fail(); ++i) {
        ifpB.readRow(i); 
        string beyond_max = "";
        if (i > max_rank) {
          beyond_max = "Beyond p = " + convertToString<double>(CUTOFF_P) + " cutoff";
          if (i <= max_i) {
            beyond_max += " - ";
          }
        }
        string jenks_class = "";
        if (i <= max_i) jenks_class = rank_to_class[i];
        ofp2 << ifpB.line << "\t"
             << beyond_max << jenks_class << endl;
      }
      ofp2.close();
      ifpB.close();
      system("rm \"" + OUT + ".temp_before_rank_classes\"");
    } else {
      if (fileSize(CUTOFF) != static_cast<ifstream::pos_type>(-1)) {
        system("rm \"" + CUTOFF + "\"");
      } 
      if (fileSize(STRATAS + ".temp_rank_classes") != 
          static_cast<ifstream::pos_type>(-1)) {
        system("rm \"" + STRATAS + ".temp_rank_classes\"");
      } 
    }
  }


  if (ADD_Z_SCORES && CUTOFF_Z != "") {
    ReadOnlyTSVFile ifpZ1(OUT, true, true);
    ifpZ1.readRow(0);
    size_t rank_col = numeric_limits<size_t>::max();
    size_t nonparamZ_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpZ1.tokens.size(); ++i) {
      if (ifpZ1.tokens[i] == "Rank") {
        rank_col = i;
      } else if (ifpZ1.tokens[i] == "Negative non-parameteric Z-score") {
        nonparamZ_col = i;
      }
    }
    ifpZ1.close();
    // The non-parametric Z-score rank is nonparamZ_col+1

    // First, sort by the rank with the above being 0-base 
    // and sort being 1-base index. 
    string nonpz_col = convertToString<size_t>(nonparamZ_col+2);
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_zsorting\"");
    system("sync");
    system("(head -n 1 \"" + OUT + ".temp_before_zsorting\" && tail -n +2 \"" + OUT + ".temp_before_zsorting\" | sort -t\"	\" -k" + nonpz_col + "," + nonpz_col + " -n) > \"" + OUT + "\"");
    system("rm -f \"" + OUT + ".temp_before_zsorting\"");
    system("sync");
    
    // Next, do stratas as usual.
    string strata = "";
    if (STRATAS_Z != "") {
      strata = "\"" + STRATAS_Z + "\" " + 
               convertToString<size_t>(NUM_STRATA) + " " +
               convertToString<size_t>(STRATA_ADD) + " \"" +
               STRATAS_Z + ".temp_rank_classes\"";
    }
    system("Rscript --vanilla \"" + EVAL_BETA_Z + "\" \"" + OUT + "\" \"" +
            CUTOFF_Z + "\" " + convertToString<size_t>(BETA_TOP_N) + " " +
            convertToString<double>(CUTOFF_P) + " TRUE " + strata);
    if (STRATAS_Z != "" &&
        fileSize(CUTOFF_Z) != static_cast<ifstream::pos_type>(-1) &&
        fileSize(STRATAS_Z + ".temp_rank_classes") !=
        static_cast<ifstream::pos_type>(-1)) {
      ReadOnlyTSVFile ifpCZ(CUTOFF_Z, true, true);
      ifpCZ.readRow(1);
      size_t max_rank = convertFromString<size_t>(ifpCZ.tokens[2]);
      ifpCZ.close();
      unordered_map<size_t, string> rank_to_class;
      ReadOnlyTSVFile ifpTZ(STRATAS_Z + ".temp_rank_classes", true, true);
      size_t max_i = 0;
      for (size_t i = 1; !ifpTZ.fail(); ++i) {
        ifpTZ.readRow(i);
        rank_to_class[i] = ifpTZ.tokens[2];
        max_i = i;
      }
      ifpTZ.close();
      system("rm \"" + STRATAS_Z + ".temp_rank_classes\"");
      system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_rank_classes\"");
      system("sync");

      ReadOnlyTSVFile ifpBZ(OUT + ".temp_before_rank_classes", true, true);
      ifpBZ.readRow(0);
      Ofstream ofp2Z(OUT);
      ofp2Z << ifpBZ.line << "\tJenks class for non-parametric Z-scores" << endl;
      for (size_t i = 1; !ifpBZ.fail(); ++i) {
        ifpBZ.readRow(i); 
        string beyond_max = "";
        if (i > max_rank) {
          beyond_max = "Beyond p = " + convertToString<double>(CUTOFF_P) + " cutoff";
          if (i <= max_i) {
            beyond_max += " - ";
          }
        }
        string jenks_class = "";
        if (i <= max_i) jenks_class = rank_to_class[i];
        ofp2Z << ifpBZ.line << "\t"
              << beyond_max << jenks_class << endl;
      }
      ofp2Z.close();
      ifpBZ.close();
      system("rm \"" + OUT + ".temp_before_rank_classes\"");
    } else {
      if (fileSize(CUTOFF_Z) != static_cast<ifstream::pos_type>(-1)) {
        system("rm \"" + CUTOFF_Z + "\"");
      }
      if (fileSize(STRATAS_Z + ".temp_rank_classes") !=
          static_cast<ifstream::pos_type>(-1)) {
        system("rm \"" + STRATAS_Z + ".temp_rank_classes\"");
      }
    }

    // Finally, re-sort on rank.
    string str_rank_col = convertToString<size_t>(rank_col+1);
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_zsorting\"");
    system("sync");
    system("(head -n 1 \"" + OUT + ".temp_before_zsorting\" && tail -n +2 \"" + OUT + ".temp_before_zsorting\" | sort -t\"	\" -k" + str_rank_col + "," + str_rank_col + " -n) > \"" + OUT + "\"");
    system("rm -f \"" + OUT + ".temp_before_zsorting\"");
  }


  if (PRECISION_RECALL != "" && POSITIVE_CONTROLS != "" &&
      ((ADD_Z_SCORES && CHEMO && CHEMO_PAIRED && CORRELATION_PREFIX != "") ||
       (ADD_Z_SCORES && 
        (!CHEMO || !CHEMO_PAIRED) && CORRELATION_PREFIX != ""))) {
    ReadOnlyTSVFile ifpZ1(OUT, true, true);
    ifpZ1.readRow(0);
    size_t rank_col = numeric_limits<size_t>::max();
    size_t nonparamZ_col = numeric_limits<size_t>::max();
    for (size_t i = 0; i < ifpZ1.tokens.size(); ++i) {
      if (ifpZ1.tokens[i] == "Rank") {
        rank_col = i;
      } else if (ifpZ1.tokens[i] == "Negative non-parameteric Z-score") {
        nonparamZ_col = i;
      }
    }
    ifpZ1.close();
    // The non-parametric Z-score rank is nonparamZ_col+1
 
    // First, sort by the rank with the above being 0-base 
    // and sort being 1-base index. 
    string nonpz_col = convertToString<size_t>(nonparamZ_col+2);
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_zsorting\"");
    system("sync");
    system("(head -n 1 \"" + OUT + ".temp_before_zsorting\" && tail -n +2 \"" + OUT + ".temp_before_zsorting\" | sort -t\"	\" -k" + nonpz_col + "," + nonpz_col + " -n) > \"" + OUT + "\"");
    system("rm -f \"" + OUT + ".temp_before_zsorting\"");
  
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_precision_recall\"");
    system("sync");
    ReadOnlyTSVFile ifpPR(OUT + ".temp_before_precision_recall", true, true);
    ifpPR.readRow(0);
    if (ifpPR.tokens[0] != "Gene") {
      cerr << "Error!  The first column is no longer the Gene column.  This cannot happen.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    Ofstream ofpPR(OUT);
    ofpPR << ifpPR.line << "\tZ-score recall\tZ-score precision" << endl;
    double totNumPos = static_cast<double>(pr_pos.size());
    double cumulative_tp = 0.0, cumulative_fp = 0.0, 
           precision = 1.0, recall = 0.0;
    for (size_t i = 1; !ifpPR.fail(); ++i) {
      ifpPR.readRow(i);
      string& curr_gene = ifpPR.tokens[0];
      if (pr_pos.find(curr_gene) != pr_pos.end()) {
        cumulative_tp += 1.0;
      } else if (noness.find(curr_gene) != noness.end()) {
        cumulative_fp += 1.0;
      }
      recall = cumulative_tp / totNumPos;
      if (cumulative_tp > 0 || cumulative_fp > 0) {
        precision = cumulative_tp/(cumulative_tp +  cumulative_fp);
      }
      ofpPR << ifpPR.line << "\t" << recall << "\t" << precision << endl;
    }
    ofpPR.close();
    ifpPR.close();
    system("rm \"" + OUT + ".temp_before_precision_recall\"");
 
    system("Rscript --vanilla \"" + EVAL_PR + "\" \"" + OUT + "\" \"" +
            PRECISION_RECALL + ".Z_Score\" 1");
   
    // Finally, re-sort on rank.
    string str_rank_col = convertToString<size_t>(rank_col+1);
    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_zsorting\"");
    system("sync");
    system("(head -n 1 \"" + OUT + ".temp_before_zsorting\" && tail -n +2 \"" + OUT + ".temp_before_zsorting\" | sort -t\"	\" -k" + str_rank_col + "," + str_rank_col + " -n) > \"" + OUT + "\"");
    system("rm -f \"" + OUT + ".temp_before_zsorting\"");
  }

  if (ADD_UNIPROT) {
    ReadOnlyTSVFile ifpU(UNIPROT, true, true);
    ifpU.readRow(0);
    bool gene_found = false;
    size_t gene_col = 0;
    for (size_t i = 0; !gene_found && i < ifpU.tokens.size(); ++i) {
      if (ifpU.tokens[i] == "Gene names") {
        gene_col = i;
        gene_found = true;
      }
    }
    if (!gene_found) {
      cerr << "Error!  Gene names column not found in " << UNIPROT << ".  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    vector<string> tokens, blanks;
    for (size_t j = gene_col + 1; j < ifpU.tokens.size(); ++j) {
      tokens.push_back(ifpU.tokens[j]);
      blanks.push_back("");
    }
    string header = joinStrings(tokens);
    string blank_uni = joinStrings(blanks);
    unordered_map<string, string> gene_uniprot;
    for (size_t i = 1; !ifpU.fail(); ++i) {
      ifpU.readRow(i);
      tokens.clear();
      for (size_t j = gene_col + 1; j < ifpU.tokens.size(); ++j) {
        tokens.push_back(ifpU.tokens[j]);
      }
      string curr_uniprot = joinStrings(tokens);
      tokens.clear();
      splitStringOnSpaces(ifpU.tokens[gene_col], tokens);
      for (size_t j = 0; j < tokens.size(); ++j) {
        if ((j == 0) || 
            (gene_uniprot.find(trimSpacesString(tokens[j])) == 
             gene_uniprot.end())) {
          gene_uniprot[trimSpacesString(tokens[j])] = curr_uniprot;
        }
      }
    }
    ifpU.close();

    system("mv \"" + OUT + "\" \"" + OUT + ".temp_before_uniprot\"");
    system("sync");
    ReadOnlyTSVFile ifpUT(OUT + ".temp_before_uniprot", true, true);
    ifpUT.readRow(0);
    if (ifpUT.tokens[0] != "Gene") {
      cerr << "Error!  CCA output file first column is not Gene.  Exiting." << endl;
      exit(EXIT_FAILURE);
    }
    Ofstream ofpUT(OUT);
    ofpUT << ifpUT.line << "\t" << header << endl;
    for (size_t i = 1; !ifpUT.fail(); ++i) {
      ifpUT.readRow(i); 
      ofpUT << ifpUT.line << "\t";
      if (gene_uniprot.find(trimSpacesString(ifpUT.tokens[0])) != gene_uniprot.end()) {
        ofpUT << gene_uniprot[trimSpacesString(ifpUT.tokens[0])] << endl;
      } else {
        ofpUT << blank_uni << endl;
      }
    }
    ofpUT.close();
    ifpUT.close();
    system("rm -f \"" + OUT + ".temp_before_uniprot\"");
  }


  string OUT_DIR = "./";
  if (OUT.find("/") != string::npos) {
    // false below finds the last instance of / in the OUT path.
    OUT_DIR = beforeString(OUT, "/", false);
  }
  // Below is part of the Imagemagick package.
  // mogrify trim will remove all white borders from the *.png files
  system("mogrify -quiet -trim " + OUT_DIR + "/*.png");


  if (RUN_GATHER_QC && CORRELATION_PREFIX != "" &&
      GATHER_QC_HTML != "" && GATHER_QC_TSV != "") {
    string at_home = AT_HOME ? "AT_HOME=true" : "";
    system("ccaGatherQC COUNTS=\"" + COUNTS + "\" CORREL_NORMED_COUNTS=\"" + CORRELATION_PREFIX + ".normed_counts.tsv\" TSV=\"" + GATHER_QC_TSV + "\" HTML=\"" + GATHER_QC_HTML + "\" PREFIX=\"cca_files/\" " + at_home);
  }

  system("mv " + OUT + " " + OUT + ".full.do_not_use");
  system("sync");
  string cca_outdir = OUT_DIR + "/cca_files";
  system("mkdir -p " + cca_outdir);
  cerr << "It is OK, if there are some no such file errors below.  Not all files are always generated, depending on the data and input parameters." << endl;
  if (INITQC_INPUT_PREFIX != "") {
    for (size_t k = 0; k < count_files.size(); ++k) {
      system("mv " + count_files[k] + INITQC_INPUT_PREFIX + "* " + cca_outdir);
      system("sync");
    }
  }
  if (CORRELATION_PREFIX != "") system("mv " + CORRELATION_PREFIX + "* " + cca_outdir);
  if (BARPLOT_HEAD != "") system("mv " + BARPLOT_HEAD + "* " + cca_outdir);
  if (INITQC_PREFIX != "") system("mv " + INITQC_PREFIX + "* " + cca_outdir);
  // if (DETREND_OUT != "") system("mv " + DETREND_OUT + "* " + cca_outdir);
  if (KILLING_PLOTS != "") system("mv " + KILLING_PLOTS + "* " + cca_outdir);
  if (CUTOFF != "") {
    // Below so can catch *tsv and *...png files.
    CUTOFF = CUTOFF.substr(0, CUTOFF.size()-3);
    system("mv " + CUTOFF + "* " + cca_outdir);
  }
  if (CUTOFF_Z != "") {
    // Below so can catch *tsv and *...png files.
    CUTOFF_Z = CUTOFF_Z.substr(0, CUTOFF_Z.size()-3);
    system("mv " + CUTOFF_Z + "* " + cca_outdir);
  }
  if (STRATAS != "") system("mv " + STRATAS + "* " + cca_outdir);
  if (STRATAS == "" && STRATAS_Z != "") system("mv " + STRATAS_Z + "* " + cca_outdir);
  if (PRECISION_RECALL != "") system("mv " + PRECISION_RECALL + "* " + cca_outdir);
  system("sync");
  if (OUT != "") {
    system("cleanUpCca IN=\"" + OUT + ".full.do_not_use\" OUT=\"" + OUT + "\"");
    system("sync");
    system("mv " + OUT + ".full.do_not_use " + cca_outdir);
    system("sync");
  }
  if (GATHER_QC_TSV != "") system("mv " + GATHER_QC_TSV + " " + cca_outdir);

  return 0;
}

