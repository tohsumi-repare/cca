#ifndef MOLBIOLIB_REPARE
#define MOLBIOLIB_REPARE

/** \file repare.hpp
 * Has routines specific to applications in Projects/Repare.
 */

#include "statistics.h"
#include "alglib.h"
#include <math.h>


float vecMean(vector<float>& input) {
  if (input.size() == 0) return 0.0;
  float sum = 0.0;
  for (size_t i = 0; i < input.size(); ++i) sum += input[i];
  return (sum/static_cast<float>(input.size()));
}
double vecMean(vector<double>& input) {
  if (input.size() == 0) return 0.0;
  double sum = 0.0;
  for (size_t i = 0; i < input.size(); ++i) sum += input[i];
  return (sum/static_cast<double>(input.size()));
}

float vecStdev(vector<float>& input, float& mean) {
  if (input.size() < 2) return 0.0;
  float sqrs = 0.0;
  for (size_t i = 0; i < input.size(); ++i) {
    sqrs += (input[i] - mean)*(input[i] - mean);
  }
  return sqrtf(sqrs/static_cast<float>(input.size() - 1));
}
double vecStdev(vector<double>& input, double& mean) {
  if (input.size() < 2) return 0.0;
  double sqrs = 0.0;
  for (size_t i = 0; i < input.size(); ++i) {
    sqrs += (input[i] - mean)*(input[i] - mean);
  }
  return sqrt(sqrs/static_cast<float>(input.size() - 1));
}

float vecMedian(vector<float>& input) {
  size_t n = input.size();
  if (n == 0) return 0.0;
  sort(input.begin(), input.end());
  if ((n % 2) == 1) {
    return input[(n-1)/2];
  } else {
    size_t index = (n + 1)/2;
    return (input[index] + input[index-1])/2.0;
  }
}
double vecMedian(vector<double>& input) {
  size_t n = input.size();
  if (n == 0) return 0.0;
  sort(input.begin(), input.end());
  if ((n % 2) == 1) {
    return input[(n-1)/2];
  } else {
    size_t index = (n + 1)/2;
    return (input[index] + input[index-1])/2.0;
  }
}

template <typename T>
T tableMedian(Table<T>& input) {
  size_t n = input.size();
  if (n == 0) return 0.0;
  sortTable<0>(input);
  if ((n % 2) == 1) {
    return input[(n-1)/2];
  } else {
    size_t index = (n + 1)/2;
    return (input[index] + input[index-1])/2.0;
  }
}

// Below is from
// https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
float vecMad(vector<float>& input, float& median) {
  size_t n = input.size();
  if (n == 0) return 0.0;
  vector<float> temp;
  for (size_t i = 0; i < n; ++i) {
    temp.push_back(fabsf(input[i] - median));
  }
  return vecMedian(temp);
}
double vecMad(vector<double>& input, double& median) {
  size_t n = input.size();
  if (n == 0) return 0.0;
  vector<double> temp;
  for (size_t i = 0; i < n; ++i) {
    temp.push_back(fabs(input[i] - median));
  }
  return vecMedian(temp);
}


template <typename T>
T tableMad(Table<T>& input, T& median) {
  size_t n = input.size();
  if (n == 0) return 0.0;
  Table<T> temp;
  for (size_t i = 0; i < n; ++i) {
    temp.push_back(fabs(input[i] - median));
  }
  return tableMedian<T>(temp);
}

// Below from https://en.wikipedia.org/wiki/Median_absolute_deviation
// and https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
float zScore(float score, float mean, float stdev, bool median) {
  if (!median) {
    return (score - mean)/stdev;
  } else {
    // Of course, in this case, mean = median and stdev = mad.
    return 0.67449*(score - mean)/stdev;
  }
}



void mannWhitneyUTest(vector<double>& sample1, vector<double>& sample2,
                      size_t& n1, size_t& n2, 
                      double& U1, double& U2, 
                      double& z1, double& z2,
                      double& mU, double& sigma) {
  n1 = sample1.size();
  n2 = sample2.size();

  //    value   sample1? 
  Table<double, bool> data;
  for (size_t i = 0; i < n1; ++i) data.push_back(sample1[i], true);
  for (size_t i = 0; i < n2; ++i) data.push_back(sample2[i], false);
  sortTable<0>(data);
  size_t n = data.size();

  // Use below to hold ranks.  Easy to get mean ranks then.
  //             value          ranks
  unordered_map< double, vector<double> > ranks;
  // Initially set rank to be the row ignoring ties.
  for (size_t i = 0; i < n; ++i) {
    ranks[data.getElement<0>(i)].push_back(static_cast<double>(i+1));
  }
  unordered_map<double, double> mean_rank;
  double sum_diff_ti3_min_ti = 0.0;
  for (auto i : ranks) {
    if (i.second.size() == 1) {
      mean_rank[i.first] = i.second[0];
    } else {
      size_t curr_ti = i.second.size();
      sum_diff_ti3_min_ti += (curr_ti*curr_ti*curr_ti) - curr_ti;
      mean_rank[i.first] = vecMean(i.second);
    }
  }
  // Now correct ties of values by using mean rank.
  // Note mean rank = median rank too in this case since it is sequential ranks.
  double R1 = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double curr_rank = mean_rank[data.getElement<0>(i)];
    if (data.getElement<1>(i)) R1 += curr_rank;
  }


  U1 = R1 - static_cast<double>(n1*(n1 + 1))/2.0;
  mU = static_cast<double>(n1*n2)/2.0;
  U2 = (2.0*mU) - U1;   // 2*mU = n1*n2.  Also U1 + U2 = n1*n2.
  sigma = sqrt((mU/6.0)*
                      (static_cast<double>(n+1) - 
                       sum_diff_ti3_min_ti/static_cast<double>(n*(n-1))));
  z1 = (U1 - mU)/sigma;
  z2 = (U2 - mU)/sigma;
} 



// We're using this test based on
// https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
void mannWhitneyU(Table<double>& mut, Table<double>& wt,
                  double& both_tails,
                  double& left_tail, double& right_tail)
{
  alglib::real_1d_array mut_arr, wt_arr;
  int mut_size = static_cast<int>(mut.size());
  int wt_size = static_cast<int>(wt.size());
  if (mut_size < 5 || wt_size < 5) {
    both_tails = -2.0;
    left_tail = -2.0;
    right_tail = -2.0;
    return;
  }
  mut_arr.setlength(mut_size);
  wt_arr.setlength(wt_size);
  for (int i = 0; i < mut_size; ++i) {
    mut_arr[i] = mut[i];
  }
  for (int i = 0; i < wt_size; ++i) {
    wt_arr[i] = wt[i];
  }
  alglib::mannwhitneyutest(mut_arr, mut_size, wt_arr, wt_size,
                           both_tails, left_tail, right_tail);
}

void mannWhitneyU(vector<float>& mut, vector<float>& wt,
                  double& both_tails,
                  double& left_tail, double& right_tail)
{
  alglib::real_1d_array mut_arr, wt_arr;
  int mut_size = static_cast<int>(mut.size());
  int wt_size = static_cast<int>(wt.size());
  if (mut_size < 5 || wt_size < 5) {
    both_tails = -2.0;
    left_tail = -2.0;
    right_tail = -2.0;
    return;
  }
  mut_arr.setlength(mut_size);
  wt_arr.setlength(wt_size);
  for (int i = 0; i < mut_size; ++i) {
    mut_arr[i] = mut[i];
  }
  for (int i = 0; i < wt_size; ++i) {
    wt_arr[i] = wt[i];
  }
  alglib::mannwhitneyutest(mut_arr, mut_size, wt_arr, wt_size,
                           both_tails, left_tail, right_tail);
}


static inline double lerp(double v0, double v1, double t)
{
    return (1 - t)*v0 + t*v1;
}

// Below based on
// https://stackoverflow.com/questions/11964552/finding-quartiles
// Using the Matlab version.
void tableQ123(vector<double>& quantile,
               double& median, double& firstQuantile, double& thirdQuantile) {

  if (quantile.size() == 0) {
    // cerr << "Warning!  Quantile size is zero.  Returning all zeros." << endl;
    median = 0.0;
    firstQuantile = 0.0;
    thirdQuantile = 0.0;
    return;
  } else if (quantile.size() == 1) {
    median = quantile[0];
    firstQuantile = quantile[0];
    thirdQuantile = quantile[0];
    return;
  } else {
    size_t q_size = quantile.size();
    // firstQuantile
    double poi = lerp(-0.5, q_size - 0.5, 0.25);
    size_t left = max(static_cast<size_t>(floor(poi)), static_cast<size_t>(0));
    size_t right = min(static_cast<size_t>(ceil(poi)), q_size - 1);
    double dataLeft = quantile[left];
    double dataRight = quantile[right];
    firstQuantile = lerp(dataLeft, dataRight, poi - left);
    // median
    poi = lerp(-0.5, q_size - 0.5, 0.5);
    left = max(static_cast<size_t>(floor(poi)), static_cast<size_t>(0));
    right = min(static_cast<size_t>(ceil(poi)), q_size - 1);
    dataLeft = quantile[left];
    dataRight = quantile[right];
    median = lerp(dataLeft, dataRight, poi - left);
    // thirdQuantile
    poi = lerp(-0.5, q_size - 0.5, 0.75);
    left = max(static_cast<size_t>(floor(poi)), static_cast<size_t>(0));
    right = min(static_cast<size_t>(ceil(poi)), q_size - 1);
    dataLeft = quantile[left];
    dataRight = quantile[right];
    thirdQuantile = lerp(dataLeft, dataRight, poi - left);
    return;
  }
}
void tableQ123(vector<float>& theVec,
               double& median, double& firstQuantile, double& thirdQuantile) {
  vector<double> quantile;
  for (size_t i = 0; i < theVec.size(); ++i) {
    quantile.push_back(static_cast<double>(theVec[i]));
  }
  tableQ123(quantile, median, firstQuantile, thirdQuantile);
}
void tableQ123(Table<double>& theTable,
               double& median, double& firstQuantile, double& thirdQuantile) {
  vector<double>& quantile = get<0>(theTable.theData);
  tableQ123(quantile, median, firstQuantile, thirdQuantile);
}



// Below based on
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
template<typename... T>
double tableColumnSkew(Table<T...>& theTable) {
  typedef typename Table<T...>::cols_type cols_type;
  typename tuple_element<0, cols_type>::type& theColumn = get<0>(theTable.theData);
  double N = static_cast<double>(theColumn.size());
  if (N == 0.0) return 0.0;
  double mean = tableColumnMean<0>(theTable),
         s = tableColumnStdDev<0>(theTable);
  double result = 0.0;
  for (size_t i = 0; i < theColumn.size(); ++i) {
    result += (theColumn[i] - mean)*(theColumn[i] - mean)*(theColumn[i] - mean)/N;
  }
  result = result/(s*s*s);
  result = sqrt(N*(N-1.0))*result/(N - 2.0);
  return result;
}


double galtonSkew(double& Q1, double& Q2, double& Q3) {
  if (Q1 == Q3) return 0.0;
  return (Q1 + Q3 - 2.0*Q2)/(Q3 - Q1);
}





#endif
