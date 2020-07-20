#include <RcppArmadillo.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <math.h>
using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
 * 
 *                 Mathematical Functions::
 *
 * ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

/*
 std::ostream& operator<<(std::ostream& os, const std::vector<int> &input)
{
for (auto const& i: input) {
os << i << " ";
}
return os;
}

template<typename T>
static inline double Lerp(T v0, T v1, T t)
{
  return (1 - t)*v0 + t*v1;
}

template<typename T>
static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs)
{
  if (inData.empty())
  {
    return std::vector<T>();
  }
  
  if (1 == inData.size())
  {
    return std::vector<T>(1, inData[0]);
  }
  
  std::vector<T> data = inData;
  std::sort(data.begin(), data.end());
  std::vector<T> quantiles;
  
  for (size_t i = 0; i < probs.size(); ++i)
  {
    T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);
    
    size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
    size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));
    
    T datLeft = data.at(left);
    T datRight = data.at(right);
    
    T quantile = Lerp<T>(datLeft, datRight, poi - left);
    
    quantiles.push_back(quantile);
  }
  
  return quantiles;
}
 */


void print(std::vector<double> const &input)
{
  for (int i = 0; i < input.size(); i++) {
    std::cout << input.at(i) << ' ';
  }
}

double sum(const vector<double>& a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += a[i];
  }
  return s;
}

double median(vector<double> vec, bool isSorted)
{
  typedef vector<double>::size_type vec_sz;
  
  vec_sz size = vec.size();
  if (size == 0) throw domain_error("[median]: median of an empty vector");
  // if (size == 0) return(-1.0);
  if (!isSorted) sort(vec.begin(), vec.end());
  
  vec_sz mid = size/2;
  
  return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];
}

double mean(const vector<double>& a)
{
  return sum(a) / a.size();
}

double sqsum(const vector<double>& a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(const vector<double>& nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

vector<double> stats5(vector<double> vec, bool isSorted) {
  typedef vector<double>::size_type vec_sz;
  
  vec_sz size = vec.size();
  if (size == 0) throw domain_error("[vector5]: median of an empty vector");
  if (!isSorted) sort(vec.begin(), vec.end());
  
  vec_sz mid = size/2;
  
  vector<double> Q1_vec(vec.cbegin(), vec.cbegin() + mid);
  vector<double> Q3_vec(vec.cbegin() + mid + 1, vec.cend());
  
  vector<double> stats;
  double Q1  = median(Q1_vec, true);
  double Q2  = median(vec, true);
  double Q3  = median(Q3_vec, true);
  double avg = mean(vec);
  double std = stdev(vec);
  
  stats.push_back(Q1);
  stats.push_back(Q2);
  stats.push_back(Q3);
  stats.push_back(avg);
  stats.push_back(std);
  
  return stats;
}

vector<double> stats4(vector<double> vec, bool isSorted) {
  typedef vector<double>::size_type vec_sz;
  
  double Q1=-1.0, Q3=-1.0, Q31=-1.0, Q2=-1.0, avg=-1.0, std=-1.0;
  vector<double> stats;
  
  vec_sz size = vec.size();
  if (size==0) {
    stats.push_back(Q31);
    stats.push_back(Q2);
    stats.push_back(avg);
    stats.push_back(std);
    
    return stats;
  }
  if (size == 0) throw domain_error("[stats4]: median of an empty vector");
  if (!isSorted) sort(vec.begin(), vec.end());
  
  vec_sz mid = size/2;
  
  vector<double> Q1_vec(vec.cbegin(), vec.cbegin() + mid);
  vector<double> Q3_vec(vec.cbegin() + mid + 1, vec.cend());
  
  if (!Q1_vec.empty() && !Q3_vec.empty()) {
    Q1  = median(Q1_vec, true);
    Q3  = median(Q3_vec, true);
    Q31 = Q3-Q1;
  }
  if (!vec.empty()) {
    Q2  = median(vec, true);
    avg = mean(vec);
    std = stdev(vec);
  }
  
  stats.push_back(Q31);
  stats.push_back(Q2);
  stats.push_back(avg);
  stats.push_back(std);
  
  return stats;
}

vector<double> quantile3(vector<double> vec, bool isSorted) {
  typedef vector<double>::size_type vec_sz;
  
  vec_sz size = vec.size();
  if (size == 0) throw domain_error("[quantile3]: median of an empty vector");
  if (!isSorted) sort(vec.begin(), vec.end());
  
  vec_sz mid = size/2;
  
  vector<double> Q1_vec(vec.cbegin(), vec.cbegin() + mid);
  vector<double> Q3_vec(vec.cbegin() + mid + 1, vec.cend());
  
  vector<double> quantiles;
  double Q1 = median(Q1_vec, true);
  double Q2 = median(vec, true);
  double Q3 = median(Q3_vec, true);
  
  quantiles.push_back(Q1);
  quantiles.push_back(Q2);
  quantiles.push_back(Q3);
  
  return quantiles;
}

vector<double> operator-(const vector<double>& a, double b)
{
  vector<double> retvect;
  for (int i = 0; i < a.size(); i++)
  {
    retvect.push_back(a[i] - b);
  }
  return retvect;
}

vector<double> operator*(const vector<double>& a, const vector<double>& b)
{
  vector<double> retvect;
  for (int i = 0; i < a.size() ; i++)
  {
    retvect.push_back(a[i] * b[i]);
  }
  return retvect;
}

// [[Rcpp::export]]
double pearsoncoeff(const vector<double>& X, const vector<double>& Y)
{
  return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

// [[Rcpp::export]]
int C_n_choose_m(int n, int r)
{
  if (r == 0) return 1;
  
  // Extra computation saving for large R,
  // using property:
  // N choose R = N choose (N-R)
  if (r > n / 2) return C_n_choose_m(n, n - r); 
  
  int res = 1; 
  // long res = 1; 
  
  for (int k = 1; k <= r; ++k)
  {
    res *= n - k + 1;
    res /= k;
  }
  return res;
}

// [[Rcpp::export]]
vector<double> C_getDeltaBetaVariationVecAbs(const vector<double>& a, const vector<double>& b,
                                             double minDelta2=-1.0, double minDelta1=-1.0, double minDelta0=-1.0) {
  vector<double> dB_vec;
  vector<double> ret_vec;
  
  // if (a.empty() || b.empty()) throw domain_error("vectors must be of equal length!");
  
  double dB_median=-1.0, dB_mean=-1.0, dB_stdev=-1.0;
  int db_cnt2 = 0;
  int db_cnt1 = 0;
  int db_cnt0 = 0;
  for (int ii=0; ii<a.size(); ii++) {
    double dB = fabs(a[ii]-b[ii]);
    dB_vec.push_back(dB);
    
    if (minDelta2>0.0 && fabs(dB)>minDelta2) db_cnt2++;
    if (minDelta1>0.0 && fabs(dB)>minDelta1) db_cnt1++;
    if (minDelta0>0.0 && fabs(dB)>minDelta0) db_cnt0++;
  }
  if (!dB_vec.empty()) {
    dB_median = median(dB_vec, false);
    dB_mean   = mean(dB_vec);
    dB_stdev  = stdev(dB_vec);
  }
  
  ret_vec.push_back( dB_median );
  ret_vec.push_back( dB_mean );
  ret_vec.push_back( dB_stdev );
  
  if (minDelta2>0.0) ret_vec.push_back( db_cnt2 );
  if (minDelta1>0.0) ret_vec.push_back( db_cnt1 );
  if (minDelta0>0.0) ret_vec.push_back( db_cnt0 );

  return ret_vec;
}

// [[Rcpp::export]]
vector<double> C_getDeltaBetaVariationVecDif(const vector<double>& a, const vector<double>& b) {
  
  vector<double> dB_vec;
  vector<double> ret_vec;
  
  double dB_median=-1.0, dB_mean=-1.0, db_pos_per=0.0;
  double dB_stdev=-1.0;
  int dB_pos_cnt=0, dB_tot_cnt=0;
  
  // if (a.empty() || b.empty()) throw domain_error("vectors must be of equal length!");
  
  for (int ii=0; ii<a.size(); ii++) {
    double dB = b[ii]-a[ii];
    if (dB>0.0) {
      dB_pos_cnt++;
      dB_vec.push_back(dB); 
    }
    dB_tot_cnt++;
  }
  
  if (!dB_vec.empty()) {
    dB_median = median(dB_vec, false);
    dB_mean   = mean(dB_vec);
    
    // Adding back starndard deviation (sd)
    dB_stdev  = stdev(dB_vec);
    
    db_pos_per = (double)(1000*dB_pos_cnt/dB_tot_cnt)/10;
  }
  
  ret_vec.push_back( dB_median );
  ret_vec.push_back( dB_mean );
  
  // Adding back starndard deviation (sd)
  ret_vec.push_back( dB_stdev );
  
  ret_vec.push_back( db_pos_per );
  ret_vec.push_back( dB_tot_cnt );
  
  return ret_vec;
}

// [[Rcpp::export]]
vector<double> C_getDeltaBetaVariationVecDifNP(const vector<double>& a, const vector<double>& b) {
  
  vector<double> dBN_vec;
  vector<double> dBP_vec;
  vector<double> ret_vec;
  
  double dBN_median=-1.0, dBN_mean=-1.0, dBP_median=-1.0, dBP_mean=-1.0;
  // double dBN_stdev=-1.0, dBP_stdev=-1.0;
  int dBN_cnt=0, dBP_cnt=0;
  
  if (a.empty() || b.empty()) throw domain_error("vectors must be of equal length!");
  
  for (int ii=0; ii<a.size(); ii++) {
    double dB = b[ii]-a[ii];
    if (dB<0.0) {
      dBN_cnt++;
      dBN_vec.push_back(dB); 
    } else {
      dBP_cnt++;
      dBP_vec.push_back(dB); 
    }
  }
  
  if (dBN_vec.size() != 0) {
    dBN_median = median(dBN_vec, false);
    dBN_mean   = mean(dBN_vec);
    // dBN_stdev  = stdev(dBN_vec);
  }
  
  if (dBP_vec.size() != 0) {
    dBP_median = median(dBP_vec, false);
    dBP_mean   = mean(dBP_vec);
    // dBP_stdev  = stdev(dBP_vec);
  }
  
  ret_vec.push_back( dBN_median );
  ret_vec.push_back( dBN_mean );
  // ret_vec.push_back( dBN_stdev );
  ret_vec.push_back( dBN_cnt );
  
  ret_vec.push_back( dBP_median );
  ret_vec.push_back( dBP_mean );
  // ret_vec.push_back( dBP_stdev );
  ret_vec.push_back( dBP_cnt );
  
  return ret_vec;
}

// [[Rcpp::export]]
NumericMatrix na_matrix(int n, int m){
  NumericMatrix mat(n,m) ;
  std::fill( mat.begin(), mat.end(), NumericVector::get_na() ) ;
  return mat ;
}

// [[Rcpp::export]]
vector<double> NumericVectorToStlVector(NumericVector nv, int start, int end) {
  vector<double> v;
  
  for (int ii = start; ii < end; ii++) {
    // Rcout << "\tii=" << ii << ", nv[ii]=" << nv[ii] << "\n";
    if (!NumericMatrix::is_na(nv[ii]))
      v.push_back(nv[ii]);
  }
  sort(v.begin(), v.end());

  return v;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
 * 
 *                 Export R Functions::
 *
 * ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

/*
 * Under Construction:: Basic RSquared Test
 *
 */
// [[Rcpp::export]]
NumericMatrix C_crossSampleLociRSquared(NumericMatrix x, NumericVector idxVec, CharacterVector samVec, 
                                        Nullable<CharacterVector> cpgVec=R_NilValue, 
                                        double minDelta2=0.2, double minDelta1=0.1, double minDelta0=0.05,
                                        bool all=false, bool cmb=false,
                                        int verbose=0, int vt=1) {
  
  string funcTag = "C_crossSampleLociRSquared";  
  int sampleCnt = idxVec.size();
  int dat_nrow = x.nrow(), dat_ncol = x.ncol();
  // int out_nrows_per_sample = 11; // Previously 9 without delta-beta threshold counts, NOW 11 (not 12), combined Q3/Q1 into Q3-Q1
  int out_nrows_per_sample = 12; // Adding back NA%, Previously 9 without delta-beta threshold counts, NOW 11 (not 12), combined Q3/Q1 into Q3-Q1
  
  // Initialize Numeric Output Matrix
  int out_ncols = (out_nrows_per_sample*sampleCnt) + 1;
  if (verbose>=vt) Rcout << "Matrix Dim: nrow=" << dat_nrow << ", ncol=" << dat_ncol << ", sampleCnt=" << sampleCnt
                         << ", out_nrows_per_sample=" << out_nrows_per_sample << ", out_ncols=" << out_ncols << "\n";
  
  // int cmb_out_cnt = 4; // Previously it was * 3 now its * 8, actually * 6 remove SD, NOW just 4 (q2,mu,pass,tot)
  int cmb_out_cnt = 5; // Previously it was * 3 now its * 8, actually * 6 retain SD, NOW back to 5 (q2,mu,sd,pass,tot)
  
  if (all) {
    out_ncols += ((sampleCnt-0) * (sampleCnt-1)) * cmb_out_cnt;
  } else if (cmb) {
    out_ncols += C_n_choose_m(sampleCnt, 2) * cmb_out_cnt;
  } else {
    out_ncols += (sampleCnt-1) * cmb_out_cnt;
  }
  if (verbose>=vt) Rcout << "Matrix Dim: nrow=" << dat_nrow << ", ncol=" << dat_ncol << ", sampleCnt=" << sampleCnt
                         << ", out_nrows_per_sample=" << out_nrows_per_sample << ", out_ncols=" << out_ncols << "\n";

  NumericMatrix out_mat = na_matrix(dat_nrow, out_ncols);
  
  // We only need to calculate this once (for the first locus)
  vector<string> ssamples_key_vec;
  vector<string> csamples_key_vec;
  
  //  Loop Over Each Loci Sequentiallly::
  for (int lociIdx = 0; lociIdx < dat_nrow; lociIdx++) {
    vector<vector<double> >  ssamples_dat_vecs;
    vector<vector<double> >  ssamples_sum_vecs;
    
    vector<double> csample_combX_vec;
    vector<double> csample_combY_vec;
    
    if (verbose>=vt+1)
      Rcout << "LociIdx=" << lociIdx << ", samVec.size=" << samVec.size() << "\n";

    int col_int_idx = 0;
    int out_col_idx = 1;
    for (int sampleIdx = 0; sampleIdx < sampleCnt; sampleIdx++) {
      int col_beg_idx = col_int_idx;
      int col_end_idx = col_beg_idx + idxVec[sampleIdx];
      int tot_rep_cnt = idxVec[sampleIdx];
      int beg_col_idx = out_col_idx;
      string sampleName = (string)samVec[sampleIdx];

      if (verbose>=vt+2)
        Rcout << "\tSampleIdx=" << sampleIdx << ", Sample=" << sampleName << ", beg/end=" << col_beg_idx << "/" << col_end_idx << "\n";
      
      // Local summary and data vectors::
      //  NOTE:: We don't need the following because we'll just directly save it to the stack
      //         vector<double> ssamples_dat_vec;
      // ssamples_dat_vecs.push_back( NumericVectorToStlVector( x( lociIdx, _ ), col_beg_idx, col_end_idx) );
      
      // Really need to remove all NA's before moving forward
      //   This next bit of code replaces the above code:: ssamples_dat_vecs.push_back( NumericVectorToStlVector( x( lociIdx, _ ), col_beg_idx, col_end_idx) );
      //   Ideallly this removes all the NA's from any future vectors...
      //   Testing Now...
      vector<double> non_nan_vec;
      for (int col_idx = col_beg_idx; col_idx < col_end_idx; col_idx++) {
        if (NumericMatrix::is_na(x(lociIdx,col_idx))) {
        } else {
          non_nan_vec.push_back((double)x(lociIdx, col_idx));
        }
      }
      ssamples_dat_vecs.push_back( non_nan_vec );

      vector<double> ssample_combA_vec;
      vector<double> ssample_combB_vec;
      
      int nan_rep_cnt = tot_rep_cnt - ssamples_dat_vecs.back().size();

      // Build single sample summary vector
      vector<double> sample_s4_vec = stats4(ssamples_dat_vecs.back(), true);

      if (verbose>=vt+2)
        Rcout << "\tSampleIdx=" << sampleIdx << ", Sample=" << sampleName << ", beg/end=" << col_beg_idx << "/" << col_end_idx << ". stats4 Completed!\n";
      
      if (sample_s4_vec[0]==-1.0 || sample_s4_vec[1]==-1.0 || sample_s4_vec[2]==-1.0 || sample_s4_vec[3]==-1.0) {
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
      } else {
        out_mat(lociIdx,out_col_idx++) = sample_s4_vec[0];
        out_mat(lociIdx,out_col_idx++) = sample_s4_vec[1];
        out_mat(lociIdx,out_col_idx++) = sample_s4_vec[2];
        out_mat(lociIdx,out_col_idx++) = sample_s4_vec[3];
      }

      // Add beta summary keys::
      if (lociIdx==0 && sampleIdx==0) {
        ssamples_key_vec.push_back( "beta_31" );
        ssamples_key_vec.push_back( "beta_q2" );
        ssamples_key_vec.push_back( "beta_mu" );
        ssamples_key_vec.push_back( "beta_sd" );
      }

      // create ssample_combA/B_vec
      int val_rep_cnt = ssamples_dat_vecs.back().size();
      for (int ii=0; ii<val_rep_cnt; ii++) {
        for (int jj=ii+1; jj<val_rep_cnt; jj++) {
          ssample_combA_vec.push_back(ssamples_dat_vecs.back().at(ii));
          ssample_combB_vec.push_back(ssamples_dat_vecs.back().at(jj));
        }
      }
      
      // Calculate inner single sample delta-beta
      vector<double> sample_dB_vec = C_getDeltaBetaVariationVecAbs(ssample_combA_vec, ssample_combB_vec,
                                                                   minDelta2, minDelta1, minDelta0);
      
      if (verbose>=vt+2) {
        Rcout << "\tSampleIdx=" << sampleIdx << ", Sample=" << sampleName << ", beg/end=" << col_beg_idx << "/" << col_end_idx 
              << ". C_getDeltaBetaVariationVecAbs Completed!\n";
      }

      if (sample_dB_vec[0]==-1.0 || sample_dB_vec[1]==-1.0 || sample_dB_vec[2]==-1.0 ||
          sample_dB_vec[3]==-1.0 || sample_dB_vec[4]==-1.0 || sample_dB_vec[5]==-1.0) {
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
        out_mat(lociIdx,out_col_idx++) = NA_REAL;
      } else {
        out_mat(lociIdx,out_col_idx++) = sample_dB_vec[0];
        out_mat(lociIdx,out_col_idx++) = sample_dB_vec[1];
        out_mat(lociIdx,out_col_idx++) = sample_dB_vec[2];
        
        out_mat(lociIdx,out_col_idx++) = sample_dB_vec[3];
        out_mat(lociIdx,out_col_idx++) = sample_dB_vec[4];
        out_mat(lociIdx,out_col_idx++) = sample_dB_vec[5];
      }
      out_mat(lociIdx,out_col_idx++) = nan_rep_cnt;
      out_mat(lociIdx,out_col_idx++) = (double)(100*nan_rep_cnt) / (double)tot_rep_cnt;
      
      // Add delta-beta summary keys::
      if (lociIdx==0 && sampleIdx==0) {
        ssamples_key_vec.push_back( "dB_q2" );
        ssamples_key_vec.push_back( "dB_mu" );
        ssamples_key_vec.push_back( "dB_sd" );
        
        ssamples_key_vec.push_back( "dB2_cnt" );
        ssamples_key_vec.push_back( "dB1_cnt" );
        ssamples_key_vec.push_back( "dB0_cnt" );
        ssamples_key_vec.push_back( "nan_cnt" );
        ssamples_key_vec.push_back( "nan_per" );
      }
      
      // Push to cross sample combination data stack (RSquared Stack)::
      csample_combX_vec.insert(csample_combX_vec.end(), ssample_combA_vec.begin(), ssample_combA_vec.end());
      csample_combY_vec.insert(csample_combY_vec.end(), ssample_combB_vec.begin(), ssample_combB_vec.end());
        
      if (verbose>=vt+2) {
        Rcout << "\tSample=" << sampleName << ", tot_rep_cnt=" << tot_rep_cnt
              << ", csample_combX_vec.size=" << csample_combX_vec.size()
              << ", csample_combY_vec.size=" << csample_combY_vec.size() << "\n";
        Rcout << "\tSample=" << sampleName
              << ", sample_dat.size="  << ssamples_dat_vecs.size()
              << ", sample_vec.size="  << ssamples_sum_vecs.size()
              << "\n";
        Rcout << "\tSample=" << sampleName << ", Summary Stats 1:\n";
        for (int ii=beg_col_idx; ii<out_col_idx; ii++) {
          Rcout << " " << ssamples_key_vec.at(ii-beg_col_idx) << "[" <<ii<< "]=" << out_mat(lociIdx,ii);
        }
        Rcout << "\n";
        Rcout << "------------------------------------------\n\n\n";
      }
      
      col_int_idx = col_end_idx;
    }
    double crossSampleR2 = 0;
    if (csample_combX_vec.size() > 1 && csample_combY_vec.size() > 1) 
      crossSampleR2 = pearsoncoeff(csample_combX_vec, csample_combY_vec);
    if (verbose>=vt) Rcout << "crossSampleR2=" << crossSampleR2 << ", ssamples_key_vec=" << ssamples_key_vec.size() 
                           << ""
                           << "\n\n";
    
    // Storing outputs::
    out_mat(lociIdx,0) = crossSampleR2;
    
    // Generate delta betas cross all sample combinations::
    for (int sampleIdxX=0; sampleIdxX<sampleCnt; sampleIdxX++) {
      string sampleNameX = (string)samVec[sampleIdxX];
      
      int y_init = sampleIdxX+1;
      int y_next = y_init+1;
      if (all) {
        y_init = 0;
        y_next = sampleCnt;
      }
      /* *** This part below seems like it was missing... Possible doesn't work! *** */
      if (cmb) {
        y_init = sampleIdxX+1;
        y_next = sampleCnt;
      }
      if (y_next>sampleCnt) y_next=sampleCnt;
      
      for (int sampleIdxY=y_init; sampleIdxY<y_next; sampleIdxY++) {
        string sampleNameY = (string)samVec[sampleIdxY];
        
        string combSampleName = sampleNameX+"_"+sampleNameY;
        
        // TBD:: Need to remove diag, use choose function for output initialization::
        // if (true) {
        if (sampleIdxX!=sampleIdxY) {
          // create ssample_combA/B_vec
          int val_cur_cntA = ssamples_dat_vecs[sampleIdxX].size();
          int val_cur_cntB = ssamples_dat_vecs[sampleIdxY].size();
          
          vector<double> csample_combA_vec;
          vector<double> csample_combB_vec;

          for (int aa=0; aa<val_cur_cntA; aa++) {
            for (int bb=0; bb<val_cur_cntB; bb++) {
              csample_combA_vec.push_back(ssamples_dat_vecs[sampleIdxX][aa]);
              csample_combB_vec.push_back(ssamples_dat_vecs[sampleIdxY][bb]);
            }
          }
          
          // Calculate cross sample delta-beta
          vector<double> sample_dB_vec = C_getDeltaBetaVariationVecDif(csample_combA_vec, csample_combB_vec);
          if (sample_dB_vec[0]==-1.0 || sample_dB_vec[1]==-1.0 || sample_dB_vec[2]==-1.0 || sample_dB_vec[3]==-1.0) {
            out_mat(lociIdx,out_col_idx++) = NA_REAL;
            out_mat(lociIdx,out_col_idx++) = NA_REAL;
            out_mat(lociIdx,out_col_idx++) = NA_REAL;
            out_mat(lociIdx,out_col_idx++) = NA_REAL;
            out_mat(lociIdx,out_col_idx++) = NA_REAL;
          } else {
            out_mat(lociIdx,out_col_idx++) = sample_dB_vec[0];
            out_mat(lociIdx,out_col_idx++) = sample_dB_vec[1];
            out_mat(lociIdx,out_col_idx++) = sample_dB_vec[2];
            out_mat(lociIdx,out_col_idx++) = sample_dB_vec[3];
            out_mat(lociIdx,out_col_idx++) = sample_dB_vec[4];
          }

          if (lociIdx==0) {
            csamples_key_vec.push_back( combSampleName+"_CSS_q2" );
            csamples_key_vec.push_back( combSampleName+"_CSS_mu" );
            
            // Adding back standard deviation (sd)
            csamples_key_vec.push_back( combSampleName+"_CSS_sd" );
            
            csamples_key_vec.push_back( combSampleName+"_CSS_per" );
            csamples_key_vec.push_back( combSampleName+"_CSS_cnt" );
          }
        }
        
      }
    }

    if (verbose>=vt+1) {
      Rcout << "Loci=" << lociIdx << ", DONE"
            << ", crossSampleR2=" << crossSampleR2
            << ", sample_dat.size="  << ssamples_dat_vecs.size()
            << ", sample_vec.size="  << ssamples_sum_vecs.size()
            << ", samVec.size=" << samVec.size() 
            << "\n--------------------------------------------------------------\n\n\n";
    }
    ssamples_dat_vecs.clear();
    ssamples_sum_vecs.clear();
  }
  
  // Gernerate Output Column Names::
  vector<string> out_cols_vec;
  out_cols_vec.push_back("RSquared");
  for (int sampleIdx = 0; sampleIdx < sampleCnt; sampleIdx++) {
    string sampleName = (string)samVec[sampleIdx];
    for (int ii=0; ii < ssamples_key_vec.size(); ii++) {
      out_cols_vec.push_back( sampleName+"_"+ssamples_key_vec[ii] );
    }
  }
  for (int ii=0; ii<csamples_key_vec.size(); ii++) {
    out_cols_vec.push_back( csamples_key_vec[ii] );
  }
  
  Rcout << "\tsampleCnt=" << sampleCnt << "\n"
        << "\tssamples_key_vec.size=" << ssamples_key_vec.size() << "\n"
        << "\tSize out_cols_vec.size=" << out_cols_vec.size() << "\n\n";
  colnames(out_mat) = wrap(out_cols_vec);
  rownames(out_mat) = rownames(x);

  if (verbose>=vt) Rcout << "["<<funcTag<<"]: Done.\n" << "\n";
  
  return out_mat;
}





/*
 * Default Test Code Provided Below::
 * 
 */


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

timesTwo(42)

*/
