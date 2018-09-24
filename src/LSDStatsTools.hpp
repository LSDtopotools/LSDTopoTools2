//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStatsTools
// Land Surface Dynamics StatsTools
//
// A collection of statistical routines for use with the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include <map>
#include <math.h>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;


#ifndef StatsTools_H
#define StatsTools_H

// tools for getting keys from a map
vector<string> extract_keys(map<string, int> input_map);
vector<string> extract_keys(map<string, float> input_map);
vector<string> extract_keys(map<string, bool> input_map);
vector<string> extract_keys(map<string, string> input_map);
vector<string> extract_keys(map<string, double> input_map);

// tools for reversing arrays
Array2D<double> reverse_array_rows(Array2D<double>& data);
Array2D<double> reverse_array_cols(Array2D<double>& data);
Array2D<float> reverse_array_rows(Array2D<float>& data);
Array2D<float> reverse_array_cols(Array2D<float>& data);
Array2D<int> reverse_array_rows(Array2D<int>& data);
Array2D<int> reverse_array_cols(Array2D<int>& data);

// computes linear regression
// replaces data in residuals with residuals and returns a 4 element vector, which has slope, intercept, r^2 and
// the Durbin-Watson test statistic which looks for autocorrelation of the residuals
vector<float> simple_linear_regression(vector<float>& x_data, vector<float>& y_data, vector<float>& residuals);
float get_mean(vector<float>& y_data);
float get_mean_ignore_ndv(vector<float>& y_data, float ndv);
float get_mean_ignore_ndv(Array2D<float>& data, float ndv);
float get_median(vector<float> y_data);
float get_median(vector<float> y_data, float ndv);
float get_median_sorted(vector<float> sorted_y_data);
float get_median_absolute_deviation(vector<float> y_data, float median);
vector<float> get_IQR_and_median(vector<float> y_data);
float get_SST(vector<float>& y_data, float mean);
float get_variance_ignore_ndv(Array2D<float>& data, float ndv, float mean);
float get_range_ignore_ndv(Array2D<float>& data, float ndv);
float get_range_from_vector(vector<float>& y_data, float ndv);
float Get_Minimum(vector<float>& y_data, float ndv);
int Get_Minimum(vector<int>& y_data, float ndv);
vector<int> Get_Index_Minimum(vector<int>& y_data, float ndv);
float Get_Maximum(vector<float>& y_data, float ndv);
vector<int> Get_Index_Maximum(vector<float>& y_data, float ndv);
float get_durbin_watson_statistic(vector<float> residuals);
float get_standard_deviation(vector<float>& y_data, float mean);
float get_standard_deviation(vector<float>& y_data, float mean, float ndv);
float get_standard_error(vector<float>& y_data, float standard_deviation);
vector<float> get_common_statistics(vector<float>& y_data);
vector<float> calculate_descriptive_stats(vector<float>& data);
float get_percentile(vector<float>& data, float percentile);


// sort a vector of vector in regards to a first vector, they all need the same number of element
// BG - some days in Januray 2018
vector<vector<float> > sort_vectors_from_one(vector<float> to_sort, vector<vector<float> > follow_the_sort);
// reorganise a vector from a vector of new IDx
vector<float> reorganize_vector_from_new_idx(vector<float> vecval, vector<int> vecid);



// orthogonal regression
// 01/04/2017 SMM No foolin
// This comes from davegiles.blogspot.co.uk/2014/11/orthogonal-regression-first-steps.html
// NOTE: THis is more generally called Total Least Squares
//  There is a solution using matrices that is probably compuationally faster
//  Might want to implement that in the future if this is slow
//  Note R^2 from simple linear regression
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> orthogonal_linear_regression( vector<float>& x_data, vector<float>& y_data, float& intercept, float& gradient, float& R_squared);

// this function gets the difference between nearest neighbours in a vector of y data
// FJC 10/11/15
vector<float> difference(vector<float>& y_data);

// this function gets the main peaks from a vector of y data using the first order difference
// FJC 13/11/15
void get_peak_indices(vector<float>& y_data, float threshold, int distance, vector<int>& peak_indices);

// sorts data; produces quartile-quantile comparison against standard normal variate, returning
// an (evenly spaced) sorted subsample of N_points, their corresponding normal variate and the
// reference value  from the standard normal distribution.  Test for departures from normality
// within the given distribution.
void generate_q_q_plot(vector<float>& data, vector<float>& values, vector<float>& standard_normal_variates, vector<float>& mn_values, int N_points);

// declaration of the quantile_quantile analysis
void quantile_quantile_analysis(vector<float>& data, vector<float>& values, vector<float>& standard_normal_variates, vector<float>& mn_values, int N_points);

// declaration of the quantile_quantile analysis
// modified to pass in percentiles as arguments
void quantile_quantile_analysis_defined_percentiles(vector<float>& data, vector<float>& values, vector<float>& standard_normal_variates, vector<float>& mn_values, int N_points, int lower_percentile, int upper_percentile);

// Bootstrapping of linear regressions
// N_iterations is the number of bootstrap iterations
// acceptance probablility is the probability that you will accept any given data point
// in an iteration. This runs without replacement
// Returns summary statistics (see cpp code for details)
vector<float> bootstrap_linear_regression(vector<float>& x_data, vector<float>& y_data, int N_iterations, float acceptance_prob);


// calculates least squares linear regression for two datasets, returning
// gradient and intercept of regression line, alongside the R-squared value.
// DTM 07/10/2014
void least_squares_linear_regression(vector<float> x_data, vector<float> y_data, float& intercept, float& gradient, float& R_squared);
// take a slice of a vector
// DTM 30/10/2014
vector<float> slice_vector(vector<float>::iterator first,vector<float>::iterator last);


// interpolation
double interp1D_ordered(vector<double>& x, vector<double>& y, double x_interp_loc);
vector<double> interp1D_ordered(vector<double>& x, vector<double>& y, vector<double> x_interp_loc);
float interp1D_ordered(vector<float>& x, vector<float>& y, float x_interp_loc);
vector<float> interp1D_ordered(vector<float>& x, vector<float>& y, vector<float> x_interp_loc);
vector<double> interp1D_spline_ordered(vector<double>& x_data, vector<double>& y_data,
                                       vector<double>& x_interp_locs);
float interp1D_unordered(vector<float> x, vector<float> y, float x_interp_loc);
vector<float> interp1D_unordered(vector<float> x, vector<float> y, vector<float>& x_interp_loc);
double interp1D_unordered(vector<double> x, vector<double> y, double x_interp_loc);
vector<double> interp1D_unordered(vector<double> x, vector<double> y, vector<double>& x_interp_loc);
vector<double> interp1D_spline_unordered(vector<double> x_data, vector<double> y_data,
                                       vector<double>& x_interp_locs);
double interp2D_bilinear(vector<double>& x_locs, vector<double>& y_locs, Array2D<double> data,
                        double x_interp, double y_interp);
float interp2D_bilinear(vector<float>& x_locs, vector<float>& y_locs, Array2D<float> data,
                        float x_interp, float y_interp);

// Generate spline curves from X and Y vectors of floats
Array2D<float> CalculateCubicSplines(vector<float> X, vector<float> Y);
void PlotCubicSplines(vector<float> X, vector<float> Y, int SplineResolution, vector<float>& Spline_X, vector<float>& Spline_Y);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// simple cubic spline interpolation library without external
// dependencies
//
// ---------------------------------------------------------------------
// Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------
//
//
// SEE: http://kluge.in-chemnitz.de/opensource/spline/
//
// band matrix solver. This is for the spline
class band_matrix
{
  private:
     std::vector< std::vector<double> > m_upper;  // upper band
     std::vector< std::vector<double> > m_lower;  // lower band
  public:
     band_matrix() {};                             // constructor
     band_matrix(int dim, int n_u, int n_l);       // constructor
     ~band_matrix() {};                            // destructor
     void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
     int dim() const;                             // matrix dimension
     int num_upper() const
     {
       return m_upper.size()-1;
     }
     int num_lower() const
     {
       return m_lower.size()-1;
     }
     // access operator
     double & operator () (int i, int j);            // write
     double   operator () (int i, int j) const;      // read
     // we can store an additional diogonal (in m_lower)
     double& saved_diag(int i);
     double  saved_diag(int i) const;
     void lu_decompose();
     std::vector<double> r_solve(const std::vector<double>& b) const;
     std::vector<double> l_solve(const std::vector<double>& b) const;
     std::vector<double> lu_solve(const std::vector<double>& b,
                              bool is_lu_decomposed=false);
};

// spline interpolation
// -----------------------
//
// USAGE:
//
// spline s;
// s.set_points(X,Y);    // currently it is required that X is already sorted
// double x=1.5;
// cout << "spline at " << x << " is: " << s(x) <<endl;
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
class spline
{
  private:
   std::vector<double> m_x,m_y;           // x,y coordinates of points
   // interpolation parameters
   // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
   std::vector<double> m_a,m_b,m_c,m_d;
  public:
   void set_points(const std::vector<double>& x,
                   const std::vector<double>& y, bool cubic_spline=true);
   double operator() (double x) const;
};
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// calculate the imaginary error function using trapezoid rule integration
double erfi(double tau);

// these look for linear segments within a data series.
void populate_segment_matrix(int start_node, int end_node, float no_data_value,
                vector<float>& all_x_data, vector<float>& all_y_data, int maximum_segment_length,
                float sigma, Array2D<float>& like_array, Array2D<float>& m_array,
                Array2D<float>& b_array, Array2D<float>& rsquared_array,
                Array2D<float>& DW_array);
void calculate_segment_matrices(vector<float>& all_x_data, vector<float>& all_y_data, int maximum_segment_length,
                float sigma, Array2D<float>& like_array, Array2D<float>& m_array,
                Array2D<float>& b_array, Array2D<float>& rsquared_array,
                Array2D<float>& DW_array);
void find_max_like_of_segments(int minimum_segment_length, Array2D<float>& like_array,
                vector<float>& max_MLE, vector< vector<int> >& segments_for_each_n_segments);
void find_max_AIC_of_segments(int minimum_segment_length, vector<float>& all_x_data, vector<float>& all_y_data,
                Array2D<float>& like_array,
                vector<float>& max_MLE, vector<float>& AIC_of_segments,
                vector<float>& AICc_of_segments, vector< vector<int> >& segments_for_each_n_segments);
void calculate_AIC_of_segments_with_normalized_sigma(float sigma,
                vector<float>& one_sigma_max_MLE, vector<float>& all_x_data,
                vector<float>& AIC_of_segments,vector<float>& AICc_of_segments);

// the below function is the main driver for the segment fitting code.
void best_fit_driver_AIC_for_linear_segments(int minimum_segment_length, float sigma,
                      vector<float> all_x_data, vector<float> all_y_data,
                      vector<float>& max_MLE);

// this gets the number of segments for several different values of sigma bassed in the vector sigma_values
void get_n_segments_for_various_sigma(vector<float> sigma_values, vector<float> one_sig_max_MLE,
                    vector<float>& all_x_data,
                      vector<int>& best_fit_AIC, vector<int>& best_fit_AICc,
                      vector< vector<float> >& AIC_for_each_n_segments,
                      vector< vector<float> >& AICc_for_each_n_segments);

// this prints full AIC and AICc information to screen
void print_AIC_and_AICc_to_screen(vector<float> sigma_values, vector< vector<int> > segments_for_each_n_segments,
                      vector<int> best_fit_AIC, vector<int> best_fit_AICc,
                      vector< vector<float> > AIC_for_each_n_segments,
                      vector< vector<float> > AICc_for_each_n_segments);

// this prints information about the most likeley segments to screen
void print_to_screen_most_likeley_segment_lengths( vector< vector<int> > segments_for_each_n_segments,
                    vector<float> MLE_for_segments);

// this returns the m, b, r2 and DW stats of each segment
void get_properties_of_best_fit_segments(int bestfit_segments_node, vector< vector<int> >& segments_for_each_n_segments,
                     vector<float>& m_values, Array2D<float>& m_array,
                     vector<float>& b_values, Array2D<float>& b_array,
                     vector<float>& r2_values, Array2D<float>& rsquared_array,
                     vector<float>& DW_values, Array2D<float>& DW_array);

// these functions manipulate likelihood matrices and vectors for use with the segment tool
Array2D<float> normalize_like_matrix_to_sigma_one(float sigma, Array2D<float>& like_array);
vector<float> normalize_like_vector_to_sigma_one(float sigma, vector<float> like_vector);
Array2D<float> change_normalized_like_matrix_to_new_sigma(float sigma, Array2D<float>& sig1_like_array);
vector<float> change_normalized_like_vector_to_new_sigma(float sigma, vector<float> sig1_like_vector);

// this uses a moving window to find segments and is incomplete
void find_linear_segments(vector<float>& all_x_data, vector<float>& all_y_data, int segment_length);


// functions for combinations
void combinations(vector<int> v, int start, int n, int k, int maxk);
void combinations(vector<int> v, int start, int n, int k, int maxk, vector< vector<int> >& combovecvec);
vector< vector<int> > combinations (int n, int k, bool zero_indexed);

// functions for partitioning and permutation (to be used with linear segment finding
int partitions_min( int x, int y);
void partition_print(int t, vector<int>& p);
void partitions_with_minimum_length(int n, int k, int t, int min_length, vector<int>& p);
void partitions_with_minimum_length(int n, int k, int t, int min_length, vector<int>& p,
                vector< vector < vector<int> > >& partitions);
void integer_partition(int n, int k, int t, vector<int>& p);
void partition_driver_to_screen(int n, int minimum_length);
vector< vector < vector<int> > > partition_driver_to_vecvecvec(int k, int minimum_length);
void partition_assign(int t, vector<int>& p, vector< vector < vector<int> > >& partitions);
void partition_vecvecvec_print(vector< vector < vector<int> > >& partitions);
void partition_vecvecvec_print_with_permutation(vector< vector < vector<int> > >& partitions);
void permute_partitioned_integer_vector(vector<int> permute_vector);

// this generates random segments for use in testing the segment finding algorithm
void generate_random_segments(float sigma, int minimum_n_nodes, int mean_segment_length, int segment_range,
                    float dx, float offset_range, float m_range,
                    vector<float>& x_data, vector<float>& y_data,
                    vector<int>& segment_length, vector<float>& slope, vector<float>& intercept);


// maxiumum likihood estimators
float calculate_MLE(vector<float>& measured, vector<float>& modelled, vector<float>& sigma);
float calculate_MLE(vector<float>& measured, vector<float>& modelled, float sigma);
float calculate_MLE_from_residuals(vector<float>& residuals, float sigma);

// RMSE estimator
float calculate_RMSE_from_residuals(vector<float>& residuals);

// a random number generator
float ran3( long *idum );
// Randomly sample from a vector without replacement DTM 21/04/2014
vector<float> sample_without_replacement(vector<float> population_vector, int N);
vector<int> sample_without_replacement(vector<int> population_vector, int N);

// conversion from numbers to strings
string itoa(int num);
string dtoa(float num);
bool atobool(string value);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function returning a Gaussian random number
// DAV 16/10/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double Gauss_rand(int Nrand, double GaussAdd, double GaussFac);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Use the Marsaglia polar method to generate random numbers drawn from a normal distribution
// with a given mean and minimum values. Set allowNegative to false to stop output values from
// dropping below zero.
//
// Extreme values can fall below or above the boundaries in < 3 sigma of cases.
//
// Seed for random number will fail post 2038. I will instruct my firstborn to resolve this
// problem.
//
// SWDG 9/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float getGaussianRandom(float minimum, float mean, bool allowNegative);

// Log binning module
// two overloaded functions:
//    -> for data stored in a 2D array (e.g. slope-area)
void log_bin_data(Array2D<float>& InputArrayX, Array2D<float>& InputArrayY, float log_bin_width, vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>& StandardDeviationX_output, vector<float>& StandardDeviationY_output,
                      vector<float>& StandardErrorX_output, vector<float>& StandardErrorY_output, vector<int>& num_observations, float NoDataValue);

//    -> for data stored in a 1D vector (e.g. for spectral analysis)
void log_bin_data(vector<float>& InputVectorX, vector<float>& InputVectorY, float log_bin_width,
                  vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>&  StandardDeviationX_output,
                      vector<float>&  StandardDeviationY_output, int NoDataValue);

// Regular binning algoritm for data stored in a 1D vector
void bin_data(vector<float>& InputVectorX, vector<float>& InputVectorY, float bin_width,
                  vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>& MedianY_output,
                      vector<float>&  StandardDeviationX_output, vector<float>&  StandardDeviationY_output,
                      vector<float>& StandardErrorX_output, vector<float>& StandardErrorY_output,
                      vector<int>& number_observations_output, float& bin_lower_limit, float NoDataValue);

// Regular binning algoritm for data stored in a 1D vector Similar to above but spits out more stats
void bin_data(vector<float>& InputVectorX, vector<float>& InputVectorY, float bin_width,
              vector<float>& midpoints_output, vector<float>&  MeanX_output,
              vector<float>&  MedianX_output, vector<float>&  StandardDeviationX_output,
              vector<float>& StandardErrorX_output, vector<float>& MADX_output,
              vector<float>& MeanY_output, vector<float>& MinimumY_output,
              vector<float>& FirstQuartileY_output, vector<float>& MedianY_output,
              vector<float>& ThirdQuartileY_output, vector<float>& MaximumY_output,
              vector<float>&  StandardDeviationY_output, vector<float>& StandardErrorY_output,
              vector<float>& MADY_output, vector<int>& number_observations_output,
              float NoDataValue);


//look for empty bins output from the log binning function and removes them to avoid
//plotting several empty bins at 0,0 in some cases. SWDG 6/11/13
void RemoveSmallBins(vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>& StandardDeviationX_output, vector<float>& StandardDeviationY_output,
                      vector<float>& StandardErrorX_output, vector<float>& StandardErrorY_output, vector<int>& number_observations, float bin_threshold);

// Load in a vector of data and convert into a histogram with a specified bin width
// that is printed to file containing:
//    Midpoint LowerLim UpperLim Count ProbabilityDensity
void print_histogram(vector<float> input_values, float bin_width, string filename);
// improved histogram functions
void calculate_histogram(vector<float> input_values, float bin_width, vector<float>& Midpoints, vector<float>& LLims, vector<float>& ULims, vector<int>& Count, vector<float>& ProbabilityDensity);
void calculate_histogram_fixed_limits(vector<float> input_values, float bin_width, float lower_limit, float upper_limit, vector<float>& Midpoints, vector<float>& LLims, vector<float>& ULims, vector<int>& Count, vector<float>& ProbabilityDensity);


// This is a much simpler version of the binning software.  It takes two vectors, and
// sorts the values held within the first vector into bins according to their respective
// values in the second vector.  The output is a vector<vector> with the binned dataset.
// and a vector of bin midpoints.  These can then be analysed ahd plotted as desired.
// DTM 14/04/2014
void bin_data(vector<float>& vector1, vector<float>& vector2, float min, float max, float bin_width, vector<float>& mid_points, vector< vector<float> >& binned_data);
// log_bin_data
// This is a similar version for log-binning
// DTM 30/10/2014
void log_bin_data(vector<float>& vector1, vector<float>& vector2, float log_bin_width, vector<float>& bin_mid_points, vector<float>& bin_vector1_mean, vector<float>& bin_vector2_mean, vector< vector<float> >& binned_data, const float NoDataValue = -9999);

// tools for sorting
template<class T> struct index_cmp;
void matlab_double_sort(vector<double>& unsorted, vector<double>& sorted, vector<size_t>& index_map);
void matlab_double_reorder(std::vector<double> & unordered, std::vector<size_t> const & index_map, std::vector<double> & ordered);
void matlab_float_sort(vector<float>& unsorted, vector<float>& sorted, vector<size_t>& index_map);
void matlab_float_reorder(std::vector<float> & unordered, std::vector<size_t> const & index_map, std::vector<float> & ordered);
void matlab_float_sort_descending(vector<float>& unsorted, vector<float>& sorted, vector<size_t>& index_map);
void matlab_int_sort(vector<int>& unsorted, vector<int>& sorted, vector<size_t>& index_map); // added 27/11/13 SWDG
void matlab_int_reorder(std::vector<int> & unordered, std::vector<size_t> const & index_map, std::vector<int> & ordered);



//Get vector of unique values in an input array of ints
vector<int> Unique(Array2D<int> InputArray, int NoDataValue);

//Get vector of unique values in an input array of floats
vector<float> Unique(Array2D<float> InputArray, int NoDataValue);

//Return unique values from a vector of ints.
//Wrapper around the std library unique method which also resizes the output vector.
// SWDG - 22/7/16
vector<int> Unique(vector<int> InputVector);

//Return unique values from a vector of float.
//Wrapper around the std library unique method which also resizes the output vector.
// SWDG - 22/7/16
vector<float> Unique(vector<float> InputVector);

// Generate vector of evenly spaced numbers between two points
vector<float> linspace(float min, float max, int n);

// convert degree bearing from north to radians from east
float BearingToRad(float Bearing);

// conversion from degrees to radians
float rad(float degree);
double rad(double degree);

// conversion from radians to degrees
float deg(float radians);
double deg(double radians);

// Get the angle between two vectors
float angle_between_vectors(float x1, float y1, float x2, float y2);

// Get the clockwise angle between two vectors
float clockwise_angle_between_vector_and_north(float x1, float y1, float x2, float y2);

// get clockwise angle between two vectors specifying the origin
float clockwise_angle_between_two_vectors(float x0, float y0, float x1, float y1, float x2, float y2);

// Get the angle between two vectors in radians
// We need to calculate the (x1,y1) and (x2,y2) coordinates by moving
// the vectors to intercept (0,0)
// the bool vectors_point_downstream is true if the vector's first element is the
// upstream node in a channel and false if the first node is downstream.
float angle_between_two_vector_datasets(vector<float>& x1_data, vector<float>& y1_data,
                                        vector<float>& x2_data, vector<float>& y2_data,
                                        bool vectors_point_downstream);

// This function takes x and y data as vectors and returns a 2 element vector
// where the 0 element is the x1 component of a directional vector
// and the 1 element is the y1 component of a directional vector
// vector vector vector, Victor.
vector<float> get_directional_vector_coords_from_dataset(vector<float> x1_data, vector<float>& y_data,
                      bool vectors_point_downstream);

// Get the data for a boxplot from an unsorted vector of floats, which does not
// contain any NDV values.
//
// Returns a vector which contains (in this order):
//
// 2Percentile 25Percentile median mean 75Percentile 98Percentile minimum maximum
//
// SWDG 12/11/15
vector<float> BoxPlot(vector<float> data);

//Method to generate Statistical distribution. - DTM
void get_distribution_stats(vector<float>& y_data, float& mean, float& median, float& UpperQuartile, float& LowerQuartile, float& MaxValue);

// Method to calculate the quadratic mean. - DTM
double get_QuadraticMean(vector<double> input_values, double bin_width);

// basic parser for parameter files   JAJ  08/01/2014
// There may be a better place to put this, but I can't think where
void parse_line(ifstream &infile, string &parameter, string &value);

// Method to get the maximum value in a 2D array - SWDG 12/6/14
float Get_Maximum(Array2D<float> Input, float NDV);
float Get_Maximum(Array2D<int> Input, float NDV);
int Get_Maximum_Index(Array2D<float> Input, int NDV);
int Get_Maximum_Index(Array2D<int> Input, int NDV);

// Method to get the maximum value in a 2D array - MDH 27/8/14
float Get_Minimum(Array2D<float> Input, int NDV);
float Get_Minimum(Array2D<int> Input, int NDV);
int Get_Minimum_Index(Array2D<float> Input, int NDV);
int Get_Minimum_Index(Array2D<int> Input, int NDV);

//Routine to count the number of values in an array - MDH 27/8/14
int Get_Value_Count(Array2D<float> Input, int NDV);
int Get_Value_Count(Array2D<int> Input, int NDV);

//Method to flatten a 2D array into a 1D vector
//generates a vector in row major order
//SWDG 12/6/14
vector<float> Flatten(Array2D<float> Input);
vector<float> Flatten_Without_Nodata(Array2D<float> Input, float NDV);
vector<int> Flatten(Array2D<int> Input);
vector<int> Flatten_Without_Nodata(Array2D<int> Input, float NDV);

//Method to count the number of instances of a given value in an array
//SWDG 17/6/14
int CountValue(Array2D<int> Input, int Value);
int CountValue(Array2D<float> Input, float Value);


//Method used to generate a Kolmogorov-Smirnov statistic and p value
//from numerical recipes
//Data1 and Data2 must be sorted.
//d is the KS statistic value
//p is the p-value. In numerical recipes it is provided as a value subtracted from
//1, this code has been modified to present value the without this subtraction
//so that it matches the result from the scipy implementation. SWDG 26/6/14
void KStwo(vector<float> Data1, vector<float> Data2, float& d, double& p);
float PKS(float z);

// gets the value of a normal distribution at a point x
// mu is mean
// sigma is standard deviation
// x is the point you want the normal distribution evaluated
float NormalDistributionAtX(float mu, float sigma, float x);

// this gets the p value of a normal distribution for a given Z
float pValueNormalDistribution(float Z);

// This is a function to perform the Mann-Whitney U test, a nonparametric
// test that checks if two data sets have the same median
float MannWhitneyUTest(vector<float>& sampleA, vector<float>& sampleB);

// this takes a sorted vector and then finds the normalised ranks (that is
// if data elements are the same the ranks take an average rank)
// It replaces two vectors passed to it
void rank_vector_with_groups(vector<float> sorted_data,
                             vector<float>& ranks, vector<int>& number_in_groups);

// Given a filestream object, read the file into memory and return
// it as a string. From: http://www.cplusplus.com/forum/general/58945/
// No error handling.
// SWDG 16/07/14
string ReadTextFile(ifstream& File);

// This reads a csv file and takes the headers out.
// These headers can't have spaces since the spaces are removed.
// SMM 18/11/2016
vector<string> ReadCSVHeader(string path, string fname);

/// Splits a string delimited by a character, c, into a sequence of strings, here
/// stored in a vector, v.
/// @author DAV, but taken out of C++ Cookbook (Stevens, Digins, Turkanis, and Coswell. O'Reilly)
void split_delimited_string(const string& s, char c, vector<string>& v);

// THis gets the size of a file
// SMM 16/10/2015
int get_file_size(string filename);

//Takes an integer vector of data and an integer vector of key values and
//returns a map of the counts of each value tied to its key.
//
//Assumes that key_values contains all of the values in Data.
//eg Data should be flattened with NoDataValues excluded and
//Key_Values should be created using Unique(Data)
//SWDG 5/6/15
void Count_Instances(vector<int> Data, vector<int> Key_Values, map<int,int>& DataMap);


// test if a string is a float - http://stackoverflow.com/a/447307/1627162
// Added by SWDG on 18/7/16
bool isFloat(string myString);

// removes control characters from the end of strings.
// This is necessary when people use a DOS file format, which
// stupidly adds control characters to the end of lines.
string RemoveControlCharactersFromEndOfString(string toRemove);

// removes all control characters
string RemoveControlCharacters(string toRemove);

// removes spaces
string RemoveSpaces(string toRemove);

// fix the path (adds a slash to end)
string FixPath(string PathtoFix);

// Unix format path
string ReformatPath(string old_path);

// INVERSE ERROR FUNCTIONS AND INVERSE COMPLEMENTARY ERROR FUNCTIONS
// DTM, Following Press et al.,2007; Numerical Recipes, the Art of Scientific Computing, CUP
// Inverse Complementary error function.  Returns x such that erfc(x)=p within limits 0<p<2
float inverfc(float p);
// Inverse Complementary error function.  Returns x such that erf(x)=p within limits -1<p<1
float inverf(float p);

float StabilityIndex(float s, float a, float c1, float c2, float t1, float t2,
                     float x1,float x2, float r1, float r2, float fs1, float fs2);

float f2s(float x1, float x2, float y1, float y2, float z);
float f3s(float x1, float x2, float y1, float y2, float b1, float b2, float z);
float fa(float y1, float y2, float b1, float b2, float a);
float fai(float y1, float y2, float b1, float b2, float a);
float fai2(float y1, float y2, float b1, float b2, float a);
float fai3(float y1, float y2, float b1, float b2, float a);
float fai4(float y1, float y2, float b1, float b2, float a);
float fai5(float y1, float y2, float b1, float b2, float a);

// CODE FOR DISJOINT SET STRUCTURE
struct DSnode{
  DSnode *parent;
  int i_node;
  int rank;
};

class DisjointSet{
private:
  vector<DSnode *> DSnodes;
  int elements;
  int sets;
public:
  DisjointSet();
  ~DisjointSet();
  void DSMakeSet(int i);
  DSnode* Find(DSnode* node);
  void Union(DSnode* node_1, DSnode* node_2);
  int Union_return_label(DSnode* node_1, DSnode* node_2);
  DSnode* get_DSnode(int i);
  int get_parent(int i);
  int ElementCount();
  int SetCount();
  int Reduce();
  void Reset();
};


struct tm Parse_time_string(string time_string);

//Returns the distance between 2 pairs of raster indexes
//SWDG 19/1/17
float distbetween(int row1, int col1, int row2, int col2);

// Normalize the values of an array of floats to between 0 and MaxValue.
// pass in percentiles eg 98 for the 98th percentile to truncate the data
// about the median. For no truncation pass in 0 and 100.
// SWDG 25/1/17
Array2D<float> normalize_terrain_index(Array2D<float> Data, float lower_percentile, float upper_percentile, float MaxValue, float NoDataValue);

// Implementation of the Jordan Curve theorem to test if a given point is inside
// a polygon.
// returns an integer counting the number of times a ray traced from the point (XCoord,YCoord)
// crosses the border of the polygon.
// An even return value (0 is even) means the point is outside the polygon, and an odd
// value means the point is inside the polygon.
//
// Adapted from: http://stackoverflow.com/a/2922778/1627162
//SWDG - 25/1/17
int PointInPolygon(int VertexCount, float XCoords[], float YCoords[], float XCoord, float YCoord);


vector<float> get_value_from_map_and_node(vector<int> vecnode, map<int,float>& map_int_float);


// Impementation of outlier detection algorithms based on the MAD
// BG - 08/01/2018
vector<float> get_absolute_deviation(vector<float> vecval, float NDV);
float get_MAD(vector<float> vecval, float NDV);
vector<float> get_modified_z_score(vector<float> vecval,float NDV);
vector<int> is_outlier_MZS(vector<float> vecval, float NDV, float threshold);




// Implementation of the Kernel Density estimation method from a vector of float
// I am using this review paper about it for the implementation:
// Sheather 2004 - DOI 10.1214/088342304000000297
// I may try to find a recent one but this last is quite well cited and post 2000 and clear ( I don't want to be a SHEATER ahah, I am not sure if this can be consider as a joke but I am laugthing)
//
// This is the fully automated version, an attempt to provide a non parametric KDE estimation
//
// Work in progress, like a lot
// BG - 04/01/2018  - Bonne annee

pair<float,vector<float> > auto_KDE(vector<float> vpoint);
vector<float> gaussian_KDE(vector<float> vpoint, float h);


//-------------------------------------------------------------------
// The code was written by Vikas C. Raykar
// and is copyrighted under the Lessr GPL:
//
// Copyright (C) 2006 Vikas C. Raykar
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA.
//
// The author may be contacted via email at: vikas(at)cs(.)umd(.)edu
//-------------------------------------------------------------------

//-------------------------------------------------------------
// File    : UnivariateDensityDerivative.h
// Purpose : Header file for UnivariateDensityDerivative.cpp
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : September 17, 2005
//-------------------------------------------------------------
// Fast implementation of the r^{th} kernel density derivative
// estimate based on the Gaussian kernel.
// [HANDLES ONLY UNIVARIATE CASE]
//
// Data is assumed to be scaled to the unit interval [0 1].
//
// Implementation based on:
//
// V. C. Raykar and R. Duraiswami 'Very fast optimal bandwidth
// selection for univariate kernel density estimation'
// Technical Report CS-TR-4774, Dept. of Computer
// Science, University of Maryland, College Park.
// ------------------------------------------------------------
//
// INPUTS [7]
// ----------------
// NSources     --> number of sources, N.
// MTargets     --> number of targets, M.
// pSources     --> pointer to sources, px(N).
// pTargets       --> pointer to the targets, py(M).
// Bandwidth    --> the source bandwidth, h.
// Order          --> order of the derivative, r.
// epsilon        --> desired error, eps.
//
// OUTPUTS [1]
// ----------------
// pDensityDerivative --> pointer the the evaluated Density
//             Derivative, pD(M).
//-------------------------------------------------------------------


// Adapted into LSDTT by B.G. - January 2018

class UnivariateDensityDerivative{
  public:
    //constructor
    UnivariateDensityDerivative(int NSources,
      int MTargets,
      double *pSources,
      double *pTargets,
      double Bandwidth,
        int Order,
      double epsilon,
      double *pDensityDerivative);

    //destructor
    ~UnivariateDensityDerivative();

    //function to evaluate the Density Derivative
    void Evaluate();


    //function to evaluate the Hermite polynomial.
    double hermite(double x, int r);

  private:
    int N;        //number of sources.
    int M;        //number of targets.
    double *px;     //pointer to sources, (N).
    double *py;       //pointer to the targets, (M).
    double  h;      //the source bandwidth.
    int r;              //the rth density derivative.
    double eps;         //the desired error
    double *pD;         //pointer the the evaluated Density Derivative, (M).

    double rx;
    double rr;
    double ry;
    int K;
    int p;
    double h_square;
    double two_h_square;

    double *pClusterCenter;
    int *pClusterIndex;

    int num_of_a_terms;
    double *a_terms;

    int num_of_B_terms;
    double *B_terms;

    double pi;
    double q;

    int factorial(int n);
    void choose_parameters();
    void space_sub_division();
    void compute_a();
    void compute_B();



};


vector<double> TV1D_denoise_v2(vector<double> input,  double lambda);

#endif
