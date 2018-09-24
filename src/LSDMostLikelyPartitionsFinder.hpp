//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDMostLikeleyPartitionsFinder
// Land Surface Dynamics MostLikeleyPartitionsFinder
//
// An object for extracting segments from x,y data
//  developed for the University of Edinburgh
//  Land Surface Dynamics group topographic toolbox.
//
// This object is mainly used in the analysis of channel profiles
//  transformed using the integral method of channel analysis.
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

/** @file LSDMostLikelyPartitionsFinder.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 1.0.0
@brief This object looks for the most likeley partitions or segments of 2D data.
@details Is principally used to identify segments of differing channel steepness in chi-zeta space.

@date 03/01/2013
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include "LSDStatsTools.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDMostLikelyPartitionsFinder_H
#define LSDMostLikelyPartitionsFinder_H

/// @brief This object looks for the most likeley partitions or segments of 2D data.
class LSDMostLikelyPartitionsFinder
{
  public:
    /// @brief Create LSDMostLikelyPartitionsFinder object.
    /// @param this_min_seg_length int Minimum segement length required.
    /// @param this_x_data vector<float> X data.
    /// @param this_y_data vector<float> Y data.
    /// @author SMM
      /// @date 01/03/13
    LSDMostLikelyPartitionsFinder(int this_min_seg_length, vector<float> this_x_data, vector<float> this_y_data)
          { create( this_min_seg_length, this_x_data, this_y_data); }

    // some get functions
    /// @return number of nodes.
    int get_n_nodes()                    {return int(x_data.size()); }
    /// @return Maximum Likelihood Estimator of segments.
    vector<float> get_MLE_of_segments()  {return MLE_of_segments;}
    /// @return Vector of X data.
    vector<float> get_x_data()             {return x_data; }
    /// @return Vector of Y data.
    vector<float> get_y_data()             {return y_data; }

    // functions for thinning the data

    /// @brief Resets all the derived data members.
    /// @details Derived data members are AIC, AICc, fit statistics, liklihoods (in the liklihood
    ///  matrix, etc.) Basically everything except for the x and y data.
    /// @author SMM
        /// @date 01/03/13
    void reset_derived_data_members();

    /// @brief Collects data as close as possible to some target dx, but does not modify invidivual data points.
    /// @param dx Target dx for thinning.
    /// @author SMM
        /// @date 01/03/13
    void thin_data_target_dx_preserve_data(float dx);

    /// @brief Collects data as close as possible to some target dx, but does not modify invidivual data points.
    /// @param dx Target dx for thinning.
    /// @param node_ref An index vector of the data points that were selected.
    /// @author SMM
        /// @date 01/03/13
    void thin_data_target_dx_preserve_data(float dx, vector<int>& node_ref);

    /// @brief Collects data as close as possible to some target dx, but does not modify invidivual data points.
    /// @param dx Target dx for thinning.
    /// @return New thinned LSDMostLikelyPartitionsFinder object.
    /// @author SMM
        /// @date 01/03/13
    LSDMostLikelyPartitionsFinder spawn_thinned_data_target_dx_preserve_data(float dx);

    /// @brief This thins the data by skipping elements.
    ///
    /// @details A positive number N means it skips N elements after each element.\n\n
    /// A negative number -N means that after N elements it skips one element.
    /// @param N skip value.
    /// @param node_ref An index vector of the data points that were selected.
    /// @author SMM
        /// @date 01/03/13
    void thin_data_skip(int N, vector<int>& node_ref);

    /// @brief Recasts data at fixed values of dx, using linear interpoloation.
    /// @param dx Target dx for thinning.
    /// @author SMM
        /// @date 01/03/13
    void thin_data_target_dx_linear_interpolation(float dx);
    /// @brief Recasts data at fixed values of dx, using linear interpoloation.
    /// @param dx Target dx for thinning.
    /// @return New thinned LSDMostLikelyPartitionsFinder object.
    /// @author SMM
        /// @date 01/03/13
    LSDMostLikelyPartitionsFinder spawn_thinned_data_target_dx_linear_interpolation(float dx);

    /// @brief Skips nodes but using a Monte Carlo scheme that samples random points along the channel profile.
    /// @param Mean_skip
    /// @param skip_range
    /// @param node_ref An index vector of the data points that were selected.
    /// @author SMM
        /// @date 01/05/13
    void thin_data_monte_carlo_skip(int Mean_skip,int skip_range, vector<int>& node_ref);

    /// @brief Thins object based on a monte carlo approach using a mean, max and minimum dchi.
    /// @param mean_dchi
    /// @param variation_dchi
    /// @param node_ref An index vector of the data points that were selected.
    /// @author SMM
    /// @date 01/03/13
    void thin_data_monte_carlo_dchi(float mean_dchi, float variation_dchi, vector<int>& node_ref);

    /// @brief Function for looking at the x and y data.
    /// @author SMM
    /// @date 01/03/13
    void print_x_y_data_to_screen();

    // this function drives the whole shebang

    /// @brief Driver function to get best fit segments.
    /// @param sigma_values vector<float> a vector containing sigma values for each node
    /// @author SMM
    /// @date 01/03/13
    void best_fit_driver_AIC_for_linear_segments(vector<float> sigma_values);
    /// @brief Driver function to get best fit segments.
    /// @param sigma_values
    /// @author SMM
    /// @date 01/03/13
    void best_fit_driver_AIC_for_linear_segments(float sigma_values);

    /// @brief Function returns data for a given sigma value.
    ///
    /// @details this_MLE, this_n_segments and this_n_nodes are all returned
    /// so the user can combine two or more segments and get an AIC or AICc.
    /// @param node index into the sigma vector. This requires a bit of explantion:
    ///   because sigma can be extracted from the product that calculates MLE as a constant
    ///   we can have multiple values of sigma and recalculate MLE for any given sigma
    ///   When this code was written we added a function to modify the sigma values. 
    ///   In this case there are no normalisation values applied so if you use 
    ///   a random sigma here you will get the base sigma (supplied to the main
    ///   computation function)
    /// @param sigma_values
    /// @param b_values A vector containing the b value (intercept) of each segment
    /// @param m_values A vector containing the m value (slope) of each segment
    /// @param r2_values A vector containing the r^2 value of each segment
    /// @param DW_values A vector containing the durbin-watson value of each segment
    /// @param fitted_y A vector containing the y value of each segment
    /// @param seg_lengths A vector with the lengths of each segment
    /// @param this_MLE The MLE of the best fit. 
    /// @param this_n_segments The number of segments
    /// @param this_n_nodes The number of total nodes
    /// @param this_AIC The Aikake information criterion
    /// @param this_AICc The corrected Aikake information criterion (needed for finite sample size)
    /// @author SMM
    /// @date 01/03/13
    void get_data_from_best_fit_lines(int node, vector<float> sigma_values,
                  vector<float>& b_values, vector<float>& m_values,
                  vector<float>& r2_values, vector<float>& DW_values,
            vector<float>& fitted_y,vector<int>& seg_lengths,
            float& this_MLE, int& this_n_segments, int& this_n_nodes,
                                                float& this_AIC, float& this_AICc);

    /// @brief Replaces two vectors, which have the starting and ending position of the best fit segments for a given sigma.
    ///
    /// @details This gets data from (most likeley) the get_data_from_best_fit_lines function.
    /// @param start_x
    /// @param end_x
    /// @param seg_lengths
    /// @author SMM
    /// @date 01/03/13
    void get_start_and_end_x_for_segments(vector<float>& start_x, vector<float>& end_x, vector<int> seg_lengths);

    // these functions populate the arrays used for calcualting the best fit segments

    /// @brief This function is used to calculate the slope, intercept, and likelihood of all possible linear segments along a series of data points.
    ///
    /// @details The function requires the data in x and y vectors, the maximum segment length
    /// and sigma, the standard deviation of the measured data. This will be approxamately
    /// the error in the surface elevation, although it might have to be increased simply because
    /// likelihood will tend to zero if this is too small. sigma should also be considered to
    /// contain the 'noise' inherent in channel incision so perhaps 1-5 metres is appropriate
    /// the maximum segment length is an integer: it is the number of data points used.
    /// these data points from raw chi data are irregularly spaced so two segments of the same
    /// 'length' can have different lengths in chi space. One remedey for this is a preprocessor that
    /// places the zeta vs chi data along evenly spaced points.
    ///
    /// The routine generates three matrices. The row of the matrix is the starting node of the segment.
    /// The column of the matrix is the ending node of the segment. Thus the routine will generate a
    /// matrix that is dimension n x n where n is the number of data points.
    ///
    /// One potential future development is to implement this using a sparse matrix from the boost mtl
    /// library to reduce the memory usage.
    /// @param sigma Standard deviation of error.
    /// @author SMM
    /// @date 01/03/13
    void calculate_segment_matrices(float sigma);

    /// @brief This function popultes the matrices of liklihood, m and b values.
    ///
    /// @details It is a recursive algorithm so in fact it doesn't just get one row
    /// but drills down through all the possible starting nodes to complete the matrix.
    /// @param start_node
    /// @param end_node
    /// @param no_data_value No data value
    /// @param sigma Standard deviation of error.
    /// @author SMM
    /// @date 01/03/13
    void populate_segment_matrix(int start_node, int end_node, float no_data_value,float sigma);

    /// @brief Function calculates the most likeley combination of segments given the liklihood of individual segments calcualted by the calculate_segment_matrices function.
  /// @author SMM
    /// @date 01/03/13
    void find_max_like_of_segments();

    /// @brief This function drives the partitioning algorithms.
    /// @param k Number of elements in the partition.
        /// @author SMM
    /// @date 01/03/13
    void partition_driver_to_vecvecvec(int k);

    // functions that perform components of the partioning

    /// @brief An integer partition algorithm.
    ///
    /// @details Algorithm and original Pascal implementation: Frank Ruskey, 1995. \n
    /// Translation to C: Joe Sawada, 1997 \n
    /// grabbed from http://theory.cs.uvic.ca/inf/nump/NumPartition.html  \n
    /// adapted smm 21/12/2012   \n
    /// algorith described in   \n
    /// http://mathworld.wolfram.com/PartitionFunctionP.html  \n
    /// and        \n
    /// Skiena, S. Implementing Discrete Mathematics: Combinatorics and Graph Theory with Mathematica. Reading, MA: Addison-Wesley, 1990. \n
    ///
    /// This is a further adaptation that only presents solution to the partition with segments of a minimum length it stores all the partitions.\n
    /// http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5123/Permutations-in-C.htm    \n
    /// http://www.cplusplus.com/reference/algorithm/next_permutation/   \n
    /// http://mdm4u1.wetpaint.com/page/4.3+Permutations+with+Some+Identical+Elements  \n
    /// @param n
    /// @param k
    /// @param t
    /// @param p
    /// @author SMM
    /// @date 01/03/13
    void partitions_with_minimum_length(int n, int k, int t, vector<int>& p);
    /// @brief Function for use with the permutations this assigns values into the vecvecvec that contains all the partitioning information.
    /// @param t
    /// @param p
        /// @author SMM
    /// @date 01/03/13
    void partition_assign(int t, vector<int>& p);
    /// @brief Function for use with the permutations gets the mininum of two values.
    /// @param x
    /// @param y
    /// @return Minimum of the two values.
        /// @author SMM
    /// @date 01/03/13
    int LSDpartitions_min( int x, int y);

    // functions for working with liklihoods. Can transform likelihoods between different sigma values

    /// @brief This takes a likelihood array that has been calcualted with a given sigma value and
    /// normalizes the sigma values as though sigma was equal to 1.
    /// @param sigma Standard deviation of error.
    /// @return Normalized sigma matrix.
     /// @author SMM
    /// @date 01/03/13
    Array2D<float> normalize_like_matrix_to_sigma_one(float sigma);

    /// @brief Normalizes but with vector data, for use with MLE vector for segments.
    /// @param sigma Standard deviation of error.
    /// @param like_vector Likelihood vector.
    /// @return Normalized sigma vector.
    /// @author SMM
    /// @date 01/03/13
    vector<float> normalize_like_vector_to_sigma_one(float sigma, vector<float> like_vector);

    /// @brief Takes a normalized likelihood array and updates the values to a new sigma value.
    /// @param sigma Standard deviation of error.
    /// @param sig1_like_array
    /// @author SMM
    /// @date 01/03/13
    void change_normalized_like_matrix_to_new_sigma(float sigma, Array2D<float>& sig1_like_array);

    /// @brief Takes a normalized likelihood vector and updates the values to a new sigma value.
    /// @param sigma Standard deviation of error.
    /// @param sig1_like_vector
    /// @return Updated sigma vector.
    /// @author SMM
    /// @date 01/03/13
    vector<float> change_normalized_like_vector_to_new_sigma(float sigma, vector<float> sig1_like_vector);

    /// @brief Takes a normalized likelihood vector and updates the values to a new sigma value.
    /// @param sigma1 Standard deviation of error.
    /// @param sig1_like_vector
    /// @param sigma2 Standard deviation of error.
    /// @return Updated sigma vector.
    /// @author SMM
    /// @date 01/03/13
    vector<float> transform_like_from_sigma1_to_sigma2(float sigma1,
                          vector<float> sig1_like_vector, float sigma2);

    // these get the best fit number of segments for a variety of sigma values.
    // because of the way the AIC and AICc algorithms work, the minimum AIC for different number of segments
    // can vary depending on sigma, which we don't know. Therefore, we include this function to scan for best fits across
    // different values fo sigma

    /// @brief Function takes the normalized MLE values (normalized with sigma = 1) and returns the best fit number of segments from both the AIC and the AICc measures.
    ///
    /// @details It also returns two vector of vectors which are the AIC values for the varius values of sigma passed to the function in the sigma values vector.
    /// @param sigma_values vector of sigma values.
    void get_n_segments_for_various_sigma(vector<float> sigma_values);

    /// @brief Function calculates AIC and AICc of segments taking the maximum_MLE based on a sigma of one.
    /// @param sigma Standard deviation of error.
    /// @param AIC_of_segments
    /// @param AICc_of_segments
     /// @author SMM
    /// @date 01/03/13
    void calculate_AIC_of_segments_with_variable_sigma(float sigma,
                    vector<float>& AIC_of_segments,
                    vector<float>& AICc_of_segments);

    /// @brief Function extracts the m, b, r^2 and DW statistic of the most likeley segments.
    /// @param bestfit_segments_node
    /// @param m_values
    /// @param b_values
    /// @param r2_values
    /// @param DW_values
    /// @author SMM
    /// @date 01/03/13
    void get_properties_of_best_fit_segments(int bestfit_segments_node,
                     vector<float>& m_values, vector<float>& b_values,
                     vector<float>& r2_values, vector<float>& DW_values);

    // some functions for printing results to screen
    /// @brief Function prints the most likeley segment lengths to screen.
    /// @author SMM
      /// @date 01/03/13
    void print_to_screen_most_likeley_segment_lengths();
    /// @brief Function prints AIC and AICc values to screen.
    /// @author SMM
      /// @date 01/03/13
    void print_AIC_and_AICc_to_screen(vector<float> sigma_values);

  protected:
    ///Minimum segment length.
    int minimum_segment_length;
    /// X data.
    vector<float> x_data;
    /// Y data.
    vector<float> y_data;

    /// The base sigma value from which the MLE of the segments is calcluated.
    float base_sigma;

    // arrays containing the properties of the segments
    // the arrays are indexed so the row is the starting node and the column is the ending node

    /// Liklihood array. Indexed so the row is the starting node and the column is the ending node.
    Array2D<float> like_array;
    /// Slope array. Indexed so the row is the starting node and the column is the ending node.
    Array2D<float> m_array;
    /// Intercept array. Indexed so the row is the starting node and the column is the ending node.
    Array2D<float> b_array;          //
    /// R^2 array. Indexed so the row is the starting node and the column is the ending node.
    Array2D<float> rsquared_array;
    /// @brief Array of Durbin-Watson statistics to test if the residuals are autocorrelated.
    ///
    /// @details Used to determine if the segment is truly linear. Values less than 1 indicate that the segment is probably not linear
    /// values < 1.5 should arouse suspicion. Indexed so the row is the starting node and the column is the ending node.
    Array2D<float> DW_array;

    /// Maximum likelihood of the different number of segments.
    vector<float> MLE_of_segments;

    /// @brief Each element of this vector contains the most likeley segments for that number of segments.
    ///
    /// @details So for example segments_for_each_n_segments[3] is a vector containing the lengths of the most likeley
    /// segments for 4 segments (note 0 indexing).
    vector< vector<int> > segments_for_each_n_segments;

    /// @brief This is vecvecvec.
    ///
    /// @details Top index references the number of segments. Second layer
    /// loops through the possible partitions of that number of segments. Third layer
    /// is the individual segment lengths.
    vector< vector < vector<int> > > partitions;

    /// Vector of best fit AIC values.
    vector<int> best_fit_AIC;
    /// Vector of best fit AICc values.
    vector<int> best_fit_AICc;
    /// Vector of vectors of AIC values.
    vector< vector<float> > AIC_for_each_n_segments;
    /// Vector of vectors of AICc values.
    vector< vector<float> > AICc_for_each_n_segments;

  private:
    void create(int this_min_seg_length, vector<float> this_x_data, vector<float> this_y_data);
};

#endif
