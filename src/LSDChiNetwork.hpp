//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiNetwork
// Land Surface Dynamics ChiNetwork
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for analysing channels using the integral method of channel
//  analysis
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

/** @file LSDChiNetwork.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 1.0.0
@brief This object is used to examine a network of channels in chi space.

@date 28/02/2013
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDChiNetwork_H
#define LSDChiNetwork_H

/// @brief This object is used to examine a network of channels in chi space.
class LSDChiNetwork
{
  public:
    /// @brief Crate routine to make a LSDChiNetwork object.
    /// @param channel_network_fname Filename.
    LSDChiNetwork(string channel_network_fname)
                 { create( channel_network_fname ); }

    LSDChiNetwork(){create();} // empty constructor
    
    LSDChiNetwork(LSDFlowInfo& FlowInfo, int SourceNode, int OutletNode, LSDRaster& Elevation,
                           LSDRaster& FlowDistance, LSDRaster& DrainageArea)
                 {create(FlowInfo, SourceNode, OutletNode, Elevation,
                         FlowDistance, DrainageArea); }
                         
    LSDChiNetwork(LSDFlowInfo& FlowInfo, int SourceNode, int OutletNode, LSDRaster& Elevation,
                           LSDRaster& FlowDistance, LSDRaster& DrainageArea, LSDRaster& Chi)
                 {create(FlowInfo, SourceNode, OutletNode, Elevation,
                         FlowDistance, DrainageArea, Chi); }


    /// @return Number of channels.
    int get_n_channels()  { return int(node_indices.size()); }

    // get functions. These are used for interfacing with
    // the LSDRaster object (not in the standalone version)
    /// @return Number of rows as an integer.
    int get_NRows() const        { return NRows; }
    /// @return Number of columns as an integer.
    int get_NCols() const        { return NCols; }
    /// @return Minimum X coordinate as an integer.
    float get_XMinimum() const        { return XMinimum; }
    /// @return Minimum Y coordinate as an integer.
    float get_YMinimum() const        { return YMinimum; }
    /// @return Data resolution as an integer.
    float get_DataResolution() const   { return DataResolution; }
    /// @return No Data Value as an integer.
    int get_NoDataValue() const        { return NoDataValue; }

    // printing routines for bug checking
    /// @brief Print channel details to screen for bug checking.
    ///
    /// @details Format of file: \n\n channel_number << " " << receiver_channel[channel_number] << " "
    ///             << node_on_receiver_channel[channel_number] << " "
    ///             << node[i] << " " << row[i] << " " << col[i] << " " << flow_distance[i] << " "
    ///         << chi[i] << " " << elevation[i] << " " << drainage_area[i]
    /// @param channel_number Channel to examine
    /// @author SMM
    /// @date 01/04/13
    void print_channel_details_to_screen(int channel_number);

    /// @brief Print channel details to file for bug checking.
    ///
    /// @details Format of file: \n\n channel_number << " " << receiver_channel[channel_number] << " "
    ///             << node_on_receiver_channel[channel_number] << " "
    ///             << node[i] << " " << row[i] << " " << col[i] << " " << flow_distance[i] << " "
    ///         << chi[i] << " " << elevation[i] << " " << drainage_area[i]
    /// @param fname Output filename.
    /// @param A_0 A_0 value.
    /// @param m_over_n  m over n ratio.
    /// @author SMM
    /// @date 01/04/13
    void print_channel_details_to_file(string fname, float A_0, float m_over_n);

    /// @brief This function prints the details of all channels to a file.
    ///
    /// @details It includes data from monte carlo fitting. Format is: \n\n
    /// A_0 m_over_n channel_number node_on_receiver_channel node_index row col flow_distance chi elevation darainage_area.
    /// @param fname Output filename.
    /// @author SMM
    /// @date 01/04/13
    void print_channel_details_to_file_full_fitted(string fname);

    /// @brief This function prints the details of all channels to a csv file that
    /// can be ingested by ArcMap.
    /// @details It includes data from monte carlo fitting. 
    /// @param fname Output filename, which has _for_Arc automatically appended.
    /// @author SMM
    /// @date 28/02/14
    void print_channel_details_to_file_full_fitted_for_ArcMap(string fname);

    /// @brief This function prints the details of all channels to a file.
    ///
    /// @details It includes data from monte carlo fitting. Format is: \n\n
    /// A_0 m_over_n channel_number node_on_receiver_channel node_index row col flow_distance chi elevation darainage_area.
    /// @param fname Output filename.
    /// @param target_nodes
    /// @param minimum_segment_length
    /// @author SMM
    /// @date 01/04/13
    void print_channel_details_to_file_full_fitted(string fname, int target_nodes,
                                                   int minimum_segment_length);

    /// @brief This extends the tributary channels all the way to the outlet.
    ///
    /// @details In its current version this only works if the tributaries all drain
    /// to the mainstem.
    /// @author SMM
    /// @date 01/04/13
    void extend_tributaries_to_outlet();

    /// @brief Routine for returning calculated data to an array.
    ///
    /// @details It includes a switch that tells the function what data member to write to the array code for data members: \n\n
    /// 1 elevations \n
    /// 2  chis        \n
    /// 3  chi_m_means   \n
    /// 4  chi_m_standard_deviations\n
    /// 5  chi_m_standard_errors\n
    /// 6  chi_b_means            \n
    /// 7  chi_b_standard_deviations\n
    /// 8  chi_b_standard_errors      \n
    /// 9  chi_DW_means                 \n
    /// 10 chi_DW_standard_deviations     \n
    /// 11 chi_DW_standard_errors           \n
    /// 12 all_fitted_DW_means                \n
    /// 13 all_fitted_DW_standard_deviations    \n
    /// 14 all_fitted_DW_standard_errors          \n
    /// @param data_member Switch to select data to be written.
    /// @return Array of data.
    /// @author SMM
    /// @date 01/04/13
    Array2D<float> data_to_array(int data_member);

    // routines for getting slope-area data
    // these print to file at the moment.

    /// @brief Extract slope over fixed vertical intervals.
    ///
    /// This one is the vertical intervals version: it measures slope over fixed vertical
    /// intervals as reccomended by Wobus et al 2006.
    /// This function gets slope data and area data for use in making slope area plots.
    ///  It generates several data elements, which are written to the file with name fname (passed to
    ///  function). The file format is for each row:\n\n
    ///   chan << " " << start_row << " " << mp_row << " " << end_row << " "
    ///		<< start_col << " " << mp_col << " " << end_col << " "
    ///		<< start_interval_elevations << " "
    ///		<< mp_interval_elevations << " " << end_interval_elevations << " "
    ///		<< start_interval_flowdistance << " " << mp_interval_flowdistance << " "
    ///		<< end_interval_flowdistance << " "
    ///		<< start_area << " " << mp_area << " " << end_area << " " << slope
    ///		<< " " << log10(mp_area) << " " << log10(slope)\n\n
    ///
    ///  where start, mp and end denote the start of the interval over which slope is measured, the midpoint
    ///    and the end.
    ///
    /// The area thin fraction is used to thin the data so that segments with large changes in drainage
    /// area are not used in the regression (because these will affect the mean slope)
    /// the fraction is determined by (downslope_area-upslope_area)/midpoint_area.
    /// So if the fraction is 1 it means that the change is area is equal to the area at the midpoint
    /// a restictive value is 0.05, you will eliminate major tributaries with a 0.2, and
    /// 1 will catch almost all of the data.
    /// @param interval
    /// @param area_thin_fraction
    /// @param fname Output filename
    /// @author SMM
    /// @date 01/04/13
    void slope_area_extraction_vertical_intervals(float interval, float area_thin_fraction,
                                                  string fname);

    /// @brief Extract slope over fixed horizontal intervals.
    ///
    /// This one is the horizontal intervals version: it measures slope over fixed flow distance
    /// as used by many authors including DiBiasie et al 2010 and Ouimet et al 2009
    /// This function gets slope data and area data for use in making slope area plots.
    ///  It generates several data elements, which are written to the file with name fname (passed to
    ///  function). The file format is for each row:   \n\n
    ///   chan << " " << start_row << " " << mp_row << " " << end_row << " "
    ///		<< start_col << " " << mp_col << " " << end_col << " "
    ///		<< start_interval_elevations << " "
    ///		<< mp_interval_elevations << " " << end_interval_elevations << " "
    ///		<< start_interval_flowdistance << " " << mp_interval_flowdistance << " "
    ///		<< end_interval_flowdistance << " "
    ///		<< start_area << " " << mp_area << " " << end_area << " " << slope
    ///		<< " " << log10(mp_area) << " " << log10(slope) \n\n
    ///
    ///  where start, mp and end denote the start of the interval over which slope is measured, the midpoint
    ///    and the end.
    ///
    /// The area thin fraction is used to thin the data so that segments with large changes in drainage
    /// area are not used in the regression (because these will affect the mean slope)
    /// the fraction is determined by (downslope_area-upslope_area)/midpoint_area.
    /// So if the fraction is 1 it means that the change is area is equal to the area at the midpoint
    /// a restictive value is 0.05, you will eliminate major tributaries with a 0.2, and
    /// 1 will catch almost all of the data.
    /// @param interval
    /// @param area_thin_fraction
    /// @param fname Output filename
    /// @author SMM
    /// @date 01/04/13
    void slope_area_extraction_horizontal_intervals(float interval, float area_thin_fraction,
                                   string fname);

    // routines for calculating chi and maniplating chi

    /// @brief This function calculates the chi values for the channel network using the rectangle rule.
    ///
    /// @details Note: the entire network must be caluculated because the chi values of the tributaries
    /// depend on the chi values of the mainstem.
    /// @param A_0 A_0 value.
    /// @param m_over_n  m over n ratio.
    /// @author SMM
    /// @date 01/04/13
    void calculate_chi(float A_0, float m_over_n);

    /// @brief This function calucaltes the chi spacing of the main stem channel (the longest channel).
    ///
    /// @details The maximum length of the dataset will be in the main stem so this will determine the
    /// target spacing of all the tributaries.
    /// @param target_nodes Node index of the target node.
    /// @return Optimal chi spacing.
    /// @author SMM
    /// @date 01/04/13
    float calculate_optimal_chi_spacing(int target_nodes);

    /// @brief This function calucaltes the skip parameter of the main stem (the longest channel).
    ///
    /// @details The maximum length of the dataset will be in the main stem so this will determine the target spacing of all the tributaries.
    /// @param target_nodes Node index of the target node.
    /// @return Skip value.
    /// @author SMM
    /// @date 01/06/13
    int calculate_skip(int target_nodes);

    /// @brief This function calucaltes the skip parameter based on a vector of chi values.
    /// @param target_nodes Node index of the target node.
    /// @param sorted_chis Vector of chi values
    /// @return Skip value.
    /// @author SMM
    /// @date 01/06/13
    int calculate_skip(int target_nodes, vector<float>& sorted_chis);

    /// @brief This function calucaltes the skip parameter of a give channel.
    ///
    /// @details The maximum length of the dataset will be in the main stem so this will determine the target spacing of all the tributaries.
    /// @param target_nodes Node index of the target node.
    /// @param channel_number The channel to be analysed.
    /// @return Skip value.
    /// @author SMM
    /// @date 01/04/13
    int calculate_skip(int target_nodes, int channel_number);


    // routines for calcualting the most likeley segments.

    /// @brief This function gets the most likely channel segments for a particular channel.
    ///
    /// @details This function replaces the b, m, r2 and DW values of each segment into vectors
    /// it also returns the fitted elevation and the index into the original channel (since this is done
    ///		with thinned data).
    /// @param channel The index into the channel.
    /// @param minimum_segment_length is how many nodes the mimimum segment will have.
    /// @param sigma is the standard deviation of error on elevation data
    /// @param N
    /// @param b_vec
    /// @param m_vec
    /// @param r2_vec
    /// @param DW_vec
    /// @param thinned_chi
    /// @param thinned_elev
    /// @param fitted_elev
    /// @param node_reference
    /// @param these_segment_lengths
    /// @param this_MLE
    /// @param this_n_segments
    /// @param n_data_nodes
    /// @param this_AIC
    /// @param this_AICc
    /// @author SMM
    /// @date 01/04/13
    void find_most_likeley_segments(int channel,int minimum_segment_length,
               float sigma, int N, vector<float>& b_vec,
               vector<float>& m_vec, vector<float>& r2_vec,vector<float>& DW_vec,
               vector<float>& thinned_chi, vector<float>& thinned_elev,
               vector<float>& fitted_elev, vector<int>& node_reference,
               vector<int>& these_segment_lengths,
               float& this_MLE, int& this_n_segments, int& n_data_nodes,
               float& this_AIC, float& this_AICc );

    /// @brief This function gets the most likely channel segments for a particular channel.
    ///
    /// @details This function replaces the b, m, r2 and DW values of each segment into vectors
    ///   it also returns the fitted elevation and the index into the original channel (since this is done
    ///   with thinned data).
    /// @param channel The index into the channel.
    /// @param minimum_segment_length is how many nodes the mimimum segment will have.
    /// @param sigma is the standard deviation of error on elevation data
    /// @param dchi
    /// @param b_vec
    /// @param m_vec
    /// @param r2_vec
    /// @param DW_vec
    /// @param thinned_chi
    /// @param thinned_elev
    /// @param fitted_elev
    /// @param node_reference
    /// @param these_segment_lengths
    /// @param this_MLE
    /// @param this_n_segments
    /// @param n_data_nodes
    /// @param this_AIC
    /// @param this_AICc
    /// @author SMM
    /// @date 01/04/13
    void find_most_likeley_segments_dchi(int channel,int minimum_segment_length,
               float sigma, float dchi, vector<float>& b_vec,
               vector<float>& m_vec, vector<float>& r2_vec,vector<float>& DW_vec,
               vector<float>& thinned_chi, vector<float>& thinned_elev,
               vector<float>& fitted_elev, vector<int>& node_reference,
               vector<int>& these_segment_lengths,
               float& this_MLE, int& this_n_segments, int& n_data_nodes,
               float& this_AIC, float& this_AICc );

    /// @brief This gets the most likely segments but uses the monte carlo data thinning method.
    ///
    /// @details The expectation is that this will be used repeatedly on channels to generate statistics of the
    /// best fit segments by individual nodes in the channel network.
    /// @param channel The index into the channel.
    /// @param minimum_segment_length is how many nodes the mimimum segment will have.
    /// @param sigma is the standard deviation of error on elevation data
    /// @param mean_skip
    /// @param skip_range
    /// @param b_vec
    /// @param m_vec
    /// @param r2_vec
    /// @param DW_vec
    /// @param thinned_chi
    /// @param thinned_elev
    /// @param fitted_elev
    /// @param node_reference
    /// @param these_segment_lengths
    /// @param this_MLE
    /// @param this_n_segments
    /// @param n_data_nodes
    /// @param this_AIC
    /// @param this_AICc
    /// @author SMM
    /// @date 01/04/13
    void find_most_likeley_segments_monte_carlo(int channel, int minimum_segment_length,
               float sigma, int mean_skip, int skip_range, vector<float>& b_vec,
               vector<float>& m_vec, vector<float>& r2_vec,vector<float>& DW_vec,
               vector<float>& thinned_chi, vector<float>& thinned_elev,
               vector<float>& fitted_elev, vector<int>& node_reference,
               vector<int>& these_segment_lengths,
               float& this_MLE, int& this_n_segments, int& n_data_nodes,
               float& this_AIC, float& this_AICc );

    /// @brief This gets the most likely segments but uses the monte carlo data thinning method.
    ///
    /// @details The mean_dchi is the mean value of dchi, and the variation is how much the
    /// chi chan vary, such that minimum_dchi = dchi-variation_dchi.
    /// The expectation is that this will be used repeatedly on channels to generate statistics of the
    /// best fit segments by individual nodes in the channel network.
    /// @param channel The index into the channel.
    /// @param minimum_segment_length is how many nodes the mimimum segment will have.
    /// @param sigma is the standard deviation of error on elevation data
    /// @param mean_dchi
    /// @param variation_dchi
    /// @param b_vec
    /// @param m_vec
    /// @param r2_vec
    /// @param DW_vec
    /// @param thinned_chi
    /// @param thinned_elev
    /// @param fitted_elev
    /// @param node_reference
    /// @param these_segment_lengths
    /// @param this_MLE
    /// @param this_n_segments
    /// @param n_data_nodes
    /// @param this_AIC
    /// @param this_AICc
    /// @author SMM
    /// @date 01/04/13
    void find_most_likeley_segments_monte_carlo_dchi(int channel, int minimum_segment_length,
               float sigma, float mean_dchi, float variation_dchi, vector<float>& b_vec,
               vector<float>& m_vec, vector<float>& r2_vec,vector<float>& DW_vec,
               vector<float>& thinned_chi, vector<float>& thinned_elev,
               vector<float>& fitted_elev, vector<int>& node_reference,
               vector<int>& these_segment_lengths,
               float& this_MLE, int& this_n_segments, int& n_data_nodes,
               float& this_AIC, float& this_AICc );

    /// @brief The master routine for calculating the best fit m over n values for a channel network, based on a fixed value of dchi.
    /// @param A_0
    /// @param n_movern
    /// @param d_movern
    /// @param start_movern
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_nodes_mainstem
    /// @param fname Output filename
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/04/13
    float search_for_best_fit_m_over_n_dchi(float A_0, int n_movern, float d_movern, float start_movern,
                 int minimum_segment_length, float sigma, int target_nodes_mainstem, string fname);

    /// @brief The master routine for calculating the best fit m over n values for a channel network
    /// @param A_0
    /// @param n_movern
    /// @param d_movern
    /// @param start_movern
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_nodes_mainstem
    /// @param fname Output filename
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/04/13
    float search_for_best_fit_m_over_n(float A_0, int n_movern, float d_movern, float start_movern,
                 int minimum_segment_length, float sigma, int target_nodes_mainstem, string fname);

    /// @brief Routine for calculating the best fit m over n values for a channel network, but calculates the mainstem and the tributaries seperately.
    /// @param A_0
    /// @param n_movern
    /// @param d_movern
    /// @param start_movern
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_nodes_mainstem
    /// @param fname Output filename
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/04/13
    float search_for_best_fit_m_over_n_seperate_ms_and_tribs(float A_0, int n_movern, float d_movern, float start_movern,
                 int minimum_segment_length, float sigma, int target_nodes_mainstem, string fname);

    /// @brief This function looks for the best fit values of m over n by simply testing for the least variation in the tributaries.
    /// @param A_0
    /// @param n_movern
    /// @param d_movern
    /// @param start_movern
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_nodes
    /// @param n_iterations
    /// @param m_over_n_values
    /// @param AICc_mean
    /// @param AICc_sdtd
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/04/13
    float search_for_best_fit_m_over_n_colinearity_test(float A_0, int n_movern, float d_movern,
          float start_movern, int minimum_segment_length, float sigma,
          int target_nodes, int n_iterations,
          vector<float>& m_over_n_values,
          vector<float>& AICc_mean, vector<float>& AICc_sdtd);

    /// @brief This function calculeates best fit m/n using the collinearity test \n
    ///  these channels are ones with breaks
    /// @param A_0 float the reference area
    /// @param n_movern int the number of m over n ratios to iterate through
    /// @param d_movern float the change in m/n in each iterations
    /// @param start_movern float the starting value of m/n
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_skip int the mean skipping value
    /// @param target_nodes int the target number of nodes in a break
    /// @param n_iterations	int the number of iterations
    /// @param m_over_n_values	vector<float>& this gets written, it contains the m/n values for the run
    /// @param AICc_mean vector<float> gets written, the mean values of the AICc for each m/n
    /// @param AICc_sdtd vector<float> gets written, the standard deviation values of the AICc for each m/n
    /// @param Monte_Carlo_switch int if 1, run the code using the iterative Monte Carlo scheme
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/07/13
    float search_for_best_fit_m_over_n_colinearity_test_with_breaks(float A_0, int n_movern, float d_movern,
                 float start_movern, int minimum_segment_length, float sigma,
                 int target_skip, int target_nodes, int n_iterations,
                 vector<float>& m_over_n_values, vector<float>& AICc_mean, vector<float>& AICc_sdtd,
                 int Monte_Carlo_switch);

    /// @brief This function calculeates best fit m/n for each channel these channels are ones with breaks
    /// @param A_0 float the reference area
    /// @details this does not report variability of the AICc values and so should not be used, instead
    /// \n use the Monte Carlo version
    /// \n retained in case you want rapid calculation of best fit m/n
    /// @param n_movern int the number of m over n ratios to iterate through
    /// @param d_movern float the change in m/n in each iterations
    /// @param start_movern float the starting value of m/n
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_skip int the mean skipping value
    /// @param target_nodes int the target number of nodes in a break
    /// @param n_iterations	int the number of iterations
    /// @param m_over_n_values	vector<float>& this gets written, it contains the m/n values for the run
    /// @param AICc_vals
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/04/13
    float search_for_best_fit_m_over_n_individual_channels_with_breaks(float A_0, int n_movern, float d_movern,
                    float start_movern, int minimum_segment_length, float sigma,
                    int target_skip, int target_nodes, int n_iterations,
                    vector<float>& m_over_n_values, vector< vector<float> >& AICc_vals);

    /// @brief This gets the best fit m over n values of all the individual tributaries.
    ///
    /// @details It uses a monte carlo appraoach so all tributaries have both the mean and the variability
    /// of the AICc values reported.
    /// @param A_0 float the reference area
    /// @param n_movern int the number of m over n ratios to iterate through
    /// @param d_movern float the change in m/n in each iterations
    /// @param start_movern float the starting value of m/n
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_skip int the mean skipping value
    /// @param target_nodes int the target number of nodes in a break
    /// @param n_iterations	int the number of iterations
    /// @param m_over_n_values	vector<float>& this gets written, it contains the m/n values for the run
    /// @param AICc_means vector<float> gets written, the mean values of the AICc for each m/n
    /// @param AICc_stddev vector<float> gets written, the standard deviation values of the AICc for each m/n
    /// @return Best fit m over n.
    /// @author SMM
    /// @date 01/07/13
    float search_for_best_fit_m_over_n_individual_channels_with_breaks_monte_carlo(float A_0, int n_movern,
                    float d_movern, float start_movern, int minimum_segment_length, float sigma,
                    int target_skip, int target_nodes, int n_iterations,
                    vector<float>& m_over_n_values,
                    vector< vector<float> >& AICc_means, vector< vector<float> >& AICc_stddev);

    /// @brief This routine uses a monte carlo approach to repeatedly sampling all the data in the channel network.
    ///
    /// @details Uses a reduced number of data elements and then populates each channel node with
    /// a distribution of m, b and fitted elevation values. These then can be averaged and details of their variation
    /// calculated. \n
    /// This is based on a fixed value of dchi
    /// @param A_0
    /// @param m_over_n
    /// @param n_iterations
    /// @param fraction_dchi_for_variation
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param target_nodes_mainstem
    /// @author SMM
    /// @date 01/04/13
    void monte_carlo_sample_river_network_for_best_fit_dchi(float A_0, float m_over_n, int n_iterations,
                              float fraction_dchi_for_variation,
                              int minimum_segment_length, float sigma,
                              int target_nodes_mainstem);


    /// @brief Monte carlo segment fitter.
    ///
    /// @details This takes a fixed m_over_n value and then samples the indivudal nodes in the full channel profile
    /// to repeadetly get the best fit segments on thinned data. the m, b fitted elevation, r^2 and DW statistic are all
    /// stored for every iteration on every node. These can then be queried later for mean, standard deviation and
    /// standard error information  \n\n
    ///
    /// the fraction_dchi_for_variation is the fration of the optimal dchi that dchi can vary over. So for example
    /// if this = 0.4 then the variation of dchi will be 0.4*mean_dchi and the minimum dchi will be
    /// min_dchi = (1-0.4)*mean_dchi \n\n
    ///
    /// Note: this is _extremely_ computationally and data intensive.
    /// @param A_0
    /// @param m_over_n
    /// @param n_iterations
    /// @param mean_skip
    /// @param skip_range
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @author SMM
    /// @date 01/04/13
    void monte_carlo_sample_river_network_for_best_fit(float A_0, float m_over_n, int n_iterations,
                int mean_skip, int skip_range,
                int minimum_segment_length, float sigma);


    /// @brief This function samples the river network using monte carlo samplig but after breaking the channels.
    /// @param A_0
    /// @param m_over_n
    /// @param n_iterations
    /// @param skip
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @author SMM
    /// @date 01/04/13
    void monte_carlo_sample_river_network_for_best_fit_after_breaks(float A_0, float m_over_n, int n_iterations,
        int skip, int minimum_segment_length, float sigma);

    /// @brief Monte carlo segment fitter.
    ///
    /// @details This takes a fixed m_over_n value and then samples the indivudal nodes in the full channel profile
    /// to repeadetly get the best fit segments on thinned data. the m, b fitted elevation, r^2 and DW statistic are all
    /// stored for every iteration on every node. These can then be queried later for mean, standard deviation and
    /// standard error information. \n\n
    ///
    /// The break nodes vector tells the algorithm where the breaks in the channel occur
    /// this function is called repeatedly until the target skip equals the all of the this_skip values.
    ///  \n\n
    /// This function continues to split the channel into segments until the target skip is achieved.
    /// @param A_0
    /// @param m_over_n
    /// @param n_iterations
    /// @param target_skip
    /// @param target_nodes
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param chan
    /// @param break_nodes
    /// @author SMM
    /// @date 01/04/13
    void monte_carlo_split_channel(float A_0, float m_over_n, int n_iterations,
        int target_skip, int target_nodes,
        int minimum_segment_length, float sigma, int chan, vector<int>& break_nodes);

    /// @brief This function uses a monte carlo sampling approach to try and split channels.
    ///
    /// @details The channel is sampled at the target skipping interval. It does it with a colinear dataset.
    /// @param A_0
    /// @param m_over_n
    /// @param n_iterations
    /// @param target_skip
    /// @param target_nodes
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param reverse_Chi
    /// @param break_nodes
    /// @param reverse_Elevation
    /// @param break_nodes
    /// @author SMM
    /// @date 01/04/13
    void monte_carlo_split_channel_colinear(float A_0, float m_over_n, int n_iterations,
        int target_skip, int target_nodes,
        int minimum_segment_length, float sigma,
        vector<float> reverse_Chi, vector<float> reverse_Elevation, vector<int>& break_nodes);

    /// @brief This function splits all the channels in one go.
    /// @param A_0
    /// @param m_over_n
    /// @param n_iterations
    /// @param target_skip
    /// @param target_nodes
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @author SMM
    /// @date 01/06/13
    void split_all_channels(float A_0, float m_over_n, int n_iterations,
        int target_skip, int target_nodes, int minimum_segment_length, float sigma);

    /// @brief This function gets the AICc after breaking the channel.
    /// @param A_0
    /// @param m_over_n
    /// @param skip
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param chan
    /// @param break_nodes
    /// @param n_total_segments
    /// @param n_total_nodes
    /// @param cumulative_MLE
    /// @return AICc value
    /// @author SMM
    /// @date 01/06/13
    float calculate_AICc_after_breaks(float A_0, float m_over_n,
        int skip, int minimum_segment_length, float sigma, int chan, vector<int> break_nodes,
        int& n_total_segments, int& n_total_nodes, float& cumulative_MLE);

    /// @brief This function gets the AICc after breaking the channel.
    ///
    /// @details It does this for n_iterations and returns a vector with all of
    /// the AICc values for each iteration reported. This can then be used
    /// to calculate the statistics of the AICc to tell if the minimum
    /// AICc is significantly different from the other AICc values for different
    /// values of m/n
    /// @param A_0
    /// @param m_over_n
    /// @param target_skip
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param chan
    /// @param break_nodes
    /// @param n_total_segments
    /// @param n_total_nodes
    /// @param cumulative_MLE
    /// @param n_iterations
    /// @return AICc value
    /// @author SMM
    /// @date 01/06/13
    vector<float> calculate_AICc_after_breaks_monte_carlo(float A_0, float m_over_n,
        int target_skip, int minimum_segment_length, float sigma, int chan, vector<int> break_nodes,
        int& n_total_segments, int& n_total_nodes, float& cumulative_MLE,
        int n_iterations);

    /// @brief This function gets the AICc after breaking the channelwith a colinear dataset.
    ///
    /// @details The reverse_chi and reverse_elevation data has to be provided.
    /// @param A_0
    /// @param m_over_n
    /// @param skip
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param reverse_chi
    /// @param reverse_elevation
    /// @param break_nodes
    /// @param n_total_segments
    /// @param n_total_nodes
    /// @param cumulative_MLE
    /// @return AICc value
    /// @author SMM
    /// @date 01/06/13
    float calculate_AICc_after_breaks_colinear(float A_0, float m_over_n,
            int skip, int minimum_segment_length, float sigma,
            vector<float> reverse_chi, vector<float> reverse_elevation,
            vector<int> break_nodes,
            int& n_total_segments, int& n_total_nodes, float& cumulative_MLE);

    /// @brief This function gets the AICc after breaking the channelwith a colinear dataset.
    ///
    /// @details The reverse_chi and reverse_elevation data has to be provided. It uses
    /// a monte carlo scheme and returns a vector with all of the AICc values calcluated from the analyses.
    /// @param A_0
    /// @param m_over_n
    /// @param skip
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @param reverse_Chi
    /// @param reverse_Elevation
    /// @param break_nodes
    /// @param n_total_segments
    /// @param n_total_nodes
    /// @param cumulative_MLE
    /// @param n_iterations
    /// @return AICc value
    /// @author SMM
    /// @date 01/06/13
    vector<float> calculate_AICc_after_breaks_colinear_monte_carlo(float A_0, float m_over_n,
        int skip, int minimum_segment_length, float sigma,
        vector<float> reverse_Chi, vector<float> reverse_Elevation,
        vector<int> break_nodes,
        int& n_total_segments, int& n_total_nodes, float& cumulative_MLE,
        int n_iterations);

    /// @brief This routine tests to see if channels are long enough to get a decent fitting from the segment finding algorithms.
    ///
    /// @details Writes to the is_tributary_long_enough vector. If this equals 1, the channel is long enough. If it is zero the channel is not long enough.
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param N
    void is_channel_long_enough_test(int minimum_segment_length,int N);

    /// @brief this function fits 2 segments to the chi-elevation data from first order basins and calculates
    /// the most likely position of the channel head - added by FC 28/06/13
    ///
    /// @param min_seg_length_for_channel_heads Minimum number of nodes used for segment fitting
    /// @return array with channel head locations
    /// @author Fiona Clubb
    /// @date 03/09/2013
    Array2D<float> calculate_channel_heads(int min_seg_length_for_channel_heads);
    
    /// @brief This gets the m_means for the channel network
    /// @return vector of vectors with m means
    /// @ author FJC
    /// @date 04/08/14
    vector< vector<float> > get_m_means();

    /// @brief This gets the m_means for the channel network
    /// @return vector of vectors with m means
    /// @ author FJC
    /// @date 04/08/14
    vector< vector<float> > get_m_standard_deviations()  { return chi_m_standard_deviations; }

    /// @brief This gets the b_means for the channel network
    /// @return vector of vectors with b means
    /// @ author SMM
    /// @date 24/05/16
    vector< vector<float> > get_b_means()  { return chi_b_means; }

    /// @brief This gets the b_standard deviations for the channel network
    /// @return vector of vectors with b tandard deviations
    /// @ author SMM
    /// @date 24/05/16
    vector< vector<float> > get_b_standard_deviations()  { return chi_b_standard_deviations; }

    /// @brief This gets the node_indices for the channel network
    /// @return vector of vectors with m means
    /// @ author DTM
    /// @date 24/03/15
    vector< vector<int> > get_node_indices()
      { return node_indices; }
    /// @brief This gets the chi coordinates for the channel network
    /// @return vector of vectors with m means
    /// @ author DTM
    /// @date 24/03/15
    vector< vector<float> > get_chis()
      { return chis; }


  protected:

    ///Number of rows.
    int NRows;
    ///Number of columns.
    int NCols;
    ///Minimum X coordinate.
    float XMinimum;
    ///Minimum Y coordinate.
    float YMinimum;

    ///Data resolution.
    float DataResolution;
    ///No data value.
    int NoDataValue;
    
    /// This boolean lets the routine know if it is to calculate chi
    bool I_should_calculate_chi; 

    /// Node indices: used in conjunction with other LSD topographic tool objects and not necessary for standalone program.
    vector< vector<int> > node_indices;
    /// Row indices: used in conjunction with other LSD topographic tool objects and not necessary for standalone program.
    vector< vector<int> > row_indices;
    /// Column indices: used in conjunction with other LSD topographic tool objects and not necessary for standalone program.
    vector< vector<int> > col_indices;
    /// The elevations along the channels
    vector< vector<float> > elevations;
    /// Flow distances along channels. Used to integrate to arrive at chi.
    vector< vector<float> > flow_distances;
    /// Drainage areas
    vector< vector<float> > drainage_areas;
    /// The chi values for the channels. This data will be overwritten as m_over_n changes.
    vector< vector<float> > chis;
    /// This is the node on the reciever channel where the tributary enters the channel. Used to find the downstream chi value of a channel.
    vector<int> node_on_receiver_channel;
    /// This is the channel that the tributary enters.
    vector<int> receiver_channel;

    // the following data elements are fitted ci slopes and intercepts, as well
    // as best fit elevation that are generated after monte-carlo sampling
    // they are only filled with data after calling the
    // monte_carlo_sample_river_network_for_best_fit function

    /// This stores the m over n value use to generate the means and standard deviations of the network properties.
    float m_over_n_for_fitted_data;
    /// This stored the A_0 value.
    float A_0_for_fitted_data;
    /// This vector is the same size as the number of channels and is 1 if the channel analysis n_nodes > 3* minimum_segment_length.
    vector<int> is_tributary_long_enough;
    /// Vector of chi_m means.
    vector< vector<float> > chi_m_means;
    /// Vector of chi_m standard deviations.
    vector< vector<float> > chi_m_standard_deviations;
    /// Vector of chi_m standard errors.
    vector< vector<float> > chi_m_standard_errors;
    /// Vector of chi_b means.
    vector< vector<float> > chi_b_means;
    /// Vector of chi_b standard deviations.
    vector< vector<float> > chi_b_standard_deviations;
    /// Vector of chi_b standard errors.
    vector< vector<float> > chi_b_standard_errors;
    /// Vector of fitted elevation means.
    vector< vector<float> > all_fitted_elev_means;
    /// Vector of fitted elevation standard deviations.
    vector< vector<float> > all_fitted_elev_standard_deviations;
    /// Vector of fitted elevation standard errors.
    vector< vector<float> > all_fitted_elev_standard_errors;
    /// Vector of Durbin-Watson means.
    vector< vector<float> > chi_DW_means;
    /// Vector of Durbin-Watson standard deviations.
    vector< vector<float> > chi_DW_standard_deviations;
    /// Vector of Durbin-Watson standard errors.
    vector< vector<float> > chi_DW_standard_errors;
    /// The parameters are generated using a monte carlo approach and not all nodes will have the same number of data points, so the number of data points is stored.
    vector< vector<int> > n_data_points_used_in_stats;
    /// This vector holds the vectors containing the node locations of breaks in the segments.
    vector< vector<int> > break_nodes_vecvec;

  private:
    void create();
    void create(string channel_network_fname);
    void create(LSDFlowInfo& FlowInfo, int SourceNode, int OutletNode, LSDRaster& Elevation,
                           LSDRaster& FlowDistance, LSDRaster& DrainageArea);
    void create(LSDFlowInfo& FlowInfo, int SourceNode, int OutletNode, LSDRaster& Elevation,
                           LSDRaster& FlowDistance, LSDRaster& DrainageArea, 
                           LSDRaster& Chi);
};

#endif
