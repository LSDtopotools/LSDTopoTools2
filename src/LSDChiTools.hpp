//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiTools
// Land Surface Dynamics ChiTools object
//
// The header of the ChiTools object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for performing various analyses in chi space
//
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2016 Simon M. Mudd 2016
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiTools.cpp
// LSDChiTools object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDChiTools_HPP
#define LSDChiTools_HPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
using namespace std;
using namespace TNT;


/// @brief This object packages a number of tools for chi analysis
class LSDChiTools
{
  public:

    /// @brief Create an LSDChiTools from a raster.
    /// @param ThisRaster An LSDRaster object
    /// @author SMM
    /// @date 24/05/2016
    LSDChiTools(LSDRaster& ThisRaster)  { create(ThisRaster); }

    /// @brief Create an LSDChiTools from a raster.
    /// @param ThisRaster An LSDIndexRaster object
    /// @author SMM
    /// @date 24/05/2016
    LSDChiTools(LSDIndexRaster& ThisRaster)  { create(ThisRaster); }

    /// @brief Create an LSDChiTools from a LSDFlowInfo object.
    /// @param ThisFI An LSDFlowInfo object
    /// @author SMM
    /// @date 24/05/2016
    LSDChiTools(LSDFlowInfo& ThisFI)  { create(ThisFI); }

    /// @brief Create an LSDChiTools from a LSDJunctionNetwork.
    /// @param ThisJN An LSDJunctionNetwork object
    /// @author SMM
    /// @date 24/05/2016
    LSDChiTools(LSDJunctionNetwork& ThisJN)  { create(ThisJN); }

    /// @brief This resets all the data maps
    /// @author SMM
    /// :date 02/06/2016
    void reset_data_maps();

    /// @brief this gets the x and y location of a node at row and column
    /// @param row the row of the node
    /// @param col the column of the node
    /// @param x_loc the x location (Northing) of the node
    /// @param y_loc the y location (Easting) of the node
    /// @author SMM
    /// @date 22/12/2014
    void get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc);

    /// @biref This function write a file similar to the MCHISegmented one slightly different to match with the knickpoint plotting requirements
    /// @ param LSDFlowInfo A LSDFlowInfo object
    /// @param string filename the path name and extension of hte required file
    /// @author BG
    /// @date 05/12/2017
    void print_mchisegmented_knickpoint_version(LSDFlowInfo& FlowInfo, string filename);

    /// @brief this gets the x and y location of a node at row and column
    /// @param row the row of the node
    /// @param col the column of the node
    /// @param x_loc the x location (Northing) of the node
    /// @param y_loc the y location (Easting) of the node
    /// @author SMM
    /// @date 22/12/2014
    void get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc);

    /// @brief a function to get the lat and long of a node in the raster
    /// @detail Assumes WGS84 ellipsiod
    /// @param row the row of the node
    /// @param col the col of the node
    /// @param lat the latitude of the node (in decimal degrees, replaced by function)
    ///  Note: this is a double, because a float does not have sufficient precision
    ///  relative to a UTM location (which is in metres)
    /// @param long the longitude of the node (in decimal degrees, replaced by function)
    ///  Note: this is a double, because a float does not have sufficient precision
    ///  relative to a UTM location (which is in metres)
    /// @param Converter a converter object (from LSDShapeTools)
    /// @author SMM
    /// @date 24/05/2015
    void get_lat_and_long_locations(int row, int col, double& lat,
                  double& longitude, LSDCoordinateConverterLLandUTM Converter);


    /// @brief a function to get the lat and long of a coordinate point in the raster
    /// @detail Assumes WGS84 ellipsiod - does not correspond necessarly to a node
    /// @param row the row of the node
    /// @param col the col of the node
    /// @param lat the latitude of the node (in decimal degrees, replaced by function)
    ///  Note: this is a double, because a float does not have sufficient precision
    ///  relative to a UTM location (which is in metres)
    /// @param long the longitude of the node (in decimal degrees, replaced by function)
    ///  Note: this is a double, because a float does not have sufficient precision
    ///  relative to a UTM location (which is in metres)
    /// @param Converter a converter object (from LSDShapeTools)
    /// @author SMM
    /// @date 24/05/2015
    void get_lat_and_long_locations_from_coordinate(float X, float Y, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter);

    /// @brief this function gets the UTM_zone and a boolean that is true if
    /// the map is in the northern hemisphere
    /// @param UTM_zone the UTM zone. Replaced in function.
    /// @param is_North a boolean that is true if the DEM is in the northern hemisphere.
    ///  replaced in function
    /// @author SMM
    /// @date 22/12/2014
    void get_UTM_information(int& UTM_zone, bool& is_North);

    /// @brief This takes a chi raster and updates the chi data map.
    /// @detail WARNING you must use a raster derived from the topography
    ///  raster that was used to make the FlowInfo object. This function
    ///  does not check the dimensions of the raster
    /// @param FlowInfo An LSDFlowInfo object
    /// @param Chi_coordinate LSDRaster of the chi coordinate
    /// @author SMM
    /// @date 17/05/2017
    void update_chi_data_map(LSDFlowInfo& FlowInfo, LSDRaster& Chi_coord);

    /// @brief This takes a chi raster and updates the chi data map.
    ///  Overloaded from the previous function: this one calculates chi
    ///  directly from the FlowInfo so you have no problems with raster size
    /// @detail WARNING you must use a raster derived from the topography
    ///  raster that was used to make the FlowInfo object. This function
    ///  does not check the dimensions of the raster
    /// @param FlowInfo An LSDFlowInfo object
    /// @param A_0 the A_0 parameter: in metres^2 suggested value is 1
    /// @param m_over_n the m/n ratio
    /// @author SMM
    /// @date 17/05/2017
    void update_chi_data_map(LSDFlowInfo& FlowInfo, float A_0, float movern);




    /// @brief This recalcualtes chi for a single basin. Built for speed
    ///   and is used by the MCMC routines.
    /// @param FlowInfo An LSDFlowInfo object
    /// @param A_0 the A_0 parameter: in metres^2 suggested value is 1
    /// @param m_over_n the m/n ratio
    /// @param mimum_contributing_pixels This minimum number of contributing pixels needed before chi is calculated
    /// @param basin_key the key to the basin you want
    /// @param outlet_node_from_basin_key_map a map where the key is the basin key and the value is the node index of the outlet.
    ///   you generate this map by calling get_outlet_node_from_basin_key_map()
    /// @return No return, but the chi values FOR THIS BASIN ONLY are updated.
    /// @author SMM
    /// @date 14/07/2017
    void update_chi_data_map_for_single_basin(LSDFlowInfo& FlowInfo, float A_0, float movern,
                                     int minimum_contributing_pixels, int basin_key,
                                     map<int,int> outlet_node_from_basin_key_map);


    /// @brief This recalcualtes chi for a single basin. Built for speed
    ///   and is used by the MCMC routines.
    /// @detail This version uses a discharge raster
    /// @param FlowInfo An LSDFlowInfo object
    /// @param A_0 the A_0 parameter: in metres^2 suggested value is 1
    /// @param m_over_n the m/n ratio
    /// @param mimum_contributing_pixels This minimum number of contributing pixels needed before chi is calculated
    /// @param basin_key the key to the basin you want
    /// @param outlet_node_from_basin_key_map a map where the key is the basin key and the value is the node index of the outlet.
    ///   you generate this map by calling get_outlet_node_from_basin_key_map()
    /// @param Discharge an LSDRaster of discharge
    /// @return No return, but the chi values FOR THIS BASIN ONLY are updated.
    /// @author SMM
    /// @date 14/07/2017
    void update_chi_data_map_for_single_basin(LSDFlowInfo& FlowInfo, float A_0, float movern,
                                     int minimum_contributing_pixels, int basin_key,
                                     map<int,int> outlet_node_from_basin_key_map, LSDRaster& Discharge);


    /// @brief This function makes a chi map and prints to a csv file
    /// @detail the lat and long coordinates in the csv are in WGS84
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The string filename including path and extension
    /// @param A_0 the A_0 parameter
    /// @param m_over_n the m/n ratio
    /// @param area_threshold the threshold over which to print chi
    /// @author SMM
    /// @date 24/05/2016
    void chi_map_to_csv(LSDFlowInfo& FlowInfo, string filename,
                        float A_0, float m_over_n, float area_threshold);

    /// @brief This function takes a raster prints to a csv file
    /// @detail the lat and long coordinates in the csv are in WGS84
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The string filename including path and extension
    /// @param chi_coord the raster of the chi coordinate (printed elsewhere)
    /// @author SMM
    /// @date 03/06/2016
    void chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname, LSDRaster& chi_coord);

    /// @brief This function takes a raster prints to a csv file. Includes the junction number in the file
    /// @detail the lat and long coordinates in the csv are in WGS84
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The string filename including path and extension
    /// @param chi_coord the raster of the chi coordinate (printed elsewhere)
    /// @param basin_raster A raster with the basin numbers (calculated elsewhere)
    /// @author SMM
    /// @date 31/01/2017
    void chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname, LSDRaster& chi_coord, LSDIndexRaster& basin_raster);

    /// @brief This function is used to tag channels with a segment number
    ///  It decides on segments if the M_Chi value has changed so should only be used
    ///  with chi networks that have used a skip of 0 and a monte carlo itertions of 1
    ///  This data is used by other routines to look at the spatial distribution of
    ///  hillslope-channel coupling.
    /// @detail WARNING: ONLY use if you have segmented with skip 0 and iterations 1. Otherwise
    ///  you will get a new segment for every channel pixel
    /// @param FlowInfo an LSDFlowInfo object
    /// @param maximum_segment_length is the longest a segment is allowed to be before a new segment is created
    /// @author SMM
    /// @date 4/02/2017
    void segment_counter(LSDFlowInfo& FlowInfo, float maximum_segment_length);

    /// @brief This function is used to tag channels with a segment number
    ///  It decides on segments if the M_Chi value has changed so should only be used
    ///  with chi networks that have used a skip of 0 and a monte carlo itertions of 1
    ///  This data is used by other routines to look at the spatial distribution of
    ///  hillslope-channel coupling.
    /// @detail WARNING: ONLY use if you have segmented with skip 0 and iterations 1. Otherwise
    ///  you will get a new segment for every channel pixel
    /// @param FlowInfo an LSDFlowInfo object
    /// @param maximum_segment_length is the longest a segment is allowed to be before a new segment is created
    /// @return LSDIndexRaster showing stream network indexed by segment ID
    /// @author MDH
    /// @date 15/06/2017
    LSDIndexRaster segment_mapping(LSDFlowInfo& FlowInfo, float maximum_segment_length);

    /// @brief This function calculates the fitted elevations: It uses m_chi and b_chi
    ///  data to get the fitted elevation of the channel points.
    /// @param FlowInfo an LSDFlowInfo object
    /// @author SMM
    /// @date 4/02/2017
    void segment_counter_knickpoint(LSDFlowInfo& FlowInfo, float threshold_knickpoint, float threshold_knickpoint_length);

    /// @brief This function extract the difference,ratio,sign between each segments of the M_segmented_chi analysis
    /// @param FlowInfo an LSDFlowInfo object
    /// @author BG
    /// @date 4/02/2017
    void ksn_knickpoint_detection(LSDFlowInfo& FlowInfo);

    /// @brief This function extract the difference,ratio,sign between each segments of the M_segmented_chi analysis
    /// @param FlowInfo an LSDFlowInfo object
    /// @author BG
    /// @date 13/11/2017
    void knickzone_weighting_completion(map<pair<int,int>, vector<int> > mapofnode);

    /// @brief Development function based on segment_counter to help
    ///  knickpoint detection. More description will be added when it will be
    ///  functional.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param float threshold_knickpoint the knickpoints detection threshold
    /// @author BG
    /// @date 10/02/2017
    void calculate_segmented_elevation(LSDFlowInfo& FlowInfo);


    /// @brief This splits all the sources from the baselevels so that
    ///  individual baselevel catchemnts can be compared in sequence.
    ///  It produces a map where the sources for each baselelvel are
    ///   split into incremetally numberered (0,1,2) channels.
    /// @param n_sources_for_baselevel The number of sources for each baselelvel node
    ///   Replaced in function.
    /// @param index_into_sources_vec The index into the ordered sources vector
    ///   that is the starting index for each baselevel. Replaced in function.
    /// @author SMM
    /// @date 26/05/2017
    void baselevel_and_source_splitter(vector<int>& n_sources_for_baselevel,
                                                vector<int>& index_into_sources_vec);

    /// @brief This prints all the indexing and keys to screen for bug checking
    /// @author SMM
    /// @date 28/05/2017
    void print_basin_and_source_indexing_to_screen();

    /// @brief This returns an maximum liklihood estiamtor by comparing
    ///  a channel (with a particular source number) against a reference channel
    /// @param FlowInfo an LSDFlowInfo object
    /// @param reference_channel the source key of the reference channel
    /// @param test_channel the source key of the test channel
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    ///     If you have many nodes this number needs to be large
    /// @return The maximum likelihood estimator
    /// @author SMM
    /// @date 04/05/2017
    float test_segment_collinearity(LSDFlowInfo& FlowInfo, int reference_channel, int test_channel, float sigma);

    /// @brief This returns an maximum liklihood estiamtor by comparing
    ///  a channel (with a particular source number) against a reference channel, using specific points
    ///  on the test channel
    /// @param FlowInfo an LSDFlowInfo object
    /// @param reference_channel the source key of the reference channel
    /// @param test_channel the source key of the test channel
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    ///     If you have many nodes this number needs to be large
    /// @param chi_distances_to_test The distances in chi space to test
    /// @return The maximum likelihood estimator
    /// @author SMM
    /// @date 20/07/2017
    float test_segment_collinearity_using_points(LSDFlowInfo& FlowInfo, int reference_channel, int test_channel,
                                             float sigma, vector<float> chi_distances_to_test);

    /// @brief This takes a basin key and returns all the residuals of the
    ///  channels.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param only_use_mainstem_as_reference If true, only use the mainstem as a reference channel
    /// @param test_channel the source key of the test channel
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    ///     If you have many nodes this number needs to be large
    /// @param baselevel_key The key to the basin you want
    /// @return A vector containing all the residuals
    /// @author SMM
    /// @date 21/07/2017
    vector<float> retrieve_all_residuals_by_basin(LSDFlowInfo& FlowInfo, bool only_use_mainstem_as_reference,
                                                 int baselevel_key);


    /// @brief This computes a collinearity metric for all combinations of
    ///  channels for a given basin
    /// @detail It takes all the combinations of sources and gets the goodness of fit between each pair
    ///  of sources.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param only_use_mainstem_as_reference True if you only want to use the mainstem
    /// @param basin_key The key into the basin you want to test all collinearity of.
    /// @param reference_source integer vector replaced in function that has the reference vector for each comparison
    /// @param test_source integer vector replaced in function that has the test vector for each comparison
    /// @param MLE_values the MLE for each comparison. Replaced in function.
    /// @param RMSE_values the RMSE for each comparison (i.e. between source 0 1, 0 2, 0 3, etc.). Replaced in function.
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    ///     If you have many nodes this number needs to be large
    /// @author SMM
    /// @date 08/05/2017
    float test_all_segment_collinearity_by_basin(LSDFlowInfo& FlowInfo, bool only_use_mainstem_as_reference,
                                        int basin_key,
                                        vector<int>& reference_source, vector<int>& test_source,
                                        vector<float>& MLE_values, vector<float>& RMSE_values,
                                        float sigma);

    /// @brief This computes a collinearity metric for all combinations of
    ///  channels for a given basin. This version uses specified points in chi space
    ///  to compare against other channels rather than the entire channel
    /// @detail It takes all the combinations of sources and gets the goodness of fit between each pair
    ///  of sources.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param only_use_mainstem_as_reference True if you only want to use the mainstem
    /// @param basin_key The key into the basin you want to test all collinearity of.
    /// @param reference_source integer vector replaced in function that has the reference vector for each comparison
    /// @param test_source integer vector replaced in function that has the test vector for each comparison
    /// @param MLE_values the MLE for each comparison. Replaced in function.
    /// @param RMSE_values the RMSE for each comparison (i.e. between source 0 1, 0 2, 0 3, etc.). Replaced in function.
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    /// @param chi_distances_to_test The distances in chi space to test
    /// @author SMM
    /// @date 20/07/2017
    float test_all_segment_collinearity_by_basin_using_points(LSDFlowInfo& FlowInfo, bool only_use_mainstem_as_reference,
                                        int basin_key,
                                        vector<int>& reference_source, vector<int>& test_source,
                                        vector<float>& MLE_values, vector<float>& RMSE_values,
                                        float sigma, vector<float> chi_distances_to_test);

    /// @brief This computes a the disorder metric of Hergarten et al 2016 by basin
    /// @param FlowInfo an LSDFlowInfo object
    /// @param basin_key The key into the basin you want to test all collinearity of.
    /// @author SMM
    /// @date 24/03/2018
    float test_collinearity_by_basin_disorder(LSDFlowInfo& FlowInfo, 
                                        int basin_key);


    /// @brief This computes a the disorder metric of Hergarten et al 2016 by basin.
    /// It uses a permutation algorithm to find all combinations of tributary channels
    /// and computes the disorder statistic of each of these
    /// @param FlowInfo an LSDFlowInfo object
    /// @param basin_key The key into the basin you want to test all collinearity of.
    /// @author SMM
    /// @date 13/04/2018
    vector<float> test_collinearity_by_basin_disorder_with_uncert(LSDFlowInfo& FlowInfo, 
                                        int basin_key);

    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics.
    /// @detail This gets the median residual. The best fit will have a median residual
    ///  closest to zero
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @author SMM
    /// @date 21/07/2017
     void calculate_goodness_of_fit_collinearity_fxn_movern_using_median_residuals(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix);



    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    ///     If you have many nodes this number needs to be large
    /// @author SMM
    /// @date 16/05/2017
    /// MODIFIED FJC 17/06/17 to take a junction network as an argument - need to print out the outlet
    /// junction of each basin to match to the basin key for visualisation
    void calculate_goodness_of_fit_collinearity_fxn_movern(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma);




    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics.
    ///  Same as above but can use a discharge raster to calculate chi
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @param Discharge and LSDRaster of discharge
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE
    ///     If you have many nodes this number needs to be large
    /// @author SMM
    /// @date 16/05/2017
    void calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge(LSDFlowInfo& FlowInfo,
                        LSDJunctionNetwork& JN, float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix,
                        LSDRaster& Discharge, float sigma);





    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics.
    /// @detail Uses discrete points rather than all the tributary data
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE.
    /// @param chi_distance_fractions These are a vector of fractions of distance of the
    ///   length of the mainstem that you want to sample the tributaries. For example if the
    ///   mainstem is 10 m long in chi space and you have {0.1,0.2,0.3} as the fraction you
    ///   sample the tributaries 1, 2 and 3 metres upstream of their confluence.
    /// @author SMM
    /// @date 20/07/2017
    void calculate_goodness_of_fit_collinearity_fxn_movern_using_points(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma,
                        vector<float> chi_distance_fractions);

    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics.
    ///  Same as above but can use a discharge raster to calculate chi
    /// @detail Uses discrete points rather than all the tributary data
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @param Discharge and LSDRaster of discharge
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE.
    /// @param chi_distance_fractions These are a vector of fractions of distance of the
    ///   length of the mainstem that you want to sample the tributaries. For example if the
    ///   mainstem is 10 m long in chi space and you have {0.1,0.2,0.3} as the fraction you
    ///   sample the tributaries 1, 2 and 3 metres upstream of their confluence.
    /// @author SMM
    /// @date 20/07/2017
    void calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_points(LSDFlowInfo& FlowInfo,
                        LSDJunctionNetwork& JN, float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix,
                        LSDRaster& Discharge, float sigma,
                        vector<float> chi_distance_fractions);



    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics. 
    /// @detail Uses discrete points rather than all the tributary data. It uses monte carlo
    ///   sampling to get the points, so one can repeatedly sample the MLE values for
    ///   a fixed number of points
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE.
    /// @param n_fracs the number of chi distance fractions you want to use
    /// @param MC_iteration The number of iteration you want to use
    /// @param max_frac the maximum chi fraction you want to examine.
    /// @author SMM
    /// @date 24/03/2017
    void calculate_goodness_of_fit_collinearity_fxn_movern_using_points_MC(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma,
                        int n_fracs,
                        int MC_iterations,
                        float max_frac);


    /// @brief This wraps the collinearity tester, looping through different m over n
    ///  values and calculating goodness of fit statistics. This one uses discharge
    /// @detail Uses discrete points rather than all the tributary data. It uses monte carlo
    ///   sampling to get the points, so one can repeatedly sample the MLE values for
    ///   a fixed number of points
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN an LSDJunctionNetwork object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param only_use_mainstem_as_reference a boolean, if true only compare channels to mainstem .
    /// @param The file prefix for the data files
    /// @param sigma The uncertainty for the MLE calculation. In practice this simply scales MLE.
    /// @param n_fracs the number of chi distance fractions you want to use
    /// @param MC_iteration The number of iteration you want to use
    /// @param max_frac the maximum chi fraction you want to examine.
    /// @param Discharge A discharge raster
    /// @author SMM
    /// @date 27/04/2018
    void calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_points_MC(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma,
                        int n_fracs,
                        int MC_iterations,
                        float max_frac, LSDRaster& Discharge);


    /// @brief This wraps the disorder collinearity tester, looping through different concavity
    ///  values and calculating disorder statistics.
    /// @detail Uses the disorder metric of Hergarten et al ESURF 2016
    /// @param FlowInfo an LSDFlowInfo object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param The file prefix for the data files
    /// @param use_uncert a bool that if true triggers the uncertainty algorithms
    /// @author SMM
    /// @date 24/03/2018
    void calculate_goodness_of_fit_collinearity_fxn_movern_using_disorder(LSDFlowInfo& FlowInfo,
                        LSDJunctionNetwork& JN, float start_movern, float delta_movern, int n_movern,
                        string file_prefix, bool use_uncert);

    /// @brief This wraps the disorder collinearity tester, looping through different concavity
    ///  values and calculating disorder statistics.
    /// @detail Uses the disorder metric of Hergarten et al ESURF 2016. Uses a discharge raster to calculate chi
    /// @param FlowInfo an LSDFlowInfo object
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param The file prefix for the data files
    /// @param use_uncert a bool that if true triggers the uncertainty algorithms
    /// @param Discharge and LSDRaster with the discharge. 
    /// @author SMM
    /// @date 27/04/2018
    void calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_disorder(LSDFlowInfo& FlowInfo,
                        LSDJunctionNetwork& JN, float start_movern, float delta_movern, int n_movern,
                        string file_prefix, bool use_uncert, LSDRaster& Discharge);

    /// @brief This function drives a Monte Carlo-Markov chain model It wraps
    /// the dmovern tuning and the main chain.
    /// @param FlowInfo An LSDFlowInfo object
    /// @param minimum_contributing_pixel chi is only calculated if the contributing pixels are bigger than this
    /// @param sigma The sigma value for checking the MLE of chi
    /// @param movern_minimum The minimum movern value to be tested
    /// @param movern_maximum The maximum movern value to be tested
    /// @param N_chain_links The number of iterations in the chain you want
    /// @param OUT_DIR the output directory where you want the file
    /// @param OUT_ID prefix of the output file
    /// @param use_points a bool that if true means you use the point version of the collinearity test
    /// @return No return but makes chaing files with extension _chain.csv
    /// @author SMM
    /// @date 17/07/2017
    void MCMC_driver(LSDFlowInfo& FlowInfo, int minimum_contributing_pixels, float sigma,
                                 float movern_minimum, float movern_maximum,
                                 int N_chain_links,
                                 string OUT_DIR, string OUT_ID, bool use_points);

    /// @brief This function drives a Monte Carlo-Markov chain model and tries to tune
    ///  the dmovern_stddev value to get between 20-30% acceptance rate
    /// @param FlowInfo An LSDFlowInfo object
    /// @param minimum_contributing_pixel chi is only calculated if the contributing pixels are bigger than this
    /// @param NIterations The number of iterations in the chain you want
    /// @param sigma The sigma value for checking the MLE of chi
    /// @param movern_minimum The minimum movern value to be tested
    /// @param movern_maximum The maximum movern value to be tested
    /// @param basin key The key of the basin to be tested
    /// @param use_points a bool that if true means you use the point version of the collinearity test
    /// @return The tuned dmovern_stddev.
    /// @author SMM
    /// @date 13/07/2017
    float MCMC_for_movern_tune_dmovern(LSDFlowInfo& FlowInfo, int minimum_contributing_pixels, float sigma,
                                 float movern_minimum, float movern_maximum,
                                 int basin_key, bool use_points);

    /// @brief This function drives a Monte Carlo-Markov chain model and tries to tune
    ///  the sigma value to get between 20-30% acceptance rate
    /// @param FlowInfo An LSDFlowInfo object
    /// @param minimum_contributing_pixel chi is only calculated if the contributing pixels are bigger than this
    /// @param dmovernstddev The desired stddev of m/n changes
    /// @param movern_minimum The minimum movern value to be tested
    /// @param movern_maximum The maximum movern value to be tested
    /// @param basin key The key of the basin to be tested
    /// @param use_points a bool that if true means you use the point version of the collinearity test
    /// @return sigma The tuned sigma value for checking the MLE of chi
    /// @author SMM
    /// @date 17/07/2017
    float MCMC_for_movern_tune_sigma(LSDFlowInfo& FlowInfo, int minimum_contributing_pixels,
                                     float dmovernstddev,
                                     float movern_minimum, float movern_maximum,
                                     int basin_key, bool use_points);


    /// @brief This function drives a Monte Carlo-Markov chain model for getting
    ///  the confidence intervals on the m/n value.
    /// @detail You can turn the chain file printing off at first to tune the
    ///  dmovern_sigma so that it arrices at a 25% acceptance rate
    /// @param ChainFname The name of the chain file (with path)
    /// @param printChain  If true, the chain file is printed
    /// @param FlowInfo An LSDFlowInfo object
    /// @param minimum_contributing_pixel chi is only calculated if the contributing pixels are bigger than this
    /// @param NIterations The number of iterations in the chain you want
    /// @param sigma The sigma value for checking the MLE of chi
    /// @param dmovernstddev The standard deviation in the change in m/n test values.
    ///    This needs to be tuned so the acceptance rate is ~25%
    /// @param movern_minimum The minimum movern value to be tested
    /// @param movern_maximum The maximum movern value to be tested
    /// @param basin key The key of the basin to be tested
    /// @param use_points a bool that if true means you use the point version of the collinearity test
    /// @return acceptanc probability.
    /// @author SMM
    /// @date 13/07/2017
    float MCMC_for_movern(string ChainFname, bool printChain, LSDFlowInfo& FlowInfo,
                          int minimum_contributing_pixels, int NIterations, float sigma, float dmovern_stddev,
                          float movern_minimum, float movern_maximum, int basin_key, bool use_points);


    /// @brief This prints a series of chi profiles as a function of m over n
    ///  for visualisation
    /// @param FlowInfo an LSDFlowInfo object
    /// @param file_prefix THe path and name of file without extension
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @author SMM
    /// @date 17/05/2017
    void print_profiles_as_fxn_movern(LSDFlowInfo& FlowInfo,string file_prefix, float start_movern, float delta_movern, int n_movern);

    /// @brief This prints a series of chi profiles as a function of m over n
    ///  for visualisation. It also burns a raster value to each of the nodes
    /// @detail The raster burning is useful for adding information like geology or K values
    /// @param FlowInfo an LSDFlowInfo object
    /// @param file_prefix THe path and name of file without extension
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param BurnRaster the raster to burn to the csv
    /// @param column_name the name of the column to which the data will be burned
    /// @author SMM
    /// @date 07/09/2017
    void print_profiles_as_fxn_movern_with_burned_raster(LSDFlowInfo& FlowInfo,
                          string filename, float start_movern, float delta_movern,
                          int n_movern, LSDRaster& BurnRaster, string column_name);
    
    void print_profiles_as_fxn_movern_with_secondary_raster(LSDFlowInfo& FlowInfo,
                          string filename, float start_movern, float delta_movern,
                          int n_movern, LSDRaster& BurnRaster, LSDRaster& SecondaryBurnRaster,string column_name, string secondary_column_name);

    /// @brief This prints a series of chi profiles as a function of mover
    ///  for visualisation
    /// @param FlowInfo an LSDFlowInfo object
    /// @param file_prefix THe path and name of file without extension
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param Discharge an LSDRaster of discharge
    /// @author SMM
    /// @date 17/05/2017
    void print_profiles_as_fxn_movern_with_discharge(LSDFlowInfo& FlowInfo,string file_prefix,
                                   float start_movern, float delta_movern,
                                   int n_movern, LSDRaster& Discharge);

    /// @brief This prints a series of chi profiles as a function of mover
    ///  for visualisation. It also burns a raster value to each of the nodes
    /// @detail The raster burning is useful for adding information like geology or K values
    /// @param FlowInfo an LSDFlowInfo object
    /// @param file_prefix THe path and name of file without extension
    /// @param start_movern the starting m/n ratio
    /// @param delta_movern the change in m/n
    /// @param n_novern the number of m/n values to use
    /// @param Discharge an LSDRaster of discharge
    /// @param BurnRaster the raster to burn to the csv
    /// @param column_name the name of the column to which the data will be burned
    /// @author SMM
    /// @date 08/09/2017
    void print_profiles_as_fxn_movern_with_discharge_and_burned_raster(LSDFlowInfo& FlowInfo, string filename,
           float start_movern, float delta_movern, int n_movern, LSDRaster& Discharge,
           LSDRaster& BurnRaster, string burned_column_name);
    
    //adding support for a secondary raster
    void print_profiles_as_fxn_movern_with_discharge_and_secondary_raster(LSDFlowInfo& FlowInfo, string filename,
           float start_movern, float delta_movern, int n_movern, LSDRaster& Discharge,
           LSDRaster& BurnRaster,LSDRaster& SecondaryBurnRaster, string burned_column_name,string secondary_burned_column_name);

    /// @brief This inverst the key_to_baselevel map
    ///  so that the key is the baselevel key and the value is the node index of the outlet node
    /// @return An <int,int> map where key is baselevel key and value is node index of outlet node
    /// @author SMM
    /// @date 14/07/2017
    map<int,int> get_outlet_node_from_basin_key_map();

    /// @brief This gets the node index (the reference into LSDFlowInfo) of a source
    ///  based on a source key
    /// @param source_key the source key of the reference channel
    /// @return the node index of the source node
    /// @author SMM
    /// @date 06/05/2017
    int get_source_from_source_key(int source_key);

    /// @brief This gets the index into the node_sequence vector of the first
    ///  node in a channel identified by its source key
    /// @param source_key the source key of the reference channel
    /// @return the index into the node_sequence vector of the source node of the channel with source_key
    /// @author SMM
    /// @date 04/05/2017
    int get_starting_node_of_source(int source_key);

    /// @brief Gets the number of channels in the DEM
    /// @return number of channels
    /// @author SMM
    /// @date 05/05/2017
    int get_number_of_channels();

    /// @brief This takes a source key and a flow info object and overwrites vectors
    ///  containing chi and elevation data from a channel tagged by a source
    ///  key. The idea is to use this in the MLE comparison between two channels
    ///  to check for collinearity
    /// @param FlowInfo and LSDFlowInfo object
    /// @param source_key The key of the source
    /// @param chi_data A vector holding chi data of the channel. Will be overwritten
    /// @param elevation_data A vector holding elevation data of the channel. Will be overwritten
    /// @author SMM
    /// @date 06/05/2017
    void get_chi_elevation_data_of_channel(LSDFlowInfo& FlowInfo, int source_key,
                                vector<float>& chi_data, vector<float>& elevation_data);

    /// @brief This takes the chi locations of a tributarry vector and then uses
    ///  linear interpolation to determine the elevation on a reference channel
    ///  at those chi values
    /// @param reference_chi the chi coordiantes of the reference channel
    /// @param reference_elevation the elevations on the reference channel
    /// @param trib_chi the chi coordiantes of the tributary channel
    /// @param trib_elevation the elevations on the tributary channel
    /// @return A vector of the elevations on the chi locations of the tributary channel
    /// @author SMM
    /// @date 07/05/2017
    vector<float> project_data_onto_reference_channel(vector<float>& reference_chi,
                                 vector<float>& reference_elevation, vector<float>& trib_chi,
                                 vector<float>& trib_elevation);


    /// @brief This takes the chi locations of a tributarry vector and then uses
    ///  linear interpolation to determine the elevation on a reference channel
    ///  of points at a vector of fixed distances upstream the tributary channel
    /// @param reference_chi the chi coordiantes of the reference channel
    /// @param reference_elevation the elevations on the reference channel
    /// @param trib_chi the chi coordiantes of the tributary channel
    /// @param trib_elevation the elevations on the tributary channel
    /// @param chi_distances_to_test the distances upstream of the confluence on the tributary channel to test the residuals
    /// @return A vector of the elevations on the chi locations of the tributary channel
    /// @author SMM
    /// @date 20/07/2017
    vector<float> project_points_onto_reference_channel(vector<float>& reference_chi,
                                 vector<float>& reference_elevation, vector<float>& trib_chi,
                                 vector<float>& trib_elevation, vector<float> chi_distances_to_test);




    /// @brief This performs slope area analysis. It goes down through each
    ///  source node and collects S-A data along these channels.
    ///  It uses the suggested appraoch of Wobus et al. 2006 in that it uses
    ///  a drop interval to measure slope.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param vertical_interval the mean intervale over which slope is measured
    /// @param midpoint_nodes The node indices of the places where slope is calculated.
    ///  This is replaced in the function.
    /// @param Slopes the slopes. This is replaced in the function.
    /// @author SMM
    /// @date 31/05/2017
    void get_slope_area_data(LSDFlowInfo& FlowInfo, float vertical_interval,
                             vector<int>& midpoint_nodes, vector<float>& slopes);


    /// @brief This function takes raw S-A data (generated by get_slope_area_data)
    ///  and runs a bootstrapping procedure to find the confidence intervals
    ///  of the regression coeffcients and the intercepts.
    /// @brief Note that it converts to log space before running the regression
    /// @param FlowInfo The LSDFlowInfo object
    /// @param SA_midpoint_node The node index of the midpoints of all the data
    /// @param SA_slope: The slope of the data at the midpoints (the slope is)
    ///  averaged over several pixels
    /// @param N_iterations: The number of bootstrap iterations.
    /// @param keep_data_prob The probability that you will keep any given data point in the data set.
    /// @param filename The name of file (with extension and directory) for printing
    /// @author SMM
    /// @date 17/07/2017
    void bootstrap_slope_area_data(LSDFlowInfo& FlowInfo,
                                          vector<int>& SA_midpoint_node,
                                          vector<float>& SA_slope,
                                          int N_iterations,
                                          float bootstrap_keep_data_prob,
                                          string filename);

    /// @detail This takes slope area data and bins the data so that we can
    ///  pretend horrible, noisy S-A data is adequate for understanding
    ///  channel behaviour.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param vertical_interval the mean intervale over which slope is measured
    /// @param midpoint_nodes The node indices of the places where slope is calculated.
    ///  This is replaced in the function.
    /// @param Slopes the slopes. This is replaced in the function.
    /// @param log_bin_width The width of the bins (in log A)
    /// @param filename The name of the output file (with path and extension)
    /// @author SMM
    /// @date 31/05/2017
    void bin_slope_area_data(LSDFlowInfo& FlowInfo, vector<int>& SA_midpoint_node,
                             vector<float>& SA_slope, float log_bin_width, string filename);



    /// @detail This takes slope area data and bins the data so that we can
    ///  pretend horrible, noisy S-A data is adequate for understanding
    ///  channel behaviour. It then segments these horrible data using the
    /// segmentation algorithm.
    /// @detail Happy 4th of July everyone!
    /// @param FlowInfo an LSDFlowInfo object
    /// @param vertical_interval the mean intervale over which slope is measured
    /// @param midpoint_nodes The node indices of the places where slope is calculated.
    ///  This is replaced in the function.
    /// @param Slopes the slopes. This is replaced in the function.
    /// @param log_bin_width The width of the bins (in log A)
    /// @param minimum_segment_length Minimum segment length for segmentation algorithm
    /// @param filename The name of the output file (with path and extension)
    /// @author SMM
    /// @date 04/07/2017
    void segment_binned_slope_area_data(LSDFlowInfo& FlowInfo,
                                          vector<int>& SA_midpoint_node,
                                          vector<float>& SA_slope,
                                          float log_bin_width,
                                          int minimum_segment_length,
                                          string filename);

    /// @brief This takes the midpoint node and slope vectors produced by the slope_area_analysis
    ///  and prints them to a csv
    /// @param SA_midpoint_node the node index of the midpoints used in to caluclate slope
    /// @param SA_slope The slope data
    /// @param filename The name (including path and extension) of the file for printing
    /// @author SMM
    /// @date 31/05/2017
    void print_slope_area_data_to_csv(LSDFlowInfo& FlowInfo,
                                              vector<int>& SA_midpoint_node,
                                              vector<float>& SA_slope,
                                              string filename);

    /// @brief This function burns the chi coordinate (and area, flow distance and elevation)
    ///  onto the data maps in the chitools object. It does not do any segmentation.
    /// @detail The purpose of this function is to get the chi coordinate without
    ///   calculating m_chi or segmenting, and is used for m/n calculations. Can
    ///   also be used for maps of chi coordinate
    /// @param FlowInfo an LSDFlowInfo object
    /// @param source_nodes a vector containing the sorted sorce nodes (by flow distance)
    /// @param outlet_nodes a vector continaing the outlet nodes
    /// @param baselevel_node_of_each_basin a vector continaing the baselelve node of the basin for each channel
    /// @param Elevation an LSDRaster containing elevation info
    /// @param DistanceFromOutlet an LSDRaster with the flow distance
    /// @param DrainageArea an LSDRaster with the drainage area
    /// @author SMM
    /// @date 16/05/2017
    void chi_map_automator_chi_only(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate);

    /// @brief This function maps out the chi steepness and other channel
    ///  metrics in chi space from all the sources supplied in the
    ///  source_nodes vector. The source and outlet nodes vector is
    ///  generated by LSDJunctionNetwork.get_overlapping_channels
    /// @detail Takes vector so source and outlet nodes and performs the segment
    ///  fitting routine on them
    /// @param FlowInfo an LSDFlowInfo object
    /// @param source_nodes a vector containing the sorted sorce nodes (by flow distance)
    /// @param outlet_nodes a vector continaing the outlet nodes
    /// @param baselevel_node_of_each_basin a vector continaing the baselelve node of the basin for each channel
    /// @param Elevation an LSDRaster containing elevation info
    /// @param DistanceFromOutlet an LSDRaster with the flow distance
    /// @param DrainageArea an LSDRaster with the drainage area
    /// @param target_nodes int the target number of nodes in a break
    /// @param n_iterations  int the number of iterations
    /// @param target_skip int the mean skipping value
    /// @param minimum_segment_length How many nodes the mimimum segment will have.
    /// @param sigma Standard deviation of error on elevation data
    /// @author SMM
    /// @date 23/05/2016
    void chi_map_automator(LSDFlowInfo& FlowInfo, vector<int> source_nodes,
                           vector<int> outlet_nodes, vector<int> baselevel_node_of_each_basin,
                           LSDRaster& Elevation, LSDRaster& FlowDistance,
                           LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                           int target_nodes, int n_iterations, int skip,
                           int minimum_segment_length, float sigma);

    /// @brief This function maps out the chi steepness and other channel
    ///  metrics in chi space from all the sources supplied in the
    ///  source_nodes vector. The source and outlet nodes vector is
    ///  generated by LSDJunctionNetwork.get_overlapping_channels
    /// @detail This is simpler than the above function: it simply performs
    ///   a linear regression ofve a fixed number of data points: no effort is
    ///   made to segment the data. It is thus much closer to a k_sn plot
    /// @param FlowInfo an LSDFlowInfo object
    /// @param source_nodes a vector continaing the sorted sorce nodes (by flow distance)
    /// @param outlet_nodes a vector continaing the outlet nodes
    /// @param baselevel_node_of_each_basin a vector continaing the baselelve node of the basin for each channel
    /// @param Elevation an LSDRaster containing elevation info
    /// @param DistanceFromOutlet an LSDRaster with the flow distance
    /// @param DrainageArea an LSDRaster with the drainage area
    /// @param regression_nodes the number of nodes in each segment over which
    ///   to perform a linear regression. This number should be odd to it has
    ///   a clear midpoint
    /// @author SMM
    /// @date 02/06/2016
    void chi_map_automator_rudimentary(LSDFlowInfo& FlowInfo, vector<int> source_nodes, vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int regression_nodes);

    /// @brief This returns an LSDIndexRaster with basins numbered by outlet junction
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN the junction network object
    /// @param Junctions The baselevel junctions to be printed
    /// @return The basin raster
    /// @author SMM
    /// @date 19/01/2017
    LSDIndexRaster get_basin_raster(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Juntions);

    /// @brief This return a map with the basin ID as key, and map of the count of the different lithology as key/value.
    /// @param FlowInfo
    /// @param JunctionNetwork
    /// @param LSDIndexRaster as lithologic or geologic map
    /// @parma vector of baselevel junctions
    /// @author BG
    /// @date 15/09/2017
    map<int,map<int,int> > get_basin_lithocount(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork, LSDIndexRaster& litho,
                                                              vector<int> Junctions);


    /// @brief This prints a csv file that has the locations of the sources and their keys
    ///  latitude,longitude,source_node, source_key
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author SMM
    /// @date 16/01/2017
    void print_source_keys(LSDFlowInfo& FlowInfo, string filename);

    /// @brief This prints a csv file that has the locations of the baselevels and their keys
    ///  latitude,longitude,baselevel_junctione, baselevel_key
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN the junction network object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author SMM
    /// @date 16/01/2017
    void print_baselevel_keys(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN, string filename);

    /// @brief This prints a basin LSDIndexRaster with basins numbered by outlet junction
    ///  and a csv file that has the latitude and longitude of both the outlet and the centroid
    /// @param FlowInfo an LSDFlowInfo object
    /// @param JN the junction network object
    /// @param Junctions The baselevel junctions to be printed
    /// @param base_filename The name of the filename to print to (should have full
    ///   path but no extension. The "_AllBasins" will be added
    /// @author SMM
    /// @date 19/01/2017
    void print_basins(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Juntions, string base_filename);

    /// @brief This prints a csv file with chi data from the data maps
    ///  the columns are:
    ///  latitude,longitude,chi,elevation,flow distance,drainage area,
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author SMM
    /// @date 05/06/2017
    void print_chi_data_map_to_csv(LSDFlowInfo& FlowInfo, string filename);


    /// @brief This prints a csv file with chi data from the data maps for a specific basin
    ///  the columns are:
    ///  latitude,longitude,chi,elevation,flow distance,drainage area,
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @param basin_key the basin key.
    /// @author SMM
    /// @date 14/07/2017
    void print_chi_data_map_to_csv_for_single_basin(LSDFlowInfo& FlowInfo, string filename, int basin_key);

    /// @brief Print a csv file with basin_key and the number of lithology pixels per lithology ID. You need the files from the rasterisation to decrypt it.
    /// @param FlowInfo
    /// @param path+name of the file
    /// @param map of map of litho obtain from a ChiTool.count_unique_values_from_litho_raster function for example
    /// @author BG
    /// @date 15/09/2017
    void simple_litho_basin_to_csv(LSDFlowInfo& FlowInfo, string csv_slbc_fname,
                                      map<int,map<int,int> > map_slbc);

    /// @brief Print a csv file with basin_key and the number/percentages of lithology pixels per lithology ID. You need the files from the rasterisation to decrypt it.
    /// @param FlowInfo
    /// @param path+name of the file
    /// @param map of map of litho obtain from a ChiTool.count_unique_values_from_litho_raster function for example
    /// @author BG
    /// @date 15/09/2017
    void extended_litho_basin_to_csv(LSDFlowInfo& FlowInfo, string csv_slbc_fname,
                                      map<int,map<int,int> > map_slbc);

    /// @brief This prints a csv file with all the data from the data maps
    ///  the columns are:
    ///  latitude,longitude,chi,elevation,flow distance,drainage area,m_chi,b_chi
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author SMM
    /// @date 02/06/2016
    void print_data_maps_to_file_full(LSDFlowInfo& FlowInfo, string filename);

    /// @brief This prints a csv file with all the knickpoint data
    ///  the columns are:
    ///  latitude,longitude,elevation,flow distance,drainage area,ratio,diff,sign
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author BG
    /// @date 06/06/2017
    void print_knickpoint_to_csv(LSDFlowInfo& FlowInfo, string filename);

    /// @brief This prints a csv file with all the segmented gradient data (adaptation of Mudd et al., for river gradient )
    ///  the columns are what they are
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author BG
    /// @date Today
    void print_segmented_gradient_to_csv(LSDFlowInfo& FlowInfo, string filename);


    /// @brief This prints a csv file with all the knickzones raw data
    ///  the columns are:
    ///  latitude,longitude,elevation,flow distance,drainage area,ratio,diff,sign
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author BG
    /// @date 06/06/2017
    void print_knickzone_to_csv(LSDFlowInfo& FlowInfo, string filename);

    /// @brief This prints a csv file with a subset of the data from the data maps
    ///  the columns are:
    ///  latitude,longitude,m_chi,b_chi
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author SMM
    /// @date 02/06/2016
    void print_data_maps_to_file_full_knickpoints(LSDFlowInfo& FlowInfo, string filename);

    /// @brief This prints a csv file with a subset of the data from the data maps
    ///  the columns are:
    ///  latitude,longitude,m_chi,b_chi,knickpoint
    /// Development function
    /// @param FlowInfo an LSDFlowInfo object
    /// @param filename The name of the filename to print to (should have full
    ///   path and the extension .csv
    /// @author SMM/BG
    /// @date 02/06/2016
    void print_data_maps_to_file_basic(LSDFlowInfo& FlowInfo, string filename);


    /// @brief fill maps containing information of receiving rivers
    /// @param FlowInfo: a FlowInfo object
    /// @author BG
    /// @date 30/11/2017 
    void get_previous_mchi_for_all_sources(LSDFlowInfo& Flowinfo);

    /// @brief It should find the ending (aka downstair) node 
    /// @param FlowInfo: a FlowInfo object
    /// @author BG
    /// @date 30/11/2017 
    int get_ending_node_of_source(LSDFlowInfo& FlowInfo, int source_key);

    /// @brief print a csv file with the receiver of each source and the corresponding source with the m_chi
    /// That is barely understable, however I have a cold so I am tired as F. Just ask me if you need more info about that
    /// @param string filename: the path/name.csv of your file
    /// @author BG
    /// @date 30/11/2017
    void print_intersources_mchi_map(string filename);


    /// @brief set a map of each source_key with the corresponding vector fo node INCLUDING the first node of the 
    /// receaving river if this abovementioned one does exist.
    /// @param Flowinfo, a LSDFlowInfo object
    /// @param int n_nodlump, the half lumping window
    /// @author BG
    /// @date 05/01/2018
    void set_map_of_source_and_node(LSDFlowInfo& FlowInfo, int n_nodlump);

    /// @brief Main function for the knickpoint analysis that control the call of other function to keep it clear and up to date
    /// @param FlowiInfo: a LSDFlowInfo object
    /// @param OUT_DIR: string containing the output directory path
    /// @param OUT_ID: string containing the output prefix
    /// @author BG
    /// @date 05/01/2018
    void ksn_knickpoint_automator(LSDFlowInfo& FlowInfo, string OUT_DIR, string OUT_ID, float MZS_th, float lambda_TVD, float lambda_TVD_b_chi,int stepped_combining_window,int window_stepped, float n_std_dev, int kp_node_search);

    void ksn_knickpoint_outlier_automator(LSDFlowInfo& FlowInfo, float MZS_th);

    /// @brief Dealing with composite knickpoints
    /// @param FlowiInfo: a LSDFlowInfo object
    /// @author BG
    /// @date 17/01/2018
    void ksn_knickpoints_combining(LSDFlowInfo& Flowinfo, int kp_node_search);

    /// @brief Detection of the knickpoints by looping through the source keys
    /// @param FlowiInfo: a LSDFlowInfo object
    /// @param OUT_DIR: string containing the output directory path
    /// @param OUT_ID: string containing the output prefix
    /// @author BG
    /// @date 05/01/2018
    void ksn_knickpoint_detection_new(LSDFlowInfo& FlowInfo);

    /// @brief increment the knickpoints for one river
    /// @param SK: the source key
    /// @param vecnode: a vector of the rive nodes
    /// @author BG
    /// @date 05/01/2018
    void ksn_knickpoint_raw_river(int SK, vector<int> vecnode);

    /// @brief write a file with the raw knickpoint informations
    /// @param FlowiInfo: a LSDFlowInfo object
    /// @param filename: string with path+name+.csv
    /// @author BG
    /// @date 05/01/2018
    void print_raw_ksn_knickpoint(LSDFlowInfo& FlowInfo, string filename);

    /// @brief Calculate KDE over the river system using the dksn/dchi map previously calculated through ksn_knickpoint_automator
    /// @author BG
    /// @date 05/01/2018
    void ksn_kp_KDE();

    /// @brief communicate with LSDStatTools to get the KDE oer river, also register the bandwidth automatically calculated
    /// @param vecnode: a vector of node index containing the data
    /// @param SK: source key
    /// @author BG
    /// @date 05/01/2018
    void KDE_vec_node_mchi(vector<int> vecnode, int SK);

    /// @brief write a file containing source key and bandwith calculated from the KDE
    /// @param vecnode: a vector of node index containing the data
    /// @param SK: source key
    /// @author BG
    /// @date 05/01/2018
    void print_bandwidth_ksn_knickpoint(string filename);

    /// @brief compute basics metrics per source keys: flow length, chi length plus more to come probably
    /// @param FlowInfo: a LSDFlowInfo Object
    /// @author BG
    /// @date 08/01/2018
    void compute_basic_matrics_per_source_keys(LSDFlowInfo& FlowInfo);


    /// @brief lump the m_chi to average the MC noise
    /// @param n_nodlump: the lumping half_window (number of nodes)
    /// @author BG
    /// @date 08/01/2018
    void lump_my_ksn(int n_nodlump);


    /// @brief lump m_chi for a specific vector of node index
    /// @param this_vec: vector of int containing the node index
    /// @param n_nodlump: the lumping half_window (number of nodes)
    /// @author BG
    /// @date 08/01/2018
    void lump_this_vec(vector<int> this_vec, int n_nodlump);

    /// @brief Apply the TVD filter () L.Condat 2013 on the m_chi and b_chi signals
    /// @param this_vec: vector of int containing the node index
    /// @param n_nodlump: the lumping half_window (number of nodes)
    /// @author BG
    /// @date 08/01/2018
    void TVD_on_my_ksn(const float lambda, float lambda_TVD_b_chi);

    vector<float> TVD_this_vec(vector<int> this_vec, const float lambda, float lambda_TVD_b_chi);
    vector<double> correct_TVD_vec(vector<double> this_val);
    float get_dksn_from_composite_kp(vector<int> vecnode);
    float get_dseg_drop_from_composite_kp(vector<int> vecnode);
    float get_kp_sharpness_length(vector<int> vecnode, LSDFlowInfo& Flowinfo);
    int get_ksn_centroid_coordinates(LSDFlowInfo& Flowinfo, vector<int> vecnode,vector<int> vecnode_river, float total_dksn);
    vector<vector<int> > group_local_kp(vector<int> vecnode_kp, vector<int> vecnode_river,LSDFlowInfo& Flowinfo, int kp_node_search);
    vector<vector<int> > old_group_local_kp(vector<int> vecnode_kp, vector<int> vecnode_river,LSDFlowInfo& Flowinfo);
    void derive_the_segmented_elevation();
    vector<int> get_vecnode_river_from_extent(int first_node, int last_node, vector<int> vecnode_river);
    
    
    void print_final_ksn_knickpoint(LSDFlowInfo& FlowInfo, string filename);
    void raw_stepped_knickpoint_detection(int SK, vector<int> vecnode);
    void stepped_knickpoints_combining(LSDFlowInfo& Flowinfo, int kp_node_search);
    void TVD_this_vec_v2(vector<int> this_vec, float lambda, float lambda_TVD_b_chi, int max_node, string type);
    void stepped_knickpoints_detection_v2(LSDFlowInfo& Flowinfo, int window, float n_std_dev);
    map<string,vector<float> > get_windowed_stats_for_knickpoints(vector<int> vecnode,int HW);








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

    /// A general incrementer for knickpoints. It has to be global for some reason
    int id_kp;

    ///A map of strings for holding georeferencing information
    map<string,string> GeoReferencingStrings;

    // Some maps to store the data
    /// A map of the M_chi values. The indices are node numbers from FlowInfo
    map<int,float> M_chi_data_map;
    /// A map of the M_chi values. The indices are node numbers from FlowInfo
    map<int,float> b_chi_data_map;
    /// A map of the M_chi values. The indices are node numbers from FlowInfo
    map<int,float> elev_data_map;
    /// A map of the M_chi values. The indices are node numbers from FlowInfo
    map<int,float> chi_data_map;
    /// A map of the M_chi values. The indices are node numbers from FlowInfo
    map<int,float>  flow_distance_data_map;
    /// A map of the M_chi values. The indices are node numbers from FlowInfo
    map<int,float>  drainage_area_data_map;
    /// A map that holds elevations regressed from fitted sections.
    map<int,float> segmented_elevation_map;
    /// A map that holds segment numbers: used with skip = 0. Can be used to map
    /// distinct segments
    map<int,int> segment_counter_map;
    /// A map that holds knickpoints information
    map<int,float> segment_counter_knickpoint_map;
    /// A map that holds knickpoints signs
    map<int,int> segment_knickpoint_sign_map;
    /// A map that holds knickpoints signs
    map<int,int> segment_length_map;
    /// A map that holds knickpoints ratio
    map<int,float> ksn_ratio_knickpoint_map;
    /// A map that holds knickpoints difference_between_segments
    map<int,float> ksn_diff_knickpoint_map;
    /// A map that holds knickpoints signs
    map<int,int> ksn_sign_knickpoint_map;
    /// Map of the knickpoints value in radian
    map<int,float> ksn_rad_knickpoint_map;
    /// Map of the knickzone by cumulative variations
    map<int,float> ksn_cumul_knickpoint_map;

    /// Map of the knickpoint by cumulative variations -  ratio
    map<int,float> rksn_cumul_knickpoint_map;

    /// Map of the knickpoint by cumulative variations
    map<int,float> rad_cumul_knickpoint_map;

    /// map of nickzones for ksn variations
    map<pair<int,int>, float> knickzone_raw_cumul_ksn;
    /// map of nickzones for rksn variations
    map<pair<int,int>, float> knickzone_raw_cumul_rksn;
    /// map of nickzones for rad variations
    map<pair<int,int>, float> knickzone_raw_cumul_rad;
    /// map of nickzones for ksn variations WEIGHTED VERSION
    map<pair<int,int>, float> knickzone_WP_ksn;
    /// map of nickzones for rksn variations WEIGHTED VERSION
    map<pair<int,int>, float> knickzone_WP_rksn;
    /// map of nickzones for rad variations WEIGHTED VERSION
    map<pair<int,int>, float> knickzone_WP_rad;
    /// map of knickzone ID to identify all the knickzones from a same base one
    map<pair<int,int>,int>  knickzone_ID;
    /// map of source_keys of receiving rivers<source_key_of_river,source_key_of_receiving_river>
    map<int,int> map_source_key_receiver;
    /// map of m_chi precedant the source_key <source_key_of_river,previous_m_chi>
    map<int,float> map_source_key_receiver_mchi;

    /// map of source key and associated vector of nodes
    map<int,vector<int> > map_node_source_key;
        /// map of source key and associated vector of nodes
    map<int,vector<int> > map_source_key_vecnode_of_receiver;
    /// map of source key and associated vector of nodes containing knickpoints used for the ksnkp calculation
    map<int,vector<int> > map_node_source_key_kp;
    /// map of raw changes in ksn, key is node and value is delta ksn from bottom to top
    map<int,float> raw_ksn_kp_map;
    /// map of raw derivative for the ksn value calculated per rivers. map[nodeindex] = dksn/dchi
    map<int,float> raw_dksndchi_kp_map;
    /// map of raw KDE, calculated using method/binning depending on the parameter file
    map<int,float> raw_KDE_kp_map;
    /// map of the automatically calculated bandwidth per source key
    map<int,float> KDE_bandwidth_per_source_key;
    /// Map[source_key] = flow_length_of_river (from source to base junction, not baselevel)
    map<int,float> map_flow_length_source_key;
    /// Map[source_key] = chi_length_of_river (from source to base junction, not baselevel)
    map<int,float> map_chi_length_source_key;
    // Map[node_index] = 0 if not outlier, 1 if outlier according to a simple Modified z score on dksn/dchi
    map<int,int> map_outlier_MZS_dksndchi;
    /// Map[node_index] = lumped m_chi
    map<int,float> lumped_m_chi_map;
        /// Map[node_index] = TVDed m_chi
    map<int,float> TVD_m_chi_map;
    /// Map[node_index] = TVDed m_chi
    map<int,float> TVD_b_chi_map;
    /// Debugging map to check the TVD correctin (deprecated - I'll clean my code when I'll be sure I'll need it)
    map<int,float>TVD_m_chi_map_non_corrected;
    /// Grouped and processed knickpoints
    map<int,float> ksn_kp_map;
    /// diffusion of a knickpoint, or sharpness of a knickpoint, I cannot decide which word I prefer
    map<int,float> sharpness_ksn_length;
    /// diffusion of a knickpoint, or sharpness of a stepped knickpoint, I cannot decide which word I prefer
    map<int,float> sharpness_stepped_length;
    /// cextent of composite knickpoints
    map<int,pair<int,int> > ksn_extent;
    map<int,pair<int,int> > stepped_extent;
    /// contains the centroid of each knickpoints
    map<int,pair<float,float> > ksn_centroid;
    /// contains the centroid of each stepped knickpoints
    map<int,pair<float,float> > stepped_centroid;
    /// Unique ID for each knickpoints
    map<int,int> ksn_kp_ID;
    map<int,float> flow_distance_kp_centroid_map;
    map<int,float> flow_distance_stepped_kp_centroid_map;
    map<int,int> nearest_node_centroid_kp;
    map<int,int> nearest_node_centroid_kp_stepped;

    map<int,int> map_outlier_MZS_combined;
    map<int,float> kp_segdrop;
    map<int,float> raw_segchange;
    map<int,float> TVD_segelev_diff;
    map<int,float> segelev_diff;
    map<int,float> segelev_diff_second;
    map<int,float> raw_delta_segelev_from_TVDb_chi;
    map<int,vector<int> > map_node_source_key_kp_stepped;

    /// Map of intermediate values for each node during TVD segmentation
    map<int,float> intermediate_TVD_m_chi;
    map<int,float> intermediate_TVD_b_chi;
    map<int,float> mean_for_kp;
    map<int,float> std_for_kp;





    /// A vector to hold the order of the nodes. Starts from longest channel
    /// and then works through sources in descending order of channel lenght
    vector<int> node_sequence;

    /// vectors to hold the source nodes and the outlet nodes
    /// The source keys are indicies into the source_to_key_map.
    /// In big DEMs the node numbers become huge so for printing efficiency we
    /// run a key that starts at 0

    /// This map contains all the nodes. The key is the node index and the value
    ///  is the source key (sorry I know this is confusing). It means if you
    ///  have the node index you can look up the source key. Used for
    ///  visualisation.
    map<int,int> source_keys_map;

    /// This has all the nodes. The key (in the map) is the node index, and the
    ///  value is the baselevel key. Again used for visualisation
    map<int,int> baselevel_keys_map;

    /// This has as many elements as there are sources. The key in the map is the
    ///  node index of the source, and the value is the source key.
    map<int,int> key_to_source_map;

    /// This has as many elements as there are baselelvels. The key is the
    /// node index and the value is the baselevel key.
    map<int,int> key_to_baselevel_map;

    /// this is an ordered list of the source nodes (from first source to last)
    vector<int> ordered_source_nodes;

    /// this is an ordered list of baselelvel nodes (from first source to last)
    vector<int> ordered_baselevel_nodes;

    /// This vector contains the rank of each source node in each basin, so the
    /// main stem in each basin is 0, the second is 1, the 3rd is 2, etc. Counting starts
    /// again when a new baselevel node starts.
    vector<int> source_nodes_ranked_by_basin;


    /// In this map the key is the basin key and the value is the node of the
    /// baselevel source
    map<int,int> source_node_of_mainstem_map;


  private:
    void create(LSDRaster& Raster);
    void create(LSDIndexRaster& Raster);
    void create(LSDFlowInfo& FlowInfo);
    void create(LSDJunctionNetwork& JN);
};

#endif
