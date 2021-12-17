//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDJunctionNetwork
// Land Surface Dynamics ChannelNetwork
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for organizing channel routing under the Fastscape algorithm
//  (see Braun and Willett, Geomorphology 2013, v180, p 170-179)
//  It uses the algorithm to create channel junction networks
//  that can be searched for network connectivity
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

/** @file LSDJunctionNetwork.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 0.0.1
@brief Object to create a channel network from an LSDFlowInfo object.
@details This object is built around Braun and Willett's fastscape algorithm and
contains a number of analysis tools built around drainage networks.

@date 30/08/2012
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDJunctionNetwork_H
#define LSDJunctionNetwork_H

#include <vector>
#include <string>
#include <map>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
using namespace std;
using namespace TNT;

///@brief Object to create a channel network from an LSDFlowInfo object.
class LSDJunctionNetwork
{
  public:
  /// @brief This defines a channel network, is empty
  /// @author SMM
  /// @date 30/07/14
  LSDJunctionNetwork()   { create(); }

  /// @brief This defines a channel network based on a FlowInfo object and a list of source nodes.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param Sources vector of source nodes.
  /// @author SMM
  /// @date 01/09/12
  LSDJunctionNetwork(vector<int> Sources, LSDFlowInfo& FlowInfo)
                  { create(Sources, FlowInfo); }


  /// @brief Assignment operator.
  LSDJunctionNetwork& operator=(const LSDJunctionNetwork& LSDR);

  /// @brief this function gets the UTM_zone and a boolean that is true if
  /// the map is in the northern hemisphere
  /// @param UTM_zone the UTM zone. Replaced in function.
  /// @param is_North a boolean that is true if the DEM is in the northern hemisphere.
  ///  replaced in function
  /// @author SMM
  /// @date 22/12/2014
  void get_UTM_information(int& UTM_zone, bool& is_North);

  /// @brief this gets the x and y location of a node at row and column
  /// @param row the row of the node
  /// @param col the column of the node
  /// @param x_loc the x location (Northing) of the node
  /// @param y_loc the y location (Easting) of the node
  /// @author SMM
  /// @date 22/12/2014
  void get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc);

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
  /// @date 22/12/2014
  void get_lat_and_long_locations(int row, int col, double& lat,
                  double& longitude, LSDCoordinateConverterLLandUTM Converter);

  /// @brief This takes latitude and longitude (in WGS 84) and converts to vectors
  ///  of easting and northing in UTM
  /// @param latitude a vector of latitudes in UTM84
  /// @param longitude a vector of longitudes in WGS84
  /// @param UTME The easting coordinate (is overwritten)
  /// @param UTMN The northing coordinate (is overwritten)
  /// @author SMM
  /// @date 13/02/2017
  void get_x_and_y_from_latlong(vector<float> latitude, vector<float> longitude,
                                                   vector<float>& UTME,vector<float>& UTMN);

  ///@brief Recursive add_to_stack routine to build the junction tree, from Braun and Willett (2012)
  ///equations 12 and 13.
  ///@param lm_index Integer
  ///@param j_index Integer
  ///@param bl_node Integer
  /// @author SMM
  /// @date 01/09/12
  void add_to_stack(int lm_index, int& j_index, int bl_node);

  // this returns all the upstream junction of a junction_number_outlet
  /// @brief This returns all the upstream junction of a junction_number_outlet.
  /// @param junction_number_outlet Integer of junction of interest.
  /// @return integer vector containing all the junction numbers upslope
  /// of the chosen junction.
  /// @author SMM
  /// @date 01/09/12
  vector<int> get_upslope_junctions(int junction_number_outlet);

  // this returns all the upstream junction of a given order junction_number_outlet
  /// @brief This returns all the upstream junction of a junction_number_outlet.
  /// @param junction_number_outlet Integer of junction of interest.
  /// @param stream order
  /// @return integer vector containing all the junction numbers upslope
  /// of the chosen junction.
  /// @author FJC
  /// @date 13/08/19
  vector<int> get_upslope_junctions_by_order(int junction_number_outlet, int stream_order);

  /// @brief This finds all the junctions that are source junctions upslope of a
  ///  given junction
  /// @param junction_number_outlet The junction number of the outlet
  /// @return source_junctions a vector of junction numbers: these are the sources
  /// @author SMM
  /// @date 18/05/2016
  vector<int> get_all_source_junctions_of_an_outlet_junction(int junction_number_outlet);

  /// @brief This finds all the nodes that are source nodes upslope of a
  ///  given junction
  /// @param junction_number_outlet The junction number of the outlet
  /// @return source_nodes a vector of node numbers: these are the sources
  ///   the nodes are node indices from the FlowInof object
  /// @author SMM
  /// @date 19/05/2016
  vector<int> get_all_source_nodes_of_an_outlet_junction(int junction_number_outlet);

  /// @brief this function gets a list of the junction indices of the donors to a particular junction
  /// @detail IMPORTANT: this has only retained the string "node" to keep equivalence
  ///  with the FlowInfo object. It takes junctions and returns junctions!!
  ///  Also note that base level nodes have themselves as a donor
  /// @param node this is the nodeindex of the node for which you want to find the donors
  /// @return a vector of the donor nodes
  /// @author SMM
  /// @date 15/06/2015
  vector<int> get_donor_nodes(int node);

  /// @brief This function maps a junction onto the indexing of the upslope junction list.
  ///
  /// @details If you get an upslope junction list the indexing starts at the furthest downslope
  /// junction. All of the junction pointing refers to the master junction list however.
  /// @param upslope_junctions Vector of upslope junctions of interest.
  /// @param junction Integer of junction of interest.
  /// @return Vector of mapped junctions.
  /// @author SMM
  /// @date 01/09/12
  int map_junction_to_upslope_junction_list(vector<int> upslope_junctions, int junction);

  // functions for finding specific basins
  /// @brief This function returns the maximum stream order in the DEM.
  /// @return Maximum stream order as an integer.
  /// @author SMM
  /// @date 01/09/12
  int get_maximum_stream_order();

  /// @brief This function returns the number of streams of a given stream order
  /// @param FlowInfo LSDFlowInfo object
  /// @param stream_order Stream order of interest
  /// @return integer with number of streams.
  /// @author FJC
  /// @date 15/03/16
  int get_number_of_streams(LSDFlowInfo& FlowInfo, int stream_order);

  /// @brief This calculates the junction angles based on a number of junctions
  /// @param JunctionList a list of junctions
  /// @param FlowInfo an LSDFlowInfo object
  /// @return A vector of junction angles
  /// @author SMM
  /// @date 21/04/2017
  map<int, vector<float> > calculate_junction_angles(vector<int> JunctionList, LSDFlowInfo& FlowInfo);

  /// @brief This is a more complete function for junction angles
  ///  It overwirtes two maps, containing all sorts of information about
  ///   the junction angles.
  /// @param JunctionList a list of junctions
  /// @param FlowInfo an LSDFlowInfo object
  /// @param JA_int_info A map where the key is the junction number
  ///    and the vector is a series of integer data about that juction
  /// @param JA_float_info A map where the key is the junction number
  ///    and the vector is a series of float data about that juction
  /// @author SMM
  /// @date 17/11/2019
  void calculate_junction_angles_complete(vector<int> JunctionList,
                                          LSDFlowInfo& FlowInfo,
                                          map<int , vector<int> >& JA_int_info,
                                          map<int, vector<float> >& JA_float_info );


  /// @brief This is a more complete function for junction angles
  ///  It overwirtes two maps, containing all sorts of information about
  ///   the junction angles.
  ///  This is an overloaded version that gives some slope and flow distance information
  /// @param JunctionList a list of junctions
  /// @param FlowInfo an LSDFlowInfo object
  /// @param JA_int_info A map where the key is the junction number
  ///    and the vector is a series of integer data about that juction
  /// @param JA_float_info A map where the key is the junction number
  ///    and the vector is a series of float data about that juction
  /// @author SMM
  /// @date 17/11/2019
  void calculate_junction_angles_complete(vector<int> JunctionList,
                                          LSDFlowInfo& FlowInfo, LSDRaster& Elevation,
                                          LSDRaster& FlowDistance, float vertical_interval,
                                          map<int , vector<int> >& JA_int_info,
                                          map<int, vector<float> >& JA_float_info );


  /// @brief This function gets the mean and standard error of every junction angle
  ///   upslope of a given junction
  /// @param target_junction The target junction
  /// @param FlowInfo an LSDFlowInfo object
  /// @return A vector of that has the mean and the standard error of the upslope junction angles
  /// @author SMM
  /// @date 23/04/2017
  vector<float> calculate_junction_angle_statistics_upstream_of_junction(int target_junction, LSDFlowInfo& FlowInfo);

  /// @brief Overloaded function similar to above but removes any junctions not greater than
  /// threshold SO
  /// @param target_junction The target junction
  /// @param FlowInfo an LSDFlowInfo object
  /// @param threshold_SO threshold stream order to keep junctions (greater than this)
  /// @return A vector of that has the stats of the upslope junction angles
  /// @author FJC
  /// @date 08/03/18
  vector<float> calculate_junction_angle_statistics_upstream_of_junction(int target_junction, LSDFlowInfo& FlowInfo, int threshold_SO);

  /// @brief This takes the junction angle statistics for all basins of a given order
  /// @param FlowInfo the LSDFlowInfo object
  /// @param BasinOrder the basin order of interest
  /// @param junction_list a vector of ints holding the junctions of interest
  ///  is replaced in the function
  /// @param junction_angle_averages Average junction angles
  ///  is replaced in the function
  /// @param junction_angle_stder a vector junction angle standard errors
  ///  is replaced in the function
  /// @param N_junctions a vector of ints holding the numer of junctions in each larger basin
  ///  is replaced in the function
  /// @author SMM
  /// @date 24/04/2017
  void calculate_junction_angle_statistics_for_order(LSDFlowInfo& FlowInfo, int BasinOrder,
                             vector<int>& junction_list,
                             vector<float>& junction_angle_averages,
                             vector<float>& junction_angle_stderr,
                             vector<int>& N_junctions);

  /// @brief This function takes a vector of basin junctions and prints statistics of all the junctions
  /// upstream of each basin junction to a CSV.  The statstics are separated by stream order.
  /// @param JunctionList list of basin junctions
  /// @param FlowInfo LSDFlowInfo object
  /// @param csv_outname name of output csv
  /// @author FJC
  /// @date 08/03/18
  void print_junction_angles_from_basin_list(vector<int> JunctionList, LSDFlowInfo& FlowInfo, string csv_outname);


  /// @brief This prints the junction angles to a csv file
  /// @param JunctionList The list of junctions to analyze. If this is an empty vector,
  ///  the code analyses all junctions in the DEM
  /// @param FlowInfo The LSDFlowInfo object
  /// @param csv_name The name of the file. Needs full path and csv extension
  /// @author SMM
  /// @date 23/04/2017
  void print_junction_angles_to_csv(vector<int> JunctionList, LSDFlowInfo& FlowInfo,
                                                       string csv_name);

  /// @brief This prints the junction angles to a csv file
  ///  It uses the much more complete junction angle code
  /// @param JunctionList The list of junctions to analyze. If this is an empty vector,
  ///  the code analyses all junctions in the DEM
  /// @param FlowInfo The LSDFlowInfo object
  /// @param csv_name The name of the file. Needs full path and csv extension
  /// @author SMM
  /// @date 17/11/2019
  void print_complete_junction_angles_to_csv(vector<int> JunctionList, LSDFlowInfo& FlowInfo,
                                                       string csv_name);

  /// @brief This prints the junction angles to a csv file
  ///  It uses the much more complete junction angle code. In addition it give slope and elevation
  ///  data for the junctions
  /// @param JunctionList The list of junctions to analyze. If this is an empty vector,
  ///  the code analyses all junctions in the DEM
  /// @param FlowInfo The LSDFlowInfo object
  /// @param Elevations An elevation raster
  /// @param FlowDistance A flow distance raster
  /// @param vertical_interval: the vertical interval around which you want the slope measured
  /// @param csv_name The name of the file. Needs full path and csv extension
  /// @author SMM
  /// @date 07/12/2019
  void print_complete_junction_angles_to_csv(vector<int> JunctionList,
                                             LSDFlowInfo& FlowInfo, LSDRaster& Elevations,
                                             LSDRaster& FlowDistance, float vertical_interval,
                                             string csv_name);

  /// @brief This gets the junction number of a given node.
  /// @param Node
  /// @param FlowInfo Flow Info object
  /// @return JunctionNumber
  /// @author FC
  /// @date 31/10/13
  int get_Junction_of_Node(int Node, LSDFlowInfo& FlowInfo);

  /// @brief This gets the junction number all the sources
  /// @param FlowInfo Flow Info object
  /// @return A vector of junctions from all the sources
  /// @author SMM
  /// @date 08/05/15
  vector<int> get_Junctions_of_Sources(LSDFlowInfo& FlowInfo);

  /// @brief returns the penultimate node of the stream link below given junction
  /// @param upstream junction of desired stream link
  /// @param FlowInfo object
  /// @return node index (for FlowInfo) of penultimate node in stream link
  /// @author DTM
  /// @date 04/06/14
  int get_penultimate_node_from_stream_link(int upstream_junction, LSDFlowInfo& FlowInfo);

  // this prints the link array to raster
  /// @brief This sends the StreamOrderArray to a LSDIndexRaster.
  /// @return LSDIndexRaster of StreamOrderArray.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexRaster StreamOrderArray_to_LSDIndexRaster();

  /// @brief Method to flatten a stream order array and place the non NDV values in a csv file.
  /// @detail Each value is placed on its own line, so that it can be read more quickly in python etc.
  ///   It includes the lat long coordinates in CSV, in WGS84 coordinate system EPSG:4326
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @author SMM
  /// @date 12/11/16
  void StreamOrderArray_to_WGS84CSV(string FileName);


  /// @brief This prints a stream network to a csv in WGS84
  /// @detail This function prints a network that is ordered by sources, channels
  ///  have stream orders and junction numbers attached
  /// @param FlowInfo the flow info object which translates node indices to actual points
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @author SMM
  /// @date 14/11/16
  void PrintChannelNetworkToCSV(LSDFlowInfo& flowinfo, string fname_prefix);

  /// @brief This prints a stream network to a csv in XY coordinates
  /// @param FlowInfo the flow info object which translates node indices to actual points
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @param elevation The elevation raster
  /// @author BG
  /// @date 28/01/2019
  void PrintChannelNetworkToCSV_nolatlon(LSDFlowInfo& flowinfo, LSDRaster& elevation, string fname_prefix);

  /// @brief This prints a stream network to a csv in WGS84. It includes the elevation data
  /// @param FlowInfo the flow info object which translates node indices to actual points
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @param elevation The elevation raster
  /// @author SMM
  /// @date 04/03/2020
  void PrintChannelNetworkToCSV_WithElevation(LSDFlowInfo& flowinfo, string fname_prefix,LSDRaster& elevation);

  /// @brief This prints a stream network to a csv in WGS84. It includes the elevation data and the donor junction
  /// @param FlowInfo the flow info object which translates node indices to actual points
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @param elevation The elevation raster
  /// @author SMM
  /// @date 08/11/2020
  void PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(LSDFlowInfo& flowinfo, string fname_prefix,LSDRaster& elevation);


  /// @brief This sends the JunctionArray to a LSDIndexRaster.
  /// @return LSDIndexRaster of JunctionArray.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexRaster JunctionArray_to_LSDIndexRaster();

  /// @brief This sends the JunctionIndexArray to a LSDIndexRaster.
  /// @return LSDIndexRaster of JunctionIndexArray.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexRaster JunctionIndexArray_to_LSDIndexRaster();

  /// @brief Turns the StreamOrderArray into a binary rastser where 1 is channel and 0 is hillslope.
  /// @return Binary LSDIndexRaster of the channel network.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexRaster StreamOrderArray_to_BinaryNetwork_LSDIndexRaster();


  /// @brief This gets the largest donor junction to the baselevel nodes so that you can
  /// automate basin selection. (e.g. for use with chi analysis)
  ///
  /// @details This function returns a integer vector containing the junction number of the largest
  /// donor catchment (i.e. donor junction with greatest drainage area) upslope of each
  /// baselevel node. These can then be used as the starting locations for performing chi
  /// analysis.
  ///
  /// IMPORTANT: the junctions always point downstream since they can have one and only
  /// one receiver. However, for a basin of given order, this starts just upstream of the
  /// confluence to the next basin order. So the basin <b>INCLUDES</b> the channel flowing
  /// downstream to the penultamite node.
  /// @return Integer vector containing the junction number of the largest donor catchment.
  /// @author MDH
  /// @date 19/6/13
  vector<int> get_BaseLevel_DonorJunctions();

  /// @brief This function takes a list of junctions and then prunes
  ///  junctions based whether they drain from the edge. This attempts to
  ///  remove junctions that are through-flowing and thus do not have the
  ///  correct drainage area
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 16/05/16
  vector<int> Prune_Junctions_Edge(vector<int>& BaseLevelJunctions_Initial,LSDFlowInfo& FlowInfo);

  /// @brief This function takes a list of junctions and then prunes
  ///  junctions based whether they drain from the edge. This attempts to
  ///  remove junctions that are through-flowing and thus do not have the
  ///  correct drainage area
  /// @detail Only gets the donor of the baselelve donor to ignore the nodes
  ///  near the outlet, which often intersect nodata in cut DEMs
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param TestRaster A raster that is just used to look for nodata
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 29/05/17
  vector<int> Prune_Junctions_Edge_Ignore_Outlet_Reach(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& TestRaster);

  /// @brief This function looks through all baselevel nodes and then
  ///  looks for the largest basin that is not influenced by the edge.
  ///  It returns a vector of these junctions.
  /// @detail Note that it only returns one basin per baselevel node at most
  ///  so might not do a great job of space filling.
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param TestRaster A raster that is just used to look for nodata
  /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 21/06/17
  vector<int> Prune_To_Largest_Complete_Basins(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& TestRaster,
                                              LSDIndexRaster& FlowAcc);


    /// @brief This function removes basins that fall outside a contributing pixel
    ///  Window
    /// @param Junctions_Initial a vector of integers containg an inital
    ///  list of junctions
    /// @param FlowInfo The LSDFlowInfo object
    /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
    /// @param lower_limit The minimum number of contributing pixels
    /// @param upper_limit The maximum number of contributing pixels
    /// @return a pruned list of base level nodes
    /// @author SMM
    /// @date 26/06/17
    vector<int> Prune_Junctions_By_Contributing_Pixel_Window(vector<int>& Junctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc,
                                              int lower_limit, int upper_limit);


    /// @brief This function removes basins that fall outside a contributing pixel
    ///  Window, those that are bounded by nodata, and those that
    ///  are nested. A rather intensive pruning process that hopeuflly results
    ///  in a number of basins that are a similar size
    /// @detail This doesn't just look for baselevel junctions: it goes through
    ///  all junctions in the DEM. Warning: computationally expensive!
    /// @param FlowInfo The LSDFlowInfo object
    /// @param TestRaster A raster that is just used to look for nodata
    /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
    /// @param lower_limit The minimum number of contributing pixels
    /// @param upper_limit The maximum number of contributing pixels
    /// @return a pruned list of base level nodes
    /// @author SMM
    /// @date 26/06/17
    vector<int>  Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(LSDFlowInfo& FlowInfo,
                                              LSDRaster& TestRaster, LSDIndexRaster& FlowAcc,
                                              int lower_limit, int upper_limit);

    /// @brief This function removes basins that are nested within any other
    ///  basin in the list
    /// @param Junctions_Initial a vector of integers containg an inital
    ///  list of junctions
    /// @param FlowInfo The LSDFlowInfo object
    /// @param TestRaster A raster that is just used to look for nodata
    /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
    /// @return a pruned list of base level nodes
    /// @author SMM
    /// @date 26/06/17
    vector<int> Prune_Junctions_If_Nested(vector<int>& Junctions_Initial,
                                      LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc);

  /// @brief This function takes a list of junctions and then prunes
  ///  junctions based on their number of contributing pixels
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
  /// @param Threshold The minimum number of accumulated pixels needed to keep
  ///   a base level node.
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 16/05/16
  vector<int> Prune_Junctions_Area(vector<int>& BaseLevelJunctions_Initial,LSDFlowInfo& FlowInfo,
                              LSDIndexRaster& FlowAcc, int Threshold);

  /// @brief This function takes a list of junctions retains ONLY the larges bains
  ///  The junction is returned as an int vector so that it can be passed to other functions
  ///  requiring junction lists.
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 03/06/16
  vector<int> Prune_Junctions_Largest(vector<int>& BaseLevelJunctions_Initial,LSDFlowInfo& FlowInfo,
                              LSDIndexRaster& FlowAcc);

  /// @brief This function takes a list of junctions retains ONLY the junctions
  ///  that have an outlet elevation greater or less than the threshold elevation
  ///  Selection of greater or lower is determined by bool keep_junctions_below_threshold
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param Elev an LSDRaster of elevation
  /// @param threshold_elevation the threshold elevation to kepp
  /// @param keep_junctions_below_threshold if true keep junctions below threshold
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 18/01/18
  vector<int> Prune_Junctions_Threshold_Elevation(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& Elev,
                                              float threshold_elevation, bool keep_junctions_below_threshold);

  /// @brief This function takes a list of junctions retains ONLY the junctions
  ///  that have an outlet elevation with an elevation window
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param Elev an LSDRaster of elevation
  /// @param lower_threshold the lower threshold elevation
  /// @param upper_threshold the lower threshold elevation
  /// @return a pruned list of base level nodes
  /// @author SMM
  /// @date 19/01/18
  vector<int> Prune_Junctions_Elevation_Window(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& Elev,
                                              float lower_threshold, float upper_threshold);


  /// @brief You give this a list of junction numbers and it returns the
  ///  number of upslope pixels
  /// @param BaseLevelJunctions_Initial a vector of integers containg an inital
  ///  list of base level nodes
  /// @param FlowInfo The LSDFlowInfo object
  /// @param FlowAcc an LSDIndexRaster with the number of pixels for flow accumulation
  /// @param Threshold The minimum number of accumulated pixels needed to keep
  ///   a base level node.
  /// @return a vector with the N contributing pixels for the junctions specified
  /// @author SMM
  /// @date 21/06/17
  vector<int> get_contributing_pixels_from_specified_junctions(vector<int>& JunctionList,
                                              LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc);

  /// @brief Get Junction number at a location.
  /// @param row Integer row index.
  /// @param col Integer column index.
  /// @return Junction number at location row,col.
  /// @author SMM
  /// @date 01/09/12
  int retrieve_junction_number_at_row_and_column(int row,int col)
                       { return JunctionIndexArray[ row ][ col ]; }

  /// @brief Function for printing out the longest channel upstream of a point.
  /// @param outlet_junction
  /// @param FInfo LSDFlowInfo object.
  /// @param dist_code
  /// @param dist_from_outlet
  /// @author SMM
  /// @date 01/09/12
  void print_longest_channel(int outlet_junction, LSDFlowInfo& FInfo, LSDIndexRaster& dist_code,
                             LSDRaster& dist_from_outlet);

  /// @brief Prints the information about the junctions to file.
  /// @param filename Output filename to be appended with '.txt'.
  /// @author SMM
  /// @date 01/09/12
  void print_junction_info_vectors(string filename);

  /// @brief This generates an LSDChannelIndex object given a junction.
  ///
  /// @details NOTE: junctions start at the upstream end of the channel section.
  /// @param start_junction Junction to extract the channel from.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return The channel for the given junction.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexChannel generate_link_index_channel_from_junction(int start_junction,LSDFlowInfo& FlowInfo);

  /// @brief This function extracts the longest channel originating from a junction number
  /// outlet_junction.
  /// @param outlet_junction Outlet of junction.
  /// @param FInfo LSDFlowInfo object.
  /// @param dist_from_outlet Distance from outlet junction.
  /// @return LSDIndexRaster of the longest channel.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexChannel generate_longest_index_channel_from_junction(int outlet_junction,LSDFlowInfo& FInfo,LSDRaster& dist_from_outlet);

  // this generates the longest channel in a basin. The basin starts where a channel of
  // some order intersects with a channel of higher order. So the bain includes the
  // basin junction, but also the channel flowing downstream from this basin
  // junction
  // It starts from the node of the receiver junction, so if one were to extract
  // the basin from this node one would get a basin that starts one node upstream from
  // the lowest node in this
  /// @brief This generates the longest channel in a basin.
  ///
  /// @details The basin starts where a channel of some order intersects with a
  /// channel of higher order. So the bain includes the basin junction, but also
  /// the channel flowing downstream from this basin junction. It starts from the
  /// node of the receiver junction, so if one were to extract the basin from
  /// this node one would get a basin that starts one node upstream from the lowest node in this.
  /// @param basin_junction
  /// @param FInfo LSDFlowInfo object.
  /// @param dist_from_outlet
  /// @return LSDIndexRaster of the longest channel.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexChannel generate_longest_index_channel_in_basin(int basin_junction, LSDFlowInfo& FInfo,
            LSDRaster& dist_from_outlet);

  /// @brief This generates the upstream source nodes from a vector of basin junctions
  ///
  /// @details The basin starts where a channel of some order intersects with a
  /// channel of higher order. So the bain includes the basin junction, but also
  /// the channel flowing downstream from this basin junction. It starts from the
  /// node of the receiver junction, so if one were to extract the basin from
  /// this node one would get a basin that starts one node upstream from the lowest node in this.
  /// @param basin_junction
  /// @param FInfo LSDFlowInfo object.
  /// @param dist_from_outlet
  /// @return LSDIndexRaster of the longest channel.
  /// @author FJC
  /// @date 21/03/17
  vector<int> get_basin_sources_from_outlet_vector(vector<int> basin_junctions, LSDFlowInfo& FlowInfo,
                             LSDRaster& dist_from_outlet);

  /// @brief This extracts the junction numbers, in a vector of integers, of all basins of a
  /// given order.
  ///
  /// @details For basins, the basin includes nodes downstream of the basinJunction,
  /// until the penulatmite node in this downstream channel.
  ///
  /// IMPORTANT: the junctions always point downstream since they can have one and only
  /// one receiver. However, for a basin of given order, this starts just upstream of the
  /// confluence to the next basin order. So the basin <b>INCLUDES</b> the channel flowing
  /// downstream to the penultamite node.
  /// @param BasinOrder Integer of the basin order.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Vector of junctions of basins of given order.
  /// @author SMM
  /// @date 01/09/12
  vector<int> extract_basins_order_outlet_junctions(int BasinOrder, LSDFlowInfo& FlowInfo);

  // this function gets the outlet node of a list of basins, contined within the
  // BasinOutletJunctions parameter. The basin outlet node is _DOWNSTREAM_ from
  // the outlet junction, it is the penultamite node of the channel index.

  /// @brief this function gets the outlet node of a list of basins.
  ///
  /// @details The basin outlet node is _DOWNSTREAM_ from the outlet junction,
  /// it is the penultamite node of the channel index.
  /// @param BasinOutletJunctions
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Vector of outlet nodes of basins.
  /// @author SMM
  /// @date 01/09/12
  vector<int> extract_basins_order_outlet_nodes(vector<int>& BasinOutletJunctions, LSDFlowInfo& FlowInfo);

  /// @brief This function gets tributaries along a continous channel.
  ///
  /// @details What it does is goes down the index channel looking at the JunctionIndexArray
  /// to see if there is a junction. If it hits a junction then all the contributing junction
  /// it overwrites two vectors: \n
  /// tributary_junctions, which lists all junctions whose receiver is the main
  /// stem and nodes_on_main_stem_of_tributaries, which are the njodes on the
  /// main_stem LSDIndexChannel where the tributaries intersect the main stem
  /// this second vector is used to calcualte the chi values of the downstream node
  /// of the tributaries.
  /// @param MainStem LSDIndexChannel of the main stem.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param tributary_junctions
  /// @param nodes_on_main_stem_of_tributaries
  /// @author SMM
  /// @date 01/09/12
  void extract_tributary_junctions_to_main_stem(LSDIndexChannel& MainStem, LSDFlowInfo& FlowInfo,
                        vector<int>& tributary_junctions,
                        vector<int>& nodes_on_main_stem_of_tributaries);

  /// @brief this function gets the tributary junctions upstream of the starting_junction based on
  /// pruning criteria.
  ///
  /// @details This function extracts tributaries juncions to the main stem of the
  /// channel, then selects a sample based on various criteria set by an integer
  /// called pruning switch \n\n
  /// pruning_switch == 0  channels are only added if they exceed a threshold drainage area \n
  /// pruning_switch == 1  channels are only added if the ratio between them
  ///   and the mainstem exceeds a certain value (pruning_threshold)\n
  /// pruning_switch == 2  channels are only added if the ratio between them
  ///   and the area of the  mainstem _at the junction_ exceeds a certain value\n
  /// pruning_switch == 3 channels are only added if the channel order is >= threshold.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param ChannelNetwork LSDJunctionNetwork object.
  /// @param starting_junction
  /// @param DistanceFromOutlet LSDIndexRaster of outlet distances.
  /// @param pruning_switch
  /// @param pruning_threshold
  /// @return Pruned tributary junctions.
  /// @author DTM
  /// @date 30/04/2013
  vector<int> get_pruned_tributaries_from_main_stem(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChannelNetwork,
                        int starting_junction, LSDRaster& DistanceFromOutlet,
                        int pruning_switch, float pruning_threshold);

  /// @brief This function extracts basin nodes according to their accumulated drainage area.
  /// @param Threshold Threshold drainage area.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Vector of basin nodes.
  /// @author DTM
  /// @date 07/05/2013
  vector<int> extract_basin_nodes_by_drainage_area(float DrainageAreaThreshold, LSDFlowInfo& FlowInfo);

	/// @brief This function extracts nodes where the basins of both tributaries are greater
	/// than a certain drainage area threshold.  Moves downstream from sources to baselevel so that
	/// nested catchments will be selected
	/// @param FlowInfo LSDFlowInfo object.
  /// @param DrainageAreaThreshold Threshold drainage area.
  /// @return Vector of basin nodes. These are the nodes just upstream of the outlet junction at
	/// the confluence of the basins.
  /// @author FJC
  /// @date 10/01/17
	vector<int> extract_basin_nodes_above_drainage_area_threshold(LSDFlowInfo& FlowInfo, float DrainageAreaThreshold);

  /// @brief This function checks all of the basin nodes to check if they fall within a mask
  /// (input raster). If they fall within the mask raster then the first node upstream
  /// not in the mask is selected.
	/// @param basin_nodes vector of basin nodes
  /// @param FlowInfo LSDFlowInfo object
  /// @param MaskRaster raster to use as mask
  /// @return vector with the modified basin nodes
  /// @author FJC
  /// @date 31/01/17
  vector<int> modify_basin_nodes_from_mask(vector<int> basin_nodes, LSDFlowInfo& FlowInfo, LSDRaster& MaskRaster);

  /// @brief This function extracts basin junctions from a list of basin outlet nodes.
  /// @param basin_nodes list of basin outlet nodes
  /// @param FlowInfo LSDFlowInfo object
  /// @return vector of basin junctions
  /// @author FJC
  /// @date 15/01/2014
  vector<int> extract_basin_junctions_from_nodes(vector<int> basin_nodes, LSDFlowInfo& FlowInfo);

  /// @brief This function gets the node indices of outlets of basins of a certain order
  ///
  /// @details IMPORTANT: The junctions always point downstream since they can have one and only
  /// one receiver. However, for a basin of given order, this starts just upstream of the
  /// confluence to the next basin order. So the basin <b>INCLUDES</b> the channel flowing
  /// downstream to the penultamite node.
  ///
  /// @param basin_junction Junction of basin to be extracted.
  /// @param basin_reference_number Reference number for printing to the IndexRaster.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of extracted basin.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexRaster extract_basin_from_junction(int basin_junction, int basin_reference_number, LSDFlowInfo& FlowInfo);

  /// @brief This function converts a list of sources used to generate the initial channel network
  /// into a list of junction indexes of channel heads which can be used to extract hollows.
  /// @param Sources A vector of source nodes that correspond to channel heads.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return A vector of junction indexes for each channel head.
  /// @author SWDG
  /// @date 05/12/13
  vector<int> Get_Channel_Head_Junctions(vector<int> Sources, LSDFlowInfo& FlowInfo);

  /// @brief This function extracts a single hollow from a given channel head junction.
  ///
  /// @details The junction index of channel heads can be extracted using LSDJunctionNetwork.Get_Channel_Head_Junctions.
  /// @param CH_junction Junction index to extract.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of the extracted hollow, coded with junction number.
  /// @author SWDG
  /// @date 05/12/13
  LSDIndexRaster extract_hollow(int CH_junction, LSDFlowInfo& FlowInfo);

  /// @brief This function extracts a series of hollows from a vector of channel head junctions.
  ///
  /// @details The junction index of channel heads can be extracted using LSDJunctionNetwork.Get_Channel_Head_Junctions.
  /// @param CH_junctions Vector of juntions to extract.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of the extracted hollows, coded with junction numbers.
  /// @author SWDG
  /// @date 05/12/13
  LSDIndexRaster extract_hollow(vector<int> CH_junctions, LSDFlowInfo& FlowInfo);

  /// @brief This function gets the an LSDIndexRaster of basins draining from a vector of junctions.
  ///
  /// @details IMPORTANT: The junctions always point downstream since they can have one and only
  /// one receiver. However, for a basin of given order, this starts just upstream of the
  /// confluence to the next basin order. So the basin <b>INCLUDES</b> the channel flowing
  /// downstream to the penultamite node.
  ///
  /// @param basin_junctions Vector of junction numbers of basins to be extracted.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of extracted basin.
  /// @author SMM
  /// @date 01/09/12
  LSDIndexRaster extract_basins_from_junction_vector(vector<int> basin_junctions, LSDFlowInfo& FlowInfo);


  /// @brief This function gets an LSDIndexRaster of basins draining from a vector of junctions.
  ///
  /// @details IMPORTANT: The junctions always point downstream since they can have one and only
  /// one receiver. However, for a basin of given order, this starts just upstream of the
  /// confluence to the next basin order. So the basin <b>INCLUDES</b> the channel flowing
  /// downstream to the penultamite node.
	/// UPDATED so that if basins are nested, they don't overwrite each other - basins are
	/// sorted by the number of contributing pixels, and the smaller basins are written
	/// first.
  ///
  /// @param basin_junctions Vector of junction numbers of basins to be extracted.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of extracted basin.
  /// @author FJC
  /// @date 10/01/17
	LSDIndexRaster extract_basins_from_junction_vector_nested(vector<int> basin_junctions, LSDFlowInfo& FlowInfo);

  /// @brief This function gets the an LSDIndexRaster of basins draining from a vector of junctions.
  /// @details IThis is a highly rudimentary version, which just collects
  ///  all the upslope nodes.
  /// @param basin_junctions Vector of junction numbers of basins to be extracted.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of extracted basin.
  /// @author SMM
  /// @date 08/05/15
  LSDIndexRaster extract_basins_from_junctions_rudimentary(vector<int> junctions, LSDFlowInfo& FlowInfo);

  /// @brief Basin extraction - extracts all drainage basins of specified stream order.
  /// @param BasinOrder Integer basin order to extract.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDIndexRaster of extracted basins.
  /// @author DTM
  /// @date 17/10/2012
  LSDIndexRaster ExtractBasinsOrder(int BasinOrder, LSDFlowInfo& FlowInfo);

  /// @brief This function extracts the juctions of all non-beheaded drainage basins of a given order, n.
  /// @param BasinOrder Integer basin order to extract.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Vector of junction indexes.
  /// @author SWDG
  /// @date 24/10/2013
  vector<int> ExtractBasinJunctionOrder(int BasinOrder, LSDFlowInfo& FlowInfo);

  /// @brief This function extracts the juctions of all non-beheaded drainage basins of a given order, n.
  ///  Like the previous version but in this case includes basins at the edge (abutting nodata)
  /// @param BasinOrder Integer basin order to extract.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Vector of junction indexes.
  /// @author SMM
  /// @date 29/04/2017
  vector<int> ExtractBasinJunctionOrderKeepEdgeBasins(int BasinOrder, LSDFlowInfo& FlowInfo);

  /// @brief Get farthest upslope hilltops.
  ///
  /// @details This function looks at all the source junctions in a network
  ///  upstream of a given junction and returns the node index of the
  ///  hilltop node that is the farthest upstream from the source junction
  /// @param JunctionNumber the junction number upstream of which you want to search for sources
  /// @param FlowInfo the flow info object
  /// @param FlowDistance distance upslope
  /// @return vector<int> a vector of node indices to the ridge nodes that are farthest upslope
  /// of the sources
  /// @author SMM
  /// @date 26/09/2013
  vector<int> FindFarthestUpslopeHilltopsFromSources(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance);

  /// @brief This function generates LSDChannels that run from the hilltops above
  /// all the sources of the junction JunctionNumber
  /// @author SMM
  /// @date 26/09/2013
  int GetChannelHeadsChiMethodFromNode(int NodeNumber,
                              int MinSegLength, float A_0, float m_over_n,
            LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& ElevationRaster);

  /// @brief This function generates LSDChannels that run from the hilltops above
  /// all the sources from the valley network down to a specified number of downstream junctions below
  /// the sources
  /// @author FJC
  /// @date 10/09/15
  int GetChannelHeadsChiMethodFromSourceNode(int NodeNumber,
                        int MinSegLength, float A_0, float m_over_n,
                        LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& ElevationRaster, int NJunctions);

  /// @brief This function generates LSDChannels that run from the hilltops above
  /// all the sources from the valley network down to a specified number of downstream junctions below
  /// the sources and writes the profile to csv
  /// @author FJC
  /// @date 23/12/16
	void write_valley_hilltop_chi_profiles_to_csv(vector<int> sources, float A_0, float m_over_n, LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& ElevationRaster, int NJunctions, string output_path, string DEM_ID);

  /// @brief This function generates an LSDIndexRaster of the channel that runs from
  /// the hilltop above the furthest upslope source of the junction JunctionNumber
  /// @param BasinOrder
  /// @param MinSegLength
  /// @param A_0
  /// @param m_over_n
  /// @param FlowInfo
  /// @param FlowDistance
  /// @param ElevationRaster
  /// @return LSDIndexRaster with channel
  /// @author FJC
  /// @date 21/08/15
  LSDIndexRaster GetChannelfromDreich(int NodeNumber, int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& ElevationRaster, string path_name, int NJunctions);



  /// @brief This function returns all potential channel heads in a DEM. It looks for
  /// channel heads organized by a basin order which is fed to the code
  /// The basin order just determines how far downstream the algorithm looks for the 'fluvial'
  /// section.
  /// It returns a vector<int> of nodeindices where the channel heads are
  /// @return vector<int> a vector of node_indices of potential channel heads
  /// @author SMM
  /// @date 26/09/2013
//   vector<int> GetChannelHeadsChiMethodBasinOrder(int BasinOrder,
//                        int MinSegLength, float A_0, float m_over_n,
//            LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
//            LSDRaster& ElevationRaster);

  /// @brief This function returns all potential channel heads in a DEM. It looks for
  /// channel heads based on the outlet junctions of the valleys (which are identified by looking
  /// for portions of the landscape with 10 or more nodes with a high curvature that are linked)
  /// @param ValleyJunctions
  /// @param MinSegLength
  /// @param A_0
  /// @param m_over_n
  /// @param FlowInfo
  /// @param FlowDistance
  /// @param ElevationRaster
  /// @return vector<int> a vector of node_indices of potential channel heads
  /// @author FC
  /// @date 31/10/2013
vector<int> GetChannelHeadsChiMethodFromValleys(vector<int> ValleyNodes,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster);

   /// @brief This function returns all potential channel heads in a DEM. It looks for
  /// channel heads based on the valley source nodes identified as concave parts of the landscape
  /// @param ValleyJunctions
  /// @param MinSegLength
  /// @param A_0
  /// @param m_over_n
  /// @param FlowInfo
  /// @param FlowDistance
  /// @param ElevationRaster
  /// @param NJunctions number of downstream junctions to run the channel profiles from
  /// @return vector<int> a vector of node_indices of potential channel heads
  /// @author FC
  /// @date 10/09/15

  vector<int> GetChannelHeadsChiMethodFromSources(vector<int> ValleySources,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster, int NJunctions);

  /// @brief This function returns all channels in the DEM that the DrEICH algorithm uses for segiment fitting.
  /// It looks for channels based on the outlet junctions of valleys.
  /// It returns a LSDIndexRaster with the channels.
  /// @param ValleyNodes
  /// @param MinSegLength
  /// @param A_0
  /// @param m_over_n
  /// @param FlowInfo
  /// @param FlowDistance
  /// @param ElevationRaster
  /// @return LSDIndexRaster with all channels
  /// @author FJC
  /// @date 21/08/15
  LSDIndexRaster GetChannelsDreich(vector<int> ValleySources, int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& ElevationRaster, string path_name, int NJunctions);


  /// @brief This function returns a 2D array containing the locations of all pixels identified
  /// as being part of the channel using chi profiles.  It calculates the chi and elevation value
  /// of every pixel upstream of the given junction, then bins this data and calculates the pixels
  /// in the 95th percentile of each bin.  Any pixels above the 95th percentile are considered part
  /// of the channel, and any below are considered to be hillslopes.  This is the first part of the
  /// channel head prediction using chi profiles.
  /// @param JunctionNumber
  /// @param A_0
  /// @param m_over_n
  /// @param bin_width
  /// @param FlowInfo Flow Info object
  /// @param ElevationRaster
  /// @return Array2D<float> with channel pixels
  /// @author FC
  /// @date 01/10/2013
  Array2D<int> GetChannelHeadsChiMethodAllPixels(int JunctionNumber,
                                      float A_0, float m_over_n, float bin_width, LSDFlowInfo& FlowInfo,
                                      LSDRaster& ElevationRaster);


  /// @brief This function returns an integer vector with the node indexes of the furthest upstream
  /// pixels identified as being part of the channel using chi profiles.  It calculates the chi and
  /// elevation value of every pixel upstream of the given junction, then bins this data and calculates
  /// the pixels in the 95th percentile of each bin.  Any pixels above the 95th percentile are considered
  /// part of the channel, and any below are considered to be hillslopes.  This is the first part of the
  /// channel head prediction using chi profiles.
  /// @param JunctionNumber
  /// @param A_0
  /// @param m_over_n
  /// @param bin_width
  /// @param FlowInfo Flow Info object
  /// @param ElevationRaster
  /// @return vector<int> with source nodes
  /// @author FC
  /// @date 04/10/2013
  vector<int> GetSourceNodesChiMethodAllPixels(int JunctionNumber,
                                      float A_0, float m_over_n, float bin_width, LSDFlowInfo& FlowInfo,
                                      LSDRaster& ElevationRaster);

  // channel head identification
  /// @brief This function is used to predict channel head locations based on the method proposed by Pelletier (2013).
  ///
  /// @details It creates a contour curvature map and identifies channel heads as pixels greater
  /// than a user defined contour curvature threshold value, set by default at 0.1.  The threshold curvature
  /// can also be defined as a multiple of the standard deviation of the curvature.  Before this function is called
  /// the DEM must be filtered using the wiener filter in the LSDRasterSpectral object in order to remove high frequency
  /// noise.
  ///
  /// Reference: Pelletier (2013) A robust, two-parameter method for the extraction of drainage
  /// networks from high-resolution digital elevation models (DEMs): Evaluation using synthetic and real-world
  /// DEMs, Water Resources Research 49: 1-15
  ///
  /// @param tan_curv_threshold Double curvature threshold value.
  /// @param FlowInfo Flow Info object
  /// @param tan_curv_array 2D array of tangential curvature.
  /// @return 2D array of predicted channel head locations.
  /// @author FC
  /// @date 16/07/13
  vector<int> calculate_pelletier_channel_heads(float tan_curv_threshold, LSDFlowInfo& FlowInfo,
                                                Array2D<float>& tan_curv_array);

  /// @brief This function predicts channel head locations based on a tangential threshold
  /// as proposed by Pelletier (2013).
  ///
  /// @detail This function is used to predict channel head locations based on the method
  /// proposed by Pelletier (2013).  First it creates a contour curvature map and
  /// identifies channel heads as pixels with tangential curvature greater than a user
  /// defined threshold value. This map needed to be reduced to give the source
  /// pixels only.  This is done by i) sorting all the possible sources by
  /// elevation and ii) routing flow from each potential source using an adaption
  /// of Freeman MD flow.  Any potential sources that are located on ANY down-slope
  /// pathway within convergent part of the topography from previously visited
  /// source pixels are excluded.
  ///
  /// Reference: Pelletier (2013) A robust, two-parameter method for the extraction of
  /// drainage networks from high-resolution digital elevation models (DEMs): Evaluation
  /// using synthetic and real-world DEMs, Water Resources Research 49: 1-15
  ///
  /// @param FlowInfo object
  /// @param raster containing elevation data
  /// @param a threshold value of tangential curvature
  /// @param an array of tangential curvature
  /// @author DTM
  /// @date 03/06/2014
  vector<int> calculate_pelletier_channel_heads_DTM(LSDFlowInfo& FlowInfo, Array2D<float> topography, float tan_curv_threshold, Array2D<float>& tan_curv_array, Array2D<float>& tan_curv_array_LW);

  /// @brief This function identifies upstream limit of channel network
  ///
  /// @detail This is used to reduce a map of channel pixels down to a vector of sources
  /// for channel extraction.  It finds the upstream limit of each channel and then
  /// removes channelised pixels that are on ANY downslope pathway, within convergent part
  /// of the topography, from previous sources.  It uses a similar algorithm to the
  /// Freeman multi-directional flow routing algorithm in the LSDRaster object.
  ///
  /// @param FlowInfo object
  /// @param raster containing elevation data
  /// @param a vector of row coordinates for possible source pixels
  /// @param a vector of column coordinates for possible source pixels
  /// @param an array of tangential curvature
  /// @author DTM
  /// @date 03/06/2014
  vector<int> identify_upstream_limits(LSDFlowInfo& FlowInfo, Array2D<float>& topography,
                  vector<int> source_row_vec,vector<int> source_col_vec, Array2D<float>& tan_curv);

  /// @brief This function is used to identify concave portions of the landscape using a tangential curvature threshold.
  ///
  /// @details It defines the threshold based on a multiple of the standard deviation
  /// of the curvature.  It then identifies valleys in which there are a linked series of pixels
  /// which have a curvature value greater than the threshold, and finds the outlet junction number
  /// of this valley.  This can be passed to the channel head prediction algorithm using the chi
  /// method.
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param tan_curv_array 2D array with curvature
  /// @param sources vector with sources of channel network
  /// @param no_connecting_nodes number of nodes that need to be above the threshold before the valley is identified
  /// @return Array2D<int> with nodes at the base of each of the valleys
  /// @author FC
  /// @date 29/10/2013
  Array2D<int> find_valleys(LSDFlowInfo& FlowInfo, Array2D<float>& tan_curv_array,
                            vector<int> sources, int no_connecting_nodes, float tan_curv_threshold = 0.1);

  /// @brief This function is used to get the outlet nodes from a vector of input source nodes
  ///
  /// @details It is used to get a list of valley nodes that can be used in the DrEICH algorithm.
  /// The function goes downstream from each source node until the stream order of the downstream node is greater than
  /// that of the current node
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param sources vector with sources of channel network
  /// @return vector<int> with node at the base of each of the valleys
  /// @author FJC
  /// @date 19/08/2015
  vector<int> get_outlet_nodes_from_sources(LSDFlowInfo& FlowInfo, vector<int> sources);

  /// @brief This function is used to identify concave portions of the landscape using a tangential curvature threshold
  /// which is adaptive for each portion of the landscape
  ///
  /// @details It defines the threshold based on the standard deviation
  /// of the curvature.  It then identifies valleys in which there are a linked series of pixels
  /// which have a curvature value greater than the threshold, and finds the outlet junction number
  /// of this valley.  This can be passed to the channel head prediction algorithm using the chi
  /// method.
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param tan_curv_array 2D array with curvature
  /// @param sources vector with sources of channel network
  /// @param no_connecting_nodes number of nodes that need to be above the threshold before the valley is identified
  /// @param tan_curv_threshold array with the curvature thresholds for each row and col
  /// @return Array2D<int> with nodes at the base of each of the valleys
  /// @author FC
  /// @date 29/10/2013
  Array2D<int> find_valleys_adaptive_threshold(LSDFlowInfo& FlowInfo, Array2D<float>& tan_curv_array, vector<int> sources, int no_connecting_nodes, Array2D<float>& tan_curv_threshold);

  /// @brief This function uses a predefined channel mask to locate valley junctions
  ///
  /// @details This uses the same approach as the find_valleys function, but allows
  /// greater flexibility in how the valley network is defined
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param channel_mask 2D binary array with multi-pixel channel network marked by 1s
  /// @param sources vector with sources of channel network
  /// @param no_connecting_nodes number of nodes that need to be above the threshold before the valley is identified
  /// @return Array2D<int> with nodes at the base of each of the valleys
  Array2D<int> find_valleys_using_channel_mask(LSDFlowInfo& FlowInfo, Array2D<int>& channel_mask, vector<int> sources, int no_connecting_nodes);


  /// @brief Ridge network extraction - extracts ridge network, defined as boundaries
  /// between two basins of the same stream order.
  ///
  /// @details This function extracts the ridge network by defining it as the co-boundaries
  /// of basins of equivalent order, for all basin orders within the landscape.
  /// This is relatively trivial since in each array containing the basins of the
  /// same order, each basin is labelled with a unique identifier, thus co-
  /// boundaries are found by locating pixels that neighbour pixels from another
  /// basin of the same order.
  ///
  /// Updated to return an LSDRaster object as ridges can now be assigned CHT values,
  /// using LSDRaster::RidgeSample, which are not integers. - SWDG
  /// @param FlowInfo LSDFlowInfo object.
  /// @return LSDRaster of ridges.
  /// @author DTM
  /// @date 18/10/2012
  LSDRaster ExtractRidges(LSDFlowInfo& FlowInfo);


  /// @brief
  ///
  /// @details This overloaded function extracts the ridge network for a defined stream
  /// order, passed in by the user.
  ///
  /// Updated to return an LSDRaster object as ridges can now be assigned CHT values,
  /// using LSDRaster::RidgeSample, which are not integers. - SWDG
  /// @param FlowInfo LSDFlowInfo object.
  /// @param min_order Lowest order of ridges to extract.
  /// @param max_order Highest order of ridges to extract.
  /// @return LSDRaster of ridges.
  /// @author DTM, SWDG
  /// @date 18/10/2012, 28/03/2013
  LSDRaster ExtractRidges(LSDFlowInfo& FlowInfo, int& min_order, int& max_order);

  /// @brief This extracts all ridges defined by every junction. It takes the penultamate node
  ///  before the next downstream junction, and then extract that basin. It then gets the outline
  ///  of that basin for tagging as ridgeline nodes. Therefore every junction, except for baselevel junctions,
  ///  will have ridgelines.
  /// @detail Overlapping is based on the stack so higher order junctions will overwrite the ridgelines of
  ///  lower order junctions.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return A map with the key being the node index and the value being the largest stream order of the ridge
  /// @author SMM
  /// @date 11/01/2021
  map<int,int> ExtractAllRidges(LSDFlowInfo& FlowInfo);

  /// @brief This extracts all ridges defined by every junction. It takes the penultamate node
  ///  before the next downstream junction, and then extract that basin. It then gets the outline
  ///  of that basin for tagging as ridgeline nodes. Therefore every junction, except for baselevel junctions,
  ///  will have ridgelines.
  /// @detail Overlapping is based on the stack so higher order junctions will overwrite the ridgelines of
  ///  lower order junctions.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param csv_name the name of the csv (with extension) to which to print the data
  /// @return A map with the key being the node index and the value being the largest stream order of the ridge
  /// @author SMM
  /// @date 11/01/2021
  map<int,int> ExtractAllRidges(LSDFlowInfo& FlowInfo, string csv_name);

  /// @brief This last function gets the hilltops: ridges limited by a maximum slope.
  ///
  /// @details Resticts ridgeline to part of ridge network where the slope is less than a
  /// threshold value.
  ///
  /// Now outputs an LSDRaster to fall in line with other hilltop tools - SWDG 29/8/13
  /// @param RidgeRaster LSDIndexRaster of extracted ridges.
  /// @param SlopeRaster LSDRaster of slope.
  /// @param MaxSlope Maximum threshold slope value.
  /// @return LSDIndexRaster of hilltops.
  /// @author DTM
  /// @date 01/04/2013
  LSDRaster ExtractHilltops(LSDRaster& RidgeRaster, LSDRaster& SlopeRaster, float MaxSlope);

  /// @brief This function iterates through the junction nodes and assigns the unique
  /// junction ID to every stream pixel.
  ///
  /// @details This can be used with the LSDRaster
  /// hilltop_flow_routing function to assign a unique ID to each hilltop
  /// section tying it to a specific section of the channel network.
  /// @param flowinfo LSDFlowInfo object.
  /// @return LSDIndexRaster of the indexed channel newtork.
  /// @author SWDG
  /// @date 04/04/13
  LSDIndexRaster ChannelIndexer(LSDFlowInfo& flowinfo);

  /// @brief This extracts vectors containing node indices, junction indices
  ///  and stream orders of pixels in the channel network.
  /// @detail The vectors are replaced by the method
  /// @author SMM
  /// @date 14/11/2016
  void GetChannelNodesAndJunctions(LSDFlowInfo& flowinfo, vector<int>& NIvec, vector<int>& JIvec, vector<int>& SOvec);


  /// SplitChannel
  /// @brief This function splits the channel into a series of segments, providing a
  /// convenient unit with which to analyse landscapes.  The user provides the
  /// TargetSegmentLength, which specifies how many nodes should be in each
  /// segment, and a MinimumSegmentLength, which specifies the fewest permissable
  /// number of nodes.  Segments smaller than this are amalgamated into the
  /// upstream segment.
  /// The algorithm loops through the sources and traces downstream, stopping a
  /// segment after the target segment length, when the stream order increases (to
  /// preserve structure of drainage network), or when a channel pixel has already
  /// been visited.
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param Sources a vector of sources
  /// @param TargetSegmentLength (suggest 200 for 1m DEM)
  /// @return LSDIndexRaster with channel segments labelled by unique ID
  /// @author DTM
  /// @date 23/10/2013
  LSDIndexRaster SplitChannel(LSDFlowInfo& FlowInfo, vector<int> Sources, int TargetSegmentLength);


  /// TypologyModel
  /// @details This function splits the channel into a series of segments,
  /// providing a convenient unit with which to analyse landscapes.
  /// Function modified from original SplitChannel function so that the
  /// segment length varies with the drainage area of the cahtchment.
  /// Length (m) is calculated based on:
  /// L = Min_reach_length * sqrt(Drainage Area (km))
  /// User must pass in the minimum reach length in metres
  /// The algorithm starts a new segment either after the target length,
  /// when the stream order increases, or when a channel pixel has already
  /// been visited.
  /// User must pass in an empty IndexRaster which will be populated with the channel
  /// segments data, and two vector of vectors which will be populated:
  /// vector< vector<int> > SegmentInfoInts has the following layout:
  /// 0 - segment IDS
  /// 1 - start node of each segment (upstream)
  /// 2 - end nodes (downstream)
  /// vector< vector<float> > SegmentInfoFloats has the following layout:
  /// 0 - segment lengths
  /// 1 - elevation of the start nodes
  /// 2 - slope of the segment
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param Sources a vector of sources
  /// @param BaselineSources vector of baseline DRN sources
  /// @param CatchIDs vector of catchment IDs from DRN
  /// @param HydroCodes vector of hydrocodes from DRN
  /// @param MinReachLength in metres
  /// @param search_radius search radius for snapping rasters to the channel segments (pixels)
  /// @param ElevationRaster raster with elevation values
  /// @param DischargeRaster raster with discharge values (CEH one is in l/second, weirdly)
  /// @param ChannelSegments empty LSDIndexRaster, returned with channel segments labelled by unique ID
  /// @param SegmentInfoInt vec<vec> with integer segment info
  /// @param SegmentInfoFloat vec<vec> with floating segment info
  /// @author FJC
  /// @date 06/02/17
  void TypologyModel(LSDFlowInfo& FlowInfo, vector<int> Sources, vector<int> BaselineSources, vector<int> CatchIDs, vector<int> HydroCodes, int MinReachLength, int search_radius, LSDRaster& ElevationRaster, LSDRaster& DischargeRaster, LSDIndexRaster& ChannelSegments, vector< vector<int> >& SegmentInfoInts, vector< vector<float> >& SegmentInfoFloats);

  /// @brief This function removes channel segments from the typology model which are not downstream of a given
  /// list of source nodes
  /// @param FlowInfo LSDFlowInfo object
  /// @param Sources vector of source nodes
  /// @param SegmentInfoInts vec<vec> of segment info (integer)
  /// @param SegmentInfoFloats vec<vec> of segment info (floating point)
  void remove_tributary_segments(LSDFlowInfo& FlowInfo, vector<int> Sources, vector <vector <int> >& SegmentInfoInts, vector <vector <float> >& SegmentInfoFloats);

  /// @brief This function prints information about the channel segments from the
  /// TypologyModel function to a csv file so it can be read by a GIS
  /// @param FlowInfo LSDFlowInfo object
  /// @param SegmentInfoInts vec<vec> of segment info (integer)
  /// @param SegmentInfoFloats vec<vec> of segment info (floating point)
  /// @param outfilename string, csv filename
  void print_channel_segments_to_csv(LSDFlowInfo& FlowInfo, vector <vector <int> > SegmentInfoInts, vector <vector <float> > SegmentInfoFloats, string outfilename);

  /// SplitHillslopes
  /// @brief This function is intended to follow the SplitChannel function.  It traces
  /// through the receiver nodes from every hillslope pixel and then assigns them
  /// an integer value that matches the index of the section of channel that is
  /// setting the base level of that hillslope.
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param ChannelSegmentsRaster a raster of channel segments, produced by the SplitChannel function
  /// @return LSDIndexRaster hillslope segments labelled by ID of channel segments
  /// @author DTM
  /// @date 29/10/2013
  LSDIndexRaster SplitHillslopes(LSDFlowInfo& FlowInfo, LSDIndexRaster& ChannelSegmentsRaster);

  /// SplitHillslopes
  /// @brief This is an overloaded function doing the same as the previous version to
  /// segment hillslopes according to the channel index of the channel setting its
  /// base level.  However, this has been adapted to include an additional input
  /// raster - MultiThreadChannelRaster - which recognises that real channels may
  /// be multithreaded and/or have widths greater than or equal to one pixel.
  /// To be rigourous, these should be removed from analyses of hillslope
  /// properties.
  ///
  /// @param FlowInfo LSDFlowInfo object
  /// @param ChannelSegmentsRaster a raster of channel segments, produced by the SplitChannel function
  /// @param MultiThreadChannelRaster a binary raster with the full channel extent
  /// @return LSDIndexRaster hillslope segments labelled by ID of channel segments
  /// @author DTM
  /// @date 29/10/2013
  LSDIndexRaster SplitHillslopes(LSDFlowInfo& FlowInfo, LSDIndexRaster& ChannelSegmentsRaster,
                                 LSDIndexRaster& MultiThreadChannelRaster);


  // simple functions for getting streams. These do not return channel data elements but
  // instead return an LSDIndexRaster with the streams of a given order retained

  /// @brief Quick and dirty way to get channels of a defined stream order.
  ///
  /// @details No input error handling, will return an LSDIndexRaster of NoDataValues if an erroneous order is passed in.
  /// @param order Integer of the required stream order.
  /// @return LSDIndexRaster of the desired channels.
  /// @author SWDG
  /// @date 04/13
  LSDIndexRaster GetStreams(int order);
  /// @brief Quick and dirty way to get channels of a defined range of stream orders.
  ///
  /// @details No input error handling, will return an LSDIndexRaster of NoDataValues if an erroneous order is passed in.
  /// @param min_order Integer of the miniumum required stream order.
  /// @param max_order Integer of the max required stream order.
  /// @return LSDIndexRaster of the desired channels.
  /// @author SWDG
  /// @date 04/13
  LSDIndexRaster GetStreams(int min_order, int max_order);

  /// @brief Function to test whether a junction's upstream nodes border areas of No Data
  /// important to ensure basins are not being artificially truncated.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param input_junction Junction to be tested.
  /// @return Boolean indicating if no data values are present or not: \n
  /// false (0) = only good data values \n
  /// true (1) = no data values present \n
  /// Updated 24/10/13 to handle junction numbers in the same way that the basin extraction code does,
  /// by searching one junction downstream of the given junction and then going back up by one node - SWDG.
  /// @author SWDG
  /// @date 27/06/2013
  bool node_tester(LSDFlowInfo& FlowInfo, int input_junction);

  /// @brief Function to snap input coordinates to the nearest junction. This
  /// enables easy extraction of a particular catchment for analysis.
  /// @param X_coordinate of point. In coordiantes of DEM (usually UTM).
  /// @param Y_coordinate of point. In coordiantes of DEM (usually UTM).
  /// @param FlowInfo LSDFlowInfo object.
  /// @author DTM
  /// @date 17/10/2013
  int get_receiver_junction_for_specified_coordinates(float X_coordinate, float Y_coordinate, LSDFlowInfo& FlowInfo);


  /// @brief Function to snap input coordinates to the nearest channel node. This
  /// enables easy extraction of a particular catchment for analysis.
  /// @param X_coordinate of point. In coordiantes of DEM (usually UTM).
  /// @param Y_coordinate of point. In coordiantes of DEM (usually UTM).
  /// @param threshold_stream_order The minimum stream order that will be considers a 'channel' by the algorithm
  /// @param search_radius_nodes the radius of the kernal used to check if there is a nearby channel. A radius of 0 only
  /// includes the centre point, a radius of 1 has a kernal diameter of 3,  radius of 2 has a kernal diameter of 5
  /// and so on
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Returns the NodeIndex of the nearest channel node.
  /// @author SMM
  /// @date 21/10/2013
  int get_nodeindex_of_nearest_channel_for_specified_coordinates(float X_coordinate,float Y_coordinate, int threshold_stream_order, int search_radius_nodes,LSDFlowInfo& FlowInfo);

  /// @brief Function to get the nodeindex of the nearest channel node for a given nodeindex
  /// @param starting_nodeindex does what it says on the tin
  /// @param threshold_stream_order The minimum stream order that will be considers a 'channel' by the algorithm
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Returns the NodeIndex of the nearest channel node.
  /// @author SMM
  /// @date 02/07/2021
  int get_nodeindex_of_nearest_channel(int starting_nodeindex, int threshold_stream_order, LSDFlowInfo& FlowInfo);

  /// @brief Function to get the distance to the nearest channel
  /// @param starting_nodeindex does what it says on the tin
  /// @param threshold_stream_order The minimum stream order that will be considers a 'channel' by the algorithm
  /// @param FlowInfo LSDFlowInfo object.
  /// @param FlowDist the distance from outlet raster
  /// @return Returns the distance to the nearest channel node following a flow path
  /// @author SMM
  /// @date 02/07/2021
  float get_distance_to_nearest_channel(int starting_nodeindex, int threshold_stream_order, LSDFlowInfo& FlowInfo, LSDRaster& FlowDist);

  /// @brief Function to get the distance to the nearest channel
  /// @param threshold_stream_order The minimum stream order that will be considers a 'channel' by the algorithm
  /// @param FlowInfo LSDFlowInfo object.
  /// @param FlowDist the flow distance raster
  /// @return Returns the distance of all pixels to the nearest channel node as a raster
  /// @author SMM
  /// @date 02/07/2021
  LSDRaster get_distance_to_nearest_channel_raster(int threshold_stream_order, LSDFlowInfo& FlowInfo, LSDRaster& FlowDist);

  /// @brief Function to get the releif to the nearest channel
  /// @param starting_nodeindex does what it says on the tin
  /// @param threshold_stream_order The minimum stream order that will be considers a 'channel' by the algorithm
  /// @param FlowInfo LSDFlowInfo object.
  /// @param Elevation the elevation raster
  /// @return Returns the relief of this pixel to the nearest channel node following a flow path
  /// @author SMM
  /// @date 02/07/2021
  float get_relief_to_nearest_channel(int starting_nodeindex, int threshold_stream_order, LSDFlowInfo& FlowInfo, LSDRaster& Elevation);

  /// @brief Function to get the distance to the nearest channel
  /// @param threshold_stream_order The minimum stream order that will be considers a 'channel' by the algorithm
  /// @param FlowInfo LSDFlowInfo object.
  /// @param Elevation the elevation raster
  /// @return Returns the relief of this pixel to the nearest channel node as a raster
  /// @author SMM
  /// @date 02/07/2021
  LSDRaster get_relief_to_nearest_channel_raster(int threshold_stream_order, LSDFlowInfo& FlowInfo, LSDRaster& Elevation);

  /// @brief Function to snap input coordinates to the nearest channel node from latitude and longitude
  /// @param latitude
  /// @param longitude
  /// @param FlowInfo LSDFlowInfo object.
  /// @param Converter LSDCoordinateConverterLLandUTM object
  /// @return Returns the NodeIndex of the nearest channel node.
  /// @author FJC
  /// @date 20/11/17
  int get_junction_of_nearest_channel_from_lat_long(double latitude, double longitude, LSDFlowInfo& FlowInfo, LSDCoordinateConverterLLandUTM Converter);
  int get_upstream_junction_from_lat_long(double latitude, double longitude, LSDFlowInfo& FlowInfo, LSDCoordinateConverterLLandUTM Converter);

  /// @brief Function to snap input coordinates to the nearest channel node from latitude and longitude
  /// @param X_coordinate in local coordinates
  /// @param Y_coordinate in local coordinates
  /// @param threshold_SO The threshold stream order
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Returns the NodeIndex of the nearest channel node that has a stream order greater than the threshold stream order.
  /// @author SMM
  /// @date 01/08/19
  int find_nearest_downslope_channel(float X_coordinate, float Y_coordinate, int threshold_SO, LSDFlowInfo& FlowInfo);

  /// @brief Function to snap input coordinates to the nearest channel node from latitude and longitude
  /// @param starting node: the node index of the starting node
  /// @param threshold_SO The threshold stream order
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Returns the NodeIndex of the nearest channel node that has a stream order greater than the threshold stream order.
  /// @author SMM
  /// @date 01/08/19
  int find_nearest_downslope_channel(int starting_node, int threshold_SO, LSDFlowInfo& FlowInfo);

	/// @brief Function to get info about the nearest channel node of a given node.
  /// @param StartingNode index of node of interest
  /// @param threshold_SO threshold stream order for finding the nearest channel
  /// @param FlowInfo LSDFlowInfo object
	/// @param DistFromOutlet LSDRaster of flow lengths
	/// @param ChannelNode int to store the NI of the nearest channel
	/// @param FlowLength float to store the flow length to the nearest channel
	/// @param DistanceUpstream float to store the distance upstream of the nearest channel
  /// @author FJC
  /// @date 29/09/16
	void get_info_nearest_channel_to_node(int& StartingNode, int& threshold_SO, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet, int& ChannelNode, float& FlowLength, float& DistanceUpstream);

	/// @brief Function to get info about the nearest channel node on the main stem of a given node.
  /// @param StartingNode index of node of interest
  /// @param FlowInfo LSDFlowInfo object
	/// @param ElevationRaster LSDRaster of elevations
	/// @param DistFromOutlet LSDRaster of flow lengths
	/// @param ChannelNode int to store the NI of the nearest channel
	/// @param FlowLength float to store the flow length to the nearest channel
	/// @param DistanceUpstream float to store the distance upstream of the nearest channel
	/// @param Relief float to store relief compared to nearest channel
  /// @author FJC
  /// @date 05/10/16
  void get_info_nearest_channel_to_node_main_stem(int& StartingNode, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, LSDIndexChannel& MainStem, int& ChannelNode, float& FlowLength, float& DistanceUpstream, float& Relief);


  /// @brief This function takes a node index, checks to see if it is on a channel,
  /// and then works its way up the channel to find the upstream junction
  /// @param ChannelNodeIndex is the node index of the channel node (if it isn't a channel
  /// the function returns NoDataValue
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Returns the Junction Number of the nearest upslope junction
  /// @author SMM
  /// @date 21/10/2013
  int find_upstream_junction_from_channel_nodeindex(int ChannelNodeIndex, LSDFlowInfo& FlowInfo);

  /// @brief This function checks whether any of the upstream nodes of a given junction are the same
  /// steam order as the junction itself. It returns an integer value which is 1 if the SO is the same
  /// and 0 if it is not the same.
  /// @param junction junction of interest
  /// @param FlowInfo LSDFlowInfo object.
  /// @return integer value 0 or 1
  /// @author FJC and MAH
  /// @date 18/03/16
  int check_stream_order_of_upstream_nodes(int junction, LSDFlowInfo& FlowInfo);

  /// @brief This function returns the node index of the donor node of a given node
  /// with the highest stream order
  /// @param current_node The current node index
  /// @param FlowInfo LSDFlowInfo object.
  /// @return node index of the donor node with the highest stream order
  /// @author FJC
  /// @date 31/01/17
  int get_upstream_node_max_stream_order(int current_node, LSDFlowInfo& FlowInfo);

  /// @brief this function is a wrapper that takes a list of x and y locations,
  ///  filters them to make sure they are in the data bounds,
  ///  and then calculates the nearest channel and junction.
  ///  It is primarily used to snap cosmo data to the channel network
  /// @param x_locs the x locations of the points
  /// @param y_locs the y locations of the points
  /// @param search_radius_nodes the number of nodes around the point to search
  ///  for a channel
  /// @param threshold_stream_order the minimum stream order to which the point
  ///  will snap
  /// @param FlowInfo the LSDFlowInfo object
  /// @param valid_cosmo_points a vector<int> of indices into the x and y vectors.
  ///  for example if the only valid points were at x_loc[12] and x_loc[34] this
  ///  would return a vector with two elements, 12 and 34. This vector is overwritten
  ///  by this function
  /// @param snapped_node_indices a vector containing the node indices of the
  ///  points snapped to the nearest channel (within search radius and over the
  ///  drainage order threshold). This is overwritten by this method.
  /// @param snapped_junction_indices a vector<int> continaing the junction numbers
  ///  downstream of the nearest channel node. This is overwritten by this method.
  /// @author SMM
  /// @date 14/11/2014
  void snap_point_locations_to_channels(vector<float> x_locs,
                vector<float> y_locs,
                int search_radius_nodes, int threshold_stream_order,
                LSDFlowInfo& FlowInfo, vector<int>& valid_cosmo_points,
                vector<int>& snapped_node_indices, vector<int>& snapped_junction_indices);

  /// @brief this function is a wrapper that takes a list of x and y locations,
  ///  filters them to make sure they are in the data bounds,
  ///  and then calculates the nearest channel node index
  /// @detail Unlike the other version of this function, the code snaps to a node index
  ///  rather than the nearest upslope junction.
  /// @param x_locs the x locations of the points
  /// @param y_locs the y locations of the points
  /// @param search_radius_nodes the number of nodes around the point to search
  ///  for a channel
  /// @param threshold_stream_order the minimum stream order to which the point
  ///  will snap
  /// @param FlowInfo the LSDFlowInfo object
  /// @param valid_cosmo_points a vector<int> of indices into the x and y vectors.
  ///  for example if the only valid points were at x_loc[12] and x_loc[34] this
  ///  would return a vector with two elements, 12 and 34. This vector is overwritten
  ///  by this function
  /// @param snapped_node_indices a vector containing the node indices of the
  ///  points snapped to the nearest channel.
  /// @author SMM
  /// @date 01/08/2019
  void snap_point_locations_to_nearest_channel_node_index(vector<float> x_locs,
                vector<float> y_locs,
                int search_radius_nodes, int threshold_stream_order,
                LSDFlowInfo& FlowInfo, vector<int>& valid_cosmo_points,
                vector<int>& snapped_node_indices);


  /// @brief This function takes the name of a csv file that contains
  ///  latitude and longitude information and then finds the junction
  ///  numbers upslope.
  /// @param csv_filename The name of the file from which the lat-long will be read.
  /// @param search_radius_nodes the number of nodes around the point to search
  ///  for a channel
  /// @param threshold_stream_order the minimum stream order to which the point
  ///  will snap
  /// @param FlowInfo the LSDFlowInfo object
  /// @author SMM
  /// @date 23/01/2020
  vector<int> snap_point_locations_to_upstream_junctions_from_latlong_csv(string csv_filename,
                int search_radius_nodes, int threshold_stream_order,
                LSDFlowInfo& FlowInfo, LSDRasterInfo& RI);

  ///@brief This function extract basins from list of nodes and overwrite given sources,
  /// baselevel_nodes, outlet_nodes and baselevel junctions.
  /// Be careful, the nodes have to be the EXACT outlet location for the basin. This
  /// intend to be used as a low-level function to extract watersheds AFTER preselecting nodes from other criterias.
  ///@param  input_nodes: node indices of the actual outlet location
  ///@param  sources: source nodes for the river: WILL BE TRIMMED TO THE SOURCES WITHIN THE BASINS
  ///@param  baselevel_nodes: WILL BE OVERWRITTEN with the new base_level nodes
  ///@param  baselevel_junctions: WILL BE OVERWRITTEN with the new baselevel_junctions nodes
  ///@param  outlet_nodes: WILL BE OVERWRITTEN with the new outlet_nodes nodes
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  LSDRaster: DistanceFromOutlet, a raster with the distance from oulet
  ///@param  check_edges: Check the influence of NoData over the basins. WARNING: this is a hard pass and just ignore the basins influence by nodata without finding equivalents.
  ///@author B.G.
  ///@date 09/12/2019
  void select_basin_from_nodes(vector<int>& input_nodes, vector<int>& sources, vector<int>& baselevel_nodes,
    vector<int>& baselevel_junctions,vector<int>& outlet_nodes,LSDFlowInfo& FlowInfo,
    LSDRaster& DistanceFromOutlet, bool check_edges);

  ///@brief This function extract basins from list of nodes and overwrite given sources,
  /// baselevel_nodes, outlet_nodes and baselevel junctions.
  /// It snaps the inputted node to the largest drainage area in a pixel radius.
  ///@param  input_nodes: node indices of the actual outlet location
  ///@param  sources: source nodes for the river: WILL BE TRIMMED TO THE SOURCES WITHIN THE BASINS
  ///@param  baselevel_nodes: WILL BE OVERWRITTEN with the new base_level nodes
  ///@param  baselevel_junctions: WILL BE OVERWRITTEN with the new baselevel_junctions nodes
  ///@param  outlet_nodes: WILL BE OVERWRITTEN with the new outlet_nodes nodes
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  DrainageArea: DA or Discharge raster to tell the code where to snap
  ///@param  n_pixels: N pixels to look for around the original point
  ///@param  check_edges: Check the influence of NoData over the basins. WARNING: this is a hard pass and just ignore the basins influence by nodata without finding equivalents.
  ///@param  LSDRaster: DistanceFromOutlet, a raster with the distance from oulet
  ///@author B.G.
  ///@date 09/12/2019
  void basin_from_node_snap_to_largest_surrounding_DA(vector<int>& input_nodes, vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  int n_pixels, bool check_edges);

  ///@brief Extract all the basins larger than a certain minimum drainage area in metre square
  ///@param  sources: source nodes for the river: WILL BE TRIMMED TO THE SOURCES WITHIN THE BASINS
  ///@param  baselevel_nodes: WILL BE OVERWRITTEN with the new base_level nodes
  ///@param  baselevel_junctions: WILL BE OVERWRITTEN with the new baselevel_junctions nodes
  ///@param  outlet_nodes: WILL BE OVERWRITTEN with the new outlet_nodes nodes
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  DrainageArea: DA or Discharge raster to tell the code where to snap
  ///@param  LSDRaster: DistanceFromOutlet, a raster with the distance from oulet
  ///@param  min_DA: the minimum drainage area
  ///@param  check_edges: Check the influence of NoData over the basins. WARNING: this is a hard pass and just ignore the basins influence by nodata without finding equivalents.
  ///@return a vector with the selected baselevel nodes
  ///@author B.G
  ///@date December 2019
  vector<int> basin_from_node_minimum_DA(vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  double min_DA, bool check_edges);


  ///@brief Takes a list of nodes, for example a river or mountain front, and extract all the basins draining to that line
  ///@param  input_nodes: node indices of the actual outlet location
  ///@param  sources: source nodes for the river: WILL BE TRIMMED TO THE SOURCES WITHIN THE BASINS
  ///@param  baselevel_nodes: WILL BE OVERWRITTEN with the new base_level nodes
  ///@param  baselevel_junctions: WILL BE OVERWRITTEN with the new baselevel_junctions nodes
  ///@param  outlet_nodes: WILL BE OVERWRITTEN with the new outlet_nodes nodes
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  DrainageArea: DA or Discharge raster to tell the code where to snap
  ///@param  LSDRaster: DistanceFromOutlet, a raster with the distance from oulet
  ///@param  min_DA: the minimum drainage area
  ///@return a vector with the selected baselevel nodes
  ///@author B.G
  ///@date December 2019
  std::vector<int> basin_from_node_minimum_DA_draining_to_list_of_nodes(vector<int>& input_nodes, vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  double min_DA);

  ///@brief  Takes an outlet node and extract all the subbasins draing to that mother-basin, whithin a range of size.
  ///@param  input_nodes: node indices of the actual outlet location
  ///@param  sources: source nodes for the river: WILL BE TRIMMED TO THE SOURCES WITHIN THE BASINS
  ///@param  baselevel_nodes: WILL BE OVERWRITTEN with the new base_level nodes
  ///@param  baselevel_junctions: WILL BE OVERWRITTEN with the new baselevel_junctions nodes
  ///@param  outlet_nodes: WILL BE OVERWRITTEN with the new outlet_nodes nodes
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  DrainageArea: DA or Discharge raster to tell the code where to snap
  ///@param  LSDRaster: DistanceFromOutlet, a raster with the distance from oulet
  ///@param  min_DA: the minimum drainage area
  ///@param  max_DA: the maximum drainage area
  ///@return a vector with the selected baselevel nodes
  ///@author B.G
  ///@date December 2019
  std::vector<int> basin_from_node_all_minimum_DA_for_one_watershed(int outlet_node_of_the_watershed, vector<int>& sources, vector<int>& baselevel_nodes,
vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
double min_DA, double max_DA);

  ///@brief  Extract all basins within a range of drainage area.
  ///@param  input_nodes: node indices of the actual outlet location
  ///@param  sources: source nodes for the river: WILL BE TRIMMED TO THE SOURCES WITHIN THE BASINS
  ///@param  baselevel_nodes: WILL BE OVERWRITTEN with the new base_level nodes
  ///@param  baselevel_junctions: WILL BE OVERWRITTEN with the new baselevel_junctions nodes
  ///@param  outlet_nodes: WILL BE OVERWRITTEN with the new outlet_nodes nodes
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  DrainageArea: DA or Discharge raster to tell the code where to snap
  ///@param  LSDRaster: DistanceFromOutlet, a raster with the distance from oulet
  ///@param  min_DA: the minimum drainage area
  ///@param  max_DA: the maximum drainage area
  ///@param  check_edges: Check the influence of NoData over the basins. WARNING: this is a hard pass and just ignore the basins influence by nodata without finding equivalents.
  ///@return a vector with the selected baselevel nodes
  ///@author B.G
  ///@date December 2019
  vector<int> basin_from_node_range_DA(vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  double min_DA, double max_DA, bool check_edges);

  ///@brief Takes a reference raster and its flowinfo and uses the stack to create a vector of bool telling which baselevel node would be influenced by nodata.
  ///@brief The index of the vector is the node index
  ///@param  FlowInfo: A fully-grown and weaned FlowInfo object
  ///@param  testRaster: a reference raster with the nodata
  vector<bool> check_nodata_influence(LSDFlowInfo& FlowInfo, LSDRaster& testrast);

  /// @brief This functions takes a junction number and then follwos the receiver
  /// junctions until it hits a baselevel junction.
  /// @param a junction node to start from.
  /// @return The base level junction to which the starting junction drains
  /// @author SMM
  /// @date 21/02/2014
  int find_base_level_node_of_junction(int StartingJunction);

  /// @brief Prints a list of junctions, with their locations in both UTM and
  ///  in lat long WGS1984 to file
  /// @detail The format of the file is:
  ///  junction,node,x,y,latitude,longitude
  /// @param FlowInfo an LSDFlowInfo object
  /// @param JunctionList A list of junctions in an integer vector
  /// @param fname The filename of the csv file
  /// @author SMM
  /// @date 18/05/2016
  void print_junctions_to_csv(LSDFlowInfo& FlowInfo, vector<int> JunctionList, string fname);

  /// @brief Prints all junctions, with their locations in both UTM and
  ///  in lat long WGS1984 to file
  /// @detail The format of the file is:
  ///  junction,node,x,y,latitude,longitude
  /// @param FlowInfo an LSDFlowInfo object
  /// @param fname The filename of the csv file
  /// @author SMM
  /// @date 21/04/2017
  void print_junctions_to_csv(LSDFlowInfo& FlowInfo, string fname);


  // Get functions

  /// @return Number of rows as an integer.
  int get_NRows() const        { return NRows; }

  /// @return Number of columns as an integer.
  int get_NCols() const        { return NCols; }

  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const      { return XMinimum; }

  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const      { return YMinimum; }

  /// @return Data resolution as an integer.
  float get_DataResolution() const  { return DataResolution; }

  /// @return No Data Value as an integer.
  int get_NoDataValue() const      { return NoDataValue; }

  /// @return Georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }

  /// @details Gets the node of a junction
  /// @param junction integer node index.
  /// @return Integer node of junction.
  /// @author SMM
  /// @date 01/01/2014
  int get_Node_of_Junction(int junction) const;

  /// @details Gets the receiver of a junction
  /// @param junction integer receiver index.
  /// @return Integer receiver of junction.
  /// @author SMM
  /// @date 01/01/2014
  int get_Receiver_of_Junction(int junction) const;

  /// @details Gets a node list from a junction list
  /// @param junction_list a vector of junctions
  /// @return vector of nodes
  /// @author SMM
  /// @date 20/01/2018
  vector<int> get_node_list_from_junction_list(vector<int> junction_list);


  /// @details Gets a node list from a junction list. The nodes are the penultimate nodes before the receiver juctions
  /// @param junction_list a vector of junctions
  /// @param FlowInfo the LSDFlowInfo object
  /// @return vector of nodes
  /// @author SMM
  /// @date 20/01/2018
  vector<int> get_node_list_of_penultimate_node_from_junction_list(vector<int> junction_list, LSDFlowInfo& FlowInfo);

  /// @details Get downstream junction
  /// @param starting_junction starting junction
  /// @param FlowInfo LSDFlowInfo object
  /// @return integer with downstream junction number
  /// @author FJC
  /// @date 08/10/15
  int get_downstream_junction(int starting_junction, LSDFlowInfo& FlowInfo);

  /// @details Get downstream junction for a specified node
  /// @param CurrentNode starting node
  /// @param FlowInfo LSDFlowInfo object
  /// @return integer with downstream junction number
  /// @author FJC
  /// @date 29/01/21
  int get_junction_downstream_of_node(int CurrentNode, LSDFlowInfo& FlowInfo);

  /// @details Gets the stream order of a node
  /// @param FlowInfo LSDFlowInfo object
  /// @param node node of interest
  /// @return integer with stream order of junction
  /// @author FJC
  /// @date 29/09/16
	int get_StreamOrder_of_Node(LSDFlowInfo& FlowInfo, int node);

  /// @details Gets the stream order of a junction
  /// @param FlowInfo LSDFlowInfo object
  /// @param junction the junction of interest
  /// @return integer with stream order of junction
  /// @author FJC
  /// @date 20/03/14
  int get_StreamOrder_of_Junction(LSDFlowInfo&FlowInfo, int junction);

  /// @details Gets the stream order of a junction
  /// @param junction the junction of interest
  /// @return integer with stream order of junction
  /// @author SMM
  /// @date 26/10/14
  int get_StreamOrder_of_Junction(int junction);

  /// @details This gets the junction that is at the next
  /// Strahler stream order from the current junction.
  /// @param junction the current junction
  /// @return the junction that is at the next stream order
  /// if the next stream order is not reached before baselevel
  /// it returns a NoDataValue
  /// @author SMM
  /// @date 26/10/2014
  int get_Next_StreamOrder_Junction(int junction);

  /// @details Returns an bool to check whether junction
  ///  is upstream of another base level
  /// @param current_junction the junction of interest
  /// @param test_junction the junction to see if it is upstream
  /// @ return true or false
  /// @author SMM
  /// @date 23/06/2017
  bool is_junction_upstream(int current_junction, int test_junction);

  /// @details Returns an integer to check whether junction
  /// is at base level.
  /// @param junction the junction of interest
  /// @ return int  1 = base level, 0 = not base level
  /// @author FJC
  /// @date 11/01/2017
  int is_Junction_BaseLevel(int junction);

  /// @return The number of junctions
  int get_NJunctions() const { return int(JunctionVector.size()); }

  /// @return The Vector of Junctions. Note that these are the node indices of the
  ///  junctions. The junction numbers just go from 0 to NJunctions
  vector<int> get_JunctionVector() const { return JunctionVector; }

  /// @return Get the baselevel junstions
  vector<int> get_BaseLevelJunctions() const { return BaseLevelJunctions; }

  /// @return The Vector of receivers.
  vector<int> get_ReceiverVector() const { return ReceiverVector; }

  /// @return The Vector of stream orders.
  vector<int> get_StreamOrderVector() const { return StreamOrderVector; }

  /// @return The vector of sources. The vector is composed of node indices
  vector<int> get_SourcesVector() const { return SourcesVector; }

  /// @return The SVector
  vector<int> get_SVector() const { return SVector; }

	/// @return the stream order array
	Array2D<int> get_StreamOrderArray() const { return StreamOrderArray; }

  void couple_hillslope_nodes_to_channel_nodes(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& D_inf_Flowdir, LSDIndexRaster& ChannelNodeNetwork, int OutletJunction, vector<int>& hillslope_nodes, vector<int>& baselevel_channel_nodes);

  /// @details This function calculates the relief of each pixel compared to the nearest downstream
  /// channel pixel equal or greater to the threshold stream order
  /// @param ElevationRaster LSDRaster with elevations
  /// @param FlowInfo LSDFlowInfo object
  /// @param threshold_SO threshold stream order to calculate relief from
  /// @return LSDRaster with channel relief
  /// @author FJC
  /// @date 17/11/15
  LSDRaster calculate_relief_from_channel(LSDRaster& ElevationRaster, LSDFlowInfo& FlowInfo, int threshold_SO);

	/// @details This function calculates the relief of each pixel compared to the nearest
	/// downstream channel pixel equal or greater to the threshold stream order for that
	/// connected components patch
  /// @param ElevationRaster LSDRaster with elevations
	/// @param ConnectedComponents connected components raster
	/// @param DistFromOutlet raster of flow lengths
  /// @param FlowInfo LSDFlowInfo object
  /// @param threshold_SO original threshold stream order to calculate relief from
  /// @return LSDRaster with channel relief
  /// @author FJC
  /// @date 29/09/16
	LSDRaster calculate_relief_from_channel_connected_components(LSDRaster& ElevationRaster, LSDIndexRaster& ConnectedComponents, LSDRaster& DistFromOutlet, LSDFlowInfo& FlowInfo, int threshold_SO, int search_distance);

	/// @details This function takes in a raster of connected component patches. It finds
	/// the elevation of the nearest channel for the patch.
  /// @param ConnectedComponents connected components raster
	/// @param ElevationRaster LSDRaster of elevations
  /// @param FlowInfo LSDFlowInfo object
	/// @param DistFromOutlet LSDRaster of flow lengths
  /// @param threshold_SO threshold stream order to calculate relief from
	/// @param search_distance length of channel reach to get elevation from
  /// @return 2D array with the elevation of nearest channel for each patch
  /// @author FJC
  /// @date 29/09/16
	Array2D<int> Get_Elevation_of_Nearest_Channel_for_Connected_Components(LSDIndexRaster& ConnectedComponents, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, LSDFlowInfo& FlowInfo, int threshold_SO, int search_distance);

	/// @details This function finds the mean elevation of the channel reach given a node on the channel network
  /// @param StartingNode node to check
  /// @param search_distance reach distance - will check both upstream and downstream this distance
  /// @param ElevationRaster elevation raster
	/// @param FlowInfo LSDFlowInfo object
  /// @return mean elevation of reach
  /// @author FJC
  /// @date 29/09/16
	float find_mean_elevation_of_channel_reach(int StartingNode, int search_distance, LSDRaster& ElevationRaster, LSDFlowInfo& FlowInfo);

		/// @details This function returns the node index of the nearest FIP to a
	/// specified node.
  /// @param point_node Node to start with
  /// @param search_distance Distance to search upstream and downstream for a FIP
  /// @param FloodplainRaster Raster with binary floodplain
	/// @param FlowInfo LSDFlowInfo object
  /// @return node index of nearest FIP
  /// @author FJC
  /// @date 09/09/16
	float find_distance_to_nearest_floodplain_pixel(int point_node, int search_distance, LSDRaster& FloodplainRaster, LSDFlowInfo& FlowInfo);

  /// @detail This overwrites two vecotrs that give all of the starting and
  ///  finishing nodes of channels in a basin
  /// @param FlowInfo an LSDFlowInfo object
  /// @param BaseLevel_Junctions an integer vector that contains the base level junctions
  /// @param DistanceFromOutlet an LSDRaster with the flow distance
  /// @param source_nodes a vector continaing the sorted sorce nodes (by flow distance)
  ///  THIS GETS OVERWRITTEN
  /// @param outlet_nodes a vector continaing the outlet nodes
  ///  THIS GETS OVERWRITTEN
  /// @author SMM
  /// @date 20/05/2016
  void get_overlapping_channels(LSDFlowInfo& FlowInfo, vector<int> BaseLevel_Junctions,
                                LSDRaster& DistanceFromOutlet,
                                vector<int>& source_nodes, vector<int>& outlet_nodes,
                                int n_nodes_to_visit);

  /// @detail This overwrites two vecotrs that give all of the starting and
  ///  finishing nodes of channels in a basin
  /// @param FlowInfo an LSDFlowInfo object
  /// @param BaseLevel_Junctions an integer vector that contains the base level junctions
  /// @param DistanceFromOutlet an LSDRaster with the flow distance
  /// @param source_nodes a vector continaing the sorted sorce nodes (by flow distance)
  ///  THIS GETS OVERWRITTEN
  /// @param outlet_nodes a vector continaing the outlet nodes
  ///  THIS GETS OVERWRITTEN
  /// @param baselevel_nodes a vector continaing the baselevel nodes (i.e. the node of the outlet of the basin)
  ///  THIS GETS OVERWRITTEN
  /// @author SMM
  /// @date 21/06/2017
  void get_overlapping_channels(LSDFlowInfo& FlowInfo, vector<int> BaseLevel_Junctions,
                                LSDRaster& DistanceFromOutlet,
                                vector<int>& source_nodes, vector<int>& outlet_nodes,
                                vector<int>& baselevel_nodes,
                                int n_nodes_to_visit);

  /// @detail This overwrites two vectors that give all of the starting and
  ///  finishing nodes of channels in a basin continuing downstream from the selected junction to its outlet
  /// @param FlowInfo an LSDFlowInfo object
  /// @param BaseLevel_Junctions an integer vector that contains the base level junctions
  /// @param DistanceFromOutlet an LSDRaster with the flow distance
  /// @param source_nodes a vector continaing the sorted sorce nodes (by flow distance)
  ///  THIS GETS OVERWRITTEN
  /// @param outlet_nodes a vector continaing the outlet nodes
  ///  THIS GETS OVERWRITTEN
  /// @author MDH
  /// @date 16/6/2017
  void get_overlapping_channels_to_downstream_outlets(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    int n_nodes_to_visit);

  /// @detail This overwrites two vectors that give all of the starting and
  ///  finishing nodes of channels in a basin continuing downstream from the selected junction to its outlet
  /// @brief I THINK THIS MIGHT CAUSE A SEG FAULT: NEED TO UPDATE!!!!!
  /// @param FlowInfo an LSDFlowInfo object
  /// @param BaseLevel_Junctions an integer vector that contains the base level junctions
  /// @param DistanceFromOutlet an LSDRaster with the flow distance
  /// @param source_nodes a vector continaing the sorted sorce nodes (by flow distance)
  ///  THIS GETS OVERWRITTEN
  /// @param outlet_nodes a vector continaing the outlet nodes
  ///  THIS GETS OVERWRITTEN
  /// @param baselevel_nodes a vector continaing the baselevel nodes (i.e. the node of the outlet of the basin)
  ///  THIS GETS OVERWRITTEN
  /// @author SMM
  /// @date 21/6/2017
  void get_overlapping_channels_to_downstream_outlets(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    vector<int>& baselevel_nodes,
                                    int n_nodes_to_visit);

  ///
  void get_overlapping_channels_to_downstream_outlets(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    vector<int>& baselevel_nodes,
                                    int n_nodes_to_visit,
                                    vector<int>& input_nodes);

/// @detail This function gets all the pixels along a line defined by a series of points and finds the pixels greater than a specified stream order.
/// @param Points PointData object with the points
/// @param ElevationRaster raster of elevations
/// @param threshold_SO threshold stream order
/// @param FlowInfo LSDFlowInfo object
/// @author FJC
/// @date 17/04/17
vector<int> get_channel_pixels_along_line(vector<int> line_rows, vector<int> line_cols, int threshold_SO, LSDFlowInfo& FlowInfo);

/// @brief function to take a vector of basin outlet junctions and write data about the longest channel in each to csv.
/// @param BasinJunctions vector of basin junctions
/// @param FlowInfo
/// @param DistanceFromOutlet
/// @param Elevation elev raster
/// @param csv_filename the output csv file name
/// @param window_size the total window size (in channel nodes) for calculating the channel slopes over.
/// @author FJC
/// @date 06/04/18
void write_river_profiles_to_csv(vector<int>& BasinJunctions, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& Elevation, string csv_filename, int window_size);

/// @brief function to take a vector of basin outlet junctions and write data about all tribs to csv
/// @param BasinJunctions vector of basin junctions
/// @param FlowInfo
/// @param DistanceFromOutlet
/// @param Elevation elev raster
/// @param csv_filename the output csv file name
/// @author FJC
/// @date 02/05/18
void write_river_profiles_to_csv_all_tributaries(vector<int>& BasinJunctions, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& Elevation, string csv_filename);

/// @brief function get the total length of channels upstream of a node
/// @param this_node node of interest
/// @param FlowInfo LSDFlowInfo object
/// @author FJC
/// @date 30/04/18
float GetTotalChannelLengthUpstream(int this_node, LSDFlowInfo& FlowInfo);

/// @brief function to write data about channels downstream of all channel heads for a specified length
/// @param channel_length length downstream to stop writing info
/// @param FlowInfo
/// @param Elevation elev raster
/// @param csv_filename the output csv file name
/// @author FJC
/// @date 02/05/18
void write_river_profiles_to_csv_all_sources(float channel_length, int slope_window_size, LSDFlowInfo& FlowInfo, LSDRaster& Elevation, string csv_filename);


/// @brief Return a map of all nodes in the channel, useful to just chek if me node is a CNode or not
/// @param FlowInfo
/// @author B.G.
/// @date 12/11/2018
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
map<int,bool> GetMapOfChannelNodes(LSDFlowInfo& flowinfo);

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

  ///A map of strings for holding georeferencing information
  map<string,string> GeoReferencingStrings;

  // Junction information
  /// The number of junctions
  int NJunctions;

  /// @brief A list of the node indices that are sources.
  ///
  /// @details Note: are not the node indices not junctions.
  /// To find the junctions you need to use the get_Junction_of_Node member function
  vector<int> SourcesVector;

  /// This vector lists the node index of the base level nodes that have a source within their catchements.
  vector<int> BaseLevelJunctions;

  /// A list of the junctions. It is an index into the nodevector from the FlowInfo object.
  vector<int> JunctionVector;

  /// @brief The stream order of the junction node/link.
  ///
  ///@details Note that each junction has one and only one receiver junction so
  ///the stream order of a junction node will apply to all nodes along the path
  ///to the next junction.
  vector<int> StreamOrderVector;

  /// A vector that give the baselevel index of each junction node.
  vector<int> BLBasinVector;

  /// Stores the number of donors to each junction.
  vector<int> NDonorsVector;

  /// Stores the node index of the receiving junction.
  vector<int> ReceiverVector;

  /// Stores the delta vector which is used to index into the donor stack
  /// and order contributing junction see Braun and Willett [2012].
  vector<int> DeltaVector;

  /// This is a vector that stores the donor junction of of the junction and is indexed by the DeltaVector.
  vector<int> DonorStackVector;

  ///@brief This vector is used to calculate flow accumulation.
  ///
  ///@details For each base level junction it progresses from a hilltop to a confluence
  ///and then jumps to the next hilltop so that by cascading down through
  ///the node indices in this list one can quickly calculate drainage area,
  ///discharge, sediment flux, etc.
  vector<int> SVector;

  /// This points to the starting point in the S vector of each node.
  vector<int> SVectorIndex;

  /// @brief The number of contributing junctions !!INCULDING SELF!! to a current pixel.
  ///
  ///@details It is used in conjunction with the SVectorIndex to build basins
  /// upslope of any and all nodes in the junction list.
  vector<int> NContributingJunctions;

  // the following arrays are for keeping track of the junctions. For large DEMs this will be quite memory intensive
  // it might be sensible to try to devise a less data intensive method in the future.
  // one could do it with much less memory but that would involve searching

  /// This array stores the stream indices of all the channels.
  Array2D<int> StreamOrderArray;

  /// @brief This array stores a junction counter.
  ///
  /// @details If zero there is no junction \n
  /// if 1 it is a junction unvisted by the junction gathering algorithm \n
  /// if 2 or more it is a previously visited junction
  Array2D<int> JunctionArray;

  /// This is an array where the elements are nodata if there is no junction
  /// and an integer indicating the junction number.
  Array2D<int> JunctionIndexArray;

  private:
  void create( void );
  void create(vector<int> Sources, LSDFlowInfo& FlowInfo);
};

#endif
