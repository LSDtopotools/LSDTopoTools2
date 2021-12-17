//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo
// Land Surface Dynamics FlowInfo
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for organizing flow routing under the Fastscape algorithm
//  (see Braun and Willett, Geomorphology 2013, v180, p 170-179)
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


/** @file LSDFlowInfo.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 1.0
@brief Object to perform flow routing.
@details This is a data object which generates and then
stores information about flow routing.
It is the object that is used to generate contributing area, etc.

@date 29/08/2012
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDFlowInfo_H
#define LSDFlowInfo_H

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
using namespace std;
using namespace TNT;


// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif



/// @brief Object to perform flow routing.
class LSDFlowInfo
{
  public:
  /// @brief The create function. This is default and throws an error.
  /// @author SMM
  /// @date 01/016/12
  LSDFlowInfo()        { create(); }

  /// @brief Creates a FlowInfo object from a binary flowinfo data.
  /// @param fname String of the binary flowinfo data file to be read.
  /// @author SMM
  /// @date 01/016/12
  LSDFlowInfo(string fname)    { create(fname); }

  /// @brief Creates a FlowInfo object from a binary flowinfo data.
  /// @detail This assumes no flux boundaries
  /// @param TopoRaster LSDRaster object containing the topographic data.
  /// @author SMM
  /// @date 3/07/2015
  LSDFlowInfo(LSDRaster& TopoRaster)
                   { create(TopoRaster); }

  /// @brief Creates a FlowInfo object from topography.
  /// @param BoundaryConditions Vector<string> of the boundary conditions at each edge of the
  /// DEM file. Boundary conditions can start with 'P' or 'p' for periodic,
  /// 'B' or 'b' for base level, or anything else for no flux.
  /// the vector shold have 4 elements, 0 is north, 1 is east, 2 is south and 3 is west
  /// @param TopoRaster LSDRaster object containing the topographic data.
  /// @author SMM
  /// @date 01/016/12
  LSDFlowInfo(vector<string>& BoundaryConditions, LSDRaster& TopoRaster)
                   { create(BoundaryConditions, TopoRaster); }

  /// @brief Copy of the LSDJunctionNetwork description here when written.
  friend class LSDJunctionNetwork;

  // some functions for retrieving information out of the data vectors
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

  /// @brief a function to get the lat and long of a location provided as northing and easting
  /// @detail Assumes WGS84 ellipsiod
  /// @param X the Easting of the location
  /// @param Y the Northing of the location
  /// @param lat the latitude of the node (in decimal degrees, replaced by function)
  ///  Note: this is a double, because a float does not have sufficient precision
  ///  relative to a UTM location (which is in metres)
  /// @param long the longitude of the node (in decimal degrees, replaced by function)
  ///  Note: this is a double, because a float does not have sufficient precision
  ///  relative to a UTM location (which is in metres)
  /// @param Converter a converter object (from LSDShapeTools)
  /// @author MDH
  /// @date 27/7/2017
  void  get_lat_and_long_locations(double X, double Y, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter);

  /// @brief a function to get the northing and easting provided as lat and long
  /// @detail Assumes WGS84 ellipsiod
  /// @param lat the latitude of the node (in decimal degrees, replaced by function)
  ///  Note: this is a double, because a float does not have sufficient precision
  ///  relative to a UTM location (which is in metres)
  /// @param long the longitude of the node (in decimal degrees, replaced by function)
  ///  Note: this is a double, because a float does not have sufficient precision
  /// @param X the Easting of the location
  /// @param Y the Northing of the location
  /// @author SMM
  /// @date 30/09/2019
  void get_x_and_y_from_latlong(double lat, double longitude, double& X, double& Y );

  /// @brief this check to see if a point is within the raster
  /// @param X_coordinate the x location of the point
  /// @param Y_coordinate the y location of the point
  /// @return is_in_raster a boolean telling if the point is in the raster
  /// @author SMM
  /// @date 13/11/2014
  bool check_if_point_is_in_raster(float X_coordinate,float Y_coordinate);

  ///@brief Gives the receiver information for a given node.
  ///@param current_node Integer
  ///@param reveiver_node Empty integer to be assigned the index of the receiver
  ///node.
  ///@param receiver_row Empty integer to be assigned the row index of the
  ///receiver node.
  ///@param receiver_col Empty integer to be assigned the column index of the
  ///receiver node.
  /// @author SMM
  /// @date 01/016/12
  void retrieve_receiver_information(int current_node,int& reveiver_node, int& receiver_row,
                                             int& receiver_col);
  ///@brief Gives the receiver information for a given node.
  ///@param current_node Integer
  ///@param reveiver_node Empty integer to be assigned the index of the receiver
  ///node.
  /// @author BG
  /// @date 06/01/2018
  void retrieve_receiver_information(int current_node,int& reveiver_node);

  ///@brief Get the row and column indices of a given node.
  ///@param current_node Integer index of a given node.
  ///@param curr_row Empty integer to be assigned the row index of the given
  ///node.
  ///@param curr_col Empty integer to be assigned the column index of the given
  ///node.
  /// @author SMM
  /// @date 01/016/12
  void retrieve_current_row_and_col(int current_node,int& curr_row,
                                             int& curr_col);

  ///@brief Get the index from row/col
  ///@param int row/col of the nodeindex 2DArray
  ///@author BG
  ///@date 27/12/2017
  int get_NodeIndex_from_row_col(int row, int col);

  ///@brief Get the X and Y coordinates of a given node.
  ///@param current_node Integer index of a given node.
  ///@param current_X Empty integer to be assigned the X coordinate of the given
  ///node.
  ///@param current_Y Empty integer to be assigned the Y coordinate of the given
  ///node.
  /// @author BG
  /// @date 20/02/2017
  void get_x_and_y_from_current_node(int current_node,float& current_X, float& current_Y);


  ///@brief Get the lat and longitude coordinates of a given node.
  ///@param current_node Integer index of a given node.
  ///@param current_lat latitude. Will be replaced by function
  ///@param current_long longitude. Will be replaced by function
  ///@param Converter A coordinate converter object
  /// @author SMM
  /// @date 26/04/2017
  void get_lat_and_long_from_current_node(int current_node, double& current_lat, double& current_long, LSDCoordinateConverterLLandUTM Converter);



  ///@brief This function takes a vector of node indices and prints a csv
  ///file that can be read by arcmap
  ///@param nodeindex vec is a vector of nodeindices (which are ints)
  ///@param outfilename is a string of the filename
  /// @author SMM
  /// @date 03/06/14
  void print_vector_of_nodeindices_to_csv_file(vector<int>& nodeindex_vec, string outfilename);

  ///@brief This function takes a vector of node indices and prints a csv
  ///file that can be read by arcmap: similar to above but also prints lat and long
  ///@param nodeindex vec is a vector of nodeindices (which are ints)
  ///@param outfilename is a string of the filename
  /// @author SMM
  /// @date 20/05/16
  void print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int>& nodeindex_vec, string outfilename);


  /// @brief This function takes a vector of node indices and prints a csv
  ///  file that can be read by arcmap: similar to above but also prints lat and long
  /// @detail Adds an option for a path
  /// @param nodeindex vec is a vector of nodeindices (which are ints)
  /// @param path to the outfile
  /// @param outfilename is a string of the filename
  /// @author SMM
  /// @date 20/05/16
  void print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int> node_list,string path, string filename);


  /// @brief This function takes a vector of node indices and prints a csv
  ///  file that can be read by arcmap: similar to above but also prints lat and long
  /// @detail Adds an option for a path
  /// @param nodeindex vec is a vector of nodeindices (which are ints)
  /// @param path to the outfile
  /// @param outfilename is a string of the filename
  /// @param Elevation an elevation raster
  /// @param FlowDistance
  /// @param drainage_area
  /// @author SMM
  /// @date 30/09/19
  void print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int> node_list,string path, string filename, LSDRaster& Elevation, LSDRaster& FlowDistance,
                                                            LSDRaster& drainage_area);

  ///@brief This function takes a vector of node indices and prints a csv
  ///file that can be read by arcmap, adding in a unique id to each row, independent of the nodeindex.
  ///
  ///@details The unique ID is used to tie triplets of channel heads together for hollow analysis.
  ///@param nodeindex vec is a vector of nodeindices (which are ints)
  ///@param outfilename is a string of the filename
  ///@author SWDG after SMM
  ///@date 2/2/16
  void print_vector_of_nodeindices_to_csv_file_Unique(vector<int>& nodeindex_vec, string outfilename);

  ///@brief This function takes a vector of node indices and moves
  ///  down the list, getting the bearings of each segment subject to an offset
  ///@param nodeindex vec is a vector of nodeindices (which are ints)
  ///@param node_spacing how many nodes between the two nodevecs
  ///@param print_bearings if true, prints to a geojson
  ///@param outfileprefix the prefix of the outfile (.geojson gets added)
  ///@param bearing_nodes The list of nodes from which the bearing is taken
  ///@param bearings the actual bearings
  ///@author SMM
  ///@date 4/3/21
  void calculate_bearings_from_nodelist(vector<int> nodeindex_vec, int node_spacing,
                                                 bool print_bearings, string outfileprefix,
                                                 vector<int>& bearing_nodes, vector<float>& bearing_vec);


  ///@brief Get the number of pixels flowing into a node.
  ///@param node Integer of node index value.
  ///@return Integer of the number of contributing pixels.
  /// @author SMM
  /// @date 01/016/12
  int retrieve_contributing_pixels_of_node(int node)
                  { return NContributingNodes[node]; }

  ///@brief Get the FlowLengthCode of a given node.
  ///@param node Integer of node index value.
  ///@return Integer of the FlowLengthCode.
  /// @author SMM
  /// @date 01/016/12
  int retrieve_flow_length_code_of_node(int node)
               { return FlowLengthCode[ RowIndex[node] ][ ColIndex[node] ]; }

  ///@brief Get the FlowDirection of a row and column pair.
  ///@param row Integer of row index.
  ///@param col Integer of col index.
  ///@return Integer of the flow direction.
  ///@author SWDG
  ///@date 04/02/14
  int get_LocalFlowDirection(int row, int col)
               { return FlowDirection[row][col]; }

  /// @brief get the number of donors to a given node
  /// @param current_node the node index from which to get n donors
  /// @return the number of donors to this cell
  /// @author SMM
  /// @date 19/9/2014
  int retrieve_ndonors_to_node(int current_node)
    { return NDonorsVector[current_node]; }

  ///@brief Get the node for a cell at a given row and column
  ///@param row index
  ///@param column index
  ///@return node index
  ///@author DTM
  ///@date 08/11/2013
  int retrieve_node_from_row_and_column(int row, int column);

  /// @brief gets a vector of all the donors to a given node
  /// @param current_node the node index from which to get n donors
  /// @return a vector containing all the donors to this node
  /// @author SMM
  /// @date 19/9/2014
  vector<int> retrieve_donors_to_node(int current_node);

  /// @brief returns the draiange area of a node in square km
  /// @param this_node node of interest
  /// @return draiange area in square km
  /// @author FJC
  /// @date 06/02/17
  float get_DrainageArea_square_km(int this_node);

  /// @brief returns the draiange area of a node in square m
  /// @param this_node node of interest
  /// @return draiange area in square km
  /// @author FJC
  /// @date 01/05/18
  float get_DrainageArea_square_m(int this_node);

  // get functions

  /// @return Number of rows as an integer.
  int get_NRows() const                   { return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const                   { return NCols; }
  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const              { return XMinimum; }
  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const              { return YMinimum; }
  /// @return Data resolution as an integer.
  float get_DataResolution() const        { return DataResolution; }
  /// @return No Data Value as an integer.
  int get_NoDataValue() const             { return NoDataValue; }
  /// @return Georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }

  /// @return Number of nodes with data as an integer.
  int get_NDataNodes () const          { return NDataNodes; }
  /// @return Vector of all base level nodes.
  vector<int> get_BaseLevelNodeList () { return BaseLevelNodeList; }
  /// @return get the number of baselevel nodes
  int get_NBaseLevelNodes () { return int( BaseLevelNodeList.size() );}

  /// @return donor stack vector (depth first search sequence of nodes)
  vector <int> get_donorStack() const { return DonorStackVector; }
  /// @return the S vector, which is a sorted list of nodes (see Braun and Willett 2012)
  vector <int> get_SVector() const { return SVector; }
  /// @return FlowDirection values as a 2D Array.
  Array2D<int> get_FlowDirection() const { return FlowDirection; }
  int get_index_SVector(int node){return SVectorIndex[node];}

  ///@brief Recursive add_to_stack routine, from Braun and Willett (2012)
  ///equations 12 and 13.
  ///@param lm_index Integer
  ///@param j_index Integer
  ///@param bl_node Integer
  void add_to_stack(int lm_index, int& j_index, int bl_node);

  // some functions that print out indices to rasters
  ///@brief Write NodeIndex to an LSDIndexRaster.
  ///@return LSDIndexRaster of node index data.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster write_NodeIndex_to_LSDIndexRaster();


  ///@brief Write FlowDirection to an LSDIndexRaster.
  ///@return LSDIndexRaster of flow directions.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster write_FlowDirection_to_LSDIndexRaster();


  ///@brief Write FlowLengthCode to an LSDIndexRaster.
  ///@return LSDIndexRaster of flow lengths.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster write_FlowLengthCode_to_LSDIndexRaster();

  /// @brief This function writes and LSDIndexRaster containing the location of nodes in the nodeindexvector.
  /// @param nodeindexvec a vector containing node indices one use is to export
  /// the LSDIndexRaster of pixels that are in the node index vector.
  /// @return LSDIndexRaster of pixels that are in the node index vector.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster write_NodeIndexVector_to_LSDIndexRaster(vector<int>& nodeindexvec);

  /// @brief This function writes an LSDIndesxRaster given a list of node indices, and give every pixel its nodeindex value, which is unique.
  /// @param nodeindexvec a vector containing node indices one use is to export
  /// the LSDIndexRaster of pixels that are in the node index vector.
  /// @return LSDIndexRaster of pixels that are in the node index vector.
  /// @author SWDG after SMM
  /// @date 2/2/16
  LSDIndexRaster write_NodeIndexVector_to_LSDIndexRaster_Unique(vector<int>& nodeindexvec);

  ///@brief Write NContributingNodes to an LSDIndexRaster.
  ///@return LSDIndexRaster of number of contributing nodes for each cell.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster write_NContributingNodes_to_LSDIndexRaster();

  ///@brief Writes flow directions to an LSDIndexRaster.
  ///Flow direction in arcmap format is: \n\n
  /// 32  64  128 \n
  /// 16   0  1    \n
  /// 8    4  2     \n
  ///
  ///@return LSDIndexRaster of flow directions in arcgis format.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster write_FlowDirection_to_LSDIndexRaster_Arcformat();

  ///@brief
  ///@return
  ///@author Fiona Clubb
  ///@date 15/11/12
  LSDRaster write_DrainageArea_to_LSDRaster();

  ///@brief Prints the flow information to file.
  ///@param filename String of the output file to be written.
  /// @author SMM
  /// @date 01/016/12
  void print_flow_info_vectors(string filename);

  ///@brief Unpickles flow information data from a binary file.
  ///@param filename String of the binary file to be read.
  /// @author SMM
  /// @date 01/016/12
  void unpickle(string filename);

  ///@brief Pickles flow information data from a binary file.
  ///@details WARNING!!! This creates HUGE files (sometimes 10x bigger than
  /// original file). Testing indicates reading this file takes
  /// almost as long as recalculating the flowinfo object so
  /// is probably not worth doing
  ///@param filename String of the binary file to be written.
  /// @author SMM
  /// @date 01/016/12
  void pickle(string filename);

  /// @brief This loads a csv file, putting the data into a data map
  /// @param filename The name of the csv file including path and extension
  /// @author SMM (ported into FlowInfo FJC 23/03/17)
  /// @date 16/02/2017
  map<string, vector<string> > load_csv_data(string filename);

  /// @brief This gets a data column from the csv file
  /// @param column_name a string that holds the column name
  /// @return a vector of strings: this holds the data.
  /// @author SMM (ported into FlowInfo FJC 23/03/17)
  /// @date 17/02/2017
  vector<string> get_data_column(string column_name, map<string, vector<string> > data_map);

  /// @brief This gets a data column from the csv file, and converts it to a
  ///   float vector
  /// @param column_name a string that holds the column name
  /// @return a vector of floats: this holds the data.
  /// @author SMM (ported into FlowInfo FJC 23/03/17)
  /// @date 17/02/2017
  vector<float> data_column_to_float(string column_name, map<string, vector<string> > data_map);

  /// @brief This gets a data column from the csv file, and converts it to an
  ///   int vector
  /// @param column_name a string that holds the column name
  /// @return a vector of ints: this holds the data.
  /// @author SMM (ported into FlowInfo FJC 23/03/17)
  /// @date 17/02/2017
  vector<int> data_column_to_int(string column_name, map<string, vector<string> > data_map);

  /// @brief Method to ingest the channel heads raster generated using channel_heads_driver.cpp
  /// into a vector of source nodes so that an LSDJunctionNetwork can be created easily
  /// from them.  **UPDATE** if the extension is a csv file it reads the node indices directly
  /// **UPDATE, FJC 20/01/16 - changed default input switch to 2**
  /// @details Assumes the FlowInfo object has the same dimensions as the channel heads raster.
  /// @param filename of the channel heads raster.
  /// @param extension of the channel heads raster.
  /// @param (optional) input_switch, ONLY NEEDED FOR LOADING .csv FILES!
  ///  An integer to determine whether to use the node index (0 -> default), row and column indices (1),
  ///  or point coordinates from .csv file (2) to locate the channel heads
  /// @return Vector of source nodes.
  /// @author SWDG updated SMM updated DTM
  /// @date 6/6/14 Happy 3rd birthday Skye!!
  vector<int> Ingest_Channel_Heads(string filename, string extension, int input_switch = 2);

  /// @brief Method to ingest the channel heads raster generated using channel_heads_driver.cpp
  /// into a vector of source nodes so that an LSDJunctionNetwork can be created easily
  /// from them.  **UPDATE** if the extension is a csv file it reads the node indices directly
  /// **UPDATE, FJC 20/01/16 - changed default input switch to 2**
  /// ***********UPDATE - overloaded function to read in the column headers using the csv logic
  /// rather than assume that the columns have to be in a specific order. ONLY WORKS WITH CSV
  /// EXTENSIONS. FJC 23/03/17*************************
  /// @details Assumes the FlowInfo object has the same dimensions as the channel heads raster.
  /// @param filename of the channel heads raster.
  /// @param extension of the channel heads raster.
  /// @param (optional) input_switch, ONLY NEEDED FOR LOADING .csv FILES!
  ///   An integer to determine whether to use the node index (0 -> default),
  ///   row and column indices (1), or point coordinates from .csv file (2) to locate the channel heads
  /// @return Vector of source nodes.
  /// @author SWDG updated SMM updated DTM
  /// @date 6/6/14 Happy 3rd birthday Skye!!
  vector<int> Ingest_Channel_Heads(string filename, int input_switch = 2);

	/// @brief Method to ingest sources from OS MasterMap Water Network Layer (csv) into a vector
	/// of source nodes so that an LSDJunctionNetwork can be easily created from them.
	/// @param csv_filename CSV file name
	/// @return vector of source nodes
  /// @author FJC
  /// @date 28/11/16
  vector<int> Ingest_Channel_Heads_OS(string csv_filename);

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Minimalistic method to ingest the channel heads from vectors of x y coordinates
  // Using xy allows a "universal" method that can ingest external or internal data
  // B.G. 11/11/2018
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  vector<int> Ingest_Channel_Heads(vector<float>& x_coord, vector<float>& y_coord);

  // functions for getting flow, discharge, sediment flux, etc

  ///@brief This function calculates the contributing pixels.
  ///It can be converted to contributing area by multiplying by the
  ///DataResolution^2. In this function a pixel that has no donors has a
  ///contributing pixel value of 0.
  ///@return LSDIndexRaster of upslope contributing pixels.
  /// @author SMM
  /// @date 01/016/12
  LSDIndexRaster calculate_n_pixels_contributing_from_upslope();

  ///@brief This calculates area and makes an index into the s vector for
  ///efficient calculation of the basin upslope of a given node.
  /// @author SMM
  /// @date 01/016/12
  void calculate_upslope_reference_indices();

  // algorithms for basin collection
  ///@brief This function returns the base level node with the greatest
  ///drainage area.
  ///@return Integer node index.
  int retrieve_largest_base_level();

  /// @brief This gets the base level node for any given node
  /// @param node the starting node
  /// @return the base level node
  /// @author SMM
  /// @date 18/12/2016
  int retrieve_base_level_node(int node);

  /// @brief This creates a new raster that has any nodes influenced by the edge removed
  /// @detail uses D8 flow routing to do this
  /// @param topography A topography raster (must be same dimension as that used to make flowinfo)
  /// @author SMM
  /// @date 08/10/2019
  LSDRaster remove_nodes_influneced_by_edge(LSDRaster& topography);

  /// @brief This creates a new raster. It takes a node list (usually a channel) and then
  ///  looks at all the donors to that channel and sees if any are influenced upslope
  ///  by nodata. It then removes all these pixels. This way you can isolate both
  ///  nodes drainaing to a local base level but eliminate nodes draining from the edge.
  /// @detail uses D8 flow routing to do this
  /// @param topography A topography raster (must be same dimension as that used to make flowinfo)
  /// @param node_list a list of nodes (often read by an LSDSpatialCSVReader object)
  /// @author SMM
  /// @date 08/10/2019
  LSDRaster find_nodes_not_influenced_by_edge_draining_to_nodelist(vector<int> node_list,LSDRaster& topography);

  ///@brief This function returns an integer vector containing all the node
  ///indexes upslope of of the node with number node_number_outlet.
  ///@param node_number_outlet Integer of the target node.
  ///@return Integer vector of upslope node indexes.
  /// @author SMM
  /// @date 01/016/12
  vector<int> get_upslope_nodes(int node_number_outlet);

  ///@brief This function returns an integer vector containing all the node
  ///indexes upslope of of the node with number node_number_outlet, and then includes the outlet
  ///@param node_number_outlet Integer of the target node.
  ///@return Integer vector of upslope node indexes plus outlet node index.
  /// @author SMM
  /// @date 30/07/19
  vector<int> get_upslope_nodes_include_outlet(int node_number_outlet);

  /// @brief This function takes a list of sources and then creates a raster
  ///  with nodata values where points are not upslope of the sources
  ///  and 1.0 if they are upslope
  /// @param source_nodes a vector of node indicies into the sources
  /// @author SMM
  /// @date 11/11/2015
  LSDRaster get_upslope_node_mask(vector<int> source_nodes);

  /// @brief This function takes a list of sources and then creates a raster
  ///  with nodata values where points are not upslope of the sources
  ///  and node_values if they are upslope
  /// @param source_nodes a vector of node indicies into the sources
  /// @param node_values a vector of the value of the nodes upslope of the sources
  ///  this vector needs to be the same size as the source_nodes vector
  /// @author SMM
  /// @date 12/11/2015
  LSDRaster get_upslope_node_mask(vector<int> source_nodes, vector<float> node_values);

  ///@brief This function accumulates some variable from an LSDRaster
  ///The most probably use is to accumulate precipitation in order
  ///to get a discharge raster
  ///@param A raster that contains the variable to be accumulated (e.g., precipitation)
  ///@return A raster containing the accumulated variable: NOTE the accumulation
  ///does not include the node itself
  ///@author SMM
  ///@date 09/06/2014
  LSDRaster upslope_variable_accumulator(LSDRaster& accum_raster);

  ///@brief This function accumulates some variable from an LSDRaster
  ///The most probably use is to accumulate precipitation in order
  ///to get a discharge raster
  ///@param A raster that contains the variable to be accumulated (e.g., precipitation)
  ///@param bool Determinde if the current node is accumulated to itself or not. Default to yes.
  ///@return A raster containing the accumulated variable
  ///@author BG, SMM
  ///@date 17/04/2019
  LSDRaster upslope_variable_accumulator_v2(LSDRaster& accum_raster);
  LSDRaster upslope_variable_accumulator_v2(LSDRaster& accum_raster, bool accum_current_node);

  ///@brief This function accumulates some variable from an LSDRaster
  ///The most probably use is to accumulate precipitation in order
  ///to get a discharge raster
  ///@detail Uses the stack vector to get the accumulation
  ///@param A raster that contains the variable to be accumulated (e.g., precipitation)
  ///@return A raster containing the accumulated variable
  ///@author SMM
  ///@date 19/11/2019
  LSDRaster upslope_variable_accumulator_v3(LSDRaster& accum_raster);

  ///@brief This function returns the upslope average value of 
  /// some raster bassed on the flow info
  ///@detail Uses the stack vector to get the accumulation
  ///@param A raster that contains the variable to be averaged (e.g., precipitation)
  ///@return A raster containing the uplsope average
  ///@author SMM
  ///@date 25/06/2021
  LSDRaster upslope_average(LSDRaster& accum_raster);


  ///@brief This function tests whether one node is upstream of another node
  ///@param current_node
  ///@param test_node The node which we want to test whether it is upstream or not
  ///@return Boolean indicating whether node is upstream or not
  ///@author FC
  ///@date 08/10/13
  int is_node_upstream(int current_node, int test_node);

  ///@brief This function tests whether a node is a base level node
  ///@param node
  ///@return int which is 1 if node is base level, 0 if not
  ///@author FJC
  ///@date 26/08/15
  int is_node_base_level(int node);

  /// @brief this function gets a list of the node indices of the donors to a particular node
  /// @param node this is the nodeindex of the node for which you want to find the donors
  /// @return a vector of the donor nodes
  /// @author SMM
  /// @date 21/10/2013
  vector<int> get_donor_nodes(int node);


// algorithms for stream profile analysis

  /// @brief This function calculates the chi function for all the nodes upslope
  /// of a given node.
  /// @param starting_node Integer index of node to analyse upslope of.
  /// @param m_over_n
  /// @param A_0
  /// @return Vector of chi values. The node indices of these values are those
  /// that would be retured from get_uplope_nodes
  /// @author SMM
  /// @date 01/16/2012
  vector<float> get_upslope_chi(int starting_node, float m_over_n, float A_0);

  /// @brief This function calculates the chi function for all the nodes upslope
  /// of a given node.
  /// @detail this version uses discharge rather than area
  /// @param starting_node Integer index of node to analyse upslope of.
  /// @param m_over_n
  /// @param A_0  the referecen discharge
  /// @param Discharge a raster of the discharge
  /// @return Vector of chi values. The node indices of these values are those
  /// that would be retured from get_uplope_nodes
  /// @author SMM
  /// @date 16/10/2015
  vector<float> get_upslope_chi(int starting_node, float m_over_n, float A_0, LSDRaster& Discharge);

  /// @brief This function calculates the chi function for a list of nodes
  /// it isn't really a standalone modules, but is only called from get_upslope_chi
  /// above
  /// @param upslope_pixel_list Vector of nodes to analyse.
  /// @param m_over_n
  /// @param A_0
  /// @return Vector of chi values.
  /// @author SMM
  /// @date 16/01/12
  vector<float> get_upslope_chi(vector<int>& upslope_pixel_list, float m_over_n, float A_0);

  /// @brief This function calculates the chi function for a list of nodes
  /// it isn't really a standalone modules, but is only called from get_upslope_chi
  /// above
  /// @detail this version uses discharge rather than area
  /// @param upslope_pixel_list Vector of nodes to analyse.
  /// @param m_over_n
  /// @param A_0
  /// @param Discharge and LSDRaster of the discharge
  /// @return Vector of chi values.
  /// @author SMM
  /// @date 16/10/15
  vector<float> get_upslope_chi(vector<int>& upslope_pixel_list, float m_over_n, float A_0, LSDRaster& Discharge);




  /// @brief This function calculates the chi function for a list of nodes
  /// it isn't really a standalone modules, but is only called from get_upslope_chi
  /// above. It returns a map, which is used to speed up computation
  /// @param upslope_pixel_list Vector of nodes to analyse.
  /// @param m_over_n
  /// @param A_0
  /// @param mimum_pixels This minimum number of contributing pixels needed before chi is calculated
  /// @return A map where the key is the node index and the value is chi
  /// @author SMM
  /// @date 13/07/17
  map<int,float> get_upslope_chi_return_map(vector<int>& upslope_pixel_list, float m_over_n, float A_0, int minimum_pixels);

  /// @brief This function calculates the chi function for a list of nodes
  /// it isn't really a standalone modules, but is only called from get_upslope_chi
  /// above. It returns a map, which is used to speed up computation
  /// @detail this version uses discharge rather than area
  /// @param upslope_pixel_list Vector of nodes to analyse.
  /// @param m_over_n
  /// @param A_0
  /// @param mimum_pixels This minimum number of contributing pixels needed before chi is calculated
  /// @param Discharge and LSDRaster of the discharge
  /// @return A map where the key is the node index and the value is chi
  /// @author SMM
  /// @date 13/07/17
  map<int,float> get_upslope_chi_return_map(vector<int>& upslope_pixel_list, float m_over_n, float A_0, int minimum_pixels, LSDRaster& Discharge);






  /// @brief this function takes a vector that contains the node indices of
  /// starting nodes, and then calculates chi upslope of these nodes
  /// to produce a chi map. A threshold drainage area can be used
  /// to only map chi where nodes have greater than the threshold drainage area
  /// @details this function is meant to mimic the function of the Willett et al
  /// (2014) Science paper. You do need to extract the wanted node indices
  /// for your starting nodes from a node index map
  /// @param starting_nodes an integer vector containing all the node indices
  /// of the node from which you want to start the chi analysis. All of these
  /// nodes will be considered to have a starting chi of 0
  /// @param m_over_n the m/n ratio. Chi is quite sensitive to this
  /// @param A_0 the reference drainage area. This is a but arbitrary. We usually use
  /// 1000 m^2. Willet et al(2014, Science) used 1m^2. As of 28 July 2014 we've not
  /// done any detailed sensitivity analysis on this parameter
  /// @param area_threshold the threshold area (in m^2) that sets the area above
  /// which chi is recorded in the chi raster
  /// @return this returns an LSDRaster for the chi values upslope of all of the
  /// nodes indicated in starting_nodes
  /// @author SMM
  /// @date 28/14/2014
  LSDRaster get_upslope_chi_from_multiple_starting_nodes(vector<int>& starting_nodes,
   float m_over_n, float A_0, float area_threshold);

  /// @brief this function takes a vector that contains the node indices of
  /// starting nodes, and then calculates chi upslope of these nodes
  /// to produce a chi map. A threshold drainage area can be used
  /// to only map chi where nodes have greater than the threshold drainage area
  /// @details this function is meant to mimic the function of the Willett et al
  /// (2014) Science paper. You do need to extract the wanted node indices
  /// for your starting nodes from a node index map
  /// This version allows computation with discharge
  /// @param starting_nodes an integer vector containing all the node indices
  /// of the node from which you want to start the chi analysis. All of these
  /// nodes will be considered to have a starting chi of 0
  /// @param m_over_n the m/n ratio. Chi is quite sensitive to this
  /// @param A_0 the reference discharge. This is arbitrary. As of 28 July 2014 we've not
  /// done any detailed sensitivity analysis on this parameter
  /// @param area_threshold the threshold area (in m^2) that sets the area above
  /// which chi is recorded in the chi raster
  /// @return this returns an LSDRaster for the chi values upslope of all of the
  /// nodes indicated in starting_nodes
  /// @author SMM
  /// @date 16/10/2015
  LSDRaster get_upslope_chi_from_multiple_starting_nodes(vector<int>& starting_nodes,
   float m_over_n, float A_0, float area_threshold, LSDRaster& Discharge);

  /// @brief Same than above, but ignore the threshold. Discharge can be any raster custom drainage area/discharge-esc raster
  /// @author BG/SMM
  /// @date 16/10/2015
  LSDRaster get_upslope_chi_from_multiple_starting_nodes_custom(vector<int>& starting_nodes,
   float m_over_n, float A_0,  LSDRaster& Discharge);

  /// @brief This funtion gets all the upslope chi of a starting node (assuming
  ///  chi at starting node is 0) and returns a map
  /// @param starting_nodes an integer containing the node index
  /// of the node from which you want to start the chi analysis.
  /// @param m_over_n the m/n ratio. Chi is quite sensitive to this
  /// @param A_0 the reference discharge.
  /// @param mimum_pixels This minimum number of contributing pixels needed before chi is calculated
  /// @return Returns a map where the key is the node index and the value is chi
  /// @author SMM
  /// @date 13/07/2017
    map<int,float> get_upslope_chi_from_single_starting_node(int starting_node,
                                 float m_over_n, float A_0, int minimum_pixels);


  /// @brief This funtion gets all the upslope chi of a starting node (assuming
  ///  chi at starting node is 0) and returns a map
  /// @details This version allows computation with discharge
  /// @param starting_nodes an integer containing the node index
  /// of the node from which you want to start the chi analysis.
  /// @param m_over_n the m/n ratio. Chi is quite sensitive to this
  /// @param A_0 the reference discharge.
  /// @param mimum_pixels This minimum number of contributing pixels needed before chi is calculated
  /// @param Discharge The discharge raster
  /// @return Returns a map where the key is the node index and the value is chi
  /// @author SMM
  /// @date 13/07/2017
    map<int,float> get_upslope_chi_from_single_starting_node(int starting_node,
                                 float m_over_n, float A_0, int minimum_pixels, LSDRaster& Discharge);


  /// @brief This function gets the chi upslope of every base level node
  /// that is, it gets the chi values of the entire DEM, assuming all
  /// base level nodes have a chi of 0
  /// @detail because this assumes all base level nodes have a chi of 0,
  /// this function is probably only appropriate for numerical models.
  /// @param m_over_n the m/n ratio. Chi is quite sensitive to this
  /// @param A_0 the reference drainage area. This is a but arbitrary. We usually use
  /// 1000 m^2. Willet et al(2014, Science) used 1m^2. As of 28 July 2014 we've not
  /// done any detailed sensitivity analysis on this parameter
  /// @param area_threshold the threshold area (in m^2) that sets the area above
  /// which chi is recorded in the chi raster
  /// @return this returns an LSDRaster for the chi values of the entire raster,
  /// with base level nodes assumed to have chi = 0
  /// @author SMM
  /// @date 28/14/2014
  LSDRaster get_upslope_chi_from_all_baselevel_nodes(float m_over_n, float A_0,
                                                float area_threshold);

  /// @brief This function gets the chi upslope of every base level node
  /// that is, it gets the chi values of the entire DEM, assuming all
  /// base level nodes have a chi of 0
  /// @detail because this assumes all base level nodes have a chi of 0,
  /// this function is probably only appropriate for numerical models.
  /// This version of the function allows computation with a discharge raster
  /// @param m_over_n the m/n ratio. Chi is quite sensitive to this
  /// @param A_0 the reference discharge (same units as discharge. As of 28 July 2014 we've not
  /// done any detailed sensitivity analysis on this parameter
  /// @param area_threshold the threshold area (in m^2) that sets the area above
  /// which chi is recorded in the chi raster
  /// @param Discharge a raster of the discharge
  /// @return this returns an LSDRaster for the chi values of the entire raster,
  /// with base level nodes assumed to have chi = 0
  /// @author SMM
  /// @date 16/10/2015
  LSDRaster get_upslope_chi_from_all_baselevel_nodes(float m_over_n, float Q_0,
                                                float area_threshold,
                                                LSDRaster& Discharge);

  /// @brief Calculates the distance from outlet of all the base level nodes.
  /// Distance is given in spatial units, not in pixels.
  /// @return LSDRaster of the distance to the outlet for all baselevel nodes.
  /// @author SMM
  /// @date 01/016/12
  LSDRaster distance_from_outlet();


  /// @biref calculates the slope measured in the d8 flow direction
  /// base level nodes have slope of 0; slope is measured from node to receiver
  /// @param the elevation raster
  /// @return  A raster of the d8 slope
  /// @author SMM
  /// @date 21/09/2014
  LSDRaster calculate_d8_slope(LSDRaster& Elevation);

  /// @brief This returns the node index of the pixel farthest upslope from the input node.
  /// @param node the node from which you want to find the farthest upslope pixel.
  /// @param DistFromOutlet an LSDRaster containing the distance from the outlet.
  /// @return This returns the node index of the pixel farthest upslope
  /// from the input node.
  /// @author SMM
  /// @date 25/19/13
  int find_farthest_upslope_node(int node, LSDRaster& DistFromOutlet);

  /// @brief This takes a list of nodes and sorts them according to a
  ///  sorting raster (it could be anything) fin ascending order
  ///  nodes are then reordered to reflect the sorting of the raster
  /// @param node_vec A vector of nodes
  /// @param SortingRaster a raster that contains some values that will
  ///  be used to sort the nodes.
  /// @return sorted_nodes a vector containing the sorted nodes
  /// @author SMM
  /// @date 20/05/2016
  vector<int> sort_node_list_based_on_raster(vector<int> node_vec, LSDRaster& SortingRaster);

  /// @brief This takes a list of nodes and sorts them according to a
  ///  sorting raster (it could be anything) fin ascending order
  ///  nodes are then reordered to reflect the sorting of the raster
  /// @param node_vec A vector of nodes
  /// @param SortingRaster a raster that contains some values that will
  ///  be used to sort the nodes.
  /// @return sorted_nodes a vector containing the sorted nodes
  /// @author SMM
  /// @date 20/05/2016
  vector<int> sort_node_list_based_on_raster(vector<int> node_vec, LSDIndexRaster& SortingRaster);


  /// @brief Function to get the node index for a point using its X and Y coordinates
  /// @param X_coordinate X_coord of point
  /// @param Y_coordinate Y_coord of point
  /// @return int with node index of point
  /// @author FJC
  /// @date 11/02/14
  int get_node_index_of_coordinate_point(float X_coordinate, float Y_coordinate);

  /// @brief Get vector of nodeindices from csv file
  /// @param csv_filename input csv file
  /// @return vector<int> with nodeindices. Ignores nodatavalues
  /// @author FJC
  /// @date 14/02/17
void get_nodeindices_from_csv(string csv_filename, vector<int>& NIs, vector<float>& X_coords, vector<float>& Y_coords);

  /// @brief Function to return the closest value of a raster to a specified node index
  /// @param NodeIndex of interest
  /// @param InputRaster raster to return value of
  /// @param search_radius rectangular window to search for, in n_pixels
  /// @return float value of raster closest to the node index
  /// @author FJC
  /// @date 08/02/17
  float snap_RasterData_to_Node(int NodeIndex, LSDRaster& InputRaster, int search_radius);

  ///@brief A get sources version that uses the flow accumulation pixels.
  ///@param FlowPixels LSDIndexRaster of flow accumulation in pixels.
  ///@param threshold Integer flow accumulation threshold.
  ///@return Vector of source integers: these refer to the node indices of the sources.
  /// @author SMM
  /// @date 01/016/12
  vector<int> get_sources_index_threshold(LSDIndexRaster& FlowPixels, int threshold);

  /// @brief A get sources version that uses AS^2 (area and slope).
  /// @param FlowPixels LSDIndexRaster of flow accumulation in pixels.
  /// @param Slope LSDRaster of slope values
  /// @param threshold Integer AS^2 threshold
  /// @return Vector of source integers: these refer to the node indices of the sources.
  /// @author FJC
  /// @date 11/02/14
  vector<int> get_sources_slope_area(LSDIndexRaster& FlowPixels, LSDRaster& Slope, int threshold);

  /// @brief Gets a vector of source nodes based on the X and Y coordinates of mapped channel
  /// heads.  Can be used if all the channel heads in a basin were mapped to get the stream
  /// network and calculate the drainage density
  /// @param X_coords X coordinates of channel heads
  /// @param Y_coords Y coordinates of channel heads
  /// @return Vector of source nodes
  /// @author FJC
  /// @date 17/02/14
  vector<int> get_sources_from_mapped_channel_heads(vector<float>& X_coords, vector<float>& Y_coords);

  /// @brief Perform a downslope trace using D8 from a given point source (i,j).
  ///
  /// @details Overwrites input parameters to return a raster of the path, the length of the
  /// trace and the final pixel coordinates of the trace.
  /// @param i Row index of starting point for trace.
  /// @param j Column index of starting point for trace.
  /// @param StreamNetwork An LSDIndexRaster of the stream network.
  /// @param length Length of trace in spatial units.
  /// @param receiver_row Row index of ending point for trace.
  /// @param receiver_col Column index of ending point for trace.
  /// @param Path Empty raster to store the final trace path.
  /// @author SWDG
  /// @date 20/1/14
  void D8_Trace(int i, int j, LSDIndexRaster StreamNetwork, float& length,
                   int& receiver_row, int& receiver_col, Array2D<int>& Path);

  /// @brief Perform a downslope trace using D8 from a given point source.
  /// @details This version uses lat-long
  /// @param latitude The latitude of the source.
  /// @param longitude The longitude of the source
  /// @return A vector of nodeindices that contain the nodes from start to finish that are
  ///  on the flow path
  /// @author SMM
  /// @date 30/09/2019
  vector<int> get_flow_path(float latitude, float longitude);


  /// @brief Perform a downslope trace using D8 from a given point source.
  /// @details This version uses a nodeindex
  /// @param ni the nodeindex of the starting point
  /// @return A vector of nodeindices that contain the nodes from start to finish that are
  ///  on the flow path
  /// @author SMM
  /// @date 12/03/2021
  vector<int> get_flow_path(int ni);


  /// @brief Move the location of the channel head downslope by a user defined distance.
  /// @param Sources a vector of node indexes of the channel heads to be moved.
  /// @param MoveDist The distance in spatial units the head is to be moved.
  /// @param DownslopeSources A vector used to contain the node indexes of the moved channel heads.
  /// @param FinalHeads A vector containing a subset of the original channel heads which corresponds to the moved heads.
  /// @author SWDG
  /// @date 27/11/15
  void MoveChannelHeadDown(vector<int> Sources, float MoveDist, vector<int>& DownslopeSources, vector<int>& FinalHeads);

  /// @brief Move the location of the channel head upslope by a user defined distance.
  /// @param Sources a vector of node indexes of the channel heads to be moved.
  /// @param MoveDist The distance in spatial units the head is to be moved.
  /// @param DEM the elevation data.
  /// @param UpslopeSources A vector used to contain the node indexes of the moved channel heads.
  /// @param FinalHeads A vector containing a subset of the original channel heads which corresponds to the moved heads.
  /// @author SWDG
  /// @date 27/11/15
  void MoveChannelHeadUp(vector<int> Sources, float MoveDist, LSDRaster DEM, vector<int>& UpslopeSources, vector<int>& FinalHeads);

  void HilltopFlowRoutingOriginal(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope, LSDRaster Aspect, LSDIndexRaster StreamNetwork);




  /// @brief Hilltop flow routing.
  ///
  /// @details Hilltop flow routing code built around original code from Martin Hurst. Based on
  /// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
  /// problem of looping flow paths implemented.
  ///
  /// This is the **REFACTORED** version!! It is still slow but hopefully not as slow as before.
  ///
  /// The algorithm now checks for local uphill flows and in the case of identifying one,
  /// D8 flow path is used to push the flow into the centre of the steepest downslope
  /// cell, at which point the trace is restarted. The same technique is used to cope
  /// with self intersections of the flow path. These problems are not solved in the
  /// original paper and I think they are caused at least in part by the high resolution
  /// topogrpahy we are using.
  ///
  /// The code is also now built to take a d infinity flow direction raster instead of an
  /// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
  ///
  /// The Basins input raster is used to code each hilltop into a basin to allow basin
  /// averaging to take place.
  ///
  /// The final 5 parameters are used to set up printing flow paths to files for visualisation,
  /// if this is not needed simply pass in false to the two boolean switches and empty variables for the
  /// others, and the code will run as normal.
  ///
  /// The structure of the returned vector< Array2D<float> > is as follows: \n\n
  /// [0] Hilltop Network coded with stream ID \n
  /// [1] Hillslope Lengths \n
  /// [2] Slope \n
  /// [3] Relief \n
  ///
  /// @param Elevation LSDRaster of elevation values.
  /// @param Slope LSDRaster of slope values.
  /// @param Hilltops LSDRaster of hilltops.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param Aspect LSDRaster of Aspect.
  /// @param Prefix String Prefix for output data filename.
  /// @param Basins LSDIndexRaster of basin outlines.
  /// @param PlanCurvature LSDRaster of planform curvature.
  /// @param print_paths_switch If true paths will be printed.
  /// @param thinning Thinning factor, value used to skip hilltops being printed, use 1 to print every hilltop.
  /// @param trace_path The file path to be used to write the path files to, must end with a slash.
  /// @param basin_filter_switch If this switch is true only basins in Target_Basin_Vector will have their paths printed.
  /// @param Target_Basin_Vector Vector of Basin IDs that the user wants to print traces for.
  /// @return Vector of Array2D<float> containing hillslope metrics.
  /// @author SWDG, SMM
  /// @date 12/2/14, 20/10/2021
  vector< Array2D<float> > HilltopFlowRouting_Refactored(LSDRaster& Elevation, LSDRaster& Hilltop_ID, LSDRaster& Slope,
               LSDRaster& Aspect, LSDRaster& HilltopCurv, LSDRaster& PlanCurvature, LSDIndexRaster& StreamNetwork, LSDIndexRaster& Basins,
               string Prefix, bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
               vector<int> Target_Basin_Vector);




  /// @ brief a refactored HFR routine that is much more memory efficient
  void HilltopFlowRouting_TerrifyingRefactored(LSDRaster& Elevation, LSDRaster& Aspect, LSDRaster& PlanCurvature,
                                                LSDIndexRaster& StreamNetwork,
                                                string ridges_csv_name,string Prefix,
                                                bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch);





  /// @brief Hilltop flow routing.
  ///
  /// @details Hilltop flow routing code built around original code from Martin Hurst. Based on
  /// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
  /// problem of looping flow paths implemented.
  ///
  /// This code is SLOW but robust, a refactored version may appear, but there may not be
  /// enough whisky in Scotland to support that endeavour.
  ///
  /// The algorithm now checks for local uphill flows and in the case of identifying one,
  /// D8 flow path is used to push the flow into the centre of the steepest downslope
  /// cell, at which point the trace is restarted. The same technique is used to cope
  /// with self intersections of the flow path. These problems are not solved in the
  /// original paper and I think they are caused at least in part by the high resolution
  /// topogrpahy we are using.
  ///
  /// The code is also now built to take a d infinity flow direction raster instead of an
  /// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
  ///
  /// The Basins input raster is used to code each hilltop into a basin to allow basin
  /// averaging to take place.
  ///
  /// The final 5 parameters are used to set up printing flow paths to files for visualisation,
  /// if this is not needed simply pass in false to the two boolean switches and empty variables for the
  /// others, and the code will run as normal.
  ///
  /// The structure of the returned vector< Array2D<float> > is as follows: \n\n
  /// [0] Hilltop Network coded with stream ID \n
  /// [1] Hillslope Lengths \n
  /// [2] Slope \n
  /// [3] Relief \n
  ///
  /// @param Elevation LSDRaster of elevation values.
  /// @param Slope LSDRaster of slope values.
  /// @param Hilltops LSDRaster of hilltops.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param Aspect LSDRaster of Aspect.
  /// @param Prefix String Prefix for output data filename.
  /// @param Basins LSDIndexRaster of basin outlines.
  /// @param PlanCurvature LSDRaster of planform curvature.
  /// @param print_paths_switch If true paths will be printed.
  /// @param thinning Thinning factor, value used to skip hilltops being printed, use 1 to print every hilltop.
  /// @param trace_path The file path to be used to write the path files to, must end with a slash.
  /// @param basin_filter_switch If this switch is true only basins in Target_Basin_Vector will have their paths printed.
  /// @param Target_Basin_Vector Vector of Basin IDs that the user wants to print traces for.
  /// @return Vector of Array2D<float> containing hillslope metrics.
  /// @author SWDG
  /// @date 12/2/14
  vector< Array2D<float> > HilltopFlowRouting(LSDRaster Elevation, LSDRaster Hilltop_ID, LSDRaster Slope,
               LSDRaster Aspect, LSDRaster HilltopCurv, LSDRaster PlanCurvature, LSDIndexRaster StreamNetwork, LSDIndexRaster Basins,
               string Prefix, bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
               vector<int> Target_Basin_Vector);

  /// @brief Hilltop flow routing which runs on unsmoothed topography.
  ///
  /// @details Hilltop flow routing code built around original code from Martin Hurst. Based on
  /// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
  /// problem of looping flow paths implemented.
  ///
  /// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
  /// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN
  /// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
  ///
  /// This code is SLOW but robust, a refactored version may appear, but there may not be
  /// enough whisky in Scotland to support that endeavour.
  ///
  /// The algorithm now checks for local uphill flows and in the case of identifying one,
  /// D8 flow path is used to push the flow into the centre of the steepest downslope
  /// cell, at which point the trace is restarted. The same technique is used to cope
  /// with self intersections of the flow path. These problems are not solved in the
  /// original paper and I think they are caused at least in part by the high resolution
  /// topogrpahy we are using.
  ///
  /// The code is also now built to take a d infinity flow direction raster instead of an
  /// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
  ///
  /// The Basins input raster is used to code each hilltop into a basin to allow basin
  /// averaging to take place.
  ///
  /// The final 5 parameters are used to set up printing flow paths to files for visualisation,
  /// if this is not needed simply pass in false to the two boolean switches and empty variables for the
  /// others, and the code will run as normal.
  ///
  /// The structure of the returned vector< Array2D<float> > is as follows: \n\n
  /// [0] Hilltop Network coded with stream ID \n
  /// [1] Hillslope Lengths \n
  /// [2] Slope \n
  /// [3] Relief \n
  ///
  /// @param Elevation LSDRaster of elevation values.
  /// @param Slope LSDRaster of slope values.
  /// @param Hilltops LSDRaster of hilltops.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param D_inf_Flowdir LSDRaster of flow directions.
  /// @param Prefix String Prefix for output data filename.
  /// @param Basins LSDIndexRaster of basin outlines.
  /// @param PlanCurvature LSDRaster of planform curvature.
  /// @param print_paths_switch If true paths will be printed.
  /// @param thinning Thinning factor, value used to skip hilltops being printed, use 1 to print every hilltop.
  /// @param trace_path The file path to be used to write the path files to, must end with a slash.
  /// @param basin_filter_switch If this switch is true only basins in Target_Basin_Vector will have their paths printed.
  /// @param Target_Basin_Vector Vector of Basin IDs that the user wants to print traces for.
  /// @return Vector of Array2D<float> containing hillslope metrics.
  /// @author SWDG
  /// @date 12/2/14
  vector< Array2D<float> > HilltopFlowRouting_RAW(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
               LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
               bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
               vector<int> Target_Basin_Vector);

  /// @brief Hilltop flow routing which generates elevation profiles.
  ///
  /// @details Hilltop flow routing code built around original code from Martin Hurst. Based on
  /// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
  /// problem of looping flow paths implemented.
  ///
  /// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
  /// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN
  /// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
  ///
  /// This code is SLOW but robust, a refactored version may appear, but there may not be
  /// enough whisky in Scotland to support that endeavour.
  ///
  /// The algorithm now checks for local uphill flows and in the case of identifying one,
  /// D8 flow path is used to push the flow into the centre of the steepest downslope
  /// cell, at which point the trace is restarted. The same technique is used to cope
  /// with self intersections of the flow path. These problems are not solved in the
  /// original paper and I think they are caused at least in part by the high resolution
  /// topogrpahy we are using.
  ///
  /// The code is also now built to take a d infinity flow direction raster instead of an
  /// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
  ///
  /// The Basins input raster is used to code each hilltop into a basin to allow basin
  /// averaging to take place.
  ///
  /// The final 5 parameters are used to set up printing flow paths to files for visualisation,
  /// if this is not needed simply pass in false to the two boolean switches and empty variables for the
  /// others, and the code will run as normal.
  ///
  /// The structure of the returned vector< Array2D<float> > is as follows: \n\n
  /// [0] Hilltop Network coded with stream ID \n
  /// [1] Hillslope Lengths \n
  /// [2] Slope \n
  /// [3] Relief \n
  ///
  /// @param Elevation LSDRaster of elevation values.
  /// @param Slope LSDRaster of slope values.
  /// @param Hilltops LSDRaster of hilltops.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param D_inf_Flowdir LSDRaster of flow directions.
  /// @param Prefix String Prefix for output data filename.
  /// @param Basins LSDIndexRaster of basin outlines.
  /// @param print_paths_switch If true paths will be printed.
  /// @param thinning Thinning factor, value used to skip hilltops being printed, use 1 to print every hilltop.
  /// @param trace_path The file path to be used to write the path files to, must end with a slash.
  /// @param basin_filter_switch If this switch is true only basins in Target_Basin_Vector will have their paths printed.
  /// @param Target_Basin_Vector Vector of Basin IDs that the user wants to print traces for.
  /// @return Vector of Array2D<float> containing hillslope metrics.
  /// @author SWDG
  /// @date 25/03/15
  vector< Array2D<float> > HilltopFlowRouting_Profile(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
                                                         LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir, string Prefix, LSDIndexRaster Basins,
                                                         bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                                                         vector<int> Target_Basin_Vector);

  /// @brief This function creates a mask depicting all cells that are influenced
  /// by a pixel that is either on the edge of the DEM or adjacent to a
  /// NoData node
  /// @param Bordered_mask and LSDIndexRaster that is created by the function
  /// find_cells_bordered_by_nodata() in LSDRaster
  /// @param Topography this is the LSDRaster containing topographic data
  /// @return Influenced_mask the LSDIndexRaster that has values 1 where pixels
  /// are influenced by a pixel on the border or next to nodata. They have 0 otherwise
  /// @author SMM
  /// @date 31/10/14
  LSDIndexRaster find_cells_influenced_by_nodata(LSDIndexRaster& Bordered_mask,
                                 LSDRaster& Topography);


  /// @brief This function looks at all uplope nodes and sees if any are bordered
  ///  by nodata.
  /// @param nodeindex The node index of the node in question
  /// @param test_raster and LSDRaster that is to be tested
  /// @return true if the node is influenced by nodata
  /// @author SMM
  /// @date 29/05/2017
  bool is_upstream_influenced_by_nodata(int nodeindex, LSDRaster& test_raster);

  /// @brief This function gets nodes that are possibly on basin edge by
  ///  removing those that do not border NoData. Intended to be passed
  ///  to function for finding concave hull of basin
  /// @param outlet node The node of the outlet
  /// @param Topography this is the LSDRaster containing topographic data
  /// @return A vector with the node indices of nodes that are adjacent to
  ///  nodata within the basin
  /// @author SMM
  /// @date 25/04/2017
  vector<int> basin_edge_extractor(int outlet_node, LSDRaster& Topography);

  /// @brief This function takes an integer raster that has some tags
  ///  and then tags all upstream nodes with that tag
  ///  designed to help extract the main drainage divide of landscapes
  /// @param tagged_raster an LSDIndexRaster containing the pixels from which we search upstream
  /// @return An index raster that has all pixels upslope of the tagged regions taking the tagged value
  /// @author SMM
  /// @date 11/11/2020
  LSDIndexRaster tag_nodes_upstream_of_baselevel(LSDIndexRaster& tagged_raster);

  /// @brief This function takes an integer raster that has some tags
  ///  and then tags all upstream nodes with that tag
  ///  designed to help extract the main drainage divide of landscapes
  /// @param tagged_raster an LSDIndexRaster containing the pixels from which we search upstream
  /// @param crit_downslope_distance The distance beyond which we don't tag the raster.
  /// @return An index raster that has all pixels upslope of the tagged regions taking the tagged value
  /// @author SMM
  /// @date 11/11/2020
  LSDIndexRaster tag_upstream_nodes(LSDIndexRaster& tagged_raster, float crit_downslope_distance);


  /// @brief This function takes an integer raster that has some tags
  ///  and then tags all downstream nodes with that tag. It stops when you have got too far downstream
  /// @param tagged_raster an LSDIndexRaster containing the pixels from which we search upstream
  /// @param crit_downslope_distance The distance beyond which we don't tag the raster.
  /// @return An index raster that has all pixels downstream of the tagged regions taking the tagged value
  ///   within a flow distance window
  /// @author SMM
  /// @date 11/12/2020
  LSDIndexRaster tag_downstream_nodes(LSDIndexRaster& tagged_raster, float crit_downslope_distance);


  /// @brief This function returns all the values from a raster for a corresponding
  /// input vector of node indices.
  /// @param An LSDRaster - must have same dimensions as the LSDFlowInfo object
  /// @param vector<float> - the node indices for which you want the values
  vector<float> get_raster_values_for_nodes(LSDRaster& Raster, vector<int>& node_indices);

  void D_Inf_single_trace_to_channel(LSDRaster Elevation, int start_node, LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir,
                                                          vector< vector<float> >& output_trace_coordinates, vector<float>& output_trace_metrics,
                                                          int& output_channel_node, bool& skip_trace);



  vector< Array2D<float> > HilltopFlowRoutingBedrock(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
               LSDIndexRaster StreamNetwork, LSDRaster Aspect, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
               bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
               vector<int> Target_Basin_Vector, LSDRaster RockExposure);

  /// @brief This method removes end nodes which are not the uppermost extent of the channel network.
  /// @param Ends an LSDIndexRaster of the end points to be processed.
  /// @return A vector of source nodes
  /// @author SWDG
  /// @date 23/7/15
  vector<int> ProcessEndPointsToChannelHeads(LSDIndexRaster Ends);

  /// @brief This method removes single pixel channels from a channel network.
  /// @param StreamNetwork an LSDIndexRaster of the channel network generated from Sources.
  /// @param Sources a vetcor of integer node indices which need cleaned
  /// @return A vector of source nodes
  /// @author SWDG
  /// @date 23/7/15
  vector<int> RemoveSinglePxChannels(LSDIndexRaster StreamNetwork, vector<int> Sources);

  /// @brief This function starts from a source and goes downstream until it
  ///  either accumulates n_nodes_to_visit or hits a base level node
  /// @param source_node The starting node
  /// @param outlet_node A node that serves as an end to the channel before
  ///  the base level. If this is set to a node not on the channel (e.g., -9999)
  ///  then the node looks for accuulation or a base level node only
  /// @param n_nodes_to_visit the number of visited pixels the flow function will
  ///   travese before it stops
  /// @VisitedRaster A raster that has pixels indicating if they have been visited
  ///  or not. A visited pixel is denoted by 1, non visted by 0
  /// @return outlet_nde the node at the end of the flow path
  /// @author SMM
  /// @date 19/05/2016
  int get_downslope_node_after_fixed_visited_nodes(int source_node,
                 int outlet_node, int n_nodes_to_visit, LSDIndexRaster& VisitedRaster);

	/// @brief This function gets the flow distance between two nodes
  /// @param UpstreamNode the upstream node
  /// @param Downstreamnode the downstream node
  /// @return flow length between the nodes
  /// @author FJC
  /// @date 29/09/2016
	float get_flow_length_between_nodes(int UpstreamNode, int DownstreamNode);

  /// @brief This function gets the slope between two nodes based on flow length
  /// @param upslope_node the upstream node
  /// @param downslope_node the downstream node
  /// @param Elevation elevation raster
  /// @return slope between the nodes
  /// @author FJC
  /// @date 04/05/18
  float get_slope_between_nodes(int upslope_node, int downslope_node, LSDRaster& Elevation);

  /// @brief This function gets the Euclidian distance between two nodes
  /// @param node_A the first node
  /// @param node_B the second node
  /// @return distance between the two nodes
  /// @author FJC
  /// @date 17/02/17
  float get_Euclidian_distance(int node_A, int node_B);

  /// @brief Method to snap a point, given as raster coordinates, to a cell in
  /// a raster of hilltops.
  /// @param a Integer row index.
  /// @param b Integer row index.
  /// @param search_radius The radius of the search window in pixels.
  /// @param Hilltops LSDRaster of the hilltops to be snapped to.
  ///@ return The nodeindex of the snapped point.
  /// @author SWDG
  /// @date 23/1/17
  int snap_to_hilltop(int a, int b, int search_radius, LSDRaster& Hilltops);

  /// @brief Wrapper around snap_to_hilltop function to process a collection of utm points.
  ///
  /// @details Writes the the nodeindex of each snapped point to SnappedNodes and the
  /// coordinate count (first coordinate pair is 0, second is 1 and so on) is written
  /// to Valid_node_IDs.
  /// @param x_locs UTM x coordinates.
  /// @param y_locs UTM y coordinates.
  /// @param search_radius The radius of the search window in pixels.
  /// @param Hilltops LSDRaster of the hilltops to be snapped to.
  /// @param SnappedNodes Empty vector where the snapped nodeindexes are written.
  /// @param Valid_node_IDs Empty vector where the coordinate count is written.
  /// @author SWDG
  /// @date 23/1/17
  void snap_to_hilltops(vector<float> x_locs, vector<float> y_locs, int search_radius, LSDRaster& Hilltops, vector<int>& SnappedNodes, vector<int>& Valid_node_IDs);

  /// @brief Extract distinct vector of flows for the entire landscape.
  /// @return Avector of vector of node index
  /// @author BG
  /// @date 29/10/2018
  vector<vector<int> > get_vectors_of_flow(LSDRaster& topo);

  /// Accessor for the DonorStackVector()
  vector<int> get_DonorStackVector() {return DonorStackVector;}
  /// Accessor for the DonorStackVector()
  vector<int> get_RowIndex() {return RowIndex;}
  /// Accessor for the DonorStackVector()
  vector<int> get_ColIndex() {return ColIndex;}


  /// @brief This is a function that produces a **HUGE** data map, but
  ///  speeds up a number of computations (i.e., memory unfriendly, but computation speed friendly)
  ///  It produces a map that has vectors of vecors where each vector is a baselevel node
  ///  and the map strings are the various data elements such as stack order, inverted stack order, rows, columns, etc
  /// @return A map of vecvecs where the key is the name of the data member and the vector into the vectors is indexed by the baselevel node.
  /// @author BG
  /// @date 01/04/2019
  map<string, vector< vector<int> > > get_map_of_vectors();


  /// @brief This return an LSDRaster draining to a specific node, this only fill it with the node ctually in the basin
  /// @brief And fill the rest with nodata. To be used with test_edge deactivated
  /// @param int node, the node index of the outlet
  /// @param LSDRaster the base elevation raster
  /// @return a LSDRater trimmed to the closest extents
  LSDRaster get_raster_draining_to_node(int node, LSDRaster& elevation_raster);

  /// @brief This takes a node list and returns true if the node is surrounded by
  ///  other nodes in the node list.
  /// @param node_list, a vector of nodes. This will almost always be a basin derived from get_upslope_nodes
  /// @return a vector of bools. True if an internal node
  /// @author SMM
  /// @date 11/01/2021
  vector<bool> internal_nodes(vector<int> node_list);

  /// @brief This takes a node and a vector of nodes and finds the nearest node in the vector to the target node.
  /// @author FJC
  /// @date 04/06/21
  int find_nearest_node_in_vector(int target_node, vector<int> node_vector);

  /// @brief Find the elevation difference between the channel and surrounding pixels using the
  /// swath tool.
  /// @return raster of channel relief.
  /// @author FJC
  /// @date 09/03/2021
  // LSDRaster calculate_channel_relief(vector <LSDRaster>& swath_rasters, LSDRaster& topography_raster, vector<int> channel_node_list);


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

  /// The number of nodes in the raster that have data.
  int NDataNodes;

  /// An array that says what node number is at a given row and column.
  Array2D<int> NodeIndex;

  /// @brief A raster of flow direction information.
  ///
  /// In the format:
  ///
  /// 7  0 1 \n
  /// 6 -1 2 \n
  /// 5  4 3 \n
  ///
  /// Nodes with flow direction of -1 drain to themselvs and are base level/sink nodes.
  Array2D<int> FlowDirection;

  /// @brief A code to denote the flow length from the node to its receiver node.
  /// <b>Each node has one and only one receiver.</b>
  /// \n\n
  /// 0 == no receiver/self receiver (base level) \n
  /// 1 == cardinal direction, flow length = DataResolution \n
  /// 2 == diagonal, flow length = DataResolution*(1/sqrt(2)) \n
  Array2D<int> FlowLengthCode;

  /// @brief This stores the row of a node in the vectorized
  /// node index. It, combined with ColIndex, is the
  /// inverse of NodeIndex.
  vector<int> RowIndex;

  /// @brief This stores the column of a node in the vectorized
  /// node index. It, combined with RowIndex, is the
  /// inverse of NodeIndex.
  vector<int> ColIndex;

  /// A list of base level nodes.
  vector<int> BaseLevelNodeList;

  /// Stores the number of donors to each node.
  vector<int> NDonorsVector;

  /// Stores the node index of the receiving node.
  vector<int> ReceiverVector;

  /// @brief Stores the delta vector which is used to index into the donor stack
  /// and order contributing nodes. See Braun and Willett (2012).
  vector<int> DeltaVector;

  /// This is a vector that stores the donor nodes of of the nodes and is
  /// indexed by the DeltaVector.
  vector<int> DonorStackVector;

  /// @brief This vector is used to caluculate flow accumulation. For each base
  /// level node it progresses from a hilltop to a confluence and then jumps to
  /// the next hilltop so that by cascading down through the node indices in
  /// this list one can quickly calculate drainage area, discharge, sediment
  /// flux, etc.
  vector<int> SVector;

  /// This stores the base level node for all of the nodes in the DEM.
  vector<int> BLBasinVector;

  /// This points to the starting point in the S vector of each node.
  vector<int> SVectorIndex;

  /// @brief The number of contributing nodes <b>INCULDING SELF</b> to a current
  /// pixel. It is used in conjunction with the SVectorIndex to build
  /// basins upslope of any and all nodes in the node list.
  vector<int> NContributingNodes;

  /// @brief Boundary conditions stored in a vector of four strings.
  /// The conditions are North[0] East[1] South[2] West[3].
  ///
  /// There are 3 kinds of edge boundaries: no flux, base level and periodic.
  ///
  /// The strings can be any length, as long as the first letter corresponds to the
  /// first letter of the boundary condition. It is not case sensitive.
  vector<string> BoundaryConditions;

  /// @brief map of donorless nodes. Key is node index, val is pair of row,col.
  map<int, pair<int,int> > No_Donor_Map;

  private:
    void create();
    void create(string fname);
    void create(LSDRaster& TopoRaster);
    void create(vector<string>& temp_BoundaryConditions, LSDRaster& TopoRaster);
};

#endif
