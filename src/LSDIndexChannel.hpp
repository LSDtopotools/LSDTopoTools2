//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexChannel
// Land Surface Dynamics IndexChannel
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for retaining information of channels
//
// Data is mostly pointers to indeces withing the LSDFlowInfo object
//  For data about the actuall channel details such as elevation,
//  drainage area, etc. once needs to use the derivative class
//  LSDChannel
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

/** @file LSDIndexChannel.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 0.0.1
@brief This object contains the node indices as well as the
row and col indices for individual channel segments.
@details These indices could be arranged arbitrailiy according to channel
junctions or simply nodes downstream of a given node and upstream
 of another arbitrary node EndNode.

@date 30/08/2012
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannel_H
#define LSDIndexChannel_H

///@brief This object contains the node indexes as well as the
///row and col indices for individual channel segments.
class LSDIndexChannel
{
  public:
  /// @brief The create function. This is default and throws an error.
  LSDIndexChannel()		{ create(); }
  /// @brief Create LSDIndexChannel object between a start and end node.
  /// @param StartNode Integer starting node.
  /// @param EndNode Integer ending node.
  /// @param FlowInfo LSDFlowInfo object.
  LSDIndexChannel(int StartNode, int EndNode, LSDFlowInfo& FlowInfo)
                 { create(StartNode, EndNode, FlowInfo); }
  /// @brief Create LSDIndexChannel object between starting junction and node and ending junction and node.
  /// @param StartJunction Integer starting junction.
  /// @param StartNode Integer starting node.
  /// @param EndJunction Integer ending junction.
  /// @param EndNode Integer ending node.
  /// @param FlowInfo LSDFlowInfo object.
  LSDIndexChannel(int StartJunction, int StartNode,
                  int EndJunction, int EndNode, LSDFlowInfo& FlowInfo)
              { create(StartJunction,StartNode,EndJunction,EndNode, FlowInfo); }

  LSDIndexChannel(vector<float>& X_coords, vector<float>& Y_coords, LSDFlowInfo& FlowInfo, float threshold_area, float threshold_distance)
            { create(X_coords, Y_coords, FlowInfo, threshold_area, threshold_distance); }

  // get functions

  /// @return Starting junction ID.
  int get_StartJunction() const			{ return StartJunction; }
  /// @return Ending junciton ID.
  int get_EndJunction() const				{ return EndJunction; }
  /// @return Starting node ID.
  int get_StartNode() const				{ return StartNode; }
  /// @return Ending Node ID.
  int get_EndNode() const					{ return EndNode; }

  /// @return Number of rows as an integer.
  int get_NRows() const				{ return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const				{ return NCols; }
  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const			{ return XMinimum; }
  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const			{ return YMinimum; }
	/// @return Data resolution as an integer.
	float get_DataResolution() const	{ return DataResolution; }
	/// @return No Data Value as an integer.
	int get_NoDataValue() const			{ return NoDataValue; }
  /// @return Georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }

  /// @return Get vector of row indexes.
	vector<int> get_RowSequence() const		{ return RowSequence; }
	/// @return Get vector of column indexes.
  vector<int> get_ColSequence() const		{ return ColSequence; }
	/// @return Get vector of node indexes.
  vector<int> get_NodeSequence() const	{ return NodeSequence; }

	/// @return This tells how many nodes are in the channel.
	int get_n_nodes_in_channel() const		{return int(NodeSequence.size()); }

  /// @brief gets the row and column of the end node
  /// @param FlowInfo the LSDFlowInfo object
  /// @param row the row that will be replaced
  /// @param col the col that will be replaced with the correct value
  /// @author SMM
  /// @date 06/05/2015
  void get_row_column_of_end_node(LSDFlowInfo& FlowInfo, int& row, int& col);

  /// @brief This gets the node index at a given node in the index channel.
  /// @param n_node index of target node.
	/// @return Node index.
	/// @author SMM
  /// @date 01/01/12
	int get_node_in_channel(int n_node);

  /// @brief This removes the final node in a channel
  /// @details Used when channel ends in a receiver junction but user wants penultamite node
  /// @author SMM
  /// @date 05/11/2013
  void truncate_final_node();

  /// @brief Get the number of contributing pixels at a given node in the channel.
	/// @param n_node index of target node.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Contributing pixels.
	/// @author SMM
  /// @date 01/01/12
	int get_contributing_pixels_at_node(int n_node, LSDFlowInfo& FlowInfo);

	/// @brief This gets the node, row, and column index.
	/// @param n_node index of target node.
	/// @param node Blank variable for output node.
	/// @param row Blank variable for output row index.
	/// @param col Blank variable for output column index.
  /// @return Node, row, and column index.
	/// @author SMM
  /// @date 01/01/12
	void get_node_row_col_in_channel(int n_node, int& node, int& row, int& col);

  /// @brief This gets the contriubting pixels at the outlet of the channel.
  /// @param FlowInfo LSDFlowInfo object.
	/// @return Contriubting pixels at the outlet of the channel.
	/// @author SMM
  /// @date 01/01/12
	int get_contributing_pixels_at_outlet(LSDFlowInfo& FlowInfo)
          { return FlowInfo.retrieve_contributing_pixels_of_node(EndNode); }
	/// @brief Gets the pixels at the penultimate node.
	///
  /// @details Useful when calculating area of basins where the tributary junction is at the EndNode.
  /// @param FlowInfo LSDFlowInfo object.
  /// @return Pixels at the penultimate node.
	/// @author SMM
  /// @date 01/01/12
	int get_contributing_pixels_at_penultimate_node(LSDFlowInfo& FlowInfo);

	/// @return This prints the channel to an LSDIndexRaster in order to see where the channel is.
	/// @author SMM
  /// @date 01/01/12
	LSDIndexRaster print_index_channel_to_index_raster();

	/// @brief This function appends a channel onto an existing LSDIndexRaster
	/// @author SMM
    /// @date 26/09/2013
	void append_index_channel_to_index_raster(LSDIndexRaster& old_raster);

  /// @brief Function to get vectors with the X and Y coordinates of nodes in the channel
  /// @param X_coordinates vector to write X_coords
  /// @param Y_coordinates vector to write Y_coords
  /// @author FJC
  /// @date 17/02/17
  void get_coordinates_of_channel_nodes(vector<double>& X_coordinates, vector<double>& Y_coordinates, LSDFlowInfo& FlowInfo);

  /// @brief Function to get vectors with the X and Y coordinates of nodes in the channel
  /// @param X_coordinates vector to write X_coords
  /// @param Y_coordinates vector to write Y_coords
  /// @author FJC
  /// @date 17/02/17
  void get_coordinates_of_channel_nodes(vector<float>& X_coordinates, vector<float>& Y_coordinates);


  /// @brief Function to get write index channel to csv
  /// @param path the path name
  /// @param filename prefix for the csv file
  /// @author FJC
  /// @date 17/02/17
  void write_channel_to_csv(string path, string filename);

  /// @brief Function to get write index channel to csv
  /// @detail Similar to above but mimics code in other objects to ensure you get the same channels
  /// @param path the path name
  /// @param filename prefix for the csv file
  /// @param FlowInfo
  /// @param FlowDistance
  /// @author SMM
  /// @date 20/03/18
  void write_channel_to_csv(string path, string filename, LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& Elevation);

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

  /// The starting junction (numbered within LSDJunctionNetwork object).
	int StartJunction;
	/// The node index of the starting Junction (numbered in the FlowInfo object).
  int StartNode;
	/// The ending junction (numbered within LSDJunctionNetwork object).
  int EndJunction;
	/// The node index of the ending Junction (numbered in the FlowInfo object).
  int EndNode;

  ///Vector of row indices.
	vector<int> RowSequence;
	///Vector of column indices.
  vector<int> ColSequence;
	///Vector of node indices.
  vector<int> NodeSequence;

	private:
	void create();
	void create(int StartNode, int EndNode, LSDFlowInfo& FlowInfo);
	void create(int StartJunction, int StartNode,
	            int EndJunction, int EndNode, LSDFlowInfo& FlowInfo);
	void create(vector<float>& X_coords, vector<float>& Y_coords, LSDFlowInfo& FlowInfo, float threshold_area, float threshold_distance);


};



#endif
