//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStrahlerLinks
// Land Surface Dynamics StrahlerLinks
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  that keeps track of Strahler ordered stream links and computes
//  various statistics on them
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
// Copyright (C) 2014 Simon M. Mudd 2014
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
// LSDStrahlerLinks.hpp
// LSDStrahlerLinks object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.1.0		26/10/2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDStrahlerLinks_HPP
#define LSDStrahlerLinks_HPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDStrahlerLinks.hpp"
using namespace std;

///@brief This object holds information on the Strahler links in a channel network
/// whereas the LSDJunctionNetwork object stores every junction, this object
/// stores information about the links that connect different strahelr orders: 
/// e.g. every 1st order channel, every 2nd order channel, etc. 
/// @author SMM
/// @date 28/10/2014
class LSDStrahlerLinks
{
  public:
    /// @brief This defines a Strahler links object, is empty
    /// @author SMM
    /// @date 26/10/14
    LSDStrahlerLinks()   { create(); }
   
    /// @brief This creates the LSDStrahlerLinks object
    /// @param JNetwork a LSDJunctionNetwork object
    /// @param FlowInfo LSDFlowInfo object.
    /// @author SMM
    /// @date 26/10/14
    LSDStrahlerLinks(LSDJunctionNetwork& JNetwork, LSDFlowInfo& FlowInfo)
               { create(JNetwork, FlowInfo); }

    /// @brief this function is called during the create process 
    /// it populates the node, row and col vectors with information
    /// about the location of the source and receiver node, row and column
    /// @param JNetwork a LSDJunctionNetwork object
    /// @param FlowInfo LSDFlowInfo object.
    /// @author SMM
    /// @date 28/10/14
    void populate_NodeRowCol_vecvecs(LSDJunctionNetwork& JNetwork, 
                                                   LSDFlowInfo& FlowInfo);

    /// @brief this function calculates the drops for each link
    /// @param FlowInfo LSDFlowInfo object
    /// @param topo_raster LSDRaster object that contains the elevations. 
    /// @author SMM
    /// @date 28/10/14
    void calculate_drops(LSDFlowInfo& FlowInfo, LSDRaster& topo_raster);
    
    /// @brief this function calculates drainage area of each link
    /// @param FlowInfo LSDFlowInfo object
    /// @author SMM
    /// @date 28/10/14
    void calculate_link_area(LSDFlowInfo& FlowInfo);    

    /// @brief this function prints drops. Modified FJC 25/03/16.
    /// @param data_directory a string containing the data dierctory. Should be
    ///  terminated with a slash
    /// @param DEM_name a string that is used to identify the file
    ///  (typically this will be the name of the DEM)
    /// @author SMM
    /// @date 28/10/14
    void print_drops(string data_directory, string DEM_name);

    /// @brief this function calculates calcualtes which basins contain nodes
    /// that receive flow from nodes on edge or adjacent to nodata
    /// and masks these basins. 
    /// @param FI the LSDFlowInfo object
    /// @param InfluenceMask LSDIndexRaster a mask raster that holds the cells
    /// that receive flow from the edge. This is generated using the LSDFlowInfo
    /// member function find_cells_influenced_by_nodata
    /// @return NotIfluencedByEdgeOrNoData an LSDIndexRaster that has values
    ///  0 for cells that receive flow from an edge or nodata cell, and 1 for cells
    ///  that do not receive flow from edge or nodata-adjacent cells 
    /// @author SMM
    /// @date 01/11/2014
    LSDIndexRaster get_no_edge_influence_mask(LSDFlowInfo& FI, 
                               LSDIndexRaster& Influence_Mask);

    /// @brief This is a one stop function that masks out all pixels
    /// that are in basins receiving flow from pixels either on the edge or 
    /// adjacent to nodata
    /// @param FI the LSDFlowInfo object
    /// @param InfluenceMask LSDIndexRaster a mask raster that holds the cells
    /// that receive flow from the edge. This is generated using the LSDFlowInfo
    /// member function find_cells_influenced_by_nodata
    /// @return NotIfluencedByEdgeOrNoData an LSDIndexRaster that has values
    ///  0 for cells that receive flow from an edge or nodata cell, and 1 for cells
    ///  that do not receive flow from edge or nodata-adjacent cells 
    /// @author SMM
    /// @date 01/11/2014
    LSDRaster get_no_edge_influence_raster(LSDFlowInfo& FI, LSDRaster& topography);
	
	  /// @brief Function to print the number of streams of each order
    /// @param data_directory directory to print file to
		/// @param DEM_name string to identify the file (e.g. the name of the DEM)
    /// @author FJC and MAH
    /// @date 17/03/16
    void print_number_of_streams(string data_directory, string DEM_name);
    
    /// @brief Function to calculate the length of each link of each order
    /// @param FlowInfo LSDFlowInfo object 
    /// @author FJC and MAH
    /// @date 24/03/16
    void calculate_lengths(LSDFlowInfo& FlowInfo);
	
	  /// @brief this function prints the lengths. Creates a different file for each stream order.
    /// @param data_directory a string containing the data dierctory. Should be
    ///  terminated with a slash
    /// @param DEM_name a string that is used to identify the file
    ///  (typically this will be the name of the DEM)
    /// @author FJC
    /// @date 25/03/16
    void print_lengths(string data_directory, string DEM_name);


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
    
    /// a vec vec containing the sources of all the Strahler links
    vector< vector<int> > SourceJunctions;
    
    /// a vec vec containing the end junctions of the Strahler links
    vector< vector<int> > ReceiverJunctions;
    
    /// a vec vec containing the node indices of the sources
    vector< vector<int> > SourceNodes;
    
    /// a vec vec containing the rows of the sources
    vector< vector<int> > SourceRows;    

    /// a vec vec containing the cols of the sources
    vector< vector<int> > SourceCols;  
 
    /// a vec vec containing the node indices of the receivers
    /// note: this does not extend to the junction of higher order:
    /// it stops on the last pixel of the channel of this order
    vector< vector<int> > ReceiverNodes;
    
    /// a vec vec containing the rows of the receivers
    /// note: this does not extend to the junction of higher order:
    /// it stops on the last pixel of the channel of this order
    vector< vector<int> > ReceiverRows;    

    /// a vec vec containing the cols of the receivers
    /// note: this does not extend to the junction of higher order:
    /// it stops on the last pixel of the channel of this order
    vector< vector<int> > ReceiverCols;  
    
    /// a vec vec containing drops of every link
    vector< vector<float> > DropData;   
    
    /// a vec vec containing lengths of every link - added FJC 24/03/16
    vector< vector<float> > LengthData;    
 
    
  private:
    void create();
    void create(LSDJunctionNetwork& JN, LSDFlowInfo& FI);
};

#endif