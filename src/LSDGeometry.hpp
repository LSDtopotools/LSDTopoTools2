//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDGeometry
// Land Surface Dynamics LSDGeometry objects
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for dealing with geometric data
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
// either version 3 of the License, or (at your option) any later version.
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
// LSDGeometry.cpp
// LSDGeometry object
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



#ifndef LSDGeometry_HPP
#define LSDGeometry_HPP

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
#include "LSDRasterInfo.hpp"
using namespace std;
using namespace TNT;


/// @brief This object packages a number of tools for chi analysis
class LSDGeometry
{
  public:

    /// @brief Empty create function. Leads to some empty vectors
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry()  { create(); }

    /// @brief Create with two vectors: no UTM provided so assumed lat/long
    /// @detail The X should be Longitude
    ///  the Y vector should be latitude
    /// @param x The longitdue data in a vector
    /// @param y the latitude data in a vector
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<double> x, vector<double> y)  { create(x,y); }

    /// @brief Create with two vectors: no UTM provided so assumed lat/long
    /// @detail The X should be Longitude
    ///  the Y vector should be latitude
    /// @param x The longitdue data in a vector
    /// @param y the latitude data in a vector
    /// @author SMM
    /// @date 10/06/2016
    LSDGeometry(vector<float> x, vector<float> y)  { create(x,y); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    ///  This version assumes northern hemisphere
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<double> x, vector<double> y, int UTMZone)  { create(x,y,UTMZone); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    ///  This version assumes northern hemisphere
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<float> x, vector<float> y, int UTMZone)  { create(x,y,UTMZone); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<double> x, vector<double> y, int UTMZone, bool isNorth)  { create(x,y,UTMZone, isNorth); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDGeometry(vector<float> x, vector<float> y, int UTMZone, bool isNorth)  { create(x,y,UTMZone, isNorth); }

    /// @brief This function converts points from Lat/Long to UTM.
    ///  The UTM zone is set as the zone of the first data point
    ///  If there is no data in the Lat/Long data vectors no action is taken
    /// @author SMM
    /// @date 10/06/2016
    void convert_points_to_UTM();

    /// @brief This function converts points from UTM to Lat/Long.
    ///  The UTM zone is set as the zone of the first data point
    ///  If there is no data in the Lat/Long data vectors no action is taken
    /// @author SMM
    /// @date 13/06/2016
    void convert_points_to_LatLong();

    /// @brief This prints the points to a csv file. It will contain both UTM and
    ///  lat-long coordinates. The UTM zone is the zone of the first point,
    ///  The lat long coordinates are in WGS84
    /// @param  path The path to the outfile. Needs the trailing slash
    /// @param file_prefix The prefix of the file **before extension**. That is,
    ///  this function will add the .csv to the end of the filename
    /// @author SMM
    /// @date 13/06/2016
    void print_points_to_csv(string path, string file_prefix);

    /// @brief This gets vectors continaing the row and columns of
    ///  the points from an LSDRasterInfo object
    /// @param RI and LSDRasterInfo object
    /// @param row_vec a vector containing the rows of the points (this is replaced by the function)
    /// @param col_vec a vector containing the cols of the points (this is replaced by the function)
    /// @author SMM
    /// @date 14/06/2106
    void find_row_and_col_of_points(LSDRasterInfo& RI, vector<int>& RowOfNodes, vector<int>& ColOfNodes);

    /// @brief This gets the row and column of a point based on an LSDRasterInfo
    ///  It will return negative and out of bounds indices: used to work
    ///  with functions for determining affected pixels.
    /// @param RI and LSDRasterInfo object
    /// @param point_index an index into the UTM data vectors
    /// @param RowOfNode the row of the point (this is replaced by the function)
    /// @param ColOfNode the col of the point (this is replaced by the function)
    /// @param IsOutOfBounds a boolean that is true if the point is out of the bounds of the raster
    /// @author SMM
    /// @date 15/06/2106
    void find_row_and_col_of_point_inc_out_of_bounds(LSDRasterInfo& RI,
                      int point_index, int& RowOfNode, int& ColOfNode,
                      bool& IsOutOfBounds);

    /// @brief This gets the row and column of a point based on an LSDRasterInfo
    ///  It will return negative and out of bounds indices: used to work
    ///  with finctions for determining affected pixels.
    /// @param RI and LSDRasterInfo object
    /// @param UTM_Easting an easting location
    /// @param UTM_Northing a northing location
    /// @param RowOfNode the row of the point (this is replaced by the function)
    /// @param ColOfNode the col of the point (this is replaced by the function)
    /// @param IsOutOfBounds a boolean that is true if the point is out of the bounds of the raster
    /// @author SMM
    /// @date 15/06/2106
    void find_row_and_col_of_point_inc_out_of_bounds(LSDRasterInfo& RI,
                      double UTM_Easting, double UTM_northing, int& RowOfNode,
                      int& ColOfNode, bool& IsOutOfBounds);

    /// @brief This function checks to see if the data has been converted to UTM
    ///  and if not updates it.
    /// @author SMM
    /// @date 14/06/2016
    void check_and_update_UTM();

    /// @brief This gets the maximum northing value
    /// @return maximum UTM northing
    /// @author SMM
    /// @date 14/06/2016
    double get_max_UTM_Northing();

    /// @brief This gets the minimum northing value
    /// @return minimum UTM northing
    /// @author SMM
    /// @date 14/06/2016
    double get_min_UTM_Northing();

    /// @brief This gets the maximum Easting value
    /// @return maximum UTM Easting
    /// @author SMM
    /// @date 14/06/2016
    double get_max_UTM_Easting();

    /// @brief This gets the minimum Easting value
    /// @return minimum UTM Easting
    /// @author SMM
    /// @date 14/06/2016
    double get_min_UTM_Easting();


    // the getter functions
    int get_UTMZone() { return UTMZone; }
    bool get_isNorth() { return isNorth; }
    vector<double> get_UTMPoints_Easting() { return UTMPoints_Easting; }
    vector<double> get_UTMPoints_Northing() { return UTMPoints_Northing; }
    vector<double> get_WGS84Points_latitude() { return WGS84Points_latitude; }
    vector<double> get_WGS84Points_longitude() { return WGS84Points_longitude; }


  protected:

    int UTMZone;
    bool isNorth;

    vector<double> UTMPoints_Easting;
    vector<double> UTMPoints_Northing;

    vector<double> WGS84Points_latitude;
    vector<double> WGS84Points_longitude;

    void create(vector<double> x, vector<double> y);
    void create(vector<float> x, vector<float> y);
    void create(vector<double> x, vector<double> y, int UTMZone);
    void create(vector<float> x, vector<float> y, int UTMZone);
    void create(vector<double> x, vector<double> y, int UTMZone, bool isNorth);
    void create(vector<float> x, vector<float> y, int UTMZone, bool isNorth);


  private:
    void create();

};


/// @brief This object packages a number of tools for chi analysis
class LSDPolyline: public LSDGeometry
{
  public:

    /// @brief Create with two vectors: no UTM provided so assumed lat/long
    /// @detail The X should be Longitude
    ///  the Y vector should be latitude
    /// @param x The longitdue data in a vector
    /// @param y the latitude data in a vector
    /// @author SMM
    /// @date 10/06/2016
    LSDPolyline(vector<double> x, vector<double> y)  { create(x,y); }

    /// @brief Create with two vectors: no UTM provided so assumed lat/long
    /// @detail The X should be Longitude
    ///  the Y vector should be latitude
    /// @param x The longitdue data in a vector
    /// @param y the latitude data in a vector
    /// @author SMM
    /// @date 10/06/2016
    LSDPolyline(vector<float> x, vector<float> y)  { create(x,y); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    ///  This version assumes northern hemisphere
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDPolyline(vector<double> x, vector<double> y, int UTMZone)  { create(x,y,UTMZone); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    ///  This version assumes northern hemisphere
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDPolyline(vector<float> x, vector<float> y, int UTMZone)  { create(x,y,UTMZone); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDPolyline(vector<double> x, vector<double> y, int UTMZone, bool isNorth)  { create(x,y,UTMZone, isNorth); }

    /// @brief Create with two vectors. UTM info provided
    /// @detail The X should be Easting or Longitude
    ///  the Y vector should be Northing
    /// @param x The Easting data in a vector
    /// @param y the Northing data in a vector
    /// @author SMM
    /// @date 13/06/2016
    LSDPolyline(vector<float> x, vector<float> y, int UTMZone, bool isNorth)  { create(x,y,UTMZone, isNorth); }


    /// @brief This sets the polyline node_order to simply be the order of
    ///  the points. This only makes the simple polyline if the node_order vector is empty
    /// @author SMM
    /// @date 14/06/2016
    void make_simple_polyline();

    /// @brief This sets the polyline node_order to simply be the order of
    ///  the points. This overwrites any existing node_order data
    /// @author SMM
    /// @date 14/06/2016
    void force_simple_polyline();

    /// @brief this  Gets the pixels of all the nodes along the path of a line
    /// @param RI an LSDRasterInfo object
    /// @param affected_rows A vector containing the row numbers of the affected pixels
    /// @param affected_cols: A vector containing the col numbers of the affected pixels
    /// @param start_index The index into the UTM coordinate vector of the starting node
    /// @param end_index the index into the UTM coodinate index of the ending node
    /// @author SMM
    /// @date 14/06/2016
    void get_affected_pixels_in_line_segment(LSDRasterInfo& RI,
                                             vector<int>& affected_rows, vector<int>& affected_cols,
                                             int start_index, int end_index);

    /// @brief this  Gets the pixels of all the nodes along a segment of a line
    ///  The strategy here is just to increment along a line so this has the possiblity
    ///  of missing nodes if the incrementing is too coarse.
    /// @param RI an LSDRasterInfo object
    /// @param affected_rows A vector containing the row numbers of the affected pixels
    /// @param affected_cols: A vector containing the col numbers of the affected pixels
    /// @param start_index The index into the UTM coordinate vector of the starting node
    /// @param end_index the index into the UTM coodinate index of the ending node
    /// @author SMM
    /// @date 15/06/2016
    void get_affected_pixels_in_line_segment_brute_force(LSDRasterInfo& RI,
                                             vector<int>& affected_rows, vector<int>& affected_cols,
                                             int start_index, int end_index);

    /// @brief this  Gets the pixels of all the nodes along the path of a line
    /// @param RI an LSDRasterInfo object
    /// @param affected_rows A vector containing the row numbers of the affected pixels
    /// @param affected_cols: A vector containing the col numbers of the affected pixels
    /// @author SMM
    /// @date 16/06/2016
    void get_affected_pixels_in_line(LSDRasterInfo& RI,vector<int>& affected_rows, vector<int>& affected_cols);

    /// @brief Method to return a vector of node indexes of every cell intersected by the polyline. Calls get_affected_pixels_in_line().
    /// @param RI an LSDRasterInfo object.
    /// @param FlowInfo an LSDFlowInfo object.
    /// @return Vector of integer node indexes.
    /// @author SWDG
    /// @date 22/7/16
    vector<int> get_flowinfo_nodes_of_line(LSDRasterInfo& RI, LSDFlowInfo& FlowInfo);

/*
    /// @detail THis traces to the next pixel ensureing no nodes are missed
    void trace_to_next_pixel(LSDRasterInfo& RI, double StartEasting,double StartNorthing, int start_row, int start_col,
                         int end_row, int end_col, int current_row, int current_col,
                         double slope,
                         double& PixelEasting, double& PixelNorthing);
*/
  protected:

    /// A vector containing the order in which poins are connected
    vector<int> node_order;

  private:

};


#endif
