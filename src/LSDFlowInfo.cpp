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

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo.cpp
// cpp file for the LSDFlowInfo object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.1.0    21/10/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDFlowInfo_CPP
#define LSDFlowInfo_CPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <limits>
#include <cstring>
#include <queue>
#include <algorithm>
#include <math.h>
#include <unordered_set>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDStatsTools.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDSpatialCSVReader.hpp"
//#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create function, this is empty, you need to include a filename
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::create()
{
  //cout << endl << "-------------------------------" << endl;
  //cout << "I am an empty flow info object. " << endl;
  //cout << "Did you forget to give me a DEM?" << endl;
  //cout << endl << "-------------------------------" << endl;
  //exit(EXIT_FAILURE);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Create function, this creates from a pickled file
// fname is the name of the pickled flow info file
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::create(string fname)
{
  unpickle(fname);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This defaults to no flux boundary conditions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::create(LSDRaster& TopoRaster)
{
  vector<string> BoundaryConditions(4, "No Flux");
  create(BoundaryConditions, TopoRaster);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function calcualtes the receiver nodes
// it returns the receiver vector r_i
// it also returns a flow direction array in this ordering:
//
// 7  0 1
// 6 -1 2
// 5  4 3
//
// note this is different from ArcMap flowdirection
//  int Arc_flowdir;      // flow direction in arcmap format
//                // 32  64  128
//                // 16  --  1
//                // 8    4  2
// one can convert nthese indices using the LSDIndexRaster object
// note in arc the row index increases down (to the south)
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::create(vector<string>& temp_BoundaryConditions,
                         LSDRaster& TopoRaster)
{

  // initialize several data members
  BoundaryConditions = temp_BoundaryConditions;
  //cout << "TBC" << endl;
  NRows = TopoRaster.get_NRows();
  //cout << "Rows: " << NRows << endl;
  NCols = TopoRaster.get_NCols();
  //cout << "Cols: " << NCols << endl;
  XMinimum = TopoRaster.get_XMinimum();
  //cout << "Xmin: " << XMinimum << endl;
  YMinimum = TopoRaster.get_YMinimum();
  //cout << "Ymin: " << YMinimum << endl;
  NoDataValue = int(TopoRaster.get_NoDataValue());
  //cout << "NDV: " << NoDataValue << endl;
  DataResolution = TopoRaster.get_DataResolution();
  //cout << "Data resolution: " <<DataResolution << endl;

  GeoReferencingStrings = TopoRaster.get_GeoReferencingStrings();
  //cout << "GRS" << endl;

  //cout << "1" << endl;

  // Declare matrices for calculating flow routing
  float one_ov_root2 = 0.707106781;
  float target_elev;        // a placeholder for the elevation of the potential receiver
  float slope;
  float max_slope;        // the maximum slope away from a node
  int max_slope_index;      // index into the maximum slope

  int row, col;            // index for the rows and column
  int receive_row,receive_col;
  string::iterator string_iterator;  // used to get characters from string

  // we need logic for all of the boundaries.
  // there are 3 kinds of edge boundaries:
  // no flux
  // base level
  // periodic
  // These are denoted in a vector of strings.
  // the vector has four elements
  // North boundary, East boundary, South bondary and West boundary
  // the strings can be any length, as long as the first letter corresponds to the
  // first letter of the boundary condition. It is not case sensitive.

  // go through the boundaries
  // A NOTE ON CARDINAL DIRECTIONS
  // If one looks at the raster data, the top row in the data corresponds to the NORTH boundary
  // This is the row that is first read into the code
  // so row 0 is the NORTH boundary
  // row NRows-1 is the SOUTH boundary
  // column 0 is the WEST boundary
  // column NCols-1 is the EAST boundary
  vector<float> slopes(8,NoDataValue);
  vector<int> row_kernal(8);
  vector<int> col_kernal(8);
  int ndv = NoDataValue;
  NDataNodes = 0;       // the number of nodes in the raster that have data
  int one_if_a_baselevel_node;  // this is a switch used to tag baseleve nodes

  // the first thing you need to do is construct a topoglogy matrix
  // the donor, receiver, etc lists are as long as the number of nodes.
  // these are made of vectors that are exactly dimension n, which is the number of nodes with
  // data. Each of these nodes has several index vectors, that point the program to where the node is
  // we construct these index vectors first
  // we need to loop through all the data before we calcualte slopes because the
  // receiver node indices must be known before the slope calculations are run
  vector<int> empty_vec;
  RowIndex = empty_vec;
  ColIndex = empty_vec;
  BaseLevelNodeList = empty_vec;
  ReceiverVector = empty_vec;
  Array2D<int> ndv_raster(NRows,NCols,ndv);

  NodeIndex = ndv_raster.copy();
  FlowDirection = ndv_raster.copy();
  FlowLengthCode = ndv_raster.copy();


  //cout << "2" << endl;

  // loop through the topo data finding places where there is actually data
  for (row = 0; row<NRows; row++)
  {
    for (col = 0; col<NCols; col++)
    {
      // only do calcualtions if there is data
      if(TopoRaster.RasterData[row][col] != NoDataValue)
      {
        RowIndex.push_back(row);
        ColIndex.push_back(col);
        NodeIndex[row][col] = NDataNodes;
        NDataNodes++;
      }
    }
  }

  //cout << "3" << endl;

  // now the row and col index are populated by the row and col of the node in row i
  // and the node index has the indeces into the row and col vectors
  // next up, make d, delta, and D vectors
  vector<int> ndn_vec(NDataNodes,0);
  vector<int> ndn_nodata_vec(NDataNodes,ndv);
  vector<int> ndn_plusone_vec(NDataNodes+1,0);
  vector<int> w_vector(NDataNodes,0);

  NDonorsVector = ndn_vec;
  DonorStackVector = ndn_vec;
  DeltaVector = ndn_plusone_vec;

  SVector = ndn_nodata_vec;
  BLBasinVector = ndn_nodata_vec;

  // this vector starts out empty and then base level nodes are added to it
  for (row = 0; row<NRows; row++)
  {
    for (col = 0; col<NCols; col++)
    {
      // only do calcualtions if there is data
      if(TopoRaster.RasterData[row][col] != NoDataValue)
      {
        // calcualte 8 slopes
        // no slopes mean get NoDataValue entries
        // the algorithm loops through the neighbors to the cells, collecting
        // receiver indices. The order is
        // 7 0 1
        // 6 - 2
        // 5 4 3
        // where the above directions are cardinal directions
        // do slope 0
        row_kernal[0] = row-1;
        row_kernal[1] = row-1;
        row_kernal[2] = row;
        row_kernal[3] = row+1;
        row_kernal[4] = row+1;
        row_kernal[5] = row+1;
        row_kernal[6] = row;
        row_kernal[7] = row-1;

        col_kernal[0] = col;
        col_kernal[1] = col+1;
        col_kernal[2] = col+1;
        col_kernal[3] = col+1;
        col_kernal[4] = col;
        col_kernal[5] = col-1;
        col_kernal[6] = col-1;
        col_kernal[7] = col-1;

        // check for periodic boundary conditions
        if( BoundaryConditions[0].find("P") == 0 || BoundaryConditions[0].find("p") == 0 )
        {
          if( BoundaryConditions[2].find("P") != 0 && BoundaryConditions[2].find("p") != 0 )
          {
            cout << "WARNING!!! North boundary is periodic! Changing South boundary to periodic" << endl;
            BoundaryConditions[2] = "P";
          }
        }
        if( BoundaryConditions[1].find("P") == 0 || BoundaryConditions[1].find("p") == 0 )
        {
          if( BoundaryConditions[3].find("P") != 0 && BoundaryConditions[3].find("p") != 0 )
          {
            cout << "WARNING!!! East boundary is periodic! Changing West boundary to periodic" << endl;
            BoundaryConditions[3] = "P";
          }
        }
        if( BoundaryConditions[2].find("P") == 0 || BoundaryConditions[2].find("p") == 0 )
        {
          if( BoundaryConditions[0].find("P") != 0 && BoundaryConditions[0].find("p") != 0 )
          {
            cout << "WARNING!!! South boundary is periodic! Changing North boundary to periodic" << endl;
            BoundaryConditions[0] = "P";
          }
        }
        if( BoundaryConditions[3].find("P") == 0 || BoundaryConditions[3].find("p") == 0 )
        {
          if( BoundaryConditions[1].find("P") != 0 && BoundaryConditions[1].find("p") != 0 )
          {
            cout << "WARNING!!! West boundary is periodic! Changing East boundary to periodic" << endl;
            BoundaryConditions[1] = "P";
          }
        }

        // reset baselevel switch for boundaries
        one_if_a_baselevel_node = 0;

        // NORTH BOUNDARY
        if (row == 0)
        {
          if( BoundaryConditions[0].find("B") == 0 || BoundaryConditions[0].find("b") == 0 )
          {
            one_if_a_baselevel_node = 1;
          }
          else
          {
            // if periodic, reflect across to south boundary
            if( BoundaryConditions[0].find("P") == 0 || BoundaryConditions[0].find("p") == 0 )
            {
              row_kernal[0] = NRows-1;
              row_kernal[1] = NRows-1;
              row_kernal[7] = NRows-1;
            }
            else
            {
              row_kernal[0] = ndv;
              row_kernal[1] = ndv;
              row_kernal[7] = ndv;
            }
          }
        }
        // EAST BOUNDAY
        if (col == NCols-1)
        {
          if( BoundaryConditions[1].find("B") == 0 || BoundaryConditions[1].find("b") == 0 )
          {
            one_if_a_baselevel_node = 1;
          }
          else
          {
            if( BoundaryConditions[1].find("P") == 0 || BoundaryConditions[1].find("p") == 0)
            {
              col_kernal[1] = 0;
              col_kernal[2] = 0;
              col_kernal[3] = 0;
            }
            else
            {
              col_kernal[1] = ndv;
              col_kernal[2] = ndv;
              col_kernal[3] = ndv;
            }
          }
        }
        // SOUTH BOUNDARY
        if (row == NRows-1)
        {
          if( BoundaryConditions[2].find("B") == 0 || BoundaryConditions[2].find("b") == 0 )
          {
            one_if_a_baselevel_node = 1;
          }
          else
          {
            if( BoundaryConditions[2].find("P") == 0 || BoundaryConditions[2].find("p") == 0)
            {
              row_kernal[3] = 0;
              row_kernal[4] = 0;
              row_kernal[5] = 0;
            }
            else
            {
              row_kernal[3] = ndv;
              row_kernal[4] = ndv;
              row_kernal[5] = ndv;
            }
          }
        }
        // WEST BOUNDARY
        if (col == 0)
        {
          if( BoundaryConditions[3].find("B") == 0 || BoundaryConditions[3].find("b") == 0 )
          {
            one_if_a_baselevel_node = 1;
          }
          else
          {
            if( BoundaryConditions[3].find("P") == 0 || BoundaryConditions[3].find("p") == 0)
            {
              col_kernal[5] = NCols-1;
              col_kernal[6] = NCols-1;
              col_kernal[7] = NCols-1;
            }
            else
            {
              col_kernal[5] = ndv;
              col_kernal[6] = ndv;
              col_kernal[7] = ndv;
            }
          }
        }

        // now loop through the surrounding nodes, calculating the slopes
        // slopes with NoData get NoData slopes
        // reminder of ordering:
        // 7 0 1
        // 6 - 2
        // 5 4 3
        // first logic for baselevel node
        if (one_if_a_baselevel_node == 1)
        {
          // get receiver index
          FlowDirection[row][col] = -1;
          ReceiverVector.push_back(NodeIndex[row][col]);
          FlowLengthCode[row][col] = 0;
        }
        // now the rest of the nodes
        else
        {
          FlowLengthCode[row][col] = 0;    // set flow length code to 0, this gets reset
          // if there is a maximum slope
          max_slope = 0;
          max_slope_index = -1;
          receive_row = row;
          receive_col = col;
          for (int slope_iter = 0; slope_iter<8; slope_iter++)
          {
            if (row_kernal[slope_iter] == ndv || col_kernal[slope_iter] == ndv)
            {
              slopes[slope_iter] = NoDataValue;
            }
            else
            {
              target_elev = TopoRaster.RasterData[ row_kernal[slope_iter] ][ col_kernal[slope_iter] ];
              if(target_elev == NoDataValue)
              {
                slopes[slope_iter] = NoDataValue;
              }
              else
              {
                if(slope_iter%2 == 0)
                {
                  //cout << "LINE 988, cardinal direction, slope iter = " << slope_iter << endl;
                  slope = TopoRaster.RasterData[row][col]-target_elev;
                }
                else
                {
                  slope = one_ov_root2*(TopoRaster.RasterData[row][col]-target_elev);
                }

                if (slope > max_slope)
                {
                  max_slope_index = slope_iter;
                  receive_row = row_kernal[slope_iter];
                  receive_col = col_kernal[slope_iter];
                  max_slope = slope;
                  if(slope_iter%2 == 0)
                  {
                    FlowLengthCode[row][col] = 1;
                  }
                  else
                  {
                    FlowLengthCode[row][col] = 2;
                  }
                }
              }
            }
          }
          // get receiver index
          FlowDirection[row][col] = max_slope_index;
          ReceiverVector.push_back(NodeIndex[receive_row][receive_col]);
        }    // end if baselevel boundary  conditional

        // if the node is a base level node, add it to the base level node list
        if (FlowLengthCode[row][col] == 0)
        {
          BaseLevelNodeList.push_back(NodeIndex[row][col]);
        }
      }      // end if there is data conditional
    }        // end col loop
  }          // end row loop



  // first create the number of donors vector
  // from braun and willett eq. 5
  for(int i = 0; i<NDataNodes; i++)
  {
    NDonorsVector[ ReceiverVector[i] ]++;
  }


  // now create the delta vector
  // this starts on the last element and works its way backwards
  // from Braun and Willett eq 7 and 8
  DeltaVector[NDataNodes] = NDataNodes;
  for(int i = NDataNodes; i>0; i--)
  {
    DeltaVector[i-1] = DeltaVector[i] -  NDonorsVector[i-1];
  }


  // now the DonorStack and the r vectors. These come from Braun and Willett
  // equation 9.
  // Note that in the manscript I have there is a typo in eqaution 9
  // (Jean Braun's code is correct)
  // it should be w_{r_i} = w_{r_i}+1
  int r_index;
  int w_index;
  int delta_index;
  for (int i = 0; i<NDataNodes; i++)
  {
    r_index = ReceiverVector[i];
    delta_index = DeltaVector[ r_index ];
    w_index = w_vector[ r_index ];
    DonorStackVector[  delta_index+w_index ] = i;
    w_vector[r_index] += 1;
    //cout << "i: " << i << " r_i: " << r_index << " delta_i: " << delta_index << " w_index: " << w_index << endl;
  }


  // now go through the base level node list, building the drainage tree for each of these nodes as one goes along
  int n_base_level_nodes;
  n_base_level_nodes = BaseLevelNodeList.size();

  int k;
  int j_index;
  int begin_delta_index, end_delta_index;
  int l_index;

  j_index = 0;
  for (int i = 0; i<n_base_level_nodes; i++)
  {
    k = BaseLevelNodeList[i];      // set k to the base level node

    // This doesn't seem to be in Braun and Willet but to get the ordering correct you
    // need to make sure that the base level node appears first in the donorstack
    // of nodes contributing to the baselevel node.
    // For example, if base level node is 4, with 4 donors
    // and the donor stack has 3 4 8 9
    // the code has to put the 4 first.
    if (DonorStackVector[ DeltaVector[k] ] != k)
    {
      int this_index = DonorStackVector[ DeltaVector[k] ];
      int bs_node = k;

      for(int ds_node = 1; ds_node < NDonorsVector[k]; ds_node++)
      {
        if( DonorStackVector[ DeltaVector[k] + ds_node ] == bs_node )
        {
          DonorStackVector[ DeltaVector[k] ] = k;
          DonorStackVector[ DeltaVector[k] + ds_node ] = this_index;
        }
      }
    }

    // now run recursive algorithm
    begin_delta_index = DeltaVector[k];
    end_delta_index = DeltaVector[k+1];

    //cout << "base_level_node is: " << k << " begin_index: " << begin_delta_index << " end: " << end_delta_index << endl;

    for (int delta_index = begin_delta_index; delta_index<end_delta_index; delta_index++)
    {
      l_index = DonorStackVector[delta_index];
      add_to_stack(l_index, j_index, k);
    }
  }

  // now calcualte the indices
  calculate_upslope_reference_indices();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDFlowInfo::get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc)
{

  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;

  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
// Same as above but with floats
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDFlowInfo::get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc)
{

  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;

  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to convert a node position with a row and column to a lat
// and long coordinate
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDFlowInfo::get_lat_and_long_locations(int row, int col, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter)
{
  // get the x and y locations of the node
  double x_loc,y_loc;
  get_x_and_y_locations(row, col, x_loc, y_loc);

  // get the UTM zone of the node
  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);
  //cout << endl << endl << "Line 1034, UTM zone is: " << UTM_zone << endl;


  if(UTM_zone == NoDataValue)
  {
    lat = NoDataValue;
    longitude = NoDataValue;
  }
  else
  {
    // set the default ellipsoid to WGS84
    int eId = 22;

    double xld = double(x_loc);
    double yld = double(y_loc);

    // use the converter to convert to lat and long
    double Lat,Long;
    Converter.UTMtoLL(eId, yld, xld, UTM_zone, is_North, Lat, Long);


    lat = Lat;
    longitude = Long;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to convert X and Y data to a lat and long coordinate
// For use with hillslope traces where X and Y may not coincide with nodes
// MDH 27/7/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDFlowInfo::get_lat_and_long_locations(double X, double Y, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter)
{
  // get the UTM zone of the node
  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);
  //cout << endl << endl << "Line 714 in LSDFlowInfo, UTM zone is: " << UTM_zone;
  //if(is_North)
  //{
  //  cout << " North" << endl;
  //}
  //else
  //{
  //  cout << " South" << endl;
  //}

  if(UTM_zone == NoDataValue)
  {
    lat = NoDataValue;
    longitude = NoDataValue;
  }
  else
  {
    // set the default ellipsoid to WGS84
    int eId = 22;

    // use the converter to convert to lat and long
    double Lat,Long;
    Converter.UTMtoLL(eId, Y, X, UTM_zone, is_North, Lat, Long);

    lat = Lat;
    longitude = Long;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



void LSDFlowInfo::get_x_and_y_from_latlong(double lat, double longitude, double& X, double& Y )
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);


  // loop through the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  Converter.LLtoUTM_ForceZone(eId, lat, longitude,X, Y, UTM_zone);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDFlowInfo::get_UTM_information(int& UTM_zone, bool& is_North)
{

  // set up strings and iterators
  map<string,string>::iterator iter;

  //check to see if there is already a map info string
  string mi_key = "ENVI_map_info";
  iter = GeoReferencingStrings.find(mi_key);
  if (iter != GeoReferencingStrings.end() )
  {
    string info_str = GeoReferencingStrings[mi_key] ;

    // now parse the string
    vector<string> mapinfo_strings;
    istringstream iss(info_str);
    while( iss.good() )
    {
      string substr;
      getline( iss, substr, ',' );
      mapinfo_strings.push_back( substr );
    }
    UTM_zone = atoi(mapinfo_strings[7].c_str());
    //cout << "Line 1041, UTM zone: " << UTM_zone << endl;
    //cout << "LINE 1042 LSDRaster, N or S: " << mapinfo_strings[7] << endl;

    // find if the zone is in the north
    string n_str = "n";
    string N_str = "N";
    is_North = false;
    size_t found = mapinfo_strings[8].find(N_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    //cout << "is_North is: " << is_North << endl;

  }
  else
  {
    UTM_zone = NoDataValue;
    is_North = false;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Checks to see is a point is in the raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDFlowInfo::check_if_point_is_in_raster(float X_coordinate, float Y_coordinate)
{
  bool is_in_raster = true;

  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum - DataResolution*0.5;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum - DataResolution*0.5;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1)
    {
      is_in_raster = false;
    }

  return is_in_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the receiver of current_node (its node, row, and column)
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_receiver_information(int current_node,
                     int& receiver_node, int& receiver_row,
                     int& receiver_col)
{
  if (current_node >= NDataNodes)
  {
    cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
    exit(0);
  }
  
  int rn, rr, rc;
  rn = ReceiverVector[current_node];
  rr = RowIndex[rn];
  rc = ColIndex[rn];
  receiver_node = rn;
  receiver_row = rr;
  receiver_col = rc;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the receiver of current_node, just its node version
//
// BG 05/01/2018
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_receiver_information(int current_node,
                     int& receiver_node)
{
  if (current_node >= NDataNodes)
  {
    cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
    exit(0);
  }

  int rn;
  rn = ReceiverVector[current_node];
  receiver_node = rn;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the row and column of the current node
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_current_row_and_col(int current_node,int& curr_row,
                                               int& curr_col)
{
  if (current_node >= NDataNodes)
  {
    cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
    exit(0);
  }

  int cr, cc;
  cr = RowIndex[current_node];
  cc = ColIndex[current_node];
  curr_row = cr;
  curr_col = cc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

int LSDFlowInfo::get_NodeIndex_from_row_col(int row, int col)
{
  return NodeIndex[row][col];
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the rows and columns of a vector of node
//
// SMM 15/06/2022
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_rows_and_cols_from_node_list(vector<int> current_nodes,vector<int>& curr_rows,
                                               vector<int>& curr_cols)
{
  
  int n_nodes = int(current_nodes.size());
  int cr, cc;
  vector<int> crs;
  vector<int> ccs;
  int current_node;
  for(int i = 0; i< n_nodes; i++)
  {
    current_node = current_nodes[i];
    if (current_node >= NDataNodes)
    {
      cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
      exit(0);
    }

    crs.push_back(RowIndex[current_node]);
    ccs.push_back(ColIndex[current_node]);
  }
  
  curr_rows = crs;
  curr_cols = ccs;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



vector<int> LSDFlowInfo::retrieve_contributing_pixels_from_node_list(vector<int> current_nodes)
{
  vector<int> NContrib;
  int current_node;
  int n_nodes = int(current_nodes.size());
  for (int i = 0; i<n_nodes; i++)
  {
    current_node = current_nodes[i];
    if (current_node >= NDataNodes)
    {
      cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
      exit(0);
    } 
    NContrib.push_back( NContributingNodes[current_node] );   
  }  
  return NContrib;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the X and Y coordinates of the current node
//
// BG 20/02/2017
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::get_x_and_y_from_current_node(int current_node, float& current_X, float& current_Y)
{
  if (current_node >= NDataNodes)
  {
    cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
    exit(0);
  }

  int cr,cc;
  retrieve_current_row_and_col(current_node, cr,cc);
  get_x_and_y_locations(cr, cc, current_X, current_Y);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// algorithms for searching the vectors
// This gets the lat and long coordinates of the current node
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::get_lat_and_long_from_current_node(int current_node, double& current_lat, double& current_long, LSDCoordinateConverterLLandUTM Converter)
{
  if (current_node >= NDataNodes)
  {
    cout << "Fatal error, you are looking for a flowinfo node that doesn't exist" << endl;
    exit(0);
  }

  int cr,cc;
  retrieve_current_row_and_col(current_node, cr,cc);
  double latitude;
  double longitude;
  get_lat_and_long_locations(cr, cc, latitude, longitude, Converter);

  current_lat = latitude;
  current_long = longitude;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// algorithms for searching the vectors
// This gets the row and column of the current node
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file(vector<int>& nodeindex_vec, string outfilename)
{

  // fid the last '.' in the filename to use in the scv filename
  unsigned dot = outfilename.find_last_of(".");

  string prefix = outfilename.substr(0,dot);
  //string suffix = str.substr(dot);
  string insert = "_nodeindices_for_Arc.csv";
  string outfname = prefix+insert;

  cout << "the Arc filename is: " << outfname << endl;

  int n_nodes = nodeindex_vec.size();
  int n_nodeindeces = RowIndex.size();

  // open the outfile
  ofstream csv_out;
  csv_out.open(outfname.c_str());
  csv_out.precision(8);

  csv_out << "x,y,node,row,col" << endl;

  int current_row, current_col;
  float x,y;

  // loop through node indices in vector
  for (int i = 0; i<n_nodes; i++)
    {
      int current_node = nodeindex_vec[i];

      // make sure the nodeindex isn't out of bounds
      if (current_node < n_nodeindeces)
  {
    // get the row and column
    retrieve_current_row_and_col(current_node,current_row,
                                 current_col);

    // get the x and y location of the node
    // the last 0.0001*DataResolution is to make sure there are no integer data points
    x = XMinimum + float(current_col)*DataResolution + 0.5*DataResolution + 0.0001*DataResolution;

    // the last 0.0001*DataResolution is to make sure there are no integer data points
    // y coord a bit different since the DEM starts from the top corner
    y = YMinimum + float(NRows-current_row)*DataResolution - 0.5*DataResolution + 0.0001*DataResolution;;
    csv_out << x << "," << y << "," << current_node << "," << current_row << "," << current_col << endl;
  }
    }

  csv_out.close();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function takes a list of junctions and prints them to a csv file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int>& nodeindex_vec, string outfilename)
{
  int n_nodes = (nodeindex_vec.size());
  int this_node;
  int row,col;
  double x_loc,y_loc;
  double latitude,longitude;

  // open the outfile
  ofstream sources_out;
  sources_out.open(outfilename.c_str());
  sources_out.precision(9);

  sources_out << "id,x,y,latitude,longitude" << endl;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  for (int i = 0; i<n_nodes; i++)
  {
    this_node = nodeindex_vec[i];

    // get the row and column
    retrieve_current_row_and_col(this_node,row,col);

    // get the x and y locations
    get_x_and_y_locations(row, col, x_loc, y_loc);

    // get the lat and long locations
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    // print to file
    sources_out << this_node << "," << x_loc << ","
                << y_loc << "," << latitude << "," << longitude << endl;

  }

  sources_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function takes a list of junctions and prints them to a csv file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int>& nodeindex_vec, string outfilename, LSDRaster& add_value, string add_column_name)
{
  int n_nodes = (nodeindex_vec.size());
  int this_node;
  int row,col;
  double x_loc,y_loc;
  double latitude,longitude;

  // open the outfile
  ofstream sources_out;
  sources_out.open(outfilename.c_str());
  sources_out.precision(9);

  sources_out << "id,x,y,latitude,longitude," << add_column_name << endl;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  for (int i = 0; i<n_nodes; i++)
  {
    this_node = nodeindex_vec[i];

    // get the row and column
    retrieve_current_row_and_col(this_node,row,col);

    // get the x and y locations
    get_x_and_y_locations(row, col, x_loc, y_loc);

    // get the lat and long locations
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    // print to file
    sources_out << this_node << "," << x_loc << ","
                << y_loc << "," << latitude << "," << longitude << "," << add_value.get_data_element(row,col) << endl;

  }

  sources_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Write nodeindex vector to csv file, and give each row a unique ID
//
// SWDG after SMM 2/2/2016
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_Unique(vector<int>& nodeindex_vec, string outfilename)
{

  // fid the last '.' in the filename to use in the scv filename
  unsigned dot = outfilename.find_last_of(".");

  string prefix = outfilename.substr(0,dot);
  //string suffix = str.substr(dot);
  string insert = "_nodeindices_for_Arc.csv";
  string outfname = prefix+insert;

  cout << "the Arc filename is: " << outfname << endl;

  int n_nodes = nodeindex_vec.size();
  int n_nodeindeces = RowIndex.size();

  // open the outfile
  ofstream csv_out;
  csv_out.open(outfname.c_str());
  csv_out.precision(8);

  csv_out << "x,y,node,row,col,unique_ID" << endl;

  int current_row, current_col;
  float x,y;

  // loop through node indices in vector
  for (int i = 0; i<n_nodes; i++)
  {
    int current_node = nodeindex_vec[i];

    // make sure the nodeindex isn't out of bounds
    if (current_node < n_nodeindeces)
    {
      // get the row and column
      retrieve_current_row_and_col(current_node,current_row,
                                 current_col);

      // get the x and y location of the node
      // the last 0.0001*DataResolution is to make sure there are no integer data points
      x = XMinimum + float(current_col)*DataResolution + 0.5*DataResolution + 0.0001*DataResolution;

      // the last 0.0001*DataResolution is to make sure there are no integer data points
      // y coord a bit different since the DEM starts from the top corner
      y = YMinimum + float(NRows-current_row)*DataResolution - 0.5*DataResolution + 0.0001*DataResolution;;
      csv_out << x << "," << y << "," << current_node << "," << current_row << "," << current_col << "," << i << endl;
    }
  }

  csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function returns the base level node with the greatest drainage area
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::retrieve_largest_base_level()
{
  int n_bl = BaseLevelNodeList.size();    // get the number of baselevel nodes
  int max_bl = 0;
  for (int i = 0; i<n_bl; i++)
    {
      if(NContributingNodes[ BaseLevelNodeList[i] ] > max_bl)
  {
    max_bl = NContributingNodes[ BaseLevelNodeList[i] ];
  }
    }
  return max_bl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function returns the base level node with the greatest drainage area
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::retrieve_base_level_node(int node)
{

  int CurrentNode = node;
  int ReceiverNode, ReceiverRow, ReceiverCol;

  // set initial ReceiverNode so that you can enter while loop
  retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);

  // now follow the nodes down to the receiver.
  while(ReceiverNode != CurrentNode)
  {
    CurrentNode = ReceiverNode;
    retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
  }

  return ReceiverNode;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function removes all pixels that drain from edge nodes
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::remove_nodes_influneced_by_edge(LSDRaster& topography)
{
  int n_nodes = int(SVector.size());

  int row,col;
  int move_row,move_col;

  vector<int> edge_node_list;

  // get a list of nodes influenced by edge
  for (int node = 0; node< n_nodes; node++)
  {
    retrieve_current_row_and_col(node,row,col);

    // if the node is on an edge, it is an edge node
    if (row == 0 || row == NRows-1 || col == 0 || col == NCols-1)
    {
      edge_node_list.push_back(node);
    }
    else
    {
      bool border_nodata = false;

      // now got through 8 neighbours, looking for nodata
      move_row = row+1;
      move_col = col+1;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row+1;
      move_col = col;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row+1;
      move_col = col-1;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row;
      move_col = col+1;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row;
      move_col = col-1;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row-1;
      move_col = col+1;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row-1;
      move_col = col;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }
      move_row = row-1;
      move_col = col-1;
      if (topography.get_data_element(move_row,move_col) == NoDataValue)
      {
        border_nodata = true;
      }


      // if any of these pixels borders nodata then you add to the list
      if (border_nodata)
      {
        edge_node_list.push_back(node);
      }
    }
  }

  // Okay, now we have an edge node list. We need to make a raster from this
  LSDRaster new_topography = topography;

  int n_nodata_nodes = int(edge_node_list.size());
  int receiver_node, current_node;
  //int receiver_row,receiver_col;
  for (int node = 0; node< n_nodata_nodes; node++)
  {

    current_node = edge_node_list[node];
    retrieve_receiver_information(current_node,receiver_node);
    retrieve_current_row_and_col(current_node,row,col);
    new_topography.set_data_element(row,col,NoDataValue);

    // This does a trace downslope of the node, setting nodatavalues all the way down.
    while (receiver_node != current_node)
    {  // need to do edge checking

      current_node = receiver_node;
      retrieve_receiver_information(current_node,receiver_node);
      retrieve_current_row_and_col(current_node,row,col);
      new_topography.set_data_element(row,col,NoDataValue);

    }

  }

  return new_topography;

}





LSDRaster LSDFlowInfo::find_nodes_not_influenced_by_edge_draining_to_nodelist(vector<int> node_list,LSDRaster& topography)
{


  int raster_NRows, raster_NCols;
  float raster_XMin, raster_YMin, raster_DataRes, raster_NoDataValue;
  int current_node, row,col;

  // first check to make sure the raster dimensions match that of the
  // raster upon which LSDFlowInfo is based
  raster_NRows =  topography.get_NRows();
  raster_NCols =  topography.get_NCols();
  raster_XMin  =  topography.get_XMinimum();
  raster_YMin  =  topography.get_YMinimum();
  raster_NoDataValue = topography.get_NoDataValue();
  raster_DataRes  =  topography.get_DataResolution();

  if (raster_NRows != NRows || raster_NCols != NCols ||
      raster_XMin != XMinimum || raster_YMin != YMinimum ||
      raster_DataRes != DataResolution)
  {
    cout << "Warning!!, LSDFlowInfo::find_nodes_not_influenced_by_edge_draining_to_nodelist\n"
      << "Accumulation raster does not match dimensions of original raster\n"
      << "The resampling need to be done before that step! I will now crash... \n"
      << "I have heard that lsdtopytools does that for you automatically now!" << endl;
    exit(EXIT_FAILURE);
  }

  // create the data array
  Array2D<float> new_topography_data(NRows,NCols,raster_NoDataValue);


  // The node list is loaded from another function. This can be a channel network or a
  // main stream channel. That channel can be generated using the LSDFlowInfo method LSDFlowInfo::get_flow_path
  int n_nodes = int(node_list.size());
  bool problem_node = false;        // for debugging
  cout << "I am trying to isolate nodes draining to a fixed channel that are not influenced by nodata." << endl;
  cout << "I am programmed in a rather rudimentary, brute force way. It may take some time. " << endl;
  cout << "I need to process " << n_nodes << " nodes." << endl;

  if(problem_node)
  {
    cout << "I am in debugging mode! To switch me back you need to go into FlowInfo line 1319" << endl;
  }

  for(int node = 0; node < n_nodes; node++)
  {
    if ( node_list[node] != NoDataValue)
    {
      if (node % 500 == 0)
      {
        cout << "Node " << node << " of " << n_nodes << endl;
      }

      // first get the current node
      retrieve_current_row_and_col(node_list[node],row,col);
      new_topography_data[row][col] = topography.get_data_element(row,col);




      // now get the donors
      vector<int> donors = retrieve_donors_to_node(  node_list[node]);
      int n_donors = int(donors.size());

      if( problem_node)
      {
        cout << "=====================================================================" << endl;
        cout << "Inspecting problem nodes: " << row << ", " << col << "," << new_topography_data[row][col] << endl;
      }

      // Loop through all the donors to this particular node. We will look through the donors to
      // find upslope influence of nodata.
      for (int donor = 0; donor<n_donors; donor++)
      {
        int this_d_node = donors[donor];

        if (problem_node)
        {
          retrieve_current_row_and_col(this_d_node,row,col);
          cout << "donor " << this_d_node << " r: " << row << ", " << col << " ";
        }

        // check to see if the donor is in the fixed channel. If it is, get rid of it.
        bool donor_in_node_list = false;
        for (int nl = 0; nl<n_nodes; nl++)
        {
          if (node_list[nl] == this_d_node)
          {
            donor_in_node_list = true;

            if (problem_node)
            {
              cout << "This donor is in the channel" << endl;
            }
          }
        }

        if (donor_in_node_list == false)
        {
          bool is_influenced = is_upstream_influenced_by_nodata(this_d_node, topography);

          if (problem_node)
          {
            cout << "This donor is not in the channel, Is it influenced by nodata? " << is_influenced << endl;
          }

          // If this donor is not influenced by nodata anywhere upstream, add all the pixels
          // draining to this node to the data.
          if (is_influenced == false)
          {
            vector<int> influenced_nodes = get_upslope_nodes_include_outlet(this_d_node);
            int n_in_basin = int(influenced_nodes.size());
            for (int in = 0; in < n_in_basin; in++)
            {
              current_node = influenced_nodes[in];
              retrieve_current_row_and_col(current_node,row,col);
              new_topography_data[row][col] = topography.get_data_element(row,col);
            }
          }
        }
      }
    }
    else
    {
      cout << "This nodeindex has a nodata value" << endl;
    }

  }

  LSDRaster temp_topo(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,new_topography_data,GeoReferencingStrings);
  return temp_topo;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get the node for a cell at a given row and column
//@author DTM
//@date 08/11/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDFlowInfo::retrieve_node_from_row_and_column(int row, int column)
{
  int Node = NodeIndex[row][column];
  return Node;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// gets a vector of all the donors to a given node
// @author SMM
// @date 19/09/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::retrieve_donors_to_node(int current_node)
{
  // get the numver of donors
  int NDonors = NDonorsVector[current_node];

  // create the vecotr of donating nodes
  vector<int> donorvec(NDonors,NoDataValue);

  // loop over the donating nodes, getting their nodeindicies
  for(int dnode = 0; dnode<NDonors; dnode++)
    {
      donorvec[dnode] = DonorStackVector[ DeltaVector[current_node]+dnode];
    }
  return donorvec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// get the drainage area of a node in km^2
// FJC 06/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDFlowInfo::get_DrainageArea_square_km(int this_node)
{
  int NContributingPixels = NContributingNodes[this_node];
  float DrainageArea = NContributingPixels*DataResolution*DataResolution;
  float DrainageAreaKm = DrainageArea/1000000;

  return DrainageAreaKm;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// get the drainage area of a node in m^2
// FJC 01/05/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDFlowInfo::get_DrainageArea_square_m(int this_node)
{
  int NContributingPixels = NContributingNodes[this_node];
  float DrainageArea = NContributingPixels*DataResolution*DataResolution;

  return DrainageArea;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// recursive add_to_stack routine, from Braun and Willett eq. 12 and 13
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::add_to_stack(int lm_index, int& j_index, int bl_node)
{
  //cout << "j_index: " << j_index << " and s_vec: " << lm_index << endl;

  SVector[j_index] = lm_index;
  BLBasinVector[j_index] = bl_node;
  j_index++;


  int begin_m,end_m;
  int l_index;
  // if donating to itself, need escape hatch
  if ( lm_index == bl_node)
    {
      begin_m = 0;
      end_m = 0;
    }
  else
    {
      begin_m = DeltaVector[lm_index];
      end_m =  DeltaVector[ lm_index+1];
    }
  //cout << "lm_index: " << lm_index << " begin_m: " << begin_m << " end m: " << end_m << endl;
  for( int m_index = begin_m; m_index<end_m; m_index++)
    {
      //cout << "recursion, begin_m: " << begin_m << " and end_m: " << end_m << endl;
      l_index = DonorStackVector[m_index];
      add_to_stack(l_index, j_index, bl_node);
    }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function pickles the data from the flowInfo object into a binary format
// which can be read by the unpickle function later
// the filename DOES NOT include and extension: this is added by the
// function
//
// WARNING: These files are HUGE and testing indicates they don't save much time
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::pickle(string filename)
{
  string ext = ".FIpickle";
  string hdr_ext = ".FIpickle.hdr";

  string hdr_fname = filename+hdr_ext;
  string data_fname = filename+ext;

  ofstream header_out;
  header_out.open(hdr_fname.c_str());

  int contributing_nodes = int(NContributingNodes.size());
  int BLNodes = int(BaseLevelNodeList.size());

  // print the header file
  header_out <<  "ncols             " << NCols
       << "\nnrows             " << NRows
       << "\nxllcorner         " << setprecision(14) << XMinimum
       << "\nyllcorner         " << setprecision(14) << YMinimum
       << "\ncellsize          " << DataResolution
       << "\nNODATA_value      " << NoDataValue
       << "\nNDataNodes        " << NDataNodes
       << "\nNBaseLevelNodes    " << BLNodes
       << "\nNContributingNodes " << contributing_nodes
       << "\nBoundaryConditions ";
  for(int i = 0; i<4; i++)
    {
      header_out << " " << BoundaryConditions[i];
    }
  header_out << endl;
  header_out.close();


  cout << "sizes RC indices: " << RowIndex.size() << " " << ColIndex.size() << endl;
  cout << "BLNL size: " << BaseLevelNodeList.size() << endl;
  cout << "donors: " << NDonorsVector.size() << " Reciev: " << ReceiverVector.size() << endl;
  cout << "delta: " << DeltaVector.size() << " S: " << SVector.size() << endl;
  cout << "donorstack: " << DonorStackVector.size() << " BBasin: " << BLBasinVector.size() << endl;
  cout << "SVectorIndex " << SVectorIndex.size() << " NContrib: " << NContributingNodes.size() << endl;


  // now do the main data
  ofstream data_ofs(data_fname.c_str(), ios::out | ios::binary);
  int temp;
  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      temp = NodeIndex[i][j];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  }

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      temp = FlowDirection[i][j];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  }

  for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
    {
      temp = FlowLengthCode[i][j];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = RowIndex[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = ColIndex[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<BLNodes; i++)
    {
      temp = BaseLevelNodeList[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = NDonorsVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = ReceiverVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes+1; i++)
    {
      temp = DeltaVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = DonorStackVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = SVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = BLBasinVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = SVectorIndex[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<contributing_nodes; i++)
    {
      temp = NContributingNodes[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }

  data_ofs.close();


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this unpickles a pickled flow info object. It is folded into a create function
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::unpickle(string filename)
{

  string ext = ".FIpickle";
  string hdr_ext = ".FIpickle.hdr";

  string hdr_fname = filename+hdr_ext;
  string data_fname = filename+ext;

  ifstream header_in;
  header_in.open(hdr_fname.c_str());

  string temp_str;
  int contributing_nodes;
  vector<string> bc(4);
  int BLNodes;

  header_in >> temp_str >> NCols >> temp_str >> NRows >> temp_str >> XMinimum
      >> temp_str >> YMinimum >> temp_str >> DataResolution
      >> temp_str >> NoDataValue >> temp_str >> NDataNodes
      >> temp_str >> BLNodes
      >> temp_str >> contributing_nodes
      >> temp_str >> bc[0] >> bc[1] >> bc[2] >> bc[3];
  header_in.close();
  BoundaryConditions = bc;


  // now read the data, using the binary stream option
  ifstream ifs_data(data_fname.c_str(), ios::in | ios::binary);
  if( ifs_data.fail() )
    {
      cout << "\nFATAL ERROR: the data file \"" << data_fname
     << "\" doesn't exist" << endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      // initialze the arrays
      Array2D<int> data_array(NRows,NCols,NoDataValue);
      NodeIndex = data_array.copy();
      FlowDirection = data_array.copy();
      FlowLengthCode = data_array.copy();

      vector<int> data_vector(NDataNodes,NoDataValue);
      vector<int> BLvector(BLNodes,NoDataValue);
      vector<int> deltaV(NDataNodes+1,NoDataValue);
      vector<int> CNvec(contributing_nodes,NoDataValue);

      int temp;
      for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
      {
        ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        NodeIndex[i][j] =temp;
      }
  }
      for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
      {
        ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        FlowDirection[i][j] =temp;
      }
  }
      for (int i=0; i<NRows; ++i)
  {
    for (int j=0; j<NCols; ++j)
      {
        ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        FlowLengthCode[i][j] =temp;
      }
  }
      RowIndex = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    RowIndex[i] =temp;

  }
      ColIndex = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    ColIndex[i] =temp;

  }
      BaseLevelNodeList = BLvector;
      for (int i=0; i<BLNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    BaseLevelNodeList[i] =temp;

  }
      NDonorsVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    NDonorsVector[i] =temp;

  }
      ReceiverVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    ReceiverVector[i] =temp;

  }
      DeltaVector = deltaV;
      for (int i=0; i<NDataNodes+1; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    DeltaVector[i] =temp;

  }
      DonorStackVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    DonorStackVector[i] =temp;

  }
      SVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    SVector[i] =temp;

  }
      BLBasinVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    BLBasinVector[i] =temp;

  }
      SVectorIndex = data_vector;
      for (int i=0; i<NDataNodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    SVectorIndex[i] =temp;

  }
      NContributingNodes = CNvec;
      for (int i=0; i<contributing_nodes; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    NContributingNodes[i] =temp;

  }


    }
  ifs_data.close();

  cout << "sizes RC indices: " << RowIndex.size() << " " << ColIndex.size() << endl;
  cout << "BLNL size: " << BaseLevelNodeList.size() << endl;
  cout << "donors: " << NDonorsVector.size() << " Reciev: " << ReceiverVector.size() << endl;
  cout << "delta: " << DeltaVector.size() << " S: " << SVector.size() << endl;
  cout << "donorstack: " << DonorStackVector.size() << " BBasin: " << BLBasinVector.size() << endl;
  cout << "SVectorIndex " << SVectorIndex.size() << " NContrib: " << NContributingNodes.size() << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// SCRIPTS FOR LOADING CSV DATA
// ported from LSDSpatialCSVReader
// FJC 23/03/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a csv file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<string, vector<string> > LSDFlowInfo::load_csv_data(string filename)
{
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv file, but the file" << filename
         << "doesn't exist; check your filename" << endl;
    exit(EXIT_FAILURE);
  }

  // Initiate the data map
  map<string, vector<string> > data_map;

  map<string, int > temp_vec_vec_key;
  vector< vector<string> > temp_vec_vec;
  map<string, vector<string> > temp_data_map;

  // << "Data map size is: " << data_map.size() << endl;
  //cout << "longitude size is: " << longitude.size() << endl;

  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  // get the headers from the first line
  getline(ifs, line_from_file);

  // reset the string vec
  this_string_vec = empty_string_vec;

  // create a stringstream
  stringstream ss(line_from_file);
  ss.precision(9);

  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_string_vec.push_back( substr );
  }
  // now check the data map
  int n_headers = int(this_string_vec.size());
  vector<string> header_vector = this_string_vec;
  for (int i = 0; i<n_headers; i++)
  {
    temp_data_map[header_vector[i]] = empty_string_vec;
  }

  // now loop through the rest of the lines, getting the data.
  while( getline(ifs, line_from_file))
  {
    //cout << "Getting line, it is: " << line_from_file << endl;
    // reset the string vec
    this_string_vec = empty_string_vec;

    // create a stringstream
    stringstream ss(line_from_file);

    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );

      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

      // add the string to the string vec
      this_string_vec.push_back( substr );
    }

    //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
    if ( int(this_string_vec.size()) <= 0)
    {
      cout << "Hey there, I am trying to load your csv data but you seem not to have" << endl;
      cout << "enough columns in your file. I am ignoring a line" << endl;
    }
    else
    {
      int n_cols = int(this_string_vec.size());
      //cout << "N cols is: " << n_cols << endl;
      for (int i = 0; i<n_cols; i++)
      {
        temp_data_map[header_vector[i]].push_back(this_string_vec[i]);
      }
      //cout << "Done with this line." << endl;
    }

  }



  data_map = temp_data_map;

//  cout << "I loaded a csv with the keys: " << endl;
//  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
//  {
//    cout << "Key is: " <<it->first << "\n";
//  }

  return data_map;

}
//==============================================================================
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns the string vector of data from a given column name
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDFlowInfo::get_data_column(string column_name, map<string, vector<string> > data_map)
{
  vector<string> data_vector;
  if ( data_map.find(column_name) == data_map.end() )
  {
    // not found
    cout << "I'm afraid the column "<< column_name << " is not in this dataset" << endl;
  }
  else
  {
    data_vector = data_map[column_name];
  }
  return data_vector;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Converts a data column to a float vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::data_column_to_float(string column_name, map<string, vector<string> > data_map)
{
  vector<string> string_vec = get_data_column(column_name, data_map);
  vector<float> float_vec;
  int N_data_elements = string_vec.size();
  for(int i = 0; i<N_data_elements; i++)
  {
    float_vec.push_back( atof(string_vec[i].c_str()));
  }
  return float_vec;
}

// Converts a data column to a float vector
vector<int> LSDFlowInfo::data_column_to_int(string column_name, map<string, vector<string> > data_map)
{
  vector<string> string_vec = get_data_column(column_name, data_map);
  vector<int> int_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {
    int_vec.push_back( atoi(string_vec[i].c_str()));
  }
  return int_vec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// END OF CSV FUNCTIONS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to ingest the channel heads raster generated using channel_heads_driver.cpp
// into a vector of source nodes so that an LSDJunctionNetwork can be created easily
// from them. Assumes the FlowInfo object has the same dimensions as the channel
// heads raster.
//
// Takes the filename and extension of the channel heads raster.
//
// SWDG 05/12/12
//
// Update: 6/6/14 Happy 3rd birthday Skye!!!!
// SMM
// Now if the file extension is "csv" then the script reads a csv channel heads
// file
//
// Update 30/09/14 Altered structure of function, but key difference is that it
// is now much better in how it goes about reading in channel heads using
// the coordinates, so that channel heads for a region determined usiong one DEM
// can be loaded in to another covering a subsample of the area, or a different
// resolution, which was impossible before.
// DTM
//
// Update 31/03/2017 Now reads the file using csv reader so that columns appear in any order.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::Ingest_Channel_Heads(string filename, string extension, int input_switch)
{

  vector<int> Sources;
  int CH_node;
  map<string, vector<string> > data_map;

  // if this is a csv file, read its contents directly into the node index vector
  if(extension == "csv")
  {
    if(input_switch != 0 && input_switch != 1 && input_switch != 2)
    {
      cout << "\t Note, you have specified an unsupported value for the input switch.  Note: \n\t\t 0=take node index\n\t\t 1=take row and column indices\n\t\t 2=take x and y coordinates"  << endl;
      cout << "\t ...taking node index by default" << endl;
    }

      // load the csv file
    data_map = load_csv_data(filename+".csv");

    vector<int> nodeindex,rowindex,colindex;
    vector<float> x_coord,y_coord;

    nodeindex = data_column_to_int("node", data_map);
    rowindex = data_column_to_int("row", data_map);
    colindex = data_column_to_int("col", data_map);
    x_coord = data_column_to_float("x", data_map);
    y_coord = data_column_to_float("y", data_map);

    int node;
    // use row and column indices to locate source nodes.
    if(input_switch == 1)
    {
      for(int i = 0; i < int(rowindex.size()); ++i)
      {
        if(rowindex[i]<NRows && rowindex[i]>=0 && colindex[i]<NCols && colindex[i] >=0 && NodeIndex[rowindex[i]][colindex[i]]!=NoDataValue)
        {
          node = retrieve_node_from_row_and_column(rowindex[i],colindex[i]);
          Sources.push_back(node);
        }
      }
    }
    // Use coordinates to locate source nodes. Note that this enables the use
    // of LiDAR derived channel heads in coarser DEMs of the same area or
    // subsets of the original DEM for more efficient processing.
    else if(input_switch == 2)
    {
      vector<int> Sources_temp;
      int N_coords = x_coord.size();
      int N_sources_1 = 0;
      for(int i = 0; i < N_coords; ++i)
      {
        node = get_node_index_of_coordinate_point(x_coord[i], y_coord[i]);
        if (node != NoDataValue)
        {
          // Test 1 - Check for channel heads that fall in same pixel
          int test1 = 0;
          N_sources_1 = Sources_temp.size();
          for(int i_test=0; i_test<N_sources_1;++i_test)
          {
            if(node==Sources_temp[i_test]) test1 = 1;
          }
          if(test1==0) Sources_temp.push_back(node);
          //else cout << "\t\t ! removed node from sources list - coincident with another source node" << endl;
        }
      }
      // Test 2 - Need to do some extra checks to load sources correctly.
      int N_sources_2 = Sources_temp.size();
      for(int i = 0; i<N_sources_2; ++i)
      {
        int test2 = 0;
        for(int i_test = 0; i_test<int(Sources_temp.size()); ++i_test)
        {
          if(i!=i_test)
          {
            if(is_node_upstream(Sources_temp[i],Sources_temp[i_test])==true) test2 = 1;
          }
        }
        if(test2 ==0) Sources.push_back(Sources_temp[i]);
        //else cout << "\t\t ! removed node from sources list - other sources upstream" << endl;
      }
    }
    // Using Node Index directly (default)
    else Sources = nodeindex;
  }
  // if not the code assumes a sources raster.
  else
  {
    LSDIndexRaster CHeads(filename, extension);

    for (int i = 0; i < NRows; ++i)
    {
      for (int j = 0; j < NCols; ++j)
      {
        if (CHeads.get_data_element(i,j) != NoDataValue)
        {
          CH_node = retrieve_node_from_row_and_column(i,j);
          if (CH_node != NoDataValue)
          {
            Sources.push_back(CH_node);
          }
        }
      }
    }
  }
  return Sources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Minimalistic method to ingest the channel heads from vectors of x y coordinates
// Using xy allows a "universal" method that can ingest external or internal data
// B.G. 11/11/2018
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::Ingest_Channel_Heads(vector<float>& x_coord, vector<float>& y_coord)
{

  vector<int> Sources, TempSources;
  int node;
  vector<int> Sources_temp;
  int N_coords = x_coord.size();
  int N_sources_1 = 0;
  for(int i = 0; i < N_coords; ++i)
  {
    node = get_node_index_of_coordinate_point(x_coord[i], y_coord[i]);
    if (node != NoDataValue)
    {
      // Test 1 - Check for channel heads that fall in same pixel
      int test1 = 0;
      N_sources_1 = Sources_temp.size();
      for(int i_test=0; i_test<N_sources_1;++i_test)
      {
        if(node==Sources_temp[i_test]) test1 = 1;
      }
      if(test1==0) Sources_temp.push_back(node);
      //else cout << "\t\t ! removed node from sources list - coincident with another source node" << endl;
    }
  }
  // Test 2 - Need to do some extra checks to load sources correctly.
  int N_sources_2 = Sources_temp.size();
  for(int i = 0; i<N_sources_2; ++i)
  {
    int test2 = 0;
    for(int i_test = 0; i_test<int(Sources_temp.size()); ++i_test)
    {
      if(i!=i_test)
      {
        if(is_node_upstream(Sources_temp[i],Sources_temp[i_test])==true) test2 = 1;
      }
    }
    if(test2 ==0) Sources.push_back(Sources_temp[i]);
    //else cout << "\t\t ! removed node from sources list - other sources upstream" << endl;
  }
  return Sources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to ingest the channel heads raster generated using channel_heads_driver.cpp
// into a vector of source nodes so that an LSDJunctionNetwork can be created easily
// from them. Assumes the FlowInfo object has the same dimensions as the channel
// heads raster.
//
// Takes the filename and extension of the channel heads raster.
//
// SWDG 05/12/12
//
// Update: 6/6/14 Happy 3rd birthday Skye!!!!
// SMM
// Now if the file extension is "csv" then the script reads a csv channel heads
// file
//
// Update 30/09/14 Altered structure of function, but key difference is that it
// is now much better in how it goes about reading in channel heads using
// the coordinates, so that channel heads for a region determined usiong one DEM
// can be loaded in to another covering a subsample of the area, or a different
// resolution, which was impossible before.
// DTM
//
// *******************************************************************************************
// UPDATE 23/03/17 - NEW OVERLOADED CHANNEL HEADS INGESTION ROUTINE.  This ONLY works with the
// csv file as this seems to be the best way of reading in the channel heads. Other formats
// should now be obsolete.
// Finds the appropriate column from the csv based on the string of the heading rather
// than by column number, so should work with different versions of the output sources
// csv file.
//
// Input switch tells what the columns the code should be looking for:
// 0 - use the node index
// 1 - use rows and columns
// 2 - use x and y (UTM coordinates)

// Could add in a 3rd switch for lat long but this requires a bunch of extra porting that
// I can't be bothered to do right now.
// FJC
// *****************************************************************************************
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::Ingest_Channel_Heads(string filename, int input_switch)
{
  vector<int> Sources;
  //int CH_node;

  // load the csv file
  map<string, vector<string> > data_map = load_csv_data(filename+".csv");


  // now check the input switch to search for the various columns
  if (input_switch == 0)
  {
    // use the node index
    vector<int> NodeIndices = data_column_to_int("node", data_map);
    Sources = NodeIndices;
  }
  else if (input_switch == 1)
  {
    // use the rows and columns
    vector<int> rows = data_column_to_int("row", data_map);
    vector<int> cols = data_column_to_int("col", data_map);
    for (int i = 0; i < int(rows.size()); i++)
    {
      int NI = retrieve_node_from_row_and_column(rows[i], cols[i]);
      Sources.push_back(NI);
    }
  }
  else if (input_switch == 2)
  {
    // use x and y (UTM coordinates)
    vector<float> x_coord = data_column_to_float("x", data_map);
    vector<float> y_coord = data_column_to_float("y", data_map);
    int N_coords = x_coord.size();

    vector<int> Sources_temp;
    int N_sources_1 = 0;
    for(int i = 0; i < N_coords; ++i)
    {
      int node = get_node_index_of_coordinate_point(x_coord[i], y_coord[i]);
      if (node != NoDataValue)
      {
        // Test 1 - Check for channel heads that fall in same pixel
        int test1 = 0;
        N_sources_1 = Sources_temp.size();
        for(int i_test=0; i_test<N_sources_1;++i_test)
        {
          if(node==Sources_temp[i_test]) test1 = 1;
        }
        if(test1==0) Sources_temp.push_back(node);
        else cout << "\t\t ! removed node from sources list - coincident with another source node" << endl;
      }
    }
    // Test 2 - Need to do some extra checks to load sources correctly.
    int N_sources_2 = Sources_temp.size();
    for(int i = 0; i<N_sources_2; ++i)
    {
      int test2 = 0;
      for(int i_test = 0; i_test<int(Sources_temp.size()); ++i_test)
      {
        if(i!=i_test)
        {
          if(is_node_upstream(Sources_temp[i],Sources_temp[i_test])==true) test2 = 1;
        }
      }
      if(test2 ==0) Sources.push_back(Sources_temp[i]);
      else cout << "\t\t ! removed node from sources list - other sources upstream" << endl;
    }
  }
  else
  {
    cout << "You have not supplied a valid input switch! Please supply either 0, 1, or 2." << endl;
  }
  return Sources;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to ingest sources from OS MasterMap Water Network Layer (csv)
// into a vector of source nodes so that an LSDJunctionNetwork can be created easily
// from them.
//
// Takes the filename and extension of the channel heads raster.
//
// FJC 28/11/2016
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::Ingest_Channel_Heads_OS(string csv_filename)
{
  vector<int> Sources;
  // read in the CSV file
  ifstream input_csv;
  string dot = ".";
  string extension = "csv";
  string fname = csv_filename+dot+extension;

  cout << "The CSV filename is: " << fname << endl;

  input_csv.open(fname.c_str());
  // check for correct input
  if (input_csv.good() == false)
  {
    cout << "I can't read the CSV file! Check your filename." << endl;
  }

  //int object_ID, name, source_ID, PosAlong;
  //float X,Y;
  vector<float> X_coords, Y_coords;

  // read in the file
  while(!input_csv.eof())
  {
    string line;
    getline(input_csv,line);

    // get the x and y coords to vectors
    istringstream ss(line);
    string param;
    int i=0;
    while(getline(ss, param, ','))
    {
      if (i == 4)
      {
        X_coords.push_back(atof(param.c_str()));
      }
      if (i == 5)
      {
        Y_coords.push_back(atof(param.c_str()));
      }
      i++;
    }
  }

  vector<int> Sources_temp;
  int N_coords = X_coords.size();
  int N_sources_1 = 0;
  for(int i = 0; i < N_coords; ++i)
  {
		int node = get_node_index_of_coordinate_point(X_coords[i], Y_coords[i]);
    if (node != NoDataValue)
    {
      // Test 1 - Check for channel heads that fall in same pixel
      int test1 = 0;
      N_sources_1 = Sources_temp.size();
      for(int i_test=0; i_test<N_sources_1;++i_test)
      {
        if(node==Sources_temp[i_test]) test1 = 1;
      }
      if(test1==0) Sources_temp.push_back(node);
      else cout << "\t\t ! removed node from sources list - coincident with another source node" << endl;
    }
  }
  // Test 2 - Need to do some extra checks to load sources correctly.
  int N_sources_2 = Sources_temp.size();
  for(int i = 0; i<N_sources_2; ++i)
  {
		cout << flush << "\t Source: " << i << " of " << N_sources_2 << "\r";
    int test2 = 0;
    for(int i_test = 0; i_test<int(Sources_temp.size()); ++i_test)
    {
			if(i!=i_test)
      {
      	if(is_node_upstream(Sources_temp[i],Sources_temp[i_test])==true) test2 = 1;
      }
    }
    if(test2 ==0) Sources.push_back(Sources_temp[i]);
    //else cout << "\t\t ! removed node from sources list - other sources upstream" << endl;
  }
	cout << "Returning sources..." << endl;

  return Sources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Get the tips of the channel network
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::get_channel_tip_raster(vector<int> sources,int n_pixels_in_tips)
{
  Array2D<int> tips(NRows,NCols,NoDataValue);

  int curr_row, curr_col, this_node;
  int receiver_node;
  int n_sources = int(sources.size());
  for (int i = 0; i<n_sources; i++)
  {
    this_node = sources[i];
    retrieve_current_row_and_col(sources[i],curr_row,curr_col);
    // tag the source
    tips[curr_row][curr_col]= 1;

    // now go down the tip the appropriate number of nodes
    for (int t = 0 ; t<n_pixels_in_tips-1; t++)
    {
      retrieve_receiver_information(this_node,receiver_node,curr_row,curr_col);  
      tips[curr_row][curr_col]= 1;
      this_node = receiver_node;    
    }

  }
  LSDIndexRaster temp_tips(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,tips,GeoReferencingStrings);
  return temp_tips;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the flownet information
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_flow_info_vectors(string filename)
{
  string string_filename;
  string dot = ".";
  string extension = "txt";
  string_filename = filename+dot+extension;
  cout << "The filename is " << string_filename << endl;

  // print out all the donor, receiver and stack info
  ofstream donor_info_out;
  donor_info_out.open(string_filename.c_str());
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << i << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << ReceiverVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << NDonorsVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes+1; i++)
    {
      donor_info_out << DeltaVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << DonorStackVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << SVector[i] << " ";
    }
  donor_info_out << endl;

  if( int(SVectorIndex.size()) == NDataNodes)
    {
      for(int i = 0; i<NDataNodes; i++)
  {
    donor_info_out << SVectorIndex[i] << " ";
  }
      donor_info_out << endl;
      for(int i = 0; i<NDataNodes; i++)
  {
    donor_info_out << NContributingNodes[i] << " ";
  }
      donor_info_out << endl;
    }

  donor_info_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// these functions write the index arrays to index rasters
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndex_to_LSDIndexRaster()
{
  cout << "NRows: " << NRows << " and NCols: " << NCols << endl;
  LSDIndexRaster temp_nodeindex(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,NodeIndex,GeoReferencingStrings);
  return temp_nodeindex;
}

LSDIndexRaster LSDFlowInfo::write_FlowDirection_to_LSDIndexRaster()
{
  LSDIndexRaster temp_flowdir(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowDirection,GeoReferencingStrings);
  return temp_flowdir;
}

LSDIndexRaster LSDFlowInfo::write_FlowLengthCode_to_LSDIndexRaster()
{
  LSDIndexRaster temp_flc(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowLengthCode,GeoReferencingStrings);
  return temp_flc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes an LSDIndexRaster given a list of node indices
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster(vector<int>& nodeindexvec)
{
  int n_node_indices = nodeindexvec.size();
  //cout << "The number of nodeindices is: " << n_node_indices << endl;
  Array2D<int> chan(NRows,NCols,NoDataValue);

  //cout << "Raster nr: " << chan.dim1() << " nc: " << chan.dim2() << endl;

  int curr_row, curr_col;

  for(int i = 0; i<n_node_indices; i++)
    {
      // make sure there is no segmentation fault for bad data
      // Note: bad data is ignored
      if(nodeindexvec[i] <= NDataNodes)
  {
    retrieve_current_row_and_col(nodeindexvec[i],curr_row,
               curr_col);

    if(chan[curr_row][curr_col] == NoDataValue)
      {
        chan[curr_row][curr_col] = 1;
      }
    else
      {
        chan[curr_row][curr_col]++;
      }
  }
      else
  {
    cout << "WARNING: LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster"
         << " node index does not exist!"<< endl;
  }
    }


  LSDIndexRaster temp_chan(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,chan,GeoReferencingStrings);
  return temp_chan;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes an LSDIndesxRaster given a list of node indices, and give every
// pixel its nodeindex value, which is unique.
//
// SWDG after SMM 2/2/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster_Unique(vector<int>& nodeindexvec)
{
  int n_node_indices = nodeindexvec.size();
  Array2D<int> chan(NRows,NCols,NoDataValue);

  int curr_row, curr_col;

  for(int i = 0; i<n_node_indices; i++){
    // make sure there is no segmentation fault for bad data
    // Note: bad data is ignored
    if(nodeindexvec[i] <= NDataNodes){
      retrieve_current_row_and_col(nodeindexvec[i], curr_row, curr_col);

      if(chan[curr_row][curr_col] == NoDataValue){
        chan[curr_row][curr_col] = i;
      }
      else{
        //revisted points will be overwritten, most recent id will be pres
        chan[curr_row][curr_col] = i;
      }
    }
    else
    {
      cout << "WARNING: LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster"
         << " node index does not exist!"<< endl;
    }
  }

  LSDIndexRaster temp_chan(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,chan,GeoReferencingStrings);
  return temp_chan;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function calculates the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
//
// This function used the S vector index and NContributing nodes to calculate the area
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NContributingNodes_to_LSDIndexRaster()
{
  Array2D<int> contributing_pixels(NRows,NCols,NoDataValue);
  int row,col;

  // loop through the node vector, adding pixels to receiver nodes
  for(int node = 0; node<NDataNodes; node++)
  {
    row = RowIndex[node];
    col = ColIndex[node];
    contributing_pixels[row][col] = NContributingNodes[node];
  }

  LSDIndexRaster temp_cp(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,contributing_pixels,GeoReferencingStrings);
  return temp_cp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// writes a index raster in arc format
// LSD format:
// 7  0 1
// 6 -1 2
// 5  4 3
//  int Arc_flowdir;      // flow direction in arcmap format
//                // 32  64  128
//                // 16   0  1
//                // 8    4  2
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_FlowDirection_to_LSDIndexRaster_Arcformat()
{

  Array2D<int> FlowDirectionArc(NRows,NCols,NoDataValue);
  for(int row = 0; row<NRows; row++)
    {
      for (int col = 0; col<NCols; col++)
  {
    if ( FlowDirection[row][col] == -1)
      {
        FlowDirectionArc[row][col] = 0;
      }
    else if ( FlowDirection[row][col] == 0)
      {
        FlowDirectionArc[row][col] = 64;
      }
    else if ( FlowDirection[row][col] == 1)
      {
        FlowDirectionArc[row][col] = 128;
      }
    else if ( FlowDirection[row][col] == 2)
      {
        FlowDirectionArc[row][col] = 1;
      }
    else if ( FlowDirection[row][col] == 3)
      {
        FlowDirectionArc[row][col] = 2;
      }
    else if ( FlowDirection[row][col] == 4)
      {
        FlowDirectionArc[row][col] = 4;
      }
    else if ( FlowDirection[row][col] == 5)
      {
        FlowDirectionArc[row][col] = 8;
      }
    else if ( FlowDirection[row][col] == 6)
      {
        FlowDirectionArc[row][col] = 16;
      }
    else if ( FlowDirection[row][col] == 7)
      {
        FlowDirectionArc[row][col] = 32;
      }
  }
    }

  LSDIndexRaster temp_fd(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowDirectionArc,GeoReferencingStrings);
  return temp_fd;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes the drainage area (number of contributing nodes
// * DataResolution^2) to an LSDRaster object
// Added by FC 15/11/12
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::write_DrainageArea_to_LSDRaster()
{
  // initialise the 2D array
  int row,col;                // node index
  float ndv = float(NoDataValue);
  float this_DA;
  Array2D<float> DrainageArea_local(NRows,NCols,ndv);

  //cout << "HELLO WHT THE FUCK" << endl;

  for(int node = 0; node<NDataNodes; node++)
    {
      if (RowIndex[node] != ndv && ColIndex[node] != ndv)
      {
        row = RowIndex[node];
        col = ColIndex[node];
        //cout << row << " " << col << endl;
        //cout << NContributingNodes[node] << endl;
        this_DA = float(NContributingNodes[node])*DataResolution*DataResolution;
        //cout << "this da" << this_DA << endl;
        DrainageArea_local[row][col] = this_DA;
      }
    }

  // create the LSDRaster object
  LSDRaster DrainageArea(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,DrainageArea_local,GeoReferencingStrings);
  return DrainageArea;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
// In this function a pixel that has no donors has a contributing pixel value of 0
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::calculate_n_pixels_contributing_from_upslope()
{
  Array2D<int> contributing_pixels(NRows,NCols,NoDataValue);
  int row,col;
  int receive_row, receive_col;
  int receiver_node;

  // loop through the s vector, adding pixels to receiver nodes
  for(int node = NDataNodes-1; node>=0; node--)
  {

    row = RowIndex[SVector[node]];
    col = ColIndex[SVector[node]];
    // if the pixel exists and has no contributing pixels,
    // change from nodata to zero

    if(contributing_pixels[row][col] == NoDataValue)
    {
      contributing_pixels[row][col] = 0;
    }

    receiver_node = ReceiverVector[ SVector[node] ] ;
    receive_row = RowIndex[ receiver_node ];
    receive_col = ColIndex[ receiver_node ];

    //cout << "node " << node << " pixel: " << SVector[node] << " receiver: " << receiver_node << endl;
    //cout << "contributing: " << contributing_pixels[row][col] << endl;

    if ( receiver_node  == SVector[node])
    {
      // do nothing
    }
    else if ( contributing_pixels[receive_row][receive_col] == NoDataValue)
    {
      contributing_pixels[receive_row][receive_col] =
      contributing_pixels[row][col]+1;
    }
    else
    {
      contributing_pixels[receive_row][receive_col] +=
      contributing_pixels[row][col]+1;
    }

    //cout << "recieving: " << contributing_pixels[receive_row][receive_col] << endl;
  }

  LSDIndexRaster temp_cp(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,contributing_pixels,GeoReferencingStrings);
  return temp_cp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
//
// In this function a pixel that has no donors contributes its own flow
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::calculate_upslope_reference_indices()
{
  vector<int> vectorized_area(NDataNodes,1);
  SVectorIndex = vectorized_area;

  int receiver_node;
  int donor_node;

  // loop through the s vector, adding pixels to receiver nodes
  for(int node = NDataNodes-1; node>=0; node--)
  {
    donor_node = SVector[node];
    receiver_node = ReceiverVector[ donor_node ];

    // every node is visited once and only once so we can map the
    // unique positions of the nodes to the SVector
    SVectorIndex[donor_node] = node;

    // add the upslope area (note no action is taken
    // for base level nodes since they donate to themselves and
    // we must avoid float counting
    if (donor_node != receiver_node)
    {
      vectorized_area[ receiver_node ] +=  vectorized_area[ donor_node ];
    }
  }

  NContributingNodes = vectorized_area;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function returns a integer vector containing all the node numbers upslope
// of of the node with number node_number_outlet
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::get_upslope_nodes(int node_number_outlet)
{
  vector<int> us_nodes;

  if(node_number_outlet < 0 || node_number_outlet > NDataNodes-1)
    {
      cout << "the node index does not exist" << endl;
      exit(EXIT_FAILURE);
    }

  int start_SVector_node = SVectorIndex[node_number_outlet];
  int end_SVector_node = start_SVector_node+NContributingNodes[node_number_outlet];

  for(int node = start_SVector_node; node < end_SVector_node; node++)
  {
    us_nodes.push_back(SVector[node]);
  }

  return us_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function returns a integer vector containing all the node numbers upslope
// of of the node with number node_number_outlet, and includes the outlet
//
// SMM 30/07/2019
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::get_upslope_nodes_include_outlet(int node_number_outlet)
{
  vector<int> us_nodes;
  us_nodes.push_back(node_number_outlet);

  if(node_number_outlet < 0 || node_number_outlet > NDataNodes-1)
  {
    cout << "the node index does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  int start_SVector_node = SVectorIndex[node_number_outlet];
  int end_SVector_node = start_SVector_node+NContributingNodes[node_number_outlet];

  for(int node = start_SVector_node; node < end_SVector_node; node++)
  {
    us_nodes.push_back(SVector[node]);
  }

  return us_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of source nodes and creates a raster where
// the pixels have a value of 1 where there are upslope nodes and nodata otherwise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_node_mask(vector<int> source_nodes)
{
  // initiate the data array
  Array2D<float> this_raster(NRows,NCols,NoDataValue);

  // loop through the nodes, collecting upslope nodes
  int n_nodes = int(source_nodes.size());
  int this_node;
  //int us_node;
  int curr_row,curr_col;
  float is_us = 1.0;

  // go through all the source nodes, find their upslope nodes
  // and set the value of these nodes to 1.0 on the data array
  for (int n = 0; n<n_nodes; n++)
  {
    this_node = source_nodes[n];

    // check if it is in the DEM
    if (this_node < NDataNodes)
    {
      vector<int> upslope_nodes = get_upslope_nodes(this_node);

      int n_us_nodes =  int(upslope_nodes.size());
      for(int us = 0; us<n_us_nodes; us++)
      {
        retrieve_current_row_and_col(upslope_nodes[us],curr_row,curr_col);
        this_raster[curr_row][curr_col]=is_us;
      }
    }
  }

  // now create the raster
  LSDRaster temp_us(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,this_raster,GeoReferencingStrings);
  return temp_us;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of source nodes and creates a raster where
// the pixels have a value of upslope_value
// where there are upslope nodes and NoData otherwise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_node_mask(vector<int> source_nodes, vector<float> upslope_values)
{
  // initiate the data array
  Array2D<float> this_raster(NRows,NCols,NoDataValue);

  if (source_nodes.size() == upslope_values.size())
  {
    // loop through the nodes, collecting upslope nodes
    int n_nodes = int(source_nodes.size());
    int this_node;
    //int us_node;
    int curr_row,curr_col;

    // go through all the source nodes, find their upslope nodes
    // and set the value of these nodes to 1.0 on the data array
    for (int n = 0; n<n_nodes; n++)
    {
      this_node = source_nodes[n];

      // check if it is in the DEM
      if (this_node < NDataNodes)
      {
        vector<int> upslope_nodes = get_upslope_nodes(this_node);

        int n_us_nodes =  int(upslope_nodes.size());
        for(int us = 0; us<n_us_nodes; us++)
        {
          retrieve_current_row_and_col(upslope_nodes[us],curr_row,curr_col);
          this_raster[curr_row][curr_col]= upslope_values[n];
        }
      }
    }
  }
  else
  {
    cout << "The uplsope vlaues vector needs to be the same lengths as the sources vector!" << endl;
    cout << "Returning an nodata raster" << endl;
  }

  // now create the raster
  LSDRaster temp_us(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,this_raster,GeoReferencingStrings);
  return temp_us;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//  Accumulate some variable (such a precipitation) from an accumulation raster
//
//  This requires summing all upslope nodes for every node. It seems a bit inefficient
//  but the other simple alternative is to do a sort() operation initially and then
//  move from upslope node down. There is probably a more efficient way to do this
//  and this algorithm should be revisited later to see if we can speed it up.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::upslope_variable_accumulator(LSDRaster& accum_raster)
{
  cout << "You are using the old and slow upslope accumulator. " << endl;
  cout << "You should probably use upslope_variable_accumulator_v2" << endl;
  cout << "This version does not accumulate at the base level nodes" << endl;

  int raster_NRows, raster_NCols;
  float raster_XMin, raster_YMin, raster_DataRes;

  // first check to make sure the raster dimensions match that of the
  // raster upon which LSDFlowInfo is based
  raster_NRows =  accum_raster.get_NRows();
  raster_NCols =  accum_raster.get_NCols();
  raster_XMin  =  accum_raster.get_XMinimum();
  raster_YMin  =  accum_raster.get_YMinimum();
  raster_DataRes  =  accum_raster.get_DataResolution();

  if (raster_NRows != NRows || raster_NCols != NCols ||
      raster_XMin != XMinimum || raster_YMin != YMinimum ||
      raster_DataRes != DataResolution)
  {
    cout << "Warning!!, LSDFlowInfo::upslope_area_accumulator\n"
         << "Accumulation raster does not match dimensions of original raster" << endl;
    exit(EXIT_FAILURE);
    // return accum_raster;
  }
  else
  {
    // create the data array
    Array2D<float> accumulated_data_array(NRows,NCols,NoDataValue);

    // loop through all the nodes, accumulating the areas
    for(int this_node = 0; this_node <NDataNodes; this_node++)
    {
      // get the upslope nodes
      vector<int> node_vec = get_upslope_nodes(this_node);

      // loop through these nodes, adding them to the accumulator
      double this_node_accumulated = 0;
      int this_row, this_col;
      for (int ni = 0; ni<int(node_vec.size()); ni++)
      {
        retrieve_current_row_and_col(node_vec[ni],this_row,this_col);
        this_node_accumulated += double(accum_raster.get_data_element(this_row, this_col));
      }

      // write the accumulated variable to the array
      int curr_row, curr_col;
      retrieve_current_row_and_col(this_node,curr_row,curr_col);
      accumulated_data_array[curr_row][curr_col] = float(this_node_accumulated);
  }
  // create the raster
  LSDRaster accumulated_flow(NRows, NCols, XMinimum, YMinimum,
            DataResolution, NoDataValue, accumulated_data_array,GeoReferencingStrings);
  return accumulated_flow;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


void LSDFlowInfo::upslope_variable_accumulator_v2(LSDRaster& accum_raster)
{
  this->upslope_variable_accumulator_v2(accum_raster, true);
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//  Accumulate some variable (such a precipitation) from an accumulation raster
//
//  v2, I am trying to optimise it and remove the dependencies of having the exact same dimensions!
//
// SMM update 30/07/2019 Making sure this reflects nodata and getting rid of the
//  topography raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::upslope_variable_accumulator_v2(LSDRaster& accum_raster, bool accum_current_node)
{
  int raster_NRows, raster_NCols;
  float raster_XMin, raster_YMin, raster_DataRes, raster_NoDataValue;

  // first check to make sure the raster dimensions match that of the
  // raster upon which LSDFlowInfo is based
  raster_NRows =  accum_raster.get_NRows();
  raster_NCols =  accum_raster.get_NCols();
  raster_XMin  =  accum_raster.get_XMinimum();
  raster_YMin  =  accum_raster.get_YMinimum();
  raster_NoDataValue = accum_raster.get_NoDataValue();
  raster_DataRes  =  accum_raster.get_DataResolution();

  if (raster_NRows != NRows || raster_NCols != NCols ||
      raster_XMin != XMinimum || raster_YMin != YMinimum ||
      raster_DataRes != DataResolution)
  {
    cout << "Warning!!, LSDFlowInfo::upslope_area_accumulator\n"
      << "Accumulation raster does not match dimensions of original raster\n"
      << "The resampling need to be done before that step! I will now crash... \n"
      << "I have heard that lsdtopytools does that for you automatically now!" << endl;
    exit(EXIT_FAILURE);
  }

  // getting my clean vectors
  map<string, vector< vector<int> > > these_vec = get_map_of_vectors();

  // create the data array
  Array2D<float> accumulated_data_array(NRows,NCols,raster_NoDataValue);
  // Array2D<float> checker(NRows,NCols,goulg);
  for(int row=0;row<NRows;row++)
  {
    for(int col=0;col<NCols;col++)
    {
      if(accum_raster.get_data_element(row,col)!=raster_NoDataValue)
      {
        accumulated_data_array[row][col] = 0;
      }
    }
  }
  cout << "Initiated the raster" << endl;

  // I will comment that code to make sure that I am not misunderstanding things
  // I am getting the number of nodes here
  size_t n_vec = these_vec["inverted_stack_order"].size();
  int tot_n = 0;
  int tot_n_check = 0;
  // Iterating through the inverted stack: the top of the hills will be processed before their downslope nodes
  // Note that my these_vec is a map of vector of vectors where each vectors is a stack/inverted stack for a watershed
  for(int i=0; i<int(n_vec); i++)
  {
    // Getting the current nodes of that watershed
    vector<int> vec_of_node = these_vec["inverted_stack_order"][i];
    // its size
    size_t dat_vec_size = vec_of_node.size();
    // Getting the numer of node I have processed in total
    tot_n += int(dat_vec_size);
    // Iterating through the nodes of that particular watersed (inverted stack)
    for(size_t j=0; j<dat_vec_size; j++)
    {
      // current node
      int this_node = vec_of_node[j];
      // Getting the row and col indices
      int row,col;
      tot_n_check++;
      row =  these_vec["inverted_rows"][i][j];
      col = these_vec["inverted_cols"][i][j];

      if (row <0 || row > NRows-1 || col < 0 || col > NCols-1)
      {
        cout << "Problem in accumulator, I'm about to crash" << endl;
        cout << "row: " << row << " col: " << col << endl;
      }

      // First getting the accumulation of that particular current node <--- THIS MIGHT BE WHERE THE DIFFEREENCE IS SIMON
      if(accum_current_node)
      {
        float increment_val = accum_raster.get_data_element(row,col);
        accumulated_data_array[row][col] += increment_val;
      }
      // if my receiving node is not the same (i.e. if I am not a base level), I will transmit that accumulation to the receiving node
      // first getting the receiving node
      int recnode;
      retrieve_receiver_information(this_node,recnode);
      // Checking if I am not a base level
      if(recnode != this_node)
      {
        int rrow,rcol;
        // Getting the location of my receiving node
        retrieve_current_row_and_col(recnode,rrow,rcol);
        // Adding my accumulated data to that next node
        accumulated_data_array[rrow][rcol] += accumulated_data_array[row][col];
      }

      // Going to the next node
    }
    // going to the next basin
  }
  // // Done

  // create the raster
  LSDRaster accumulated_flow(NRows, NCols, XMinimum, YMinimum,
        DataResolution, NoDataValue, accumulated_data_array,GeoReferencingStrings);
  return accumulated_flow;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//  Accumulate some variable (such a precipitation) from an accumulation raster
//
//  v3, trying to simplify the code (SMM, 19/11/2019)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::upslope_variable_accumulator_v3(LSDRaster& accum_raster)
{
  int raster_NRows, raster_NCols;
  float raster_XMin, raster_YMin, raster_DataRes, raster_NoDataValue;

  // first check to make sure the raster dimensions match that of the
  // raster upon which LSDFlowInfo is based
  raster_NRows =  accum_raster.get_NRows();
  raster_NCols =  accum_raster.get_NCols();
  raster_XMin  =  accum_raster.get_XMinimum();
  raster_YMin  =  accum_raster.get_YMinimum();
  raster_NoDataValue = accum_raster.get_NoDataValue();
  raster_DataRes  =  accum_raster.get_DataResolution();

  if (raster_NRows != NRows || raster_NCols != NCols ||
      raster_XMin != XMinimum || raster_YMin != YMinimum ||
      raster_DataRes != DataResolution)
  {
    cout << "Warning!!, LSDFlowInfo::upslope_area_accumulator\n"
      << "Accumulation raster does not match dimensions of original raster\n"
      << "The resampling need to be done before that step! I will now crash... \n"
      << "I have heard that lsdtopytools does that for you automatically now!" << endl;
    exit(EXIT_FAILURE);
  }

  // We first loop through the the nodes, adding the value on that particular pixel
  int row,col;
  vector<float> vectorized_area(NDataNodes,raster_NoDataValue);
  for (int node = 0; node<NDataNodes; node++)
  {
    row = RowIndex[node];
    col = ColIndex[node];
    vectorized_area[node] = accum_raster.get_data_element(row,col);
  }

  int receiver_node;
  int donor_node;

  // loop through the s vector, adding pixels to receiver nodes
  for(int node = NDataNodes-1; node>=0; node--)
  {
    donor_node = SVector[node];
    receiver_node = ReceiverVector[ donor_node ];

    // add the upslope area (note no action is taken
    // for base level nodes since they donate to themselves and
    // we must avoid float counting
    if (donor_node != receiver_node)
    {
      vectorized_area[ receiver_node ] +=  vectorized_area[ donor_node ];
    }
  }

  // create the data array
  Array2D<float> accumulated_data_array(NRows,NCols,raster_NoDataValue);
  for (int node = 0; node<NDataNodes; node++)
  {
    row = RowIndex[node];
    col = ColIndex[node];
    accumulated_data_array[row][col] = vectorized_area[node];
  }

  // create the raster
  LSDRaster accumulated_flow(NRows, NCols, XMinimum, YMinimum,
        DataResolution, NoDataValue, accumulated_data_array,GeoReferencingStrings);
  return accumulated_flow;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//  This averages an upslope quantity
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::upslope_average(LSDRaster& accum_raster)
{
  LSDRaster Accumulated_value = upslope_variable_accumulator_v3(accum_raster);
  float raster_NoDataValue = accum_raster.get_NoDataValue();
  int row,col;
  // Now loop through nodes getting the average values  
  // create the data array
  Array2D<float> average_data_array(NRows,NCols,raster_NoDataValue);
  for (int node = 0; node<NDataNodes; node++)
  {
    row = RowIndex[node];
    col = ColIndex[node];

    average_data_array[row][col] = Accumulated_value.get_data_element(row,col) / float( NContributingNodes[node] );
  }  


  // create the raster
  LSDRaster average_value(NRows, NCols, XMinimum, YMinimum,
        DataResolution, NoDataValue, average_data_array,GeoReferencingStrings);
  return average_value;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function tests whether the test is upstream of the current node
//
// FC 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::is_node_upstream(int current_node, int test_node)
{
  int i = 0;

  int start_SVector_node = SVectorIndex[current_node];
  int end_SVector_node = start_SVector_node+NContributingNodes[current_node];

  int SVector_test_node = SVectorIndex[test_node];

  for(int node = start_SVector_node; node < end_SVector_node; node++)
  {
    if (node == SVector_test_node)
    {
      i = 1;
    }
  }

  return i;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function tests whether a node is a base level node
//
// FC 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::is_node_base_level(int node)
{
  int i = 0;

  for (int j = 0; j < int(BaseLevelNodeList.size()); j++)
  {
    if (node == BaseLevelNodeList[j])
    {
      i = 1;
    }
  }

  return i;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function redurns a vector of node indices to all the donor
// nodes of a particular node
//
// SMM 21/10/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_donor_nodes(int current_node)
{
  int start_D = DeltaVector[current_node];
  int end_D = DeltaVector[current_node+1];

  vector<int> donor_nodes;
  for(int this_node = start_D; this_node<end_D; this_node++)
  {
    //cout << "node " << current_node << " and donor: " << DonorStackVector[ this_node ] << endl;
    donor_nodes.push_back( DonorStackVector[ this_node ] );
  }

  return donor_nodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the chi function for all the nodes upslope a given node
// it takes a node list
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::get_upslope_chi(int starting_node, float m_over_n, float A_0)
{
  vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);
  vector<float> chi_vec = get_upslope_chi(upslope_pixel_list, m_over_n, A_0);
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the chi function for all the nodes upslope a given node
// it takes a node list
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::get_upslope_chi(int starting_node, float m_over_n, float A_0,
                                           LSDRaster& Discharge)
{
  vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);
  vector<float> chi_vec = get_upslope_chi(upslope_pixel_list, m_over_n, A_0, Discharge);
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is called from the the get_upslope_chi that only has an integer
// it returns the acutal chi values in a vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::get_upslope_chi(vector<int>& upslope_pixel_list,
                                           float m_over_n, float A_0)
{

  int receiver_node;
  int IndexOfReceiverInUplsopePList;
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  float pixel_area = DataResolution*DataResolution;
  int node,row,col;
  // get the number of nodes upslope
  int n_nodes_upslope = upslope_pixel_list.size();
  vector<float> chi_vec(n_nodes_upslope,0.0);

  if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
  {
    cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
    exit(EXIT_FAILURE);
  }

  int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];

  for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
  {
    node = upslope_pixel_list[n_index];
    receiver_node = ReceiverVector[ node ];
    IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
    row = RowIndex[node];
    col = ColIndex[node];

    if (FlowLengthCode[row][col] == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }


    chi_vec[n_index] = dx*(pow( (A_0/ (float(NContributingNodes[node])*pixel_area) ),m_over_n))
                          + chi_vec[IndexOfReceiverInUplsopePList];

    //  cout << "node: " << upslope_pixel_list[n_index] << " receiver: " << receiver_node
    //       << " SIndexReciever: " << IndexOfReceiverInUplsopePList
    //       << " and checked: " << upslope_pixel_list[IndexOfReceiverInUplsopePList]
    //       << " and chi: " << chi_vec[n_index] << endl;

  }
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is called from the the get_upslope_chi that only has an integer
// it returns the actual chi values in a vector
// same as above but uses a discharge raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::get_upslope_chi(vector<int>& upslope_pixel_list,
                                           float m_over_n, float A_0,
                                           LSDRaster& Discharge)
{

  int receiver_node;
  int IndexOfReceiverInUplsopePList;
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  int node,row,col;
  // get the number of nodes upslope
  int n_nodes_upslope = upslope_pixel_list.size();
  vector<float> chi_vec(n_nodes_upslope,0.0);

  if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
  {
    cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
    exit(EXIT_FAILURE);
  }

  int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];

  for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
  {
    node = upslope_pixel_list[n_index];
    receiver_node = ReceiverVector[ node ];
    IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
    row = RowIndex[node];
    col = ColIndex[node];

    if (FlowLengthCode[row][col] == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }


    chi_vec[n_index] = dx*(pow( (A_0/ ( Discharge.get_data_element(row, col) ) ),m_over_n))
                          + chi_vec[IndexOfReceiverInUplsopePList];

    //  cout << "node: " << upslope_pixel_list[n_index] << " receiver: " << receiver_node
    //       << " SIndexReciever: " << IndexOfReceiverInUplsopePList
    //       << " and checked: " << upslope_pixel_list[IndexOfReceiverInUplsopePList]
    //       << " and chi: " << chi_vec[n_index] << endl;

  }
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is called from the the get_upslope_chi that only has an integer
// it returns the acutal chi values in a map where the key is the node and
// the value is chi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,float> LSDFlowInfo::get_upslope_chi_return_map(vector<int>& upslope_pixel_list,
                                           float m_over_n, float A_0, int minimum_pixels)
{

  map<int,float> map_of_chi;


  int receiver_node;
  //int IndexOfReceiverInUplsopePList;
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  float pixel_area = DataResolution*DataResolution;
  int node,row,col;
  // get the number of nodes upslope
  int n_nodes_upslope = upslope_pixel_list.size();
  //cout << "n_nodes upslope are: " << n_nodes_upslope << endl;
  //vector<float> chi_vec(n_nodes_upslope,0.0);

  if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
  {
    cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
    exit(EXIT_FAILURE);
  }

  //int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];
  map_of_chi[ upslope_pixel_list[0] ] = 0.0;


  for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
  {
    // We only get the data if the number of contributing pixels is greater than
    // the minimum pixels
    node = upslope_pixel_list[n_index];
    if (NContributingNodes[node] >= minimum_pixels)
    {
      receiver_node = ReceiverVector[ node ];
      //IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
      row = RowIndex[node];
      col = ColIndex[node];

      if (FlowLengthCode[row][col] == 2)
      {
        dx = diag_length;
      }
      else
      {
        dx = DataResolution;
      }
      map_of_chi[ node ] = dx*(pow( (A_0/ (float(NContributingNodes[node])*pixel_area) ),m_over_n))
                          + map_of_chi[ receiver_node ];
    }
  }
  return map_of_chi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is called from the the get_upslope_chi that only has an integer
// it returns the actual chi values in a map where the key is the node and
// the value is chi
// same as above but uses a discharge raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,float> LSDFlowInfo::get_upslope_chi_return_map(vector<int>& upslope_pixel_list,
                                           float m_over_n, float A_0, int minimum_pixels,
                                           LSDRaster& Discharge)
{

  map<int,float> map_of_chi;


  int receiver_node;
  //int IndexOfReceiverInUplsopePList;
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  int node,row,col;
  // get the number of nodes upslope
  int n_nodes_upslope = upslope_pixel_list.size();
  //vector<float> chi_vec(n_nodes_upslope,0.0);

  if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
  {
    cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
    exit(EXIT_FAILURE);
  }

  //int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];
  map_of_chi[ upslope_pixel_list[0] ] = 0.0;



  for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
  {
    // We only get the data if the number of contributing pixels is greater than
    // the minimum pixels
    node = upslope_pixel_list[n_index];
    if (NContributingNodes[node] >= minimum_pixels)
    {
      receiver_node = ReceiverVector[ node ];
      //IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
      row = RowIndex[node];
      col = ColIndex[node];

      if (FlowLengthCode[row][col] == 2)
      {
        dx = diag_length;
      }
      else
      {
        dx = DataResolution;
      }
      map_of_chi[node] = dx*(pow( (A_0/ ( Discharge.get_data_element(row, col) ) ),m_over_n))
                          + map_of_chi[ receiver_node ];
    }

  }
  return map_of_chi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-









//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function calculates chi upslope of a given starting node
// It returns a map with the node index as the key and the chi value as
// the value
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,float> LSDFlowInfo::get_upslope_chi_from_single_starting_node(int starting_node, float m_over_n, float A_0, int minimum_pixels)
{
  // get the pixel list
  vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);

  //cout << "Number of upslope nodes is: " << upslope_pixel_list.size() << endl;

  // Now get the upslope chi
  map<int,float> upslope_chi_map = get_upslope_chi_return_map(upslope_pixel_list,
                                                         m_over_n, A_0, minimum_pixels);

  return upslope_chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function calculates chi upslope of a given starting node
// It returns a map with the node index as the key and the chi value as
// the value
// Same as above but uses discharge
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,float> LSDFlowInfo::get_upslope_chi_from_single_starting_node(int starting_node,
                                 float m_over_n, float A_0, int minimum_pixels, LSDRaster& Discharge)
{
  // get the pixel list
  vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);

  // Now get the upslope chi
  map<int,float> upslope_chi_map = get_upslope_chi_return_map(upslope_pixel_list,
                                                         m_over_n, A_0, minimum_pixels, Discharge);

  return upslope_chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of starting nodes and calculates chi
// it assumes each chi value has the same base level.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_multiple_starting_nodes(vector<int>& starting_nodes,
                              float m_over_n, float A_0, float area_threshold)
{
  // some variables for writing to the raster
  int curr_row;
  int curr_col;
  int curr_node;

  float PixelArea = DataResolution*DataResolution;
  float DrainArea;

  // an array to hold the chi values
  Array2D<float> new_chi(NRows,NCols,NoDataValue);

  int n_starting_nodes = int(starting_nodes.size());

  for(int sn = 0; sn<n_starting_nodes; sn++)
  {
    // first check to see if this node has already been visited. If it has
    // all upslope nodes have also been visited so there is no point continuing
    // with this node
    retrieve_current_row_and_col(starting_nodes[sn],curr_row,curr_col);
    if(new_chi[curr_row][curr_col] == NoDataValue)
    {
      vector<float> us_chi = get_upslope_chi(starting_nodes[sn], m_over_n, A_0);
      vector<int> upslope_pixel_list = get_upslope_nodes(starting_nodes[sn]);

      int n_chi_nodes = int(us_chi.size());
      for (int cn = 0; cn<n_chi_nodes; cn++)
      {
        // get the current row and column
        curr_node =  upslope_pixel_list[cn];
        retrieve_current_row_and_col(curr_node,curr_row,curr_col);

        // check to see if the drainage area is greater than the threshold
        // if so, calcualte chi
        DrainArea = PixelArea*NContributingNodes[curr_node];
        if(DrainArea > area_threshold)
        {
          new_chi[curr_row][curr_col]= us_chi[cn];
        }
      }
   }
  }

  LSDRaster chi_map(NRows, NCols, XMinimum, YMinimum,
                    DataResolution, NoDataValue, new_chi,GeoReferencingStrings);
  return chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of starting nodes and calculates chi
// it assumes each chi value has the same base level.
// Same as above but calculates using discharge
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_multiple_starting_nodes(vector<int>& starting_nodes,
                              float m_over_n, float A_0, float area_threshold,
                              LSDRaster& Discharge )
{
  // some variables for writing to the raster
  int curr_row;
  int curr_col;
  int curr_node;

  float PixelArea = DataResolution*DataResolution;
  float DrainArea;

  // an array to hold the chi values
  Array2D<float> new_chi(NRows,NCols,NoDataValue);

  int n_starting_nodes = int(starting_nodes.size());

  for(int sn = 0; sn<n_starting_nodes; sn++)
  {
    // first check to see if this node has already been visited. If it has
    // all upslope nodes have also been visited so there is no point continuing
    // with this node
    retrieve_current_row_and_col(starting_nodes[sn],curr_row,curr_col);
    if(new_chi[curr_row][curr_col] == NoDataValue)
    {
      vector<float> us_chi = get_upslope_chi(starting_nodes[sn], m_over_n, A_0,Discharge);
      vector<int> upslope_pixel_list = get_upslope_nodes(starting_nodes[sn]);

      int n_chi_nodes = int(us_chi.size());
      for (int cn = 0; cn<n_chi_nodes; cn++)
      {
        // get the current row and column
        curr_node =  upslope_pixel_list[cn];
        retrieve_current_row_and_col(curr_node,curr_row,curr_col);

        // check to see if the drainage area is greater than the threshold
        // if so, calcualte chi
        DrainArea = PixelArea*NContributingNodes[curr_node];
        if(DrainArea > area_threshold)
        {
          new_chi[curr_row][curr_col]= us_chi[cn];
        }
      }
   }
  }

  LSDRaster chi_map(NRows, NCols, XMinimum, YMinimum,
                    DataResolution, NoDataValue, new_chi,GeoReferencingStrings);
  return chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of starting nodes and calculates chi
// it assumes each chi value has the same base level.
// Same as above but calculates using custom drainage_area/discharge, modified to ignore the threshold
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_multiple_starting_nodes_custom(vector<int>& starting_nodes,
                              float m_over_n, float A_0, LSDRaster& Discharge )
{
  // some variables for writing to the raster
  int curr_row;
  int curr_col;
  int curr_node;
  // an array to hold the chi values
  Array2D<float> new_chi(NRows,NCols,NoDataValue);

  int n_starting_nodes = int(starting_nodes.size());

  for(int sn = 0; sn<n_starting_nodes; sn++)
  {
    // first check to see if this node has already been visited. If it has
    // all upslope nodes have also been visited so there is no point continuing
    // with this node
    retrieve_current_row_and_col(starting_nodes[sn],curr_row,curr_col);
    if(new_chi[curr_row][curr_col] == NoDataValue)
    {
      vector<float> us_chi = get_upslope_chi(starting_nodes[sn], m_over_n, A_0,Discharge);
      vector<int> upslope_pixel_list = get_upslope_nodes(starting_nodes[sn]);

      int n_chi_nodes = int(us_chi.size());
      for (int cn = 0; cn<n_chi_nodes; cn++)
      {
        // get the current row and column
        curr_node =  upslope_pixel_list[cn];
        retrieve_current_row_and_col(curr_node,curr_row,curr_col);
        new_chi[curr_row][curr_col]= us_chi[cn];
      }
    }
  }

  LSDRaster chi_map(NRows, NCols, XMinimum, YMinimum,
                    DataResolution, NoDataValue, new_chi,GeoReferencingStrings);
  return chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function assumes all base level nodes are at the same base level
// and calculates chi for them. Essentially it covers the entire map in
// chi values.
// This function is probably most appropriate for looking at numerical
// model results
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_all_baselevel_nodes(float m_over_n, float A_0,
                  float area_threshold)
{
  LSDRaster all_chi = get_upslope_chi_from_multiple_starting_nodes(BaseLevelNodeList,
                                      m_over_n, A_0, area_threshold);
  return all_chi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function assumes all base level nodes are at the same base level
// and calculates chi for them. Essentially it covers the entire map in
// chi values.
// This function is probably most appropriate for looking at numerical
// model results
// same as above but calculates with a discharge
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_all_baselevel_nodes(float m_over_n, float A_0,
                  float area_threshold, LSDRaster& Discharge)
{
  LSDRaster all_chi = get_upslope_chi_from_multiple_starting_nodes(BaseLevelNodeList,
                                m_over_n, A_0, area_threshold, Discharge);
  return all_chi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// distance from outlet function
// this is overloaded.
// if it isn't provided any argument, it calculates the distance from outlet
// of all the base level nodes
// if it is given a node index number or a row and column, then
// the distance from outlet includes all the distances upstream of that node
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::distance_from_outlet()
{
  // initialize the array2d that will become the LSDRaster
  float ndv = float(NoDataValue);
  Array2D<float> flow_distance(NRows,NCols,ndv);

  // initialize the 1/root(2)
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;

  int row,col,bl_row,bl_col,receive_row,receive_col;

  int start_node = 0;
  int end_node;
  int nodes_in_bl_tree;
  int baselevel_node;
  // loop through the base level node list
  int n_base_level_nodes = BaseLevelNodeList.size();
  for(int bl = 0; bl<n_base_level_nodes; bl++)
  {
    baselevel_node = BaseLevelNodeList[bl];

    bl_row = RowIndex[baselevel_node];
    bl_col = ColIndex[baselevel_node];
    // get the number of nodes upslope and including this node
    nodes_in_bl_tree = NContributingNodes[baselevel_node];
    //cout << "LINE 938, FlowInfo, base level: " << bl << " with " << nodes_in_bl_tree << " nodes upstream" << endl;

    end_node = start_node+nodes_in_bl_tree;

    // set the distance of the outlet to zero
    flow_distance[bl_row][bl_col] = 0;

    // now loop through stack
    for(int s_node = start_node; s_node < end_node; s_node++)
    {
      //cout << "Line 953 flow info, s_node is: " << s_node << endl;

      //cout << SVector.size() << " " << ReceiverVector.size() << " " << RowIndex.size() << " " << ColIndex.size() << endl;
      row = RowIndex[ SVector[ s_node]  ];
      col = ColIndex[ SVector[ s_node]  ];
      //cout << "got rows and columns " << row << " " << col << endl;
      receive_row = RowIndex[ ReceiverVector[SVector[s_node] ]];
      receive_col = ColIndex[ ReceiverVector[SVector[s_node] ]];
      //cout <<  "get receive " << receive_row << " " << receive_col << endl;

      if ( FlowLengthCode[row][col] == 1)
      {
        flow_distance[row][col] = flow_distance[receive_row][receive_col]+DataResolution;
      }
      else if ( FlowLengthCode[row][col] == 2 )
      {
        flow_distance[row][col] = flow_distance[receive_row][receive_col]
                                  + diag_length;
      }
      //cout << "Flow distance: " << flow_distance << endl;
    }
    start_node = end_node;
  }
  //cout << "LINE 971 FlowInfo Flow distance complete, flow_distance is: " << endl;
  //cout << flow_distance << endl;
  LSDRaster FlowLength(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,
                       flow_distance,GeoReferencingStrings);
  return FlowLength;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function get the d8 slope. It points downslope from each node
//
// SMM 22/09/2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::calculate_d8_slope(LSDRaster& Elevation)
{
  float ndv = float(NoDataValue);
  Array2D<float> d8_slope(NRows,NCols,ndv);

  // these save a bit of computational expense.
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;

  int this_row;
  int this_col;
  int r_row;
  int r_col;
  int r_node;
  float dx;

  for (int node = 0; node<NDataNodes; node++)
    {
      // get the row and column
      retrieve_current_row_and_col(node,this_row,this_col);

      // get the distance between nodes. Depends on flow direction
      switch (retrieve_flow_length_code_of_node(node))
  {
  case 0:
    dx = -99;
    break;
  case 1:
    dx = DataResolution;
    break;
  case 2:
    dx = dx_root2;
    break;
  default:
    dx = -99;
    break;
  }

      // get the receiver information
      retrieve_receiver_information(node,r_node, r_row, r_col);

      // now calculate the slope
      if (r_node == node)
  {
    d8_slope[this_row][this_col] = 0;
  }
      else
  {
    d8_slope[this_row][this_col] = (1/dx)*
      (Elevation.get_data_element(this_row,this_col)
       -Elevation.get_data_element(r_row,r_col));
  }

    }

  LSDRaster d8_slope_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,
          d8_slope,GeoReferencingStrings);
  return d8_slope_raster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//
// this finds the node that is farthest upstream from a given node
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::find_farthest_upslope_node(int node, LSDRaster& DistFromOutlet)
{
  // set the farthest node to the current node; if the node has no contributing pixels
  // the function will just return itself
  int farthest_upslope_node = node;

  // first get the nodes that are upslope
  vector<int> upslope_node_list = get_upslope_nodes(node);

  int row, col;
  float this_flow_distance;

  // now loop through these, looking for the farthest upstream node
  float farthest = 0.0;
  int n_upslope_nodes = upslope_node_list.size();
  for (int i = 0; i<n_upslope_nodes; i++)
  {
    // get the row and col of upslope nodes
    row = RowIndex[ upslope_node_list[i] ];
    col = ColIndex[ upslope_node_list[i] ];

    // get the flow distance
    this_flow_distance = DistFromOutlet.get_data_element(row, col);

    // test if it is the farthest
    if (this_flow_distance > farthest)
    {
      farthest = this_flow_distance;
      farthest_upslope_node =  upslope_node_list[i];
    }
  }

  return farthest_upslope_node;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a list of nodes and sorts them according to a sorting
// raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::sort_node_list_based_on_raster(vector<int> node_vec, LSDRaster& SortingRaster)
{
  // First make a vector of the data from the sorting vector
  vector<float> VectorisedData;
  vector<float> SortedData;
  vector<int> SortedNodes;

  vector<size_t> index_map;

  int row,col;

  // loop through node vec, getting the value of the raster at each node
  int n_nodes = int(node_vec.size());
  for (int n = 0; n<n_nodes; n++)
  {
    retrieve_current_row_and_col(node_vec[n],row,col);

    // now get the data element
    VectorisedData.push_back(SortingRaster.get_data_element(row,col));
  }

  // now sort that data using the matlab sort
  matlab_float_sort(VectorisedData, SortedData, index_map);

  // now sort the nodes based on this sorting
  matlab_int_reorder(node_vec, index_map, SortedNodes);

  return SortedNodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a list of nodes and sorts them according to a sorting
// raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::sort_node_list_based_on_raster(vector<int> node_vec, LSDIndexRaster& SortingRaster)
{
  // First make a vector of the data from the sorting vector
  vector<int> VectorisedData;
  vector<int> SortedData;
  vector<int> SortedNodes;

  vector<size_t> index_map;

  int row,col;

  // loop through node vec, getting the value of the raster at each node
  int n_nodes = int(node_vec.size());
  for (int n = 0; n<n_nodes; n++)
  {
    retrieve_current_row_and_col(node_vec[n],row,col);

    // now get the data element
    VectorisedData.push_back(SortingRaster.get_data_element(row,col));
  }

  // now sort that data using the matlab sort
  matlab_int_sort(VectorisedData, SortedData, index_map);

  // now sort the nodes based on this sorting
  matlab_int_reorder(node_vec, index_map, SortedNodes);

  return SortedNodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// get node index of point from X and Y coordinates
// this is different from the above function in that it does not snap to the nearest channel
//
// FJC 11/02/14
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::get_node_index_of_coordinate_point(float X_coordinate, float Y_coordinate)
{
  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(ceil(Y_coordinate_shifted_origin/DataResolution)-0.5);


  // Get node of point
  int CurrentNode;
  if(col_point>=0 && col_point<NCols && row_point>=0 && row_point<NRows)
  {
    CurrentNode = retrieve_node_from_row_and_column(row_point, col_point);
  }
  else
  {
    CurrentNode = NoDataValue;
  }
  return CurrentNode;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get vector of nodeindices from a csv file with X and Y coordinates
// The csv file must have the format: ID, X, Y
//
// FJC 14/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::get_nodeindices_from_csv(string csv_filename, vector<int>& NIs, vector<float>& X_coords, vector<float>& Y_coords)
{
  ifstream csv_input(csv_filename.c_str());
  vector<int> temp_NIs;
  vector<float> temp_X_Coords;
  vector<float> temp_Y_Coords;
  //initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  // discard the first line
  getline(csv_input, line_from_file);

  // now loop through the rest of the lines
  while (getline(csv_input,line_from_file))
  {
    // reset the string vec
    this_string_vec = empty_string_vec;

    // create a stringstream
    stringstream ss(line_from_file);

    while (ss.good())
    {
      string substr;
      getline(ss,substr,',');

      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

      // add the string to the string vec
      this_string_vec.push_back( substr );
    }

    // for some reason our compiler can't deal with stof so converting to doubles
    double X_coordinate = atof(this_string_vec[1].c_str());
    double Y_coordinate = atof(this_string_vec[2].c_str());

    // Shift origin to that of dataset
    float X_coordinate_shifted_origin = X_coordinate - XMinimum - DataResolution*0.5;
    float Y_coordinate_shifted_origin = Y_coordinate - YMinimum - DataResolution*0.5;

    // Get row and column of point
    int col_point = int(X_coordinate_shifted_origin/DataResolution);
    int row_point = (NRows-1) - int(round(Y_coordinate_shifted_origin/DataResolution));

    // Get node of point
    int CurrentNode = NoDataValue;
    if(col_point>=0 && col_point<NCols && row_point>=0 && row_point<NRows)
    {
      CurrentNode = retrieve_node_from_row_and_column(row_point, col_point);
    }
    if (CurrentNode != NoDataValue)
    {
      temp_X_Coords.push_back(float(X_coordinate));
      temp_Y_Coords.push_back(float(Y_coordinate));
      temp_NIs.push_back(CurrentNode);
    }

  }

  //copy to output vectors
  NIs = temp_NIs;
  X_coords = temp_X_Coords;
  Y_coords = temp_Y_Coords;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// For a given node, find the nearest point in the raster that is not a NoDataValue and
// snap it to the node.  User must specify the search radius in pixels.
// FJC 08/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float LSDFlowInfo::snap_RasterData_to_Node(int NodeIndex, LSDRaster& InputRaster, int search_radius)
{
  float RasterValue_at_Node;
  Array2D<float> InputRasterData = InputRaster.get_RasterData();
  int i,j;
  retrieve_current_row_and_col(NodeIndex,i,j);
  if (InputRasterData[i][j] != NoDataValue)
  {
    // The node index already has a raster value!
    RasterValue_at_Node = InputRasterData[i][j];
  }
  else
  {

    vector<float> Data_in_window;
    vector<float> Dists_in_window;
    vector<float> Dists_in_window_sorted;
    vector<size_t> index_map;
    Data_in_window.reserve(4 * search_radius);
    Dists_in_window.reserve(4 * search_radius);
    Dists_in_window_sorted.reserve(4 * search_radius);
    index_map.reserve(4 * search_radius);

    //set up the bounding box
    int i_min = i - search_radius;
    int i_max = i + search_radius;
    int j_min = j - search_radius;
    int j_max = j + search_radius;

    //out of bounds checking
    if (i_min < 0){i_min = 0;}
    if (j_min < 0){j_min = 0;}
    if (i_max > (NRows - 1)){i_max = (NRows - 1);}
    if (j_max > (NCols - 1)){j_max = (NCols - 1);}

    // only iterate over the search area.
    for (int row = i_min; row < i_max; ++row){
      for (int col = j_min; col < j_max; ++col){

        if (InputRasterData[row][col] != NoDataValue)
        {
          //get the  raster data and distance at each point in the window
          Data_in_window.push_back(InputRasterData[row][col]);

          float Dist = distbetween(i,j,row,col);
          Dists_in_window.push_back(Dist);
        }
      }
    }
    if (int(Data_in_window.size()) != 0)
    {
      matlab_float_sort(Dists_in_window, Dists_in_window_sorted, index_map);

      // return the raster value at the smallest distance from the point
      RasterValue_at_Node = Data_in_window[index_map[0]];
    }
    else
    {
      // if no raster values in the window then return nodatavalue
      RasterValue_at_Node = NoDataValue;
    }
  }

  return RasterValue_at_Node;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses the flow accumulation pixels
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_index_threshold(LSDIndexRaster& FlowPixels, int threshold)
{
  vector<int> sources;
  int row,col;
  //int n_donors;
  int donor_row,donor_col;
  int thresh_switch;
  int donor_node;

  // drop down through the stack
  // if the node is greater than or equal to the threshold
  // check all the donors
  // if there are no donors, then the node is a source
  // if none of the donors are greater than the threshold, then it also is a source
  for (int node = 0; node<NDataNodes; node++)
    {
      row = RowIndex[node];
      col = ColIndex[node];

      // see if node is greater than threshold
      if(FlowPixels.get_data_element(row,col)>=threshold)
  {
    //cout << "node " << node << " is a potential source, it has a value of "
    //     << FlowPixels.get_data_element(row,col)
    //     << "and it has " << NDonorsVector[node] <<" donors " << endl;

    // if it doesn't have donors, it is a source
    if(NDonorsVector[node] == 0)
      {
        sources.push_back(node);
      }
    else
      {
        thresh_switch = 1;
        // figure out where the donor nodes are, and if
        // the donor node is greater than the threshold
        for(int dnode = 0; dnode<NDonorsVector[node]; dnode++)
    {
      donor_node = DonorStackVector[ DeltaVector[node]+dnode];
      donor_row = RowIndex[ donor_node ];
      donor_col = ColIndex[ donor_node ];

      // we don't float count base level nodes, which donate to themselves
      if (donor_node != node)
      {
        // if the donor node is greater than the threshold,
        // then this node is not a threhold
        if(FlowPixels.get_data_element(donor_row,donor_col)>=threshold)
        {
          thresh_switch = 0;
        }
      }

      //cout << "thresh_switch is: " << thresh_switch << endl;
    }
        // if all of the donors are below the threhold, this is a source
        if (thresh_switch == 1)
    {
      sources.push_back(node);
    }
      }
  }
    }
  return sources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses a threshold of drainage area * slope^2
//
//
// FJC 11/02/14
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_slope_area(LSDIndexRaster& FlowPixels, LSDRaster& Slope, int threshold)
{
  vector<int> sources;
  int row,col;
  //int n_donors;
  int donor_row,donor_col;
  int thresh_switch;
  int donor_node;

  // drop down through the stack
  // if the node is greater than or equal to the threshold
  // check all the donors
  // if there are no donors, then the node is a source
  // if none of the donors are greater than the threshold, then it also is a source
  for (int node = 0; node<NDataNodes; node++)
    {
      row = RowIndex[node];
      col = ColIndex[node];

      float area = FlowPixels.get_data_element(row,col);
      float slope = Slope.get_data_element(row,col);
      float SA_product = area * (slope*slope);
      // see if node is greater than threshold
      if(SA_product >= threshold)
  {
    //cout << "node " << node << " is a potential source, it has a value of "
    //     << SA_product
    //     << "and it has " << NDonorsVector[node] <<" donors " << endl;

    // if it doesn't have donors, it is a source
    if(NDonorsVector[node] == 0)
      {
        sources.push_back(node);
      }
    else
      {
        thresh_switch = 1;
        // figure out where the donor nodes are, and if
        // the donor node is greater than the threshold
        for(int dnode = 0; dnode<NDonorsVector[node]; dnode++)
    {
      donor_node = DonorStackVector[ DeltaVector[node]+dnode];
      donor_row = RowIndex[ donor_node ];
      donor_col = ColIndex[ donor_node ];

      // we don't float count base level nodes, which donate to themselves
      if (donor_node != node)
        {
          // if the donor node is greater than the threshold,
          // then this node is not a threhold
          float area_donor = FlowPixels.get_data_element(donor_row,donor_col);
          float slope_donor = Slope.get_data_element(donor_row,donor_col);
          float SA_product_donor = area_donor * (slope_donor*slope_donor);
          if(SA_product_donor >= threshold)
      {
        thresh_switch = 0;
      }
        }

      //cout << "thresh_switch is: " << thresh_switch << endl;
    }
        // if all of the donors are below the threhold, this is a source
        if (thresh_switch == 1)
    {
      sources.push_back(node);
    }
      }
  }
    }
  return sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses the X and Y coordinates of mapped channel heads
//
//
// FJC 17/02/14
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_from_mapped_channel_heads(vector<float>& X_coords, vector<float>& Y_coords)
{
  vector<int> SourceNodes;
  cout << "N of channel heads: " << X_coords.size() << endl;
  for (unsigned int i = 0; i < X_coords.size(); i++)
    {
      int NI = get_node_index_of_coordinate_point(X_coords[i], Y_coords[i]);
      if (NI != NoDataValue)
  {
    SourceNodes.push_back(NI);
  }
    }

  return SourceNodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Perform a downslope trace using D8 from a given point source (i,j).
//Overwrites input parameters to return a raster of the path, the length of the
//trace and the final pixel coordinates of the trace.
// SWDG 20/1/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::D8_Trace(int i, int j, LSDIndexRaster StreamNetwork, float& length, int& receiver_row, int& receiver_col, Array2D<int>& Path){

  float root_2 = 1.4142135623;

  Array2D<int> stnet = StreamNetwork.get_RasterData();

  length = 0;

  int node;

  int receiver_node = retrieve_node_from_row_and_column(i, j);
  receiver_row = i;
  receiver_col = j;

  Path[receiver_row][receiver_col] = 1;

  while (StreamNetwork.get_data_element(receiver_row, receiver_col) == NoDataValue){  // need to do edge checking

    retrieve_receiver_information(receiver_node, node, receiver_row, receiver_col);

    Path[receiver_row][receiver_col] = 1;

    //update length
    if (retrieve_flow_length_code_of_node(receiver_node) == 1){ length += DataResolution; }
    else if (retrieve_flow_length_code_of_node(receiver_node) == 2){ length += (DataResolution * root_2); }

    receiver_node = node;

  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets a flow path. Like D8 trace, but just get a vector of nodeindices rather than anything
// more fancy
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_flow_path(float x_loc, float y_loc)
{
  vector<int> node_list;

  // check to see if the node is in the raster
  bool is_in_raster = check_if_point_is_in_raster(x_loc,y_loc);

  int receiver_node, current_node;

  // check that the point is not nodata and set is_in_raster to false if it is
  if (is_in_raster)
  {
    current_node = get_node_index_of_coordinate_point(x_loc,y_loc);
    node_list.push_back(current_node);

    retrieve_receiver_information(current_node,receiver_node);

    while (receiver_node != current_node)
    {  // need to do edge checking

      current_node = receiver_node;
      retrieve_receiver_information(current_node,receiver_node);

      node_list.push_back(current_node);
    }
  }
  return node_list;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets a flow path. Like D8 trace, but just get a vector of nodeindices rather than anything
// more fancy
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_flow_path(int ni)
{
  vector<int> node_list;

  int receiver_node, current_node;

  current_node = ni;
  node_list.push_back(current_node);

  retrieve_receiver_information(current_node,receiver_node);

  while (receiver_node != current_node)
  {  // need to do edge checking

    current_node = receiver_node;
    retrieve_receiver_information(current_node,receiver_node);

    node_list.push_back(current_node);
  }
  return node_list;
}


float LSDFlowInfo::calculate_nci(int start_node, int end_node, LSDRaster& elevation, LSDRaster& distance_from_outlet)
{
  vector<float> elevations;
  vector<float> distances;

  int this_node = start_node;
  int receiver_node;
  int this_row, this_col;
  float this_elevation; 
  float this_distance;
  
  this_row = RowIndex[this_node];
  this_col = ColIndex[this_node];
  this_elevation = elevation.get_data_element(this_row,this_col);
  this_distance = distance_from_outlet.get_data_element(this_row,this_col);

  elevations.push_back(this_elevation);
  distances.push_back(this_distance);

  bool not_found_baselevel = true;
  do
  {
    retrieve_receiver_information(this_node,receiver_node);
    if(this_node == receiver_node)
    {
      not_found_baselevel = false;
    }
    else
    {
      this_node = receiver_node;

      this_row = RowIndex[this_node];
      this_col = ColIndex[this_node]; 
      this_elevation = elevation.get_data_element(this_row,this_col);
      this_distance = distance_from_outlet.get_data_element(this_row,this_col);
      elevations.push_back(this_elevation);
      distances.push_back(this_distance);
    }
  } while ( this_node != end_node && not_found_baselevel);

  int n_elevations = int(elevations.size());
  
  float min_elevation = elevations[n_elevations-1];
  float max_elevation = elevations[0];
  float min_distance = distances[n_elevations-1];
  float max_distance = distances[0];
  float slope = (max_elevation-min_elevation)/(max_distance-min_distance);
  float intercept = max_elevation-slope*max_distance;

  vector<float> norm_delta_z;
  float dz;
  for(int i = 0; i< n_elevations; i++)
  {
    dz =  elevations[i] - (slope*distances[i]+intercept);
    norm_delta_z.push_back( dz/(max_elevation - min_elevation)  );
  } 
  float NDV = -9999;
  float nci = get_median(norm_delta_z,NDV);
  return nci;
}


void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int> rows,vector<int> cols,string path, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  double x_loc,y_loc;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_rows = int(rows.size());

  if (rows.size() != cols.size())
  {
    cout << "You are trying to print a csv from a row and column list but the two vectors are not the same size" << endl;
    cout << "Fatal error, exiting. " << endl;
    exit(0);
  }

  // open the data file
  ofstream chan_out;
  string chan_fname = path+filename+"_nodes.csv";
  chan_out.open(chan_fname.c_str());
  chan_out << "nodeindex,row,column,latitude,longitude,x,y" << endl;

  for (int n = 0; n<n_rows; n++)
  {
    this_node = get_NodeIndex_from_row_col(rows[n],cols[n]);
    get_lat_and_long_locations(rows[n],cols[n], latitude, longitude, Converter);
    get_x_and_y_locations(rows[n],cols[n], x_loc, y_loc);

    chan_out << this_node << ","
                   << rows[n] << ","
                   << cols[n] << ",";
    chan_out.precision(9);
    chan_out << latitude << ","
             << longitude << ",";
    chan_out.precision(9);
    chan_out << x_loc << "," << y_loc << endl;
  }

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Prints a node list
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int> node_list,string path, string filename)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  double x_loc,y_loc;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = int(node_list.size());

  // open the data file
  ofstream chan_out;
  string chan_fname = path+filename+"_nodes.csv";
  chan_out.open(chan_fname.c_str());
  chan_out << "nodeindex,row,column,latitude,longitude,x,y" << endl;
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }

  for (int n = 0; n<n_nodes; n++)
  {
    this_node = node_list[n];
    retrieve_current_row_and_col(this_node,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);
    get_x_and_y_locations(row, col, x_loc, y_loc);

    chan_out << this_node << ","
                   << row << ","
                   << col << ",";
    chan_out.precision(9);
    chan_out << latitude << ","
             << longitude << ",";
    chan_out.precision(9);
    chan_out << x_loc << "," << y_loc << endl;
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Prints a node list, with some extra bells and whistles
// Basically you get extra data members so that you can print a channel profile
// in excel (so undergraduate friendly)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_with_latlong(vector<int> node_list,string path, string filename,
                                                                       LSDRaster& Elevation, LSDRaster& FlowDistance,
                                                                       LSDRaster& drainage_area)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  double x_loc,y_loc;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = int(node_list.size());

  // open the data file
  ofstream chan_out;
  string chan_fname = path+filename+"_nodes.csv";
  chan_out.open(chan_fname.c_str());
  chan_out << "id,row,column,latitude,longitude,x,y, flow distance(m), drainage area (m^2), elevation(m)" << endl;
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }

  for (int n = 0; n<n_nodes; n++)
  {
    this_node = node_list[n];
    retrieve_current_row_and_col(this_node,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);
    get_x_and_y_locations(row, col, x_loc, y_loc);

    chan_out << this_node << ","
                   << row << ","
                   << col << ",";
    chan_out.precision(9);
    chan_out << latitude << ","
             << longitude << ",";
    chan_out.precision(9);
    chan_out << x_loc << "," << y_loc << ",";
    chan_out.precision(6);
    chan_out << FlowDistance.get_data_element(row,col) << "," << drainage_area.get_data_element(row,col) << "," << Elevation.get_data_element(row,col) << endl;
  }
}



void LSDFlowInfo::calculate_bearings_from_nodelist(vector<int> nodeindex_vec, int node_spacing,
                                                 bool print_bearings, string outfile_name,
                                                 vector<int>& bearing_nodes, vector<float>& bearing_vec)
{

  ofstream bearings;
  if(print_bearings)
  {
    bearings.open(outfile_name.c_str());
    bearings << "latitude,longitude,bearing,bearing_degree" << endl;
    bearings.precision(9);
  }

  vector<int> nodes;
  vector<float> bearings_vec;

  int top_node;
  int middle_node;
  int bottom_node;

  int t_row,t_col,m_row,m_col,b_row,b_col;
  float t_x_loc,t_y_loc,m_x_loc,m_y_loc,b_x_loc,b_y_loc;

  float bearing_u, bearing_d, bearing;

  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  int n_nodes = int( nodeindex_vec.size());

  top_node = 0;
  middle_node = node_spacing;
  bottom_node = node_spacing*2;

  while (bottom_node < n_nodes)
  {
    retrieve_current_row_and_col(nodeindex_vec[top_node],t_row,t_col);
    get_x_and_y_locations(t_row, t_col, t_x_loc, t_y_loc);

    retrieve_current_row_and_col(nodeindex_vec[middle_node],m_row,m_col);
    get_x_and_y_locations(m_row, m_col, m_x_loc, m_y_loc);

    retrieve_current_row_and_col(nodeindex_vec[bottom_node],b_row,b_col);
    get_x_and_y_locations(b_row, b_col, b_x_loc, b_y_loc);

    bearing_u = clockwise_angle_between_vector_and_north( t_x_loc,  t_y_loc, m_x_loc,  m_y_loc);
    bearing_d = clockwise_angle_between_vector_and_north( m_x_loc,  m_y_loc, b_x_loc,  b_y_loc);

    bearing = 0.5*(bearing_u-bearing_d)+bearing_d;

    nodes.push_back(nodeindex_vec[middle_node]);
    bearings_vec.push_back(bearing);

    if(print_bearings)
    {
      get_lat_and_long_locations(m_row, m_col, latitude, longitude, Converter);
      bearings << latitude << "," << longitude << "," << bearing << "," << deg(bearing) << endl;
    }

    top_node++;
    middle_node++;
    bottom_node++;

  }

  if(print_bearings)
  {
    bearings.close();
  }

  bearing_nodes = nodes;
  bearing_vec = bearings_vec;


}





// Move the location of the channel head downslope by a user defined distance.
// Returns A vector of node indexes pointing to the moved heads.
// SWDG 27/11/15

void LSDFlowInfo::MoveChannelHeadDown(vector<int> Sources, float MoveDist, vector<int>& DownslopeSources, vector<int>& FinalHeads){

  float root_2 = 1.4142135623;

  float length;

  int receiver_row;
  int receiver_col;
  int node;

  //for loop goes here to iterate over the Sources vector
  for(int q=0; q<int(Sources.size());++q){
    length = 0;

    int receiver_node = Sources[q];

    while (length < MoveDist){

      retrieve_receiver_information(receiver_node, node, receiver_row, receiver_col);

      //update length
      if (retrieve_flow_length_code_of_node(receiver_node) == 1){ length += DataResolution; }
      else if (retrieve_flow_length_code_of_node(receiver_node) == 2){ length += (DataResolution * root_2); }
      else if (retrieve_flow_length_code_of_node(receiver_node) == 0){break;}

      receiver_node = node;


    }

    DownslopeSources.push_back(receiver_node);
    FinalHeads.push_back(Sources[q]);

  }

  //end of for loop

}

// Move the location of the channel head upslope by a user defined distance.
// Returns A vector of node indexes pointing to the moved heads.
// SWDG 27/11/15
void LSDFlowInfo::MoveChannelHeadUp(vector<int> Sources, float MoveDist, LSDRaster DEM, vector<int>& UpslopeSources, vector<int>& FinalHeads){

  float root_2 = 1.4142135623;

  float length;

  Array2D<float> Elevation = DEM.get_RasterData();

  int new_node;

  int new_i;
  int new_j;

  //for loop goes here to iterate over the Sources vector
  for(int q=0; q<int(Sources.size());++q){
    length = 0;

    int i;
    int j;

    retrieve_current_row_and_col(Sources[q], i, j);

    //test for channel heads at edges
    if (i == 0 || i == NRows - 1 || j == 0 || j == NCols - 1){
      cout << "Hit an edge, skipping" << endl;
    }

    else{

    while (length < MoveDist){

      float currentElevation;
      int Direction = 1; //1 is cardinal 2 is diagonal 1 is default
      new_j = 1;    // default
      new_i = 1;    // default

      //find the neighbour with the maximum Elevation

      currentElevation = Elevation[i][j];

      //top left
      if(currentElevation < Elevation[i-1][j-1]){
        currentElevation = Elevation[i-1][j-1];
        new_i = i-1;
        new_j = j-1;
        Direction = 2;
      }
      //top
      if(currentElevation < Elevation[i][j-1]){
        currentElevation = Elevation[i][j-1];
        new_i = i;
        new_j = j-1;
        Direction = 1;
      }
      //top right
      if(currentElevation < Elevation[i+1][j-1]){
        currentElevation = Elevation[i+1][j-1];
        new_i = i+1;
        new_j = j-1;
        Direction = 2;
      }
      //right
      if(currentElevation < Elevation[i+1][j]){
        currentElevation = Elevation[i+1][j];
        new_i = i+1;
        new_j = j;
        Direction = 1;
      }
      //botom right
      if(currentElevation < Elevation[i+1][j+1]){
        currentElevation = Elevation[i+1][j+1];
        new_i = i+1;
        new_j = j+1;
        Direction = 2;
      }
      //bottom
      if(currentElevation < Elevation[i][j+1]){
        currentElevation = Elevation[i][j+1];
        new_i = i;
        new_j = j+1;
        Direction = 1;
      }
      //bottom left
      if(currentElevation < Elevation[i-1][j+1]){
        currentElevation = Elevation[i-1][j+1];
        new_i = i-1;
        new_j = j+1;
        Direction = 2;
      }
      //left
      if(currentElevation < Elevation[i-1][j]){
        currentElevation = Elevation[i-1][j];
        new_i = i-1;
        new_j = j;
        Direction = 1;
      }
      //test that we have not hit the ridge
      //this will exit the loop and add the final visited
      //node to the upper soureces vector
      if (currentElevation == Elevation[i][j]){
        cout << "Warning, unable to move channel head " << Sources[q] << " up by user defined distance" << endl;
        break;
      }

      //update length
      if (Direction == 1)
      {
        length += DataResolution;
      }
      else if (Direction == 2)
      {
        length += (DataResolution * root_2);
      }

      i = new_i;
      j = new_j;

    }
}
    new_node = retrieve_node_from_row_and_column(i,j);

    UpslopeSources.push_back(new_node);
    FinalHeads.push_back(Sources[q]);

  } //end of for loop

}

void LSDFlowInfo::HilltopFlowRoutingOriginal(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope, LSDRaster Aspect, LSDIndexRaster StreamNetwork)
{
  //Declare parameters
  int i,j,a,b;
  //double X,Y;
  float dem_res = DataResolution;
  //float mean_slope, relief;
  int ndv = NoDataValue;
  double slope_total, length, d;
  int flag, count;
  double PI = 3.14159265;
  double degs, degs_old, degs_new, theta;
  double s_local;
  double s_edge;
  double xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;
  double xmin = XMinimum;
  double ymin = YMinimum;
  double ymax = ymin + NRows*dem_res;
  //double xmax = xmin + NCols*dem_res;

  //Declare Arrays
  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols);
  Array2D<float> blank(NRows,NCols,NoDataValue);

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  //for creating file names for hillslope profiles
  string file_part_1 = "prof_";
  string file_part_2;
  string file_part_3 = ".txt";
  string filename;
  char filename_c[256];

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*dem_res - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = xmin + j*dem_res + 0.5;
  }

  int ht_count = 0;

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    for (j=1; j<NCols-1; ++j) {
      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != ndv) {

  //reset slope, length, hillslope flag and pixel counter
  slope_total = 0;
  length = 0;
  flag = true;
  count = 1;

  //copt blank raster to map hillslope trace
  path = blank.copy();

  //update hilltop counter
  ++ht_count;

  //get aspect in radians
  degs = aspect[i][j];
  theta = (M_PI/180.)*((-1*degs)+90.);

  //setup indices
  a = i;
  b = j;
  path[a][b] = 1;

  //add first pixel to easting northing vectors
  east_vec[0] = easting[b];
  north_vec[0] = northing[a];

  //get local slope
  s_local = slope[a][b];

  //test direction, calculate outlet coordinates and update indicies
  // easterly, dir == 1
  if (degs >= 45 && degs < 135) {
    //find location where trace exits the cell and distance
    xo = 1., yo = (1.+tan(theta))/2.;
    d = abs(1./(2.*cos(theta)));
    //transmit to next cell over
    xi = 0., yi = yo;
    dir = 1;
    //add to vector
    east_vec[count] = easting[b] + 0.5*dem_res;
    north_vec[count] = northing[a] + yo - 0.5*dem_res;
    //increment indices
    ++b;
    //check we're not right in the corner!
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //southerly
  else if (degs >= 135 && degs < 225) {
    //find location where trace exits the cell and distance
    xo = (1-(1/tan(theta)))/2, yo = 0;
    d = abs(1/(2*cos((PI/2)-theta)));
    //transmit to next cell over
    xi = xo, yi = 1;
    dir = 2;
    //add to vector
    east_vec[count] = easting[b] + xo - 0.5*dem_res;
    north_vec[count] = northing[a] - 0.5*dem_res;
    //increment indices
    ++a;
    //check we're not right in the corner!
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  // westerly
  else if (degs >= 225 && degs < 315) {
    //find location where trace exits the cell and distance
    xo = 0, yo = (1-tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    //transmit to next cell over
    xi = 1,  yi = yo;
    dir = 3;
    //add to vector
    east_vec[count] = easting[b] -0.5*dem_res;
    north_vec[count] = northing[a] + yo - 0.5*dem_res;
    //increment indices
    --b;
    //check we're not right in the corner!
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //northerly
  else if (degs >= 315 || degs < 45) {
    //find location where trace exits the cell and distance
    xo = (1+(1/tan(theta)))/2, yo = 1;
    d = abs(1/(2*cos((PI/2) - theta)));
    //transmit to next cell over
    xi = xo, yi = 0;
    dir = 4;
    //add to vector
    east_vec[count] = easting[b] + xo - 0.5*dem_res;
    north_vec[count] = northing[a] + 0.5*dem_res;
    //increment indices
    --a;
    //check we're not right in the corner!
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  else {
    cout << "FATAL ERROR, Kinematic routing algorithm encountered null aspect value" << endl;
    exit(EXIT_FAILURE);
  }

  //collect slopes and totals weighted by path length
  //slope_total += s_local*d;
  //length += d;

  s_local = slope[a][b];


  //continue trace until a stream node is encountered
  while (flag == true) {

    path[a][b] = 1;

    degs_old = degs;
    degs_new = aspect[a][b];
    theta = (M_PI/180.)*((-1*degs_new)+90.);
    ++ count;

    //Test for perimeter flow paths
    if ((dir == 1 && degs_new > 0 && degs_new < 180)
        ||  (dir == 2 && degs_new > 90 && degs_new < 270)
        ||  (dir == 3 && degs_new > 180 && degs_new < 360)
        ||  ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

      //DO NORMAL FLOW PATH
      //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
      temp_yo1 = yi + (1-xi)*tan(theta);     // xo = 1
      temp_xo1 = xi + (1-yi)*(1/tan(theta));   // yo = 1
      temp_yo2 = yi - xi*tan(theta);      // xo = 0
      temp_xo2 = xi - yi*(1/tan(theta));    // yo = 0

      // can't outlet at same point as inlet
      if (dir == 1) temp_yo2 = -1;
      else if (dir == 2) temp_xo1 = -1;
      else if (dir == 3) temp_yo1 = -1;
      else if (dir == 4) temp_xo2 = -1;

      s_local = slope[a][b];

      if (temp_yo1 <= 1 && temp_yo1 > 0) {
        xo = 1, yo = temp_yo1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 0, yi = yo,
    dir = 1;
        east_vec[count] = easting[b] + 0.5*dem_res;
        north_vec[count] = northing[a] + yo - 0.5*dem_res;
        ++b;
        if (xi== 0 && yi == 0) yi = 0.00001;
        else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
      }
      else if (temp_xo2 <= 1 && temp_xo2 > 0) {
        xo = temp_xo2, yo = 0;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 1,
    dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*dem_res;
        north_vec[count] = northing[a] - 0.5*dem_res;
        ++a;
        if (xi== 0 && yi == 1) xi = 0.00001;
        else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
      }
      else if (temp_yo2 <= 1 && temp_yo2 > 0) {
        xo = 0, yo = temp_yo2;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 1, yi = yo,
    dir = 3;
        east_vec[count] = easting[b] -0.5*dem_res;
        north_vec[count] = northing[a] + yo - 0.5*dem_res;
        --b;
        if (xi== 1 && yi == 0) yi = 0.00001;
        else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
      }

      else if (temp_xo1 <= 1 && temp_xo1 > 0) {
        xo = temp_xo1, yo = 1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 0,
    dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*dem_res;
        north_vec[count] = northing[a] + 0.5*dem_res;
        --a;

        if (xi == 0 && yi == 0) xi = 0.00001;
        else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
      }
      slope_total += s_local*d;
    }

    else
    {

      // ROUTE ALONG EDGES
      if (dir  == 1)
      {
        if   (degs_old <= 90 || degs_new >= 270)
        {
          xo = 0.00001, yo = 1;
          s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo, yi = 1-yo;
          dir = 4;
          east_vec[count] = easting[b] + xo - 0.5*dem_res;
          north_vec[count] = northing[a] + 0.5*dem_res;
          --a;
        }
        else if (degs_old > 90 && degs_new < 270)
        {
          xo = 0.00001, yo = 0;
          s_edge = abs(s_local*sin((PI/2)-theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo, yi = 1-yo;
          dir = 2;
          east_vec[count] = easting[b] + xo - 0.5*dem_res;
          north_vec[count] = northing[a] - 0.5*dem_res;
          ++a;
        }
        else
        {
          cout << "Flow unable to route N or S" << endl;
          exit(EXIT_FAILURE);
        }
      }
      else if (dir == 2)
      {
        if (degs_old <= 180 && degs_new >= 0)
        {
          xo = 1, yo = 1-0.00001;
          s_edge = abs(s_local*sin((2/PI)-theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 1;
          east_vec[count] = easting[b] + 0.5*dem_res;
          north_vec[count] = northing[a] + yo - 0.5*dem_res;
          ++b;
        }
        else if (degs_old > 180 && degs_new < 360)
        {
          xo = 0, yo = 1-0.00001;
          s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 3;
          east_vec[count] = easting[b] -0.5*dem_res;
          north_vec[count] = northing[a] + yo - 0.5*dem_res;
          --b;
        }
        else
         {
          cout << "Flow unable to route E or W" << endl;
          exit(EXIT_FAILURE);
        }
      }
      else if (dir == 3)
      {
        if   (degs_old <= 270 && degs_new >= 90)
        {
          xo = 1-0.00001, yo = 0;
          s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo, yi = 1-yo;
          dir = 2;
          east_vec[count] = easting[b] + xo - 0.5*dem_res;
          north_vec[count] = northing[a] - 0.5*dem_res;
          ++a;
        }
        else if (degs_old > 270 || degs_new < 90)
        {
          xo = 1-0.00001, yo = 1;
          s_edge = abs(s_local*sin((2/PI) - theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo, yi = 1- yo;
          dir = 4;
          east_vec[count] = easting[b] + xo - 0.5*dem_res;
          north_vec[count] = northing[a] + 0.5*dem_res;
          --a;
        }
        else
        {
          cout << "Flow unable to route N or S" << endl;
          exit(EXIT_FAILURE);
        }
      }
      else if (dir == 4) {
        if   (degs_old <= 360 && degs_new >= 180)
        {
          xo = 0, yo = 0.00001;
          s_edge = abs(s_local*sin((PI/2) - theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 3;
          east_vec[count] = easting[b] -0.5*dem_res;
          north_vec[count] = northing[a] + yo - 0.5*dem_res;
          --b;
        }
        else if (degs_old > 0 && degs_new < 180)
        {
          xo = 1, yo = 0.00001;
          s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 1;
          east_vec[count] = easting[b] + 0.5*dem_res;
          north_vec[count] = northing[a] + yo - 0.5*dem_res;
          ++b;
        }
        else
        {
          cout << "Flow unable to route E or W" << endl;
          exit(EXIT_FAILURE);
        }
      }
      slope_total += s_edge*d;

    }
    length += d;
    degs = degs_new;

    //cout << "[a][b]: " << a << " " << b << endl;

    if (a <= 0 || b <= 0 ||  a >= NRows-1 || b >= NCols-1) flag = false;
    else if (stnet[a][b] != ndv || path[a][b] == 1) flag = false;
  }

  //if trace finished at a stream, print hillslope info.
  if (a <= 0 || b <= 0 ||  a >= NRows-1 || b >= NCols-1) continue;
  else
    {
      if (path[a][b] == 1)
        {
          cout << "Didn't make it to a channel!" << endl;
        }

      //          // PRINT TO FILE Cht Sbar Relief Lh
      //          X = xmin + j*dem_res;
      //          Y = ymin + (NRows-i)*dem_res;
      //          relief = zeta[i][j] - zeta[a][b];
      //          length = length*dem_res;
      //          mean_slope = slope_total/(length/dem_res);

      //          ofs << X << " " << Y << " " << seg[i][j] << " "
      //            << cht[i][j] << " " << mean_slope << " "
      //            << relief << " " << length << " " << "/n"; //area[a][b] << "\n";

      //PRINT FILE OF PATH NODES FOR EACH HILLSLOPE VECTOR
      stringstream s;
      s << ht_count;
      file_part_2 = s.str();
      filename = file_part_1;
      filename.append(file_part_2), filename.append(file_part_3);
      strcpy(filename_c,filename.c_str());

      ofstream prof_out;
      prof_out.open(filename_c);
      prof_out << "Easting " << "Northing" << endl;
      prof_out.precision(10);

      for (int c=0;c<count;++c) {
        prof_out   << east_vec[c] << " "
      << north_vec[c] << endl;
      }
      prof_out.close();
    }
      }
      //return condition for debugging purposes only
      //if (ht_count > 50) return;
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting(LSDRaster Elevation, LSDRaster Hilltop_ID, LSDRaster Slope, LSDRaster Aspect, LSDRaster HilltopCurv, LSDRaster PlanCurvature,
                                                         LSDIndexRaster StreamNetwork, LSDIndexRaster Basins,
                                                         string Prefix,
                                                         bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                                                         vector<int> Target_Basin_Vector)
{
  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float s_local;
  //float s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltop_ID = Hilltop_ID.get_RasterData(); //hilltops
  Array2D<float> CHT = HilltopCurv.get_RasterData(); //hilltip curv
  Array2D<float> slope = Slope.get_RasterData(); //hilltop slope
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);

  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;
  vector<double> easting, northing, east_vec, north_vec, zeta_vec, length_vec, empty_vec;

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "easting,northing,i,j,hilltop_id,Cht,S,R,Lh,BasinID,a,b,StreamID,HilltopSlope,DivergentCount,PlanarCountFlag,E_Star,R_Star,EucDist\n";

  //calculate northing and easting
  cout << "XMinimum is " << XMinimum << endl;
  cout << "YMinimum is " << YMinimum << endl;
  cout << "ymax is " << ymax << endl;

  for (i=0;i<NRows;++i) northing.push_back(ymax - DataResolution*(i - 0.5));
  for (j=0;j<NCols;++j) easting.push_back(XMinimum + DataResolution*(j + 0.5));


  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i)
  {
    for (j=0; j<NCols; ++j)
    {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i)
  {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j)
    {

      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltop_ID[i][j] != NoDataValue)
      {

        length = 0;
        flag = true;
        count = 1;
        path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        PlanarCountFlag = 0;
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

        E_Star = 0;
        R_Star = 0;
        EucDist = 0;

        ++ht_count;

        degs = aspect[i][j];
        theta = rads[i][j];
        a = i;
        b = j;
        path[a][b] += 1;
        east_vec = empty_vec;
        north_vec = empty_vec;
        east_vec.push_back(easting[b]);
        north_vec.push_back(northing[a]);
        zeta_vec.push_back(zeta[a][b]);
        length_vec.push_back(length);

        //s_local = slope[a][b];

        //test direction, calculate outlet coordinates and update indicies
        // easterly
        if (degs >= 45 && degs < 135)
        {
          //cout << "\neasterly" << endl;
          xo = 1, yo = (1+tan(theta))/2;
          d = abs(1/(2*cos(theta)));
          xi = 0, yi = yo;
          dir = 1;
          east_vec.push_back(easting[b] + 0.5*DataResolution);
          north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
          ++b;
          if (yi == 0) yi = 0.00001;
          else if (yi == 1) yi = 1 - 0.00001;
        }
        //southerly
        else if (degs >= 135 && degs < 225)
        {
          //cout << "\nsoutherly" << endl;
          xo = (1-(1/tan(theta)))/2, yo = 0;
          d = abs(1/(2*cos((PI/2)-theta)));
          xi = xo, yi = 1;
          dir = 2;
          east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
          north_vec.push_back(northing[a] - 0.5*DataResolution);
          ++a;
          if (xi == 0) xi = 0.00001;
          else if (xi == 1) xi = 1 - 0.00001;
        }
        // westerly
        else if (degs >= 225 && degs < 315)
        {
          xo = 0, yo = (1-tan(theta))/2;
          d = abs(1/(2*cos(theta)));
          xi = 1,  yi = yo;
          dir = 3;
          east_vec.push_back(easting[b] -0.5*DataResolution);
          north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
          --b;
          if (yi == 0) yi = 0.00001;
          else if (yi == 1) yi = 1 - 0.00001;
        }
        //northerly
        else if (degs >= 315 || degs < 45)
        {
          xo = (1+(1/tan(theta)))/2, yo = 1;
          d = abs(1/(2*cos((PI/2) - theta)));
          xi = xo, yi = 0;
          dir = 4;
          east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
          north_vec.push_back(northing[a] + 0.5*DataResolution);
          --a;
          if (xi == 0) xi = 0.00001;
          else if (xi == 1) xi = 1 - 0.00001;
        }
        else
        {
          cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
          exit(EXIT_FAILURE);
        }

        //collect slopes and totals weighted by path length
        length += d;
        //s_local = slope[a][b];

        //update elevation length vectors
        zeta_vec.push_back(zeta[a][b]);
        length_vec.push_back(length*DataResolution);

        //continue trace until a stream node is encountered
        while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1)
        {   //added boudary checking to catch cells which flow off the  edge of the DEM tile.

          path[a][b] += 1;

          degs_new = aspect[a][b];
          theta = rads[a][b];
          ++count;

          //Test for perimeter flow paths
          if ((dir == 1 && degs_new > 0 && degs_new < 180)
              || (dir == 2 && degs_new > 90 && degs_new < 270)
              || (dir == 3 && degs_new > 180 && degs_new < 360)
              || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90)))
          {

            //DO NORMAL FLOW PATH
            //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
            temp_yo1 = yi + (1-xi)*tan(theta);     // xo = 1
            temp_xo1 = xi + (1-yi)*(1/tan(theta));   // yo = 1
            temp_yo2 = yi - xi*tan(theta);      // xo = 0
            temp_xo2 = xi - yi*(1/tan(theta));    // yo = 0

            // can't outlet at same point as inlet
            if (dir == 1) temp_yo2 = -1;
            else if (dir == 2) temp_xo1 = -1;
            else if (dir == 3) temp_yo1 = -1;
            else if (dir == 4) temp_xo2 = -1;

            //s_local = slope[a][b];

            if (temp_yo1 <= 1 && temp_yo1 > 0)
            {
              xo = 1, yo = temp_yo1;
              d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
              xi = 0, yi = yo,
              dir = 1;
              east_vec.push_back(easting[b] + 0.5*DataResolution);
              north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
              ++b;
              if (xi== 0 && yi == 0) yi = 0.00001;
              else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
            }
            else if (temp_xo2 <= 1 && temp_xo2 > 0)
            {
              xo = temp_xo2, yo = 0;
              d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
              xi = xo, yi = 1,
              dir = 2;
              east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
              north_vec.push_back(northing[a] - 0.5*DataResolution);
              ++a;
              if (xi== 0 && yi == 1) xi = 0.00001;
              else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
            }
            else if (temp_yo2 <= 1 && temp_yo2 > 0)
            {
              xo = 0, yo = temp_yo2;
              d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
              xi = 1, yi = yo,
              dir = 3;
              east_vec.push_back(easting[b] -0.5*DataResolution);
              north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
              --b;
              if (xi== 1 && yi == 0) yi = 0.00001;
              else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
            }
            else if (temp_xo1 <= 1 && temp_xo1 > 0)
            {
              xo = temp_xo1, yo = 1;
              d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
              xi = xo, yi = 0,
              dir = 4;
              east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
              north_vec.push_back(northing[a] + 0.5*DataResolution);
              --a;
              if (xi == 0 && yi == 0) xi = 0.00001;
              else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
            }
          }
          else
          {
            // ROUTE ALONG EDGES
            if (dir  == 1)
            {
              if (degs_new <= 90 || degs_new >= 270)
              {
                //secondary compenent of flow is north
                xo = 0.00001, yo = 1;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = xo, yi = 1-yo;
                dir = 4;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] + 0.5*DataResolution);
                --a;
              }
              else if (degs_new > 90 && degs_new < 270) {  //secondary component is south
                xo = 0.00001, yo = 0;
                //s_edge = abs(s_local*sin((PI/2)-theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = xo, yi = 1-yo;
                dir = 2;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] - 0.5*DataResolution);
                ++a;
              }
              else
              {
                cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                //exit(EXIT_FAILURE);
              }
            }
            else if (dir == 2)
            {
              if (degs_new >= 0 && degs_new <= 180)
              {
                //secondary component is East
                xo = 1, yo = 1-0.00001;
                //s_edge = abs(s_local*sin((2/PI)-theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = 1-xo, yi = yo;
                dir = 1;
                east_vec.push_back(easting[b] + 0.5*DataResolution);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                ++b;
              }
              else if (degs_new > 180 && degs_new <= 360)
              {
                //secondary component is West
                xo = 0, yo = 1-0.00001;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = 1-xo, yi = yo;
                dir = 3;
                east_vec.push_back(easting[b] -0.5*DataResolution);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                --b;
              }
              else
              {
                cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                            //exit(EXIT_FAILURE);
              }
            }
            else if (dir == 3)
            {
              if   (degs_new >= 90 && degs_new <= 270)
              {
                //secondary component is South
                xo = 1-0.00001, yo = 0;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = xo, yi = 1-yo;
                dir = 2;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] - 0.5*DataResolution);
                ++a;
              }
              else if (degs_new > 270 || degs_new < 90)
              {
                //secondary component is North
                xo = 1-0.00001, yo = 1;
                //s_edge = abs(s_local*sin((2/PI) - theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = xo, yi = 1- yo;
                dir = 4;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] + 0.5*DataResolution);
                --a;
              }
              else
              {
                cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                          //exit(EXIT_FAILURE);
              }
            }
            else if (dir == 4)
            {
              if   (degs_new >= 180 && degs_new <= 360)
              {
                //secondary component is West
                xo = 0, yo = 0.00001;
                //s_edge = abs(s_local*sin((PI/2) - theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = 1-xo, yi = yo;
                dir = 3;
                east_vec.push_back(easting[b] -0.5*DataResolution);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                --b;
              }
              else if (degs_new >= 0 && degs_new < 180)
              {
                //secondary component is East
                xo = 1, yo = 0.00001;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
                xi = 1-xo, yi = yo;
                dir = 1;
                east_vec.push_back(easting[b] + 0.5*DataResolution);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                ++b;
              }
              else
              {
                cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                            //exit(EXIT_FAILURE);
              }
            }
          }

          if (path[a][b] < 1)
          {  // only update length on 'first slosh'
            length += d;
          }
          else if (path[a][b] >= 3)
          { //update the skip trace flag so we can categorise each trace
             skip_trace = true;
          }
          degs = degs_new;

          zeta_vec.push_back(zeta[a][b]);
          length_vec.push_back(length*DataResolution);

          // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
          // The larger the counter the more convergent or divergent the trace is
          if (abs(PlanCurvature.get_data_element(a,b)) > (0.001))
          {
            ++DivergentCountFlag;
          }
          else
          {
            ++PlanarCountFlag;
          }

          if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
        }

        if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 )
        {
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
        {
        //if trace finished at a stream, print hillslope info.
          if (stnet[a][b] != NoDataValue) // || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
          {
            path[a][b] = 1;

            ++s_count;

            X = easting[j];
            Y = northing[i];
            relief = zeta[i][j] - zeta[a][b];
            mean_slope = relief/(length * DataResolution);

            // update arrays with the current metrics
            RoutedHilltops[i][j] = 1;
            HillslopeLength_Array[i][j] = (length * DataResolution);
            Slope_Array[i][j] = mean_slope;
            Relief_Array[i][j] = relief;

            //calculate an E* and R* Value assuming S_c of 0.8
            E_Star = (2.0 * abs(CHT[i][j])*(length*DataResolution))/0.8;
            R_Star = relief/((length*DataResolution)*0.8);

            //calulate the Euclidean distance between the start and end points of the trace
            EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;

            if (relief > 0)
            {
              ofs << X << "," << Y << "," << i << "," << j << "," << hilltop_ID[i][j] << "," << CHT[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[a][b] << "," << a << "," << b << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "\n";
            }
            else
            {
              ++neg_count;
            }
          }
          else
          {
            //unable to route using aspects
            //this will encompass skipped traces
            //ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
            ++ns_count;
          }
        }

        //This block checks the various path printing options and writes the data out accordingly
        if (print_paths_switch == true)
        {
          if (ht_count % thinning == 0)
          {
            if (hilltop_ID[i][j] != NoDataValue) // && skip_trace == false)
            { //check that the current i,j tuple corresponds to a hilltop, ie there is actually a trace to write to file, and check that the trace was valid.

              //declare some params for lat long conversion
              double latitude,longitude;
              LSDCoordinateConverterLLandUTM Converter;

	            //create the output filename from the user supplied path
              ofstream pathwriter;
              string OutputFileName = Prefix+"_hillslope_traces.csv";
	            ifstream oftest(OutputFileName.c_str());
	            bool FileExists = false;
	            if (oftest) FileExists = true;
	            oftest.close();

	            //open the output filestream and write headers
	            ofstream WriteTracesFile;
	            if (FileExists == 0)
	            {
		            WriteTracesFile.open(OutputFileName.c_str());
                //write headers to allow pandas dataframe reading
                if (WriteTracesFile.is_open()) WriteTracesFile << "HilltopID,Easting,Northing,Latitude,Longitude,Distance,Elevation" << endl;
	            }
	            WriteTracesFile.close();

              //open output filestream again to  coastline data
	            WriteTracesFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	            //Check if file exists if not open a new one and write headers
	            if (WriteTracesFile.is_open())
	            {
		            for (int v = 0; v < count+1; ++v)
		            {
                  //get lat long for printing to file
                  get_lat_and_long_locations(east_vec[v], north_vec[v], latitude, longitude, Converter);

                  if (basin_filter_switch == false)
                  {
                    WriteTracesFile << ht_count << "," << setiosflags(ios::fixed) << setprecision(10) << east_vec[v] << "," << north_vec[v] << "," << latitude << "," << longitude << "," << length_vec[v] << "," << zeta_vec[v] << endl;
                  }
                  else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end())
                  {  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
                    WriteTracesFile << ht_count << "," << setiosflags(ios::fixed) << setprecision(10) << east_vec[v] << "," << north_vec[v] << "," << latitude << "," << longitude << "," << length_vec[v] << "," << zeta_vec[v] << endl;
                  }
                }
              }
              WriteTracesFile.close();
            }
          }
        } // End of path printing logic
      }
    } //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;

  return OutputArrays;
}





















//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting_Refactored(LSDRaster& Elevation, LSDRaster& Hilltop_ID, LSDRaster& Slope,
                                                         LSDRaster& Aspect, LSDRaster& HilltopCurv, LSDRaster& PlanCurvature,
                                                         LSDIndexRaster& StreamNetwork, LSDIndexRaster& Basins,
                                                         string Prefix,
                                                         bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                                                         vector<int> Target_Basin_Vector)
{
  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float s_local;
  //float s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;


  // The array of slopes. Hopefully this will save some sime
  Array2D<float> tangent_raster = Aspect.get_RasterData();

  //Get data arrays from LSDRasters
  // SMM Not sure why we need these
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltop_ID = Hilltop_ID.get_RasterData(); //hilltops
  Array2D<float> CHT = HilltopCurv.get_RasterData(); //hilltop curv
  Array2D<float> slope = Slope.get_RasterData(); //hilltop slope
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols,NoDataValue);
  Array2D<float> tangents(NRows,NCols,NoDataValue);
  Array2D<float> one_over_tangents(NRows,NCols,NoDataValue);



  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);

  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;
  vector<double> easting, northing, east_vec, north_vec, zeta_vec, length_vec, empty_vec;

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "easting,northing,i,j,hilltop_id,Cht,S,R,Lh,BasinID,a,b,StreamID,HilltopSlope,DivergentCount,PlanarCountFlag,E_Star,R_Star,EucDist\n";

  //calculate northing and easting
  cout << "XMinimum is " << XMinimum << endl;
  cout << "YMinimum is " << YMinimum << endl;
  cout << "ymax is " << ymax << endl;

  for (i=0;i<NRows;++i) northing.push_back(ymax - DataResolution*(i - 0.5));
  for (j=0;j<NCols;++j) easting.push_back(XMinimum + DataResolution*(j + 0.5));


  //convert aspects to radians with east as theta = 0/2*pi
  // also converts tangents. This saves time later
  for (i=0; i<NRows; ++i)
  {
    for (j=0; j<NCols; ++j)
    {
      //convert aspects to radians with east as theta = 0/2*pi
      if (aspect[i][j] != NoDataValue)
      {
        rads[i][j] = BearingToRad(aspect[i][j]);
        tangents[i][j] = tan(rads[i][j]);
        one_over_tangents[i][j] = 1/tan(rads[i][j]);
      }
    }
  }

  // some time saving parameters
  float half_datares = 0.5*DataResolution;
  float pi_over_two = M_PI*0.5;


  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i)
  {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j)
    {

      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltop_ID[i][j] != NoDataValue)
      {

        length = 0;
        flag = true;
        count = 1;
        path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        PlanarCountFlag = 0;
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

        E_Star = 0;
        R_Star = 0;
        EucDist = 0;

        ++ht_count;

        degs = aspect[i][j];
        theta = rads[i][j];
        a = i;
        b = j;
        path[a][b] += 1;
        east_vec = empty_vec;
        north_vec = empty_vec;
        east_vec.push_back(easting[b]);
        north_vec.push_back(northing[a]);
        zeta_vec.push_back(zeta[a][b]);
        length_vec.push_back(length);

        //s_local = slope[a][b];

        //test direction, calculate outlet coordinates and update indicies
        // easterly
        if (degs >= 45 && degs < 135)
        {
          //cout << "\neasterly" << endl;
          xo = 1, yo = (1+tangents[i][j])*0.5;
          d = abs(1/(2*cos(theta)));
          xi = 0, yi = yo;
          dir = 1;
          east_vec.push_back(easting[b] + half_datares);
          north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
          ++b;
          if (yi == 0) yi = 0.00001;
          else if (yi == 1) yi = 1 - 0.00001;
        }
        //southerly
        else if (degs >= 135 && degs < 225)
        {
          //cout << "\nsoutherly" << endl;
          xo = (1-one_over_tangents[i][j])*0.5, yo = 0;
          d = abs(1/(2*cos((pi_over_two)-theta)));
          xi = xo, yi = 1;
          dir = 2;
          east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
          north_vec.push_back(northing[a] - half_datares);
          ++a;
          if (xi == 0) xi = 0.00001;
          else if (xi == 1) xi = 1 - 0.00001;
        }
        // westerly
        else if (degs >= 225 && degs < 315)
        {
          xo = 0, yo = (1-tangents[i][j])*0.5;
          d = abs(1/(2*cos(theta)));
          xi = 1,  yi = yo;
          dir = 3;
          east_vec.push_back(easting[b] -half_datares);
          north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
          --b;
          if (yi == 0) yi = 0.00001;
          else if (yi == 1) yi = 1 - 0.00001;
        }
        //northerly
        else if (degs >= 315 || degs < 45)
        {
          xo = (1+one_over_tangents[i][j])*0.5, yo = 1;
          d = abs(1/(2*cos((pi_over_two) - theta)));
          xi = xo, yi = 0;
          dir = 4;
          east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
          north_vec.push_back(northing[a] + half_datares);
          --a;
          if (xi == 0) xi = 0.00001;
          else if (xi == 1) xi = 1 - 0.00001;
        }
        else
        {
          cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
          exit(EXIT_FAILURE);
        }

        //collect slopes and totals weighted by path length
        length += d;
        //s_local = slope[a][b];

        //update elevation length vectors
        zeta_vec.push_back(zeta[a][b]);
        length_vec.push_back(length*DataResolution);

        //continue trace until a stream node is encountered
        while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1)
        {   //added boudary checking to catch cells which flow off the  edge of the DEM tile.

          path[a][b] += 1;

          degs_new = aspect[a][b];
          theta = rads[a][b];
          ++count;

          //Test for perimeter flow paths
          if ((dir == 1 && degs_new > 0 && degs_new < 180)
              || (dir == 2 && degs_new > 90 && degs_new < 270)
              || (dir == 3 && degs_new > 180 && degs_new < 360)
              || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90)))
          {

            //DO NORMAL FLOW PATH
            //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
            temp_yo1 = yi + (1-xi)*tangents[a][b];     // xo = 1
            temp_xo1 = xi + (1-yi)*(one_over_tangents[a][b]);   // yo = 1
            temp_yo2 = yi - xi*tangents[a][b];      // xo = 0
            temp_xo2 = xi - yi*(one_over_tangents[a][b]);    // yo = 0

            // can't outlet at same point as inlet
            if (dir == 1) temp_yo2 = -1;
            else if (dir == 2) temp_xo1 = -1;
            else if (dir == 3) temp_yo1 = -1;
            else if (dir == 4) temp_xo2 = -1;

            //s_local = slope[a][b];

            if (temp_yo1 <= 1 && temp_yo1 > 0)
            {
              xo = 1, yo = temp_yo1;
              d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
              xi = 0, yi = yo,
              dir = 1;
              east_vec.push_back(easting[b] + half_datares);
              north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
              ++b;
              if (xi== 0 && yi == 0) yi = 0.00001;
              else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
            }
            else if (temp_xo2 <= 1 && temp_xo2 > 0)
            {
              xo = temp_xo2, yo = 0;
              d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
              xi = xo, yi = 1,
              dir = 2;
              east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
              north_vec.push_back(northing[a] - half_datares);
              ++a;
              if (xi== 0 && yi == 1) xi = 0.00001;
              else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
            }
            else if (temp_yo2 <= 1 && temp_yo2 > 0)
            {
              xo = 0, yo = temp_yo2;
              d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
              xi = 1, yi = yo,
              dir = 3;
              east_vec.push_back(easting[b] -half_datares);
              north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
              --b;
              if (xi== 1 && yi == 0) yi = 0.00001;
              else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
            }
            else if (temp_xo1 <= 1 && temp_xo1 > 0)
            {
              xo = temp_xo1, yo = 1;
              d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
              xi = xo, yi = 0,
              dir = 4;
              east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
              north_vec.push_back(northing[a] + half_datares);
              --a;
              if (xi == 0 && yi == 0) xi = 0.00001;
              else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
            }
          }
          else
          {
            // ROUTE ALONG EDGES
            if (dir  == 1)
            {
              if (degs_new <= 90 || degs_new >= 270)
              {
                //secondary compenent of flow is north
                xo = 0.00001, yo = 1;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = xo, yi = 1-yo;
                dir = 4;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] + half_datares);
                --a;
              }
              else if (degs_new > 90 && degs_new < 270) {  //secondary component is south
                xo = 0.00001, yo = 0;
                //s_edge = abs(s_local*sin((PI/2)-theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = xo, yi = 1-yo;
                dir = 2;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] - half_datares);
                ++a;
              }
              else
              {
                cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                //exit(EXIT_FAILURE);
              }
            }
            else if (dir == 2)
            {
              if (degs_new >= 0 && degs_new <= 180)
              {
                //secondary component is East
                xo = 1, yo = 1-0.00001;
                //s_edge = abs(s_local*sin((2/PI)-theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = 1-xo, yi = yo;
                dir = 1;
                east_vec.push_back(easting[b] + half_datares);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                ++b;
              }
              else if (degs_new > 180 && degs_new <= 360)
              {
                //secondary component is West
                xo = 0, yo = 1-0.00001;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = 1-xo, yi = yo;
                dir = 3;
                east_vec.push_back(easting[b] -half_datares);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                --b;
              }
              else
              {
                cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                            //exit(EXIT_FAILURE);
              }
            }
            else if (dir == 3)
            {
              if   (degs_new >= 90 && degs_new <= 270)
              {
                //secondary component is South
                xo = 1-0.00001, yo = 0;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = xo, yi = 1-yo;
                dir = 2;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] - half_datares);
                ++a;
              }
              else if (degs_new > 270 || degs_new < 90)
              {
                //secondary component is North
                xo = 1-0.00001, yo = 1;
                //s_edge = abs(s_local*sin((2/PI) - theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = xo, yi = 1- yo;
                dir = 4;
                east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
                north_vec.push_back(northing[a] + half_datares);
                --a;
              }
              else
              {
                cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                          //exit(EXIT_FAILURE);
              }
            }
            else if (dir == 4)
            {
              if   (degs_new >= 180 && degs_new <= 360)
              {
                //secondary component is West
                xo = 0, yo = 0.00001;
                //s_edge = abs(s_local*sin((PI/2) - theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = 1-xo, yi = yo;
                dir = 3;
                east_vec.push_back(easting[b] -half_datares);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                --b;
              }
              else if (degs_new >= 0 && degs_new < 180)
              {
                //secondary component is East
                xo = 1, yo = 0.00001;
                //s_edge = abs(s_local*sin(theta));
                d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = 1-xo, yi = yo;
                dir = 1;
                east_vec.push_back(easting[b] + half_datares);
                north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
                ++b;
              }
              else
              {
                cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
                skip_trace = true;
                            //exit(EXIT_FAILURE);
              }
            }
          }

          if (path[a][b] < 1)
          {  // only update length on 'first slosh'
            length += d;
          }
          else if (path[a][b] >= 3)
          { //update the skip trace flag so we can categorise each trace
             skip_trace = true;
          }
          degs = degs_new;

          zeta_vec.push_back(zeta[a][b]);
          length_vec.push_back(length*DataResolution);

          // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
          // The larger the counter the more convergent or divergent the trace is
          if (abs(PlanCurvature.get_data_element(a,b)) > (0.001))
          {
            ++DivergentCountFlag;
          }
          else
          {
            ++PlanarCountFlag;
          }

          if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
        }

        if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 )
        {
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
        {
        //if trace finished at a stream, print hillslope info.
          if (stnet[a][b] != NoDataValue) // || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
          {
            path[a][b] = 1;

            ++s_count;

            X = easting[j];
            Y = northing[i];
            relief = zeta[i][j] - zeta[a][b];
            mean_slope = relief/(length * DataResolution);

            // update arrays with the current metrics
            RoutedHilltops[i][j] = 1;
            HillslopeLength_Array[i][j] = (length * DataResolution);
            Slope_Array[i][j] = mean_slope;
            Relief_Array[i][j] = relief;

            //calculate an E* and R* Value assuming S_c of 0.8
            E_Star = (2.0 * abs(CHT[i][j])*(length*DataResolution))/0.8;
            R_Star = relief/((length*DataResolution)*0.8);

            //calulate the Euclidean distance between the start and end points of the trace
            EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;

            if (relief > 0)
            {
              ofs << X << "," << Y << "," << i << "," << j << "," << hilltop_ID[i][j] << "," << CHT[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[a][b] << "," << a << "," << b << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "\n";
            }
            else
            {
              ++neg_count;
            }
          }
          else
          {
            //unable to route using aspects
            //this will encompass skipped traces
            //ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
            ++ns_count;
          }
        }

        //This block checks the various path printing options and writes the data out accordingly
        if (print_paths_switch == true)
        {
          if (ht_count % thinning == 0)
          {
            if (hilltop_ID[i][j] != NoDataValue) // && skip_trace == false)
            { //check that the current i,j tuple corresponds to a hilltop, ie there is actually a trace to write to file, and check that the trace was valid.

              //declare some params for lat long conversion
              double latitude,longitude;
              LSDCoordinateConverterLLandUTM Converter;

	            //create the output filename from the user supplied path
              ofstream pathwriter;
              string OutputFileName = Prefix+"_hillslope_traces.csv";
	            ifstream oftest(OutputFileName.c_str());
	            bool FileExists = false;
	            if (oftest) FileExists = true;
	            oftest.close();

	            //open the output filestream and write headers
	            ofstream WriteTracesFile;
	            if (FileExists == 0)
	            {
		            WriteTracesFile.open(OutputFileName.c_str());
                //write headers to allow pandas dataframe reading
                if (WriteTracesFile.is_open()) WriteTracesFile << "HilltopID,Easting,Northing,Latitude,Longitude,Distance,Elevation" << endl;
	            }
	            WriteTracesFile.close();

              //open output filestream again to  coastline data
	            WriteTracesFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	            //Check if file exists if not open a new one and write headers
	            if (WriteTracesFile.is_open())
	            {
		            for (int v = 0; v < count+1; ++v)
		            {
                  //get lat long for printing to file
                  get_lat_and_long_locations(east_vec[v], north_vec[v], latitude, longitude, Converter);

                  if (basin_filter_switch == false)
                  {
                    WriteTracesFile << ht_count << "," << setiosflags(ios::fixed) << setprecision(10) << east_vec[v] << "," << north_vec[v] << "," << latitude << "," << longitude << "," << length_vec[v] << "," << zeta_vec[v] << endl;
                  }
                  else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end())
                  {  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
                    WriteTracesFile << ht_count << "," << setiosflags(ios::fixed) << setprecision(10) << east_vec[v] << "," << north_vec[v] << "," << latitude << "," << longitude << "," << length_vec[v] << "," << zeta_vec[v] << endl;
                  }
                }
              }
              WriteTracesFile.close();
            }
          }
        } // End of path printing logic
      }
    } //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;

  return OutputArrays;
}



















//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
//
// SWDG 12/2/14
// SMM refactored to be more memory efficient 28/01/2021
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::HilltopFlowRouting_TerrifyingRefactored(LSDRaster& Elevation, LSDRaster& Aspect, LSDRaster& PlanCurvature,
                                                         LSDIndexRaster& StreamNetwork, string ridges_csv_name,string Prefix,
                                                         bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch)
{

  // Load the ridges csv
  LSDSpatialCSVReader RidgesData(Elevation, ridges_csv_name);

  // Grab the data elements from the Ridges data
  vector<int> ht_row_vec = RidgesData.data_column_to_int("row");
  vector<int> ht_col_vec = RidgesData.data_column_to_int("col");
  vector<int> basin_vec = RidgesData.data_column_to_int("basin_id");
  vector<float> ht_slope_vec = RidgesData.data_column_to_float("slope");
  vector<float> ht_curv_vec = RidgesData.data_column_to_float("curvature");

  vector<double> ht_lat_vec = RidgesData.get_latitude();
  vector<double> ht_long_vec = RidgesData.get_longitude();


  LSDCoordinateConverterLLandUTM Converter;

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float s_local;
  //float s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  float CHT,HT_slope,trace_length;
  double HT_latitude, HT_longitude;
  int basin_ID;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;
  float ymax = YMinimum + NRows*DataResolution;


  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols,NoDataValue);
  Array2D<float> tangents(NRows,NCols,NoDataValue);
  Array2D<float> one_over_tangents(NRows,NCols,NoDataValue);


  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);

  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;
  vector<double> easting, northing, east_vec, north_vec, zeta_vec, length_vec, empty_vec;

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData_TN.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "latitude,longitude,easting,northing,row,col,Cht,S,R,Lh,BasinID,channel_row,channel_col,StreamID,HilltopSlope,DivergentCount,PlanarCountFlag,E_Star,R_Star,EucDist\n";

  //calculate northing and easting
  cout << "XMinimum is " << XMinimum << endl;
  cout << "YMinimum is " << YMinimum << endl;
  cout << "ymax is " << ymax << endl;

  for (i=0;i<NRows;++i) northing.push_back(ymax - DataResolution*(i - 0.5));
  for (j=0;j<NCols;++j) easting.push_back(XMinimum + DataResolution*(j + 0.5));


  //convert aspects to radians with east as theta = 0/2*pi
  // also converts tangents. This saves time later
  for (i=0; i<NRows; ++i)
  {
    for (j=0; j<NCols; ++j)
    {
      //convert aspects to radians with east as theta = 0/2*pi
      if (Aspect.get_data_element(i,j) != NoDataValue)
      {
        rads[i][j] = BearingToRad(Aspect.get_data_element(i,j));
        tangents[i][j] = tan(rads[i][j]);
        one_over_tangents[i][j] = 1/tan(rads[i][j]);
      }
    }
  }

  // some time saving parameters
  float half_datares = 0.5*DataResolution;
  float pi_over_two = M_PI*0.5;


  // cycle through the hilltops, find hilltops and trace downstream
  bool trace_sucess;
  for (int ht_node = 0; ht_node < int(ht_row_vec.size()); ht_node++)
  {
    cout << flush <<  "\tHilltop pixel: " << ht_node << " of = " << int(ht_row_vec.size()) << "              \r";

    // resest the trace sucess. This only gets set to true if the trace worked.

    trace_sucess = false;
    // Get the hilltop data
    i =  ht_row_vec[ht_node];
    j =  ht_col_vec[ht_node];
    basin_ID = basin_vec[ht_node];

    HT_slope = ht_slope_vec[ht_node];
    CHT = ht_curv_vec[ht_node];

    HT_latitude = ht_lat_vec[ht_node];
    HT_longitude = ht_long_vec[ht_node];

    length = 0;
    flag = true;
    count = 1;
    path = blank.copy();
    DivergentCountFlag = 0; //initialise count of divergent cells in trace
    PlanarCountFlag = 0;
    skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

    E_Star = 0;
    R_Star = 0;
    EucDist = 0;

    ++ht_count;

    degs = Aspect.get_data_element(i,j);
    theta = rads[i][j];
    a = i;
    b = j;
    path[a][b] += 1;
    east_vec = empty_vec;
    north_vec = empty_vec;
    east_vec.push_back(easting[b]);
    north_vec.push_back(northing[a]);
    zeta_vec.push_back(Elevation.get_data_element(a,b));
    length_vec.push_back(length);

    //s_local = slope[a][b];

    //test direction, calculate outlet coordinates and update indicies
    // easterly
    if (degs >= 45 && degs < 135)
    {
      //cout << "\neasterly" << endl;
      xo = 1, yo = (1+tangents[i][j])*0.5;
      d = abs(1/(2*cos(theta)));
      xi = 0, yi = yo;
      dir = 1;
      east_vec.push_back(easting[b] + half_datares);
      north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
      ++b;
      if (yi == 0) yi = 0.00001;
      else if (yi == 1) yi = 1 - 0.00001;
    }
    //southerly
    else if (degs >= 135 && degs < 225)
    {
      //cout << "\nsoutherly" << endl;
      xo = (1-one_over_tangents[i][j])*0.5, yo = 0;
      d = abs(1/(2*cos((pi_over_two)-theta)));
      xi = xo, yi = 1;
      dir = 2;
      east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
      north_vec.push_back(northing[a] - half_datares);
      ++a;
      if (xi == 0) xi = 0.00001;
      else if (xi == 1) xi = 1 - 0.00001;
    }
    // westerly
    else if (degs >= 225 && degs < 315)
    {
      xo = 0, yo = (1-tangents[i][j])*0.5;
      d = abs(1/(2*cos(theta)));
      xi = 1,  yi = yo;
      dir = 3;
      east_vec.push_back(easting[b] -half_datares);
      north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
      --b;
      if (yi == 0) yi = 0.00001;
      else if (yi == 1) yi = 1 - 0.00001;
    }
    //northerly
    else if (degs >= 315 || degs < 45)
    {
      xo = (1+one_over_tangents[i][j])*0.5, yo = 1;
      d = abs(1/(2*cos((pi_over_two) - theta)));
      xi = xo, yi = 0;
      dir = 4;
      east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
      north_vec.push_back(northing[a] + half_datares);
      --a;
      if (xi == 0) xi = 0.00001;
      else if (xi == 1) xi = 1 - 0.00001;
    }
    else
    {
      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
      exit(EXIT_FAILURE);
    }

    //collect slopes and totals weighted by path length
    length += d;
    //s_local = slope[a][b];

    //update elevation length vectors
    zeta_vec.push_back(Elevation.get_data_element(a,b));
    length_vec.push_back(length*DataResolution);

    //continue trace until a stream node is encountered
    while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1)
    {   //added boudary checking to catch cells which flow off the  edge of the DEM tile.

      path[a][b] += 1;

      degs_new = Aspect.get_data_element(a,b);
      theta = rads[a][b];
      ++count;

      //Test for perimeter flow paths
      if ((dir == 1 && degs_new > 0 && degs_new < 180)
          || (dir == 2 && degs_new > 90 && degs_new < 270)
          || (dir == 3 && degs_new > 180 && degs_new < 360)
          || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90)))
      {

        //DO NORMAL FLOW PATH
        //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
        temp_yo1 = yi + (1-xi)*tangents[a][b];     // xo = 1
        temp_xo1 = xi + (1-yi)*(one_over_tangents[a][b]);   // yo = 1
        temp_yo2 = yi - xi*tangents[a][b];      // xo = 0
        temp_xo2 = xi - yi*(one_over_tangents[a][b]);    // yo = 0

        // can't outlet at same point as inlet
        if (dir == 1) temp_yo2 = -1;
        else if (dir == 2) temp_xo1 = -1;
        else if (dir == 3) temp_yo1 = -1;
        else if (dir == 4) temp_xo2 = -1;

        //s_local = slope[a][b];

        if (temp_yo1 <= 1 && temp_yo1 > 0)
        {
          xo = 1, yo = temp_yo1;
          d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
          xi = 0, yi = yo,
          dir = 1;
          east_vec.push_back(easting[b] + half_datares);
          north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
          ++b;
          if (xi== 0 && yi == 0) yi = 0.00001;
          else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
        }
        else if (temp_xo2 <= 1 && temp_xo2 > 0)
        {
          xo = temp_xo2, yo = 0;
          d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
          xi = xo, yi = 1,
          dir = 2;
          east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
          north_vec.push_back(northing[a] - half_datares);
          ++a;
          if (xi== 0 && yi == 1) xi = 0.00001;
          else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
        }
        else if (temp_yo2 <= 1 && temp_yo2 > 0)
        {
          xo = 0, yo = temp_yo2;
          d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
          xi = 1, yi = yo,
          dir = 3;
          east_vec.push_back(easting[b] -half_datares);
          north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
          --b;
          if (xi== 1 && yi == 0) yi = 0.00001;
          else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
        }
        else if (temp_xo1 <= 1 && temp_xo1 > 0)
        {
          xo = temp_xo1, yo = 1;
          d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
          xi = xo, yi = 0,
          dir = 4;
          east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
          north_vec.push_back(northing[a] + half_datares);
          --a;
          if (xi == 0 && yi == 0) xi = 0.00001;
          else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
        }
      }
      else
      {
        // ROUTE ALONG EDGES
        if (dir  == 1)
        {
          if (degs_new <= 90 || degs_new >= 270)
          {
            //secondary compenent of flow is north
            xo = 0.00001, yo = 1;
            //s_edge = abs(s_local*sin(theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = xo, yi = 1-yo;
            dir = 4;
            east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
            north_vec.push_back(northing[a] + half_datares);
            --a;
          }
          else if (degs_new > 90 && degs_new < 270) {  //secondary component is south
            xo = 0.00001, yo = 0;
            //s_edge = abs(s_local*sin((PI/2)-theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = xo, yi = 1-yo;
            dir = 2;
            east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
            north_vec.push_back(northing[a] - half_datares);
            ++a;
          }
          else
          {
            cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
            cout << "Trace skipped.\n" << endl;
            skip_trace = true;
            //exit(EXIT_FAILURE);
          }
        }
        else if (dir == 2)
        {
          if (degs_new >= 0 && degs_new <= 180)
          {
            //secondary component is East
            xo = 1, yo = 1-0.00001;
            //s_edge = abs(s_local*sin((2/PI)-theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = 1-xo, yi = yo;
            dir = 1;
            east_vec.push_back(easting[b] + half_datares);
            north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
            ++b;
          }
          else if (degs_new > 180 && degs_new <= 360)
          {
            //secondary component is West
            xo = 0, yo = 1-0.00001;
            //s_edge = abs(s_local*sin(theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = 1-xo, yi = yo;
            dir = 3;
            east_vec.push_back(easting[b] -half_datares);
            north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
            --b;
          }
          else
          {
            cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
            cout << "Trace skipped.\n" << endl;
            skip_trace = true;
                        //exit(EXIT_FAILURE);
          }
        }
        else if (dir == 3)
        {
          if   (degs_new >= 90 && degs_new <= 270)
          {
            //secondary component is South
            xo = 1-0.00001, yo = 0;
            //s_edge = abs(s_local*sin(theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = xo, yi = 1-yo;
            dir = 2;
            east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
            north_vec.push_back(northing[a] - half_datares);
            ++a;
          }
          else if (degs_new > 270 || degs_new < 90)
          {
            //secondary component is North
            xo = 1-0.00001, yo = 1;
            //s_edge = abs(s_local*sin((2/PI) - theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = xo, yi = 1- yo;
            dir = 4;
            east_vec.push_back(easting[b] + (xo - 0.5)*DataResolution);
            north_vec.push_back(northing[a] + half_datares);
            --a;
          }
          else
          {
            cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
            cout << "Trace skipped.\n" << endl;
            skip_trace = true;
                      //exit(EXIT_FAILURE);
          }
        }
        else if (dir == 4)
        {
          if   (degs_new >= 180 && degs_new <= 360)
          {
            //secondary component is West
            xo = 0, yo = 0.00001;
            //s_edge = abs(s_local*sin((PI/2) - theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = 1-xo, yi = yo;
            dir = 3;
            east_vec.push_back(easting[b] -half_datares);
            north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
            --b;
          }
          else if (degs_new >= 0 && degs_new < 180)
          {
            //secondary component is East
            xo = 1, yo = 0.00001;
            //s_edge = abs(s_local*sin(theta));
            d = sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
            xi = 1-xo, yi = yo;
            dir = 1;
            east_vec.push_back(easting[b] + half_datares);
            north_vec.push_back(northing[a] + (yo - 0.5)*DataResolution);
            ++b;
          }
          else
          {
            cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
            cout << "Trace skipped.\n" << endl;
            skip_trace = true;
                        //exit(EXIT_FAILURE);
          }
        }
      }

      if (path[a][b] < 1)
      {  // only update length on 'first slosh'
        length += d;
      }
      else if (path[a][b] >= 3)
      { //update the skip trace flag so we can categorise each trace
          skip_trace = true;
      }
      degs = degs_new;

      zeta_vec.push_back(Elevation.get_data_element(a,b));
      length_vec.push_back(length*DataResolution);

      // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
      // The larger the counter the more convergent or divergent the trace is
      if (abs(PlanCurvature.get_data_element(a,b)) > (0.001))
      {
        ++DivergentCountFlag;
      }
      else
      {
        ++PlanarCountFlag;
      }

      if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1
                 || StreamNetwork.get_data_element(a,b) != NoDataValue
                 || path[a][b] >= 3 || skip_trace == true)
      {
        flag = false;
      }
    }

    if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 )
    {
      // avoid going out of bounds.

      // this is caused by having a hilltop on the first row or col away from the border
      // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
      // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

      ++edge_count;
    }
    else
    {
    //if trace finished at a stream, print hillslope info.
      if (StreamNetwork.get_data_element(a,b) != NoDataValue)
      {
        path[a][b] = 1;

        ++s_count;

        X = easting[j];
        Y = northing[i];
        relief = Elevation.get_data_element(i,j) - Elevation.get_data_element(a,b);
        trace_length = length * DataResolution;
        mean_slope = relief/(trace_length);


        //calculate an E* and R* Value assuming S_c of 0.8
        E_Star = (2.0 * abs(CHT)*(length*DataResolution))/0.8;
        R_Star = relief/((length*DataResolution)*0.8);

        //calulate the Euclidean distance between the start and end points of the trace
        EucDist = sqrt( ((i+0.5)-(a+yo))*((i+0.5)-(a+yo)) + ((j+0.5)-(b+xo))*((j+0.5)-(b+xo)) ) * DataResolution;

        if (relief > 0)
        {
          trace_sucess = true;
          ofs.precision(9);
          ofs << HT_latitude << "," << HT_longitude << ",";
          ofs.precision(6);
          ofs << X << "," << Y << "," << i << "," << j << "," << CHT
              << "," << mean_slope << "," << relief << "," << length*DataResolution << ","
              << basin_ID << "," << a << "," << b << "," << StreamNetwork.get_data_element(a,b)
              << "," << HT_slope << "," << DivergentCountFlag << "," << PlanarCountFlag << ","
              << E_Star << "," << R_Star << "," << EucDist << "\n";
        }
        else
        {
          ++neg_count;
        }
      }
      else
      {
        //unable to route using aspects
        //this will encompass skipped traces
        //ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
        ++ns_count;
      }
    }

    //This block checks the various path printing options and writes the data out accordingly
    if (print_paths_switch == true)
    {
      if (ht_count % thinning == 0)
      {
        if (trace_sucess) // && skip_trace == false)
        { //check that the current i,j tuple corresponds to a hilltop, ie there is actually a trace to write to file, and check that the trace was valid.

          //declare some params for lat long conversion
          double latitude,longitude;

          //create the output filename from the user supplied path
          ofstream pathwriter;
          string OutputFileName = Prefix+"_hillslope_traces.csv";
          ifstream oftest(OutputFileName.c_str());
          bool FileExists = false;
          if (oftest) FileExists = true;
          oftest.close();

          //open the output filestream and write headers
          ofstream WriteTracesFile;
          if (FileExists == 0)
          {
            WriteTracesFile.open(OutputFileName.c_str());
            //write headers to allow pandas dataframe reading
            if (WriteTracesFile.is_open()) WriteTracesFile << "HilltopID,Easting,Northing,Latitude,Longitude,Distance,Elevation" << endl;
          }
          WriteTracesFile.close();

          //open output filestream again to  coastline data
          WriteTracesFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

          //Check if file exists if not open a new one and write headers
          if (WriteTracesFile.is_open())
          {
            for (int v = 0; v < count+1; ++v)
            {
              //get lat long for printing to file
              get_lat_and_long_locations(east_vec[v], north_vec[v], latitude, longitude, Converter);

              WriteTracesFile << ht_count << "," << setiosflags(ios::fixed) << setprecision(10) << east_vec[v] << ","
                              << north_vec[v] << "," << latitude << "," << longitude << "," << length_vec[v] << ","
                              << zeta_vec[v] << endl;

            }
          }
          WriteTracesFile.close();
        }
      }
    }
  }

  ofs.close();


  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;

}





















































//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN
// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting_RAW(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
                   LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
                   bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                   vector<int> Target_Basin_Vector){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag = 0;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float degs_old;
  //float s_local;
  //float s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = D_inf_Flowdir.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);

  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData_RAW.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope,DivergentCount\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {

      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {

  length = 0;
  flag = true;
  count = 1;
  path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

  ++ht_count;

  degs = aspect[i][j];
  theta = rads[i][j];
  a = i;
  b = j;
  path[a][b] += 1;
  east_vec[0] = easting[b];
  north_vec[0] = northing[a];
  //s_local = slope[a][b];

  //test direction, calculate outlet coordinates and update indicies
  // easterly
  if (degs >= 45 && degs < 135) {
    //cout << "\neasterly" << endl;
    xo = 1, yo = (1+tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 0, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //southerly
  else if (degs >= 135 && degs < 225) {
    //cout << "\nsoutherly" << endl;
    xo = (1-(1/tan(theta)))/2, yo = 0;
    d = abs(1/(2*cos((PI/2)-theta)));
    xi = xo, yi = 1;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  // westerly
  else if (degs >= 225 && degs < 315) {
    xo = 0, yo = (1-tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 1,  yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //northerly
  else if (degs >= 315 || degs < 45) {
    xo = (1+(1/tan(theta)))/2, yo = 1;
    d = abs(1/(2*cos((PI/2) - theta)));
    xi = xo, yi = 0;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  else {
    cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
    exit(EXIT_FAILURE);
  }

  //collect slopes and totals weighted by path length
  length += d;
  //s_local = slope[a][b];

  //continue trace until a stream node is encountered
  while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the  edge of the DEM tile.
    int a_2 = a;
          int b_2 = b;

    path[a][b] += 1;

    //degs_old = degs;
    degs_new = aspect[a][b];
    theta = rads[a][b];
          ++count;

    //Test for perimeter flow paths
    if ((dir == 1 && degs_new > 0 && degs_new < 180)
        || (dir == 2 && degs_new > 90 && degs_new < 270)
        || (dir == 3 && degs_new > 180 && degs_new < 360)
        || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

      //DO NORMAL FLOW PATH
      //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
      temp_yo1 = yi + (1-xi)*tan(theta);     // xo = 1
      temp_xo1 = xi + (1-yi)*(1/tan(theta));   // yo = 1
      temp_yo2 = yi - xi*tan(theta);      // xo = 0
      temp_xo2 = xi - yi*(1/tan(theta));    // yo = 0

      // can't outlet at same point as inlet
      if (dir == 1) temp_yo2 = -1;
      else if (dir == 2) temp_xo1 = -1;
      else if (dir == 3) temp_yo1 = -1;
      else if (dir == 4) temp_xo2 = -1;

      //s_local = slope[a][b];

      if (temp_yo1 <= 1 && temp_yo1 > 0) {
        xo = 1, yo = temp_yo1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 0, yi = yo,
    dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
        if (xi== 0 && yi == 0) yi = 0.00001;
        else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
      }
      else if (temp_xo2 <= 1 && temp_xo2 > 0) {
        xo = temp_xo2, yo = 0;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 1,
    dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
        if (xi== 0 && yi == 1) xi = 0.00001;
        else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
      }
      else if (temp_yo2 <= 1 && temp_yo2 > 0) {
        xo = 0, yo = temp_yo2;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 1, yi = yo,
    dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
        if (xi== 1 && yi == 0) yi = 0.00001;
        else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
      }

      else if (temp_xo1 <= 1 && temp_xo1 > 0) {
        xo = temp_xo1, yo = 1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 0,
    dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
        if (xi == 0 && yi == 0) xi = 0.00001;
        else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
      }

          }
    else {

      // ROUTE ALONG EDGES
      if (dir == 1) {
        if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north
    xo = 0.00001, yo = 1;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
        }
        else if (degs_new > 90 && degs_new < 270) {  //secondary component is south
    xo = 0.00001, yo = 0;
    //s_edge = abs(s_local*sin((PI/2)-theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
        }
        else {
    cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 2) {
        if   (degs_new >= 0 && degs_new <= 180) { //secondary component is East
    xo = 1, yo = 1-0.00001;
    //s_edge = abs(s_local*sin((2/PI)-theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
        }
        else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West
    xo = 0, yo = 1-0.00001;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
        }
        else {
    cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 3) {
        if   (degs_new >= 90 && degs_new <= 270) {  //secondary component is South
    xo = 1-0.00001, yo = 0;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
        }
        else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
    xo = 1-0.00001, yo = 1;
    //s_edge = abs(s_local*sin((2/PI) - theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1- yo;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
        }
        else {
    cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 4) {
        if   (degs_new >= 180 && degs_new <= 360) { //secondary component is West
    xo = 0, yo = 0.00001;
    //s_edge = abs(s_local*sin((PI/2) - theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
        }
        else if (degs_new >= 0 && degs_new < 180) { //secondary component is East
    xo = 1, yo = 0.00001;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
        }
        else {
    cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }

    }

    if (path[a][b] < 1){  // only update length on 'first slosh'
      length += d;
          }

          degs = degs_new;

          if(zeta[a][b] - zeta[a_2][b_2] > 0){

            length -= d;    //remove uphill length from trace

            a = a_2;
            b = b_2;

            //restart trace
            degs = aspect[a][b];
      theta = rads[a][b];
      path[a][b] += 1;
      //s_local = slope[a][b];

            length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

      //test direction, calculate outlet coordinates and update indicies
      // easterly
      if (degs >= 45 && degs < 135) {
        xo = 1, yo = (1+tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 0, yi = yo;
        dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
      }
      //southerly
      else if (degs >= 135 && degs < 225) {
        xo = (1-(1/tan(theta)))/2, yo = 0;
        d = abs(1/(2*cos((PI/2)-theta)));
        xi = xo, yi = 1;
        dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
      }
      // westerly
      else if (degs >= 225 && degs < 315) {
        xo = 0, yo = (1-tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 1,  yi = yo;
        dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
      }
      //northerly
      else if (degs >= 315 || degs < 45) {
        xo = (1+(1/tan(theta)))/2, yo = 1;
        d = abs(1/(2*cos((PI/2) - theta)));
        xi = xo, yi = 0;
        dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
      }
      else {
        cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
        exit(EXIT_FAILURE);
      }

      //collect slopes and totals weighted by path length

      length += d;
      //s_local = slope[a][b];

          }

    if (path[a][b] >= 1){  //self intersect/'slosh'

            degs = aspect[a][b];
      theta = rads[a][b];
      path[a][b] += 1;
      //s_local = slope[a][b];

            a_2 = a;
            b_2 = b;

      length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

      //test direction, calculate outlet coordinates and update indicies
      // easterly
      if (degs >= 45 && degs < 135) {
        xo = 1, yo = (1+tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 0, yi = yo;
        dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
      }
      //southerly
      else if (degs >= 135 && degs < 225) {
        xo = (1-(1/tan(theta)))/2, yo = 0;
        d = abs(1/(2*cos((PI/2)-theta)));
        xi = xo, yi = 1;
        dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
      }
      // westerly
      else if (degs >= 225 && degs < 315) {
        xo = 0, yo = (1-tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 1,  yi = yo;
        dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
      }
      //northerly
      else if (degs >= 315 || degs < 45) {
        xo = (1+(1/tan(theta)))/2, yo = 1;
        d = abs(1/(2*cos((PI/2) - theta)));
        xi = xo, yi = 0;
        dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
      }
      else {
        cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
        exit(EXIT_FAILURE);
      }

      //collect slopes and totals weighted by path length
      length += d;
      //s_local = slope[a][b];

    }

    // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
    // The larger the counter the more convergent or divergent the trace is
    if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)){
      ++DivergentCountFlag;
    }
    else {
      ++PlanarCountFlag;
    }

    if (path[a][b] >=3){ //update flag if a trace cannot complete, so that we can track errors.
      skip_trace = true;
    }

    if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
  }

        if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 ){
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
    {
      //if trace finished at a stream, print hillslope info.
      if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
        {
    path[a][b] = 1;

    ++s_count;

    X = XMinimum + j*DataResolution;
    Y = YMinimum - (NRows-i)*DataResolution;
    relief = zeta[i][j] - zeta[a][b];
    mean_slope = relief/(length * DataResolution);

    // update arrays with the current metrics
    RoutedHilltops[i][j] = 1;
    HillslopeLength_Array[i][j] = (length * DataResolution);
    Slope_Array[i][j] = mean_slope;
    Relief_Array[i][j] = relief;

    //calculate an E* and R* Value assuming S_c of 0.8
    E_Star = (2.0 * abs(hilltops[i][j])*(length*DataResolution))/0.8;
    R_Star = relief/((length*DataResolution)*0.8);

    //calulate the Euclidean distance between the start and end points of the trace
    EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;

    if (relief > 0){
      ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "\n";
    }
    else {
      ++neg_count;
    }
        }
      else{  //unable to route using aspects
        //this will encompass the skipped traces
        ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
        ++ns_count;
      }
    }

  //This block checks the various path printing options and writes the data out accordingly
  if (print_paths_switch == true){
    if (ht_count % thinning == 0){
      if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop and has a valid trace, ie there is actually a trace to write to file.

        //create stringstream object to create filename
        ofstream pathwriter;

        //create the output filename from the user supplied path
        stringstream ss_path;
        ss_path << trace_path << i << "_" << j << "_trace.txt";

        pathwriter.open(ss_path.str().c_str());

        if(pathwriter.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
    exit(EXIT_FAILURE);
        }

        for (int v = 0; v < count+1; ++v){
    if (basin_filter_switch == false){
      pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << endl;
    }
    else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
      pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << endl;
    }
        }
        pathwriter.close();
      }
    }
  }
  // End of path printing logic
      }
    }   //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;

  return OutputArrays;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN
// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// This version is used to generate elevation profiles of the traces. The elevation data is encoded within
// the trace files.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 25/3/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting_Profile(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
                 LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir, string Prefix, LSDIndexRaster Basins,
                 bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                 vector<int> Target_Basin_Vector){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;

  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float degs_old;
  //float s_local;
  //float s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = D_inf_Flowdir.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);

  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);
  vector<float> ZetaList;
  vector<float> LengthList;

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData_RAW.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {

      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {

  length = 0;
  flag = true;
  count = 1;
  path = blank.copy();
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

  ++ht_count;

  degs = aspect[i][j];
  theta = rads[i][j];
  a = i;
  b = j;
  path[a][b] += 1;
  east_vec[0] = easting[b];
  north_vec[0] = northing[a];
  //s_local = slope[a][b];
  ZetaList.clear();
  LengthList.clear();

  //test direction, calculate outlet coordinates and update indicies
  // easterly
  if (degs >= 45 && degs < 135) {
    //cout << "\neasterly" << endl;
    xo = 1, yo = (1+tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 0, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //southerly
  else if (degs >= 135 && degs < 225) {
    //cout << "\nsoutherly" << endl;
    xo = (1-(1/tan(theta)))/2, yo = 0;
    d = abs(1/(2*cos((PI/2)-theta)));
    xi = xo, yi = 1;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  // westerly
  else if (degs >= 225 && degs < 315) {
    xo = 0, yo = (1-tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 1,  yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //northerly
  else if (degs >= 315 || degs < 45) {
    xo = (1+(1/tan(theta)))/2, yo = 1;
    d = abs(1/(2*cos((PI/2) - theta)));
    xi = xo, yi = 0;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  else {
    cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
    exit(EXIT_FAILURE);
  }

  //collect slopes and totals weighted by path length
  length += d;
  //s_local = slope[a][b];

  //continue trace until a stream node is encountered
  while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the  edge of the DEM tile.
    int a_2 = a;
          int b_2 = b;

    path[a][b] += 1;

    //degs_old = degs;
    degs_new = aspect[a][b];
    theta = rads[a][b];
          ++count;

    //Test for perimeter flow paths
    if ((dir == 1 && degs_new > 0 && degs_new < 180)
        || (dir == 2 && degs_new > 90 && degs_new < 270)
        || (dir == 3 && degs_new > 180 && degs_new < 360)
        || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

      //DO NORMAL FLOW PATH
      //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
      temp_yo1 = yi + (1-xi)*tan(theta);     // xo = 1
      temp_xo1 = xi + (1-yi)*(1/tan(theta));   // yo = 1
      temp_yo2 = yi - xi*tan(theta);      // xo = 0
      temp_xo2 = xi - yi*(1/tan(theta));    // yo = 0

      // can't outlet at same point as inlet
      if (dir == 1) temp_yo2 = -1;
      else if (dir == 2) temp_xo1 = -1;
      else if (dir == 3) temp_yo1 = -1;
      else if (dir == 4) temp_xo2 = -1;

      //s_local = slope[a][b];

      if (temp_yo1 <= 1 && temp_yo1 > 0) {
        xo = 1, yo = temp_yo1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 0, yi = yo,
    dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
        if (xi== 0 && yi == 0) yi = 0.00001;
        else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
      }
      else if (temp_xo2 <= 1 && temp_xo2 > 0) {
        xo = temp_xo2, yo = 0;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 1,
    dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
        if (xi== 0 && yi == 1) xi = 0.00001;
        else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
      }
      else if (temp_yo2 <= 1 && temp_yo2 > 0) {
        xo = 0, yo = temp_yo2;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 1, yi = yo,
    dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
        if (xi== 1 && yi == 0) yi = 0.00001;
        else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
      }

      else if (temp_xo1 <= 1 && temp_xo1 > 0) {
        xo = temp_xo1, yo = 1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 0,
    dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
        if (xi == 0 && yi == 0) xi = 0.00001;
        else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
      }

          }
    else {

      // ROUTE ALONG EDGES
      if (dir  == 1) {
        if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north
    xo = 0.00001, yo = 1;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
        }
        else if (degs_new > 90 && degs_new < 270) {  //secondary component is south
    xo = 0.00001, yo = 0;
    //s_edge = abs(s_local*sin((PI/2)-theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
        }
        else {
    cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl;
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 2) {
        if   (degs_new >= 0 && degs_new <= 180) { //secondary component is East
    xo = 1, yo = 1-0.00001;
    //s_edge = abs(s_local*sin((2/PI)-theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
        }
        else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West
    xo = 0, yo = 1-0.00001;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
        }
        else {
    cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 3) {
        if   (degs_new >= 90 && degs_new <= 270) {  //secondary component is South
    xo = 1-0.00001, yo = 0;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
        }
        else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
    xo = 1-0.00001, yo = 1;
    //s_edge = abs(s_local*sin((2/PI) - theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1- yo;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
        }
        else {
    cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 4) {
        if   (degs_new >= 180 && degs_new <= 360) { //secondary component is West
    xo = 0, yo = 0.00001;
    //s_edge = abs(s_local*sin((PI/2) - theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
        }
        else if (degs_new >= 0 && degs_new < 180) { //secondary component is East
    xo = 1, yo = 0.00001;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
        }
        else {
    cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
    skip_trace = true;
                //exit(EXIT_FAILURE);
        }
      }

    }

    if (path[a][b] < 1){  // only update length on 'first slosh'
      length += d;
          }

          degs = degs_new;

          if(zeta[a][b] - zeta[a_2][b_2] > 0){

            length -= d;    //remove uphill length from trace

            a = a_2;
            b = b_2;

            //restart trace
            degs = aspect[a][b];
      theta = rads[a][b];
      path[a][b] += 1;
      //s_local = slope[a][b];

            length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

      //test direction, calculate outlet coordinates and update indicies
      // easterly
      if (degs >= 45 && degs < 135) {
        xo = 1, yo = (1+tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 0, yi = yo;
        dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
      }
      //southerly
      else if (degs >= 135 && degs < 225) {
        xo = (1-(1/tan(theta)))/2, yo = 0;
        d = abs(1/(2*cos((PI/2)-theta)));
        xi = xo, yi = 1;
        dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
      }
      // westerly
      else if (degs >= 225 && degs < 315) {
        xo = 0, yo = (1-tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 1,  yi = yo;
        dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
      }
      //northerly
      else if (degs >= 315 || degs < 45) {
        xo = (1+(1/tan(theta)))/2, yo = 1;
        d = abs(1/(2*cos((PI/2) - theta)));
        xi = xo, yi = 0;
        dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
      }
      else {
        cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
        exit(EXIT_FAILURE);
      }

      //collect slopes and totals weighted by path length

      length += d;
      //s_local = slope[a][b];

          }

    if (path[a][b] >= 1){  //self intersect/'slosh'

            degs = aspect[a][b];
      theta = rads[a][b];
      path[a][b] += 1;
      //s_local = slope[a][b];

            a_2 = a;
            b_2 = b;

      length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

      //test direction, calculate outlet coordinates and update indicies
      // easterly
      if (degs >= 45 && degs < 135) {
        xo = 1, yo = (1+tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 0, yi = yo;
        dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
      }
      //southerly
      else if (degs >= 135 && degs < 225) {
        xo = (1-(1/tan(theta)))/2, yo = 0;
        d = abs(1/(2*cos((PI/2)-theta)));
        xi = xo, yi = 1;
        dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
      }
      // westerly
      else if (degs >= 225 && degs < 315) {
        xo = 0, yo = (1-tan(theta))/2;
        d = abs(1/(2*cos(theta)));
        xi = 1,  yi = yo;
        dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
      }
      //northerly
      else if (degs >= 315 || degs < 45) {
        xo = (1+(1/tan(theta)))/2, yo = 1;
        d = abs(1/(2*cos((PI/2) - theta)));
        xi = xo, yi = 0;
        dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
      }
      else {
        cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
        exit(EXIT_FAILURE);
      }

      //collect slopes and totals weighted by path length
      length += d;
      //s_local = slope[a][b];

    }

    if (path[a][b] >=3){ //update flag if a trace cannot complete, so that we can track errors.
      skip_trace = true;
    }

    ZetaList.push_back(zeta[a][b]);
    LengthList.push_back(length*DataResolution);

    if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
  }

        if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 ){
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
    {
      //if trace finished at a stream, print hillslope info.
      if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
        {
    path[a][b] = 1;

    ++s_count;

    X = XMinimum + j*DataResolution;
    Y = YMinimum - (NRows-i)*DataResolution;
    relief = zeta[i][j] - zeta[a][b];
    mean_slope = relief/(length * DataResolution);

    // update arrays with the current metrics
    RoutedHilltops[i][j] = 1;
    HillslopeLength_Array[i][j] = (length * DataResolution);
    Slope_Array[i][j] = mean_slope;
    Relief_Array[i][j] = relief;

    if (relief > 0){
      ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "\n";
    }
    else {
      ++neg_count;
    }
        }
      else{  //unable to route using aspects
        //this will encompass the skipped traces
        ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
        ++ns_count;
      }
    }

  //This block checks the various path printing options and writes the data out accordingly
  if (print_paths_switch == true){
    if (ht_count % thinning == 0){
      if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop and has a valid trace, ie there is actually a trace to write to file.

        //create stringstream object to create filename
        ofstream pathwriter;

        //create the output filename from the user supplied path
        stringstream ss_path;
        ss_path << trace_path << i << "_" << j << "_trace.txt";

        pathwriter.open(ss_path.str().c_str());

        if(pathwriter.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
    exit(EXIT_FAILURE);
        }

        for (int v = 0; v < count+1; ++v){
    if (basin_filter_switch == false){
      pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << ZetaList[v] << " " << LengthList[v] << endl;
    }
    else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
      pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << ZetaList[v] << " " << LengthList[v] << endl;
    }
        }
        pathwriter.close();
      }
    }
  }
  // End of path printing logic
      }
    }   //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;

  return OutputArrays;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes an integer raster that has nodata and classification
// codes. These are then propagated upslope, so all upslope nodes have
// this code. Nodes that have recieves a classification do not reclassify:
// the downstream node classification "wins"
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDFlowInfo::tag_nodes_upstream_of_baselevel(LSDIndexRaster& tagged_raster)
{
  // We just double check if the dimensions are the same
  int NRows_ir = tagged_raster.get_NRows();
  int NCols_ir = tagged_raster.get_NCols();
  int NDV = tagged_raster.get_NoDataValue();

  Array2D<int> output_array(NRows_ir,NCols_ir,NDV);

  if ( (NRows_ir != NRows) || (NCols_ir != NCols) )
  {
    cout << "Fatal error. You are trying to make a tagged raster using LSDFlowInfo::tag_nodes_upstream_of_baselevel" << endl;
    cout << "but the raster and the flow info object do not have the same rows and columns." << endl;
    exit(0);
  }

  // Okay, now we go through the stack and tag the upstream nodes
  int row,col;
  int this_tag;
  int s_node = 0;
  int this_node = SVector[s_node];


  vector<int> bln = get_BaseLevelNodeList();
  vector<int> upslope_nodes;

  for (int bl = 0; bl< int(bln.size()); bl++)
  {
    this_node = bln[bl];
    retrieve_current_row_and_col(this_node, row, col);
    this_tag = tagged_raster.get_data_element(row,col);

    cout << "baselevel: " << bl << endl;

    if(this_tag != NDV)
    {
      upslope_nodes = get_upslope_nodes(this_node);
      for (int n = 0; n< int(upslope_nodes.size()); n++)
      {
        retrieve_current_row_and_col(upslope_nodes[n], row, col);
        output_array[row][col] = this_tag;
      }
    }
  }


  // now write the mask as an LSDIndexRaster
  LSDIndexRaster final_tagged(NRows,NCols,XMinimum,YMinimum,
          DataResolution,int(NoDataValue),output_array,GeoReferencingStrings);
  return final_tagged;


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes an integer raster that has nodata and classification
// codes. These are then propagated upslope, so all upslope nodes have
// this code. Nodes that have recieved a classification do not reclassify:
// the downstream node classification "wins"
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDFlowInfo::tag_upstream_nodes(LSDIndexRaster& tagged_raster, float crit_upslope_distance)
{
  // We just double check if the dimensions are the same
  int NRows_ir = tagged_raster.get_NRows();
  int NCols_ir = tagged_raster.get_NCols();
  int NDV = tagged_raster.get_NoDataValue();

  Array2D<int> output_array(NRows_ir,NCols_ir,NDV);

  if ( (NRows_ir != NRows) || (NCols_ir != NCols) )
  {
    cout << "Fatal error. You are trying to make a tagged raster using LSDFlowInfo::tag_upslope_nodes" << endl;
    cout << "but the raster and the flow info object do not have the same rows and columns." << endl;
    exit(0);
  }

  // get the flow distance and contributing pixel rasters
  LSDRaster fd = distance_from_outlet();

  // Loop through the stack, collecting tagged pixels.
  // This should have the effect of ordering the pixels
  vector<int> tagged_node_list;
  int row,col;
  for (int si = 0; si<NDataNodes; si++)
  {
    retrieve_current_row_and_col(SVector[si], row, col);

    if( tagged_raster.get_data_element(row,col) != NDV)
    {
      tagged_node_list.push_back(SVector[si]);
    }
  }

  // Now loop through the tagged node list.
  int NTN = int(tagged_node_list.size());
  float tagged_node_FD;
  float node_FD;
  float FD_diff;
  int this_tag;
  vector<int> upslope_nodes;
  for(int tn = 0; tn< NTN; tn++)
  {
    // get the flow distance of this node
    retrieve_current_row_and_col(tagged_node_list[tn], row, col);
    tagged_node_FD = fd.get_data_element(row,col);

    // get the tag
    this_tag = tagged_raster.get_data_element(row,col);

    upslope_nodes = get_upslope_nodes(tagged_node_list[tn]);
    for (int n = 0; n< int(upslope_nodes.size()); n++)
    {
      retrieve_current_row_and_col(upslope_nodes[n], row, col);

      node_FD = fd.get_data_element(row,col);

      FD_diff = node_FD - tagged_node_FD;

      if (FD_diff < crit_upslope_distance)
      {
        output_array[row][col] = this_tag;
      }
    }
  }

  // now write the mask as an LSDIndexRaster
  LSDIndexRaster final_tagged(NRows,NCols,XMinimum,YMinimum,
          DataResolution,int(NoDataValue),output_array,GeoReferencingStrings);
  return final_tagged;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This is for replacing all values upstream (and including) and outlet node
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::replace_upstream_pixels_with_value(LSDRaster& value_raster, int outlet_node, float value, bool overwrite)
{
  // We just double check if the dimensions are the same
  int NRows_ir = value_raster.get_NRows();
  int NCols_ir = value_raster.get_NCols();
  float NDV = value_raster.get_NoDataValue();

  int curr_row,curr_col,curr_node;
  float this_value;

  if ( (NRows_ir != NRows) || (NCols_ir != NCols) )
  {
    cout << "Fatal error. You are trying to make a tagged raster using LSDFlowInfo::tag_upslope_nodes" << endl;
    cout << "but the raster and the flow info object do not have the same rows and columns." << endl;
    exit(0);
  }


  vector<int> upslope_nodes;
  upslope_nodes = get_upslope_nodes_include_outlet(outlet_node);

  int n_nodes = int(upslope_nodes.size());
  for(int n= 0; n<n_nodes; n++)
  {
    curr_node = upslope_nodes[n];
    retrieve_current_row_and_col(curr_node,curr_row,curr_col);

    if(overwrite)
    {
      value_raster.set_data_element(curr_row,curr_col,value);
    }
    else
    {
      this_value = value_raster.get_data_element(curr_row,curr_col);
      if(this_value == NDV)
      {
        value_raster.set_data_element(curr_row,curr_col,value);
      }
    }   
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes an integer raster that has nodata and classification
// codes. These are then propagated downslope, so all downslope nodes have
// this code. Nodes that have receives a classification do not reclassify:
// the downstream node classification "wins"
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDFlowInfo::tag_downstream_nodes(LSDIndexRaster& tagged_raster, float crit_downslope_distance)
{
  // We just double check if the dimensions are the same
  int NRows_ir = tagged_raster.get_NRows();
  int NCols_ir = tagged_raster.get_NCols();
  int NDV = tagged_raster.get_NoDataValue();

  Array2D<int> output_array(NRows_ir,NCols_ir,NDV);

  if ( (NRows_ir != NRows) || (NCols_ir != NCols) )
  {
    cout << "Fatal error. You are trying to make a tagged raster using LSDFlowInfo::tag_downstream_nodes" << endl;
    cout << "but the raster and the flow info object do not have the same rows and columns." << endl;
    exit(0);
  }

  // get the flow distance and contributing pixel rasters
  LSDRaster fd = distance_from_outlet();

  // Loop through the stack, collecting tagged pixels.
  // This should have the effect of ordering the pixels
  vector<int> tagged_node_list;
  int row,col;
  for (int si = 0; si<NDataNodes; si++)
  {
    retrieve_current_row_and_col(SVector[si], row, col);

    if( tagged_raster.get_data_element(row,col) != NDV)
    {
      tagged_node_list.push_back(SVector[si]);
    }
  }

  // Now loop through the tagged node list.
  int NTN = int(tagged_node_list.size());
  int this_node;
  float tagged_node_FD;
  float receiver_FD;
  int reciever_node;
  float FD_diff;
  int this_tag;
  for(int tn = 0; tn< NTN; tn++)
  {
    // get the flow distance of this node
    retrieve_current_row_and_col(tagged_node_list[tn], row, col);
    tagged_node_FD = fd.get_data_element(row,col);

    // get the tag
    this_tag = tagged_raster.get_data_element(row,col);

    this_node = tagged_node_list[tn];
    retrieve_receiver_information(this_node,reciever_node,row,col);
    receiver_FD = fd.get_data_element(row,col);

    FD_diff = tagged_node_FD-receiver_FD;
    while( this_node != reciever_node && FD_diff < crit_downslope_distance )
    {
      this_node = reciever_node;
      retrieve_receiver_information(this_node,reciever_node,row,col);
      receiver_FD = fd.get_data_element(row,col);

      output_array[row][col] = this_tag;



      FD_diff = tagged_node_FD-receiver_FD;
    }
  }

  // now write the mask as an LSDIndexRaster
  LSDIndexRaster final_tagged(NRows,NCols,XMinimum,YMinimum,
          DataResolution,int(NoDataValue),output_array,GeoReferencingStrings);
  return final_tagged;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This takes a node list and a category list and tags all upstream pixes of the 
// *donors* to those tagged pixels. 
// Used to tag pixels draining to a categorised channel
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDFlowInfo::tag_donor_pixels_exclude_nodelist(vector<int> node_list, vector<int> categories)
{

  int current_node, row,col;
  vector<int> donor_list;
  int this_donor;

  // create the data array
  Array2D<int> new_tags(NRows,NCols,int(NoDataValue));

  for(int node = 0; node < int(node_list.size()); node++)
  {
    donor_list = get_donor_nodes( node_list[node]);

    for(int d = 0; d< int(donor_list.size()); d++)
    {
      this_donor = donor_list[d];
      if(find(node_list.begin(), node_list.end(), this_donor) != node_list.end()) 
      {
        // The donor is in the channel. Do nothing
      } 
      else 
      {
        // The donor is not in the channel. Apply tag. 
        retrieve_current_row_and_col(this_donor,row,col);
        new_tags[row][col] = categories[node];
      }
    }
  }

  // Now you need to tag the upstream nodes
  LSDIndexRaster temp_tagged(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,new_tags,GeoReferencingStrings);

  temp_tagged.write_raster("temp_tagged","bil");

  // Now tag these values
  float crit_upslope_distance = 1000000000000;      // big number to get all pixels
  LSDIndexRaster tagged_raster = tag_upstream_nodes(temp_tagged, crit_upslope_distance);
  
  return tagged_raster;  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a starting node, gets all upslope nodes, and determines
// if they are bounded by noddata. Those that are not are eliminated from the
// list so that what remains are nodes that are fully within the
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::basin_edge_extractor(int outlet_node, LSDRaster& Topography)
{

  int n_nodes = (RowIndex.size());
  int i,j;

  vector<int> upslope_nodes;
  if (outlet_node < n_nodes)
  {
    // first get the nodes upstream of the source node
    upslope_nodes = get_upslope_nodes(outlet_node);
  }
  else
  {
    cout << "Fatal error, outlet node doesn't exist." << endl;
    exit(EXIT_FAILURE);
  }

  // stupidly memory intensive but I can't think of a better way to do it without
  // masses of logic statments.
  Array2D<float> BasinData(NRows, NCols, NoDataValue);

  //create subset arrays for just the basin data - this should be rolled into its own method.
  cout << "The number of nodes in this basin is: " <<  upslope_nodes.size() << endl;
  for (int q = 0; q < int(upslope_nodes.size()); ++q)
  {

    retrieve_current_row_and_col(upslope_nodes[q], i, j);
    BasinData[i][j] = upslope_nodes[q];
  }

  vector<int> perim;
  int NDVCount;
  for (int q = 0; q < int(upslope_nodes.size()); ++q)
  {

    retrieve_current_row_and_col(upslope_nodes[q], i, j);
    NDVCount = 0;

    if (i == 0 || j == 0 || i == NRows-1 || j == NRows-1)
    {
      // We are not going to worry about corners since anything
      // with NDVCount > 1 is classed as a potential boundary.
      NDVCount = 3;
    }
    else
    {
      //count border cells that are NDV
      if (BasinData[i-1][j-1] == NoDataValue){ ++NDVCount; }
      if (BasinData[i][j-1] == NoDataValue){ ++NDVCount; }
      if (BasinData[i+1][j-1] == NoDataValue){ ++NDVCount; }
      if (BasinData[i-1][j] == NoDataValue){ ++NDVCount; }
      if (BasinData[i+1][j] == NoDataValue){ ++NDVCount; }
      if (BasinData[i-1][j+1] == NoDataValue){ ++NDVCount; }
      if (BasinData[i][j+1] == NoDataValue){ ++NDVCount; }
      if (BasinData[i+1][j+1] == NoDataValue){ ++NDVCount; }
    }
    if (NDVCount >= 1 && NDVCount < 8)
    {
      perim.push_back(upslope_nodes[q]);
    }
  }

  return perim;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function makes a mask of all the pixels that recieve flow (d8)
// from a pixel that is either nodata or is on the boundary of the DEM
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDFlowInfo::find_cells_influenced_by_nodata(LSDIndexRaster& Bordered_mask,
                  LSDRaster& Topography)
{

  // set up the array
  Array2D<int> influenced_mask(NRows,NCols,int(NoDataValue));
  for(int row = 0; row <NRows; row++)
    {
      for(int col = 0; col<NCols; col++)
  {
    if(Topography.get_data_element(row,col) != NoDataValue)
      {
        influenced_mask[row][col] = 0;
      }
  }
    }


  int curr_node;
  int next_node;
  int next_row,next_col;

  // now loop through every node in the array
  for(int row = 0; row <NRows; row++)
    {
      for(int col = 0; col<NCols; col++)
  {
    if(Topography.get_data_element(row,col) != NoDataValue)
      {
        // this node has data.
        // first see if it has already been tagged
        if(influenced_mask[row][col] != 1)
    {

      //See if it is borderd by a NDV
      if(Bordered_mask.get_data_element(row,col) == 1)
        {
          // it is bordered by nodata. Work your way down the node list
          curr_node = retrieve_node_from_row_and_column(row, col);
          next_node = ReceiverVector[curr_node];

          influenced_mask[row][col] = 1;
          retrieve_current_row_and_col(next_node, next_row, next_col);

          //cout << "I am bordered by NDV, entering search loop" << endl;
          //cout << "Row: " << row <<  " col: " << col << " node: " << curr_node
          //     << " receiver: " << next_node << " next infl mask: "
          //     << influenced_mask[next_row][next_col] << endl;

          // loop until you hit another influenced node or a baselevel node
          while(next_node != curr_node && influenced_mask[next_row][next_col] != 1 )
      {
        curr_node = next_node;
        next_node = ReceiverVector[curr_node];

        // the index here say next row and column but actually this is
        // preserved from the previous loop so is the current node.
        influenced_mask[next_row][next_col] = 1;

        // get the row and column of the receiver
        retrieve_current_row_and_col(next_node, next_row, next_col);

        //cout << "Looping thought influence, next influenced is: "
        //     << influenced_mask[next_row][next_col] << endl;

      }
        }
    }
      }
  }
    }

  // now write the mask as an LSDIndexRaster
  LSDIndexRaster Influence_by_NDV(NRows,NCols,XMinimum,YMinimum,
          DataResolution,int(NoDataValue),influenced_mask,GeoReferencingStrings);
  return Influence_by_NDV;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDFlowInfo::is_upstream_influenced_by_nodata(int nodeindex, LSDRaster& test_raster)
{
  // get all the upslope nodes of this node.
  vector<int> upslope_node_list = get_upslope_nodes(nodeindex);

  int i,j;

  bool flag = false;
  float raster_value;
  float NDV = test_raster.get_NoDataValue();

  // now loop through all these nodes, seeing if any of them is bounded by nodata
  for (int node = 0; node < int(upslope_node_list.size()); node++)
  {
    retrieve_current_row_and_col( upslope_node_list[node] ,i,j);

    //check for edges of the file
    if (i == 0 || i == (NRows - 1) || j == 0 || j == (NCols - 1))
    {
      flag = true;
      return flag;
    }
    else
    {
      for(int ii = -1; ii<=1; ii++)
      {
        for(int jj = -1; jj<=1; jj++)
        {
          raster_value = test_raster.get_data_element(i+ii,j+jj);
          if (raster_value == NDV)
          {
            flag = true;
            return flag;
          }

        }
      }
    }
  }

  return flag;
}


//----------------------------------------------------------------------------------------
// get_raster_values_for_nodes
//----------------------------------------------------------------------------------------
// This function gets the values from a raster corresponding to the given nodes.
vector<float> LSDFlowInfo::get_raster_values_for_nodes(LSDRaster& Raster, vector<int>& node_indices)
{
  int N_nodes = node_indices.size();
  vector<float> return_values(N_nodes,float(NoDataValue));
  int row=0;
  int col=0;
  for(int i = 0; i < N_nodes; ++i)
    {
      if(node_indices[i] == NoDataValue)
  {
    return_values[i] = NoDataValue;
  }
      else
  {
    retrieve_current_row_and_col(node_indices[i],row,col);
    return_values[i] = Raster.get_data_element(row,col);
  }
    }
  return return_values;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This version performs a single trace from a specified node, and routes down until it
// reaches a channel pixel
//
// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN
// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// There are 4 outputs:
// output_trace_coordinates - the output coordinates tracing the flow path
// output_trace_metrics - the metrics derived from the flow routing
//                        (i) X
//                        (i) Y
//                        (i) mean slope
//                        (i) hillslope relief
//                        (i) hillslope length
//                        (i) channel ID @ lower boundary
// output_channel_node -the nodeindex at the bounding channel
// skip_trace - a bool object that specifies whether this trace has routed
// successfully to the channel.
//
// SWDG (adapted by DTM) 23/3/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::D_Inf_single_trace_to_channel(LSDRaster Elevation, int start_node, LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir,
            vector< vector<float> >& output_trace_coordinates, vector<float>& output_trace_metrics,
            int& output_channel_node, bool& skip_trace)
{

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  //   int flag;
  int count = 0;
  //   int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float degs_old;
  //   float s_local, s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  //   bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  //   int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = D_inf_Flowdir.get_RasterData(); //aspect
  //   Array2D<float> slope = Slope.get_RasterData(); //slope

  Array2D<float> rads(NRows,NCols,NoDataValue);
  Array2D<float> path(NRows, NCols,0.0);
  Array2D<float> blank(NRows,NCols,0.0);

  int channel_node = int(NoDataValue);
  vector<float> trace_metrics;
  vector< vector<float> > trace_coordinates;
  vector<float> empty;
  trace_coordinates.push_back(empty);
  trace_coordinates.push_back(empty);

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  //calculate northing and easting
  for (i=0;i<NRows;++i)
    {
      northing[i] = ymax - i*DataResolution - 0.5;
    }
  for (j=0;j<NCols;++j)
    {
      easting[j] = XMinimum + j*DataResolution + 0.5;
    }
  // find node and trace downstream

  // ignore edge cells and non-hilltop cells
  // route initial node by aspect and get outlet coordinates
  int start_row, start_col;
  retrieve_current_row_and_col(start_node,start_row,start_col);
  bool flag = false;
  if (zeta[start_row][start_col] != NoDataValue)
    {
      length = 0;
      flag = true;
      count = 1;
      path = blank.copy();
      //     DivergentCountFlag = 0; //initialise count of divergent cells in trace
      skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

      //     ++ht_count;

      degs = aspect[start_row][start_col];
      theta = BearingToRad(aspect[start_row][start_col]);
      a = start_row;
      b = start_col;
      path[a][b] += 1;
      east_vec[0] = easting[b];
      north_vec[0] = northing[a];
      //     s_local = slope[a][b];

      //test direction, calculate outlet coordinates and update indicies
      // easterly
      if (degs >= 45 && degs < 135)
  {
    //cout << "\neasterly" << endl;
    xo = 1, yo = (1+tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 0, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
      //southerly
      else if (degs >= 135 && degs < 225)
  {
    //cout << "\nsoutherly" << endl;
    xo = (1-(1/tan(theta)))/2, yo = 0;
    d = abs(1/(2*cos((PI/2)-theta)));
    xi = xo, yi = 1;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
      // westerly
      else if (degs >= 225 && degs < 315)
  {
    xo = 0, yo = (1-tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 1,  yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
      //northerly
      else if (degs >= 315 || degs < 45)
  {
    xo = (1+(1/tan(theta)))/2, yo = 1;
    d = abs(1/(2*cos((PI/2) - theta)));
    xi = xo, yi = 0;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
      else
  {
    cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
    exit(EXIT_FAILURE);
  }

      //collect slopes and totals weighted by path length
      length += d;
      //     s_local = slope[a][b];
      // place coordinates into output vector
      trace_coordinates[0].push_back(east_vec[count]);
      trace_coordinates[1].push_back(north_vec[count]);
      //continue trace until a stream node is encountered
      while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1)   //added boudary checking to catch cells which flow off the edge of the DEM tile.
  {
    int a_2 = a;
    int b_2 = b;

    path[a][b] += 1;

    //degs_old = degs;
    degs_new = aspect[a][b];
    theta = BearingToRad(aspect[a][b]);
    ++count;

    //       cout << "TEST1" << endl;
    //Test for perimeter flow paths
    if ((dir == 1 && degs_new > 0 && degs_new < 180)
        || (dir == 2 && degs_new > 90 && degs_new < 270)
        || (dir == 3 && degs_new > 180 && degs_new < 360)
        || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90)))
      {

        //       cout << "TEST1a" << endl;
        //DO NORMAL FLOW PATH
        //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
        temp_yo1 = yi + (1-xi)*tan(theta);      // xo = 1
        temp_xo1 = xi + (1-yi)*(1/tan(theta));  // yo = 1
        temp_yo2 = yi - xi*tan(theta);          // xo = 0
        temp_xo2 = xi - yi*(1/tan(theta));      // yo = 0

        // can't outlet at same point as inlet
        if (dir == 1) temp_yo2 = -1;
        else if (dir == 2) temp_xo1 = -1;
        else if (dir == 3) temp_yo1 = -1;
        else if (dir == 4) temp_xo2 = -1;

        //         s_local = slope[a][b];
        if (temp_yo1 <= 1 && temp_yo1 > 0)
    {
      xo = 1, yo = temp_yo1;
      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
      xi = 0;
      yi = yo;
      dir = 1;
      east_vec[count] = easting[b] + 0.5*DataResolution;
      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
      ++b;
      if (xi== 0 && yi == 0) yi = 0.00001;
      else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
    }
        else if (temp_xo2 <= 1 && temp_xo2 > 0)
    {
      xo = temp_xo2, yo = 0;
      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
      xi = xo;
      yi = 1;
      dir = 2;
      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
      north_vec[count] = northing[a] - 0.5*DataResolution;
      ++a;
      if (xi== 0 && yi == 1) xi = 0.00001;
      else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
    }
        else if (temp_yo2 <= 1 && temp_yo2 > 0)
    {
      xo = 0, yo = temp_yo2;
      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
      xi = 1;
      yi = yo;
      dir = 3;
      east_vec[count] = easting[b] -0.5*DataResolution;
      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
      --b;
      if (xi== 1 && yi == 0) yi = 0.00001;
      else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
    }
        else if (temp_xo1 <= 1 && temp_xo1 > 0)
    {
      xo = temp_xo1, yo = 1;
      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
      xi = xo;
      yi = 0;
      dir = 4;
      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
      north_vec[count] = northing[a] + 0.5*DataResolution;
      --a;
      if (xi == 0 && yi == 0) xi = 0.00001;
      else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
    }

        //       cout << "TEST1_end" << endl;
      }
    else
      {
        // ROUTE ALONG EDGES

        //       cout << "TEST-" << endl;
        if (dir == 1)
    {
      //       cout << "TEST2" << endl;
      if (degs_new <= 90 || degs_new >= 270) //secondary compenent of flow is north
        {
          xo = 0.00001;
          yo = 1;
          //             s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo;
          yi = 1-yo;
          dir = 4;
          east_vec[count] = easting[b] + xo - 0.5*DataResolution;
          north_vec[count] = northing[a] + 0.5*DataResolution;
          --a;
        }
      else if (degs_new > 90 && degs_new < 270)  //secondary component is south
        {
          xo = 0.00001;
          yo = 0;
          //             s_edge = abs(s_local*sin((PI/2)-theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo;
          yi = 1-yo;
          dir = 2;
          east_vec[count] = easting[b] + xo - 0.5*DataResolution;
          north_vec[count] = northing[a] - 0.5*DataResolution;
          ++a;
        }
      else
        {
          cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
          cout << "Trace skipped.\n" << endl;
          skip_trace = true;
          //exit(EXIT_FAILURE);
        }
    }
        else if (dir == 2)
    {
      //       cout << "TEST3" << endl;
      if (degs_new >= 0 && degs_new <= 180) //secondary component is East
        {
          xo = 1, yo = 1-0.00001;
          //             s_edge = abs(s_local*sin((2/PI)-theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 1;
          east_vec[count] = easting[b] + 0.5*DataResolution;
          north_vec[count] = northing[a] + yo - 0.5*DataResolution;
          ++b;
        }
      else if (degs_new > 180 && degs_new <= 360)  //secondary component is West
        {
          xo = 0, yo = 1-0.00001;
          //             s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 3;
          east_vec[count] = easting[b] -0.5*DataResolution;
          north_vec[count] = northing[a] + yo - 0.5*DataResolution;
          --b;
        }
      else
        {
          cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
          cout << "Trace skipped.\n" << endl; //something has gone very wrong...
          skip_trace = true;
          //exit(EXIT_FAILURE);
        }
    }
        else if (dir == 3)
    {

      //           cout << "TEST4" << endl;
      if(degs_new >= 90 && degs_new <= 270)   //secondary component is South
        {
          xo = 1-0.00001;
          yo = 0;
          //             s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo;
          yi = 1-yo;
          dir = 2;
          east_vec[count] = easting[b] + xo - 0.5*DataResolution;
          north_vec[count] = northing[a] - 0.5*DataResolution;
          ++a;
        }
      else if (degs_new > 270 || degs_new < 90)   //secondary component is North
        {
          xo = 1-0.00001, yo = 1;
          //             s_edge = abs(s_local*sin((2/PI) - theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = xo;
          yi = 1- yo;
          dir = 4;
          east_vec[count] = easting[b] + xo - 0.5*DataResolution;
          north_vec[count] = northing[a] + 0.5*DataResolution;
          --a;
        }
      else
        {
          cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
          cout << "Trace skipped.\n" << endl; //something has gone very wrong...
          skip_trace = true;
          //exit(EXIT_FAILURE);
        }
    }
        else if (dir == 4)
    {
      //       cout << "TEST5" << endl;
      if(degs_new >= 180 && degs_new <= 360) //secondary component is West
        {
          xo = 0, yo = 0.00001;
          //             s_edge = abs(s_local*sin((PI/2) - theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo;
          yi = yo;
          dir = 3;
          east_vec[count] = easting[b] -0.5*DataResolution;
          north_vec[count] = northing[a] + yo - 0.5*DataResolution;
          --b;
        }
      else if (degs_new >= 0 && degs_new < 180) //secondary component is East
        {
          xo = 1, yo = 0.00001;
          //             s_edge = abs(s_local*sin(theta));
          d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
          xi = 1-xo, yi = yo;
          dir = 1;
          east_vec[count] = easting[b] + 0.5*DataResolution;
          north_vec[count] = northing[a] + yo - 0.5*DataResolution;
          ++b;
        }
      else
        {
          cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
          cout << "Trace skipped.\n" << endl; //something has gone very wrong...
          skip_trace = true;
          //exit(EXIT_FAILURE);
        }
    }

        //       cout << "TEST6" << endl;
      }

    if (path[a][b] < 1) length += d; // only update length on 'first slosh'

    degs = degs_new;

    if(zeta[a][b] - zeta[a_2][b_2] > 0)
      {
        length -= d;    //remove uphill length from trace

        a = a_2;
        b = b_2;

        //         cout << "TEST7" << endl;
        //restart trace
        degs = aspect[a][b];
        theta = BearingToRad(aspect[a][b]);
        path[a][b] += 1;
        //         s_local = slope[a][b];
        length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

        //test direction, calculate outlet coordinates and update indices easterly
        if (degs >= 45 && degs < 135)
    {
      xo = 1, yo = (1+tan(theta))/2;
      d = abs(1/(2*cos(theta)));
      xi = 0, yi = yo;
      dir = 1;
      east_vec[count] = easting[b] + 0.5*DataResolution;
      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
      ++b;
    }
        //southerly
        else if (degs >= 135 && degs < 225)
    {
      xo = (1-(1/tan(theta)))/2, yo = 0;
      d = abs(1/(2*cos((PI/2)-theta)));
      xi = xo, yi = 1;
      dir = 2;
      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
      north_vec[count] = northing[a] - 0.5*DataResolution;
      ++a;
    }
        // westerly
        else if (degs >= 225 && degs < 315)
    {
      xo = 0, yo = (1-tan(theta))/2;
      d = abs(1/(2*cos(theta)));
      xi = 1,  yi = yo;
      dir = 3;
      east_vec[count] = easting[b] -0.5*DataResolution;
      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
      --b;
    }
        //northerly
        else if (degs >= 315 || degs < 45)
    {
      xo = (1+(1/tan(theta)))/2, yo = 1;
      d = abs(1/(2*cos((PI/2) - theta)));
      xi = xo, yi = 0;
      dir = 4;
      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
      north_vec[count] = northing[a] + 0.5*DataResolution;
      --a;
    }
        else
    {
      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
      exit(EXIT_FAILURE);
    }
        //collect slopes and totals weighted by path length
        length += d;
        //         s_local = slope[a][b];
      }

    if (path[a][b] >= 1)  //self intersect/'slosh'
      {
        degs = aspect[a][b];
        theta = rads[a][b];
        path[a][b] += 1;
        //         s_local = slope[a][b];

        a_2 = a;
        b_2 = b;

        length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

        //test direction, calculate outlet coordinates and update indices
        // easterly
        if (degs >= 45 && degs < 135)
    {
      xo = 1, yo = (1+tan(theta))/2;
      d = abs(1/(2*cos(theta)));
      xi = 0, yi = yo;
      dir = 1;
      east_vec[count] = easting[b] + 0.5*DataResolution;
      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
      ++b;
    }
        //southerly
        else if (degs >= 135 && degs < 225)
    {
      xo = (1-(1/tan(theta)))/2, yo = 0;
      d = abs(1/(2*cos((PI/2)-theta)));
      xi = xo, yi = 1;
      dir = 2;
      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
      north_vec[count] = northing[a] - 0.5*DataResolution;
      ++a;
    }
        // westerly
        else if (degs >= 225 && degs < 315)
    {
      xo = 0, yo = (1-tan(theta))/2;
      d = abs(1/(2*cos(theta)));
      xi = 1,  yi = yo;
      dir = 3;
      east_vec[count] = easting[b] -0.5*DataResolution;
      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
      --b;
    }
        //northerly
        else if (degs >= 315 || degs < 45)
    {
      xo = (1+(1/tan(theta)))/2, yo = 1;
      d = abs(1/(2*cos((PI/2) - theta)));
      xi = xo, yi = 0;
      dir = 4;
      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
      north_vec[count] = northing[a] + 0.5*DataResolution;
      --a;
    }
        else
    {
      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
      exit(EXIT_FAILURE);
    }
        //collect slopes and totals weighted by path length
        length += d;
        //         s_local = slope[a][b];
      }

    // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
    // The larger the counter the more convergent or divergent the trace is
    //       if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)) ++DivergentCountFlag;
    if (path[a][b] >=3) skip_trace = true;//update flag if a trace cannot complete, so that we can track errors.
    if (a == 0 || b == 0 || a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;

    // save trace coordinates for this iteration.
    trace_coordinates[0].push_back(east_vec[count]);
    trace_coordinates[1].push_back(north_vec[count]);
  }

      if (a == 0 || b == 0 || a == NRows-1 || b == NCols-1 )
  {
    // avoid going out of bounds.
    // this is caused by having a hilltop on the first row or col away from the border
    // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
    // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*
    ++edge_count;
  }
      else
  {
    //if trace finished at a stream, print hillslope info.
    if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
      {
        path[a][b] = 1;
        ++s_count;
        X = XMinimum + j*DataResolution;
        Y = YMinimum - (NRows-i)*DataResolution;
        relief = zeta[start_row][start_col] - zeta[a][b];
        mean_slope = relief/(length * DataResolution);

        trace_metrics.push_back(X);
        trace_metrics.push_back(Y);
        trace_metrics.push_back(float(start_node));
        trace_metrics.push_back(mean_slope);
        trace_metrics.push_back(relief);
        trace_metrics.push_back(length*DataResolution);

        if (stnet[a][b] != NoDataValue)
    {
      channel_node = retrieve_node_from_row_and_column(a,b);
      trace_metrics.push_back(float(channel_node));
    }
        // find nearest channel pixel within 1m buffer - if more than one, choose furthest downstream
        else
    {
      float min_elev=NoDataValue;
      if (stnet[a-1][b-1] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a-1][b-1] < min_elev)
      {
        min_elev = zeta[a-1][b-1];
        channel_node = retrieve_node_from_row_and_column(a-1,b-1);
      }
        }
      else if (stnet[a-1][b] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a-1][b] < min_elev)
      {
        min_elev = zeta[a-1][b];
        channel_node = retrieve_node_from_row_and_column(a-1,b);
      }
        }
      else if (stnet[a-1][b+1] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a-1][b+1] < min_elev)
      {
        min_elev = zeta[a-1][b+1];
        channel_node = retrieve_node_from_row_and_column(a-1,b+1);
      }
        }
      else if (stnet[a][b-1] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a][b-1] < min_elev)
      {
        min_elev = zeta[a][b-1];
        channel_node = retrieve_node_from_row_and_column(a,b-1);
      }
        }
      else if (stnet[a][b+1] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a][b+1] < min_elev)
      {
        min_elev = zeta[a][b+1];
        channel_node = retrieve_node_from_row_and_column(a,b+1);
      }
        }
      else if (stnet[a+1][b-1] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a+1][b-1] < min_elev)
      {
        min_elev = zeta[a+1][b-1];
        channel_node = retrieve_node_from_row_and_column(a+1,b-1);
      }
        }
      else if (stnet[a+1][b] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a+1][b] < min_elev)
      {
        min_elev = zeta[a+1][b];
        channel_node = retrieve_node_from_row_and_column(a+1,b);
      }
        }
      else if (stnet[a+1][b+1] != NoDataValue)
        {
          if (min_elev == NoDataValue || zeta[a+1][b+1] < min_elev)
      {
        min_elev = zeta[a+1][b+1];
        channel_node = retrieve_node_from_row_and_column(a+1,b+1);
      }
        }
      trace_metrics.push_back(float(channel_node));
    }
        //         if (relief > 0) ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "\n";
        //         else ++neg_count;
        if (relief <= 0) ++neg_count;
      }
    else
      { //unable to route using aspects
        //this will encompass the skipped traces
        //         ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
        ++ns_count;
      }
  }
    }
  output_trace_coordinates = trace_coordinates;
  output_trace_metrics = trace_metrics;
  output_channel_node = channel_node;
}
//----------------------------------------------------------------------------------------------------------------------



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
// (4) Fraction Rock Exposure
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRoutingBedrock(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
                LSDIndexRaster StreamNetwork, LSDRaster Aspect, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
                bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                vector<int> Target_Basin_Vector, LSDRaster RockExposure){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  float rock_exposure;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  //float s_local;
  //float s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins
  Array2D<float> rock = RockExposure.get_RasterData(); // Rock Exposure

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);
  Array2D<float> Rock_Array(NRows,NCols,NoDataValue);

  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope,DivergentCount\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {

      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {

  length = 0;
  rock_exposure = 0;
  flag = true;
  count = 1;
  path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        PlanarCountFlag = 0;
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

        E_Star = 0;
        R_Star = 0;
        EucDist = 0;

  ++ht_count;

  degs = aspect[i][j];
  theta = rads[i][j];
  a = i;
  b = j;
  path[a][b] += 1;
  east_vec[0] = easting[b];
  north_vec[0] = northing[a];
  //s_local = slope[a][b];

  //test direction, calculate outlet coordinates and update indicies
  // easterly
  if (degs >= 45 && degs < 135) {
    //cout << "\neasterly" << endl;
    xo = 1, yo = (1+tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 0, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //southerly
  else if (degs >= 135 && degs < 225) {
    //cout << "\nsoutherly" << endl;
    xo = (1-(1/tan(theta)))/2, yo = 0;
    d = abs(1/(2*cos((PI/2)-theta)));
    xi = xo, yi = 1;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  // westerly
  else if (degs >= 225 && degs < 315) {
    xo = 0, yo = (1-tan(theta))/2;
    d = abs(1/(2*cos(theta)));
    xi = 1,  yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
    if (yi == 0) yi = 0.00001;
    else if (yi == 1) yi = 1 - 0.00001;
  }
  //northerly
  else if (degs >= 315 || degs < 45) {
    xo = (1+(1/tan(theta)))/2, yo = 1;
    d = abs(1/(2*cos((PI/2) - theta)));
    xi = xo, yi = 0;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
    if (xi == 0) xi = 0.00001;
    else if (xi == 1) xi = 1 - 0.00001;
  }
  else {
    cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
    exit(EXIT_FAILURE);
  }

  //collect slopes and totals weighted by path length
  length += d;
  //s_local = slope[a][b];
  rock_exposure += rock[a][b]*d;

  //continue trace until a stream node is encountered
  while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the  edge of the DEM tile.

    path[a][b] += 1;

    degs_new = aspect[a][b];
    theta = rads[a][b];
    ++count;

    //Test for perimeter flow paths
    if ((dir == 1 && degs_new > 0 && degs_new < 180)
        || (dir == 2 && degs_new > 90 && degs_new < 270)
        || (dir == 3 && degs_new > 180 && degs_new < 360)
        || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

      //DO NORMAL FLOW PATH
      //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
      temp_yo1 = yi + (1-xi)*tan(theta);     // xo = 1
      temp_xo1 = xi + (1-yi)*(1/tan(theta));   // yo = 1
      temp_yo2 = yi - xi*tan(theta);      // xo = 0
      temp_xo2 = xi - yi*(1/tan(theta));    // yo = 0

      // can't outlet at same point as inlet
      if (dir == 1) temp_yo2 = -1;
      else if (dir == 2) temp_xo1 = -1;
      else if (dir == 3) temp_yo1 = -1;
      else if (dir == 4) temp_xo2 = -1;

      //s_local = slope[a][b];

      if (temp_yo1 <= 1 && temp_yo1 > 0) {
        xo = 1, yo = temp_yo1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 0, yi = yo,
    dir = 1;
        east_vec[count] = easting[b] + 0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        ++b;
        if (xi== 0 && yi == 0) yi = 0.00001;
        else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
      }
      else if (temp_xo2 <= 1 && temp_xo2 > 0) {
        xo = temp_xo2, yo = 0;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 1,
    dir = 2;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] - 0.5*DataResolution;
        ++a;
        if (xi== 0 && yi == 1) xi = 0.00001;
        else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
      }
      else if (temp_yo2 <= 1 && temp_yo2 > 0) {
        xo = 0, yo = temp_yo2;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = 1, yi = yo,
    dir = 3;
        east_vec[count] = easting[b] -0.5*DataResolution;
        north_vec[count] = northing[a] + yo - 0.5*DataResolution;
        --b;
        if (xi== 1 && yi == 0) yi = 0.00001;
        else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
      }

      else if (temp_xo1 <= 1 && temp_xo1 > 0) {
        xo = temp_xo1, yo = 1;
        d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
        xi = xo, yi = 0,
    dir = 4;
        east_vec[count] = easting[b] + xo - 0.5*DataResolution;
        north_vec[count] = northing[a] + 0.5*DataResolution;
        --a;
        if (xi == 0 && yi == 0) xi = 0.00001;
        else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
      }

    }
    else {

      // ROUTE ALONG EDGES
      if (dir  == 1) {
        if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north
    xo = 0.00001, yo = 1;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
        }
        else if (degs_new > 90 && degs_new < 270) {  //secondary component is south
    xo = 0.00001, yo = 0;
    //s_edge = abs(s_local*sin((PI/2)-theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
        }
        else {
    cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
    cout << "Trace skipped.\n" << endl;
    skip_trace = true;
    //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 2) {
        if   (degs_new >= 0 && degs_new <= 180) { //secondary component is East
    xo = 1, yo = 1-0.00001;
    //s_edge = abs(s_local*sin((2/PI)-theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
        }
        else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West
    xo = 0, yo = 1-0.00001;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
        }
        else {
    cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
    cout << "Trace skipped.\n" << endl;
    skip_trace = true;
    //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 3) {
        if   (degs_new >= 90 && degs_new <= 270) {  //secondary component is South
    xo = 1-0.00001, yo = 0;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1-yo;
    dir = 2;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] - 0.5*DataResolution;
    ++a;
        }
        else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
    xo = 1-0.00001, yo = 1;
    //s_edge = abs(s_local*sin((2/PI) - theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = xo, yi = 1- yo;
    dir = 4;
    east_vec[count] = easting[b] + xo - 0.5*DataResolution;
    north_vec[count] = northing[a] + 0.5*DataResolution;
    --a;
        }
        else {
    cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
    cout << "Trace skipped.\n" << endl;
    skip_trace = true;
    //exit(EXIT_FAILURE);
        }
      }
      else if (dir == 4) {
        if   (degs_new >= 180 && degs_new <= 360) { //secondary component is West
    xo = 0, yo = 0.00001;
    //s_edge = abs(s_local*sin((PI/2) - theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 3;
    east_vec[count] = easting[b] -0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    --b;
        }
        else if (degs_new >= 0 && degs_new < 180) { //secondary component is East
    xo = 1, yo = 0.00001;
    //s_edge = abs(s_local*sin(theta));
    d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
    xi = 1-xo, yi = yo;
    dir = 1;
    east_vec[count] = easting[b] + 0.5*DataResolution;
    north_vec[count] = northing[a] + yo - 0.5*DataResolution;
    ++b;
        }
        else {
    cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
    cout << "Trace skipped.\n" << endl;
    skip_trace = true;
    //exit(EXIT_FAILURE);
        }
      }

    }

    if (path[a][b] < 1){  // only update length on 'first slosh'
      length += d;
      rock_exposure += rock[a][b]*d;
    }
    else if (path[a][b] >= 3){ //update the skip trace flag so we can categorise each trace
      skip_trace = true;
    }

    degs = degs_new;

    // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
    // The larger the counter the more convergent or divergent the trace is
    if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)){
      ++DivergentCountFlag;
    }
    else {
      ++PlanarCountFlag;
    }

    if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
  }

  if (a == 0 || b == 0 ||  a == NRows-1 || b == NCols-1 ){
    // avoid going out of bounds.

    // this is caused by having a hilltop on the first row or col away from the border
    // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
    // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

    ++edge_count;

  }
  else
    {
      //if trace finished at a stream, print hillslope info.
      if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
        {
    path[a][b] = 1;

    ++s_count;

    X = XMinimum + j*DataResolution;
    Y = YMinimum - (NRows-i)*DataResolution;
    relief = zeta[i][j] - zeta[a][b];
    mean_slope = relief/(length * DataResolution);

    // update arrays with the current metrics
    RoutedHilltops[i][j] = 1;
    HillslopeLength_Array[i][j] = (length * DataResolution);
    Slope_Array[i][j] = mean_slope;
    Relief_Array[i][j] = relief;
    Rock_Array[i][j] = rock_exposure/length;

    //calculate an E* and R* Value assuming S_c of 0.8
    E_Star = (2.0 * abs(hilltops[i][j])*(length*DataResolution))/0.8;
    R_Star = relief/((length*DataResolution)*0.8);

    //calulate the Euclidean distance between the start and end points of the trace
    EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;


    if (relief > 0){
      ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "," << rock_exposure/length << "\n";
    }
    else {
      ++neg_count;
    }
        }
      else{  //unable to route using aspects
        //this will encompass skipped traces
        ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
        ++ns_count;
      }
    }

  //This block checks the various path printing options and writes the data out accordingly
  if (print_paths_switch == true){
    if (ht_count % thinning == 0){
      if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop, ie there is actually a trace to write to file, and check that the trace was valid.

        //create stringstream object to create filename
        ofstream pathwriter;

        //create the output filename from the user supplied path
        stringstream ss_path;
        ss_path << trace_path << i << "_" << j << "_trace.txt";

        pathwriter.open(ss_path.str().c_str());

        if(pathwriter.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
    exit(EXIT_FAILURE);
        }

        for (int v = 0; v < count+1; ++v){
    if (basin_filter_switch == false){
      pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << "," << rock_exposure/length << endl;
    }
    else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
      pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << "," << rock_exposure/length << endl;
    }
        }
        pathwriter.close();
      }
    }
  }
  // End of path printing logic
      }
    }   //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);
  OutputArrays.push_back(Rock_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;

  return OutputArrays;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This method removes end nodes which are not the uppermost extent of the channel network.
// SWDG 23/7/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::ProcessEndPointsToChannelHeads(LSDIndexRaster Ends){

  Array2D<int> EndArray = Ends.get_RasterData();
  vector<int> Sources;

  //make a map containing each nodeindex where true means it is a valid channel head, eg the top of the network
  map<int,bool> EndStatus;
  vector<int> EndNodes;

  for(int i=1; i<NRows-1; ++i){
    for(int j=1; j<NCols-1; ++j){
      if (EndArray[i][j] != NoDataValue){
        int nodeindex = retrieve_node_from_row_and_column (i,j);
        EndStatus[nodeindex] = true;
        EndNodes.push_back(nodeindex);
      }
    }
  }

  for (int q = 0; q < int(EndNodes.size());++q){
    cout << flush << q << " of " << EndNodes.size() << "\r";
    int CurrentNode = EndNodes[q];
    if (EndStatus[CurrentNode] == true){

      bool stop = false;

      while (stop == false){
        int DownslopeNode;
        int Downslopei;
        int Downslopej;

        //get steepest descent neighbour
        retrieve_receiver_information(CurrentNode,DownslopeNode,Downslopei,Downslopej);

        if (find(EndNodes.begin(), EndNodes.end(), DownslopeNode) != EndNodes.end()){
          EndStatus[DownslopeNode] = false;
          stop = true;
        }

        //check for out of bounds
        if (Downslopei == 0 || Downslopei == NRows - 1 || Downslopej == 0 || Downslopej == NCols - 1){
          stop = true;
        }

        //check for a node with no downslope neughbours
        if (CurrentNode == DownslopeNode){
          stop = true;
        }

        CurrentNode = DownslopeNode;

      }

    }
  }
  cout << endl;

  for (int q = 0; q < int(EndNodes.size());++q){
    if (EndStatus[EndNodes[q]] == true){
      if (EndNodes[q] != NoDataValue){
        Sources.push_back(EndNodes[q]);
      }
    }
  }

  return Sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This method removes single pixel channels from a channel network.
// SWDG 23/7/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::RemoveSinglePxChannels(LSDIndexRaster StreamNetwork, vector<int> Sources){

  for (int q = 0; q < int(Sources.size());++q){

    int CurrentNode = Sources[q];
    int Currenti;
    int Currentj;
    retrieve_current_row_and_col(CurrentNode,Currenti,Currentj);
    int CurrentOrder = StreamNetwork.get_data_element(Currenti,Currentj);

    //get steepest descent neighbour
    int DownslopeNode;
    int Downslopei;
    int Downslopej;
    retrieve_receiver_information(CurrentNode,DownslopeNode,Downslopei,Downslopej);
    int DownslopeOrder = StreamNetwork.get_data_element(Downslopei,Downslopej);

    if (CurrentOrder != DownslopeOrder){
      //remove the value from the list of nodes -> Sources is passed by val, so this will not change values in sources outide this method
      Sources.erase(remove(Sources.begin(), Sources.end(), Sources[q]), Sources.end());
    }

  }

  return Sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function starts from a given node index and then goes downstream
// until it either hits a baselevel node or until it has accumulated a
// number of visited pixels
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDFlowInfo::get_downslope_node_after_fixed_visited_nodes(int source_node,
                 int outlet_node, int n_nodes_to_visit, LSDIndexRaster& VisitedRaster)
{
  int n_visited = 0;
  int current_node, receiver_node,row, col;
  int bottom_node;

  bool Am_I_at_the_bottom_of_the_channel = false;

  current_node = source_node;

  // you start from the source node and work your way downstream
  while( Am_I_at_the_bottom_of_the_channel == false )
  {
    // get the receiver node
    retrieve_receiver_information(current_node,receiver_node, row,col);

    // check if this is a base level node
    if (current_node == receiver_node)
    {
      Am_I_at_the_bottom_of_the_channel = true;
      bottom_node = receiver_node;
    }
    else if (receiver_node == outlet_node)
    {
      Am_I_at_the_bottom_of_the_channel = true;
      bottom_node = receiver_node;
    }
    else
    {
      // check to see if this node has been visited, if so increment the n_visited
      // iterator
      if (VisitedRaster.get_data_element(row,col) == 1)
      {
        n_visited++;
      }
      else
      {
        VisitedRaster.set_data_element(row, col, 1);
      }

      // see if we have collected enough nodes to visit
      if (n_visited >= n_nodes_to_visit && n_nodes_to_visit != -9999)
      {
        Am_I_at_the_bottom_of_the_channel = true;
        bottom_node = receiver_node;
      }
    }

    current_node = receiver_node;
  }
  return bottom_node;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function gets the flow length between two nodes.
// FJC 29/09/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDFlowInfo::get_flow_length_between_nodes(int UpstreamNode, int DownstreamNode)
{
	float length = 0;
	float root_2 = 1.4142135623;

  if (UpstreamNode != DownstreamNode)
  {
  	int upstream_test = is_node_upstream(DownstreamNode, UpstreamNode);
  	if (upstream_test != 1)
  	{
  		cout << "FlowInfo 7430: FATAL ERROR: The selected node is not upstream" << endl;
      length = float(NoDataValue);
  	}
    else
    {
    	bool ReachedChannel = false;
    	int CurrentNode = UpstreamNode;
    	while (ReachedChannel == false)
    	{
    		//get receiver information
    		int ReceiverNode, ReceiverRow, ReceiverCol;
    		retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
    		//if node is at baselevel then exit
    		if (CurrentNode == ReceiverNode)
    		{
    			ReachedChannel = true;
    			//cout << "You reached a baselevel node, returning baselevel" << endl;
    		}
    		//if receiver is a channel > threshold then get the stream order
    		if (ReceiverNode == DownstreamNode)
    		{
    			ReachedChannel = true;
    		}
    		else
    		{
    			//move downstream
    			CurrentNode = ReceiverNode;
    			// update length
    			if (retrieve_flow_length_code_of_node(ReceiverNode) == 1){ length += DataResolution; }
          else if (retrieve_flow_length_code_of_node(ReceiverNode) == 2){ length += (DataResolution * root_2); }
    		}
    	}
    }
  }
	return length;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function gets the Euclidian distance between two nodes in metres
// FJC 17/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDFlowInfo::get_Euclidian_distance(int node_A, int node_B)
{
  int row_A, row_B, col_A, col_B;
  // get the row and cols of the nodes
  retrieve_current_row_and_col(node_A, row_A, col_A);
  retrieve_current_row_and_col(node_B, row_B, col_B);

  float row_length = (row_B - row_A)*DataResolution;
  float col_length = (col_B - col_A)*DataResolution;

  //find the distance between these nodes
  float dist = sqrt(row_length*row_length + col_length*col_length);

  return dist;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Snap a given point to the nearest hilltop pixel, within a search radius.
// Returns the nodeindex of the snapped point.
// SWDG 23/1/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDFlowInfo::snap_to_hilltop(int a, int b, int search_radius, LSDRaster& Hilltops){

  int tmpNode;

  if (Hilltops.get_data_element(a,b) != NoDataValue){
    // The point is already on a hilltop pixel!
    tmpNode = retrieve_node_from_row_and_column(a,b);
  }
  else{

    vector<int> Nodes_in_window;
    vector<float> Dists_in_window;
    vector<float> Dists_in_window_sorted;
    vector<size_t> index_map;
    Nodes_in_window.reserve(4 * search_radius);
    Dists_in_window.reserve(4 * search_radius);
    Dists_in_window_sorted.reserve(4 * search_radius);
    index_map.reserve(4 * search_radius);

    //set up the bounding box
    int a_min = a - search_radius;
    int a_max = a + search_radius;
    int b_min = b - search_radius;
    int b_max = b + search_radius;

    //out of bounds checking
    if (a_min < 0){a_min = 0;}
    if (b_min < 0){b_min = 0;}
    if (a_max > (NRows - 1)){a_max = (NRows - 1);}
    if (b_max > (NCols - 1)){b_max = (NCols - 1);}

    // only iterate over the search area.
    for (int i = a_min; i < a_max; ++i){
      for (int j = b_min; j < b_max; ++j){

        if (Hilltops.get_data_element(i, j) != NoDataValue){

          //get the nodeindex and distance from user defined point for each cell in the search window
          tmpNode = retrieve_node_from_row_and_column(i,j);
          Nodes_in_window.push_back(tmpNode);

          float Dist = distbetween(a,b,i,j);
          Dists_in_window.push_back(Dist);

        }
      }
    }

  matlab_float_sort(Dists_in_window, Dists_in_window_sorted, index_map);

  //the hilltop node with the smallest distance to the user defined point
  tmpNode = Nodes_in_window[index_map[0]];

  }

  return tmpNode;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Wrapper around snap_to_hilltop function to process a collection of utm points.
// Writes the the nodeindex of each snapped point to SnappedNodes and the
// coordinate count (first coordinate pair is 0, second is 1 and so on) is written
// to Valid_node_IDs.
// SWDG 23/1/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::snap_to_hilltops(vector<float> x_locs, vector<float> y_locs, int search_radius, LSDRaster& Hilltops, vector<int>& SnappedNodes, vector<int>& Valid_node_IDs){

  for (int q = 0; q < int(x_locs.size()); ++q){

    bool is_in_raster = check_if_point_is_in_raster(x_locs[q], y_locs[q]);

    if (is_in_raster){

      // Shift origin to that of dataset
      float X_coordinate_shifted_origin = x_locs[q] - XMinimum;
      float Y_coordinate_shifted_origin = y_locs[q] - YMinimum;

      // Get row and column of point
      int col_point = int(X_coordinate_shifted_origin/DataResolution);
      int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

      int tmpNode = snap_to_hilltop(row_point, col_point, search_radius, Hilltops);
      SnappedNodes.push_back(tmpNode);
      Valid_node_IDs.push_back(q);

    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get slope between nodes
// FJC 03/05/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDFlowInfo::get_slope_between_nodes(int upslope_node, int downslope_node, LSDRaster& Elevation)
{
  float slope = NoDataValue;
  int upslope_row, upslope_col, downslope_row, downslope_col;
  bool us_node = is_node_upstream(downslope_node, upslope_node);
  if (us_node == false)
  {
    cout << "Warning! Your downslope node is not downslope of the upslope one. Returning NDV." << endl;
  }
  else
  {
    retrieve_current_row_and_col(upslope_node, upslope_row, upslope_col);
    retrieve_current_row_and_col(downslope_node, downslope_row, downslope_col);
    float upslope_elev = Elevation.get_data_element(upslope_row, upslope_col);
    float downslope_elev = Elevation.get_data_element(downslope_row, downslope_col);
    float FlowDist = get_flow_length_between_nodes(upslope_node, downslope_node);

    slope = (upslope_elev - downslope_elev)/FlowDist;
  }
  return slope;
}



/// @brief Structure used to sort vectors of flow
struct vector_PQ
{
  float base_elevation;
  vector<int> nodes;
};

bool operator>( const vector_PQ& lhs, const vector_PQ& rhs )
{
  return lhs.base_elevation > rhs.base_elevation;
}
bool operator<( const vector_PQ& lhs, const vector_PQ& rhs )
{
  return lhs.base_elevation < rhs.base_elevation;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get separated vector of flow for the entire landscape
// BG 10/2018
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector< vector<int> > LSDFlowInfo::get_vectors_of_flow(LSDRaster& topo)
{
  vector<int> No_donors ;
  vector<vector<int> > out_vector;
  priority_queue< vector_PQ, vector<vector_PQ>, greater<vector_PQ> > processing_PQ;

  // Preparing a raster to check who has been processed
  Array2D<short> processed(NRows,NCols,NoDataValue);
  for(int i=0;i<NRows;i++)
  {
    for(int j=0;j<NCols;j++)
    {
      if(topo.get_data_element(i,j)!=NoDataValue)
      {
        processed[i][j] = 0;
      }
    }
  }
  // ofstream  data_out;
  // string filename = "/home/s1675537/PhD/LSDTopoData/knickpoint/test_location_paper/Smugglers_SC/Response_To_WS/TESTAGUL.csv";
  // data_out.open(filename.c_str());
  // data_out << "row,col,x,y,z" << endl;
  // data_out.precision(9);
  for(int i=0; i<int(RowIndex.size());i++)
  {
    int tRow = RowIndex[i];
    int tCol = ColIndex[i];
    int NDonors_here = NDonorsVector[i];
    int this_node = int(i);

    if(NDonors_here == 0)
    {
      float x,y;
      get_x_and_y_locations(tRow,tCol,x,y);
      // data_out << tRow << ",";
      // data_out << tCol << ",";
      // data_out << x << ",";
      // data_out << y << ",";
      // data_out << topo.get_data_element(tRow,tCol) << endl;
      No_donors.push_back(this_node);
    }

  }
  // data_out.close();
  int cpt = 1, ma = No_donors.size();
  for(vector<int>::iterator olga=No_donors.begin();olga!=No_donors.end(); olga++)
  {
    cout << "Processing #" << cpt <<"/"<<ma << "\r";
    cpt++;
    int this_node = *olga;
    vector<int> this_vector;
    this_vector.push_back(this_node);
    int row,col;
    retrieve_current_row_and_col(this_node,row,col);
    processed[row][col] = 1;
    bool proc = true;
    do
    {
      int this_receiver;
      retrieve_receiver_information(this_node,this_receiver,row,col);

      if(this_node != this_receiver && processed[row][col] != 1)
      {
        processed[row][col] = 1;
        this_vector.push_back(this_receiver);
        this_node = this_receiver;
      }
      else
      {
        reverse(this_vector.begin(),this_vector.end());
        out_vector.push_back(this_vector);
        proc = false;
        this_vector.clear();
      }
    }
    while(proc);


  }
  cout << endl;
  // OK now I am getting all the vector in order
  // the iterator is just the name of a book in front of me
  cout << "I got the flow vectors, I am now sorting it per base elevation to make sure it is processing first the bases" << endl;
  for(vector<vector<int> >::iterator Introduction_to_Differentail_Equations_by_Mark_h_Holmes = out_vector.begin(); Introduction_to_Differentail_Equations_by_Mark_h_Holmes != out_vector.end(); Introduction_to_Differentail_Equations_by_Mark_h_Holmes++)
  {
    vector<int> this_vector = *Introduction_to_Differentail_Equations_by_Mark_h_Holmes;
    vector_PQ this_inclusion;
    int row,col;
    int this_node = this_vector[0];
    retrieve_current_row_and_col(this_node,row,col);
    float this_elevation = topo.get_data_element(row,col);
    this_inclusion.base_elevation = this_elevation;
    this_inclusion.nodes = this_vector;
    processing_PQ.push(this_inclusion);
  }
  out_vector.clear();

  while(processing_PQ.size()>0)
  {
    vector_PQ this_inclusion;
    this_inclusion = processing_PQ.top();
    processing_PQ.pop();
    vector<int> this_vector =  this_inclusion.nodes;
    out_vector.push_back(this_vector);
  }


  return out_vector;
}


// Alright, for the following operation, I will need to create a structure:
// I want to sort my base-level nodes per ascending order to make sure I multithread the largest basins first
struct base_level_node_sorter
{
  int Node_ID;
  int n_donors_to_that_node;
};
//Overload the less than and greater than operators to consider Zeta data only
//N.B. Fill only needs greater than but less than useful for mdflow routing
//(I've coded this but not yet added to LSDRaster, it's only faster than presorting
//when applied to pretty large datasets).
bool operator>( const base_level_node_sorter& lhs, const base_level_node_sorter& rhs )
{
  return lhs.n_donors_to_that_node > rhs.n_donors_to_that_node;
}
bool operator<( const base_level_node_sorter& lhs, const base_level_node_sorter& rhs )
{
  return lhs.n_donors_to_that_node < rhs.n_donors_to_that_node;
}

// This function attemps to return a user friendly vectors of different element calculated in the FlowInfo object
// If you question the relevance of that function, try to get the stack vector of node from Braun and Willett
// just reading the comments of that file.
// the returned map has the several string keys:
// "stack order" -> the stack order sense Braun and willett, per node index
// ...
// B.G. 14/04/2019
map<string, vector< vector<int> > > LSDFlowInfo::get_map_of_vectors()
{

  // I first, I want the base levels again
  priority_queue< base_level_node_sorter, vector<base_level_node_sorter>, less<base_level_node_sorter> > PriorityQueue;

  vector<int> base_levels;
  int check_my_vectors = 0;
  for(size_t i=0; i<ReceiverVector.size(); i++)
  {
    check_my_vectors++;
    // checking if this is a baselevel AND has some receivers
    if(int(i) == ReceiverVector[i] && retrieve_ndonors_to_node(i) > 1)
    {
      base_level_node_sorter BLNS;BLNS.Node_ID = int(i);BLNS.n_donors_to_that_node = retrieve_ndonors_to_node(i);
      PriorityQueue.push( BLNS);
    }
  }

  while(PriorityQueue.size()>0)
  {
    base_levels.push_back(PriorityQueue.top().Node_ID);
    PriorityQueue.pop();
  }

  // cout << "Gougne: " << base_levels.size() << " || vs: " << check_my_vectors << endl;

  vector<vector<int> > nodes_per_basins(base_levels.size()), inverted_nodes_per_basins(base_levels.size()), n_row(base_levels.size()), n_col(base_levels.size()), i_row(base_levels.size()), i_col(base_levels.size());
  for(size_t i=0;i<base_levels.size();i++)
  {
    // First get the stack_order of that base_level
    vector<int> this_vec = get_upslope_nodes(base_levels[i]);

    if(this_vec.size() > 1)
    {
      this_vec.insert(this_vec.begin(),base_levels[i]); // Do not forget dat first node!!!
      nodes_per_basins.push_back(this_vec);
      vector<int> rev_vec;
      for(int j= int(this_vec.size()-1); j>=0; j--)
      {
        int new_val = this_vec[j];
        rev_vec.push_back(new_val);
      }

      // reverse(rev_vec.begin(),rev_vec.end());
      inverted_nodes_per_basins.push_back(rev_vec);
      // Now getting the row and col
      vector<int> this_n_row, this_n_col, this_i_row, this_i_col;
      for(size_t j = 0; j< this_vec.size(); j++)
      {
        int this_node = this_vec[j];
        int i_node = rev_vec[j];

        int row,col;
        retrieve_current_row_and_col(this_node,row,col);
        this_n_row.push_back(row);
        this_n_col.push_back(col);

        retrieve_current_row_and_col(i_node,row,col);
        this_i_row.push_back(row);
        this_i_col.push_back(col);
      }

      n_row.push_back(this_n_row);
      n_col.push_back(this_n_col);
      i_row.push_back(this_i_row);
      i_col.push_back(this_i_col);
    }

  }


  // formatting the output
  map<string, vector< vector<int> > > output;

  // adding a tracer
  vector<vector<int> > tracer;
  int cpt = 0;
  for(size_t i = 0; i<size_t(nodes_per_basins.size()); i++ )
  {
    vector<int> goris2;
    for(size_t j=0; j<nodes_per_basins[i].size();j++)
      {
        goris2.push_back(cpt);
        cpt ++;
      }
      tracer.push_back(goris2);
  }

  output["stack_order"] = nodes_per_basins;
  output["inverted_stack_order"] = inverted_nodes_per_basins;
  output["rows"] = n_row;
  output["cols"] = n_col;
  output["inverted_rows"] = i_row;
  output["inverted_cols"] = i_col;
  output["tracer"] = tracer;

  // DEBUGGING, KEEP FOR A BIT------------------------------------------------------------------------------------------------------------
  // float goulg = 0;
  // Array2D<float> checker(NRows,NCols,goulg);
  // for(size_t i=0; i<i_row.size(); i++)
  // {
  //   for(size_t j=0; j<i_row[i].size() ; j++)
  //   {
  //     int row,col;
  //     row = i_row[i][j];
  //     col = i_col[i][j];
  //     checker[row][col] = 1;
  //   }
  // }

  // // Done
  // int cpt2 = 0;
  // for(size_t i=0;i<NRows;i++)
  // {
  //   for(size_t j=0;j<NCols;j++)
  //   {
  //     if(checker[i][j] == 0)
  //       cpt2++;
  //   }
  // }

  // cout << "DEBUG::ORERING::There are " << cpt2 << " 0s out of " << NRows*NCols  <<  endl;
  // DEBUGGING, KEEP FOR A BIT------------------------------------------------------------------------------------------------------------


  return output;

}


// This function returns a raster draining with all the nodes draining to an outlet.
// rest is no data
//B.G
LSDRaster LSDFlowInfo::get_raster_draining_to_node(int node, LSDRaster& elevation_raster)
{
  // getting all the node indices draining to that point
  vector<int> nodes_in_basin = this->get_upslope_nodes_include_outlet(node);

  int min_row = std::numeric_limits<int>::max();
  int max_row = std::numeric_limits<int>::min();
  int min_col = std::numeric_limits<int>::max();
  int max_col = std::numeric_limits<int>::min();

  // first I am getting the min and max row
  for(auto n:nodes_in_basin)
  {

    int this_row,this_col; this->retrieve_current_row_and_col(n,this_row,this_col);
    if(this_row<min_row)
      min_row = this_row;
    else if(this_row > max_row)
      max_row = this_row;
    if(this_col<min_col)
      min_col = this_col;
    else if(this_col > max_col)
      max_col = this_col;
  }


  // creating an empty raster
  Array2D<float> this_raster_data(max_row - min_row + 1, max_col - min_col + 1, NoDataValue);

  for(auto n:nodes_in_basin)
  {
    // getting row col of main raster
    int global_row, global_col; this->retrieve_current_row_and_col(n,global_row,global_col);
    int new_row = global_row - min_row;
    int new_col =  global_col - min_col;

    this_raster_data[new_row][new_col] = elevation_raster.get_data_element(global_row,global_col);
  }


  float new_XLL = XMinimum + min_col * DataResolution;
  float new_YLL = YMinimum + (NRows - 1 - max_row) * DataResolution;

  LSDRaster output(max_row - min_row + 1, max_col - min_col + 1, new_XLL, new_YLL,
      DataResolution, NoDataValue, this_raster_data, GeoReferencingStrings);

  return output;


}

// This function takes a nodelist and determines which of the nodes in the
// nodelist are internal nodes
// Used to get basin outlines
vector<bool> LSDFlowInfo::internal_nodes(vector<int> node_list)
{
  int n_nodes = int(node_list.size());

  // this speeds up the searching
  unordered_set<int> node_set(node_list.begin(),node_list.end());
  vector<bool> is_internal(n_nodes,false);

  //for(int i = 0; i<n_nodes; i++)
  //{
    //node_set.push_back(node_list[i]);
  //  is_internal.push_back(false);
  //}

  // okay, now loop through all the nodes, looking for internal nodes.
  int row,col,trow,tcol,tni;
  bool this_node_internal;
  for(int i = 0; i<n_nodes; i++)
  {
    retrieve_current_row_and_col(node_list[i],row,col);
    this_node_internal = true;

    if (row != 0 && row != NRows-1 && col != 0 && col != NCols -1)
    {
      // check all the adjacent nodes
      trow = row+1;
      tcol = col+1;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row+1;
      tcol = col;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row+1;
      tcol = col-1;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row;
      tcol = col+1;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row;
      tcol = col-1;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row-1;
      tcol = col+1;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row-1;
      tcol = col;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      trow = row-1;
      tcol = col-1;
      tni = NodeIndex[trow][tcol];
      if(node_set.find(tni) == node_set.end())
      {
        this_node_internal = false;
      }

      if(this_node_internal)
      {
        is_internal[i] = true;
      }
    }
  }
  return is_internal;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Search through a vector of nodes and find the nearest node to the target one. 
// FJC
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDFlowInfo::find_nearest_node_in_vector(int target_node, vector<int> node_vector)
{
  int this_node, nearest_node;
  float this_dist;
  float min_dist = 10000000000000;
  for (int i = 0; i < int(node_vector.size()); i++)
  {
    this_node = node_vector[i];
    this_dist = get_Euclidian_distance(this_node, target_node);
    if (this_dist < min_dist) 
    { 
      min_dist = this_dist;
      nearest_node = this_node; 
    }
  }
  return nearest_node;

}

//-----------------------------------------------------------------------------//
// Function to calculate the upstream hypsometric integral from a node
// HI = (Z_mean - Z_min) / (Z_max - Z_min)
//-----------------------------------------------------------------------------//
float LSDFlowInfo::calculate_upstream_HI(int node, LSDRaster& Elevation)
{
  float HI;
  // get the vector of upstream nodes
  vector<int> upstream_nodes = get_upslope_nodes_include_outlet(node);
  // get a vector of elevation of each of these nodes
  int this_elev, this_row, this_col;
  vector<float> upstream_elevs;
  for (int i = 0; i < int(upstream_nodes.size()); i++)
  {
    retrieve_current_row_and_col(upstream_nodes[i], this_row, this_col);
    this_elev = Elevation.get_data_element(this_row, this_col);
    upstream_elevs.push_back(this_elev);
  }

  // now work out the HI from the mean, max and min elevations
  float mean_elev = get_mean_ignore_ndv(upstream_elevs, NoDataValue);
  float min_elev = Get_Minimum(upstream_elevs, NoDataValue);
  float max_elev = Get_Maximum(upstream_elevs, NoDataValue);
  HI = (mean_elev - min_elev)/(max_elev - min_elev);

  return HI;
}

//-----------------------------------------------------------------------------//
// Function to calculate the % upstream elevation within specific elevation 
// bands for a node
// FJC
//-----------------------------------------------------------------------------//
void LSDFlowInfo::calculate_upstream_elevations(int node, LSDRaster& Elevation, float bin_width, float lower_limit, float upper_limit, vector<float>& Midpoints, vector<float>& ProbabilityDensity)
{
  float HI;
  // get the vector of upstream nodes
  vector<int> upstream_nodes = get_upslope_nodes_include_outlet(node);
  // get a vector of elevation of each of these nodes
  int this_elev, this_row, this_col;
  vector<float> upstream_elevs;
  for (int i = 0; i < int(upstream_nodes.size()); i++)
  {
    retrieve_current_row_and_col(upstream_nodes[i], this_row, this_col);
    this_elev = Elevation.get_data_element(this_row, this_col);
    upstream_elevs.push_back(this_elev);
  }

  // calculate histogram from the vector of upstream elevations
  vector<float> LLims, ULims;
  vector<int> Count;
  calculate_histogram_fixed_limits(upstream_elevs, bin_width, lower_limit, upper_limit, Midpoints, LLims, ULims, Count, ProbabilityDensity);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates the elevation difference between each point in the swath and
// the channel
// FJC
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDRaster LSDFlowInfo::calculate_channel_relief(vector <LSDRaster>& swath_rasters, LSDRaster& topography_raster, vector<int> channel_node_list)
// {
//   Array2D<float> ChannelRelief_array(NRows, NCols, NoDataValue);
//
//   int channel_node, channel_row, channel_col, channel_idx;
//   float channel_elev, this_elev;
//   for (int row = 0; row < NRows; row++)
//   {
//     for (int col = 0; col < NCols; col++)
//     {
//       channel_idx = swath_rasters[2].get_data_element(row, col);
//       if (channel_idx != NoDataValue)
//       {
//         channel_node = channel_node_list[channel_idx];
//         retrieve_current_row_and_col(channel_node, channel_row, channel_col);
//         channel_elev = topography_raster.get_data_element(channel_row, channel_col);
//         this_elev = topography_raster.get_data_element(row, col);
//         ChannelRelief_array[row][col] = this_elev - channel_elev;
//       }
//     }
//   }
//
//   LSDRaster ChannelRelief(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ChannelRelief_array, GeoReferencingStrings);
//
//   return ChannelRelief;
// }


#endif
