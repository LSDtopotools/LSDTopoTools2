//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDFloodplain.cpp
//
// Land Surface Dynamics Floodplain Object
//
// This object creates and stores information about floodplain and terraces
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//
// Developed by:
//  Fiona J. Clubb
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
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
#ifndef LSDFloodplain_CPP
#define LSDFloodplain_CPP

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDFloodplain.hpp"
#include "LSDSpatialCSVReader.hpp"
using namespace std;
using namespace TNT;

// OPEN CV
#include "opencv2/opencv.hpp"

using namespace cv;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// create
// this populates the binary array and connected components array for the floodplain
// given rasters of channel relief and slope and thresholds for both. Each pixel
// must be below the slope and channel relief threshold to be classified as floodplain.
// User must set a minimum floodplain patch size (in pixels, set to 0 if all patches are kept).
//
// FJC 18/10/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFloodplain::create(LSDRaster& ChannelRelief, LSDRaster& Slope, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, float relief_thresh, float slope_thresh, int min_patch_size, int threshold_SO)
{

  /// set the protected variables
  NRows = ChannelRelief.get_NRows();
  NCols = ChannelRelief.get_NCols();
  XMinimum = ChannelRelief.get_XMinimum();
  YMinimum = ChannelRelief.get_YMinimum();
  DataResolution = ChannelRelief.get_DataResolution();
  NoDataValue = ChannelRelief.get_NoDataValue();
  GeoReferencingStrings = ChannelRelief.get_GeoReferencingStrings();

  relief_threshold = relief_thresh;
  slope_threshold = slope_thresh;

  //declare the arrays
  Array2D<int> TempBinArray(NRows,NCols,0);
  Array2D<int> TempLinkArray(NRows,NCols,NoDataValue);
  BinaryArray = TempBinArray.copy();
  BinaryMSArray = TempBinArray.copy();
  ConnectedComponents_Array = TempLinkArray.copy();
  Array2D<int> StreamOrderArray = ChanNetwork.get_StreamOrderArray();
  FloodplainNodes_array = TempLinkArray.copy();

  //declare the vectors
  vector<int> FloodplainNodes_temp, patch_ids_channel;

  //loop through every row and col and get the slope and relief values
  for (int i =0; i < NRows; i++)
  {
    for (int j = 0; j < NCols; j++)
    {
      if (ChannelRelief.get_data_element(i,j) != NoDataValue && Slope.get_data_element(i,j) != NoDataValue)
      {
        float slope = Slope.get_data_element(i,j);
        float relief = ChannelRelief.get_data_element(i,j);
        if (relief < relief_threshold && slope < slope_threshold)        //floodplain points must be lower than both the relief
        {                                                                //and the slope threshold.
          BinaryArray[i][j] = 1;
        }
        if (StreamOrderArray[i][j] != NoDataValue)                    // assign all parts of the channel network as part of the valley/floodplain network.
        {
          BinaryArray[i][j] = 2;
        }
      }
    }
  }
  LSDIndexRaster FloodplainRaster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,BinaryArray,GeoReferencingStrings);

  // run the connected components algorithm on the floodplain array
  LSDIndexRaster ConnectedComponents = FloodplainRaster.ConnectedComponents();
  if (min_patch_size > 0)
  {
    LSDIndexRaster ConnectedComponents_final = ConnectedComponents.RemoveSmallPatches(min_patch_size);
    ConnectedComponents_Array = ConnectedComponents_final.get_RasterData();
  }
  else
  {
    ConnectedComponents_Array = ConnectedComponents.get_RasterData();
  }

  // separate into floodplain and terrace patches

  //loop through the DEM and get the ID of all patches connected to the channel network
  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (ConnectedComponents_Array[row][col] != NoDataValue)
      {
      //check if the pixel is part of the channel network
        if (StreamOrderArray[row][col] >= threshold_SO)
        {
          patch_ids_channel.push_back(ConnectedComponents_Array[row][col]);
        }
      }
    }
  }

  //for each pixel, find if it is in a patch with an ID in patch_ids_channel vector
  vector<int>::iterator find_it;
  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (ConnectedComponents_Array[row][col] != NoDataValue)
      {
        int ThisNode = FlowInfo.retrieve_node_from_row_and_column(row, col);
        int patch_id = ConnectedComponents_Array[row][col];
        find_it = find(patch_ids_channel.begin(), patch_ids_channel.end(), patch_id);   //search ID vector for patch ID of pixel
        if (find_it != patch_ids_channel.end())
        {
          FloodplainNodes_temp.push_back(ThisNode);
          FloodplainNodes_array[row][col] = ThisNode;
        }
      }
    }
  }

  //copy to vectors
  FloodplainNodes = FloodplainNodes_temp;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// overloaded create function to read in an existing floodplain raster.
// FJC 04/11/21
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFloodplain::create(LSDRaster& ExistingRaster)
{

  /// set the protected variables
  NRows = ExistingRaster.get_NRows();
  NCols = ExistingRaster.get_NCols();
  XMinimum = ExistingRaster.get_XMinimum();
  YMinimum = ExistingRaster.get_YMinimum();
  DataResolution = ExistingRaster.get_DataResolution();
  NoDataValue = ExistingRaster.get_NoDataValue();
  GeoReferencingStrings = ExistingRaster.get_GeoReferencingStrings();

  //declare the arrays
  Array2D<int> TempBinArray(NRows,NCols,0);
  BinaryArray = TempBinArray.copy();

  for (int i = 0; i < NRows; i++)
  {
    for (int j = 0; j < NCols; j++)
    {
      BinaryArray[i][j] = ExistingRaster.get_data_element(i,j);
    }
  }

}

void LSDFloodplain::create()
{
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function finds the nearest channel greater than a threshold stream order
// to each patch.
// Floodplains - just finds elevation of the nearest channel > threshold SO
// FJC 21/10/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFloodplain::Get_Relief_of_Nearest_Channel(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, int threshold_SO)
{
  //set up the arrays
  Array2D<int> TempIntArray(NRows,NCols,NoDataValue);
  Array2D<float> TempFloatArray(NRows,NCols,NoDataValue);
  NearestChannelElev_array = TempFloatArray.copy();
  ChannelRelief_array = TempFloatArray.copy();
  UpstreamDist_array = TempFloatArray.copy();
  NearestChannelNode_array = TempIntArray.copy();

  for (int i = 0; i < int(FloodplainNodes.size()); i++)
  {
    int row, col, ChannelRow, ChannelCol, ChannelNode;
    FlowInfo.retrieve_current_row_and_col(FloodplainNodes[i], row, col);
    float FlowLength, DistanceUpstream;
    ChanNetwork.get_info_nearest_channel_to_node(FloodplainNodes[i], threshold_SO, FlowInfo, DistFromOutlet, ChannelNode, FlowLength, DistanceUpstream);
    FlowInfo.retrieve_current_row_and_col(ChannelNode, ChannelRow, ChannelCol);
    NearestChannelElev_array[row][col] = ElevationRaster.get_data_element(ChannelRow,ChannelCol);
    NearestChannelNode_array[row][col] = ChannelNode;
    UpstreamDist_array[row][col] = DistanceUpstream;
  }

  cout << "\t Got all the channel elevations! Now calculating the relief..." << endl;

  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (NearestChannelElev_array[row][col] != NoDataValue)
      {
        float this_elev = ElevationRaster.get_data_element(row,col);
        ChannelRelief_array[row][col] = this_elev - NearestChannelElev_array[row][col];
      }
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Takes in a junction number and generates the main stem channel from this point
// Creates 2 vectors for upstream distance and channel relief only for TERRACE pixels which are
// connected to the main channel.
// FJC 26/10/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFloodplain::get_distance_upstream_along_main_stem(int junction_number, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet)
{
  // get the main stem channel from this junction
  LSDIndexChannel MainStem = ChanNetwork.generate_longest_index_channel_in_basin(junction_number, FlowInfo, DistFromOutlet);
  // get the main stem nodes
  MainStemNodes = MainStem.get_NodeSequence();

  vector<int>::iterator find_it;

  // loop through all the nodes and find ones connected to the main stem
  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (NearestChannelNode_array[row][col] != NoDataValue)
      {
        int NearestChanNI = NearestChannelNode_array[row][col];
        //cout << "Channel node: " << NearestChanNI << endl;
        find_it = find(MainStemNodes.begin(), MainStemNodes.end(), NearestChanNI);
        if (find_it != MainStemNodes.end())
        {
          //found a pixel connected to the main stem! Get the distance upstream for this pixel.
          UpstreamDist.push_back(UpstreamDist_array[row][col]);
          ChannelRelief.push_back(ChannelRelief_array[row][col]);
          BinaryMSArray[row][col] = 1;
        }
      }
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to get valley walls - edge of floodplain raster
// FJC 29/01/21
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<LSDRaster> LSDFloodplain::get_valley_walls(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, int threshold_SO)
{
    //set up the arrays
    Array2D<float> ValleyWallArray(NRows,NCols,NoDataValue);
    Array2D<float> ValleyWallSigns_array(NRows,NCols,NoDataValue);

    // for each point in the floodplain, get the nodes of the nearest channel pixels
    Get_Relief_of_Nearest_Channel(ChanNetwork, FlowInfo, ElevationRaster, DistFromOutlet, threshold_SO);

    // create the floodplain raster - this will be used to calculate the distance to the nearest bank
    // in this raster floodplain pixels are set to No Data and the rest of the landscape is set to 1
    LSDRaster Floodplain = get_NoData_FloodplainRaster();

    // get the distance and value of the nearest bank masks for each floodplain pixel
    vector<LSDRaster> nearest_rs = Floodplain.get_nearest_distance_and_value_masks();

    Array2D<float> DistanceArray = nearest_rs[0].get_RasterData();
    int chan_node, ustream_node, dstream_node, ustream_junc, dstream_junc;
    float wall_x, wall_y, ustream_x, ustream_y, dstream_x, dstream_y;
    float det;
    for (int row = 0; row < NRows; row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        if (DistanceArray[row][col] == DataResolution)
        {
          // assign the valley wall to the array
          ValleyWallArray[row][col] = 1;

          // find the nearest channel pixel to that pixel
          chan_node = NearestChannelNode_array[row][col];

          // get the determinant of the vectors between the valley wall and the channel and between the ustream and dstream channel nodes.
          // the determinant is calculated by:
          // det = (wall_x - ustream_x)(dstream_y - ustream_y) - (wall_y - ustream_y)(dstream_x - ustream_x)
          // If det is positive = right bank, negative = left bank (while looking downstream)
          // https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
          // first get x and y locations of the wall pixel
          Floodplain.get_x_and_y_locations(row, col, wall_x, wall_y);
          // now get x and y locations of ustream and dstream channel nodes. We will just get the nearest upstream and downstream junctions to the node
          ustream_junc = ChanNetwork.find_upstream_junction_from_channel_nodeindex(chan_node, FlowInfo);
          dstream_junc = ChanNetwork.get_Receiver_of_Junction(ustream_junc);
          ustream_node = ChanNetwork.get_Node_of_Junction(ustream_junc);
          dstream_node = ChanNetwork.get_Node_of_Junction(dstream_junc);

          FlowInfo.get_x_and_y_from_current_node(ustream_node, ustream_x, ustream_y);
          FlowInfo.get_x_and_y_from_current_node(dstream_node, dstream_x, dstream_y);
          cout << "ustream junc: " << ustream_junc << ", dstream_junc: " << dstream_junc << endl;

          det = (wall_x - ustream_x)*(dstream_y - ustream_y) - (wall_y - ustream_y)*(dstream_x - ustream_x);
          cout << det << endl;
          ValleyWallSigns_array[row][col] = det;
        }
      }
    }

    LSDRaster ValleyWallRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ValleyWallArray, GeoReferencingStrings);
    LSDRaster ValleyWallSigns(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ValleyWallSigns_array, GeoReferencingStrings);
    nearest_rs.push_back(ValleyWallRaster);
    nearest_rs.push_back(ValleyWallSigns);

    return  nearest_rs;
}

//----------------------------------------------------------------------------------------
// Get the valley WIDTHS
// returns a vector of vectors (valley_widths_vec)
// valley_widths_vec[0] = left bank half widths
// valley_widths_vec[1] = right bank half widths
//
// also writes the valley widths to a csv file and geojson
// FJC 05/03/21
//---------------------------------------------------------------------------------------
vector<vector<float>> LSDFloodplain::calculate_valley_widths(LSDSpatialCSVReader channel_csv, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, int channel_bearing_node__spacing, bool print_bearings, string bearing_fname, int template_scale, string outfilename, vector<int> nodes_to_remove)
{
  // open the outfile
  ofstream out_file;
  out_file.open((outfilename+".csv").c_str());
  out_file.precision(9);

  out_file << "node,latitude,longitude,distance_from_outlet,flow_bearing,orthogonal_bearing,left_valley_width,right_valley_width,total_valley_width" << endl;

  // open the geojson for orthogonal lines
  ofstream width_lines;
  width_lines.precision(9);
  width_lines.open((outfilename+"_orthogonal.geojson").c_str());

  width_lines << "{" << endl;
  width_lines << "\"type\": \"FeatureCollection\"," << endl;
  width_lines << "\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\" } }," << endl;
  width_lines << "\"features\": [" << endl;

  vector<float> left_bank_widths, right_bank_widths;
  // get the flow direction bearings along this channel
  vector<int> node_list =  channel_csv.get_nodeindices_from_lat_long(FlowInfo);
  vector<int> bearing_nodes;
  vector<float> bearing_vec;
  FlowInfo.calculate_bearings_from_nodelist(node_list, channel_bearing_node__spacing,
                        print_bearings, bearing_fname,bearing_nodes, bearing_vec);

  cout << "got the bearings" << endl;
  cout << bearing_nodes.size() << " " << bearing_vec.size() << endl;

  // get the distance from outlet array
  Array2D<float> DistFromOutlet_array = DistanceFromOutlet.get_RasterData();

  // loop through each node along the channel
  int this_node, this_row, this_col, left_bank_node, right_bank_node;
  float flow_bearing, orthogonal_bearing, left_bank_width, right_bank_width, total_width;
  LSDCoordinateConverterLLandUTM Converter;
  double left_latitude = NoDataValue;
  double left_longitude = NoDataValue;
  double right_latitude = NoDataValue;
  double right_longitude = NoDataValue;
  int n_nodes = bearing_nodes.size();
  vector<int>::iterator find_it;

  int count = 0;
  for (int i = 0; i < n_nodes; i++)
  {
    this_node = bearing_nodes[i];
    find_it = find(nodes_to_remove.begin(), nodes_to_remove.end(), this_node); // check whether this node is in the list of nodes to remove. If not, then calculate the width.
    if (find_it == nodes_to_remove.end())
    {
      // check if the centre of the template array is floodplain. if not, then the valley width is = Data Resolution. Left bank width and right bank width are both DataRes/2.
      FlowInfo.retrieve_current_row_and_col(this_node, this_row, this_col);
      //cout << this_node << " " << this_row << " " << this_col << endl;
      // now for each node, find the orthogonal to that bearing.
      flow_bearing = deg(bearing_vec[i]);
      //cout << flow_bearing << endl;
      if (flow_bearing < 270)
      {
        orthogonal_bearing = flow_bearing + 90;
      }
      else
      {
        orthogonal_bearing = flow_bearing + 90 - 360;
      }
      //cout << "flow bearing: " << flow_bearing << " orthogonal_bearing: " << orthogonal_bearing << endl;
      if (BinaryArray[this_row][this_col] == 0) // centre node is not in the floodplain.
      {
        left_bank_width = DataResolution/2;
        right_bank_width = DataResolution/2;
        cout << "This node is not on the floodplain, setting the width to the data resolution" << endl;
      }
      else // in this case the centre node is in the floodplain.
      {
        // get a template of pixels along the line of the orthogonal bearing
        Array2D<int> template_array = make_template_for_vector_bearing(orthogonal_bearing, template_scale);

        // loop through the template array and find the pixels with the correct bearing and which are floodplain
        int array_dim = 2*template_scale + 1;
        int row_offset, col_offset, fp_check, bearing_check, check_node;
        int DEM_row, DEM_col;
        float this_dist, this_direction;
        float min_left_dist = 10000000;
        float min_right_dist = 10000000;
        for (int row = 0; row < array_dim; row++)
        {
          for (int col = 0; col < array_dim; col++)
          {
            row_offset = template_scale - row;
            col_offset = template_scale - col;
            DEM_row = this_row - row_offset;
            DEM_col = this_col - col_offset;
            //cout << DEM_row << " " << DEM_col << endl;
            // check if this row or column is off the edge of the DEM
            if (DEM_row < 0 || DEM_row > NRows || DEM_col < 0 || DEM_col > NCols)
            {
              cout << "At the edge of the DEM, setting width to NoDataValue" << endl;
              min_left_dist = 10000000;
              min_right_dist = 10000000;
              goto endloop;
            }
            else
            {
              // check if there is a floodplain pixel here (1) or no floodplain (0)
              fp_check = BinaryArray[DEM_row][DEM_col];
              // check if this pixel lies along the correct bearing (1)
              bearing_check = template_array[row][col];
              //cout << fp_check << " " << bearing_check << endl;
              if (fp_check == 0 && bearing_check == 1)
              {
                // this pixel is not in the floodplain and on the RIGHT side of the bank while looking downstream.
                // now find the distance away from the centre pixel and get the nearest one
                check_node = FlowInfo.retrieve_node_from_row_and_column(DEM_row, DEM_col);
                if (check_node != NoDataValue)
                {
                  this_dist = FlowInfo.get_Euclidian_distance(this_node, check_node);
                  //cout << "check node " << check_node << " this dist " << this_dist << endl;
                  if (this_dist < min_right_dist)
                  {
                    min_right_dist = this_dist;
                    FlowInfo.get_lat_and_long_from_current_node(check_node, right_latitude, right_longitude, Converter);
                    right_bank_node = check_node;
                  }
                }
              }
              if (fp_check == 0 && bearing_check == 2)
              {
                // this pixel is not in the floodplain and on the LEFT side of the bank while looking downstream
                // now find the distance away from the centre pixel and get the nearest one
                check_node = FlowInfo.retrieve_node_from_row_and_column(DEM_row, DEM_col);
                if (check_node != NoDataValue)
                {
                  this_dist = FlowInfo.get_Euclidian_distance(this_node, check_node);
                  //cout << "check node " << check_node << " this dist " << this_dist << endl;
                  if (this_dist < min_left_dist)
                  {
                    min_left_dist = this_dist;
                    FlowInfo.get_lat_and_long_from_current_node(check_node, left_latitude, left_longitude, Converter);
                    left_bank_node = check_node;
                  }
                }
              }
            }
          }
        }
      endloop:
        //cout << "Left half width: " << min_left_dist << ", right half width: " << min_right_dist << endl;
        if (min_left_dist > 1000000 || min_right_dist > 1000000)
        {
          cout << "Could not find the edge of the valley, setting width to NoDataValue" << endl;
          left_bank_width = NoDataValue;
          right_bank_width = NoDataValue;
          total_width = NoDataValue;
        }
        else
        {
          //left_bank_width = min_left_dist - DataResolution/2;
          //right_bank_width = min_right_dist - DataResolution/2;
          left_bank_width = min_left_dist;
          right_bank_width = min_right_dist;
          //total_width = FlowInfo.get_Euclidian_distance(left_bank_node, right_bank_node);
          total_width = left_bank_width + right_bank_width;
        }
      }
      left_bank_widths.push_back(left_bank_width);
      right_bank_widths.push_back(right_bank_width);

      // this is for latitude and longitude
      double latitude, longitude;
      FlowInfo.get_lat_and_long_from_current_node(this_node, latitude, longitude, Converter);
      //cout << "this node: " << this_node << " latitude: " << latitude << " longitude: " << longitude << endl;

      // write to a csv file along the channel
      out_file << this_node << "," << latitude << "," << longitude << "," << DistFromOutlet_array[this_row][this_col] << "," << flow_bearing << "," << orthogonal_bearing << "," << left_bank_width << "," << right_bank_width << "," << total_width << endl;

      if (total_width != NoDataValue && left_latitude != NoDataValue && right_latitude != NoDataValue)
      {
        count+=1;
        // now write the geojson with orthogonal lines
        string first_bit = "{ \"type\": \"Feature\", \"properties\": { \"total_valley_width\": ";
        //string second_bit = dtoa(latitude[i])+", \"longitude\": "+ dtoa(this_longitude);
        string second_bit = ", \"left_valley_width\": ";
        string third_bit = ", \"right_valley_width\": ";
        string fourth_bit = ", \"centreline_node\": ";
        string fifth_bit = ", \"centreline_latitude\": ";
        string sixth_bit = ", \"centreline_longitude\": ";
        string seventh_bit = ", \"flow_bearing\": ";
        string eighth_bit = ", \"orthogonal_bearing\": ";

        string ninth_bit = " }, \"geometry\": { \"type\": \"Linestring\", \"coordinates\": [ [";
        //string fifth_bit = dtoa(this_longitude) +","+ dtoa(latitude[i]) +" ] } },";

        string last_bit = " ] ] } }";
        if (count == 1)
        {
          width_lines << first_bit << total_width << second_bit << left_bank_width
                  << third_bit << right_bank_width << fourth_bit << this_node << fifth_bit << latitude << sixth_bit << longitude << seventh_bit 
                  << flow_bearing << eighth_bit << orthogonal_bearing << ninth_bit << left_longitude << "," << left_latitude << "], [" << right_longitude << "," << right_latitude << last_bit
                  << endl;
        }
        else
        {
          width_lines << "," << first_bit << total_width << second_bit << left_bank_width
                  << third_bit << right_bank_width << fourth_bit << this_node << fifth_bit << latitude << sixth_bit << longitude << seventh_bit 
                  << flow_bearing << eighth_bit << orthogonal_bearing << ninth_bit << left_longitude << "," << left_latitude << "], [" << right_longitude << "," << right_latitude << last_bit
                  << endl;
        }
      }
    }
  }

  out_file.close();

  width_lines << "]" << endl;
  width_lines << "}" << endl;

  width_lines.close();

  cout << "Finished getting valley width, printing to geojson..." << endl;

  // create the vector of vectors
  vector<vector<float>> valley_widths_vec;
  valley_widths_vec.push_back(left_bank_widths);
  valley_widths_vec.push_back(right_bank_widths);

  // print to geojson
  string gjson_name =outfilename+".geojson";
  LSDSpatialCSVReader thiscsv(outfilename+".csv");
  thiscsv.print_data_to_geojson_linestring(gjson_name);

  return valley_widths_vec;

}

//----------------------------------------------------------------------------------------
// Get the valley WIDTHS for multiple channels
// returns a vector of vectors (valley_widths_vec)
// valley_widths_vec[0] = left bank half widths
// valley_widths_vec[1] = right bank half widths
//
// also writes the valley widths to a csv file and geojson
// FJC 05/03/21
//---------------------------------------------------------------------------------------
void LSDFloodplain::calculate_valley_widths_multiple_channels(vector<int> basin_junctions, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNetwork, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea, LSDRaster& Elevation, int channel_bearing_node__spacing, bool print_bearings, string bearing_fname, int template_scale, string outfilename)
{
  // open the outfile
  ofstream out_file;
  out_file.open((outfilename+".csv").c_str());
  out_file.precision(9);

  out_file << "source_id,node,latitude,longitude,distance_from_outlet,drainage_area,elevation,flow_bearing,orthogonal_bearing,left_valley_width,right_valley_width,total_valley_width" << endl;

  // open the geojson for orthogonal lines
  ofstream width_lines;
  width_lines.precision(9);
  width_lines.open((outfilename+"_orthogonal.geojson").c_str());

  width_lines << "{" << endl;
  width_lines << "\"type\": \"FeatureCollection\"," << endl;
  width_lines << "\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\" } }," << endl;
  width_lines << "\"features\": [" << endl;

  //-----------------------------------------------------------------------------------//
  // GETTING CHANNEL SEGMENTS FOR WIDTH CALCULATION
  //-----------------------------------------------------------------------------------//
  cout << "Getting the channel segments for width calculation..." << endl;
  vector<int> source_nodes;
  vector<int> outlet_nodes;
  vector<int> baselevel_node_of_each_basin;
  int n_nodes_to_visit = 10;

  ChanNetwork.get_overlapping_channels(FlowInfo, basin_junctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

  cout << "Got the channel segments! Proceeding to width analysis..." << endl;

  int count = 0;
  for (int chan = 0; chan < int(source_nodes.size()); chan++)
  {
    //cout << "source node: " << source_nodes[chan] << " outlet node: " << outlet_nodes[chan] << endl;
    // get this particular channel (it is a chi network with only one channel)
    LSDIndexChannel ThisChannel(source_nodes[chan], outlet_nodes[chan], FlowInfo);
    vector<int> node_list = ThisChannel.get_NodeSequence();
    //cout << "This channel is " << int(node_list.size()) << " nodes long" << endl;

    //-----------------------------------------------------------------------------------//
    // VALLEY WIDTH CALCULATIONS
    //-----------------------------------------------------------------------------------//

    vector<float> left_bank_widths, right_bank_widths;
    // get the flow direction bearings along this channel
    vector<int> bearing_nodes;
    vector<float> bearing_vec;
    FlowInfo.calculate_bearings_from_nodelist(node_list, channel_bearing_node__spacing,
                          print_bearings, bearing_fname,bearing_nodes, bearing_vec);

    //cout << "got the bearings" << endl;
    //cout << bearing_nodes.size() << " " << bearing_vec.size() << endl;


    // loop through each node along the channel
    int this_node, this_row, this_col, left_bank_node, right_bank_node;
    float flow_bearing, orthogonal_bearing, left_bank_width, right_bank_width, total_width;
    LSDCoordinateConverterLLandUTM Converter;
    double left_latitude = NoDataValue;
    double left_longitude = NoDataValue;
    double right_latitude = NoDataValue;
    double right_longitude = NoDataValue;
    int n_nodes = bearing_nodes.size();

    for (int i = 0; i < n_nodes; i++)
    {
      this_node = bearing_nodes[i];
      //cout << "THIS NODE " << this_node << endl;
      //cout << "Node " << i << " of " << n_nodes << endl;
      // check if the centre of the template array is floodplain. if not, then the valley width is = Data Resolution. Left bank width and right bank width are both DataRes/2.
      FlowInfo.retrieve_current_row_and_col(this_node, this_row, this_col);
      //cout << this_node << " " << this_row << " " << this_col << endl;
      // now for each node, find the orthogonal to that bearing.
      flow_bearing = deg(bearing_vec[i]);
      //cout << flow_bearing << endl;
      if (flow_bearing < 270)
      {
        orthogonal_bearing = flow_bearing + 90;
      }
      else
      {
        orthogonal_bearing = flow_bearing + 90 - 360;
      }
      //cout << "flow bearing: " << flow_bearing << " orthogonal_bearing: " << orthogonal_bearing << endl;
      if (BinaryArray[this_row][this_col] == 0) // centre node is not in the floodplain.
      {
        left_bank_width = DataResolution/2;
        right_bank_width = DataResolution/2;
        //cout << "This node is not on the floodplain, setting the width to the data resolution" << endl;
      }
      else // in this case the centre node is in the floodplain.
      {
        // loop through the template array and find the pixels with the correct bearing and which are floodplain
        int array_dim = 2*template_scale + 1;
        int row_offset, col_offset, fp_check, bearing_check, check_node;
        int DEM_row, DEM_col;
        float this_dist, this_direction;
        float min_left_dist = 10000000;
        float min_right_dist = 10000000;

        // bounds checking for making the template array
        int min_row_offset = this_row - template_scale;
        int max_row_offset = this_row + template_scale;
        int min_col_offset = this_col - template_scale;
        int max_col_offset = this_col + template_scale;
        if (min_row_offset >= 0 || max_row_offset <= NRows || min_col_offset >= 0 || max_col_offset <= NCols)
        {
        
          // get a template of pixels along the line of the orthogonal bearing
          Array2D<int> template_array = make_template_for_vector_bearing(orthogonal_bearing, template_scale);
          //cout << "got the template array" << endl;

          for (int row = 0; row < array_dim; row++)
          {
            for (int col = 0; col < array_dim; col++)
            {
              row_offset = template_scale - row;
              col_offset = template_scale - col;
              DEM_row = this_row - row_offset;
              DEM_col = this_col - col_offset;
              //cout << "CHECK ROW AND COL " << DEM_row << " " << DEM_col << endl;
              // check if this row or column is off the edge of the DEM
              if (DEM_row <= 0 || DEM_row >= NRows || DEM_col <= 0 || DEM_col >= NCols)
              {
                cout << "At the edge of the DEM, setting width to NoDataValue" << endl;
                min_left_dist = 10000000;
                min_right_dist = 10000000;
                goto endloop;
              }
              else
              {
                // check if there is another issue with this row and col
                check_node = FlowInfo.retrieve_node_from_row_and_column(DEM_row, DEM_col);
                if (check_node != NoDataValue)
                {
                  //cout << "LINE 724" << endl;
                  // check if there is a floodplain pixel here (1) or no floodplain (0)
                  fp_check = BinaryArray[DEM_row][DEM_col];
                  //cout << "LINE 727" << endl;
                  // check if this pixel lies along the correct bearing (1)
                  bearing_check = template_array[row][col];
                  //cout << "LINE 730" << endl;
                  //cout << fp_check << " " << bearing_check << endl;
                  if (fp_check == 0 && bearing_check == 1)
                  {
                    // this pixel is not in the floodplain and on the RIGHT side of the bank while looking downstream.
                    // now find the distance away from the centre pixel and get the nearest one
                    //cout << "RIGHT check node " << DEM_row << " " << DEM_col << endl;
                    this_dist = FlowInfo.get_Euclidian_distance(this_node, check_node);
                    //cout << "RIGHT check node " << check_node << " this dist " << this_dist << " min dist " << min_right_dist << endl;
                    if (this_dist < min_right_dist)
                    {
                      min_right_dist = this_dist;
                      FlowInfo.get_lat_and_long_from_current_node(check_node, right_latitude, right_longitude, Converter);
                      //cout << "lat " << right_latitude << endl;
                      right_bank_node = check_node;
                    }
                    //cout << "ENDING THIS BIT" << endl;
                  }
                  if (fp_check == 0 && bearing_check == 2)
                  {
                    // this pixel is not in the floodplain and on the LEFT side of the bank while looking downstream
                    // now find the distance away from the centre pixel and get the nearest one
                    //cout << "LEFT check node " << check_node << endl;
                    this_dist = FlowInfo.get_Euclidian_distance(this_node, check_node);
                    //cout << "LEFT check node " << check_node << " this dist " << this_dist << endl;
                    if (this_dist < min_left_dist)
                    {
                      min_left_dist = this_dist;
                      FlowInfo.get_lat_and_long_from_current_node(check_node, left_latitude, left_longitude, Converter);
                      left_bank_node = check_node;
                    }
                  }
                }
              }
            }
          }
        }
      endloop:
        //cout << "REACHED END OF LOOP" << endl;
        //cout << "Left half width: " << min_left_dist << ", right half width: " << min_right_dist << endl;
        if (min_left_dist > 1000000 || min_right_dist > 1000000)
        {
          cout << "Could not find the edge of the valley, setting width to NoDataValue" << endl;
          left_bank_width = NoDataValue;
          right_bank_width = NoDataValue;
        }
        else
        {
          //left_bank_width = min_left_dist - DataResolution/2;
          //right_bank_width = min_right_dist - DataResolution/2;
          left_bank_width = min_left_dist;
          right_bank_width = min_right_dist;
          //total_width = FlowInfo.get_Euclidian_distance(left_bank_node, right_bank_node);
          count+=1;
          //cout << "L WIDTH " << left_bank_width << " R WIDTH" << right_bank_width << endl;
        }
      }
      if (left_bank_width == NoDataValue || right_bank_width == NoDataValue) { total_width = NoDataValue; }
      else { total_width = left_bank_width + right_bank_width; }
      left_bank_widths.push_back(left_bank_width);
      right_bank_widths.push_back(right_bank_width);

      // this is for latitude and longitude
      double latitude, longitude;
      FlowInfo.get_lat_and_long_from_current_node(this_node, latitude, longitude, Converter);
      //cout << "This node " << this_node << " latitude " << latitude << " longitude " << longitude << endl;

      // write to a csv file along the channel
      out_file << source_nodes[chan] << "," << this_node << "," << latitude << "," << longitude << "," << DistanceFromOutlet.get_data_element(this_row, this_col) 
               << "," << DrainageArea.get_data_element(this_row, this_col) << "," << Elevation.get_data_element(this_row, this_col)
               << "," << flow_bearing << "," << orthogonal_bearing << "," << left_bank_width << "," << right_bank_width << "," << total_width << endl;
      //cout << "wrote line to csv" << endl;

      if (total_width != NoDataValue && left_latitude != NoDataValue && right_latitude != NoDataValue)
      {
        // now write the geojson with orthogonal lines
        string first_bit = "{ \"type\": \"Feature\", \"properties\": { \"total_valley_width\": ";
        //string second_bit = dtoa(latitude[i])+", \"longitude\": "+ dtoa(this_longitude);
        string source_bit = ", \"source_id\": ";
        string dist_bit = ", \"distance_from_outlet\": ";
        string area_bit = ", \"drainage_area\": ";
        string elev_bit = ", \"elevation\": ";
        string second_bit = ", \"left_valley_width\": ";
        string third_bit = ", \"right_valley_width\": ";
        string fourth_bit = ", \"centreline_node\": ";
        string fifth_bit = ", \"centreline_latitude\": ";
        string sixth_bit = ", \"centreline_longitude\": ";
        string seventh_bit = ", \"flow_bearing\": ";
        string eighth_bit = ", \"orthogonal_bearing\": ";

        string ninth_bit = " }, \"geometry\": { \"type\": \"Linestring\", \"coordinates\": [ [";
        //string fifth_bit = dtoa(this_longitude) +","+ dtoa(latitude[i]) +" ] } },";

        string last_bit = " ] ] } }";
        if (count == 1)
        {
          width_lines << first_bit << total_width << source_bit << source_nodes[chan] << dist_bit << DistanceFromOutlet.get_data_element(this_row, this_col) 
                  << area_bit << DrainageArea.get_data_element(this_row, this_col) << elev_bit << Elevation.get_data_element(this_row, this_col)
                  << second_bit << left_bank_width << third_bit << right_bank_width << fourth_bit << this_node << fifth_bit << latitude << sixth_bit << longitude << seventh_bit 
                  << flow_bearing << eighth_bit << orthogonal_bearing << ninth_bit << left_longitude << "," << left_latitude << "], [" << right_longitude << "," << right_latitude << last_bit
                  << endl;
        }
        else
        {
          width_lines << "," << first_bit << total_width << source_bit << source_nodes[chan] << dist_bit << DistanceFromOutlet.get_data_element(this_row, this_col) 
                  << area_bit << DrainageArea.get_data_element(this_row, this_col) << elev_bit << Elevation.get_data_element(this_row, this_col)
                  << second_bit << left_bank_width << third_bit << right_bank_width << fourth_bit << this_node << fifth_bit << latitude << sixth_bit << longitude << seventh_bit 
                  << flow_bearing << eighth_bit << orthogonal_bearing << ninth_bit << left_longitude << "," << left_latitude << "], [" << right_longitude << "," << right_latitude << last_bit
                  << endl;
        }
      }
    }
  }

  out_file.close();

  width_lines << "]" << endl;
  width_lines << "}" << endl;

  width_lines.close();

  cout << "Finished getting valley width, printing to geojson..." << endl;

  // create the vector of vectors
  //vector<vector<float>> valley_widths_vec;
  //valley_widths_vec.push_back(left_bank_widths);
  //valley_widths_vec.push_back(right_bank_widths);

  // print to geojson
  string gjson_name =outfilename+".geojson";
  LSDSpatialCSVReader thiscsv((outfilename+".csv").c_str());
  thiscsv.print_data_to_geojson(gjson_name);

  //return valley_widths_vec;

}

//----------------------------------------------------------------------------------------
// Loop through the valley csv and get a list of nodes that are within search radius of a tributary
// junction greater or equal to the threshold stream order
// FJC 19/03/21
//----------------------------------------------------------------------------------------
vector<int> LSDFloodplain::remove_tributary_nodes(LSDSpatialCSVReader& channel_csv, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, int trib_search_radius, int trib_stream_order)
{
  vector<int> nodes_to_remove;

  vector<int> node_list = channel_csv.get_nodeindices_from_lat_long(FlowInfo);
  int dstream_junc, ustream_junc, ustream_SO, dstream_SO, ustream_node, dstream_node, n_nodes;

  for(int i = 0; i < int(node_list.size()); i++)
  {
    // get the upslope and downslope junction for this node
    ustream_junc = ChanNetwork.find_upstream_junction_from_channel_nodeindex(node_list[i], FlowInfo);
    dstream_junc = ChanNetwork.get_junction_downstream_of_node(node_list[i], FlowInfo);
    ustream_SO = ChanNetwork.get_StreamOrder_of_Junction(ustream_junc);
    dstream_SO = ChanNetwork.get_StreamOrder_of_Junction(dstream_junc);

    if (ustream_SO >= trib_stream_order)
    {
      // check the number of nodes between this node and the upstream junction
      {
        ustream_node = ChanNetwork.get_Node_of_Junction(ustream_junc);
        LSDIndexChannel this_chan(ustream_node, node_list[i], FlowInfo);
        n_nodes = this_chan.get_n_nodes_in_channel();
        if (n_nodes < trib_search_radius)
        {
          nodes_to_remove.push_back(node_list[i]);
        }
      }
    }
    if (dstream_SO >= trib_stream_order)
    {
      // check the number of nodes between this node and the downstream junction
      {
        dstream_node = ChanNetwork.get_Node_of_Junction(dstream_junc);
        LSDIndexChannel this_chan(node_list[i], dstream_node, FlowInfo);
        n_nodes = this_chan.get_n_nodes_in_channel();
        if (n_nodes < trib_search_radius)
        {
          nodes_to_remove.push_back(node_list[i]);
        }
      }
    }
  }

  return nodes_to_remove;
}

//----------------------------------------------------------------------------------------
// Fill holes in the floodplain raster using OpenCV floodfill.
// FJC 29/10/21
//----------------------------------------------------------------------------------------
void LSDFloodplain::fill_holes_in_floodplain()
{
    // use OpenCV floodfill to fill holes in the floodplain raster
    // floodfill from point (0,0)

    // create an OpenCV Mat object from the raster data
    Mat img(NRows, NCols, CV_8UC1);
    for(int i=0; i<NRows; ++i)
         for(int j=0; j<NCols; ++j)
              img.at<unsigned char>(Point(j, i)) = BinaryArray[i][j];

    cout << "got mat object" << endl;
    Mat img_floodfill = img.clone();

    cv::floodFill(img_floodfill, cv::Point(0,0), 1);
    cout << "completed flood filling" << endl;

    // // invert the array to turn 0 to 1 and 1 to 0
    Mat img_floodfill_inv = img_floodfill.clone();
    for(int i=0; i<img_floodfill.rows; ++i)
     for(int j=0; j<img_floodfill.cols; ++j)
          if (img_floodfill.at<unsigned char>(i, j) == 0) { img_floodfill_inv.at<unsigned char>(i, j) = 1; }
          else if (img_floodfill.at<unsigned char>(i, j) == 1) { img_floodfill_inv.at<unsigned char>(i, j) = 0; }

    cout << "inverted array" << endl;

    // // combine the two images to get the filled floodplain
    Mat im_out = img.clone();
    for(int i=0; i<img.rows; i++)
      for(int j=0; j<img.cols; j++)
        im_out.at<unsigned char>(i, j) = img.at<unsigned char>(i, j) + img_floodfill_inv.at<unsigned char>(i, j);
    cout << "combined images to get filled floodplain" << endl;

    // write out image to the binary array
    for(int i=0; i<im_out.rows; ++i)
     for(int j=0; j<im_out.cols; ++j)
          BinaryArray[i][j] = im_out.at<unsigned char>(i, j);
}
//----------------------------------------------------------------------------------------
// FUNCTIONS TO GENERATE RASTERS
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Get the connected components array to a raster
// FJC 20/10/16
//----------------------------------------------------------------------------------------
LSDIndexRaster LSDFloodplain::print_ConnectedComponents_to_Raster()
{
  LSDIndexRaster ConnectedComponents(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ConnectedComponents_Array, GeoReferencingStrings);
  return ConnectedComponents;
}

//----------------------------------------------------------------------------------------
// Print binary raster of floodplain locations
// FJC 24/11/16
//----------------------------------------------------------------------------------------
LSDIndexRaster LSDFloodplain::print_BinaryRaster()
{
  LSDIndexRaster BinaryRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, BinaryArray, GeoReferencingStrings);
  return BinaryRaster;
}

//----------------------------------------------------------------------------------------
// Get raster of floodplain locations where the floodplain is no data and the
// surrounding landscape has the elevation of the topography
// FJC 29/01/21
//----------------------------------------------------------------------------------------
LSDRaster LSDFloodplain::get_NoData_FloodplainRaster(LSDRaster& TopographyRaster)
{
  //Array2D<int> FloodplainArray = get_FloodplainArray();
  Array2D<float> NoDataArray(NRows, NCols, NoDataValue);

  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (BinaryArray[row][col] == 0)
      {
        float this_elev = TopographyRaster.get_data_element(row, col);
        NoDataArray[row][col] = this_elev;
      }
    }
  }
  LSDRaster NoDataRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, NoDataArray, GeoReferencingStrings);
  return NoDataRaster;

}

LSDRaster LSDFloodplain::get_NoData_FloodplainRaster()
{
  //Array2D<int> FloodplainArray = get_FloodplainArray();
  Array2D<float> NoDataArray(NRows, NCols, 1);

  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (BinaryArray[row][col] != 0)
      {
        NoDataArray[row][col] = NoDataValue;
      }
    }
  }
  LSDRaster NoDataRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, NoDataArray, GeoReferencingStrings);
  return NoDataRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to buffer the nodata floodplain by a defined number of pixels. This is to 
// ensure successful extraction of the centreline.
// FJC 06/05/21
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFloodplain::Buffer_NoData_FloodplainRaster(LSDRaster& NoData_FloodplainRaster, int n_pixels)
{
  Array2D<float> RasterData = NoData_FloodplainRaster.get_RasterData();
  Array2D<float> NewRasterData = RasterData.copy();

  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      if (RasterData[row][col] == NoDataValue)
      {
        // loop through pixels in each direction in a spiral and add the buffer
        for (int i = 1; i <= n_pixels; i++)
        {
          for (int j = -i; j <= i; j++) 
          {
            if ((row-i) > 0 && (row+i) < NRows && (col-i) > 0 && (col+i) < NCols) { 
              //cout << "i " << i << " j " << j << endl;
              NewRasterData[row-i][col+j] = NoDataValue;  // top row
              NewRasterData[row+j][col+i] = NoDataValue;  // right column
              NewRasterData[row+i][col+j] = NoDataValue;  // bottom row
              NewRasterData[row+j][col-i] = NoDataValue;  // left column
            }
          }
        } 
      }
    }
  }

  LSDRaster NoDataRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, NewRasterData, GeoReferencingStrings);
  return NoDataRaster;
}

////----------------------------------------------------------------------------------------
//// Get the raster of channel relief relative to the nearest channel reach
//// FJC 18/10/16
////----------------------------------------------------------------------------------------
LSDRaster LSDFloodplain::print_ChannelRelief_to_Raster()
{
  LSDRaster ChannelRelief(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ChannelRelief_array, GeoReferencingStrings);
  return ChannelRelief;
}

////----------------------------------------------------------------------------------------
//// Get the raster of upstream distance along main stem channel
//// FJC 24/10/16
////----------------------------------------------------------------------------------------
LSDRaster LSDFloodplain::print_UpstreamDistance_to_Raster()
{
  LSDRaster UpstreamDistance(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, UpstreamDist_array, GeoReferencingStrings);
  return UpstreamDistance;
}


////----------------------------------------------------------------------------------------
//// FUNCTIONS TO PRINT TEXT FILES
////----------------------------------------------------------------------------------------
////----------------------------------------------------------------------------------------
//// Write a text file with the distance along main stem and channel relief for each
//// CC pixel.  The format is:
//// distance_upstream channel_relief
//// FJC 19/10/16
////----------------------------------------------------------------------------------------
void LSDFloodplain::print_ChannelRelief_to_File(string filename)
{
  ofstream output_file;
  output_file.open(filename.c_str());

  for (int i = 0; i < int(UpstreamDist.size()); i++)
  {
    output_file << UpstreamDist[i] << " " << ChannelRelief[i] << endl;
  }

  output_file.close();
}

////----------------------------------------------------------------------------------------
//// Write a text file with the distance along main stem and channel relief for each
//// CC pixel, binned by distance along main stem.  The user specifies the
//// bin width. The format is:
//// distance_upstream channel_relief
//// FJC 20/10/16
////----------------------------------------------------------------------------------------
void LSDFloodplain::print_Binned_ChannelRelief_to_File(string filename, float& bin_width, float& bin_lower_limit, float& bin_threshold)
{
  // declare vectors for binning
  vector<float> MeanDistances, MeanReliefs, Midpoints_distance, MedianReliefs, StDev_distance, StDev_relief, StErr_distance, StErr_relief;
  vector<int> n_observations;

  cout << "\t Binning, there are " << UpstreamDist.size() << " observations" << endl;

  // bin the data
  bin_data(UpstreamDist, ChannelRelief, bin_width, MeanDistances, MeanReliefs, Midpoints_distance, MedianReliefs, StDev_distance, StDev_relief, StErr_distance, StErr_relief, n_observations, bin_lower_limit, NoDataValue);

  cout << "\t Binned the data" << endl;

  RemoveSmallBins(MeanDistances, MeanReliefs, Midpoints_distance, StDev_distance, StDev_relief, StErr_distance, StErr_relief, n_observations, bin_threshold);

  cout << "\t Removed small bins" << endl;

  // write to file
  ofstream output_file;
  output_file.open(filename.c_str());

  for (int i = 0; i < int(MeanDistances.size()); i++)
  {
    output_file << MeanDistances[i] << " " << StDev_distance[i] << " " << StErr_distance[i] << " " << MeanReliefs[i] << " " << StDev_relief[i] << " " << StErr_relief[i] << endl;
  }

  output_file.close();
}


#endif