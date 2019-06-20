//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChannel
// Land Surface Dynamics Channel
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for retaining information of channels
//
// This is a derivative class of LSDIndexChannel.
//  LSDIndexChannel alone holds pointers to data in
//  LSDFlowInfo and LSDRaster, whereas LSDChannel
//  contains actual data about the channel such as
//  elevation and drainage area.
//
// These two objects are seperated to save on memory overhead
//  during runtime.
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

// Source code for the LSDChannel object.
// This obect is linked to a LSDIndexChannel object.
// The LSDIndexChannel object holds the row, column and node
// index for each channel. The LSDChannel object contains
// additional information such as elevation, drainage area
// and chi (the transformed coordiante for integral analysis of
// channel profiles

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------


#include <vector>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDChannel.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDMostLikelyPartitionsFinder.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDChannel_CPP
#define LSDChannel_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// first create routine
// creates an LSDChannel by copying from an IndexChannel
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(LSDIndexChannel& InChann)
{

  vector<float> empty_vec;
  Elevation = empty_vec;
  Chi = empty_vec;
  DrainageArea = empty_vec;

  StartJunction = InChann.get_StartJunction();
  EndJunction = InChann.get_EndJunction();
  StartNode = InChann.get_StartNode();
  EndNode = InChann.get_EndNode();

  NRows = InChann.get_NRows();
  NCols = InChann.get_NCols();
  XMinimum = InChann.get_XMinimum();
  YMinimum = InChann.get_YMinimum();
  DataResolution = InChann.get_DataResolution();
  NoDataValue = InChann.get_NoDataValue();
  GeoReferencingStrings = InChann.get_GeoReferencingStrings();

  RowSequence =  InChann.get_RowSequence();
  ColSequence =  InChann.get_ColSequence();
  NodeSequence =  InChann.get_NodeSequence();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this calculates all the channel areas, elevations and chi parameters based on
// for a starting node index and ending node index
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJN, int EJN, float downslope_chi,
                             float m_over_n, float A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
{

  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  DataResolution = FlowInfo.get_DataResolution();
  NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();

  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  float pixel_area = DataResolution*DataResolution;



  StartJunction = -1;
  EndJunction = -1;
  StartNode = SJN;
  EndNode = EJN;

  vector<int> RowI;
  vector<int> ColI;
  vector<int> NdI;

  int curr_node = StartNode;

  // push back the data vecotors with the starting node
  int curr_row, curr_col;
  FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
  NdI.push_back(StartNode);
  RowI.push_back(curr_row);
  ColI.push_back(curr_col);

  int receive_node = -99;
  int receive_row, receive_col;

  // loop through receivers until you get to EndNode
  while(curr_node != EndNode)
  {
    FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

    NdI.push_back(receive_node);
    RowI.push_back(receive_row);
    ColI.push_back(receive_col);

    if (receive_node == curr_node)
    {
      EndNode = curr_node;
      cout << "Warning, the channel has come to a baselevel node before it has"
           << endl << "reached the end node" << endl;
    }
    else
    {
      curr_node = receive_node;
    }
  }

  RowSequence = RowI;
  ColSequence = ColI;
  NodeSequence = NdI;

  // get the number of nodes in the channel
  int n_nodes_in_channel = int(NodeSequence.size());

  // the bottom node is at chi of downslope_chi
  // initiate the chi vector
  vector<float> empty_vec;


  vector<float> chi_temp(n_nodes_in_channel,downslope_chi);


  vector<float> elev_temp(n_nodes_in_channel,float(NoDataValue));
  vector<float> area_temp(n_nodes_in_channel,float(NoDataValue));

  // get the first node
  float curr_area;

  curr_node = NodeSequence[n_nodes_in_channel-1];
  curr_row = RowI[n_nodes_in_channel-1];
  curr_col = ColI[n_nodes_in_channel-1];
  curr_area = float(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
  area_temp[n_nodes_in_channel-1] = curr_area;
  elev_temp[n_nodes_in_channel-1] = Elevation_Raster.get_data_element(curr_row,curr_col);

  // now loop up through the channel, adding chi values
  // note, the channel index are arranges with upstream element first, so you need to go through the channel
  // in reverse order
  for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
  {
     //cout << "ChIndex is: " << ChIndex << endl;
    curr_node = NodeSequence[ChIndex];
    FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);
    if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }
    //cout << "dx is: " << dx << endl;

    curr_area = float(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
    area_temp[ChIndex] = curr_area;
    elev_temp[ChIndex] = Elevation_Raster.get_data_element(curr_row,curr_col);
    chi_temp[ChIndex] = dx*(pow( (A_0/curr_area ),
                          m_over_n))
                           + chi_temp[ChIndex+1];
    //cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
    //     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
  }

  Chi = chi_temp;
  Elevation = elev_temp;
  DrainageArea = area_temp;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this calculates all the channel areas, elevations and chi parameters based on
// for a starting node index and ending node index
// Similar to above but you assign the drainage area so you can use dinfinity
// or some other method if you want
//
// SMM 2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJN, int EJN, float downslope_chi,
                             float m_over_n, float A_0, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster, LSDRaster& Drainage_Area_Raster)
{

  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  DataResolution = FlowInfo.get_DataResolution();
  NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();

  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  //float pixel_area = DataResolution*DataResolution;

  StartJunction = -1;
  EndJunction = -1;
  StartNode = SJN;
  EndNode = EJN;

  vector<int> RowI;
  vector<int> ColI;
  vector<int> NdI;

  int curr_node = StartNode;

  // push back the data vecotors with the starting node
  int curr_row, curr_col;
  FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
  NdI.push_back(StartNode);
  RowI.push_back(curr_row);
  ColI.push_back(curr_col);

  int receive_node = -99;
  int receive_row, receive_col;

  // loop through receivers until you get to EndNode
  while(curr_node != EndNode)
  {
    FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

    NdI.push_back(receive_node);
    RowI.push_back(receive_row);
    ColI.push_back(receive_col);

    if (receive_node == curr_node)
    {
      EndNode = curr_node;
      cout << "Warning, the channel has come to a baselevel node before it has"
           << endl << "reached the end node" << endl;
    }
    else
    {
      curr_node = receive_node;
    }
  }

  RowSequence = RowI;
  ColSequence = ColI;
  NodeSequence = NdI;

  // get the number of nodes in the channel
  int n_nodes_in_channel = int(NodeSequence.size());

  // the bottom node is at chi of downslope_chi
  // initiate the chi vector
  vector<float> empty_vec;


  vector<float> chi_temp(n_nodes_in_channel,downslope_chi);


  vector<float> elev_temp(n_nodes_in_channel,float(NoDataValue));
  vector<float> area_temp(n_nodes_in_channel,float(NoDataValue));

  // get the first node
  float curr_area;

  curr_node = NodeSequence[n_nodes_in_channel-1];
  curr_row = RowI[n_nodes_in_channel-1];
  curr_col = ColI[n_nodes_in_channel-1];
  curr_area = Drainage_Area_Raster.get_data_element(curr_row,curr_col);
  area_temp[n_nodes_in_channel-1] = curr_area;
  elev_temp[n_nodes_in_channel-1] = Elevation_Raster.get_data_element(curr_row,curr_col);

  // now loop up through the channel, adding chi values
  // note, the channel index are arranges with upstream element first, so you need to go through the channel
  // in reverse order
  for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
  {
     //cout << "ChIndex is: " << ChIndex << endl;
    curr_node = NodeSequence[ChIndex];
    FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);
    if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }
    //cout << "dx is: " << dx << endl;

    curr_area = Drainage_Area_Raster.get_data_element(curr_row,curr_col);
    area_temp[ChIndex] = curr_area;
    elev_temp[ChIndex] = Elevation_Raster.get_data_element(curr_row,curr_col);
    chi_temp[ChIndex] = dx*(pow( (A_0/curr_area ),
                          m_over_n))
                           + chi_temp[ChIndex+1];
    //cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
    //     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
  }

  Chi = chi_temp;
  Elevation = elev_temp;
  DrainageArea = area_temp;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this calculates all the channel areas, elevations and chi parameters based on
// for a given LSDChannelIndex
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(float downslope_chi,
                             float m_over_n, float A_0, LSDIndexChannel& InChann, LSDFlowInfo& FlowInfo,
                             LSDRaster& Elevation_Raster)
{


  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  DataResolution = FlowInfo.get_DataResolution();
  NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();

  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  float pixel_area = DataResolution*DataResolution;
  //cout << "data res: " << DataResolution << endl;



  StartJunction = InChann.get_StartJunction();
  EndJunction = InChann.get_EndJunction();
  StartNode = InChann.get_StartNode();
  EndNode = InChann.get_EndNode();

  RowSequence =  InChann.get_RowSequence();
  ColSequence =  InChann.get_ColSequence();
  NodeSequence =  InChann.get_NodeSequence();

  int curr_node = StartNode;

  // push back the data vectors with the starting node
  int curr_row, curr_col;

  // get the number of nodes in the channel
  int n_nodes_in_channel = int(NodeSequence.size());

  // the bottom node is at chi of downslope_chi
  // initiate the chi vector
  vector<float> empty_vec;
  vector<float> chi_temp(n_nodes_in_channel,downslope_chi);
  vector<float> elev_temp(n_nodes_in_channel,float(NoDataValue));
  vector<float> area_temp(n_nodes_in_channel,float(NoDataValue));

  // get the first node
  float curr_area;

  //cout << "downslope_chi: " << downslope_chi << endl;

  curr_node = NodeSequence[n_nodes_in_channel-1];
  FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);
  curr_area = float(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
  area_temp[n_nodes_in_channel-1] = curr_area;
  elev_temp[n_nodes_in_channel-1] = Elevation_Raster.get_data_element(curr_row,curr_col);

  // now loop up through the channel, adding chi values
  // note, the channel index are arranges with upstream element first, so you need to go through the channel
  // in reverse order
  for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
  {

    curr_node = NodeSequence[ChIndex];
    FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,
                                             curr_col);

    //cout << "ChIndex is: " << ChIndex << " curr_node: " << curr_node << " row: "
                //                         << curr_row << " curr_col: " << curr_col << endl;

    if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }
    //cout << "dx is: " << dx << endl;

    curr_area = float(FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area;
    area_temp[ChIndex] = curr_area;
    elev_temp[ChIndex] = Elevation_Raster.get_data_element(curr_row,curr_col);
    chi_temp[ChIndex] = dx*(pow( (A_0/curr_area ),
                          m_over_n))
                           + chi_temp[ChIndex+1];
    //cout << "node " << curr_node << " and chi: " << chi_temp[ChIndex]
    //     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
  }

  Chi = chi_temp;
  Elevation = elev_temp;
  DrainageArea = area_temp;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// creates an index channel with just the node index of the starting and ending nodes
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJN, int EJN, LSDFlowInfo& FlowInfo)
{

  vector<float> empty_vec;
  Elevation = empty_vec;
  Chi = empty_vec;
  DrainageArea = empty_vec;

  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  DataResolution = FlowInfo.get_DataResolution();
  NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();

  StartJunction = -1;
  EndJunction = -1;
  StartNode = SJN;
  EndNode = EJN;

  vector<int> RowI;
  vector<int> ColI;
  vector<int> NdI;

  int curr_node = StartNode;

  // push back the data vecotors with the starting node
  int curr_row, curr_col;
  FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
  NdI.push_back(StartNode);
  RowI.push_back(curr_row);
  ColI.push_back(curr_col);

  int receive_node = -99;
  int receive_row, receive_col;

  // loop through receivers until you get to EndNode
  while(curr_node != EndNode)
  {
    FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

    NdI.push_back(receive_node);
    RowI.push_back(receive_row);
    ColI.push_back(receive_col);

    if (receive_node == curr_node)
    {
      EndNode = curr_node;
      cout << "Warning, the channel has come to a baselevel node before it has"
           << endl << "reached the end node" << endl;


    }
    else
    {
      curr_node = receive_node;
    }
  }

  RowSequence = RowI;
  ColSequence = ColI;
  NodeSequence = NdI;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// third create routine
// creates an index channel with just the node index of the starting and ending nodes
// also includes junction information
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::create_LSDC(int SJ, int SJN, int EJ, int EJN, LSDFlowInfo& FlowInfo)
{

  vector<float> empty_vec;
  Elevation = empty_vec;
  Chi = empty_vec;
  DrainageArea = empty_vec;

  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  DataResolution = FlowInfo.get_DataResolution();
  NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();

  StartJunction = SJ;
  EndJunction = EJ;
  StartNode = SJN;
  EndNode = EJN;

  vector<int> RowI;
  vector<int> ColI;
  vector<int> NdI;

  int curr_node = StartNode;

  // push back the data vecotors with the starting node
  int curr_row, curr_col;
  FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
  NdI.push_back(StartNode);
  RowI.push_back(curr_row);
  ColI.push_back(curr_col);

  int receive_node = -99;
  int receive_row, receive_col;

  // loop through receivers until you get to EndNode
  while(curr_node != EndNode)
  {
    FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

    //cout << "receive_node: " << receive_node << " and Endnode: " << EndNode << endl;

    NdI.push_back(receive_node);
    RowI.push_back(receive_row);
    ColI.push_back(receive_col);

    if (receive_node == curr_node)
    {
      EndNode = curr_node;
      cout << "Warning, the channel has come to a baselevel node before it has"
           << endl << "reached the end node" << endl;

    }
    else
    {
      curr_node = receive_node;
    }
  }
  RowSequence = RowI;
  ColSequence = ColI;
  NodeSequence = NdI;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function prints the LSDChannel to an LSDIndexRaster
//
// FJC 21/08/15
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDChannel::print_channel_to_IndexRaster(LSDFlowInfo& FlowInfo)
{

  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  DataResolution = FlowInfo.get_DataResolution();
  NoDataValue = FlowInfo.get_NoDataValue();
  GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();
  Array2D<int> nodes_in_channel(NRows,NCols,NoDataValue);

  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      for (int i = 0; i < int(RowSequence.size()); i++)
      {
        if (RowSequence[i] == row && ColSequence[i] == col)
        {
          nodes_in_channel[row][col] = NodeSequence[i];
        }
      }
    }
  }

  LSDIndexRaster Channel(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, nodes_in_channel, GeoReferencingStrings);
  return Channel;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function uses a flow info object to calculate the chi values in the channel
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::calculate_chi(float downslope_chi, float m_over_n, float A_0, LSDFlowInfo& FlowInfo )
{
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  float pixel_area = DataResolution*DataResolution;
  int curr_node;

  // get the number of nodes in the channel
  int n_nodes_in_channel = int(NodeSequence.size());

  // the bottom node is at chi of downslope_chi
  // initiate the chi vector
  vector<float> empty_vec;
  vector<float> chi_temp(n_nodes_in_channel,downslope_chi);

  // now loop up through the channel, adding chi values
  // note, the channel index are arranges with upstream element first, so you need to go through the channel
  // in reverse order
  for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
  {
     //cout << "ChIndex is: " << ChIndex << endl;
    curr_node = NodeSequence[ChIndex];
    if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }
    //cout << "dx is: " << dx << endl;

    chi_temp[ChIndex] = dx*(pow( (A_0/ (float(
                          FlowInfo.retrieve_contributing_pixels_of_node(curr_node))*pixel_area) ),
                          m_over_n))
                           + chi_temp[ChIndex+1];
    //cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
    //     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
  }

  Chi = chi_temp;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates chi but with a flow accumulation raster
//
// SMM 2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::calculate_chi(float downslope_chi, float m_over_n, float A_0,
                            LSDRaster& FlowAccum, LSDFlowInfo& FlowInfo )
{
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  //float pixel_area = DataResolution*DataResolution;
  int curr_node;
  int this_row,this_col;

  // get the number of nodes in the channel
  int n_nodes_in_channel = int(NodeSequence.size());

  // the bottom node is at chi of downslope_chi
  // initiate the chi vector
  vector<float> empty_vec;
  vector<float> chi_temp(n_nodes_in_channel,downslope_chi);

  // now loop up through the channel, adding chi values
  // note, the channel index are arranges with upstream element first, so you need to go through the channel
  // in reverse order
  for (int ChIndex = n_nodes_in_channel-2; ChIndex>=0; ChIndex--)
  {
     //cout << "ChIndex is: " << ChIndex << endl;
    curr_node = NodeSequence[ChIndex];
    if (FlowInfo.retrieve_flow_length_code_of_node(curr_node) == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }
    //cout << "dx is: " << dx << endl;

    // get the row and columm
    FlowInfo.retrieve_current_row_and_col(curr_node, this_row, this_col);

    chi_temp[ChIndex] = dx*(pow( (A_0/FlowAccum.get_data_element(this_row,this_col) ),
                          m_over_n))
                           + chi_temp[ChIndex+1];
    //cout << "link 0, node " << curr_node << " and chi: " << chi_temp[ChIndex]
    //     << " and chi_temp+1: " << chi_temp[ChIndex+1] << endl;
  }

  Chi = chi_temp;
}


// this function gets the most likely channel segments
//
// SMM 2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::find_most_likeley_segments(int minimum_segment_length, float sigma, int target_nodes,
               vector<float>& b_vec, vector<float>&  m_vec,
               vector<float>&   r2_vec,vector<float>&  DW_vec,
               vector<float>& thinned_chi, vector<float>& thinned_elev,
               vector<float>& fitted_elev, vector<int>& node_ref_thinned,
               vector<int>& these_segment_lengths,
               float& this_MLE, int& this_n_segments, int& n_data_nodes,
               float& this_AIC, float& this_AICc )
{
  // first create a segment finder object
        //cout << "making MLEfinder object, " << endl;
  vector<int> empty_vec;
  node_ref_thinned = empty_vec;

  vector<float> reverse_Chi = Chi;
  reverse(reverse_Chi.begin(), reverse_Chi.end());
  vector<float> reverse_Elevation = Elevation;
  reverse(reverse_Elevation.begin(), reverse_Elevation.end());
  vector<int> this_node_sequence = NodeSequence;
  reverse(this_node_sequence.begin(), this_node_sequence.end());

  LSDMostLikelyPartitionsFinder channel_MLE_finder(minimum_segment_length, reverse_Chi, reverse_Elevation);
  //cout << "got MLE finder object" << endl;
  //cout << "rc size: " << reverse_Chi.size() << " and r_elev size: " << reverse_Elevation.size() << endl;
  //cout << "ns size: " << NodeSequence.size() << " and rns sz: " << this_node_sequence.size() << endl;
  // this needs to be thinned. Get the maximum chi value and then determine dchi
  int n_nodes = reverse_Chi.size();
  float max_chi = reverse_Chi[n_nodes-1];
  float min_chi = reverse_Chi[0];
  //cout << "LSDChannel::find_most_likeley_segments, max_chi: " << max_chi << " and min: " << min_chi << endl;
        //cout << "n_nodes is: " <<  channel_MLE_finder.get_n_nodes() << endl;

  float dchi = (max_chi-min_chi)/float(target_nodes);
  cout << "LSDChannel 533, dchi is: " << dchi << endl;


  // now thin the data, preserving the data (not interpoalting)
  vector<int> node_reference;
  channel_MLE_finder.thin_data_target_dx_preserve_data(dchi, node_reference);
  n_nodes = node_reference.size();
  //cout << "number of nodes in node reference: " << n_nodes << endl;
  for (int i = 0; i< n_nodes; i++)
    {
      //cout << " the node reference is: " << node_reference[i] << endl;
      //cout << " node sequence: " << this_node_sequence[ node_reference[i]] << endl;
      node_ref_thinned.push_back(this_node_sequence[ node_reference[i] ]);
    }

  //cout << "thinned, n_nodes is: " <<  channel_MLE_finder.get_n_nodes() << endl;

  // now create a single sigma value vector
  vector<float> sigma_values;
  sigma_values.push_back(sigma);

  // compute the best fit AIC
  channel_MLE_finder.best_fit_driver_AIC_for_linear_segments(sigma_values);


  channel_MLE_finder.get_data_from_best_fit_lines(0, sigma_values, b_vec, m_vec,
                 r2_vec, DW_vec, fitted_elev,these_segment_lengths,
                   this_MLE, this_n_segments, n_data_nodes, this_AIC, this_AICc);

  thinned_chi = channel_MLE_finder.get_x_data();
  thinned_elev = channel_MLE_finder.get_y_data();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function finds the best fit m/n ratio
// THIS DOES NOT SEEM TO DO ANYTHING SMM 11/2017
// SMM 2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChannel::find_best_fit_m_over_n_with_segments(int n_movern, float d_movern,float start_movern,
                  float downslope_chi, float A_0, LSDFlowInfo& FlowInfo,
                                          int minimum_segment_length, float sigma, float target_nodes )
{
  // now get the details of the best fit
  vector<float> m_vec;
  vector<float> b_vec;
  vector<float> r2_vec;
  vector<float> DW_vec;
  vector<float> fitted_y;
  int n_data_nodes;
  int this_n_segments;
  float this_MLE, this_AIC, this_AICc;
  vector<int> these_segment_lengths;
  vector<float> chi_thinned;
  vector<float> elev_thinned;
  vector<float> elev_fitted;
  vector<int> node_ref_thinned;

  // these are used for storing the best fits in the list of m_over_n
  vector<float> best_m_vec;
  vector<float> best_b_vec;
  vector<float> best_r2_vec;
  vector<float> best_DW_vec;
  vector<float> best_fitted_y;
  //int best_n_data_nodes;
  //int best_this_n_segments;
  //float best_this_MLE;
  //float best_this_AIC, best_this_AICc;
  vector<int> best_these_segment_lengths;
  vector<float> best_chi_thinned;
  vector<float> best_elev_thinned;
  vector<float> best_elev_fitted;
  vector<int> best_node_ref_thinned;

  float min_AICc = 9999;
  //float best_movern = start_movern;
  float m_over_n;
  for(int i = 0; i<n_movern; i++)
  {

    m_over_n = float(i)*d_movern+start_movern;

    // recalculate chi
    calculate_chi(downslope_chi, m_over_n,A_0, FlowInfo );

    find_most_likeley_segments(minimum_segment_length, sigma, target_nodes,
       b_vec, m_vec,r2_vec,DW_vec,
             chi_thinned, elev_thinned,elev_fitted, node_ref_thinned,
       these_segment_lengths, this_MLE, this_n_segments, n_data_nodes,
       this_AIC, this_AICc);

    if (this_AICc < min_AICc)
      {
  best_b_vec =b_vec;
  best_m_vec = m_vec;
  best_r2_vec = r2_vec;
  best_DW_vec = DW_vec;
  best_chi_thinned = chi_thinned;
  best_elev_thinned =  elev_thinned;
  best_elev_fitted = elev_fitted;
  best_node_ref_thinned = node_ref_thinned;
  best_these_segment_lengths = these_segment_lengths;
  //best_this_MLE = this_MLE;
  //best_this_n_segments =  this_n_segments;
  //best_n_data_nodes =  n_data_nodes;
  //best_this_AIC = this_AIC;
  //best_this_AICc = this_AICc;

  min_AICc = this_AICc;
  //best_movern = m_over_n;

  //cout << "best AICc: " << this_AICc << " and m_over_n: " << best_movern << endl;
      }
  }

  // now print the channel profile
  //cout << endl << endl << endl << "best fit m_over_n: " << best_movern << " with AICc: " << min_AICc << endl;
  //int n_nodes = chi_thinned.size();
  //for(int i = 0; i<n_nodes; i++)
  //  {
  //    cout << chi_thinned[i] << " " << elev_thinned[i] << " " << elev_fitted[i] << endl;
  //  }
}

///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/// Calculate channel head locations using chi segment fitting.
///
/// Fitting segments to the chi-elevation data of the main stem.  We assume that the profile
/// is made up of 2 segments in chi-space: a linear channel segment and a non-linear hillslope
/// segment.  We loop through the possible combinations of segment lengths, performing a linear
/// regression to calculate the r^2 and DW of each segment length.  We then calculate a test
/// value: r^2 of the channel segment - ((DW of the hillslope segment - 2)/2).  This value
/// will vary between 0 and 1.  The maximum test_value will give the best fit channel and
/// hillslope segments. Need to get the best fit m_over_n value first.
/// Parameters: min_seg_length_for_channel_heads (length used for fitting segments to the chi-
/// elevation profile, a value of 10 is suggested), A_0, m over n, FlowInfo.
/// Return value: integer with the node index of the channel head location.
/// FC 25/09/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChannel::calculate_channel_heads(int min_seg_length_for_channel_heads, float A_0,
                                            float m_over_n, LSDFlowInfo& FlowInfo)
{
    float downslope_chi = 0;
    calculate_chi(downslope_chi, m_over_n, A_0, FlowInfo);
    vector<float> channel_chi;
    vector<float> hillslope_chi;
    vector<float> channel_elev;
    vector<float> hillslope_elev;
    int end_node = Chi.size();
    float test_value;
    float max_test_value = 0;
    int start_node = 0;
    int node_index = 0;
    float elev_intersection = 0;

    vector<float>::iterator vec_iter_start;
    vector<float>::iterator vec_iter_end;

    // Looping through the combinations of hillslope and channel segment lengths
    for (int hill_seg_length = min_seg_length_for_channel_heads; hill_seg_length <= end_node-min_seg_length_for_channel_heads; hill_seg_length++)
    {
      int chan_seg_length = end_node - hill_seg_length;

      // assigning the chi values of the hillslope segment
      hillslope_chi.resize(hill_seg_length);
      vec_iter_start = Chi.begin()+start_node;
      vec_iter_end = vec_iter_start+hill_seg_length;
      hillslope_chi.assign(vec_iter_start,vec_iter_end);

      // assigning the elevation values of the hillslope segment
      hillslope_elev.resize(hill_seg_length);
      vec_iter_start = Elevation.begin()+start_node;
      vec_iter_end = vec_iter_start+hill_seg_length;
      hillslope_elev.assign(vec_iter_start,vec_iter_end);

      // assigning the chi values of the channel segment
      channel_chi.resize(chan_seg_length);
      vec_iter_start = Chi.begin()+start_node+hill_seg_length;
      vec_iter_end = vec_iter_start+chan_seg_length;
      channel_chi.assign(vec_iter_start,vec_iter_end);

      // assigning the elevation values of the channel segment
      channel_elev.resize(chan_seg_length);
      vec_iter_start = Elevation.begin()+start_node+hill_seg_length;
      vec_iter_end = vec_iter_start+chan_seg_length;
      channel_elev.assign(vec_iter_start,vec_iter_end);

      // performing linear regression on the channel segment
      vector<float> residuals_chan;
      vector<float> results_chan = simple_linear_regression(channel_chi,channel_elev, residuals_chan);

      // performing linear regression on the hillslope segment
      vector<float> residuals_hill;
      vector<float> results_hill = simple_linear_regression(hillslope_chi, hillslope_elev, residuals_hill);

      // calculating the test value
      test_value = results_chan[2] - ((results_hill[3] - 2)/2);

      // looping through test values to find the max test value

      if (test_value > max_test_value)
      {
         max_test_value = test_value;
         elev_intersection = channel_elev.front();
         //chi_intersection = channel_chi.front();
      }

    }

    //Getting the node index of the channel heads
    for (unsigned int i = 0; i < Elevation.size(); i++)
    {
      if (Elevation[i] == elev_intersection)
      {
        node_index =  NodeSequence[i];
      }
    }

    return node_index;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function writes a channel to a CSV file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChannel::write_channel_to_csv(string path, string filename, LSDRaster& flow_dist)
{

  ofstream chan_out;
  string chan_fname = path+filename+"_chan.csv";
  chan_out.open(chan_fname.c_str());

  int NNodes = NodeSequence.size();

  // make sure elevation chi and area vectors exist
  if (int(Elevation.size()) !=  NNodes)
  {
    vector<float> temp(NNodes,0);
    Elevation = temp;
  }
  if (int(Chi.size()) !=  NNodes)
  {
    vector<float> temp(NNodes,0);
    Chi = temp;
  }
  if (int(DrainageArea.size()) !=  NNodes)
  {
    vector<float> temp(NNodes,0);
    DrainageArea = temp;
  }


  int this_row,this_col;
  float x,y;

  int id = 0;

  chan_out << "id,x,y,row,column,distance_from_outlet,elevation,drainge_area,chi" << endl;
  for(int node = 0; node<NNodes; node++)
  {
    id++;
    this_row = RowSequence[node];
    this_col = ColSequence[node];
    x = XMinimum + float(this_col)*DataResolution + 0.5*DataResolution;

    // Slightly different logic for y because the DEM starts from the top corner
    y = YMinimum + float(NRows-this_row)*DataResolution - 0.5*DataResolution;


    chan_out << id << "," << x << "," << y << "," << this_row << ","
             << this_col << "," << flow_dist.get_data_element(this_row,this_col) << ","
             << Elevation[node] <<"," << DrainageArea[node] <<","<<Chi[node] << endl;
  }
  chan_out.close();
}

//----------------------------------------------------------------------------------------//
// calculate slopes along the channel using linear regression
// FJC 19/06/19
// completely forgot how to code in c++ so this is probably awful
//----------------------------------------------------------------------------------------//
vector<float> LSDChannel::calculate_channel_slopes(int window_size, LSDRaster& flow_distance)
{
  vector<float> channel_slopes;
  int NNodes = NodeSequence.size();

  // get flow distances along the channel
  vector<float> flow_distances;
  int this_row, this_col;
  float this_dist;
  for (int i = 0; i < NNodes; i++)
  {
    this_row = RowSequence[i];
    this_col = ColSequence[i];
    this_dist = flow_distance.get_data_element(this_row,this_col);
    flow_distances.push_back(this_dist);
  }
  //cout << "Got flow distances" << endl;

  int slicer = (window_size - 1)/2;

  for (int i = 0; i < NNodes; i++)
  {
    int start_idx = i - slicer;
    if (start_idx < 0) start_idx = 0;

    int end_idx = i+slicer;
    if (end_idx >= NNodes) end_idx = NNodes-1;
    //cout << "start idx: " << start_idx << " end idx: " << end_idx << endl;
    //cout << "n elevation nodes: " << Elevation.size() << endl;

    // find the rows above and below relating to the window size. We use whatever nodes
    // are available to not waste the data.
    vector<float>::iterator first = Elevation.begin() + start_idx;
    vector<float>::iterator last = Elevation.begin() + end_idx;
    vector<float> these_elevs(first, last);
    //cout << "Got elev vector" << endl;

    first = flow_distances.begin() + start_idx;
    last = flow_distances.begin() + end_idx;
    vector<float> these_dists(first,last);
    //cout << "Got the vector slices" << endl;
    // now regress this slice
    vector<float> residuals;
    vector<float> regress = simple_linear_regression(these_dists, these_elevs, residuals);
    float slope = regress[0];
    //cout << "slope: " << slope << endl;
    channel_slopes.push_back(abs(slope));
  }

  return channel_slopes;
}


#endif
