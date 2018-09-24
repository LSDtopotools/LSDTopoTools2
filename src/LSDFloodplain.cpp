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
using namespace std;
using namespace TNT;

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
				}
			}
		}
	}
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
	Array2D<int> BinaryArray (NRows,NCols,NoDataValue);

	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			if (FloodplainNodes_array[row][col] != NoDataValue)
			{
				BinaryArray[row][col] = 1;
			}
		}
	}
	LSDIndexRaster BinaryRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, BinaryArray, GeoReferencingStrings);
	return BinaryRaster;
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
