//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDTerrace.cpp
//
// Land Surface Dynamics Terrace Object
//
// This object creates and stores information about extracted terraces
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
#ifndef LSDTerrace_CPP
#define LSDTerrace_CPP

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
#include "LSDSwathProfile.hpp"
#include "LSDFloodplain.hpp"
#include "LSDTerrace.hpp"
using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// create
// this populates the binary array and connected components array for the terraces
// given rasters of channel relief and slope and thresholds for both. Each pixel
// must be below the slope and channel relief threshold to be classified as a terrace.
// User must set a minimum patch size (in pixels, set to 0 if all patches are kept).
// Any patches connected to the channel network are removed as these should represent
// the modern floodplain
// FJC 18/10/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDTerrace::create(LSDRaster& ChannelRelief, LSDRaster& Slope, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, float relief_thresh, float slope_thresh, int min_patch_size, int threshold_SO, int RemoveChannelThreshold)
{

  /// set the protected variables
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
	GeoReferencingStrings = FlowInfo.get_GeoReferencingStrings();

	relief_threshold = relief_thresh;
	slope_threshold = slope_thresh;

	//declare the arrays
	Array2D<int> TempBinArray(NRows,NCols,0);
	Array2D<int> TempLinkArray(NRows,NCols,NoDataValue);
	BinaryArray = TempBinArray.copy();
	ConnectedComponents_Array = TempLinkArray.copy();
	Array2D<int> StreamOrderArray = ChanNetwork.get_StreamOrderArray();
	TerraceNodes_array = TempLinkArray.copy();

	//declare the vectors
	vector<int> TerraceNodes_temp, patch_ids_channel;

	//loop through every row and col and get the slope and relief values
  for (int i =0; i < NRows; i++)
  {
    for (int j = 0; j < NCols; j++)
    {
      if (ChannelRelief.get_data_element(i,j) != NoDataValue && Slope.get_data_element(i,j) != NoDataValue)
      {
        float slope = Slope.get_data_element(i,j);
        float relief = ChannelRelief.get_data_element(i,j);
				//terraces must fall within the relief and slope thresholds
				// if (relief < relief_threshold && relief > RemoveChannelThreshold && slope < slope_threshold && StreamOrderArray[i][j] < 3)
				if (relief < relief_threshold && relief > RemoveChannelThreshold && slope < slope_threshold)
        {
          BinaryArray[i][j] = 1;
        }
      }
    }
  }
	LSDIndexRaster TerraceRaster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,BinaryArray,GeoReferencingStrings);

	// run the connected components algorithm on the terrace array
	LSDIndexRaster ConnectedComponents = TerraceRaster.ConnectedComponents();
	if (min_patch_size > 0)
	{
		LSDIndexRaster ConnectedComponents_final = ConnectedComponents.RemoveSmallPatches(min_patch_size);
		ConnectedComponents_Array = ConnectedComponents_final.get_RasterData();
	}
	else
	{
		ConnectedComponents_Array = ConnectedComponents.get_RasterData();
	}

	// push back the terrace IDs to vector
	vector<int> TerraceIDs_temp = Unique(ConnectedComponents_Array, NoDataValue);

	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			if (ConnectedComponents_Array[row][col] != NoDataValue)
			{
				int ThisNode = FlowInfo.retrieve_node_from_row_and_column(row, col);
				TerraceNodes_temp.push_back(ThisNode);
				TerraceNodes_array[row][col] = ThisNode;
			}
		}
	}

	//copy to vector
	TerraceNodes = TerraceNodes_temp;
	TerraceIDs = TerraceIDs_temp;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function finds the nearest channel greater than a threshold stream order
// to each patch.
// Terraces - It calculates the mean elevation of a reach defined by this channel
// and then gets the elevation of each pixel compared to this reach.
// Floodplains - just finds elevation of the nearest channel > threshold SO
// FJC 21/10/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDTerrace::Get_Relief_of_Nearest_Channel(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, int threshold_SO, int search_distance)
{
	//set up the arrays
	Array2D<int> TempIntArray(NRows,NCols,NoDataValue);
	Array2D<float> TempFloatArray(NRows,NCols,NoDataValue);
	NearestChannelElev_array = TempFloatArray.copy();
	ChannelRelief_array = TempFloatArray.copy();
	UpstreamDist_array = TempFloatArray.copy();
	NearestChannelNode_array = TempIntArray.copy();


	//set up the vectors
	vector<int> PatchIDs_vector, ChannelNodes_temp, ChannelNodes_final, Elevation_vector;
	vector<float> FlowLengths;


	cout << "\t Getting elevations of nearest channels for the terraces" << endl;

  // START WITH THE TERRACES
	for (int i = 0; i < int(TerraceNodes.size()); i++)
  {
		int row, col, ChannelNode;
		FlowInfo.retrieve_current_row_and_col(TerraceNodes[i], row, col);
		float FlowLength, DistanceUpstream;
		ChanNetwork.get_info_nearest_channel_to_node(TerraceNodes[i], threshold_SO, FlowInfo, DistFromOutlet, ChannelNode, FlowLength, DistanceUpstream);
		//cout << "Channel node: " << ChannelNode << " Flow length: " << FlowLength << endl;
		int patch_id = ConnectedComponents_Array[row][col];
		PatchIDs_vector.push_back(patch_id);
		FlowLengths.push_back(FlowLength);
		ChannelNodes_temp.push_back(ChannelNode);
		UpstreamDist_array[row][col] = DistanceUpstream;
	}

	cout << "\t Got the flow lengths and nodes for each patch, now finding nearest channel..." << endl;

	// Find the nearest channel node for each patch ID
	 //get unique patch IDs
  vector<int> Unique_Patches = Unique(PatchIDs_vector);

	for (int i =0; i < int(Unique_Patches.size()); i++)
	{
		float ShortestLength = 100000000000;		//arbitrary large number
		int NearestChannel = 0;
		//float UpstreamDist = 0;
		for (int j = 0; j < int (PatchIDs_vector.size()); j++)
		{
			// find the nearest channel node
			if (PatchIDs_vector[j] == Unique_Patches[i])
			{
				if (FlowLengths[j] < ShortestLength)
				{
					//update the flow length and node
					ShortestLength = FlowLengths[j];
					NearestChannel = ChannelNodes_temp[j];
					//UpstreamDist = DistUpstream_temp[j];
				}
			}
		}
		//cout << "Length: " << ShortestLength << " Channel node: " << NearestChannel << " Distance upstream: " << UpstreamDist << endl;

		// Get the average elevation of the reach for this patchID
		float mean_elev = ChanNetwork.find_mean_elevation_of_channel_reach(NearestChannel,search_distance,ElevationRaster,FlowInfo);
		if (mean_elev >= 0)
		{
			Elevation_vector.push_back(mean_elev);
		}
		else { Elevation_vector.push_back(0.0); }
		//cout << "Patch ID: " << Unique_Patches[i] << " Elevation of channel: " << mean_elev << endl;
		ChannelNodes_final.push_back(NearestChannel);
		//DistUpstream_final.push_back(UpstreamDist);
	}

	//cout << "Elev vector: " << Elevation_vector.size() << " Channel nodes vector: " << ChannelNodes_final.size() << " Dist upstream vector: " << DistUpstream_final.size() << endl;

	vector<int>::iterator it;
	// Get the relief relative to channel for each terrace node
	for (int i = 0; i < int(TerraceNodes.size()); i++)
	{
		int row,col;
		FlowInfo.retrieve_current_row_and_col(TerraceNodes[i], row, col);
		it = find(Unique_Patches.begin(), Unique_Patches.end(), PatchIDs_vector[i]);
		int index = it - Unique_Patches.begin();
		//cout << "Index: " << index << endl;
		// update the vector with the elevation for this patch ID
		NearestChannelElev_array[row][col] = Elevation_vector[index];
		NearestChannelNode_array[row][col] = ChannelNodes_final[index];
	}

	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			if (TerraceNodes_array[row][col] != NoDataValue)
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
void LSDTerrace::get_terraces_along_main_stem(int junction_number, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet)
{
	Array2D<float> TempFloatArray(NRows,NCols,NoDataValue);
	MainStemRelief_array = TempFloatArray.copy();
	MainStemDist_array = TempFloatArray.copy();

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
			if (NearestChannelNode_array[row][col] != NoDataValue && TerraceNodes_array[row][col] != NoDataValue)
			{
				int NearestChanNI = NearestChannelNode_array[row][col];
				//cout << "Channel node: " << NearestChanNI << endl;
				find_it = find(MainStemNodes.begin(), MainStemNodes.end(), NearestChanNI);
				if (find_it != MainStemNodes.end())
				{
					//found a pixel connected to the main stem! Get the distance upstream for this pixel.
					UpstreamDist.push_back(UpstreamDist_array[row][col]);
					ChannelRelief.push_back(ChannelRelief_array[row][col]);
					MainStemRelief_array[row][col] = ChannelRelief_array[row][col];
					MainStemDist_array[row][col] = UpstreamDist_array[row][col];
				}
			}
		}
	}
}

//----------------------------------------------------------------------------------------
// This function compiles the data along each terrace into a single bin for each point
// along the main stem. It returns a vector of vectors with
//----------------------------------------------------------------------------------------
// vector< vector<float> > LSDTerrace::Collate_TerraceData_Along_MainStem()
// {
//
// }

//----------------------------------------------------------------------------------------
// FUNCTIONS TO GENERATE RASTERS
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Get the connected components array to a raster
// FJC 20/10/16
//----------------------------------------------------------------------------------------
LSDIndexRaster LSDTerrace::print_ConnectedComponents_to_Raster()
{
	LSDIndexRaster ConnectedComponents(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ConnectedComponents_Array, GeoReferencingStrings);
	return ConnectedComponents;
}

////----------------------------------------------------------------------------------------
//// Get the raster of channel relief relative to the nearest channel reach
//// FJC 18/10/16
////----------------------------------------------------------------------------------------
LSDRaster LSDTerrace::print_ChannelRelief_to_Raster()
{
	LSDRaster ChannelRelief(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ChannelRelief_array, GeoReferencingStrings);
	return ChannelRelief;
}

////----------------------------------------------------------------------------------------
//// Get the raster of channel relief for main stem terraces only
//// FJC 28/10/16
////----------------------------------------------------------------------------------------
LSDRaster LSDTerrace::print_ChannelRelief_to_Raster_MainStem()
{
	LSDRaster MainStemRelief(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, MainStemRelief_array, GeoReferencingStrings);
	return MainStemRelief;
}

////----------------------------------------------------------------------------------------
//// Get the raster of upstream distance along main stem channel
//// FJC 24/10/16
////----------------------------------------------------------------------------------------
LSDRaster LSDTerrace::print_UpstreamDistance_to_Raster()
{
	LSDRaster UpstreamDistance(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, UpstreamDist_array, GeoReferencingStrings);
	return UpstreamDistance;
}


////----------------------------------------------------------------------------------------
//// Get the raster of upstream distance for only main stem terraces
//// FJC 28/10/16
////----------------------------------------------------------------------------------------
LSDRaster LSDTerrace::print_UpstreamDistance_to_Raster_MainStem()
{
	LSDRaster UpstreamDistance(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, MainStemDist_array, GeoReferencingStrings);
	return UpstreamDistance;
}

//----------------------------------------------------------------------------------------
// Print binary raster of terrace locations
// FJC 19/01/17
//----------------------------------------------------------------------------------------
LSDIndexRaster LSDTerrace::print_BinaryRaster()
{
	Array2D<int> BinaryArray (NRows,NCols,NoDataValue);

	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			if (ConnectedComponents_Array[row][col] != NoDataValue)
			{
				BinaryArray[row][col] = 1;
			}
		}
	}
	LSDIndexRaster BinaryRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, BinaryArray, GeoReferencingStrings);
	return BinaryRaster;
}

//----------------------------------------------------------------------------------------
// Print raster of terrace locations with values assigned by an input raster
// FJC 02/02/17
//----------------------------------------------------------------------------------------
LSDRaster LSDTerrace::get_Terraces_RasterValues(LSDRaster& InputRaster)
{
	Array2D<float> Terrace_RasterValues(NRows,NCols,NoDataValue);
	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			float RasterValue = InputRaster.get_data_element(row,col);
			if ((ConnectedComponents_Array[row][col] != NoDataValue) && (RasterValue != NoDataValue))
			{
				Terrace_RasterValues[row][col] = RasterValue;
			}
		}
	}

	LSDRaster TerraceRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Terrace_RasterValues, GeoReferencingStrings);
	return TerraceRaster;
}

//----------------------------------------------------------------------------------------
// Print index raster of terrace and floodplains
// Terraces = 1
// Floodplains = 2
// FJC 03/05/17
//----------------------------------------------------------------------------------------
LSDIndexRaster LSDTerrace::get_combined_terraces_and_floodplains_raster(LSDFloodplain& Floodplains)
{
	Array2D<int> TerraceFloodplain_Array(NRows,NCols,NoDataValue);
	Array2D<int> FloodplainNodes_array = Floodplains.get_FloodplainArray();

	for (int i = 0; i < NRows; i++)
	{
		for (int j = 0; j < NCols; j++)
		{
			if (TerraceNodes_array[i][j] != NoDataValue)
			{
				TerraceFloodplain_Array[i][j] = 1;
			}
			if (FloodplainNodes_array[i][j] != NoDataValue)
			{
				TerraceFloodplain_Array[i][j] = 2;
			}
		}
	}

	LSDIndexRaster TerraceFloodplainRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, TerraceFloodplain_Array, GeoReferencingStrings);
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
void LSDTerrace::print_ChannelRelief_to_File(string filename)
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
void LSDTerrace::print_Binned_ChannelRelief_to_File(string filename, float& bin_width, float& bin_lower_limit, float& bin_threshold)
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

////----------------------------------------------------------------------------------------
//// Write a csv file with the area of each terrace.
//// FJC 07/03/17
////----------------------------------------------------------------------------------------
void LSDTerrace::print_TerraceAreas_to_file(string filename, LSDFlowInfo& FlowInfo)
{
	ofstream output_file;
	output_file.open(filename.c_str());
	output_file << "TerraceID,Area(m)" << endl;
	for (int i = 0; i < int(TerraceIDs.size()); i++)
	{
		//get the n pixels of each CC value
		int this_CC = TerraceIDs[i];
		int n_pixels = 0;
		for (int j = 0; j < int(TerraceNodes.size()); j++)
		{
			int row, col;
			FlowInfo.retrieve_current_row_and_col(TerraceNodes[j], row, col);
			if (ConnectedComponents_Array[row][col] == this_CC)
			{
				n_pixels++;
			}
		}
		// get the area of this terrace - n pixels * DataRes^2
		float TerraceArea = n_pixels*DataResolution*DataResolution;
		output_file << this_CC << "," << TerraceArea << endl;
	}
	output_file.close();
}

//----------------------------------------------------------------------------------------
// Get the nearest channel node on the baseline channel to each terrace pixel
// FJC 30/11/17
//----------------------------------------------------------------------------------------
Array2D<int> LSDTerrace::get_ChannelNodeArray(LSDSwath& Swath, Array2D<float> BaselineDistanceArray, LSDFlowInfo& FlowInfo)
{
	Array2D<int> BaselineNodes(NRows,NCols,NoDataValue);

	// get the baseline values, rows, and cols
	vector<int> BaselineRows = Swath.get_BaselineRows();
	vector<int> BaselineCols = Swath.get_BaselineCols();
	vector<float> BaselineValue = Swath.get_BaselineValue();
	vector<float> BaselineDistance = Swath.get_DistanceAlongBaseline();

	cout << "N baseline rows: " << BaselineRows.size() << " N baseline cols: " << BaselineCols.size() << " baseline distance vector length: " << BaselineDistance.size() << endl;

	// loop through the terrace nodes, and find their position on the channel.
	for (int i = 0; i < NRows; i++)
	{
		for (int j = 0; j < NCols; j++)
		{
			// find the baseline dist of this node
			if (ConnectedComponents_Array[i][j] != NoDataValue)
			{
				float this_dist = BaselineDistanceArray[i][j];
				//cout << this_dist << endl;
				vector<float>::iterator find_it;
				find_it = find(BaselineDistance.begin(), BaselineDistance.end(), this_dist);
				if (find_it != BaselineDistance.end())
				{
					int idx = distance(BaselineDistance.begin(), find_it);
					int BaselineNode = FlowInfo.retrieve_node_from_row_and_column(BaselineRows[idx], BaselineCols[idx]);
					//cout << "baseline node: " << BaselineNode << endl;
					BaselineNodes[i][j] = BaselineNode;
				}
			}
		}
	}

	return BaselineNodes;
}

////----------------------------------------------------------------------------------------
//// Write a csv file giving elevation and distance information for each pixel in each terrace.
//// FJC 28/09/17
////----------------------------------------------------------------------------------------
void LSDTerrace::print_TerraceInfo_to_csv(string csv_filename, LSDRaster& ElevationRaster, LSDRaster& ChannelRelief, LSDFlowInfo& FlowInfo, LSDSwath& Swath)
{
	ofstream output_file;
	output_file.open(csv_filename.c_str());
  output_file.precision(8);

	if (!output_file)
 {
		 cout << "\n Error opening output csv file. Please check your filename";
		 exit(1);
 }
 cout << "Opened the csv" << endl;

	output_file << "TerraceID,Latitude,Longitude,X,Y,Elevation,DistAlongBaseline,DistToBaseline,BaselineNode,ChannelRelief" << endl;

	LSDIndexRaster ConnectedComponents(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,ConnectedComponents_Array,GeoReferencingStrings);

	cout << "Got the CC raster" << endl;

	// get the baseline distance array
	Array2D<float> BaselineDistance = Swath.get_BaselineDist_ConnectedComponents(ConnectedComponents);
	Array2D<float> ElevationArray = ElevationRaster.get_RasterData();
	Array2D<float> ReliefArray = ChannelRelief.get_RasterData();
	Array2D<float> DistToBaseline = Swath.get_DistanceToBaseline_ConnectedComponents(ConnectedComponents);
	Array2D<int> ChannelNodes = get_ChannelNodeArray(Swath, BaselineDistance, FlowInfo);

	cout << "Now writing the terrace information to the csv file..." << endl;

  // the x and y locations
	double x_loc, y_loc;
	double latitude,longitude;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

	// loop through all the rows and cols and print some information
	for (int row=0; row<NRows; row++)
  {
    for (int col=0; col<NCols; col++)
    {
			if (ConnectedComponents_Array[row][col] != NoDataValue && BaselineDistance[row][col] != NoDataValue && ElevationArray[row][col] != NoDataValue && ReliefArray[row][col] != NoDataValue && DistToBaseline[row][col] != NoDataValue)
			{
				// get the latitude and longitude of the point
				ElevationRaster.get_x_and_y_locations(row, col, x_loc, y_loc);
				//cout << "Row: " << row << " Col: " << col << " X: " << x_loc << " Y: " << y_loc << endl;
				ElevationRaster.get_lat_and_long_locations(row, col, latitude, longitude, Converter);

				float this_elev = ElevationRaster.get_data_element(row,col);

				output_file << ConnectedComponents_Array[row][col] << "," << latitude << "," << longitude << "," << x_loc << "," << y_loc << "," << this_elev << "," << BaselineDistance[row][col] << "," << DistToBaseline[row][col] << "," << ChannelNodes[row][col] << "," << ReliefArray[row][col] << endl;
			}
		}
	}
	output_file.close();
}
//-----------------------------------------------------------------------//
// function to print terrace widths to csv
// FJC 21/11/17
//-----------------------------------------------------------------------//

void LSDTerrace::print_TerraceWidths_to_csv(string csv_filename, LSDSwath& Swath)
{
	// get the cc to index raster
	LSDIndexRaster ConnectedComponents = print_ConnectedComponents_to_Raster();

	// now get the terrace widths
	vector<float> Widths = Swath.get_widths_along_swath(ConnectedComponents);

	vector<float> DistAlongBaseline = Swath.get_DistanceAlongBaseline();

	// open the csv
	ofstream output_file;
	output_file.open(csv_filename.c_str());
  output_file.precision(8);

	if (!output_file)
	{
	  cout << "\n Error opening output csv file. Please check your filename";
	  exit(1);
	}
  cout << "Opened the csv" << endl;

	// write to file
	output_file << "DistAlongBaseline,TerraceWidth" << endl;

	for (int i = 0; i < int(DistAlongBaseline.size()); i++)
	{
		output_file << DistAlongBaseline[i] << "," << Widths[i] << endl;
	}
	output_file.close();

}

#endif
