//beginning of the LSDTerrace class

#ifndef LSDTerrace_H
#define LSDTerrace_H

#include <vector>
#include <string>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDFloodplain.hpp"
#include "LSDSwathProfile.hpp"

// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>

using namespace std;
using namespace TNT;

///@brief Class to store information about floodplains and terraces...
class LSDTerrace
{
  public:

  /// @brief create function
  ///
  /// @details This populates the binary array and connected components array for the terraces
	/// given rasters of channel relief and slope and thresholds for both. Each pixel
  /// must be below the slope and channel relief threshold to be classified as a terrace pixel.
  /// @author FJC
	/// 18/10/16
  LSDTerrace(LSDRaster& ChannelRelief, LSDRaster& Slope, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, float relief_threshold, float slope_threshold, int min_patch_size, int threshold_SO, int RemoveChannelThreshold)
					{ create(ChannelRelief, Slope, ChanNetwork, FlowInfo, relief_threshold, slope_threshold, min_patch_size, threshold_SO, RemoveChannelThreshold); }

	/// @return Number of rows as an integer.
  int get_NRows() const        { return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const        { return NCols; }
  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const        { return XMinimum; }
  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const        { return YMinimum; }
  /// @return Data resolution as an integer.
  float get_DataResolution() const        { return DataResolution; }
  /// @return No Data Value as an integer.
  int get_NoDataValue() const        { return NoDataValue; }
	/// @return Georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }

	/// @brief This function gets the elevation of the nearest channel reach to each patch
	/// @details For each pixel this function finds the nearest channel to a patch greater than a
	/// threshold stream order
	/// @details Terraces - calculates the mean elevation of a reach defined by this
	/// channel and gets the elevation of each pixel compared to this reach. Floodplains - gets the
	/// elevation of the nearest channel pixel.
	/// @param ChanNetwork LSDJunctionNetwork object
	/// @param FlowInfo LSDFlow info object
	/// @param ElevationRaster LSDRaster of elevations
	/// @param DistFromOutlet LSDRaster of distance from outlet
	/// @param threshold_SO threshold stream order
	/// @param search_distance search distance for channel reach
	/// @author FJC
	/// @date 21/10/16
void Get_Relief_of_Nearest_Channel(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, int threshold_SO, int search_distance);

	/// @brief This function gets the information about all the floodplain pixels connected to the main stem channel from a junction
	/// @details Takes a junction number and generates the main stem channel from this point. THe information about each floodplain or terrace pixel is then calculated relative to the main channel.
	/// @param junction_number junction number of interest
	/// @param ChanNetwork LSDJunctionNetwork object
	/// @param FlowInfo LSDFlow info object
	/// @param DistFromOutlet LSDRaster of distances from the outlet
	/// @param ElevationRaster LSDRaster of elevations
	/// @author FJC
	/// @date 26/10/16
	void get_terraces_along_main_stem(int junction_number, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet);

  /// @brief This function prints terrace widths along the swath to csv
	/// @param csv_filename name of csv
  /// @param Swath swath object
	/// @author FJC
	/// @date 21/11/17
  void print_TerraceWidths_to_csv(string csv_filename, LSDSwath& Swath);

  /// @brief This function gets the nearest channel node on the baseline to each terrace pixel
  /// @param Swath LSDSwath object
  /// @param BaselineDistance array of distance along baseline of each channel
  /// @param FlowInfo LSDFlowInfo object
  /// @author FJC
  /// @date 30/11/17
  Array2D<int> get_ChannelNodeArray(LSDSwath& Swath, Array2D<float> BaselineDistance, LSDFlowInfo& FlowInfo);

	/// FUNCTIONS TO GENERATE RASTERS

	/// @brief This function prints the connected components array to a raster
	/// @return ConnectedComponents connected components raster
	/// @author FJC
	/// @date 20/10/16
	LSDIndexRaster print_ConnectedComponents_to_Raster();

	/// @brief This function prints the channel relief to a raster
	/// @return ChannelRelief LSDRaster of channel relief
	/// @author FJC
	/// @date 18/10/16
	LSDRaster print_ChannelRelief_to_Raster();

	/// @brief This function prints the channel relief compared to the main stem to a raster
	/// @return MainStemRelief raster of main stem relief
	/// @author FJC
	/// @date 28/10/16
	LSDRaster print_ChannelRelief_to_Raster_MainStem();

	/// @brief This function prints the upstream distance compared of the nearest main stem to a raster
	/// @return UpstreamDist LSDRaster of upstream distance
	/// @author FJC
	/// @date 19/10/16
	LSDRaster print_UpstreamDistance_to_Raster();

	/// @brief This function prints the upstream distance compared of the nearest main stem to a raster
	/// @return UpstreamDist LSDRaster of upstream distance
	/// @author FJC
	/// @date 28/10/16
	LSDRaster print_UpstreamDistance_to_Raster_MainStem();


	/// @brief This function prints the flow lengths to the nearest main stem to a raster
	/// @return Flow Lengths LSDRaster of flow length
	/// @author FJC
	/// @date 19/10/16
	LSDRaster print_FlowLengths_to_Raster();

	/// @brief This function prints a binary raster of terrace locations
	/// @return BinaryRaster binary raster
	/// @author FJC
	/// @date 19/01/17
	LSDIndexRaster print_BinaryRaster();

  /// @brief This function assigns raster values based on terrace locations
  /// @param InputRaster Input raster
  /// @return raster of terrace locations with assigned raster values
  /// @author FJC
  /// @date 02/02/17
  LSDRaster get_Terraces_RasterValues(LSDRaster& InputRaster);

  /// @brief This function combines the floodplains and terraces into one IndexRaster
  /// terraces = 1; floodplains = 2
  /// @param Floodplains LSDFloodplain object
  /// @return index raster of terrace and floodplain locations
  /// @author FJC
  /// @date 03/05/17
  LSDIndexRaster get_combined_terraces_and_floodplains_raster(LSDFloodplain& Floodplains);

	/// FUNCTIONS TO PRINT TEXT FILES

	/// @brief This function prints the upstream distance and channel relief of the floodplain pixels
	/// to a text file
	/// @author FJC
	/// @date 19/10/16
	void print_ChannelRelief_to_File(string filename);

	/// @brief This function prints the binned upstream distance and channel relief of all the CC
	/// pixels to a text file.
	/// @details The format is: mean_distance st_dev_distance st_err_distance mean_relief st_dev_relief st_err_relief
	/// @author FJC
	/// @date 20/10/16
	void print_Binned_ChannelRelief_to_File(string filename, float& bin_width, float& bin_lower_limit, float& bin_threshold);

  ////----------------------------------------------------------------------------------------
  //// Write a csv file with the area of each terrace.
  //// FJC 07/03/17
  ////----------------------------------------------------------------------------------------
  void print_TerraceAreas_to_file(string filename, LSDFlowInfo& FlowInfo);

  /// @brief This function prints information about each terrace pixel to a csv file.
  /// @param csv_filename filename of the csv
  /// @param ElevationRaster elevation raster
  /// @param FlowInfo LSDFlowInfo object
  /// @param Swath LSDSwathProfile object
	/// @author FJC
	/// @date 28/09/17
  void print_TerraceInfo_to_csv(string csv_filename, LSDRaster& ElevationRaster, LSDRaster& ChannelRelief,  LSDFlowInfo& FlowInfo, LSDSwath& Swath);

  protected:

	/// Number of rows
	int NRows;
	/// Number of columns
	int NCols;
	/// X minimum
	float XMinimum;
	/// Y minimum
	float YMinimum;
	/// Data resolution
	float DataResolution;
	/// No data value
	int NoDataValue;
	/// A map of strings for holding georeferencing information
  map<string,string> GeoReferencingStrings;

	/// Relief threshold
	float relief_threshold;
	/// Slope threshold
	float slope_threshold;

	/// The binary array of terrace data
	Array2D<int> BinaryArray;
	/// The array of connected components
	Array2D<int> ConnectedComponents_Array;
	/// Array of nearest channel elevations
	Array2D<float> NearestChannelElev_array;
	/// Array of nearest channel node indices
	Array2D<int> NearestChannelNode_array;
	/// Array of relief to the nearest channel
	Array2D<float> ChannelRelief_array;
	/// Array of distance upstream of nearest channel
	Array2D<float> UpstreamDist_array;

	/// vector of terrace nodes
	vector<int> TerraceNodes;
	/// array of terrace nodes
	Array2D<int> TerraceNodes_array;
  /// vector of terrace IDs, unique values from the connected components array
  vector<int> TerraceIDs;

	/// Store information for a specified main stem junction
	/// vector of nodes on the main stem
	vector<int> MainStemNodes;
	/// vector of distance upstream for every pixel connected to the main stem
	vector<float> UpstreamDist;
	/// vector of channel relief for every pixel connected to the main stem
	vector<float> ChannelRelief;
	/// array of main stem relief
	Array2D<float> MainStemRelief_array;
	/// array of main stem distance
	Array2D<float> MainStemDist_array;

  private:
	void create(LSDRaster& ChannelRelief, LSDRaster& Slope, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, float relief_threshold, float slope_threshold, int min_patch_size, int threshold_SO, int RemoveChannelThreshold);

};

#endif
