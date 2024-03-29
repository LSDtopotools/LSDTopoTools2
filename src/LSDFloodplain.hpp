//beginning of the LSDFloodplain class

#ifndef LSDFloodplain_H
#define LSDFloodplain_H

#include <vector>
#include <string>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDSpatialCSVReader.hpp"
using namespace std;
using namespace TNT;

///@brief Class to store information about floodplains and terraces...
class LSDFloodplain
{
  public:

  /// @brief create function
  ///
  /// @details This populates the binary array and connected components array for the floodplain
	/// given rasters of channel relief and slope and thresholds for both. Each pixel
  /// must be below the slope and channel relief threshold to be classified as floodplain.
  /// @author FJC
	/// 18/10/16
  LSDFloodplain(LSDRaster& ChannelRelief, LSDRaster& Slope, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, float relief_threshold, float slope_threshold, int min_patch_size, int threshold_SO)
					{ create(ChannelRelief, Slope, ChanNetwork, FlowInfo, relief_threshold, slope_threshold, min_patch_size, threshold_SO); }

	LSDFloodplain(LSDRaster& ExistingRaster)
					{ create (ExistingRaster); }

  LSDFloodplain()             { create(); }

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
	void Get_Relief_of_Nearest_Channel(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, int threshold_SO);

	/// @brief This function gets the information about all the floodplain pixels connected to the main stem channel from a junction
	/// @details Takes a junction number and generates the main stem channel from this point. THe information about each floodplain or terrace pixel is then calculated relative to the main channel.
	/// @param junction_number junction number of interest
	/// @param ChanNetwork LSDJunctionNetwork object
	/// @param FlowInfo LSDFlow info object
	/// @param DistFromOutlet LSDRaster of distances from the outlet
	/// @param ElevationRaster LSDRaster of elevations
	/// @author FJC
	/// @date 26/10/16
	void get_distance_upstream_along_main_stem(int junction_number, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet);

  /// @brief This function gets the information about all the floodplain pixels connected to the main stem channel from starting and ending junctions
	/// @param junction_number junction number of interest
	/// @param ChanNetwork LSDJunctionNetwork object
	/// @param FlowInfo LSDFlow info object
	/// @param DistFromOutlet LSDRaster of distances from the outlet
	/// @param ElevationRaster LSDRaster of elevations
	/// @author FJC
	/// @date 31/01/21
  void isolate_to_selected_channel(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo&FlowInfo, vector<int> nodes);

  /// @brief This function gets the valley centreline as a raster
  /// @author FJC
  /// @date 30/01/21
  //LSDRaster get_valley_centreline(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, vector<int> nodes, int threshold_SO, float window_radius, int neighbourhood_switch);

  /// @brief This function gets the information about the valley walls
	/// @author FJC
	/// @date 29/01/21
  vector<LSDRaster> get_valley_walls(LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, int threshold_SO);

  /// @brief This function gets the continuous valley width along a csv file of nodes, either the channel or the centreline. 
	/// @param channel_csv The csv file of channel nodes
	/// @param FlowInfo LSDFlow info object
	/// @param DistanceFromOutlet LSDRaster of distances from the outlet
	/// @param channel_bearing_node__spacing The spacing between nodes upstream and downstream of the node of interest used to to calculate the flow bearing. A higher spacing will mean a more smoothed bearing calculation.
	/// @param print_bearings Set to true if you want to print a geojson of the channel nodes with their flow bearings for checking
	/// @param bearing_fname The name of the file to print the bearings
	/// @param template_scale The search radius from the main channel to look for the valley banks. If you increase this it will increase the run time, but look for wider valleys.
	/// @param outfilename The name of the file to print the channel widths
	/// @param nodes_to_remove vector of nodes to remove
	/// @author FJC
	/// @date 31/01/21
  vector<vector<float>> calculate_valley_widths(LSDSpatialCSVReader channel_csv, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea, LSDRaster& Elevation,
                                                int channel_bearing_node__spacing, bool print_bearings, string bearing_fname, int template_scale, string outfilename, vector<int> nodes_to_remove);

  void calculate_valley_widths_multiple_channels(vector<int> basin_junctions, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNetwork, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea, LSDRaster& Elevation,
                                                int channel_bearing_node__spacing, bool print_bearings, string bearing_fname, int template_scale, string outfilename);


 	/// @brief This function looks through the channel csv file and gets a list of any nodes that are within a certain number of nodes upstream or downstream of a 
 	/// tributary junction greater or equal to the threshold stream order. It returns a vector of these nodes which will then be removed in the valley width code.
 	/// @param channel_csv The csv file of channel nodes
 	/// @param ChanNetwork the junction network object
	/// @param FlowInfo LSDFlow info object
	/// @param trib_search_radius The number of channel nodes upstream and downstream of a trib junction that will be removed
	/// @param trib_stream_order The threshold stream order for trib junctions to be removed. widths will be removed if they are within the search radius of any tributaries greater or equal to this stream order.
  vector<int> remove_tributary_nodes(LSDSpatialCSVReader& channel_csv, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, int trib_search_radius, int trib_stream_order);



 	/// @brief This function fills holes in the floodplain before running the valley width technique using the OpenCV flood fill algorithm
	/// @author FJC
	/// @date 29/10/21
  void fill_holes_in_floodplain();


  void fill_holes_in_floodplain_spiral(int window_width);

  /// @brief This function gets the valley centreline flow path from a trough raster and a trough FlowInfo.
  /// If it finds the end of a flow path, it first checks to see that none of the downslope neighbours are
  /// in the channel network (if they are then it forces the flow path to continue)
  /// @param start_ni starting node index of the centreline
  /// @param FlowInfo LSDFlowInfo object corresponding to the trough
  /// @param topo_raster raster of the original topography.
  /// @author FJC
  /// @date 17/06/21
 /// vector<int> get_centreline_flow_path(int start_ni, LSDFlowInfo& FlowInfo, LSDRaster& topo_raster);

	/// FUNCTIONS TO GENERATE RASTERS

	/// @brief This function prints the connected components array to a raster
	/// @return ConnectedComponents connected components raster
	/// @author FJC
	/// @date 20/10/16
	LSDIndexRaster print_ConnectedComponents_to_Raster();

	/// @brief This function prints a binary raster of floodplain locations
	/// @return BinaryRaster binary raster
	/// @author FJC
	/// @date 24/11/16
	LSDIndexRaster print_BinaryRaster();

  /// @brief This function gets a raster of floodplain locations where the floodplain is no data and the surrounding landscape has the elevation value of the topography.
	/// @return TopographyRaster 
	/// @author FJC
	/// @date 29/01/21
  LSDRaster get_NoData_FloodplainRaster(LSDRaster& TopographyRaster);
  LSDRaster get_NoData_FloodplainRaster();


  /// @brief Function to buffer the nodata floodplain by a defined number of pixels. This is to ensure successful extraction of the centreline.
  /// @author FJC 
  /// @date 06/05/21
  LSDRaster Buffer_NoData_FloodplainRaster(LSDRaster& NoData_FloodplainRaster, int n_pixels);

	/// @brief This function prints the channel relief compared to the main stem to a raster
	/// @return ChannelRelief LSDRaster of channel relief
	/// @author FJC
	/// @date 18/10/16
	LSDRaster print_ChannelRelief_to_Raster();

	/// @brief This function prints the upstream distance compared of the nearest main stem to a raster
	/// @return UpstreamDist LSDRaster of upstream distance
	/// @author FJC
	/// @date 19/10/16
	LSDRaster print_UpstreamDistance_to_Raster();

	/// @brief This function prints the nearest channel node for every floodplain pixel to a raster
	/// @return LSDRaster of nearest channel nodes
	/// @author FJC
	/// @date 19/01/23
	LSDIndexRaster print_NearestChannelNode_to_Raster();

	/// @brief This function prints the flow lengths to the nearest main stem to a raster
	/// @return Flow Lengths LSDRaster of flow length
	/// @author FJC
	/// @date 19/10/16
	LSDRaster print_FlowLengths_to_Raster();

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

  /// @brief This function returns the array of floodplain nodes
	/// @author FJC
	/// @date 03/05/17
  Array2D<int> get_FloodplainArray() const { return FloodplainNodes_array; }

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

	/// The binary array of floodplain data
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

	/// vectors to separate floodplain and terrace nodes
	/// vector of floodplain nodes
	vector<int> FloodplainNodes;
	Array2D<int> FloodplainNodes_array;

	/// Store information for a specified main stem junction
	/// vector of nodes on the main stem
	vector<int> MainStemNodes;
	/// vector of distance upstream for every pixel connected to the main stem
	vector<float> UpstreamDist;
	/// vector of channel relief for every pixel connected to the main stem
	vector<float> ChannelRelief;
  /// binary array of data on the main stem.
  Array2D<int> BinaryMSArray;

  private:
  void create();
	void create(LSDRaster& ChannelRelief, LSDRaster& Slope, LSDJunctionNetwork& ChanNetwork, LSDFlowInfo& FlowInfo, float relief_threshold, float slope_threshold, int min_patch_size, int threshold_SO);
	void create(LSDRaster& ExistingRaster);
};

#endif