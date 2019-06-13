//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// lsdtt-hillslope-channel-coupling.cpp
//
// This programme combines the chi analysis of Mudd et al. 2014 and the hillslope analysis
// of Hurst et al 2012 and Grieve et al. 2016 to allow combined topographic analysis comparing hillslope and 
// channel metrics for erosion rate in order to explore landscape morphology and transience.
//
// This program was first deployed in Hurst et al., 2019 EPSL
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate
// different kinds of chi analysis
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2019 Martin D. Hurst and Simon M. Mudd 2019
//
// Developers can be contacted:
//
//    martin.hurst _at_ ed.ac.uk
//    Martin D. Hurst
//    School of
//    Room 404, East Quad
//    University of Glasgow
//    Glasgow, G12 8QQ
//    Scotland
//    United Kindom
//
//    simon.m.mudd _at_ ed.ac.uk
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDBasin.hpp"
#include "../LSDChiTools.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the coupled channel and hillslope tool!  ||" << endl;
  cout << "|| This program has a number of options to make chi    ||" << endl;
  cout << "|| plots and to map out slopes in chi space.           ||" << endl;
  cout << "|| This program was developed by Martin D. Hurst       ||" << endl;
  cout << "||  at the University of Glasgow                       ||" << endl;
  cout << "|| and Simon M. Mudd                                   ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;    
  cout << "=========================================================" << endl;   
  cout << "|| If you use these routines please cite:              ||" << endl;   
  cout << "|| *Hurst et al., Paper in press to be updated*        ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| Documentation can be found at:                      ||" << endl;
  cout << "|| https://lsdtopotools.github.io/LSDTT_documentation/ ||" << endl;
  cout << "=========================================================" << endl;

  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);


  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // Basic DEM preprocessing
  float_default_map["minimum_elevation"] = 0.0;
  float_default_map["maximum_elevation"] = 30000;
  float_default_map["min_slope_for_fill"] = 0.0001;
  float_default_map["WindowRadius"] = 12.0;  // This is for the surface mapping (in metres)
  bool_default_map["raster_is_filled"] = false; // assume base raster is already filled
  bool_default_map["remove_seas"] = false; // elevations above minimum and maximum will be changed to nodata
  bool_default_map["only_check_parameters"] = false;

   //Defining hilltops
	int_default_map["StreamNetworkPadding"] = 0;
	int_default_map["min_stream_order_to_extract_basins"] = 0;
	int_default_map["max_stream_order_to_extract_basins"] = 100;
	bool_default_map["RemovePositiveHilltops"] = true;
	bool_default_map["RemoveSteepHilltops"] = true;
	float_default_map["Threshold_Hilltop_Gradient"] = 0.4;
	bool_default_map["MaskHilltopstoBasins"] = true;

  // Input filenames
  string_default_map["ChannelSegments_file"] = "NULL";
	string_default_map["Floodplain_file"] = "NULL";  
  string_default_map["CHeads_file"] = "NULL";
  string_default_map["BaselevelJunctions_file"] = "NULL";

  // Selecting basins
  int_default_map["threshold_contributing_pixels"] = 1000;
  int_default_map["minimum_basin_size_pixels"] = 1000;
  bool_default_map["extend_channel_to_node_before_receiver_junction"] = true;
  bool_default_map["test_drainage_boundaries"] = false;
  bool_default_map["only_take_largest_basin"] = false;
  bool_default_map["print_basin_raster"] = false;

  // printing of rasters and data before chi analysis
	bool_default_map["print_fill_raster"] = false;
	bool_default_map["write_hillshade"] = false;
  bool_default_map["print_slope"] = false;
  bool_default_map["print_aspect"] = false;
  bool_default_map["print_curvature"] = false;
  bool_default_map["print_planform_curvature"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_channels_to_csv"] = false;
  bool_default_map["print_junction_index_raster"] = false;
  bool_default_map["print_junctions_to_csv"] = false;

  bool_default_map["convert_csv_to_geojson"] = false;  // This converts all cv files to geojson (for easier loading in a GIS)

  // These are the HFR printing flags
  bool_default_map["run_HRF_analysis"] = false;
  bool_default_map["write_hilltops"] = false;
  bool_default_map["write_hilltop_curvature"] = false;
  bool_default_map["write_hillslope_length"] = false;
  bool_default_map["write_hillslope_gradient"] = false;
  bool_default_map["write_hillslope_relief"] = false;


 	// these params do not need changed during normal use of the HFR algorithm
	bool_default_map["print_hillslope_traces"] = false;
	int_default_map["hillslope_trace_thinning"] = 1;
	string_default_map["hillslope_traces_file"] = "";
	bool_default_map["hillslope_traces_basin_filter"] = false;

  // basic parameters for calculating chi
  float_default_map["A_0"] = 1;
  float_default_map["m_over_n"] = 0.5;
  int_default_map["threshold_pixels_for_chi"] = 0;
  bool_default_map["print_chi_data_maps"] = false;

    // These give unique IDs to each segment and then add this information to the
  // MChi and also semgent raster. We use this to map segments to other landscape
  // properties such as various hillslope metrics
  bool_default_map["print_segments"] = false;
  bool_default_map["print_segments_raster"] = false;

  // This burns a raster value to any csv output of chi data
  // Useful for appending geology data to chi profiles
  bool_default_map["burn_raster_to_csv"] = false;
  string_default_map["burn_raster_prefix"] = "NULL";
  string_default_map["burn_data_csv_column_header"] = "burned_data";
    
  // This burns a secondary raster value to any csv output of chi data
  // Useful when there are two datasets, e.g., precipitation data and geology data
  bool_default_map["secondary_burn_raster_to_csv"] = false;
  string_default_map["secondary_burn_raster_prefix"] = "NULL";
  string_default_map["secondary_burn_data_csv_column_header"] = "secondary_burned_data";


  // parameters for various chi calculations as well as slope-area
  int_default_map["n_iterations"] = 20;
  int_default_map["minimum_segment_length"] = 10;
  int_default_map["maximum_segment_length"] = 100000; //make super large so as not to be a factor unless user defined
  int_default_map["n_nodes_to_visit"] = 10;
  int_default_map["target_nodes"] = 80;
  int_default_map["skip"] = 2;
  float_default_map["sigma"] = 20;

  // these are routines that run segmentation
  bool_default_map["print_segmented_M_chi_map_to_csv"] = false;

  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  //cout << "PRINT THE PARAMETERS..." << endl;
  //LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  string BaselevelJunctions_file = LSDPP.get_BaselevelJunctions_file();
  string ChannelSegments_file = LSDPP.get_ChannelSegments_file();
  string Floodplain_file = LSDPP.get_Floodplain_file();

  cout << "Read filename is: " <<  DATA_DIR+DEM_ID << endl;
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;

    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);

  // check the threshold pixels for chi
  if (this_int_map["threshold_pixels_for_chi"] > this_int_map["threshold_contributing_pixels"])
  {
    cout << "WARNING: you have chosen a threshold pixels for chi which is greater" << endl;
    cout << "   the threshold contributing pixels. Defaulting so these are equal." << endl;
    this_int_map["threshold_pixels_for_chi"] = this_int_map["threshold_contributing_pixels"];
  }

  // Set the printing functions to true is you are running the HFR analysis
  if (this_bool_map["run_HRF_analysis"])
  {
    this_bool_map["write_hilltops"] = true;
    this_bool_map["write_hilltop_curvature"] = true;
    this_bool_map["write_hillslope_length"] = true;
    this_bool_map["write_hillslope_gradient"] = true;
    this_bool_map["write_hillslope_relief"] = true;
  }


  // initialise variables to be assigned from .driver file
  // These will all be assigned default values
  float A_0 = this_float_map["A_0"];
  float movern = this_float_map["m_over_n"];
  int n_iterations = this_int_map["n_iterations"];
  int minimum_segment_length = this_int_map["minimum_segment_length"];
  int n_nodes_to_visit = this_int_map["n_nodes_to_visit"];             // when constructing channel network, this
  float sigma = this_float_map["sigma"];
  int target_nodes = this_int_map["target_nodes"];
  int skip = this_int_map["skip"];
  int threshold_contributing_pixels = this_int_map["threshold_contributing_pixels"];
  int minimum_basin_size_pixels = this_int_map["minimum_basin_size_pixels"];

  //============================================================================
  //
  // RASTER BURNING
  //
  //============================================================================
  //============================================================================
  // check to see if the raster for burning exists
  LSDRaster BurnRaster;
  LSDRaster SecondaryBurnRaster;
  bool burn_raster_exists = false;
  bool secondary_burn_raster_exists = false;
  string burn_raster_header;
  string secondary_burn_raster_header;
  string burn_prefix;
  string secondary_burn_prefix;

  // Get al the filenames.
  burn_raster_header = DATA_DIR+this_string_map["burn_raster_prefix"]+".hdr";
  burn_prefix = DATA_DIR+this_string_map["burn_raster_prefix"];
  secondary_burn_raster_header = DATA_DIR+this_string_map["secondary_burn_raster_prefix"]+".hdr";
  secondary_burn_prefix = DATA_DIR+this_string_map["secondary_burn_raster_prefix"];

  
  if (this_bool_map["burn_raster_to_csv"])
  {
    cout << "I am going to burn a raster to all your csv files. The header name for this raster is: " << endl;
    cout <<  burn_raster_header << endl;
  }
  ifstream burn_head_in2;
  burn_head_in2.open(burn_raster_header.c_str());
  if( not burn_head_in2.fail() )
  {
    burn_raster_exists = true;
    string burn_fname = burn_prefix;
    cout << "The burn raster exists. It has a prefix of: " << endl;
    cout <<  burn_fname << endl;
    LSDRaster TempRaster(burn_fname,raster_ext);
    BurnRaster = TempRaster;
  }
  else
  {
    cout << "The burn raster doesn't exist! I am turning off the  burn flag" << endl;
    this_bool_map["burn_raster_to_csv"] = false;
  }

  // Adding a second burn raster    
  if (this_bool_map["burn_raster_to_csv"])
  {
    if (this_bool_map["secondary_burn_raster_to_csv"])
    {
    
      cout << "I am going to burn a secondary raster to all your csv files. The header name for this raster is: " << endl;
      cout <<  secondary_burn_raster_header << endl;
    }
    ifstream secondary_burn_head_in2;
    secondary_burn_head_in2.open(secondary_burn_raster_header.c_str());
    if( not secondary_burn_head_in2.fail() )
    {
      secondary_burn_raster_exists = true;
      string secondary_burn_fname = secondary_burn_prefix;
      cout << "The secondary burn raster exists. It has a prefix of: " << endl;
      cout <<  secondary_burn_fname << endl;
      LSDRaster TempRaster(secondary_burn_fname,raster_ext);
      SecondaryBurnRaster = TempRaster;
    }
    else
    {
      cout << "The secondary burn raster doesn't exist! I am turning off the  burn flag" << endl;
      this_bool_map["secondary_burn_raster_to_csv"] = false;
    }
  }
  else
  {
    cout << "Cannot burn a secondary raster if there isn't already a raster being burned to csv." << endl;
    cout << "I am turning off secondary raster burning!" << endl;
    this_bool_map["secondary_burn_raster_to_csv"] = false;
  }
  //============================================================================


  //============================================================================
  //
  // LOAD THE DEM
  //
  //============================================================================
  //============================================================================  
  LSDRaster topography_raster;
  if (this_bool_map["remove_seas"])
  {
    cout << "I am removing high and low values to get rid of things that should be nodata." << endl;
    LSDRaster start_raster((DATA_DIR+DEM_ID), raster_ext);
    // now get rid of the low and high values
    float lower_threshold = this_float_map["minimum_elevation"];
    float upper_threshold = this_float_map["maximum_elevation"];
    bool belowthresholdisnodata = true;
    LSDRaster Flooded = start_raster.mask_to_nodata_using_threshold(lower_threshold,belowthresholdisnodata);
    belowthresholdisnodata = false;
    topography_raster = Flooded.mask_to_nodata_using_threshold(upper_threshold,belowthresholdisnodata);
  }
  else
  {
    LSDRaster start_raster((DATA_DIR+DEM_ID), raster_ext);
    topography_raster = start_raster;
  }
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;


  float Resolution = topography_raster.get_DataResolution();
  map<string,string> GRS = topography_raster.get_GeoReferencingStrings();

  float thresh_area_for_chi = float(this_int_map["threshold_pixels_for_chi"])*Resolution*Resolution;

  if(this_bool_map["only_check_parameters"])
  {
    cout << "You set the only_check_parameters flag to true; I have only printed" << endl;
    cout << "the parameters to file and am now exiting." << endl;
    exit(0);
  }
  //============================================================================


  //============================================================================
  //
  // Start gathering necessary rasters
  //
  //============================================================================
  //============================================================================  
  LSDRaster filled_topography;
  // now get the flow info object
  if ( this_bool_map["raster_is_filled"] )
  {
    cout << "You have chosen to use a filled raster." << endl;
    filled_topography = topography_raster;
  }
  else
  {
    cout << "Let me fill that raster for you, the min slope is: "
         << this_float_map["min_slope_for_fill"] << endl;
    filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
  }

  if (this_bool_map["print_fill_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }

	//Surface fitting to get slope, aspect, curvature and planform curvature
	static const int Arr[] = {0,1,1,1,1,0,0,0};
	vector<int> RasterSelection (Arr, Arr + sizeof(Arr) / sizeof(Arr[0]));
  cout << endl << endl << "=========================================" << endl;
  cout << "I am starting the analysis. This is very memory intensive!! " << endl;
  cout << "If you get a segmentation fault the most likeley problem is that your DEM is too large." << endl;
  cout << "To fix you can either clip the DEM, reduce the resolution, or buy a more expensive computer." << endl;
  cout << "Note that if you are using docker you can change the settings to give your container more memory." << endl;
	cout << "\tCalculating surface metrics..." << endl;
	vector<LSDRaster> Surfaces = filled_topography.calculate_polyfit_surface_metrics(this_float_map["WindowRadius"], RasterSelection);
  LSDRaster Aspect = Surfaces[2];


  cout << "\t Flow routing..." << endl;
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

  // calculate the flow accumulation
  cout << "\t Calculating flow accumulation (in pixels)..." << endl;
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  cout << "\t Converting to flow area..." << endl;
  LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

  if (this_bool_map["print_DrainageArea_raster"])
  {
    string DA_raster_name = OUT_DIR+OUT_ID+"_DArea";
    DrainageArea.write_raster(DA_raster_name,raster_ext);
  }

  // calcualte the distance from outlet
  cout << "\t Calculating flow distance..." << endl;
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  cout << "\t Loading Sources..." << endl;
  cout << "\t Source file is... " << CHeads_file << endl;
  
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    cout << endl << "\t\t==================================" << endl;
    cout << "\t\tThe channel head file is null. " << endl;
    cout << "\t\tGetting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, threshold_contributing_pixels);

    cout << "\t\tThe number of sources is: " << sources.size() << endl;

  }
  else
  {
    cout << "\t\tLoading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t\tGot sources!" << endl;
  }

  // now get the junction network
  cout << "\t I am getting the junction network..." << endl;
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);


  // Print channels and junctions if you want them.
  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "\t\tprinting channels to csv..." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
    JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);

    // convert to geojson if that is what the user wants
    // It is read more easily by GIS software but has bigger file size
    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_CN.geojson";
      LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_CN.csv");
      thiscsv.print_data_to_geojson(gjson_name);
    }

  }

  // print junctions
  if( this_bool_map["print_junctions_to_csv"])
  {
    cout << "\t\tprinting junctions to csv..." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_JN.csv";
    JunctionNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_JN.geojson";
      LSDSpatialCSVReader thiscsv(channel_csv_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

  // Print sources
  if( this_bool_map["print_sources_to_csv"])
  {
    cout << "\t\tprinting sources to csv..." << endl;
    string sources_csv_name = OUT_DIR+OUT_ID+"_ATsources.csv";

    //write channel_heads to a csv file
    FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(sources, sources_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_ATsources.geojson";
      LSDSpatialCSVReader thiscsv(sources_csv_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }
  
  // Print stream order rasters if you want them
  if (this_bool_map["print_stream_order_raster"])
  {
    cout << "\t\tprinting stream order raster but this is really wasteful. Why don't you print to csv instead?" << endl;
    LSDIndexRaster SOArray = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
    string SO_raster_name = OUT_DIR+OUT_ID+"_SO";
    SOArray.write_raster(SO_raster_name,raster_ext);
  }
  if (this_bool_map["print_junction_index_raster"])
  {
    cout << "\t\tprinting junction raster but this is really wasteful. Why don't you print to csv instead?" << endl;
    LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
    string JI_raster_name = OUT_DIR+OUT_ID+"_JI";
    JIArray.write_raster(JI_raster_name,raster_ext);
  }

  //set up the stream network to get index raster for channel mapping
	LSDIndexRaster StreamNetwork;
	if (ChannelSegments_file == "NULL" || ChannelSegments_file == "Null" || ChannelSegments_file == "null")
  {
    cout << "\tGetting stream network from source nodes..." << endl;
	  StreamNetwork = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
  }
  else
  {
    cout << "\tGetting stream network from Channel segments file... " << DATA_DIR+ChannelSegments_file << endl;
    StreamNetwork = LSDIndexRaster((DATA_DIR+ChannelSegments_file), raster_ext);
  }

  //pad the stream network by a number of pixels
  if (this_int_map["StreamNetworkPadding"] > 0)
  {
    cout << "\tPadding the stream network by " << this_int_map["StreamNetworkPadding"] << " pixels..." << endl;
    StreamNetwork.PadRaster(this_int_map["StreamNetworkPadding"]);
  }

  //load floodplain and merge with the channel network if required, otherwise the
 	//floodplain mask will only contain the channel data
  if (Floodplain_file == "NULL" || Floodplain_file == "Null" || Floodplain_file == "null") {}
  else
  {
    cout << "\tCombining the channelnetwork and floodplain masks..." << endl;
		LSDIndexRaster Floodplains((Floodplain_file), raster_ext);
		StreamNetwork.MergeIndexRasters(Floodplains);
	}

  //Check to see if a list of junctions for extraction exists
  // need to get base-level nodes , otherwise these catchments will be missed!
  vector< int > BaseLevelJunctions;
  vector< int > BaseLevelJunctions_Initial;
  cout << "\tNow I am going to deal with the baselevel junctions. " << endl;
  cout << "\n\nThe BaselevelJunctions_file is: " << BaselevelJunctions_file << endl;
  if (BaselevelJunctions_file == "NULL" || BaselevelJunctions_file == "Null" || BaselevelJunctions_file == "null" || BaselevelJunctions_file.empty() == true)
  {
    cout << "\tSelecting basins..." << endl;
    // remove basins drainage from edge if that is what the user wants
    if (this_bool_map["find_complete_basins_in_window"])
    {
      cout << "\tFinding basins not influended by nodata and removing nested basins..." << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(FlowInfo, filled_topography, FlowAcc,
                                              this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
    }
    else
    {
      //Get baselevel junction nodes from the whole network
      BaseLevelJunctions_Initial = JunctionNetwork.get_BaseLevelJunctions();

      // now prune these by drainage area
      cout << "\tRemoving basins with fewer than " << this_int_map["minimum_basin_size_pixels"] << " pixels" << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions_Initial,
                                              FlowInfo, FlowAcc, this_int_map["minimum_basin_size_pixels"]);
      cout << "\tNow I have " << BaseLevelJunctions.size() << " baselelvel junctions." << endl;

      if (this_bool_map["find_largest_complete_basins"])
      {
        cout << "\tFinding largest basin not influenced by nodata..." << endl;
        BaseLevelJunctions = JunctionNetwork.Prune_To_Largest_Complete_Basins(BaseLevelJunctions,FlowInfo, filled_topography, FlowAcc);
      }
      else
      {
        if (this_bool_map["test_drainage_boundaries"])     // now check for edge effects
        {
          cout << endl << endl << "\tRemoving basins draining to the edge." << endl;
          BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions,FlowInfo, filled_topography);
          //BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge(BaseLevelJunctions,FlowInfo);
        }
      }
    }
  }
  else if (this_bool_map["extract_basins_by_stream_order"])
  {
    //Or just use a list of basins?
  	cout << "\tExtracting baselevel nodes for a strahler stream order of " << this_int_map["stream_order_to_extract_basins"] << "..." << endl;
  	BaseLevelJunctions = JunctionNetwork.ExtractBasinJunctionOrder(this_int_map["stream_order_to_extract_basins"], FlowInfo);

  }
  else
  {
    //specify junctions to work on from a list file
    string JunctionsFile = DATA_DIR+DEM_ID+"_junctions.list";

    cout << "\tReading junctions from a junction list... " << JunctionsFile << endl;

    vector<int> JunctionsList;
    ifstream infile(JunctionsFile.c_str());
    if (infile)
    {
      int n;
      while (infile >> n) BaseLevelJunctions_Initial.push_back(n);
    }
    else
    {
      cout << "Fatal Error: Junctions File " << JunctionsFile << " does not exist" << endl;
      exit(EXIT_FAILURE);
    }

    // Now make sure none of the basins drain to the edge
    cout << "\tPruning junctions that drain to the edge of the DEM..." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions_Initial, FlowInfo, filled_topography);
  }

  // Now we get channel segments for use with chi plotting. 
  vector<int> source_nodes;
  vector<int> outlet_nodes;
  vector<int> baselevel_node_of_each_basin;
  if (this_bool_map["extend_channel_to_node_before_receiver_junction"])
  {
    cout << endl << endl << "=====================================================" << endl;
    cout << "I am now getting the channels by basin." << endl;
    cout << "  These channels extend below the junction to the channel that stops" << endl;
    cout << "  just before the reciever junction. This option is used to remain" << endl;
    cout << "  consitent with basin ordering, since a 2nd order basin will begin" << endl;
    cout << "  at the channel one node upslope of the most upstream 3rd order junction." << endl;
    cout << "  If you simply want the channel starting from the selcted junction, " << endl;
    cout << "  set the option:" << endl;
    cout << "    extend_channel_to_node_before_receiver_junction" << endl;
    cout << "  to false." << endl;
    cout << "=====================================================" << endl << endl;

    JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

  }
  else
  {
    cout << endl << endl << "=====================================================" << endl;
    cout << "I am now getting the channels by basin." << endl;
    cout << "  These channels will start from the baselevel junctions selected. " << endl;
    cout << "  If you want them to extend to below the junction to the channel that stops" << endl;
    cout << "  just before the reciever junction, then set the option:" << endl;
    cout << "    extend_channel_to_node_before_receiver_junction" << endl;
    cout << "  to true." << endl;
    cout << "=====================================================" << endl << endl;

    JunctionNetwork.get_overlapping_channels(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);
  }

	//get the basins
	cout << "\tExtracting basins..." << endl;
	LSDIndexRaster BasinsRaster = JunctionNetwork.extract_basins_from_junction_vector(BaseLevelJunctions, FlowInfo);

  // Get the chi coordinate if needed
  LSDRaster chi_coordinate;
  if ( this_bool_map["print_basin_raster"] ||
        this_bool_map["print_chi_data_maps"] ||
        this_bool_map["print_segmented_M_chi_map_to_csv"])
  {
    cout << "I am calculating the chi coordinate." << endl;
    chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(this_float_map["m_over_n"],this_float_map["A_0"],this_int_map["threshold_contributing_pixels"]);
  
    // Print the basin raster if you want it:
    if(this_bool_map["print_basin_raster"])
    {
      cout << "I am going to print the basins for you. " << endl;
      LSDChiTools ChiTool_basins(FlowInfo);
      ChiTool_basins.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                              filled_topography, DistanceFromOutlet,
                              DrainageArea , chi_coordinate);
      string basin_raster_prefix = OUT_DIR+OUT_ID;
      ChiTool_basins.print_basins(FlowInfo, JunctionNetwork, BaseLevelJunctions, basin_raster_prefix);
    }


    //======================================================================
    // This is for visualisation of the channel data in the basin
    //======================================================================
    if(this_bool_map["print_chi_data_maps"])
    {
      cout << "I am going to print some simple chi data maps for visualisation." << endl;
      cout << "These data maps are also useful for visualising channel networks and making channel profiles." << endl;
      LSDChiTools ChiTool(FlowInfo);
      ChiTool.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                              filled_topography, DistanceFromOutlet,
                              DrainageArea, chi_coordinate);
  
      string chi_data_maps_string = OUT_DIR+OUT_ID+"_chi_data_map.csv";
      ChiTool.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);
  
      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_chi_data_map.geojson";
        LSDSpatialCSVReader thiscsv(chi_data_maps_string);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
    //======================================================================


    // This does all the segmenting. 
    // It uses a chi coordinate raster that has been inherited from earlier in this
    // program
    if (this_bool_map["print_segmented_M_chi_map_to_csv"])
    {
      cout << "\t\tI am going to segment your channel network." << endl;
      LSDChiTools ChiTool(FlowInfo);
      ChiTool.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                              filled_topography, DistanceFromOutlet,
                              DrainageArea, chi_coordinate);
      cout << "I am calculating the segmented channels" << endl;
      if (source_nodes.size() == 0)
      {
        cout << "I don't seem to have any source nodes!" << endl;
      }

      // check to see if we want segments. If that is the case, the
      // skip and iterations default to 0 and 1
      if (this_bool_map["print_segments"])
      {
        n_iterations = 1;
        skip = 0;
        ChiTool.chi_map_automator(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                              filled_topography, DistanceFromOutlet,
                              DrainageArea, chi_coordinate, target_nodes,
                              n_iterations, skip, minimum_segment_length, sigma);
        ChiTool.segment_counter(FlowInfo, this_int_map["maximum_segment_length"]);
        if (this_bool_map["print_segments_raster"])
        {
          LSDIndexRaster SegmentsRaster = ChiTool.segment_mapping(FlowInfo, this_int_map["maximum_segment_length"]);
          string Segments_raster_name = OUT_DIR+OUT_ID+"_Segments";
          SegmentsRaster.write_raster(Segments_raster_name,raster_ext);
        }
      }
      else
      {
        ChiTool.chi_map_automator(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                              filled_topography, DistanceFromOutlet,
                              DrainageArea, chi_coordinate, target_nodes,
                              n_iterations, skip, minimum_segment_length, sigma);
      }

      string csv_full_fname = OUT_DIR+OUT_ID+"_MChiSegmented.csv";
      cout << "Let me print all the data for you into a csv file called " << csv_full_fname << endl;
      ChiTool.print_data_maps_to_file_full(FlowInfo, csv_full_fname);
      
      //adding burned raster data to MChiSegmented output
      if (this_bool_map["burn_raster_to_csv"])  
      {  
          cout << "You asked me to burn a raster to the csv, I'm attemping to add this to your MChiSegmented output.";
          if(burn_raster_exists)
          {
              string header_for_burn_data = this_string_map["burn_data_csv_column_header"];        
              LSDSpatialCSVReader CSVFile(RI,csv_full_fname);
          
              cout << "I am burning the raster to the csv file." << endl;
              CSVFile.burn_raster_data_to_csv(BurnRaster,header_for_burn_data);
          
              string csv_full_fname_burned = OUT_DIR+OUT_ID+"_MChiSegmented_burned.csv";
              cout << "Now I'll print the data to a new file" << endl;
              CSVFile.print_data_to_csv(csv_full_fname_burned);
          }
          else
          {
              cout << "Burn raster does not exist!";
          }
      }
      if (this_bool_map["secondary_burn_raster_to_csv"])
      {
          cout << "You asked me to burn a secondary raster to the csv, attemping to add to MChiSegmented output.";
          if(secondary_burn_raster_exists)
          {
              string header_for_burn_data = this_string_map["secondary_burn_data_csv_column_header"];        
              string csv_full_fname_burned = OUT_DIR+OUT_ID+"_MChiSegmented_burned.csv";
              LSDSpatialCSVReader CSVFile(RI,csv_full_fname_burned);
          
              cout << "I am burning the secondary raster to the csv file." << endl;
              CSVFile.burn_raster_data_to_csv(SecondaryBurnRaster,header_for_burn_data);
          
              cout << "Now I'll print the data to a new file" << endl;
              CSVFile.print_data_to_csv(csv_full_fname_burned);
          }
          else
          {
              cout << "Secondary Burn raster does not exist!";
          }
      }
      cout << "That is your file printed!" << endl;

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        cout << "Now let me print your chi network to a geojson" << endl;
        string gjson_name = OUT_DIR+OUT_ID+"_MChiSegmented.geojson";
        LSDSpatialCSVReader thiscsv(csv_full_fname);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }


  }

	// Extract hilltops
	cout << "\tExtracting hilltop network..." << endl;
	LSDRaster Hilltops = JunctionNetwork.ExtractRidges(FlowInfo, this_int_map["min_stream_order_to_extract_basins"], this_int_map["max_stream_order_to_extract_basins"]);

  // Mask to only use hilltops inside basins being analysed
  if (this_bool_map["MaskHilltopstoBasins"])
  {
    cout << "\tMasking hilltops to basins being analysed..." << endl;
    Hilltops.MaskRaster(BasinsRaster);
  }

	if (this_bool_map["RemoveSteepHilltops"])
	{
	  cout << "\tFiltering out hilltops with slopes greater than " << this_float_map["Threshold_Hilltop_Gradient"] << "..." << endl;
	  string Condition = "<";
    LSDIndexRaster SteepHilltopsMask = Surfaces[1].Create_Mask(Condition,this_float_map["Threshold_Hilltop_Gradient"]);
    Hilltops.MaskRaster(SteepHilltopsMask);
  }

  if (this_bool_map["RemovePositiveHilltops"])
	{
	  cout << "\tFiltering out hilltops with positive curvature..." << endl;
	  float Zero = 0;
	  string Condition = "<";
    LSDIndexRaster NegativeCurvatureMask = Surfaces[3].Create_Mask(Condition,Zero);
    Hilltops.MaskRaster(NegativeCurvatureMask);
	}

  //get hilltop curvature using filters to remove positive curvatures and steep slopes
  cout << "\tGetting hilltop curvature..." << endl;
	LSDRaster CHT = filled_topography.get_hilltop_curvature(Surfaces[3], Hilltops);

  if (this_bool_map["run_HFR_analysis"])
  {
    // Run hilltop flow routing
    vector<int> Target_Basin_Vector;
    cout << "\tRunning hillslope flow routing. You might want to find a good book..." << endl;
    string HillslopeTracesFolder = DATA_DIR;
    vector< Array2D<float> > HFR_Arrays = FlowInfo.HilltopFlowRouting(filled_topography, Hilltops, Surfaces[1], Surfaces[2], Surfaces[3], Surfaces[4], StreamNetwork, BasinsRaster, (DATA_DIR+DEM_ID), this_bool_map["print_hillslope_traces"], this_int_map["hillslope_trace_thinning"], HillslopeTracesFolder, this_bool_map["hillslope_traces_basin_filter"], Target_Basin_Vector);
    cout << "\tDone" << endl;  

    //hilltops
    if (this_bool_map["write_hilltops"])
    {
      string hilltops_raster_name = OUT_DIR+OUT_ID+"_Hilltops";
      cout << "\tWriting hilltops to " << hilltops_raster_name << raster_ext << "..." << endl;
      Hilltops.write_raster(hilltops_raster_name, raster_ext);
      cout << "done.";
    }
    //hilltop curvature
    if (this_bool_map["write_hilltop_curvature"])
    {
      string hilltop_curvature_raster_name = OUT_DIR+OUT_ID+"_HFR_CHT";
      cout << "\tWriting hilltop curvature to " << hilltop_curvature_raster_name << raster_ext << "..." << endl;
      CHT.write_raster(hilltop_curvature_raster_name, raster_ext);
      cout << "done.";
    }
    //hillslope length raster
    if (this_bool_map["write_hillslope_length"])
    {
      string hillslope_length_raster_name = OUT_DIR+OUT_ID+"_HFR_LH";
      cout << "\tWriting hillslope_lengths to " << hillslope_length_raster_name << raster_ext << "..." << endl;
      LSDRaster HFR_LH = CHT.LSDRasterTemplate(HFR_Arrays[1]);
      HFR_LH.write_raster(hillslope_length_raster_name,raster_ext);
      cout << "done.";
    }
    //hillslope gradient raster
    if (this_bool_map["write_hillslope_gradient"])
    {
      string hillslope_gradient_raster_name = OUT_DIR+OUT_ID+"_HFR_SLP";
      cout << "\tWriting hillslope gradients to " << hillslope_gradient_raster_name << raster_ext << "..." << endl;
      LSDRaster HFR_Slope = CHT.LSDRasterTemplate(HFR_Arrays[2]);
      HFR_Slope.write_raster(hillslope_gradient_raster_name,raster_ext);
      cout << "done.";
    }
    //hillslope relief raster
    if (this_bool_map["write_hillslope_relief"])
    {
      string hillslope_relief_raster_name = OUT_DIR+OUT_ID+"_HFR_Relief";
      cout << "\tWriting hillslope relief to " << hillslope_relief_raster_name << raster_ext << "..." << endl;
      LSDRaster relief = CHT.LSDRasterTemplate(HFR_Arrays[3]);
      relief.write_raster(hillslope_relief_raster_name,raster_ext);
      cout << "done.";
    }
  }



  // Write rasters to file
  //fill
  if (this_bool_map["print_fill_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    cout << "\tWriting filled topography to "  << filled_raster_name <<raster_ext << endl;
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }

  //hillshade
  if (this_bool_map["write_hillshade"])
  {
    string hillshade_raster_name = OUT_DIR+OUT_ID+"_Hillshade";
    cout << "Writing hillshade to " << hillshade_raster_name << raster_ext << endl;
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster HillshadeRaster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);
    HillshadeRaster.write_raster(hillshade_raster_name,raster_ext);
  }

  //slope
  if (this_bool_map["print_slope"])
  {
    string slope_raster_name = OUT_DIR+OUT_ID+"_Slope";
    cout << "\tWriting slope to " << slope_raster_name << raster_ext << "..." << endl;
    Surfaces[1].write_raster((OUT_DIR+OUT_ID+"_Slope"),raster_ext);
    cout << "done.";
  }
  //aspect
  if (this_bool_map["print_aspect"])
  {
    string aspect_raster_name = OUT_DIR+OUT_ID+"_Aspect";
    cout << "\tWriting aspect to " << aspect_raster_name << raster_ext << "..." << endl;
    Surfaces[2].write_raster(aspect_raster_name,raster_ext);
    cout << "done.";
  }
  //curvature
  if (this_bool_map["print_curvature"])
  {
    string curvature_raster_name = OUT_DIR+OUT_ID+"_Curvature";
    cout << "\tWriting curvature to " << curvature_raster_name << raster_ext << "..." << endl;
    Surfaces[3].write_raster(curvature_raster_name,raster_ext);
    cout << "done.";
  }
  //planform curvature
  if (this_bool_map["print_planform_curvature"])
  {
    string plancurv_raster_name = OUT_DIR+OUT_ID+"_PlanCurv";
    cout << "\tWriting planform curvature to " << plancurv_raster_name << raster_ext << "..." << endl;
    Surfaces[4].write_raster(plancurv_raster_name, raster_ext);
    cout << "done.";
  }
  //stream network
  if (this_bool_map["print_stream_order_raster"])
  {
    string stnet_raster_name = OUT_DIR+OUT_ID+"_MaskedStreamNetwork";
    cout << "\tWriting stream newtork mask to " << stnet_raster_name << raster_ext << "..." << endl;
    StreamNetwork.write_raster(stnet_raster_name, raster_ext);
    cout << "done.";
  }
  //basins
  if (this_bool_map["print_basins"])
  {
    string basins_raster_name = OUT_DIR+OUT_ID+"_Basins";
    cout << "\tWriting basins mask to " << basins_raster_name << raster_ext << "..." << endl;
    BasinsRaster.write_raster(basins_raster_name, raster_ext);
    cout << "done.";
  }


}
