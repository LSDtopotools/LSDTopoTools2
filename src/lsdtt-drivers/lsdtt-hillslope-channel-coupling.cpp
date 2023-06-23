 //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// lsdtt-hillslope-channel-coupling.cpp
//
// This program combines the chi analysis of Mudd et al. 2014 and the hillslope analysis
// of Hurst et al 2012 and Grieve et al. 2016 to allow combined topographic analysis comparing hillslope and
// channel metrics for erosion rate in order to explore landscape morphology and transience.
//
// This program was first deployed in Hurst et al., 2019 EPSL
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate
// different kinds of chi analysis
//
// Call with -h to generate a help file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2022 Martin D. Hurst and Simon M. Mudd 2022
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
#include <chrono>
#include "../LSDStatsTools.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDRasterSpectral.hpp"
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

  string version_number = "0.9";
  string citation = "http://doi.org/10.5281/zenodo.4577879";

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the coupled channel and hillslope tool!  ||" << endl;
  cout << "|| This program has a number of options to make chi    ||" << endl;
  cout << "|| plots and to map out slopes in chi space.           ||" << endl;
  cout << "|| This program was developed by Martin D. Hurst       ||" << endl;
  cout << "||  at the University of Glasgow                       ||" << endl;
  cout << "|| and Simon M. Mudd and Stuart W.D. Grieve            ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| Citation for this code is:                          ||" << endl;
  cout << "|| " << citation << endl;
  cout << "|| If you use these routines please cite:              ||" << endl;
  cout << "|| Grieve et al. 2016. How long is a hillslope?        ||" << endl;
  cout << "|| ESPL 41, 1039–1054.                                 ||" << endl;
  cout << "|| https://doi.org/10.1002/esp.3884                    ||" << endl;
  cout << "|| and                                                 ||" << endl;
  cout << "|| Hurst et al., 2019, Detection of channel-hillslope  ||" << endl;
  cout << "|| coupling along a tectonic gradient, EPSL, 522, 30-39||" << endl;
  cout << "|| , https://doi.org/10.1016/j.epsl.2019.06.018        ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| Documentation can be found at:                      ||" << endl;
  cout << "|| https://lsdtopotools.github.io/LSDTT_documentation/ ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| This is LSDTopoTools2 version                       ||" << endl;
  cout << "|| " << version_number << endl;
  cout << "|| If the version number has a d at the end it is a    ||" << endl;
  cout << "||  development version.                               ||" << endl;
  cout << "=========================================================" << endl;
  
  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);
  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

  // Check if we are doing the version or the citation
  if(f_name == "lsdtt_citation.txt")
  {

    cout << endl << endl << endl << "==============================================" << endl;
    cout << "To cite this code, please use this citation: " << endl;
    cout << citation << endl;
    cout << "Copy this url to find the full citation." << endl;
    cout << "also see above for more detailed citation information." << endl;
    cout << "=========================================================" << endl;

    ofstream ofs;
    ofs.open("./lsdtt-hillslope-channel-coupling-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    cout << endl << endl << endl << "==============================================" << endl;    
    cout << "This is lsdtt-hillslope-channel-coupling version number " << version_number << endl;
    cout << "If the version contains a 'd' then you are using a development version." << endl;
    cout << "=========================================================" << endl;
    ofstream ofs;
    ofs.open("./lsdtt-hillslope-channel-coupling-version.txt");
    ofs << version_number << endl;
    ofs.close();

    exit(0);
  }


  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // this will contain the help file
  map< string, vector<string> > help_map;

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  //==================================================================================
  //
  // .#####....####...#####....####...##...##..######..######..######..#####....####..
  // .##..##..##..##..##..##..##..##..###.###..##........##....##......##..##..##.....
  // .#####...######..#####...######..##.#.##..####......##....####....#####....####..
  // .##......##..##..##..##..##..##..##...##..##........##....##......##..##......##.
  // .##......##..##..##..##..##..##..##...##..######....##....######..##..##...####..
  //
  //=================================================================================
  // Basic DEM preprocessing
  float_default_map["minimum_elevation"] = 0.0;
  help_map["minimum_elevation"] = { "float","0.0","All elevation values below this become nodata if remove_seas is true.","Usually 0."};

  float_default_map["maximum_elevation"] = 30000;
  help_map["maximum_elevation"] = {  "float","0.0","All elevation values above this become nodata if remove_seas is true.","Pick a big number."};

  float_default_map["min_slope_for_fill"] = 0.0001;
  help_map["min_slope_for_fill"] = {  "float","0.0001","Minimum slope between pixels for the filling algorithm.","Best not to change the default."};

  bool_default_map["raster_is_filled"] = false; // assume base raster is already filled
  help_map["raster_is_filled"] = {  "bool","false","This reads a pre-existing fill raster to save time.","You need to have printed the fill raster if you set this to true."};

  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  help_map["remove_seas"] = {  "bool","true","Slightly misleading name; it replaces both high and low DEM values with nodata.","This gets rid of low lying areas but also is handy when the nodata is not translated from the raw DEM and it is full of funny large numbers."};

  bool_default_map["only_check_parameters"] = false;
  help_map["only_check_parameters"] = {  "bool","false","This just checks parameters without running an analysis.","For bug checking."};

  float_default_map["surface_fitting_radius"] = 12;
  help_map["surface_fitting_radius"] = {  "float","12","Our surface fitting routines fit a polynomial over the points with in a radius defined by surface_fitting_radius and then differentiate this surface to get the surface metrics like gradient and curvature","If not bigger than the pixel_size*sqrt(2) then will increase to that number."};

  bool_default_map["use_spiral_trimmer_for_edge_influence"] = false;
  help_map["use_spiral_trimmer_for_edge_influence"] = {  "bool","false","Makes sure raster is rectangular before running drainage extraction.","Use this if you want to avoid edge effects of your flow routing that can cause seg faults."};


  // Channel extraction
  bool_default_map["print_wiener_channels"] = false;
  help_map["print_wiener_channels"] = {  "bool","false","Prints the channel network determined by the Wiener method which is a mashup of Passalacqua and Pelletier methods first reported by grieve et al 2016.","Output is a csv file."};

  bool_default_map["print_perona_malik_channels"] = false;
  help_map["print_perona_malik_channels"] = {  "bool","false","Prints the channel network determined by the Perona-Malik method (Passalacqua et al., 2010)","Output is a csv file."};

  float_default_map["pruning_drainage_area"] = 1000;
  help_map["pruning_drainage_area"] = {  "float","1000","In both driech and Wiener channels with less than this drainage area are removed during the pruning process.","In m^2."};

  int_default_map["threshold_contributing_pixels"] = 1000;
  help_map["threshold_contributing_pixels"] = {  "int","1000","The number of contributing pixels needed to start a channel using the threshold method.","This is in pixels not drainage area. More options are in the lsdtt-channel-extraction tool."};

  int_default_map["connected_components_threshold"] = 100;
  help_map["connected_components_threshold"] = {  "int","100","Number of connected pixels to classify as a channel.","Used in Pelletier and Wiener methods."};

   //Defining hilltops
	int_default_map["StreamNetworkPadding"] = 0;
  help_map["StreamNetworkPadding"] = {  "int","0","Distance in pixels from channel network that can't be classed as a hilltop","Use if you want to remove hilltops near channels."};

  bool_default_map["extract_ridges"] = false;
  help_map["extract_ridges"] = {  "bool","false","Does what it says on the tin.","Output is a csv file. If you set run_HFR_analysis to true this will also be set to true by default"};

	float_default_map["Threshold_Hilltop_Gradient"] = 0.4;
  help_map["Threshold_Hilltop_Gradient"] = {  "float","0.4","Hilltops with greater than this gradient are removed from the analysis.","Used to remove plunging ridgelines."};

  float_default_map["Threshold_Hilltop_Curvature"] = 0.0;
  help_map["Threshold_Hilltop_Curvature"] = {  "float","0.4","Hilltops with greater than this curvature are removed from the analysis.","Used to remove concave ridgetops."};

  int_default_map["Threshold_Hilltop_Contributing_Pixels"] = 1;
  help_map["Threshold_Hilltop_Contributing_Pixels"] = {  "int","1","Hilltops with greater than this number of contributing pixels are removed from the analysis.","Used to remove ridgetops that are getting sediment from adjacent hilltop pixels."};

  int_default_map["Threshold_Distance_From_Channel"] = 0;
  help_map["Threshold_Distance_From_Channel"] = { "int","0","Hilltops closer to the channel than this distance (in metres) will be removed from the ridge extraction.","Used to filter network to remove pixels close to channel. You can set the threshold stream order that will be used using the Threshold_Stream_Order parameter."};
	 
  int_default_map["Threshold_Stream_Order"] = 1;
  help_map["Threshold_Stream_Order"] = { "int","1","Used to determine which channels will be used to filter the ridge network to remove hilltops close to the channel network.","Default is 1 which means all channels will be considered."};


  // This is all for legacy HFR
  bool_default_map["complete_hilltop_masking"] = true;
  help_map["complete_hilltop_masking"] = {  "bool","false","Legacy hilltop flow routing. Runs all the masking routines.","Legacy code only. Use run_HFR_analysis now."};

  bool_default_map["use_legacy_HFR"] = false;
  help_map["use_legacy_HFR"] = {  "bool","false","Use the legacy hilltop flow routing.","Legacy code only. Use run_HFR_analysis now with this set to false."};

  bool_default_map["RemovePositiveHilltops"] = true;
  help_map["RemovePositiveHilltops"] = {  "bool","false","Removes positive (concave) hilltops in legacy hilltop flow routing.","Legacy code only. This defaults to true for new version."};
 
	bool_default_map["RemoveSteepHilltops"] = true;
  help_map["RemoveSteepHilltops"] = {  "bool","false","Removes steep (>0.4 gradient) hilltops in legacy hilltop flow routing.","Legacy code only. This defaults to true for new version."};
 
  bool_default_map["MaskHilltopstoBasins"] = true;
  help_map["MaskHilltopstoBasins"] = {  "bool","false","Removes hilltops not around study basins in legacy hilltop flow routing.","Legacy code only. This defaults to true for new version."};
 
  int_default_map["min_stream_order_to_extract_basins"] = 0;
  help_map["min_stream_order_to_extract_basins"] = {  "int","0","Minimum basin order in legacy hilltop flow routing.","Legacy code only."};
 
	int_default_map["max_stream_order_to_extract_basins"] = 100;
  help_map["max_stream_order_to_extract_basins"] = {  "int","100","Maximum basin order in legacy hilltop flow routing.","Legacy code only."};

  // Input filenames
  string_default_map["ChannelSegments_file"] = "NULL";
  help_map["ChannelSegments_file"] = {  "string","NULL","If you already calculated a channel segments file you can use this to load the file.","This will be a csv file that you can get by using the print_segments option."};

	string_default_map["Floodplain_file"] = "NULL";
  help_map["Floodplain_file"] = {  "string","NULL","The prefix (no extension) of the floodplain mask. Is a binary with true (1) for the floodplain. Hilltop traces will stop here","Floodplain mask can be generated using lsdtt-valley-metrics."};

  string_default_map["CHeads_file"] = "NULL";
  help_map["CHeads_file"] = {  "string","NULL","The name of a channel heads file.","You can output this csv file with the channel extraction algorithms. It contains latitude and longitude values of the channel heads."};
 
  bool_default_map["get_basins_from_outlets"] = false;
  help_map["get_basins_from_outlets"] = {  "bool","false","Switches on the outlet based basin finding.","See BaselevelJunctions_file for format of outlets csv."};

  int_default_map["search_radius_nodes"] = 8;
  help_map["search_radius_nodes"] = {  "int","8","A parameter for snapping to the nearest channel. It will search for the largest channel (by stream order) within the pixel window.","You will want smaller pixel numbers if you have a dense channel network."};
 
  int_default_map["threshold_stream_order_for_snapping"] = 1;
  help_map["threshold_stream_order_for_snapping"] = {  "int","2","If you are snapping to a channel this routine will ignore channel with lower stream order than this number.","Set this to a higher number to avoid snapping to small channels."};
  
  string_default_map["BaselevelJunctions_file"] = "NULL";
  help_map["BaselevelJunctions_file"] = {  "string","NULL","The name of a csv file with basin outlets for selecting basins using junction numbers.","An old method. You should use get_basins_from_outlets: true and basin_outlet_csv instead."};
  
  string_default_map["basin_outlet_csv"] = "NULL";
  help_map["basin_outlet_csv"] = {  "string","NULL","A csv file with the lat long of basin outlets.","csv should have latitude and longitude columns and rows with basin outlets."};
  

  // Selecting basins
  int_default_map["minimum_basin_size_pixels"] = 1000;
  help_map["minimum_basin_size_pixels"] = {  "int","1000","For basin finding algorithm this value is the minimum size of a selected basin.","Will reject basins along edge."};
  
  int_default_map["maximum_basin_size_pixels"] = 1000000000;
  help_map["maximum_basin_size_pixels"] = {  "int","1000000000","For basin finding algorithm this value is the maximum size of a selected basin.","Will reject basins along edge."};
  
  bool_default_map["extend_channel_to_node_before_receiver_junction"] = true;
  help_map["extend_channel_to_node_before_receiver_junction"] = {  "bool","true","For various basin extractions the basin snaps to the nearest junction. If this is true then the outlet of the basin is one pixel upstream of the receiver junction of the snapped channel.","If false it will pick the donor junction of the channel rather than one pixel above the receiver."};

  bool_default_map["test_drainage_boundaries"] = true;
  help_map["test_drainage_boundaries"] = {  "bool","true","Looks for basins influenced by edge and removes them if they are.","chi coordinate must be calculated using complete basins so this tests for that."};
  
  bool_default_map["only_take_largest_basin"] = false;
  help_map["only_take_largest_basin"] = {  "bool","false","This only retains the largest complete basin in the raster.","Will reject basins along edge."};
  
  bool_default_map["print_basin_raster"] = false;
  help_map["print_basin_raster"] = {  "bool","false","This prints a raster where the values are the basin number.","You can combine this with python tools to get basin shapefiles."};
    

  // printing of rasters and data before chi analysis
  bool_default_map["print_fill_raster"] = false;
  help_map["print_fill_raster"] = {  "bool","false","Prints the fill raster.","Filename includes _FILL"};

  bool_default_map["write_hillshade"] = false;
  help_map["write_hillshade"] = {  "bool","false","Write the hillshade raster.","You need this for a lot of our plotting routines. Filename includes _HS"};

  bool_default_map["print_slope"] = false;
  help_map["print_slope"] = {  "bool","false","Prints slope raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_aspect"]= false;
  help_map["print_aspect"] = {  "bool","false","Prints aspect raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_curvature"]= false;
  help_map["print_curvature"] = {  "bool","false","Prints curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_planform_curvature"]= false;
  help_map["print_planform_curvature"] = {  "bool","false","Prints planform curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_stream_order_raster"] = false;
  help_map["print_stream_order_raster"] = {  "bool","false","Prints a raster with _SO in filename with stream orders of channel in the appropriate pixel.","Generates a big file so we suggest printing the network to csv."};

  bool_default_map["print_channels_to_csv"] = false;
  help_map["print_channels_to_csv"] = {  "bool","false","Prints the channel network to a csv file.","This version produces smaller files than the raster version."};

  bool_default_map["use_extended_channel_data"] = false;
  help_map["use_extended_channel_data"] = {  "bool","false","If this is true you get more data columns in your channel network csv.","I will tell you what these columns are one day."};

  bool_default_map["print_junction_index_raster"] = false;
  help_map["print_junction_index_raster"] = {  "bool","false","Prints a raster with junctions and their number.","Makes big files. It is better to use the csv version."};

  bool_default_map["print_junctions_to_csv"] = false;
  help_map["print_junctions_to_csv"] = {  "bool","false","Prints a csv with the locations and numbers of the junctions.","This is better to use than the raster version."};
 
  bool_default_map["print_sources_to_csv"] = false;
  help_map["print_sources_to_csv"] = {  "bool","false","Prints the sources to a csv file.","Each source on its own row with latitude and longitude columns."};

  bool_default_map["convert_csv_to_geojson"] = false;
  help_map["convert_csv_to_geojson"] = {  "bool","false","Converts csv files to geojson files","Makes csv output easier to read with a GIS. Warning: these files are much bigger than csv files."};

  bool_default_map["print_channels_to_line_vector"] = false;
  help_map["print_channels_to_line_vector"] = { "bool","true","If this is true print a geojson with the channel network as a polyline.","Used for more efficient visualisation of the channel network."};


  // These are the HFR printing flags
  bool_default_map["run_HFR_analysis"] = false;
  help_map["run_HFR_analysis"] = {  "bool","false","Run the main hilltop flow routing code","Makes csv output of hilltop curvature and hillslope metrics."};

  bool_default_map["map_hilltops_to_channels"] = false;
  help_map["map_hilltops_to_channels"] = {  "bool","false","This connects channel pixels to the hilltop pixels so you can compare channel steepness to hilltops","Makes csv output of hillslope metrics with channel metrics connected to each hilltop pixel."};

  // The below are for the legacy code
  bool_default_map["write_hilltops"] = false;
  help_map["write_hilltops"] = {  "bool","false","Writes hilltops to a raster. Legacy code.","Legacy code. New code writes to csv."};

  bool_default_map["write_hilltop_curvature"] = false;
  help_map["write_hilltop_curvature"] = {  "bool","false","Writes hilltop curvature to a raster. Legacy code.","Legacy code. New code writes to csv."};

  bool_default_map["write_hillslope_length"] = false;
  help_map["write_hillslope_length"] = {  "bool","false","Writes hilltop length to a raster. Legacy code.","Legacy code. New code writes to csv."};

  bool_default_map["write_hillslope_gradient"] = false;
  help_map["write_hillslope_gradient"] = {  "bool","false","Writes hilltop gradient to a raster. Legacy code.","Legacy code. New code writes to csv."};

  bool_default_map["write_hillslope_relief"] = false;
  help_map["write_hillslope_relief"] = {  "bool","false","Writes hilltop relief to a raster. Legacy code.","Legacy code. New code writes to csv."};

  // set these if you want to read in a premade curvature raster
  bool_default_map["read_curvature_raster"] = false;
  help_map["read_curvature_raster"] = {  "bool","false","Reads a precalculated curvature raster.","Use this if you already calculated curvature."};

  string_default_map["curvature_fname"] = "NULL";
  help_map["curvature_fname"] = {  "string","NULL","The raster prefix (without extension) of a precalculated curvature raster.","You need to set read_curvature_raster to true to read the raster."};


 	// these params do not need changed during normal use of the HFR algorithm
	bool_default_map["print_hillslope_traces"] = false;
  help_map["print_hillslope_traces"] = {  "bool","false","Prints the (very large) hillslope traces files for plotting traces.","Warning: results in a huge file."};

	int_default_map["hillslope_trace_thinning"] = 1;
  help_map["hillslope_trace_thinning"] = {  "int","1","The number of traces to print. If 1 it prints every trace. If 2 prints every other trace. And so on.","Use this if you already calculated curvature."};
 
	bool_default_map["hillslope_traces_basin_filter"] = false;
  help_map["hillslope_traces_basin_filter"] = {  "bool","false","If true only traces from selected basins will be printed.","Used to thin hillslope traces."};


  // basic parameters for calculating chi
  float_default_map["A_0"] = 1.0;
  help_map["A_0"] = {  "float","1.0","The A_0 parameter for chi computation. See https://doi.org/10.1002/esp.3302","Usually set to 1 so that the slope in chi-elevation space is the same as k_sn"};
   
  float_default_map["m_over_n"] = 0.5;
  help_map["m_over_n"] = {  "float","0.5","The concavity index for chi calculations. Usually denoted as the Greek symbol theta.","Default is 0.5 but possibly 0.45 is better as Kwang and Parker suggest 0.5 leads to unrealistic behaviour in landscape evolution models."};

  bool_default_map["print_chi_data_maps"] = false;
  help_map["print_chi_data_maps"] = {  "bool","false","If true prints the chi network to csv.","csv file has chidatamaps in the filename. Has the locations of the channel pixels with their chi coordinates and other information."};
  
  int_default_map["threshold_pixels_for_chi"] = 0;
  help_map["threshold_pixels_for_chi"] = {  "int","0","Minimum number of contributing pixels for calculating chi.","You can reduce the size of the chi csv files by setting a high number."};


  // These give unique IDs to each segment and then add this information to the
  // MChi and also semgent raster. We use this to map segments to other landscape
  // properties such as various hillslope metrics
  bool_default_map["print_segments"] = false;
  help_map["print_segments"] = {  "bool","false","This will print a channel csv that includes segment numbers.","csv file will have lat-long chi and other information but also segment numbers."};

  bool_default_map["print_segments_raster"] = false;
  help_map["print_segments_raster"] = {  "bool","false","Prints a raster of segment numbers.","For visualisation. The print_segments is more data efficient."};


  // This burns a raster value to any csv output of chi data
  // Useful for appending geology data to chi profiles
  bool_default_map["burn_raster_to_csv"] = false;
  help_map["burn_raster_to_csv"] = {  "bool","false","Takes a raster with burn_raster_prefix and then samples that raster with the points in the csv file. The new column will be burn_data_csv_column_header.","Useful for adding raster data to csv file. Often used to add lithological information to csv data (you must rasterize the lithology data first."};

  string_default_map["burn_raster_prefix"] = "NULL";
  help_map["burn_raster_prefix"] = {  "string","NULL","The prefix of the raster to burn to a csv.","No extension required."};

  string_default_map["burn_data_csv_column_header"] = "burned_data";
  help_map["burn_data_csv_column_header"] = {  "string","burned_data","Column header in csv of data burned from raster.","For example lithocode."};

  // This burns a secondary raster value to any csv output of chi data
  // Useful when there are two datasets, e.g., precipitation data and geology data
  bool_default_map["secondary_burn_raster_to_csv"] = false;
  help_map["secondary_burn_raster_to_csv"] = {  "bool","false","Takes a second raster with burn_raster_prefix and then samples that raster with the points in the csv file. The new column will be burn_data_csv_column_header.","Useful for adding raster data to csv file. Often used to add lithological information to csv data (you must rasterize the lithology data first."};

  string_default_map["secondary_burn_raster_prefix"] = "NULL";
  help_map["secondary_burn_raster_prefix"] = {  "string","NULL","The prefix of the second raster to burn to a csv.","No extension required."};

  string_default_map["secondary_burn_data_csv_column_header"] = "secondary_burned_data";
  help_map["secondary_burn_data_csv_column_header"] = {  "string","secondary_burned_data","Column header in csv of data burned from second raster.","For example lithocode."};


  // parameters for various chi calculations as well as slope-area
  int_default_map["n_iterations"] = 20;
  help_map["n_iterations"] = {  "int","20","Number of iterations of random sampling of chi-elevation space to make segments. Used to constrain segmentation uncertainty.","See Mudd et al 2014 JGR-ES for details.."};

  int_default_map["minimum_segment_length"] = 10;
  help_map["minimum_segment_length"] = {  "int","10","The minimum number of pixels in a segment. If too short computation is very expensive.","See Mudd et al 2014 JGR-ES for details. Sensitivity testing suggest values between 8 and 14 are appropriate."};

  int_default_map["maximum_segment_length"] = 100000; //make super large so as not to be a factor unless user defined
  help_map["maximum_segment_length"] = {  "int","100000","The maximum number of pixels in a segment. Usually not used so set to a very high number.","See Mudd et al 2014 JGR-ES for details."};

  int_default_map["n_nodes_to_visit"] = 10;
  help_map["n_nodes_to_visit"] = {  "int","10","The number of nodes downslope of a junction to visit for segmentation across junctions.","See Mudd et al 2014 JGR-ES for details."};

  int_default_map["target_nodes"] = 80;
  help_map["target_nodes"] = {  "int","80","Segmentation breaks channels into chunks of this length. It tests all combinations within that window so computation time increases a lot when this is bigger than 120.","See Mudd et al 2014 JGR-ES for details. 80 is okay. 120 will take forever. 60 too short"};

  int_default_map["skip"] = 2;
  help_map["skip"] = {  "int","2","A parameter used in segmentation. This is the mean number of pixels skipped for each pixel selected. The actual number skipped in each iteration is uniformly distributed between 0 and 2*skip. This approach is taken to constrain uncertainty in segments.","See Mudd et al 2014 JGR-ES for details. A number between 2 and 4 is normal. Larger numbers can be used for lidar. "};

  float_default_map["sigma"] = 20;
  help_map["sigma"] = {  "float","20","A parameter used in segmentation. This should be the geomorphic noise in metres so the topographic error plus the noise in the channels from boulders local erodibility variability etc.","See Mudd et al 2014 JGR-ES for details. Will vary between 5-20 unless you are using terrible ASTER data in which case use a bigger number. Never more than 50."};

  // these are routines that run segmentation
  bool_default_map["print_segmented_M_chi_map_to_csv"] = false;
  help_map["print_segmented_M_chi_map_to_csv"] = {  "bool","false","This runs the Mudd et al 2014 JGR segmentation routine so you get all the chi slope segments.","See Mudd et al 2014 JGR-ES for details. If you want maps of k_sn this is the way to do it."};

 
  //=========================================================================
  //
  //.#####....####...#####....####...##...##..######..######..######..#####..
  //.##..##..##..##..##..##..##..##..###.###..##........##....##......##..##.
  //.#####...######..#####...######..##.#.##..####......##....####....#####..
  //.##......##..##..##..##..##..##..##...##..##........##....##......##..##.
  //.##......##..##..##..##..##..##..##...##..######....##....######..##..##.
  //
  //..####...##..##..######...####...##..##...####..                         
  //.##..##..##..##..##......##..##..##.##...##.....                         
  //.##......######..####....##......####.....####..                         
  //.##..##..##..##..##......##..##..##.##.......##.                         
  //..####...##..##..######...####...##..##...####..    
  //                     
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
  string Floodplain_file = this_string_map["Floodplain_file"];

  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./lsdtt-hillslope-channel-coupling-README.csv" << endl;
    string help_prefix = "lsdtt-hillslope-channel-coupling-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }


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

  if (this_bool_map["map_hilltops_to_channels"])
  {
    cout << "I will map hilltops to channels, so I am setting up the HFR and segmentation." << endl;
    this_bool_map["run_HFR_analysis"] = true;
    this_bool_map["print_segmented_M_chi_map_to_csv"] = true;
  }

  // Set the printing functions to true is you are running the HFR analysis
  if (this_bool_map["run_HFR_analysis"])
  {
    cout << "I am running HFR analysis so I am setting this up with ridge extraction." << endl;
    this_bool_map["write_hilltops"] = true;
    this_bool_map["write_hilltop_curvature"] = true;
    this_bool_map["write_hillslope_length"] = true;
    this_bool_map["write_hillslope_gradient"] = true;
    this_bool_map["write_hillslope_relief"] = true;
    this_bool_map["complete_hilltop_masking"] = true;
    this_bool_map["extract_ridges"] = true;
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

  //============================================================================
  //
  //.#####....####....####...######..######..#####...........#####...##..##..#####...##..##.
  //.##..##..##..##..##........##....##......##..##..........##..##..##..##..##..##..###.##.
  //.#####...######...####.....##....####....#####...........#####...##..##..#####...##.###.
  //.##..##..##..##......##....##....##......##..##..........##..##..##..##..##..##..##..##.
  //.##..##..##..##...####.....##....######..##..##..........#####....####...##..##..##..##.
  //
  //============================================================================
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


  //========================================================================
  //
  //.##.......####....####...#####...........#####....####...######...####..
  //.##......##..##..##..##..##..##..........##..##..##..##....##....##..##.
  //.##......##..##..######..##..##..........##..##..######....##....######.
  //.##......##..##..##..##..##..##..........##..##..##..##....##....##..##.
  //.######...####...##..##..#####...........#####...##..##....##....##..##.
  //
  //========================================================================
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
  //..####...##..##..#####...######...####....####...######..........######..######..######..######..######..##..##...####..
  //.##......##..##..##..##..##......##..##..##..##..##..............##........##......##......##......##....###.##..##.....
  //..####...##..##..#####...####....######..##......####............####......##......##......##......##....##.###..##.###.
  //.....##..##..##..##..##..##......##..##..##..##..##..............##........##......##......##......##....##..##..##..##.
  //..####....####...##..##..##......##..##...####...######..........##......######....##......##....######..##..##...####..
  //
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

  if( this_bool_map["use_spiral_trimmer_for_edge_influence"])
  {
    LSDRaster SpiralTrim = filled_topography.RasterTrimmerSpiral();
    filled_topography = SpiralTrim;
    cout << "Warning, I have replaced the edges of your DEM with nodata" << endl;
    cout << "This usually means I am trying to remove nodes influence by the edge" << endl;
    filled_topography.replace_edges_with_nodata();
  }

  if (this_bool_map["print_fill_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }

	static const int Arr[] = {0,1,1,1,1,0,0,0};
	vector<int> RasterSelection (Arr, Arr + sizeof(Arr) / sizeof(Arr[0]));
  cout << endl << endl << "=========================================" << endl;
  cout << "I am starting the analysis. This is very memory intensive!! " << endl;
  cout << "If you get a segmentation fault the most likeley problem is that your DEM is too large." << endl;
  cout << "To fix you can either clip the DEM, reduce the resolution, or buy a more expensive computer." << endl;
  cout << "Note that if you are using docker you can change the settings to give your container more memory." << endl;
	cout << "\tCalculating surface metrics..." << endl;
	vector<LSDRaster> Surfaces = filled_topography.calculate_polyfit_surface_metrics(this_float_map["surface_fitting_radius"], RasterSelection);
  LSDRaster Aspect = Surfaces[2];

  //Surface fitting to get slope, aspect, curvature and planform curvature
  if (this_bool_map["read_curvature_raster"])
  {
    cout << "I'm reading in a pre-calculated curvature raster." << endl;
    cout << "NOTE: This will overwrite the curvature calculated through polynomial surface fitting" << endl;
    LSDRaster curvature_raster((DATA_DIR+this_string_map["curvature_fname"]),raster_ext);
    if( this_bool_map["use_spiral_trimmer_for_edge_influence"])
    {
      LSDRaster SpiralTrim = curvature_raster.RasterTrimmerSpiral();
      curvature_raster = SpiralTrim;
      cout << "Warning, I have replaced the edges of your DEM with nodata" << endl;
      cout << "This usually means I am trying to remove nodes influence by the edge" << endl;
      curvature_raster.replace_edges_with_nodata();
    }
    Surfaces[3] = curvature_raster;
  }



  //============================================================================
  //
  //.######..##.......####...##...##..........#####....####...##..##..######..######..##..##...####..
  //.##......##......##..##..##...##..........##..##..##..##..##..##....##......##....###.##..##.....
  //.####....##......##..##..##.#.##..........#####...##..##..##..##....##......##....##.###..##.###.
  //.##......##......##..##..#######..........##..##..##..##..##..##....##......##....##..##..##..##.
  //.##......######...####....##.##...........##..##...####....####.....##....######..##..##...####..
  //
  //============================================================================
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

  // calculate the distance from outlet
  cout << "\t Calculating flow distance..." << endl;
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  cout << "\t Loading Sources..." << endl;
  cout << "\t Source file is... " << CHeads_file << endl;


  //==============================================================
  //
  //..####....####...##..##..#####....####...######...####..
  //.##......##..##..##..##..##..##..##..##..##......##.....
  //..####...##..##..##..##..#####...##......####.....####..
  //.....##..##..##..##..##..##..##..##..##..##..........##.
  //..####....####....####...##..##...####...######...####..
  //
  //==============================================================
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    cout << endl << endl << endl;
    cout << endl << "\t\t==================================" << endl;
    cout << "\t\tThe channel head file is NULL. " << endl;
    cout << "\t\tThis means I will calculate the sources." << endl;

    if (this_bool_map["print_wiener_channels"])
    {
      cout << "I am calculating channels using the Wiener algorithm (doi:10.1029/2012WR012452)." << endl;
      cout << "This algorithm was used by Clubb et al. (2016, DOI: 10.1002/2015JF003747) " << endl;
      cout << "and Grieve et al. (2016, doi:10.5194/esurf-4-627-2016) " << endl;
      cout << "and combines elements of the Pelletier and Passalacqua et al  methods: " << endl;
      cout << " doi:10.1029/2012WR012452 and doi:10.1029/2009JF001254" << endl;

      // initiate the spectral raster
      LSDRasterSpectral Spec_raster(topography_raster);

      string QQ_fname = OUT_DIR+OUT_ID+"__qq.txt";

      cout << "I am am getting the connected components using a weiner QQ filter." << endl;
      cout << "Area threshold is: " << this_float_map["pruning_drainage_area"] << " window is: " <<  this_float_map["surface_fitting_radius"] << endl;

      LSDIndexRaster connected_components = Spec_raster.IsolateChannelsWienerQQ(this_float_map["pruning_drainage_area"],
                                                        this_float_map["surface_fitting_radius"], QQ_fname);

      cout << "I am filtering by connected components" << endl;
      LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(this_int_map["connected_components_threshold"]);
      LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();

      cout << "I am thinning the network to a skeleton." << endl;
      LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();

      cout << "I am finding the finding end points" << endl;
      LSDIndexRaster Ends = skeleton_raster.find_end_points();
      Ends.remove_downstream_endpoints(CC_raster, Spec_raster);

      //this processes the end points to only keep the upper extent of the channel network
      cout << "getting channel heads" << endl;
      sources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);

    }
    if (this_bool_map["print_perona_malik_channels"])
    {
      cout << "I am calculating channels using a QQ threshold algorithm (doi:10.1029/2012WR012452)." << endl;
      cout << "and a Perona-Malik filter similar to GeoNet (Passalacqua et al, 2010): " << endl;
      cout << "doi:10.1029/2009JF001254" << endl;

      // initiate the spectral raster
      LSDRasterSpectral Spec_raster(topography_raster);

      string QQ_fname = OUT_DIR+OUT_ID+"__qq.txt";

      cout << "I am am getting the connected components using a weiner QQ filter." << endl;
      cout << "Area threshold is: " << this_float_map["pruning_drainage_area"] << " window is: " <<  this_float_map["surface_fitting_radius"] << endl;

      LSDIndexRaster connected_components = Spec_raster.IsolateChannelsPeronaMalikQQ(this_float_map["pruning_drainage_area"],
                                                         this_float_map["surface_fitting_radius"], QQ_fname);

      cout << "I am filtering by connected components" << endl;
      LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(this_int_map["connected_components_threshold"]);
      LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();

      cout << "I am thinning the network to a skeleton." << endl;
      LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();

      cout << "I am finding the finding end points" << endl;
      LSDIndexRaster Ends = skeleton_raster.find_end_points();
      Ends.remove_downstream_endpoints(CC_raster, Spec_raster);

      //this processes the end points to only keep the upper extent of the channel network
      cout << "getting channel heads" << endl;
      vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);

      cout << "got all the end points" << endl;

      // we need a temp junction network to search for single pixel channels
      LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
      LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

      cout << "removing single px channels" << endl;
      sources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);

      //Now we have the final channel heads, so we can generate a channel network from them
      //LSDJunctionNetwork ChanNetwork(FinalSources, FlowInfo);

    }
    else
    {
      cout << "\t\tI'm calculating sources using a threshold area routine." << endl;
      cout << "\t\tGetting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
      sources = FlowInfo.get_sources_index_threshold(FlowAcc, threshold_contributing_pixels);
    }
  }
  else
  {
    cout << endl << endl << endl;
    cout << endl << "\t\t==================================" << endl;
    cout << "\t\tI found a channel head filename. " << endl;
    cout << "\t\tLoading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    cout << "\t\tWarning: if you got the filename wrong this won't work." << endl;
    cout << "\t\tThe filename MUST INCLUDE THE csv extension." << endl;
    LSDRasterInfo ThisRI(filled_topography);
    string csv_filename = DATA_DIR+CHeads_file;
    LSDSpatialCSVReader CHeadCSV(ThisRI,csv_filename);
    sources = CHeadCSV.get_nodeindices_from_lat_long(FlowInfo);
    cout << "\t\tGot sources!" << endl;
  }

  cout << endl << endl << "================================" << endl;
  cout << "I've now ingested the channel source nodes. Yum yum. Sources are tasty." << endl;
  cout << "\t\tThe number of sources is: " << sources.size() << endl;
  if (sources.size() == 0)
  {
    cout << endl << endl << endl << "====================" << endl;
    cout << "FATAL ERROR: there are no sources. " << endl;
    cout << "The most likeley cause of this is that you tried to load a sources filename " << endl;
    cout << "That was incorrect. Please check the filename and try again." << endl;
    cout << "=================================" << endl;
    exit(0);
  }

  //=============================================================================================
  //  
  //.######..........##..##..######..######..##...##...####...#####...##..##.
  //.....##..........###.##..##........##....##...##..##..##..##..##..##.##..
  //.....##..........##.###..####......##....##.#.##..##..##..#####...####...
  //.##..##..........##..##..##........##....#######..##..##..##..##..##.##..
  //..####...........##..##..######....##.....##.##....####...##..##..##..##.
  //
  //=============================================================================================
  // now get the junction network
  cout << "\t I am now getting the junction network..." << endl;
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);
  cout << "\t Got the junction network" << endl;


  // Print channels and junctions if you want them.
  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "\t\tprinting channels to csv..." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
    if ( this_bool_map["use_extended_channel_data"])
    {
      cout << "I am going to use the extended channel network data outputs." << endl;
      JunctionNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, topography_raster);
    }
    else
    {
      JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
    }

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
    cout << "===================" << endl;
    cout << "SOURCES!!" << endl;
    cout << "Sources" << endl;
    cout << "===================" << endl;
    string sources_csv_name = OUT_DIR+OUT_ID+"_sources.csv";

    //write channel_heads to a csv file
    FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(sources, sources_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_sources.geojson";
      LSDSpatialCSVReader thiscsv(sources_csv_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

  // Print channel network to polyline
  if (this_bool_map["print_channels_to_line_vector"])
  {
    string channel_geojson_name = OUT_DIR+OUT_ID+"_PM_CN_line";
    cout << "I am going to print your channel network to a line geojson" << endl;
    // distance from outlet
    JunctionNetwork.PrintChannelNetworkToLineVector(FlowInfo, channel_geojson_name, filled_topography, DrainageArea, DistanceFromOutlet);
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

  //=================================================================================================================== 
  //
  //..####...######..#####...######...####...##...##..........##..##..######..######..##...##...####...#####...##..##.
  //.##........##....##..##..##......##..##..###.###..........###.##..##........##....##...##..##..##..##..##..##.##..
  //..####.....##....#####...####....######..##.#.##..........##.###..####......##....##.#.##..##..##..#####...####...
  //.....##....##....##..##..##......##..##..##...##..........##..##..##........##....#######..##..##..##..##..##.##..
  //..####.....##....##..##..######..##..##..##...##..........##..##..######....##.....##.##....####...##..##..##..##.
  //
  //=================================================================================================================== 
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
    cout << endl << endl << endl << "=============================" << endl;
    cout << "\tCombining the channelnetwork and floodplain masks..." << endl;
		LSDIndexRaster Floodplains((Floodplain_file), raster_ext);
		StreamNetwork.MergeIndexRasters(Floodplains);
	}



  //=================================================================================
  //.#####....####....####...######..##..##...####..
  //.##..##..##..##..##........##....###.##..##.....
  //.#####...######...####.....##....##.###...####..
  //.##..##..##..##......##....##....##..##......##.
  //.#####...##..##...####...######..##..##...####..
  //................................................ 
  //=================================================================================   
  //Check to see if a list of junctions for extraction exists
  // need to get base-level nodes , otherwise these catchments will be missed!
  vector< int > BaseLevelJunctions;
  vector< int > BaseLevelJunctions_Initial;
  
  cout << endl << endl << endl;
  cout << "======================================" << endl;
  cout << "\tNow I am going to deal with the baselevel junctions. " << endl;

  // deal with the baselevel junctions file
  string test_BaselevelJunctions_file = LSDPP.get_BaselevelJunctions_file();
  if(this_string_map["BaselevelJunctions_file"] == "NULL" && test_BaselevelJunctions_file == "NULL")
  {
    cout << "No baselevel junctions file found. I'll need to do a bit more checking to see what to do." << endl;
    BaselevelJunctions_file = "NULL";
  }
  else if(this_string_map["BaselevelJunctions_file"] == "NULL" && test_BaselevelJunctions_file != "NULL")
  {
    cout << "I am loading a baselevel junctions file." << endl;
    BaselevelJunctions_file = test_BaselevelJunctions_file;
  }
  else if(this_string_map["BaselevelJunctions_file"] != "NULL" && test_BaselevelJunctions_file == "NULL")
  {
    cout << "I am loading a baselevel junctions file." << endl;
    BaselevelJunctions_file = this_string_map["BaselevelJunctions_file"];
  }
  else
  {
    cout << "WARNING You have defined the baselevel junction file in two ways. " << endl;
    cout << "This is because the authors of LSDTopoTools created a dumb inheritance problem and can't fix it or it will break legacy code." << endl;
    cout << "I will use the newer version." << endl;
    BaselevelJunctions_file = this_string_map["BaselevelJunctions_file"];
    cout << "The junctions file I am using is: " <<  BaselevelJunctions_file << endl;
  }
  // now check to see if there is a full path
  cout << endl << endl << "I need to check your baselevel junctions file, to see if it is in the correct path. " << endl;
  BaselevelJunctions_file = LSDPP.check_for_path_and_add_read_path_if_required(BaselevelJunctions_file);


  cout << "\n\nThe BaselevelJunctions_file is: " << BaselevelJunctions_file << endl;
  if (BaselevelJunctions_file == "NULL" || BaselevelJunctions_file == "Null" || BaselevelJunctions_file == "null" || BaselevelJunctions_file.empty() == true)
  {
    cout << "\tI didn't find a baselevel junctions file. Checking to see if I use lat-long." << endl;  
    if(this_bool_map["get_basins_from_outlets"])
    {
      cout << "I am going to get basins lat-long coordinates" << endl;
      string full_BL_LL_name = DATA_DIR+this_string_map["basin_outlet_csv"];
      cout << "The file is: " << full_BL_LL_name << endl;
      int search_radius_nodes = this_int_map["search_radius_nodes"];
      int threshold_stream_order = this_int_map["threshold_stream_order_for_snapping"];
      BaseLevelJunctions = JunctionNetwork.snap_point_locations_to_upstream_junctions_from_latlong_csv(full_BL_LL_name,
                                                          search_radius_nodes, threshold_stream_order,FlowInfo, RI);

      BaseLevelJunctions_Initial = BaseLevelJunctions;

      cout << "My preliminary list of baselevel junctions is:" << endl;
      for (int i = 0; i<int(BaseLevelJunctions.size()); i++)
      {
        cout << BaseLevelJunctions[i] << endl;
      }
    }
    else
    {
      //Get baselevel junction nodes from the whole network
      cout << "I am going to get baselevel junctions using an algorithm." << endl;
      BaseLevelJunctions_Initial = JunctionNetwork.get_BaseLevelJunctions();

      // now prune these by drainage area
      cout << "\tRemoving basins with fewer than " << this_int_map["minimum_basin_size_pixels"] << " pixels" << endl;
      //BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions_Initial,FlowInfo, FlowAcc, this_int_map["minimum_basin_size_pixels"]);
      cout << "The pixel limits are: lower: " << this_int_map["minimum_basin_size_pixels"] << " and upper: " << this_int_map["maximum_basin_size_pixels"] << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(FlowInfo, filled_topography, FlowAcc,
                                                      this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
      cout << "\tNow I have " << BaseLevelJunctions.size() << " baselelvel junctions." << endl;

      BaseLevelJunctions_Initial = BaseLevelJunctions;
    }

  }
  else
  {
    cout << "I am attempting to read base level junctions from a base level junction list." << endl;
    cout << "If this is not a simple text file that only contains integers there will be problems!" << endl;

    //specify junctions to work on from a list file
    //string JunctionsFile = DATA_DIR+BaselevelJunctions_file;
    cout << "The junctions file is: " << BaselevelJunctions_file << endl;

    vector<int> JunctionsList;
    ifstream infile(BaselevelJunctions_file.c_str());
    if (infile)
    {
      cout << "Junctions File " << BaselevelJunctions_file << " exists" << endl;;
      int n;
      while (infile >> n) BaseLevelJunctions_Initial.push_back(n);
    }
    else
    {
      cout << "Fatal Error: Junctions File " << BaselevelJunctions_file << " does not exist" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Now make sure none of the basins drain to the edge
  cout << "\tPruning junctions that drain to the edge of the DEM..." << endl;
  BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions_Initial, FlowInfo, filled_topography);
  
  cout << "I have finished with the baselevel junctions. They are: " << endl;
  for (int i = 0; i<int(BaseLevelJunctions.size()); i++)
  {
    cout << BaseLevelJunctions[i] << "  ";
    if (i%10==0)
    {
      cout << endl;
    }
  }
  cout << endl << "======================================" << endl;
  cout << endl << endl << endl;


  // Now we get channel segments for use with chi plotting.
  vector<int> source_nodes;
  vector<int> outlet_nodes;
  vector<int> baselevel_node_of_each_basin;
  if (this_bool_map["extend_channel_to_node_before_receiver_junction"])
  {
    cout << endl << endl << "=====================================================" << endl;
    cout << "I am now getting the channels by basin." << endl;
    cout << "  These channels extend below the junction to the channel that stops" << endl;
    cout << "  just before the receiver junction. This option is used to remain" << endl;
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
    cout << "  just before the receiver junction, then set the option:" << endl;
    cout << "    extend_channel_to_node_before_receiver_junction" << endl;
    cout << "  to true." << endl;
    cout << "=====================================================" << endl << endl;

    JunctionNetwork.get_overlapping_channels(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);
  }

	//get the basins
	cout << "\tExtracting basins..." << endl;
	LSDIndexRaster BasinsRaster = JunctionNetwork.extract_basins_from_junction_vector(BaseLevelJunctions, FlowInfo);


  //=====================================================================================================
  //
  //.#####...#####...######..##..##..######..........#####....####....####...######..######..#####....####..
  //.##..##..##..##....##....###.##....##............##..##..##..##..##........##....##......##..##..##.....
  //.#####...#####.....##....##.###....##............#####...######...####.....##....####....#####....####..
  //.##......##..##....##....##..##....##............##..##..##..##......##....##....##......##..##......##.
  //.##......##..##..######..##..##....##............##..##..##..##...####.....##....######..##..##...####..
  //
  //=====================================================================================================
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
    string hillshade_raster_name = OUT_DIR+OUT_ID+"_HS";
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











  //=====================================================================================================
  //
  //.######..##..##..######..#####....####....####...######...........####...##..##..######.
  //.##.......####.....##....##..##..##..##..##..##....##............##..##..##..##....##...
  //.####......##......##....#####...######..##........##............##......######....##...
  //.##.......####.....##....##..##..##..##..##..##....##............##..##..##..##....##...
  //.######..##..##....##....##..##..##..##...####.....##.............####...##..##..######.
  //
  //=====================================================================================================
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

  
  //=====================================================================================================
  //
  //.######..##..##..######..#####....####....####...######..........#####...######..#####....####...######...####..
  //.##.......####.....##....##..##..##..##..##..##....##............##..##....##....##..##..##......##......##.....
  //.####......##......##....#####...######..##........##............#####.....##....##..##..##.###..####.....####..
  //.##.......####.....##....##..##..##..##..##..##....##............##..##....##....##..##..##..##..##..........##.
  //.######..##..##....##....##..##..##..##...####.....##............##..##..######..#####....####...######...####..
  //
  //=====================================================================================================
  if ( this_bool_map["run_HFR_analysis"] ||
       this_bool_map["extract_ridges"])
  {
    cout << "\tExtracting hilltop network the new way" << endl;
    map<int,int> hilltop_nodes = JunctionNetwork.ExtractAllRidges(FlowInfo);

    map<int,int> filtered_hilltop_nodes;
    // remove hilltops near to the channel
    if (this_int_map["Threshold_Distance_From_Channel"] > 0)
    {
      filtered_hilltop_nodes = JunctionNetwork.filter_nodes_by_distance_to_channel(hilltop_nodes, this_int_map["Threshold_Distance_From_Channel"], FlowInfo, DistanceFromOutlet, this_int_map["Threshold_Stream_Order"]);
    }
    else { filtered_hilltop_nodes = hilltop_nodes; }

    map<int, vector<int> > full_hilltop_int_map;
    map<int, vector<float> > full_hilltop_float_map;
    map<int, vector<double> > full_hilltop_latlong_map;

    ofstream ht_csv_out;
    string ht_csv_full_fname = OUT_DIR+OUT_ID+"_RidgeData.csv";
    ht_csv_out.open(ht_csv_full_fname.c_str());
    ht_csv_out.precision(9);
    ht_csv_out << "latitude,longitude,row,col,basin_id,stream_order,contributing_pixels,slope,curvature"<<endl;

    map<int, int>::iterator it;
    LSDCoordinateConverterLLandUTM Converter;
    double latitude, longitude;
    int curr_row,curr_col;
    int this_node;
    int this_SO;
    float slope, curv;
    int this_CP;
    int this_basin;
    
    vector<int> hilltop_int_vec(5,0);
    vector<float> hilltop_float_vec(2,0);
    vector<double> hilltop_latlong_vec(2,0);
    for (it = filtered_hilltop_nodes.begin(); it != filtered_hilltop_nodes.end(); it++)  
    {
      this_node = it->first;
      this_SO = it->second;

      //cout << "Node is: " << this_node << endl;

      FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);
      slope = Surfaces[1].get_data_element(curr_row,curr_col);
      curv = Surfaces[3].get_data_element(curr_row,curr_col);
      this_CP = FlowInfo.retrieve_contributing_pixels_of_node(this_node);
      this_basin = BasinsRaster.get_data_element(curr_row,curr_col);

      // check if the hilltop pixel meets all the masking criteria
      if (slope < this_float_map["Threshold_Hilltop_Gradient"] && curv < this_float_map["Threshold_Hilltop_Curvature"] &&
          this_CP <= this_int_map["Threshold_Hilltop_Contributing_Pixels"] && this_basin >= 0)
      {       
        hilltop_int_vec[0] = curr_row;
        hilltop_int_vec[1] = curr_col;
        hilltop_int_vec[2] = this_basin;
        hilltop_int_vec[3] = this_SO;
        hilltop_int_vec[4] = this_CP;

        hilltop_float_vec[0] = slope;
        hilltop_float_vec[1] = curv;

        FlowInfo.get_lat_and_long_from_current_node(this_node, latitude,longitude, Converter);
        hilltop_latlong_vec[0] = latitude;
        hilltop_latlong_vec[1] = longitude;

        full_hilltop_int_map[this_node] = hilltop_int_vec;
        full_hilltop_float_map[this_node] = hilltop_float_vec;
        full_hilltop_latlong_map[this_node] = hilltop_latlong_vec;

        ht_csv_out << latitude << "," << longitude << "," << curr_row << "," << curr_col 
                << "," << this_basin << "," << this_SO << "," << this_CP 
                << "," << slope << "," << curv << endl;

      }
    } 
    ht_csv_out.close(); 

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      cout << "Now let me print your ridge data geojson" << endl;
      string gjson_name = OUT_DIR+OUT_ID+"_RidgeData.geojson";
      LSDSpatialCSVReader thiscsv(ht_csv_full_fname);
      thiscsv.print_data_to_geojson(gjson_name);
    }

  }

  //=========================================
  //
  //.##..##..######..#####..
  //.##..##..##......##..##.
  //.######..####....#####..
  //.##..##..##......##..##.
  //.##..##..##......##..##.
  //
  //==========================================
  if (this_bool_map["run_HFR_analysis"])
  {
    // Run hilltop flow routing
    cout << "\tRunning hillslope flow routing. You might want to find a good book..." << endl;
    string HillslopeTracesFolder = DATA_DIR;


    // Do the legacy code first. 
    if (this_bool_map["use_legacy_HFR"])
    {
      cout << "This is the legacy version of the hilltop flow routing. I first need to do some preprocessing"<< endl;
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

      vector<int> Target_Basin_Vector;   // this doesn't actually do anything
      cout << "NOTE I am using the legacy hilltop flow routing." << endl;
      auto t1 = std::chrono::high_resolution_clock::now();
      vector< Array2D<float> >HFR_Arrays = FlowInfo.HilltopFlowRouting_Refactored(filled_topography, Hilltops, 
                                                            Surfaces[1], Surfaces[2], Surfaces[3], Surfaces[4], 
                                                            StreamNetwork, BasinsRaster, (DATA_DIR+DEM_ID), 
                                                            this_bool_map["print_hillslope_traces"], this_int_map["hillslope_trace_thinning"], 
                                                            HillslopeTracesFolder, this_bool_map["hillslope_traces_basin_filter"], 
                                                            Target_Basin_Vector);    
      auto t2 = std::chrono::high_resolution_clock::now();   
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();                                                                                                     
      cout << "\tDone. That took " << duration << " ms. " << endl << endl;

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
    else            // Everything below here is the new version of the code
    {
      cout << "Calculating hilltop flow tracing." << endl;
      cout << "I am using the new and refactored version of the routing code." << endl;
      auto t1 = std::chrono::high_resolution_clock::now();
      string ht_csv_full_fname = OUT_DIR+OUT_ID+"_RidgeData.csv";   // name of the ridgeline file
      FlowInfo.HilltopFlowRouting_TerrifyingRefactored(filled_topography, Surfaces[2], Surfaces[4], 
                                                            StreamNetwork, ht_csv_full_fname, (OUT_DIR+OUT_ID), 
                                                            this_bool_map["print_hillslope_traces"], this_int_map["hillslope_trace_thinning"], 
                                                            HillslopeTracesFolder, this_bool_map["hillslope_traces_basin_filter"]);    

      auto t2 = std::chrono::high_resolution_clock::now();   
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();                                                                                                     
      cout << "\tDone. That took " << duration << " ms. " << endl << endl;

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        cout << "Now let me print your ridge data geojson" << endl;
        string gjson_name = OUT_DIR+OUT_ID+"_HilltopData_TN.geojson";
        string csv_name = OUT_DIR+OUT_ID+"_HilltopData_TN.csv";
        LSDSpatialCSVReader thiscsv(csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
  }   // End HFR logic

  //=====================================================================================================
  // 
  //.##..##..######..##......##......######...####...#####............####............####...##..##...####...##..##.
  //.##..##....##....##......##........##....##..##..##..##..............##..........##..##..##..##..##..##..###.##.
  //.######....##....##......##........##....##..##..#####............####...........##......######..######..##.###.
  //.##..##....##....##......##........##....##..##..##..............##..............##..##..##..##..##..##..##..##.
  //.##..##..######..######..######....##.....####...##..............######...........####...##..##..##..##..##..##.
  //  
  //=====================================================================================================
  if (this_bool_map["map_hilltops_to_channels"])
  {
    cout << "I am going to map the hilltops to your segmented channel network" << endl;
    // first get the two csv files
    string csv_name = OUT_DIR+OUT_ID+"_HilltopData_TN.csv";
    LSDSpatialCSVReader hilltop_csv(csv_name);

    csv_name = OUT_DIR+OUT_ID+"_MChiSegmented.csv";
    LSDSpatialCSVReader segments_csv(csv_name);

    // now create a data map from the segments
    map<int, int> segment_data_map;

    vector<int> segment_node_vec = segments_csv.data_column_to_int("node");
    vector<float> segment_chi_vec = segments_csv.data_column_to_float("chi");
    vector<float> segment_z_vec = segments_csv.data_column_to_float("elevation");
    vector<float> segment_fd_vec = segments_csv.data_column_to_float("flow_distance");
    vector<float> segment_da_vec = segments_csv.data_column_to_float("drainage_area");
    vector<float> segment_m_chi_vec = segments_csv.data_column_to_float("m_chi");
    vector<float> segment_b_chi_vec = segments_csv.data_column_to_float("b_chi");
    vector<int> segment_basin_key_vec = segments_csv.data_column_to_int("basin_key");
    vector<int> segment_source_key_vec = segments_csv.data_column_to_int("source_key");

    int n_channel_nodes = int(segment_node_vec.size());
    for(int i = 0; i<n_channel_nodes; i++)
    {
      segment_data_map[ segment_node_vec[i] ] = i;
    }

    // now loop through the hilltop nodes, mapping the hilltops to channel data  
    vector<int> ht_row_vec = hilltop_csv.data_column_to_int("row");
    vector<int> ht_col_vec = hilltop_csv.data_column_to_int("col");
    vector<int> chan_row_vec = hilltop_csv.data_column_to_int("channel_row");
    vector<int> chan_col_vec = hilltop_csv.data_column_to_int("channel_col");

    vector<int> basin_vec = hilltop_csv.data_column_to_int("BasinID");
    vector<int> stream_vec = hilltop_csv.data_column_to_int("StreamID");
    vector<int> planar_count_vec = hilltop_csv.data_column_to_int("PlanarCountFlag");
    vector<int> divergent_count_vec = hilltop_csv.data_column_to_int("DivergentCount");

    vector<float> ht_slope_vec = hilltop_csv.data_column_to_float("S");
    vector<float> ht_curv_vec = hilltop_csv.data_column_to_float("Cht");
    vector<float> ht_R_vec = hilltop_csv.data_column_to_float("R");
    vector<float> ht_Lh_vec = hilltop_csv.data_column_to_float("Lh");
    vector<float> ht_easting_vec = hilltop_csv.data_column_to_float("easting");
    vector<float> ht_northing_vec = hilltop_csv.data_column_to_float("northing");
    vector<float> ht_hilltopS_vec = hilltop_csv.data_column_to_float("HilltopSlope");
    vector<float> ht_ES_vec = hilltop_csv.data_column_to_float("E_Star");
    vector<float> ht_RS_vec = hilltop_csv.data_column_to_float("R_Star");
    vector<float> ht_ED_vec = hilltop_csv.data_column_to_float("EucDist");

    vector<double> ht_lat_vec = hilltop_csv.get_latitude();
    vector<double> ht_long_vec = hilltop_csv.get_longitude();   

    ofstream ofs;

    //create the output filename from the user supplied filename prefix
    string ss_filename = OUT_DIR+OUT_ID+"_combined_HT_chan.csv";

    ofs.open(ss_filename.c_str());

    if( ofs.fail() )
    {
      cout << "\nFATAL ERROR: unable to write to " << ss_filename.c_str() << endl;
      exit(EXIT_FAILURE);
    }
    ofs << "latitude,longitude,easting,northing,row,col,Cht,S,R,Lh,BasinID,channel_row,channel_col," 
        << "StreamID,HilltopSlope,DivergentCount,PlanarCountFlag,E_Star,R_Star,EucDist,"
        << "chan_node,chi,chan_elev,flow_distance,drainage_area,m_chi,b_chi,basin_key,source_key\n";



    int this_chan_ni;
    int no_link_count = 0;
    int link_count = 0;
    int cn;
    int n_ht = int( ht_row_vec.size());
    for (int htn = 0; htn < n_ht; htn++)
    {
      // get the nodeindex, then search for the key
      this_chan_ni = FlowInfo.retrieve_node_from_row_and_column(chan_row_vec[htn],chan_col_vec[htn]);
      //cout << "This channel node: " << this_chan_ni << endl;

      if ( segment_data_map.find(this_chan_ni) == segment_data_map.end() )
      {
        // not found
        no_link_count++;
      }
      else
      {
        link_count++;
        cn = segment_data_map[this_chan_ni];
              // this would be much easier in python but we just reprint everything. 
        ofs.precision(9);
        ofs << ht_lat_vec[htn] << "," << ht_long_vec[htn] << ",";
        ofs.precision(6);
        ofs << ht_easting_vec[htn] << "," << ht_northing_vec[htn] << ","
          << ht_row_vec[htn] << "," << ht_col_vec[htn] << "," << ht_curv_vec[htn] << "," << ht_slope_vec[htn] << ","
          << ht_R_vec[htn] << "," << ht_Lh_vec[htn] << "," << basin_vec[htn] << "," << chan_row_vec[htn] << "," << chan_col_vec[htn] << ","
          << stream_vec[htn] << "," << ht_hilltopS_vec[htn] << "," << divergent_count_vec[htn] << "," << planar_count_vec[htn] << ","
          << ht_ES_vec[htn] << "," << ht_RS_vec[htn] << "," << ht_ED_vec[htn] << ",";
        ofs << segment_node_vec[cn] << "," << segment_chi_vec[cn] << "," << segment_z_vec[cn] << ","
            << segment_fd_vec[cn] << "," << segment_da_vec[cn] << "," << segment_m_chi_vec[cn] << ","
            << segment_b_chi_vec[cn] << "," << segment_basin_key_vec[cn] << "," << segment_source_key_vec[cn] << endl;
      }
    }
    cout << "I was able to link " << link_count << " hilltops to channels, and failed to link " << no_link_count << " hilltops." << endl;

    ofs.close();


    if ( this_bool_map["convert_csv_to_geojson"])
    {
      cout << "Now let me print your ridge data geojson" << endl;
      string gjson_name = OUT_DIR+OUT_ID+"_combined_HT_chan.geojson";
      string csv_name = OUT_DIR+OUT_ID+"_combined_HT_chan.csv";
      LSDSpatialCSVReader combinedcsv(csv_name);
      combinedcsv.print_data_to_geojson(gjson_name);
    }
  }
}
