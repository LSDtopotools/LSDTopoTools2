//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// lsdtt-valley-metrics
// A series of tools to analyse floodplains and terraces and get valley widths
//
// These include the floodplain and terrace algorithms developed for the paper:
// https://www.earth-surf-dynam.net/5/369/2017/
//
// and valley extraction algorithms in a forthcoming publication
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Copyright (C) 2021 Fiona J. Clubb and Simon M. Mudd 2021
//
// Developers can be contacted:
//
//    fiona.j.clubb _at_ durham.ac.uk
//    Fiona Clubb
//    Durham University
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
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDFloodplain.hpp"
#include "../LSDTerrace.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDSpatialCSVReader.hpp"


int main (int nNumberofArgs,char *argv[])
{
  //start the clock
  clock_t begin = clock();

  string version_number = "0.5";
  string citation = "http://doi.org/10.5281/zenodo.4577879";

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the LSDTopoTools valley metrics tool!    ||" << endl;
  cout << "|| This program extracts floodplains and               ||" << endl;
  cout << "|| terraces using slope and relief_thresholds.         ||" << endl;
  cout << "|| This program was developed by                       ||" << endl;
  cout << "|| Fiona J. Clubb, Durham University, and              ||" << endl;
  cout << "|| Simon M Mudd at the University of Edinburgh         ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| Citation for this code is:                          ||" << endl;
  cout << "|| " << citation << endl;
  cout << "|| And                                                 ||" << endl;
  cout << "|| https://www.earth-surf-dynam.net/5/369/2017/        ||" << endl;
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
    ofs.open("./lsdtt-valley-metrics-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    cout << endl << endl << endl << "==============================================" << endl;    
    cout << "This is lsdtt-valley-metrics version number " << version_number << endl;
    cout << "If the version contains a 'd' then you are using a development version." << endl;
    cout << "=========================================================" << endl;
    ofstream ofs;
    ofs.open("./lsdtt-valley-metrics-version.txt");
    ofs << version_number << endl;
    ofs.close();

    exit(0);
  }


  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // this will contain the help file
  map< string, vector<string> > help_map;

  // single channel
  string_default_map["channel_source_fname"] = "NULL";
  help_map["channel_source_fname"] = {  "string","NULL","Name of the csv file without extension of the single channel source. Needs column headers latitude and longitude.","If latitude and longitude are not column headers in the csv this will not work."};

  bool_default_map["extract_single_channel"] = false;
  help_map["extract_single_channel"] = {  "bool","false","This extracts a flow path from a line that starts at a point denoted in the file single_channel_source.","Used for imposing baselevel in various simulations."};

  bool_default_map["extract_multiple_channels"] = false;
  help_map["extract_multiple_channels"] = {  "bool","false","This extracts multiple flow path from upstream sources from the csv channel_source_fname.","Used for imposing baselevel in various simulations."};

  bool_default_map["extract_all_channels"] = false;
  help_map["extract_all_channels"] = {  "bool","false","This extracts the whole channel network. Use to get a width network for the whole DEM.","Used for imposing baselevel in various simulations."};

  int_default_map["search_radius_nodes"] = 8;
  help_map["search_radius_nodes"] = {  "int","8","A parameter for snapping to the nearest channel. It will search for the largest channel (by stream order) within the pixel window.","You will want smaller pixel numbers if you have a dense channel network."};
 
  int_default_map["threshold_stream_order_for_snapping"] = 2;
  help_map["threshold_stream_order_for_snapping"] = {  "int","2","If you are snapping to a channel, it will ignore channel with lower stream order than this number.","Set this to a higher number to avoid snapping to small channels."};
  

  // valley centreline
  bool_default_map["get_valley_centreline"] = false;
  help_map["get_valley_centreline"] = {  "bool","false","This gets a valley centreline.","."};

  int_default_map["centreline_loops"] = 5;
  help_map["centreline_loops"] = {  "int","5","The centreline routine iteratively digs and fills the trough and this sets the number of iterations.","We find 5 is about right but if the results are looking bad then try increasing this number. Warning: that might make things worse."};
 
  float_default_map["trough_scaling_factor"] = 0.1;
  help_map["trough_scaling_factor"] = {  "float","0.1","The centreline routine digs a trough (to route flow through the centre of the valley) with a depth that is this factor times the distance from the edge. Increase this number for a deeper trough.","A deeper trough has more sucess at routing flow but you get backwater effects at the bottom of the valley. Play with this number to get the best result."};
  
  float_default_map["minimum_bank_elevation_window_radius"] = 20;
  help_map["minimum_bank_elevation_window_radius"] = {  "float","20","The valley routine will dig a trough where the starting elevation at the side is determined by the lowest bank elevation within this radius.","For best results this needs to be about the width of the valley."};
   

  // filling
  bool_default_map["carve_before_fill"] = false; // Implements a carving algorithm
  help_map["carve_before_fill"] = {  "bool","false","This implements a breaching algorithm before filling.","Good for landscapes with DEM obstructions (like roads) across the channels."};
 
  float_default_map["min_slope_for_fill"] = 0.0001;
  help_map["min_slope_for_fill"] = {  "float","0.0001","Minimum slope between pixels for the filling algorithm.","Best not to change the default."};

  bool_default_map["remove_seas"] = false;
  help_map["remove_seas"] = {  "bool","true","Slightly misleading name; it replaces both high and low DEM values with nodata.","This gets rid of low lying areas but also is handy when the nodata is not translated from the raw DEM and it is full of funny large numbers."};

  // valley and floodplain extraction
  int_default_map["threshold_SO"] = 3;
  help_map["threshold_SO"] = {  "int","3","Threshold stream order for extracting the floodplain.","Seems to still get lower order channels so check with Fiona."};

  int_default_map["relief_lower_percentile"] = 25;
  help_map["relief_lower_percentile"] = {  "int","25","The lower relief percentile that is used to fit the QQ plot.","See Clubb et al ESURF 2017."};

  int_default_map["relief_upper_percentile"] = 75;
  help_map["relief_upper_percentile"] = {  "int","75","The upper relief percentile that is used to fit the QQ plot.","See Clubb et al ESURF 2017."};

  int_default_map["slope_lower_percentile"] = 25;
  help_map["slope_lower_percentile"] = {  "int","25","The lower slope percentile that is used to fit the QQ plot.","See Clubb et al ESURF 2017."};

  int_default_map["slope_upper_percentile"] = 75;
  help_map["slope_upper_percentile"] = {  "int","75","The upper slope percentile that is used to fit the QQ plot.","See Clubb et al ESURF 2017."};

  int_default_map["minimum_patch_size"] = 1000;
  help_map["minimum_patch_size"] = {  "int","1000","The minimum number of pixels in a patch that will be classed as a floodplain or terrace.","If this is low you will get more floodplains but they may be disconnected."};

  int_default_map["threshold_contributing_pixels"] = 1000;
  help_map["threshold_contributing_pixels"] = {  "int","1000","The number of contributing pixels needed to start a channel using the threshold method.","This is in pixels not drainage area. More options are in the lsdtt-channel-extraction tool."};

  float_default_map["surface_fitting_radius"] = 30;
  help_map["surface_fitting_radius"] = {  "float","30","Our surface fitting routines fit a polynomial over the points with in a radius defined by surface_fitting_radius and then differentiate this surface to get the surface metrics like gradient and curvature","If not bigger than the pixel_size*sqrt(2) then will increase to that number."};

  float_default_map["QQ_threshold"] = 0.005;
  help_map["QQ_threshold"] = {  "float","0.005","Once the QQ plot is fitted, this determines the level below which a floodplain is identified. Higher numbers give you more floodplain or terrace pixels.","See Clubb et al ESURF 2017."};


  // choice for absolute thresholds
  bool_default_map["use_absolute_thresholds"] = false;
  help_map["use_absolute_thresholds"] = {  "bool","false","This overrides the quantile-qunatile approach and just sets anything below the relief and slope thresholds as a floodplain.","This is not recommended but might work in very high relief landscapes."};
 
  int_default_map["relief_threshold"] = 50;
  help_map["relief_threshold"] = {  "int","50","When use_absolute_thresholds is true this is the threshold below which you get a floodplain","The slope threshold must also be lower to trigger a floodplain"};
  
  float_default_map["slope_threshold"] = 0.2;
  help_map["slope_threshold"] = {  "int","50","When use_absolute_thresholds is true this is the threshold below which you get a floodplain","The relief threshold must also be lower to trigger a floodplain"};
 
  // valley width
  bool_default_map["fill_floodplain"] = false;
  help_map["fill_floodplain"] = {  "bool","false","This looks for holes in the floodplain and then fills them to get a continuous floodplain mask","Needed to make sure holes in the floodplain don't truncate the channel width."};
   
  int_default_map["fill_loops"] = 1;
  help_map["fill_loops"] = {  "int","1","For the floodplain filling routine this controls the number of loops that are tried to fill the floodplain.","Needed to make sure holes in the floodplain don't truncate the channel width."};
 
  int_default_map["fp_search_radius"] = 2;
  help_map["fp_search_radius"] = {  "int","2","For the floodplain filling routine this controls how far the routine looks for holes.","Needed to make sure holes in the floodplain don't truncate the channel width."};
 
  bool_default_map["get_valley_widths"] = false;
  help_map["get_valley_widths"] = {  "bool","false","Runs the valley width routines. It uses the floodplain mask to get the valley widths so will be sensitive to floodplain settings.","Will result in valley width extraction."};
    
  int_default_map["width_node_spacing"]= 2;
  help_map["width_node_spacing"] = {  "int","2","Width is measured along valley centreline nodes and this sets how frequently width is sampled. 1 means using every node.","On lidar you probably want a relatively high number (e.g. 10) whereas coarse data you might want something small."};
 
  bool_default_map["print_channel_bearings"] = false;
  help_map["print_channel_bearings"] = {  "bool","false","This prints to csv the directional bearings (looking downslope) of channel pixels.","Helps to bug check and visualise the direction of the channel where width is measured orthogonal to the bearing."};
 
  int_default_map["valley_banks_search_radius"] = 5;
  help_map["valley_banks_search_radius"] = {  "int","5","A radius to look for valley edges.","Need to check exactly what this does later--SMM."};
 
  string_default_map["valley_points_csv"]  = "NULL";
  help_map["valley_points_csv"] = {  "string","NULL","If this isn't NULL it will load a valley point csv that need latitude and longitude and flow distance columns.","If this is NULL the valley centreline will be calculated algorithmically."};
   

  bool_default_map["remove_tributary_widths"] = false;
  help_map["remove_tributary_widths"] = {  "bool","false","This takes away widths that are around tributaries.","Ensures you don't get wide widths where the widths go up into the tributary valleys."};
 
  int_default_map["tributary_search_radius"] = 5;
  help_map["tributary_search_radius"] = {  "int","5","How far along the channel pixels to check if the width intersects with a tributary.","Set remove_tributary_widths to true to enable this."};
  
  int_default_map["tributary_stream_order"] = 3;
  help_map["tributary_stream_order"] = {  "int","3","The tributary search only looks for tributaries this order or more.","If this is low you will get a lot of little tributaries that don't actually modify the channel width eliminated."};
 

  // terraces
  bool_default_map["get_terraces"] = false;
  help_map["get_terraces"] = {  "bool","false","Finds terraces","Refer to Clubb et al 2017 ESURF"};
 
  int_default_map["minimum_terrace_height"] = -10;
  help_map["minimum_terrace_height"] = {  "int","-10","The minimum height above the channel for detecting a terrace","Default is low to get all terraces."};

  float_default_map["swath_width"] = 1000;
  help_map["swath_width"] = {  "float","1000","Width in metres of the swath along the channel that is used to find terraces","Make this wider if you want terraces far from the channel."};

  bool_default_map["print_swath_raster"] = false;
  help_map["print_swath_raster"] = {  "bool","false","Prints the raster of the swath used to find terraces","Good for bug checking to make sure the correct river is being sampled"};

  // printing
  bool_default_map["write_hillshade"] = false;
  help_map["write_hillshade"] = {  "bool","false","Write the hillshade raster.","You need this for a lot of our plotting routines. Filename includes _HS"};

  bool_default_map["print_channels_to_csv"] = false;
  help_map["print_channels_to_csv"] = {  "bool","false","Prints the channel network to a csv file.","This version produces smaller files than the raster version."};

  bool_default_map["print_junctions_to_csv"] = false;
  help_map["print_junctions_to_csv"] = {  "bool","false","Prints a csv with the locations and numbers of the junctions.","This is better to use than the raster version."};
 
  bool_default_map["use_extended_channel_data"] = false;
  help_map["use_extended_channel_data"] = {  "bool","false","If this is true you get more data columns in your channel network csv.","I will tell you what these columns are one day."};

  bool_default_map["convert_csv_to_geojson"] = false;
  help_map["convert_csv_to_geojson"] = {  "bool","false","Converts csv files to geojson files","Makes csv output easier to read with a GIS. Warning: these files are much bigger than csv files."};


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
  //============================================================================
  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);
  LSDPP.force_bil_extension();

  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // Parses the help file
  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./lsdtt-basic-metrics-README.csv" << endl;
    string help_prefix = "lsdtt-valley-metrics-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string DEM_extension =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();


  cout << "I've parsed all your parameters...now I'm going to get you some FLOODPLAINS :)" << endl;

  bool need_filtered_DEM = true;

  LSDRaster filled_topography;
  LSDRaster topography_raster((DATA_DIR+DEM_ID), DEM_extension);
  LSDRasterInfo RI(topography_raster);

  // Hillshade option
  if(this_bool_map["write_hillshade"])
  {
    cout << "Let me print the hillshade for you. " << endl;
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+OUT_ID+"_hs";
    hs_raster.write_raster(hs_fname,DEM_extension);
  }


  //=========================================================================
  //
  // .%%%%%%..%%%%%%..%%......%%%%%%..%%%%%%..%%%%%...%%%%%%..%%..%%...%%%%.
  // .%%........%%....%%........%%....%%......%%..%%....%%....%%%.%%..%%....
  // .%%%%......%%....%%........%%....%%%%....%%%%%.....%%....%%.%%%..%%.%%%.
  // .%%........%%....%%........%%....%%......%%..%%....%%....%%..%%..%%..%%.
  // .%%......%%%%%%..%%%%%%....%%....%%%%%%..%%..%%..%%%%%%..%%..%%...%%%%..
  //
  //=========================================================================
  // First we need to test for filtering
  cout << "====================================================" << endl;
  cout << " I need to check for filtering. " << endl;
  string filtered_prefix = DATA_DIR+DEM_ID+"_filtered";
  string filtered_test_file = DATA_DIR+DEM_ID+"_filtered.hdr";

  ifstream test_filter(filtered_test_file.c_str());
  if (test_filter)
  {
    cout << "I found a filtered file. I am going to use that as my DEM." << endl;
    LSDRaster load_DEM(filtered_prefix, DEM_extension);
    filled_topography = load_DEM;
  }
  else
  {
    cout << "I did not find a filtered DEM." << endl;
    cout << "For these routines I need to performa a pre-procesing step" << endl;
    cout << "where I use a perona-malik filter on your topographic data." << endl;
    cout << "I'm starting that now." << endl;

    // load the DEM
    cout << "Loading the DEM..." << endl;

    filled_topography = topography_raster;

    if(this_bool_map["remove_seas"])
    {
      filled_topography.remove_seas();
    }

    // filter using Perona Malik
    int timesteps = 50;
    float percentile_for_lambda = 90;
    float dt = 0.1;
    cout << "Now for the filtering step. This might take a little while." << endl;
    filled_topography = filled_topography.PeronaMalikFilter(timesteps, percentile_for_lambda, dt);

    // fill
    cout << "Let me fill that raster for you, the min slope is: "
        << this_float_map["min_slope_for_fill"] << endl;
    if(this_bool_map["carve_before_fill"])
    {
      cout << "I am using the carving algorithm before filling." << endl;
      LSDRaster carved_topography = filled_topography.Breaching_Lindsay2016();
      filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
    }
    else
    {
      filled_topography = filled_topography.fill(this_float_map["min_slope_for_fill"]);
    }

    string fill_name = "_filtered";
    filled_topography.write_raster((OUT_DIR+OUT_ID+fill_name), DEM_extension);
  }
  cout << "Done with filtering routine" << endl;
  cout << "==============================================" << endl << endl;

  //=================================================================================================
  //
  // .######..##.......####...##...##..........#####....####...##..##..######..######..##..##...####..
  // .##......##......##..##..##...##..........##..##..##..##..##..##....##......##....###.##..##.....
  // .####....##......##..##..##.#.##..........#####...##..##..##..##....##......##....##.###..##.###.
  // .##......##......##..##..#######..........##..##..##..##..##..##....##......##....##..##..##..##.
  // .##......######...####....##.##...........##..##...####....####.....##....######..##..##...####..
  // 
  //=================================================================================================
  cout << "\t I now need to do some flow routing." << endl;
  cout << "\t This is memory intensive. If it crashes it probably means" << endl;
  cout << "\t You either need to give your container more memory or you need" << endl;
  cout << "\t to downscale your DEM." << endl;
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions, filled_topography);
  // calcualte the distance from outlet
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  // calculate the flow accumulation
  cout << "\t Calculating flow accumulation (in pixels)..." << endl;
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  int threshold_contributing_pixels = this_int_map["Chan area threshold"];

  // load the sources
  cout << "\t Loading Sources, if you have them..." << endl;
  cout << "\t Source file is... " << CHeads_file << endl;
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

    cout << "The number of sources is: " << sources.size() << endl;
  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    cout << "The channel heads file *MUST* have a csv extension or this will crash!!" << endl;
    LSDRasterInfo ThisRI(filled_topography);
    string csv_filename = DATA_DIR+CHeads_file;
    LSDSpatialCSVReader CHeadCSV(ThisRI,csv_filename);
    sources = CHeadCSV.get_nodeindices_from_lat_long(FlowInfo);
    cout << "\t Got sources!" << endl;
    cout << "The number of sources is: " << sources.size() << endl;
  }

  // remove sources not in the filtered DEM
  vector<int> new_sources;
  int this_row, this_col;
  float NDV = filled_topography.get_NoDataValue();
  for (int i = 0; i < int(sources.size()); i++)
  {
    if (sources[i] != NDV) { 
      FlowInfo.retrieve_current_row_and_col(sources[i], this_row, this_col);
      if (this_row != NDV) {
        new_sources.push_back(sources[i]); 
      }
      else { cout << " No row and col at this source, removing" << endl;}
    }
    else { cout << "Not a valid source, removing" << endl; }
  }

  cout << "New number of sources is: " << new_sources.size() << endl;

  // now get the junction network
  LSDJunctionNetwork ChanNetwork(new_sources, FlowInfo);
  cout << "\t Got the channel network" << endl;

  // Print channels and junctions if you want them.
  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "I am going to print the channel network." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
    if ( this_bool_map["use_extended_channel_data"])
    {
      cout << "I am going to use the extended channel network data outputs." << endl;
      ChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, filled_topography);
    }
    else
    {
      ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
    }

    cout << "I've printed the channel network. " << endl;

    // convert to geojson if that is what the user wants
    // It is read more easily by GIS software but has bigger file size
    if ( this_bool_map["convert_csv_to_geojson"])
    {
      cout << "Let me convert that data to json." << endl;
      string gjson_name = OUT_DIR+OUT_ID+"_CN.geojson";
      LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_CN.csv");
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

  // needed for all channel routines
  LSDRaster DA_d8 = FlowInfo.write_DrainageArea_to_LSDRaster();

  // print junctions
  if( this_bool_map["print_junctions_to_csv"])
  {
    cout << "I am writing the junctions to csv." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_JN.csv";
    ChanNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_JN.geojson";
      LSDSpatialCSVReader thiscsv(channel_csv_name);
      thiscsv.print_data_to_geojson(gjson_name);
    }
  }

  //=========================================================================
  //
  // ..%%%%...%%%%%%..%%..%%...%%%%...%%......%%%%%%.
  // .%%........%%....%%%.%%..%%......%%......%%.....
  // ..%%%%.....%%....%%.%%%..%%.%%%..%%......%%%%...
  // .....%%....%%....%%..%%..%%..%%..%%......%%.....
  // ..%%%%...%%%%%%..%%..%%...%%%%...%%%%%%..%%%%%%.
  // ................................................
  // ..%%%%...%%..%%...%%%%...%%..%%..%%..%%..%%%%%%..%%.....
  // .%%..%%..%%..%%..%%..%%..%%%.%%..%%%.%%..%%......%%.....
  // .%%......%%%%%%..%%%%%%..%%.%%%..%%.%%%..%%%%....%%.....
  // .%%..%%..%%..%%..%%..%%..%%..%%..%%..%%..%%......%%.....
  // ..%%%%...%%..%%..%%..%%..%%..%%..%%..%%..%%%%%%..%%%%%%.
  //
  //=========================================================================
  if ( this_bool_map["extract_single_channel"])
  {
    string filename = OUT_DIR+OUT_ID+"_swath_channel_nodes.csv";
    ifstream test_single_chan(filename.c_str());
    if (test_single_chan)
    {
      cout << "I found a channel file already. I'll just load that" << endl;
    }
    else
    {
      cout << "I am going to extract a channel. " << endl;
      cout << "To do this I need to get the flow info object, which will take a little while." << endl;
      cout << "Also I am assuming your first point is the upstream point" << endl;
      cout << "and the second point is the downstream point." << endl;
      cout << "If you don't have a downstream point" << endl;
      cout << "I just use the first point as the source and follow it to" << endl;
      cout << "the edge of the DEM." << endl;

      cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
      LSDSpatialCSVReader CSVFile( RI, (DATA_DIR+this_string_map["channel_source_fname"]) );

      // get the x and y locations of the points
      vector<float> UTME;
      vector<float> UTMN;
      vector<double> latitudes = CSVFile.get_latitude();
      vector<double> longitudes = CSVFile.get_longitude();
      CSVFile.get_x_and_y_from_latlong(UTME,UTMN);

      int n_nodes = int(UTME.size());
      bool end_node_in_csv = false;
      if (n_nodes > 1)
      {
        end_node_in_csv = true;
        cout << "You have given me a channel source file with a second node." << endl;
        cout << "I am assuming this is an end node. " << endl;
        cout << "I will search for this node." << endl;
      }

      vector<int> NI_of_points =  CSVFile.get_nodeindices_from_lat_long(FlowInfo);

      // flow path in NI
      vector<int> final_channel_list;
      vector<int> flow_path = FlowInfo.get_flow_path(UTME[0],UTMN[0]);
      int downstream_NI;
      if(end_node_in_csv)
      {
        // snap the lower point to a channel
        cout << "The search radius is: " << this_int_map["search_radius_nodes"] << endl;
        downstream_NI = ChanNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(UTME[1],UTMN[1],
                            this_int_map["search_radius_nodes"],
                            this_int_map["threshold_stream_order_for_snapping"],
                            FlowInfo);
      }
      else
      {
        // There is no end node in the csv so we set the end node to a pixel that does not exist
        cout << "You didn't specify a downstream node, I'll follow down to the outlet" << endl;
        downstream_NI = -1;
      }

      // travel down the flow path until you either reach the outlet node or
      // get to the end
      for(int node = 0; node < int(flow_path.size()); node ++)
      {
        final_channel_list.push_back(flow_path[node]);
        //cout << "node: " << flow_path[node] << endl;

        // get the flow distance of the point, as well as the x and y locations
        int curr_row,curr_col;
        float curr_E,curr_N;
        FlowInfo.retrieve_current_row_and_col(flow_path[node],curr_row,curr_col);
        FlowInfo.get_x_and_y_from_current_node(flow_path[node],curr_E, curr_N);

        // If this is the downstream node, jump to the end of the flowpath
        if ( flow_path[node] == downstream_NI)
        {
          node = flow_path.size();
        }
      }

      cout << "writing to a csv file" << endl;

      string out_filename = OUT_ID+"_swath_channel";
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(final_channel_list,
                        OUT_DIR, out_filename,filled_topography, DistanceFromOutlet,
                        DA_d8);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string in_name = OUT_DIR+OUT_ID+"_swath_channel_nodes.csv";
        cout << "Printing the geojson of the swath channel." << endl;
        string gjson_name = OUT_DIR+OUT_ID+"_swath_channel.geojson";
        LSDSpatialCSVReader thiscsv(in_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
  }

//=========================================================================
//                 .__   __  .__       .__                  
//   _____  __ __|  |_/  |_|__|_____ |  |   ____          
//  /     \|  |  \  |\   __\  \____ \|  | _/ __ \         
// |  Y Y  \  |  /  |_|  | |  |  |_> >  |_\  ___/         
// |__|_|  /____/|____/__| |__|   __/|____/\___  >        
//       \/                   |__|             \/         
//        .__                                .__          
//   ____ |  |__ _____    ____   ____   ____ |  |   ______
// _/ ___\|  |  \\__  \  /    \ /    \_/ __ \|  |  /  ___/
// \  \___|   Y  \/ __ \|   |  \   |  \  ___/|  |__\___ \ 
//  \___  >___|  (____  /___|  /___|  /\___  >____/____  >
//      \/     \/     \/     \/     \/     \/          \/ 
//=========================================================================
  if(this_bool_map["extract_multiple_channels"])
  {
    cout << "I am going to extract multiple channels. " << endl;
    cout << "To do this I need to get the flow info object, which will take a little while." << endl;
    cout << "I am assuming that the CSV file has a list of upstream source points: I will " << endl;
    cout << "extract the channels downstream of each point." << endl; 

    cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
    LSDSpatialCSVReader CSVFile( RI, (DATA_DIR+this_string_map["channel_source_fname"]) );

    // get the x and y locations of the points
    vector<float> UTME;
    vector<float> UTMN;
    vector<double> latitudes = CSVFile.get_latitude();
    vector<double> longitudes = CSVFile.get_longitude();
    CSVFile.get_x_and_y_from_latlong(UTME,UTMN);
    cout << int(latitudes.size()) << " " << int(longitudes.size()) << endl;

    int n_nodes = int(UTME.size());
    vector<int> NI_of_points =  CSVFile.get_nodeindices_from_lat_long(FlowInfo);
    int downstream_NI = -1;

    cout << "The number of channels I will extract is: " << n_nodes << endl;

    // loop through the upstream nodes and get the channels
    for(int i = 0; i < n_nodes; i++)
    {
      cout << "I am getting the channel for node: " << NI_of_points[i] << endl;
      vector<int> final_channel_list;
      vector<int> flow_path = FlowInfo.get_flow_path(UTME[i],UTMN[i]);

      // travel down the flow path until you either reach the outlet node or
      // get to the end
      for(int node = 0; node < int(flow_path.size()); node ++)
      {
        final_channel_list.push_back(flow_path[node]);
        //cout << "node: " << flow_path[node] << endl;

        // get the flow distance of the point, as well as the x and y locations
        int curr_row,curr_col;
        float curr_E,curr_N;
        FlowInfo.retrieve_current_row_and_col(flow_path[node],curr_row,curr_col);
        FlowInfo.get_x_and_y_from_current_node(flow_path[node],curr_E, curr_N);

        // If this is the downstream node, jump to the end of the flowpath
        if ( flow_path[node] == downstream_NI)
        {
          node = flow_path.size();
        }
      } 
      cout << "writing to a csv file" << endl;

      string out_filename = OUT_ID+"_swath_channel_"+itoa(i);
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(final_channel_list,
                        OUT_DIR, out_filename,filled_topography, DistanceFromOutlet,
                        DA_d8);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string in_name = OUT_DIR+OUT_ID+"_swath_channel_"+itoa(i)+"_nodes.csv";
        cout << "Printing the geojson of the swath channel." << endl;
        string gjson_name = OUT_DIR+OUT_ID+"_swath_channel_"+itoa(i)+".geojson";
        LSDSpatialCSVReader thiscsv(in_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
  }

  //=========================================================================
  //
  // o      'O         o  o                      Oo                   o
  // O       o        O  O                      o  O                 O              o
  // o       O        o  o                     O    o                o
  // o       o        O  O                    oOooOoOo               O
  // O      O' .oOoO' o  o  .oOo. O   o       o      O 'OoOo. .oOoO' o  O   o .oOo  O  .oOo
  // `o    o   O   o  O  O  OooO' o   O       O      o  o   O O   o  O  o   O `Ooo. o  `Ooo.
  //  `o  O    o   O  o  o  O     O   o       o      O  O   o o   O  o  O   o     O O      O
  //   `o'     `OoO'o Oo Oo `OoO' `OoOO       O.     O  o   O `OoO'o Oo `OoOO `OoO' o' `OoO'
  //                                  o                                     o
  //                               OoO'                                  OoO'
  //
  //=========================================================================

  cout << endl << endl << "=========================================================" << endl;
  cout << "In this analysis (lsdtt-valley-metrics) I automatically extract floodplains" << endl;
  cout << "I am going to extract several rasters for you." << endl;
  cout << "A relief raster, a slope raster, and a valley raster" << endl;
  cout << "This will take a bit of computation. " << endl;
  cout << "WARNING if you are using 30m data and you are not in quite a big valley this will often" << endl;
  cout << "result in channel with a single pixel." << endl;
  //calculate the channel relief
  LSDRaster ChannelRelief;
  string relief_prefix = OUT_DIR+OUT_ID+"_channel_relief";
  string relief_test_file = OUT_DIR+OUT_ID+"_channel_relief.hdr"; 
  ifstream test_relief(relief_test_file.c_str());
  if (test_relief)
  {
    cout << "I found a channel relief raster already. I'll just load that" << endl;
    LSDRaster load_DEM(relief_prefix, DEM_extension);
    ChannelRelief = load_DEM;
  }
  else
  {
    cout << "\t Getting relief relative to channel" << endl;
    cout << "\t Threshold stream order = " << this_int_map["threshold_SO"] << endl;
    ChannelRelief = ChanNetwork.calculate_relief_from_channel(filled_topography, FlowInfo, this_int_map["threshold_SO"]);
    string relief_name = "_channel_relief";
    ChannelRelief.write_raster((OUT_DIR+OUT_ID+relief_name), DEM_extension);
    cout << "\t Got the relief!" << endl;
  }

  //get the slope
  LSDRaster Slope;
  string slope_prefix = OUT_DIR+OUT_ID+"_slope";
  string slope_test_file = OUT_DIR+OUT_ID+"_slope.hdr"; 
  ifstream test_slope(slope_test_file.c_str());
  if (test_slope)
  {
    cout << "I found a slope raster already. I'll just load that" << endl;
    LSDRaster load_DEM(slope_prefix, DEM_extension);
    Slope = load_DEM;
  }
  else
  {
    cout << "\t Calculating slope..." << endl;
    vector<LSDRaster> surface_fitting;
    vector<int> raster_selection(8, 0);
    raster_selection[1] = 1;             // this means you want the slope
    surface_fitting = filled_topography.calculate_polyfit_surface_metrics(this_float_map["surface_fitting_radius"], raster_selection);
    Slope = surface_fitting[1];
    cout << "\t Done!" << endl;
    string slope_name = "_slope";
    Slope.write_raster((OUT_DIR+OUT_ID+slope_name), DEM_extension);
  }

  float relief_threshold, slope_threshold;
  if(this_bool_map["use_absolute_thresholds"])
  {
    relief_threshold = this_int_map["relief_threshold"];
    slope_threshold = this_float_map["slope_threshold"];
  }
  else
  {
    // get the channel relief and slope_threshold using quantile-quantile plots
    cout << "Getting channel relief_threshold from QQ plots" << endl;
    string qq_fname = DATA_DIR+DEM_ID+"_qq_relief.txt";
    relief_threshold = ChannelRelief.get_threshold_for_floodplain_QQ(qq_fname, this_float_map["QQ_threshold"], this_int_map["relief_lower_percentile"], this_int_map["relief_upper_percentile"]);

    cout << "Getting slope_threshold from QQ plots" << endl;
    string qq_slope = path_name+DEM_ID+"_qq_slope.txt";
    slope_threshold = Slope.get_threshold_for_floodplain_QQ(qq_slope, this_float_map["QQ_threshold"], this_int_map["slope_lower_percentile"], this_int_map["slope_upper_percentile"]);
  }

  cout << "relief_threshold: " << relief_threshold << " slope_threshold: " << slope_threshold << endl;

  // get the floodplain object
  cout << "Getting the floodplain object" << endl;
  LSDFloodplain Floodplain(ChannelRelief, Slope, ChanNetwork, FlowInfo, relief_threshold, slope_threshold, this_int_map["minimum_patch_size"], this_int_map["threshold_SO"]);

  //print connected components
  //LSDIndexRaster CC = Floodplain.print_ConnectedComponents_to_Raster();
  //string cc_ext = "_CC";
  //CC.write_raster((OUT_DIR+OUT_ID+cc_ext), DEM_extension);

  // print a binary raster: valley = 1, non-valley = 0
  LSDIndexRaster BinaryRaster = Floodplain.print_BinaryRaster();
  string bin_ext = "_valley";
  BinaryRaster.write_raster((OUT_DIR+OUT_ID+bin_ext), DEM_extension);

  // fill holes in floodplain raster
  if (this_bool_map["fill_floodplain"])
  {
    cout << "I'm going to fill holes in the floodplain for you..." << endl;
    for (int i = 0; i < this_int_map["fill_loops"]; i++)
    {
      Floodplain.fill_holes_in_floodplain(this_int_map["fp_search_radius"]);
    }
    LSDIndexRaster FilledRaster = Floodplain.print_BinaryRaster();
    string bin_ext = "_valley_FILLED";
    FilledRaster.write_raster((OUT_DIR+OUT_ID+bin_ext), DEM_extension);
  }
  cout << "Finished with the floodplain extraction" << endl;
  cout << "=========================================================" << endl << endl;

  //=============================================================================================
  //                     _       _              _  _                                     _                        _       _                             
  //   __ __   __ _     | |     | |     ___    | || |    o O O   __      ___    _ _     | |_      _ _    ___     | |     (_)    _ _      ___      o O O 
  //   \ V /  / _` |    | |     | |    / -_)    \_, |   o       / _|    / -_)  | ' \    |  _|    | '_|  / -_)    | |     | |   | ' \    / -_)    o      
  //   _\_/_  \__,_|   _|_|_   _|_|_   \___|   _|__/   TS__[O]  \__|_   \___|  |_||_|   _\__|   _|_|_   \___|   _|_|_   _|_|_  |_||_|   \___|   TS__[O] 
  // _|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_| """"| {======|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""| {======| 
  // "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'./o--000'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'./o--000' 
  //
  //=============================================================================================
  if(this_bool_map["get_valley_centreline"])
  {
    cout << endl << endl << "=================================================" << endl;
	  cout << "I am looking for the valley centreline now." << endl;
	  string filename = OUT_DIR+OUT_ID+"_valley_centreline_nodes.csv";
    ifstream test_centreline(filename.c_str());
    if (test_centreline)
    {
      cout << "I found a centreline file already. I'll just load that" << endl;
    }
    else
    {

  	  cout << "I have not found a centreline file so I need to generate one." << endl;
  	  cout << "First what I will do is take the floodplain and punch a nodata raster" << endl;
      LSDRaster Punched = Floodplain.get_NoData_FloodplainRaster(filled_topography);
      string punched_name = OUT_DIR+OUT_ID+"_Punched";
      Punched.write_raster(punched_name, DEM_extension);  
  	  cout << "Finished with the punched raster." << endl;
        
      cout << "Grabbing rasters giving the distance and elevation of the nearest" << endl;
      cout << "pixels in the punched area to pixels with data." << endl;
  	  // This uses the topography raster to ignore any nodata values around the edges
      int n_pixels = 2;
      LSDRaster Buffered = Floodplain.Buffer_NoData_FloodplainRaster(Punched, n_pixels);
      vector<LSDRaster> nearest_rs = Buffered.get_nearest_distance_and_value_masks(filled_topography);

      string distance_R_name = OUT_DIR+OUT_ID+"_DistToND";
      string value_R_name = OUT_DIR+OUT_ID+"_ValueToND";    

      nearest_rs[0].write_raster(distance_R_name,DEM_extension);
      nearest_rs[1].write_raster(value_R_name,DEM_extension);

      // now get the lowest value in a window
      cout << "Now get the lowest value in a window." << endl;
      float window_radius = this_float_map["minimum_bank_elevation_window_radius"];
      int window_type = 1;
      bool find_maximum = false;

      cout << "Finding the lowest pixels along the river or valley edge." << endl;
      LSDRaster minimum_river_elev = nearest_rs[1].neighbourhood_statistics_local_min_max(window_radius, window_type, find_maximum);
      string minelev_R_name = OUT_DIR+OUT_ID+"_MinBankElev";
      minimum_river_elev.write_raster(minelev_R_name,DEM_extension);


      LSDRaster trough = nearest_rs[0];
      trough.raster_multiplier(this_float_map["trough_scaling_factor"]); 
      LSDRaster river_path = minimum_river_elev;
      LSDRaster filled_topography;

      LSDRaster add_trough = trough;
      
      // DO A FEW LOOPS
      int n_loops = this_int_map["centreline_loops"];
      for (int i = 0; i< n_loops; i++)
      {  
        cout << "Cutting and filling the river, loop number " << i << " of " << n_loops << endl;
        river_path = river_path.MapAlgebra_subtract(trough);
        LSDRaster carved_topography = river_path.Breaching_Lindsay2016();
        river_path = carved_topography.fill(this_float_map["min_slope_for_fill"]);
        
        // We keep track of the trough additions so at the end we add these back on except for the last one
        // this means that the depth of the middle of the trough should be set by the trough_scaling_factor
        // and not related to the number of loops
        // if (i < n_loops-2)
        // {
        //   add_trough = add_trough.MapAlgebra_add(trough);
        // }
      }

      //river_path = river_path.MapAlgebra_add(add_trough);
      //river_path.AdjustElevation(this_float_map["river_depth"]);


      // Now merge the rasters
      //cout << "Combining rasters." << endl;
      //LSDRaster topo_copy = topography_raster;
      //topo_copy.OverwriteRaster(river_path);

      string imposedriver_R_name = OUT_DIR+OUT_ID+"_RiverTrough";
      //topo_copy.write_raster(imposedriver_R_name,DEM_extension); 
      river_path.write_raster(imposedriver_R_name,DEM_extension);    

      if(this_bool_map["write_hillshade"])
      {
        cout << "Let me print the hillshade for you. " << endl;
        float hs_azimuth = 315;
        float hs_altitude = 45;
        float hs_z_factor = 1;
        LSDRaster hs_raster = river_path.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

        string hs_fname = OUT_DIR+OUT_ID+"_RiverTrough_hs";
        hs_raster.write_raster(hs_fname,DEM_extension);
      }

      // Now we need the flow path
      cout << "\t Flow routing for centreline. Note this is memory intensive. If your DEM is very large you may get a segmentation fault here..." << endl;
      // get a flow info object
      LSDFlowInfo FI(boundary_conditions,river_path);
      cout << "Finished flow routing." << endl; 

      // get the junction network
      int trough_threshold = 1;
      LSDIndexRaster FA = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      // now get the starting node for the flow path from your csv file
      // this is assuming you want a single channel
      int start_NI, downstream_NI;
      bool end_node_in_csv = false;
      vector<int> NI_of_points;
      if (this_bool_map["extract_single_channel"])
      {
        cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
        LSDSpatialCSVReader CSVFile( RI, (DATA_DIR+this_string_map["channel_source_fname"]) );

        // get the x and y locations of the points
        vector<float> UTME;
        vector<float> UTMN;
        vector<double> latitudes = CSVFile.get_latitude();
        vector<double> longitudes = CSVFile.get_longitude();
        CSVFile.get_x_and_y_from_latlong(UTME,UTMN);

        NI_of_points =  CSVFile.get_nodeindices_from_lat_long(FI);

        int n_nodes = NI_of_points.size();
        if (n_nodes > 1)
        {
          end_node_in_csv = true;
          cout << "You have given me a channel source file with a second node." << endl;
          cout << "I am assuming this is an end node. " << endl;
          cout << "I will search for this node." << endl;
        }
        // find the nearest channel to the starting node. This is necessary so that the start node is definitely in the centreline
        // this is a bit of a pain because of the different flow infos...
        int snapping_SO = 1;
        int snapped_NI = ChanNetwork.find_nearest_downslope_channel(UTME[0], UTMN[0], snapping_SO, FlowInfo);
        float snapped_UTME, snapped_UTMN;
        FlowInfo.get_x_and_y_from_current_node(snapped_NI, snapped_UTME, snapped_UTMN);
        start_NI = FI.get_node_index_of_coordinate_point(snapped_UTME, snapped_UTMN);
        cout << "start ni " << start_NI << endl;
        cout << snapped_UTME << " " << snapped_UTMN << endl;
      }
      else
      {
  		  cout << "I didn't find a channel source, so I am going to route from the highest elevation in the river path. " << endl;
        // get the maximum elevation
        int max_row,max_col;
        float max_elev = river_path.max_elevation(max_row,max_col);

        start_NI = FI.get_NodeIndex_from_row_col(max_row,max_col);
      }

      vector<int> flow_path = FI.get_flow_path(start_NI);
      cout << "n nodes " << int(flow_path.size()) << endl;
      // loop through the flow path and find the node nearest to the end node (if one exists)
      if(end_node_in_csv)
      {
        // snap the lower point to a channel
        cout << "The search radius is: " << this_int_map["search_radius_nodes"] << endl;
        downstream_NI = FI.find_nearest_node_in_vector(NI_of_points[1], flow_path);
      }
      else
      {
        // There is no end node in the csv so we set the end node to a pixel that does not exist
        cout << "You didn't specify a downstream node, I'll follow down to the outlet" << endl;
        downstream_NI = -1;
      }

      // travel down the flow path until you either reach the outlet node or
      // get to the end
      vector<int> final_channel_list;
      for(int node = 0; node < int(flow_path.size()); node ++)
      {
        final_channel_list.push_back(flow_path[node]);
        //cout << "node: " << flow_path[node] << endl;

        // get the flow distance of the point, as well as the x and y locations
        int curr_row,curr_col;
        float curr_E,curr_N;
        FI.retrieve_current_row_and_col(flow_path[node],curr_row,curr_col);
        FI.get_x_and_y_from_current_node(flow_path[node],curr_E, curr_N);

        // If this is the downstream node, jump to the end of the flowpath
        if ( flow_path[node] == downstream_NI)
        {
          node = flow_path.size();
        }
      }

      cout << "I need to calculate the flow distance now." << endl;
      LSDRaster FD = FI.distance_from_outlet();
      LSDRaster DA = FI.write_DrainageArea_to_LSDRaster();
      cout << "Got the dist from outlet and DA" << endl;

      string fname = OUT_ID+"_valley_centreline";
      FI.print_vector_of_nodeindices_to_csv_file_with_latlong(flow_path,OUT_DIR, fname,
                                                              river_path, FD, DA);      

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_valley_centreline.geojson";
        string full_csv = OUT_DIR+OUT_ID+"_valley_centreline_nodes.csv"; 
        LSDSpatialCSVReader thiscsv(full_csv);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
  }

  //======================================================================
  //               ____                       _     ____  __        
  //  _   ______ _/ / /__  __  __   _      __(_)___/ / /_/ /_  _____
  // | | / / __ `/ / / _ \/ / / /  | | /| / / / __  / __/ __ \/ ___/
  // | |/ / /_/ / / /  __/ /_/ /   | |/ |/ / / /_/ / /_/ / / (__  ) 
  // |___/\__,_/_/_/\___/\__, /    |__/|__/_/\__,_/\__/_/ /_/____/  
  //                    /____/                               
  //======================================================================
  if (this_bool_map["get_valley_widths"])
  {
    if (this_bool_map["extract_multiple_channels"])
    {
      cout << "You chose to extract multiple channels earlier in this routine." << endl;
      cout << "I am assuming you want to use those channel so I am overriding" << endl;
      cout << "Your valley_points_csv filename." << endl;

      LSDSpatialCSVReader CSVFile( RI, (DATA_DIR+this_string_map["channel_source_fname"]) );

      // get the x and y locations of the points
      vector<float> UTME;
      vector<float> UTMN;
      vector<double> latitudes = CSVFile.get_latitude();
      vector<double> longitudes = CSVFile.get_longitude();
      CSVFile.get_x_and_y_from_latlong(UTME,UTMN);

      int n_nodes = int(UTME.size());
      vector<int> NI_of_points =  CSVFile.get_nodeindices_from_lat_long(FlowInfo);
      for(int i = 0; i < n_nodes; i++)
      {
        cout << "This node index is: " << NI_of_points[i] << endl;
        LSDSpatialCSVReader ThisCSVFile(filled_topography, OUT_DIR+OUT_ID+"_swath_channel_"+itoa(i)+"_nodes.csv");
        cout << "Got the channel csv file" << endl;

        vector<int> nodes_to_remove;
        if(this_bool_map["remove_tributary_widths"])
        {
          cout << "I'm getting a list of tributary nodes to remove from the width calculation..." << endl;
          nodes_to_remove = Floodplain.remove_tributary_nodes(ThisCSVFile, ChanNetwork, FlowInfo, this_int_map["tributary_search_radius"], this_int_map["tributary_stream_order"]);
        }

        // now get the valley widths
        string bearing_fname =  OUT_DIR+OUT_ID+"_bearings_"+itoa(i)+".csv";
        string outfilename = OUT_DIR+OUT_ID+"_valley_widths_"+itoa(i);
        vector<vector<float>> valley_widths = Floodplain.calculate_valley_widths(ThisCSVFile, FlowInfo, DistanceFromOutlet, 
                            this_int_map["width_node_spacing"], this_bool_map["print_bearings"], bearing_fname, this_int_map["valley_banks_search_radius"], outfilename, nodes_to_remove);

      }
    }
    if (this_bool_map["extract_all_channels"])
    {
      cout << endl << endl << "=========================================================" << endl;
      cout << "I am going to calculate width for all channels in the DEM. This might be slow. " << endl;
      cout << "To do this I need to get the flow info object, which will take a little while." << endl;
      cout << "=========================================================" << endl;

      // basin junctions
      vector<int> basin_junctions = ChanNetwork.get_BaseLevelJunctions();
      cout << "The number of junctions I will analyse is: " << basin_junctions.size() << endl;

      // now get the valley widths
      string bearing_fname =  OUT_DIR+OUT_ID+"_bearings_network.csv";
      string outfilename = OUT_DIR+OUT_ID+"_valley_width_network";
      Floodplain.calculate_valley_widths_multiple_channels(basin_junctions, FlowInfo, ChanNetwork, DistanceFromOutlet, DA_d8, filled_topography,
                          this_int_map["width_node_spacing"], this_bool_map["print_bearings"], bearing_fname, this_int_map["valley_banks_search_radius"], outfilename);

    }
    else
    {
      // Read in the channel
      if (this_bool_map["extract_single_channel"])
      {
        cout << "You chose to extract a single channel earlier in this routine." << endl;
        cout << "I am assuming you want to use that channel so I am overriding" << endl;
        cout << "Your valley_points_csv filename." << endl;
        this_string_map["valley_points_csv"] = OUT_DIR+OUT_ID+"_swath_channel_nodes.csv";
      }
      if (this_bool_map["get_valley_centreline"])
      {
        cout << "You also chose to extract a valley centreline in this routine." << endl;
        cout << "I am assuming you want to use that centreline so I am overriding" << endl;
        cout << "to use that centreline instead of the river channel." << endl;
        this_string_map["valley_points_csv"] = OUT_DIR+OUT_ID+"_valley_centreline_nodes.csv"; 
      }
      LSDSpatialCSVReader CSVFile(filled_topography, this_string_map["valley_points_csv"]);

      cout << "Got the channel csv file" << endl;

      vector<int> nodes_to_remove;
      if(this_bool_map["remove_tributary_widths"])
      {
        cout << "I'm getting a list of tributary nodes to remove from the width calculation..." << endl;
        nodes_to_remove = Floodplain.remove_tributary_nodes(CSVFile, ChanNetwork, FlowInfo, this_int_map["tributary_search_radius"], this_int_map["tributary_stream_order"]);
      }

      // now get the valley widths
      string bearing_fname =  OUT_DIR+OUT_ID+"_bearings.csv";
      string outfilename = OUT_DIR+OUT_ID+"_valley_widths";
      vector<vector<float>> valley_widths = Floodplain.calculate_valley_widths(CSVFile, FlowInfo, DistanceFromOutlet, 
                          this_int_map["width_node_spacing"], this_bool_map["print_bearings"], bearing_fname, this_int_map["valley_banks_search_radius"], outfilename, nodes_to_remove);
    }
  }

  //=========================================================================
  //
  // oOoOOoOOo
  //     o
  //     o
  //     O
  //     o     .oOo. `OoOo. `OoOo. .oOoO' .oOo  .oOo. .oOo
  //     O     OooO'  o      o     O   o  O     OooO' `Ooo.
  //     O     O      O      O     o   O  o     O         O
  //     o'    `OoO'  o      o     `OoO'o `OoO' `OoO' `OoO'
  //
  //=========================================================================
  if (this_bool_map["get_terraces"])
  {
    // Read in the channel
    if (this_bool_map["extract_single_channel"])
    {
      cout << "You chose to extract a single channel earier in this routine." << endl;
      cout << "I am assuming you want to use that channel so I am overriding" << endl;
      cout << "Your valley_points_csv filename." << endl;
      this_string_map["valley_points_csv"] = OUT_DIR+OUT_ID+"_swath_channel_nodes.csv";
    }
    LSDSpatialCSVReader CSVFile(filled_topography, this_string_map["valley_points_csv"]);

    cout << "Got the channel csv file" << endl;

    // now run the swath
    vector<float> spaced_eastings, spaced_northings, spaced_distances;
    // get the easting and northings from the csv file
    CSVFile.get_x_and_y_from_latlong(spaced_eastings,spaced_northings);
    // get the flow distance along the channel
    string fd_column_name = "flowdistance(m)";
    if ( CSVFile.is_column_in_csv(fd_column_name) )
    {
      cout << "Super, I've got the distances in the file." << endl;
      spaced_distances = CSVFile.data_column_to_float(fd_column_name);
    }
    else
    {
      cout << "You are attempting to load a swath baseline file but it does not contain" << endl;
      cout << "A flowdistance(m) column (note that whitespace in column names are removed)" << endl;
      cout << "The column names in this file are: " << endl;
      CSVFile.print_data_map_keys_to_screen();
      cout << "I am afaid I will need to exit." << endl;
      exit(0);
    }

    string swath_data_prefix = OUT_DIR+OUT_ID;
    vector< LSDRaster > swath_rasters = filled_topography.find_nearest_point_from_list_of_points(spaced_eastings, spaced_northings, spaced_distances, this_float_map["swath_width"]);

    vector<int> channel_node_list = CSVFile.data_column_to_int("id");
    // now get the channel relief along the swath
    LSDRaster swath_channel_relief = ChannelRelief.isolate_to_smaller_raster(swath_rasters[0]);
    if (this_bool_map["print_swath_raster"])
    {
      string relief_fname = "_swath_channel_relief";
      swath_channel_relief.write_raster(OUT_DIR+OUT_ID+relief_fname, DEM_extension);
      }

    float mask_threshold = 1.0;
    bool below = 0;
    // remove any stupid slope values
    LSDRaster Slope_new = Slope.mask_to_nodata_using_threshold(mask_threshold, below);

    // get the terrace pixels
    LSDTerrace Terraces(swath_channel_relief, Slope_new, ChanNetwork, FlowInfo, relief_threshold, slope_threshold, this_int_map["minimum_patch_size"], this_int_map["threshold_SO"], this_int_map["minimum_terrace_height"]);
    LSDIndexRaster ConnectedComponents = Terraces.print_ConnectedComponents_to_Raster();
    string CC_ext = "_terraces";
    ConnectedComponents.write_raster((OUT_DIR+OUT_ID+CC_ext), DEM_extension);

    //cout << "\t Testing connected components" << endl;
    //vector <vector <float> > CC_vector = TestSwath.get_connected_components_along_swath(ConnectedComponents, RasterTemplate, this_int_map["NormaliseToBaseline"]);

    // print the terrace information to a csv
    //string csv_fname = "_terrace_info.csv";
    //string full_csv_name = DATA_DIR+DEM_ID+csv_fname;
    //cout << "The full csv filename is: " << full_csv_name << endl;
    //Terraces.print_TerraceInfo_to_csv(full_csv_name, filled_topography, swath_channel_relief, FlowInfo, TestSwath);
  }

  // Done, check how long it took
  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
