//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDTT_BasicMetrics.cpp
//
// This is a program for calculating basic landscape metrics.
// It includes options for slope, curvature, drainage area, and other metrics
//
// This program takes two arguments, the path name and the driver name
//
// The documentation is here:
// https://lsdtopotools.github.io/LSDTT_documentation/
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2022 Simon M. Mudd 2022
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
#include <utility>
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDChiTools.hpp"
#include "../LSDRasterMaker.hpp"
#include "../LSDBasin.hpp"

int main (int nNumberofArgs,char *argv[])
{

  string version_number = "0.8";
  string citation = "http://doi.org/10.5281/zenodo.4577879";

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the LSDTopoTools basic metrics tool!     ||" << endl;
  cout << "|| This program has a number of options for calculating||" << endl;
  cout << "|| simple landscape metrics.                           ||" << endl;
  cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| If you use these routines please cite:              ||" << endl;
  cout << "|| http://doi.org/10.5281/zenodo.4577879               ||" << endl;
  cout << "|| If you use the roughness routine please cite:       ||" << endl;
  cout << "|| https://www.doi.org/10.5194/esurf-3-483-2015        ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| Documentation can be found at:                      ||" << endl;
  cout << "|| https://lsdtopotools.github.io/LSDTT_documentation/ ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| This is LSDTopoTools2 version                       ||" << endl;
  cout << "|| " << version_number  << endl;
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
    ofs.open("./lsdtt-basic-metrics-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    cout << endl << endl << endl << "==============================================" << endl;    
    cout << "This is lsdtt-basic-metrics version number " << version_number << endl;
    cout << "If the version contains a 'd' then you are using a development version." << endl;
    cout << "=========================================================" << endl;
    ofstream ofs;
    ofs.open("./lsdtt-basic-metrics-version.txt");
    ofs << version_number << endl;
    ofs.close();

    exit(0);
  }


  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // this will contain the help file
  map< string, vector<string> > help_map;

  //==================================================================================
  //
  // .#####....####...#####....####...##...##..######..######..######..#####....####..
  // .##..##..##..##..##..##..##..##..###.###..##........##....##......##..##..##.....
  // .#####...######..#####...######..##.#.##..####......##....####....#####....####..
  // .##......##..##..##..##..##..##..##...##..##........##....##......##..##......##.
  // .##......##..##..##..##..##..##..##...##..######....##....######..##..##...####..
  //
  //=================================================================================
  
  //========================================================
  // Basic DEM preprocessing
  //========================================================
  float_default_map["minimum_elevation"] = 0.0;
  help_map["minimum_elevation"] = { "float","0.0","All elevation values below this become nodata if remove_seas is true.","Ususally 0."};

  float_default_map["maximum_elevation"] = 30000;
  help_map["maximum_elevation"] = {  "float","0.0","All elevation values above this become nodata if remove_seas is true.","Pick a big number."};

  float_default_map["min_slope_for_fill"] = 0.0001;
  help_map["min_slope_for_fill"] = {  "float","0.0001","Minimum slope between pixels for the filling algorithm.","Best not to change the default."};

  bool_default_map["raster_is_filled"] = false; // assume base raster is already filled
  help_map["raster_is_filled"] = {  "bool","false","This reads a pre-existing fill raster to save time.","You need to have printed the fill raster if you set this to true."};

  bool_default_map["carve_before_fill"] = false; // Implements a carving algorithm
  help_map["carve_before_fill"] = {  "bool","false","This implements a breaching algorithm before filling.","Good for landscapes with DEM obstructions (like roads) across the channels."};
 
  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  help_map["remove_seas"] = {  "bool","true","Slightly misleading name; it replaces both high and low DEM values with nodata.","This gets rid of low lying areas but also is handy when the nodata is not translated from the raw DEM and it is full of funny large numbers."};
 
  string_default_map["CHeads_file"] = "NULL";
  help_map["CHeads_file"] = {  "string","NULL","The name of a channel heads file.","You can output this csv file with the channel extraction algorithms. It contains latitude and longitude values of the channel heads."};
 
  bool_default_map["only_check_parameters"] = false;
  help_map["only_check_parameters"] = {  "bool","false","This just checks parameters without running an analysis.","For bug checking."};
 
  bool_default_map["replace_pixels"] = false;
  help_map["replace_pixels"] = {  "bool","false","When true the raster will replace specific pixels designated by a csv file.","Pixels to be replaced are designated in the pixels_to_replace_file"};
   
  string_default_map["pixels_to_replace_file"] = "pixels_to_replace.csv";
  help_map["pixels_to_replace_file"] = {  "string","pixels_to_replace.csv","The pixels to be replaced. Needs three columns X; Y; value. Header is ignored so it needs to be in this order","Incorrect formatting of file will crash the code."};

  float_default_map["elevation_change"] = 0;
  help_map["elevation_change"] = {  "float","0","Adjusts the elevation of the raster by this amount. Does nothing if the elevation change is 0. This prints to a new raster and exits.","For when you need to lift or drop the raster."};
  


  //========================================================
  // Various ways of trimming your raster
  //========================================================
  bool_default_map["remove_nodes_influenced_by_edge"] = false;
  help_map["remove_nodes_influenced_by_edge"] = {  "bool","false","Runs flow routing and then any node downstream of a node adjacent to a nodata pixel is turned to nodata.","Use this is you want to be sure all you data is in complete basins."};

  bool_default_map["remove_outer_nodes_for_edge_influence"] = false;
  help_map["remove_outer_nodes_for_edge_influence"] = {  "bool","false","Used in conjunction with remove_nodes_influenced_by_edge. In some cases the fill function causes an edge effect so this makes sure the edges with nodata can be propigated in the flow direction to make sure nodes influenced by the edge are removed","Use this is you want to be sure all you data is in complete basins."};

  bool_default_map["use_spiral_trimmer_for_edge_influence"] = false;
  help_map["use_spiral_trimmer_for_edge_influence"] = {  "bool","false","Used in conjunction with remove_nodes_influenced_by_edge. Makes sure raster is rectangular before running edge influence routine.","Use this is you want to be sure all you data is in complete basins."};




  bool_default_map["isolate_pixels_draining_to_fixed_channel"] = false;
  help_map["isolate_pixels_draining_to_fixed_channel"] = {  "bool","false","Runs flow routing and only take nodes that drain into a list of nodes defined in the fixed_channel_csv_name csv file.","Use this if you want to isolate a main drainage pathway."};

  //========================================================
  // raster trimming, to take care of rasters that have a bunch of nodata at the edges
  //========================================================
  bool_default_map["print_trimmed_raster"] = false;
  help_map["print_trimmed_raster"] = {  "bool","false","This is for rasters that have nodata around the edges. Trims away nodata edges and updates the raster georeferencing.","Use if your raster has lots of nodata and you want to shrink it to the area with data."};

  int_default_map["trimming_buffer_pixels"] = 0;
  help_map["trimming_buffer_pixels"] = {  "int","0","Used with the print_trimmed_raster Creates a buffer of this many nodata pixels around the edge.","You might want a buffer for use in numerical modelling or other applications."};

  //========================================================
  // Parameters for swath mapping
  //========================================================
  bool_default_map["calculate_swath_profile"] = false;
  help_map["calculate_swath_profile"] = {  "bool","false","Program sets this itself if you choose one of two options: calculate_swath_along_line or calculate_swath_along_channel. This makes a swath that prints a csv with distance along swath and the statistics.","The user should not select this option but instead use the swath_along_channel or swath_along_points flags instead. This flag is used only for completeness in the log file."};

  bool_default_map["calculate_swath_along_line"] = false;
  help_map["calculate_swath_along_line"] = {  "bool","false","Swath mapping where you give it a series of points in the csv file and it creates line segments between these points. This serves as the swath baseline.","Only one of this and calculate_swath_along_channel can be true."};

  bool_default_map["calculate_swath_along_channel"] = false;
  help_map["calculate_swath_along_channel"] = {  "bool","false","Swath mapping where you give it a starting point and a finish point (the latter is optional) and it will follow a flow path from the starting point and use that as the swath baseline.","Only one of this and calculate_swath_along_points can be true."};

  bool_default_map["print_swath_rasters"] = true;
  help_map["print_swath_rasters"] = {  "bool","true","If this is true and swath mapping is true then rasters showing the distance to the baseline and the closest node along the baseline file are printed.","These rasters can be used to visualise the pixels that go into the swath."};

  string_default_map["swath_points_csv"] = "swath.csv";
  help_map["swath_points_csv"] = {  "string","swath.csv","The name of the csv file for swath mapping. It can be the start and (optional) end point of the swath channel or a series of points that form a polyline that serves as the baseline for the swath.","Filename needs to include the csv extension and have latitude and longitude in the column headers."};

  float_default_map["swath_point_spacing"] = 500;
  help_map["swath_point_spacing"] = {  "float","500","Only used by calcualte_swath_along_line. How closely spaced the points are (in metres) along the baseline.","If this number is large (many times the DEM resolution) the swath will have an irregular edge so this should be close to the DEM resolution."};

  float_default_map["swath_bin_spacing"] = 1000;
  help_map["swath_bin_spacing"] = {  "float","1000","The spacing of the bins in metres that are used to calculate statistics along the swath.","Should be at least a few times bigger than swath_point_spacing or the pixel size of the DEM."};

  float_default_map["swath_width"] = 1000;
  help_map["swath_width"] = {  "float","1000","The width of the swath in metres.","Up to you. Don't let me tell you how wide your swath should be."};

  bool_default_map["make_categorised_swath"] = false;
  help_map["make_categorised_swath"] = {  "bool","false","Creates a swath that has binned values sorted and separated by a category. This is most often (or always) used with terrace data.","Categories are defined by a raster with prefix swath_category_raster_prefix. If you want terraces you should use the lsdtt-valley-metrics tool"};

  string_default_map["swath_category_raster_prefix"] = "terraces";
  help_map["swath_category_raster_prefix"] = {  "string","terraces","The prefix of the integer rasterthat will be used to categorise the swaths. This is usually the terrace file.","The terrace file is generated used lsdt-valley-metrics."};


  //========================================================
  // some tools for punching out polygons and then getting distance to nearest nodata
  // used for creating rivers in noisy DEMs
  //========================================================
  bool_default_map["punch_nodata"] = false;
  help_map["punch_nodata"] = {  "bool","false","This takes a raster mask and punches out nodata pixels where the mask is true based on a threshold.","The threshold is set by punch_threshold and the name of the raster to use as a mask is punch_raster_prefix."};

  bool_default_map["print_nearest_to_nodata_rasters"] = false; 
  help_map["print_nearest_to_nodata_rasters"] = {  "bool","false","This looks at interior nodata pixels and then makes rasters with values on the nodata pixels with distance to nearest nodata and the value of the nearest pixel with data.","The raster extensions are _DistToND and _ValueToND."};

  string_default_map["punch_raster_prefix"] = "DEM_punch";
  help_map["punch_raster_prefix"] = {  "string","DEM_punch","For raster punching this is the prefix of the mask.","This mask will turn values in the punched raster to nodata."};

  bool_default_map["belowthresholdisnodata"] = true;
  help_map["belowthresholdisnodata"] = {  "bool","true","For raster punching you set a threshold with punch_threshold and this controls whether values above or below will turn to nodata.","For example if the threshold is 1 and this is true then all mask values less than 1 will be turned to nodata in the punched raster."};

  float_default_map["punch_threshold"] = 1000;
  help_map["punch_threshold"] = {  "float","1000","For raster punching this is the threshold above or below which you set to nodata.","Use belowthresholdisnodata to determine if values above or below are set to nodata."};

  float_default_map["minimum_bank_elevation_window_radius"] = 20;
  help_map["minimum_bank_elevation_window_radius"] = {  "float","20","The valley routine will dig a trough where the starting elevation at the side is determined by the lowest bank elevation within this radius.","For best results this needs to be about the width of the valley."};
   
  float_default_map["river_depth"] = 1.0;
  help_map["river_depth"] = {  "float","1","If you are using create_valley_trough this sets the depth of the edges of the river below the minimum bank elevation.","Set to zero if you don't want the river inset into the banks."};

  float_default_map["trough_scaling_factor"] = 0.1;
  help_map["trough_scaling_factor"] = {  "float","0.1","The centreline routine digs a trough (to route flow through the centre of the valley) with a depth that is this factor times the distance from the edge. Increase this number for a deeper trough.","A deeper trough has more success at routing flow but you get backwater effects at the bottom of the valley. Play with this number to get the best result."};
  
  bool_default_map["create_valley_trough"] = false;
  help_map["create_valley_trough"] = {  "bool","false","If you have a valley mask this will punch a trough in that valley mask. The new river will have a depth of river_depth and will be sloping toward the middle of the valley based on trough_scaling_factor. Used to condition noisy rasters to have rivers through them.","We use this routine for flood modelling when the data over the river is very noisy."};

  bool_default_map["get_valley_centreline"] = false;
  help_map["get_valley_centreline"] = {  "bool","false","This extracts as single csv centreline for a valley mask.","Can be used to make cross sections."};

  int_default_map["centreline_loops"] = 5;
  help_map["centreline_loops"] = {  "int","5","The centreline routine iteratively digs and fills the trough and this sets the number of iterations.","We find 5 is about right but if the results are looking bad then try increasing this number. Warning: that might make things worse."};
 

  //========================================================
  // More river processing tools
  //========================================================
  bool_default_map["channel_and_valley_width_extraction"] = false;
  help_map["channel_and_valley_width_extraction"] = {  "bool","false","Some testing tools for channel and valley width. Operational routines in lsdtt-valley-metrics","For testing only. If you want valley width use lsdtt-valley-metrics."};
 
  string_default_map["channel_or_valley_raster_prefix"] = "channel";
  help_map["channel_or_valley_raster_prefix"] = {  "string","channel","Some testing tools for channel and valley width. Operational routines in lsdtt-valley-metrics","For testing only. If you want valley width use lsdtt-valley-metrics."};
 
  string_default_map["channel_or_valley_skeleton_prefix"] = "skeleton";
  help_map["channel_or_valley_skeleton_prefix"] = {  "string","skeleton","Some testing tools for channel and valley width. The prefix of the channel skeleton file.","For testing only. If you want valley width use lsdtt-valley-metrics."};
 
  int_default_map["test_scale"] = 4;
  help_map["test_scale"] = {  "int","4","Some testing tools for channel bearings. Sets the size of the bearing template.","For testing only. If you want valley width use lsdtt-valley-metrics."};
 
  float_default_map["test_bearing"] = 90;
  help_map["test_bearing"] = {  "float","90","Some testing tools for bearing","For testing only."};
 
  bool_default_map["test_bearing_template"] = false;
  help_map["test_bearing_template"] = {  "bool","false","A small debugging routine to test the size of the channel bearing template","Helps to bug check and visualise the direction of the channel where width is measured orthogonal to the bearing."};
   
  int_default_map["channel_bearing_node_spacing"]= 4;
  help_map["channel_bearing_node_spacing"] = {  "int","4","This sets how frequently the channel is sampled for its bearing.","Helps to bug check and visualise the direction of the channel where width is measured orthogonal to the bearing."};
  
  bool_default_map["print_channel_bearings"] = true; 
  help_map["print_channel_bearings"] = {  "bool","true","This prints to csv the directional bearings (looking downslope) of channel pixels.","Helps to bug check and visualise the direction of the channel where width is measured orthogonal to the bearing."};
 
  string_default_map["valley_points_csv"]  = "NULL";
  help_map["valley_points_csv"] = {  "string","NULL","If this isn't NULL it will load a valley point csv that need latitude and longitude and flow distance columns.","If this is NULL the valley centreline will be calculated algorithmically."};
   


  //========================================================
  // Calculate the basic relief within a window in a raster
  //========================================================
  bool_default_map["print_relief_raster"] = false;
  help_map["print_relief_raster"] = {  "bool","false","Prints relief in a moving window of size defined by relief_window.","A rudimentary measure of relief. If you want something more robust use the tools in lsdtt-hillslope-channel-coupling."};

  float_default_map["relief_window"] = 200;
  help_map["relief_window"] = {  "float","200","The radius of the relief window in metres.","Choose what you like."};

  int_default_map["relief_window_kernel_type"] = 1;
  help_map["relief_window_kernel_type"] = {  "int","1","A flag that defines the relief window.","1 = circular."};

  //========================================================
  // the most basic raster printing
  //========================================================
  bool_default_map["write_hillshade"] = false;
  help_map["write_hillshade"] = {  "bool","false","Write the hillshade raster.","You need this for a lot of our plotting routines. Filename includes _HS"};

  bool_default_map["print_raster_without_seas"] = false;
  help_map["print_raster_without_seas"] = {  "bool","false","Writes the raster without seas. Default is to overwrite the raster","DANGER this will replace your existing raster with the high and low points replaced by nodata. See the remove_seas flag"};

  bool_default_map["overwrite_raster_without_seas"] = true;
  help_map["overwrite_raster_without_seas"] = {  "bool","true","Overwrites the raster without seas.","DANGER if true this will replace your existing raster with the high and low points replaced by nodata. See the remove_seas flag"};

  bool_default_map["print_distance_from_outlet"] = false;
  help_map["print_distance_from_outlet"] = {  "bool","false","Prints a raster of the distance from the outlet.","Filename includes _FD"};

  bool_default_map["print_fill_raster"] = false;
  help_map["print_fill_raster"] = {  "bool","false","Prints the fill raster.","Filename includes _FILL"};

  // This converts all csv files to geojson (for easier loading in a GIS)
  bool_default_map["convert_csv_to_geojson"] = false;
  help_map["convert_csv_to_geojson"] = {  "bool","false","Converts csv files to geojson files","Makes csv output easier to read with a GIS. Warning: these files are much bigger than csv files."};

  //========================================================
  // Slope calculations
  //========================================================
  float_default_map["surface_fitting_radius"] = 30;
  help_map["surface_fitting_radius"] = {  "float","30","Our surface fitting routines fit a polynomial over the points with in a radius defined by surface_fitting_radius and then differentiate this surface to get the surface metrics like gradient and curvature","If not bigger than the pixel_size*sqrt(2) then will increase to that number."};

  bool_default_map["print_smoothed_elevation"]= false;
  help_map["print_smoothed_elevation"] = {  "bool","false","Prints smoothed surface after polynomial fitting.","Part of surface fitting metrics."};
 
  bool_default_map["print_slope"] = false;
  help_map["print_slope"] = {  "bool","false","Prints slope raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_aspect"]= false;
  help_map["print_aspect"] = {  "bool","false","Prints aspect raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_curvature"]= false;
  help_map["print_curvature"] = {  "bool","false","Prints curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_planform_curvature"]= false;
  help_map["print_planform_curvature"] = {  "bool","false","Prints planform curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_profile_curvature"]= false;
  help_map["print_profile_curvature"] = {  "bool","false","Prints profile curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_tangential_curvature"]= false;
  help_map["print_tangential_curvature"] = {  "bool","false","Prints tangential curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_point_classification"]= false;
  help_map["print_point_classification"] = {  "bool","false","Prints point classification raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_directional_gradients"] = false;
  help_map["print_directional_gradients"] = {  "bool","false","Prints two rasters: the x and y direction gradients after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["calculate_basin_statistics"] = false;
  help_map["calculate_basin_statistics"] = {  "bool","false","Prints some statistics for each basin. SMM needs to check what this does.","Part of surface fitting metrics."};

  //========================================================
  // Window size estimation
  //========================================================  
  bool_default_map["calculate_window_size"] = false;
  help_map["calculate_window_size"] = {  "bool","false","This is a routine that computes the optimal surface_fitting_radius for the surface fitting metrics.","Warning: time consuming."};


  //========================================================
  // Roughness calculations
  //========================================================
  bool_default_map["print_REI_raster"] = false;
  help_map["print_REI_raster"] = {  "bool","false","This computes the value of REI (rock exposure index) which is a roughness metric.","See DiBiase et al 2012 ESPL for details."};
 
  float_default_map["REI_critical_slope"] = 1.0;
  help_map["REI_critical_slope"] = {  "float","1.0","Critical slope value for the REI metric. Above this slope the ground is considered to be rock.","See DiBiase et al 2012 ESPL for details."};

  float_default_map["REI_window_radius"] = 10;
  help_map["REI_window_radius"] = {  "float","10","Radius of window withig which you test for pixels above the REI_critical_slope.","See DiBiase et al 2012 ESPL for details."};

  bool_default_map["print_roughness_rasters"] = false;
  help_map["print_roughness_rasters"] = {  "bool","false","Computes roughness rasters from divergence of surface normals.","See Milidowski et al 2015 ESURF for details."};
   
  float_default_map["roughness_radius"] = 3;
  help_map["roughness_radius"] = {  "float","3","Radius about centre point about which you calculate the surface normals.","See Milidowski et al 2015 ESURF for details. Should be around 3 times the pixel size."};
 
  //========================================================
  // drainage area
  //========================================================
  bool_default_map["print_dinf_drainage_area_raster"] = false;
  help_map["print_dinf_drainage_area_raster"] = {  "bool","false","Prints d-infinity drainage area raster.","Raster extension is Dinf."};
 
  bool_default_map["print_d8_drainage_area_raster"] = false;
  help_map["print_d8_drainage_area_raster"] = {  "bool","false","Prints d8 drainage area raster.","Raster extension is D8."};

  bool_default_map["print_QuinnMD_drainage_area_raster"] = false;
  help_map["print_QuinnMD_drainage_area_raster"] = {  "bool","false","Prints Quinn multidirection drainage area raster.","Raster extension is QMD."};

  bool_default_map["print_FreemanMD_drainage_area_raster"] = false;
  help_map["print_FreemanMD_drainage_area_raster"] = {  "bool","false","Prints Freeman multidirection drainage area raster.","Raster extension is FMD."};

  bool_default_map["print_MD_drainage_area_raster"] = false;
  help_map["print_MD_drainage_area_raster"] = {  "bool","false","Prints multidirection (fully divergent flow) drainage area raster.","Raster extension is MD."};


  //========================================================
  // Extracting a single channel
  //========================================================
  bool_default_map["extract_single_channel"] = false;
  help_map["extract_single_channel"] = {  "bool","false","This extracts a flow path from a line that starts at a point denoted in the file single_channel_source.","Used for imposing baselevel in various simulations."};

  bool_default_map["use_dinf_for_single_channel"] = false;
  help_map["use_dinf_for_single_channel"] = {  "bool","false","Uses dinf accumulation to print to the single channel file.","Might be useful in headwaters."};

  string_default_map["channel_source_fname"] = "single_channel_source";
  help_map["channel_source_fname"] = {  "string","single_channel_source","Name of the csv file without extension of the single channel source. Needs column headers latitude and longitude.","If latitude and longitude are not column headers in the csv this will not work."};

  // imposing a single channel
  bool_default_map["impose_single_channel"] = false;
  help_map["impose_single_channel"] = {  "bool","false","This forces the elevations of pixels read from the fixed_channel_csv_name.","Useful for enforcing a main stem flow path."};

  string_default_map["fixed_channel_csv_name"] = "NULL";
  help_map["fixed_channel_csv_name"] = {  "string","NULL","The prefix of the csv that holds the fixed channel information. csv much have latitude longitude flowdistance and elevation columns.","Obtain channel using extract_single_channel."};

  bool_default_map["buffer_single_channel"] = true;
  help_map["buffer_single_channel"] = {  "bool","true","Adds data to nodata pixels next to the single channel so that after flow routing the channel does not drain off the edge of the DEM.","Useful for when the single channel runs along the edge of the DEM."};

  bool_default_map["use_XY_for_buffer"] = false;
  help_map["use_XY_for_buffer"] = {  "bool","true","For buffering routine from buffer_single_channel this uses easting and northing instead of lat long","Preferred method is lat-long coordinates this is a legacy option"};

  bool_default_map["force_single_channel_slope"] = false;
  help_map["force_single_channel_slope"] = {  "bool","false","For the single channel this enforces a minimum slope between nodes. Need the flow distance in the single channel csv","Set the flow distance column with single_channel_fd_string"};

  bool_default_map["fixed_channel_dig_and_slope"] = false;
  help_map["fixed_channel_dig_and_slope"] = {  "bool","false","For the single channel this enforces the slope and digs the channel in a bit to ensure flow is routed through it","The digging is used to enforce flow routing"};

  float_default_map["single_channel_drop"] = 0;
  help_map["single_channel_drop"] = {  "float","0","For the single channel this drops the elevation","Use if you want a big drop in base level"};

  float_default_map["single_channel_dig"] = 0.01;
  help_map["single_channel_dig"] = {  "float","0.01","For the single channel this does the same as the drop but is used in the initial imposition of channel slope to try to enforce flow routing","Redundant with the fixed_channel_drop. Will probably remove at some point. "};

  string_default_map["single_channel_fd_string"] = "flow distance(m)";
  help_map["single_channel_fd_string"] = {  "string","flow distance(m)","The name of the flow distance column in the single channel csv.","Case sensitive"};

  string_default_map["single_channel_elev_string"] = "elevation(m)";
  help_map["single_channel_elev_string"] = {  "string","elevation(m)","The name of the elevation column in the single channel csv.","Case sensitive"};

  string_default_map["test_single_channel_name"] = "test_single_channel";
  help_map["test_single_channel_name"] = {  "string","test_single_channel","This is the filename of the single channel that prints after the slope and dig are imposed","File is used for bug checking."};


  // centre channel nodes and interpolate across gaps
  bool_default_map["centre_and_interpolate_channel_coordinates"] = false;
  help_map["centre_and_interpolate_channel_coordinates"] = {  "bool","false","This takes channel nodes that might not be exactly in the centre of pixels (possibly because they came from a different projection) and snaps them to pixel centres","Used for mapping channels from different projections onto a raster"};

  string_default_map["X_column_name"] = "X";
  help_map["X_column_name"] = {  "string","X","For channels that use projected data this is the X column name (sometimes called Easting)","Some csv files will have this of easting or some variant thereof"};

  string_default_map["Y_column_name"] = "Y";
  help_map["Y_column_name"] = {  "string","Y","For channels that use projected data this is the Y column name (sometimes called Northing)","Some csv files will have this of northing or some variant thereof"};

  // for reading in a channel csv file
  bool_default_map["use_xy_for_node_index"] = false;
  help_map["use_xy_for_node_index"] = {  "bool","true","For routine from centre_and_interpolate_channel_coordinates this uses easting and northing instead of lat long","Preferred method is lat-long coordinates this is a legacy option"};

  // this converts from files that are sent by nagra to a working single channel
  bool_default_map["prepare_single_channel_from_xy"] = false;
  help_map["prepare_single_channel_from_xy"] = {  "bool","false","This is a wrapper function that takes pixel x,y locations and attempts to create a single channel file with lat-long, area, flow distance, etc","Used to impose a channel from a different source onto a raster. This is designed for nagra work but might have future applications."};



  // Basic channel network
  int_default_map["threshold_contributing_pixels"] = 1000;
  help_map["threshold_contributing_pixels"] = {  "int","1000","The number of contributing pixels needed to start a channel using the threshold method.","This is in pixels not drainage area. More options are in the lsdtt-channel-extraction tool."};

  bool_default_map["print_stream_order_raster"] = false;
  help_map["print_stream_order_raster"] = {  "bool","false","Prints a raster with _SO in filename with stream orders of channel in the appropriate pixel.","Generates a big file so we suggest printing the network to csv."};

  bool_default_map["print_channels_to_csv"] = false;
  help_map["print_channels_to_csv"] = {  "bool","false","Prints the channel network to a csv file.","This version produces smaller files than the raster version."};

  bool_default_map["print_sources_to_csv"] = false;
  help_map["print_sources_to_csv"] = {  "bool","false","Prints the sources to a csv file.","Each source on its own row with latitude and longitude columns."};

  bool_default_map["use_extended_channel_data"] = false;
  help_map["use_extended_channel_data"] = {  "bool","false","If this is true you get more data columns in your channel network csv.","I will tell you what these columns are one day."};

  bool_default_map["print_channel_data_plus_surface_metrics"] = false;
  help_map["print_channel_data_plus_surface_metrics"] = { "bool", "false", "If this is true you get the channel network with various surface metric data like local slope relief and curvature.","This is slow because it needs to run the polyfitting routines."};

  bool_default_map["print_junction_index_raster"] = false;
  help_map["print_junction_index_raster"] = {  "bool","false","Prints a raster with junctions and their number.","Makes big files. It is better to use the csv version."};

  bool_default_map["print_junctions_to_csv"] = false;
  help_map["print_junctions_to_csv"] = {  "bool","false","Prints a csv with the locations and numbers of the junctions.","This is better to use than the raster version."};
 
  bool_default_map["print_channel_tips_raster"] = false;
  help_map["print_channel_tips_raster"] = {  "bool","false","This prints a raster with only the channel tips. It follows from the source down a number of pixels defined by channel_tips_pixel_distance","Extracts only the tips of the channel network. We use this when comparing topographically extracted channels with channels found using optical and radar satellite data."};
 
  int_default_map["channel_tips_pixel_distance"] = 10;
  help_map["channel_tips_pixel_distance"] = {  "int","10","The number of contributing pixels needed to start a channel using the threshold method.","This is in pixels not drainage area. More options are in the lsdtt-channel-extraction tool."};


  // Basin-based channel extraction
  bool_default_map["find_basins"] = false;
  help_map["find_basins"] = {  "bool","false","If true enters basin finding algorithms.","Used to try to extract basins of similar size. If you use outlets this flag is not required."};
 
  int_default_map["minimum_basin_size_pixels"] = 50000;
  help_map["minimum_basin_size_pixels"] = {  "int","50000","For basin finding algorithm this value the minimum size of a selected basin.","Will reject basins along edge."};
  
  int_default_map["maximum_basin_size_pixels"] = 1000000;
  help_map["maximum_basin_size_pixels"] = {  "int","1000000","For basin finding algorithm this value the maximum size of a selected basin.","Will reject basins along edge."};

  bool_default_map["fill_interior_nodata"] = false;
  help_map["fill_interior_nodata"] = {  "bool","false","For basin finding and other extration methods this removes interior nodata pixels.","Needs a window size that controls how big of an area is searched."};
  
  int_default_map["window_for_filling_interior_nodata"] = 10;
  help_map["window_for_filling_interior_nodata"] = {  "int","10","When you need to get rid of interior nodata this is the size of the window for filling.","Larger numbers fill more holes but slows down the computation"};




  bool_default_map["only_take_largest_basin"] = false;
  help_map["only_take_largest_basin"] = {  "bool","false","This only retains the largest complete basin in the raster.","Will reject basins along edge."};
  
  string_default_map["BaselevelJunctions_file"] = "NULL";
  help_map["BaselevelJunctions_file"] = {  "string","NULL","The name of a csv file with basin outlets for selecting basins using junction numbers.","An old method. You should use get_basins_from_outlets: true and basin_outlet_csv instead."};
  
  bool_default_map["get_basins_from_outlets"] = false;
  help_map["get_basins_from_outlets"] = {  "bool","false","Switches on the outlet based basin finding.","See BaselevelJunctions_file for format of outlets csv."};
  
  bool_default_map["prune_edge_draining_basins"] = true;
  help_map["prune_edge_draining_basins"] = {  "bool","true","If this is false then the basin finding algorithm will keep basins that drain from the edge meaning that there is a chance your chi coordinate and drainage area will be incorrect.","WARNING switch this off with extreme caution since you might get basins that drain from the edge of the DEM."};  

  int_default_map["search_radius_nodes"] = 8;
  help_map["search_radius_nodes"] = {  "int","8","A parameter for snapping to the nearest channel. It will search for the largest channel (by stream order) within the pixel window.","You will want smaller pixel numbers if you have a dense channel network."};
 
  int_default_map["threshold_stream_order_for_snapping"] = 2;
  help_map["threshold_stream_order_for_snapping"] = {  "int","2","If you are snapping to a channel the routine it will ignore channel with lower stream order than this number.","Set this to a higher number to avoid snapping to small channels."};
  
  string_default_map["basin_outlet_csv"] = "NULL";
  help_map["basin_outlet_csv"] = {  "string","NULL","A csv file with the lat long of basin outlets.","csv should have latitude and longitude columns and rows with basin outlets."};
  
  bool_default_map["extend_channel_to_node_before_receiver_junction"] = true;
  help_map["extend_channel_to_node_before_receiver_junction"] = {  "bool","true","For various basin extractions the basin snaps to the nearest junction. If this is true then the outlet of the basin is one pixel upstream of the receiver junction of the snapped channel.","If false it will pick the donor junction of the channel rather than one pixel above the receiver."};
    
  bool_default_map["print_basin_raster"] = false;
  help_map["print_basin_raster"] = {  "bool","false","This prints a raster where the values are the basin number.","You can combine this with python tools to get basin shapefiles."};

  bool_default_map["clip_raster_to_basins"] = false;
  help_map["clip_raster_to_basins"] = {  "bool","false","This the selected basins and then clips the raster to just those basins.","This will change the extent of the raster."};

  // finding major drainage divides
  bool_default_map["divide_finder"] = false;
  help_map["divide_finder"] = {  "bool","false","An experimental tool for finding the main drainage divide. This tags strips on the edge of the DEM and then follows the tagged pixels upstream.","You can choose horizontal or vertical strips."};
  
  bool_default_map["horizontal_strips"] = true;
  help_map["horizontal_strips"] = {  "bool","true","For use with the divide finder. If true strips are horizontal and vertical if not.","For now you need to hope your mountain is not diagonal across the DEM."};
 
  int_default_map["n_row_or_col_for_strips"] = 10;
  help_map["n_row_or_col_for_strips"] = {  "int","10","For use with the divide finder. Number of rows or columns to tag on the edges.","Make wider strips to ensure you capture all the valleys."};
 
  // Getting all the ridges
  bool_default_map["extract_ridges"] = false;
  help_map["extract_ridges"] = {  "bool","false","Extracts the ridges and prints the ridge information to csv.","This is a basic function. You should use lsdtt-hillslope-channel-coupling for full analysis."};
  

  // Tagging pixels
  bool_default_map["tag_nodes"] = false;
  help_map["tag_nodes"] = {  "bool","false","A routine for entering a raster with value and then propagating those values down the flow network.","The raster tag values are given by tagged_raster_input_name."};
    
  string_default_map["tagged_raster_input_name"] = "NULL";
  help_map["tagged_raster_input_name"] = {  "strig","NULL","The name of the raster with which to tag nodes.","The returned raster will have these values tagged in the output raster."};
  
  bool_default_map["tag_downslope_nodes"] = false;
  help_map["tag_downslope_nodes"] = {  "bool","false","Determines if the tagging routine will follow nodes upslope or downslope.","The default is to tag nodes upslope."};
  
  float_default_map["downslope_tagging_distance"] = 100;
  help_map["downslope_tagging_distance"] = {  "float","100","The distance through the flow network in metres cells will be tagged.","Make this number big if you want to follow through the entire drainage network."};
  
  float_default_map["upslope_tagging_distance"] = 100;
  help_map["upslope_tagging_distance"] = {  "float","100","The distance through the flow network in metres cells will be tagged.","Make this number big if you want to follow through the entire drainage network."};
 
  bool_default_map["tag_nodes_from_categorised_channel"] = false;
  help_map["tag_nodes_from_categorised_channel"] = {  "bool","false","This will read the channel file single_channel_csv and then take a column and categories it based on the boundary values in category_boundaries.","The data to be categorised is given by column_to_categorise. Most useful to get at drainage divides"};
    
  string_default_map["category_boundaries"] = "0";
  help_map["category_boundaries"] = { "string","0","This is a comma separated string into which you place category boundaries.","This will be automatically sorted to the order of the values does not matter. Used with tag_nodes_from_categorised_channel"};

  string_default_map["column_to_categorise"] = "flow_distance";
  help_map["column_to_categorise"] = {  "string","flow_distance","If you are categorising nodes in a channel this is the column that is categorised.","Used with tag_nodes_from_categorised_channel and the boundaries are set by category_boundaries"};

  string_default_map["random_seed_for_categorisation_test"] = 10;
  help_map["random_seed_for_categorisation_test"] = {  "int","10","Allows you to input a seed for the categorisation test only","This seed only works on the categorisation test and nowhere else"};


  // Some chi coordinate settings
  float_default_map["A_0"] = 1.0;
  help_map["A_0"] = {  "float","1.0","The A_0 parameter for chi computation. See https://doi.org/10.1002/esp.3302","Usually set to 1 so that the slope in chi-elevation space is the same as k_sn"};
   
  float_default_map["m_over_n"] = 0.5;
  help_map["m_over_n"] = {  "float","0.5","The concavity index for chi calculations. Usually denoted as the Greek symbol theta.","Default is 0.5 but possibly 0.45 is better as Kwang and Parker suggest 0.5 leads to unrealistic behaviour in landscape evolution models."};

  bool_default_map["print_chi_data_maps"] = false;
  help_map["print_chi_data_maps"] = {  "bool","false","If true prints the chi network to csv.","csv file has chidatamaps in the filename. Has the locations of the channel pixels with their chi coordinates and other information."};
  
  if (bool_default_map["print_chi_data_maps"] == true)
  {
    cout << "You want the chi data maps, so I am setting the basin finding to true." << endl;
    bool_default_map["find_basins"] = true;
  }

  // The wiener filter
  bool_default_map["print_wiener_filtered_raster"] = false;
  help_map["print_wiener_filtered_raster"] = {  "bool","false","Runs the raster through a wiener filter and prints the result.","Output file has _Wfilt in filename."};

  // This burns a raster value to any csv output of chi data
  // Useful for appending geology data to chi profiles
  bool_default_map["burn_raster_to_csv"] = false;
  help_map["burn_raster_to_csv"] = {  "bool","false","Takes a raster with burn_raster_prefix and then samples that raster with the points in the csv file. The new column will be burn_data_csv_column_header.","Useful for adding raster data to csv file. Often used to add lithological information to csv data (you must rasterize the lithology data first."};

  string_default_map["burn_raster_prefix"] = "NULL";
  help_map["burn_raster_prefix"] = {  "string","NULL","The prefix of the raster to burn to a csv.","No extension required."};

  string_default_map["burn_data_csv_column_header"] = "burned_data";
  help_map["burn_data_csv_column_header"] = {  "string","burned_data","Column header in csv of data burned from raster.","For example lithocode."};

  string_default_map["csv_to_burn_name"] = "NULL";
  help_map["csv_to_burn_name"] = {  "string","NULL","Name of csv file to which data will be burned.","You need to include the csv extension."};

  // This is for junction angles
  bool_default_map["print_junction_angles_to_csv"] = false;
  help_map["print_junction_angles_to_csv"] = {  "bool","false","Prints a csv with ALL the locations of the junctions and associated statistics.","csv file contains junction angles; bending angles; areas; gradients; and other information."};

  bool_default_map["print_junction_angles_to_csv_in_basins"] = false;
  help_map["print_junction_angles_to_csv_in_basins"] = {  "bool","false","Prints a csv with the locations of the junctions within selected basins and associated statistics.","csv file contains junction angles; bending angles; areas; gradients; and other information."};

  float_default_map["SA_vertical_interval"] = 10;
  help_map["SA_vertical_interval"] = {  "float","10","This is used in both slope-area routines and also junction angle routines. It sets the vertical drop over which the gradient is measured and the junction angle is measured on points between a tributary within the height and the junction.","For S-A analysis the value should be greater than the vertical uncertainty of your DEM. For junction angles if you set this to a large number it will measure angles between the junction the donor junctions and the receiver junction."};


  // Now for connectivity
  //bool_default_map["calculate_connectivity_index"] = false;
  //help_map["calculate_connectivity_index"] = {  "bool","false","Calculates the connectivity index NOT WORKING YET.","Based on https://doi.org/10.1016/j.cageo.2017.10.009"};

  // calculate the hypsometric integral for the channel network
  bool_default_map["calculate_hypsometric_integral"] = false;
  help_map["calculate_hypsometric_integral"] = {  "bool","false","Calculates the hypsometric integral for every part of the channel network.", "Useful for examining the distribution of usptream elevations."};

  float_default_map["HI_bin_width"] = 500;
  help_map["HI_bin_width"] = {  "float","500","The elevation bins used to calculate the hypsometric integral.","Should be set by examining the elevation range in your DEM to get a reasonable bin width."};
 
  float_default_map["HI_lower_limit"] = 0;
  help_map["HI_lower_limit"] = {  "float","0","The lowest elevation in the hypsometric integral.","Should be set by examining the elevation range in your DEM to get a reasonable bin width."};
 
  float_default_map["HI_upper_limit"] = 8000;
  help_map["HI_upper_limit"] = {  "float","0","The highest elevation in the hypsometric integral.","Should be set by examining the elevation range in your DEM to get a reasonable bin width."};
 

  // Steady state parameters
  bool_default_map["create_a_steady_landscape"] = false;
  help_map["create_a_steady_landscape"] = {  "bool","false","Turn on to create a steady state landscape initiated by the diamond square algorithm.","Used to quickly spin up steady landscapes."};

  float_default_map["steady_state_resolution"] = 30.0;
  help_map["steady_state_resolution"] = {  "float","30.0","Grid spacing in a synthetic steady state landscape.","Used to quickly spin up steady landscapes."};

  int_default_map["steady_state_nrows"] = 200;
  help_map["steady_state_nrows"] = {  "int","200","Number of rows in a synthetic steady state landscape.","Used to quickly spin up steady landscapes."};

  int_default_map["steady_state_ncols"] = 400;
  help_map["steady_state_ncols"] = {  "int","400","Number of column in a synthetic steady state landscape.","Used to quickly spin up steady landscapes."};

  int_default_map["steady_state_n_cycles"] = 20;
  help_map["steady_state_n_cycles"] = {  "int","20","The number of cycles to use in the steady state landscape.","Increased cycles ensure landscape is fully dissected."};

  int_default_map["steady_state_ds_feature_order"] = 3;
  help_map["steady_state_ds_feature_order"] = {  "int","3","For the steady state landscape this is the diamond square feature order which contols how repetituve the landscape is.","Values between 3 and 6 are usually fine."};

  float_default_map["random_seed"] = 1.0;
  help_map["random_seed"] = {  "float","1.0","The seed you use for random values","For routines that have randomness this sets the seed. Values will be random but you will get the same numbers if you use the same seed. Only implemented for the steady state landscape at the moment."};

  int_default_map["snap_print_rate"] = 10;     // how frequently the snapping is printed in terms of snap cycles
  help_map["snap_print_rate"] = {  "int","10","This prints the DEM evey n cycles set by this parameters.","Used to just get the last few snaps for an animation."};

  float_default_map["steady_state_n"] = 1.0;
  help_map["steady_state_n"] = {  "float","1.0","The n value for the steady state simulation.","Used only in steady state computations."};

  float_default_map["steady_state_m"] = 0.45;
  help_map["steady_state_m"] = {  "float","0.45","The m value for the steady state simulation.","Used only in steady state computations."};

  float_default_map["steady_state_threshold_hillslope_area"] = 100000.0;
  help_map["steady_state_threshold_hillslope_area"] = {  "float","100000.0","Threshold drainage area for a hillslope pixel in m^2 for steady state simulations.","Used only in steady state computations."};  

  float_default_map["steady_state_threshold_hillslope_gradient"] = 0.4;
  help_map["steady_state_threshold_hillslope_gradient"] = {  "float","0.4","Threshold gradient hillslope pixel in m^2 for steady state simulations.","Used only in steady state computations."};  

  float_default_map["steady_state_U"] = 0.0004;
  help_map["steady_state_U"] = {  "float","0.0004","Uplift rate in m/yr for steady state simulations.","Used only in steady state computations."};  

  float_default_map["steady_state_target_fluvial_relief"] = 200.0;
  help_map["steady_state_target_fluvial_relief"] = {  "float","200.0","Target relief in metres for steady state simulations.","Used only in steady state computations."};  


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
  //============================================================================
  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./lsdtt-basic-metrics-README.csv" << endl;
    string help_prefix = "lsdtt-basic-metrics-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }


  // Now print the parameters for bug checking
  cout << "PRINT THE PARAMETERS..." << endl;
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();

  cout << "Read filename is: " <<  DATA_DIR+DEM_ID << endl;
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;

  // A little switch to find basins if you supply an outlet file
  if( this_bool_map["get_basins_from_outlets"] )
  {
    this_bool_map["find_basins"] = true;
  }

  // Some switches to turn on the connectivity index
  if ( this_bool_map["calculate_connectivity_index"])
  {
    this_bool_map["print_dinf_drainage_area_raster"] = true;    
  }

  //============================================================================
  //
  //..####...######..######...####...#####...##..##...........####...######...####...######..######.........
  //.##........##....##......##..##..##..##...####...........##........##....##..##....##....##.............
  //..####.....##....####....######..##..##....##.............####.....##....######....##....####...........
  //.....##....##....##......##..##..##..##....##................##....##....##..##....##....##.............
  //..####.....##....######..##..##..#####.....##.............####.....##....##..##....##....######.........
  //
  //============================================================================
  if (this_bool_map["create_a_steady_landscape"])
  {
    cout << endl << endl << endl << endl << "<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]" << endl;
    cout << "Hello, I am going to create a steady state landscape for you." << endl;
    cout << "I will use this steady landscape for all further analyses you requested." << endl;
    cout << "<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]" << endl;    

    // Initiate the DEM
    LSDRasterMaker RM;
    float placeholder_elevation = 1.0;
    RM.resize_and_reset(this_int_map["steady_state_nrows"],this_int_map["steady_state_ncols"],
                         this_float_map["steady_state_resolution"],placeholder_elevation);


    long this_seed = long(this_float_map["random_seed"]);
    float noise_relief = 0.25; // set the noise to something small;
    float ds_relief = 10;
    float parabola_relief = 10;

    // Now run the diamond square 
    RM.create_diamond_square_surface(this_int_map["steady_state_ds_feature_order"], ds_relief, 
                    noise_relief, parabola_relief, this_float_map["min_slope_for_fill"], &this_seed);

    // Add a parabolic surface (this forces the basins to drain to the )

    // Print the DS surface with hillshade
    
    string DS_initial_fname = DATA_DIR+DEM_ID+"_ds_initial";
    LSDRaster DS_initial = RM.return_as_raster();
    DS_initial.write_raster(DS_initial_fname,"bil");

    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = DS_initial.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = DATA_DIR+DEM_ID+"_ds_initial_hs";
    hs_raster.write_raster(hs_fname,raster_ext);

    // Tune the K value for the desired relief
    float target_K = RM.tune_K_for_relief(this_float_map["steady_state_U"], 
                         this_float_map["steady_state_target_fluvial_relief"],
                         this_float_map["steady_state_m"],
                         this_float_map["steady_state_n"],
                         this_bool_map["carve_before_fill"], 
                         this_float_map["min_slope_for_fill"],
                         this_float_map["steady_state_threshold_hillslope_area"]);
    cout << "The value of K for your desired relief of " << this_float_map["steady_state_target_fluvial_relief"] << " is: " << target_K << endl;

    // now enter snapping snapping cycle
    for(int i = 1; i<= this_int_map["steady_state_n_cycles"]; i++)
    {
      cout << "\rSnap number " << i << " of "  << this_int_map["steady_state_n_cycles"] << flush;
      RM.snap_to_steady(target_K, this_float_map["steady_state_U"], 
                        this_float_map["steady_state_threshold_hillslope_gradient"],
                        this_float_map["steady_state_m"],
                        this_float_map["steady_state_n"],
                        this_bool_map["carve_before_fill"], 
                        this_float_map["min_slope_for_fill"],
                        this_float_map["steady_state_threshold_hillslope_area"]);
    }

    // Now print the final state
    // We use DATA_DIR and DEM_ID because if this switch is on it uses
    string DS_snapped_fname = OUT_DIR+OUT_ID+"_ds_snapped";
    LSDRaster DS_snapped = RM.return_as_raster();
    DS_snapped.write_raster(DS_snapped_fname,"bil");

    cout << "Updating the DEM prefix so we use the Diamond square raster as the base DEM" << endl;
    DEM_ID = OUT_ID+"_ds_snapped";
    OUT_ID = DEM_ID;
    DATA_DIR = OUT_DIR;
  }



  // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);

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

    if (this_bool_map["print_raster_without_seas"])
    {
      if (this_bool_map["overwrite_raster_without_seas"])
      {
        cout << "I'm replacing your raster with a raster without seas." << endl;
        string this_raster_name = OUT_DIR+OUT_ID;
        topography_raster.write_raster(this_raster_name,raster_ext);
      }
      else
      {
        cout << "I'm writing a remove seas raster with the suffix RSEAS" << endl;
        string this_raster_name = OUT_DIR+OUT_ID+"_RSEAS";
        topography_raster.write_raster(this_raster_name,raster_ext);
      }
    }
  }
  else
  {
    LSDRaster start_raster((DATA_DIR+DEM_ID), raster_ext);
    topography_raster = start_raster;
  }
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;

  // now for the pixel replacement. This exits upon completion
  if (this_bool_map["replace_pixels"])
  {
    string pixel_if_name = DATA_DIR+this_string_map["pixels_to_replace_file"];

    topography_raster.replace_pixels(pixel_if_name);

    string this_raster_name = OUT_DIR+OUT_ID+"_pixels_replaced";
    topography_raster.write_raster(this_raster_name,raster_ext);

    cout << "I replaced some pixels in your DEM and am now exiting" << endl;
    exit(0);
  }

  // Adjust the elevation of the raster
  if (this_float_map["elevation_change"] != 0)
  {
    topography_raster.AdjustElevation(this_float_map["elevation_change"]);
    cout << "Adjusted the elevation of your raster by " << this_float_map["elevation_change"] << endl;
    topography_raster.write_raster(DATA_DIR+DEM_ID+"_elevadjusted", "bil");
    cout << "I adjusted the elevations of your raster and printed the new raster." << endl;
    cout << "It has the extension _elevadjusted" << endl;
    cout << "I am now exiting." << endl;
    exit(0);
 }


  if(this_bool_map["only_check_parameters"])
  {
    cout << "You set the only_check_parameters flag to true; I have only printed" << endl;
    cout << "the parameters to file and am now exiting." << endl;
    exit(0);
  }

  //============================================================================
  // Raster trimming
  // This trims the raster to the smallest nodata.
  // It only prints the trimmed raster! It will not use it in computation.
  // This is because we want to avoid messed up georeferencing.
  // Print the trimmed raster and then use in raster the next steps
  // of raster processing
  if(this_bool_map["print_trimmed_raster"])
  {
    cout << "Let me trim that raster for you." << endl;
    LSDRaster trimmed_raster = topography_raster.RasterTrimmerPadded(this_int_map["trimming_buffer_pixels"]);
    string this_raster_name = OUT_DIR+OUT_ID+"_TRIM";
    trimmed_raster.write_raster(this_raster_name,raster_ext);
  }

  //=========================================================================
  // ..####...##...##...####...######..##..##.
  // .##......##...##..##..##....##....##..##.
  // ..####...##.#.##..######....##....######.
  // .....##..#######..##..##....##....##..##.
  // ..####....##.##...##..##....##....##..##.
  // 
  //=========================================================================
  // We need some logic to keep people from messing up.
  if (this_bool_map["calculate_swath_along_line"] && this_bool_map["calculate_swath_along_channel"])
  {
    cout << "You can only use one swath method." << endl;
    cout << "Choose *EITHER* calculate_swath_along_line *OR* calculate_swath_along_channel" << endl;
    cout << "Exiting so you can fix." << endl;
    exit(0);
  }
  if (this_bool_map["calculate_swath_along_line"] || this_bool_map["calculate_swath_along_channel"])
  {
    this_bool_map["calculate_swath_profile"] = true;
  }
  if (not this_bool_map["calculate_swath_along_line"] && 
      not this_bool_map["calculate_swath_along_channel"] &&
      this_bool_map["calculate_swath_profile"])
  {
    cout << "You can't set calculate_swath_profile to true without choosing a method." << endl;
    cout << "Please select EITHER calculate_swath_along_line OR calculate_swath_along_channel" << endl;
    cout << "You cannot choose both." << endl;
    cout << "Exiting for you to fix." << endl;
    exit(0);
  }


  // okay, now the calculations. 
  if(this_bool_map["calculate_swath_profile"])
  {
    cout << "This is exciting! I get to make a swath profile!!" << endl;
    LSDRasterInfo RI(topography_raster);


    LSDSpatialCSVReader CSVFile(RI,this_string_map["swath_points_csv"]);

    // get the x and y locations of the points
    vector<float> UTME;
    vector<float> UTMN;
    CSVFile.get_x_and_y_from_latlong(UTME,UTMN); 

    int n_nodes = int(UTME.size()); 
    CSVFile.check_if_points_are_in_raster();
    vector<bool> is_in_raster = CSVFile.get_if_points_are_in_raster_vector();
    for (int i= 0; i< n_nodes; i++)
    {
      //cout << "x: " << UTME[i] << ", y: " << UTMN[i] << endl;
      if(not is_in_raster[i])
      {
        cout << "There is a point in the swath baseline that is not in the raster." << endl;
        cout << "You need to either correct your raster or correct this point." << endl;
        exit(0);
      }
    }   

    vector<float> spaced_eastings;
    vector<float> spaced_northings; 
    vector<float> spaced_distances;

    if(this_bool_map["calculate_swath_along_line"])
    {
      cout << "I am going to get the swath profile along a line that is" << endl;
      cout << "made up of segments derived from your swath points file" << endl;

      // Now segment the line
      float spacing = this_float_map["swath_point_spacing"];
      cout << "Segmenting the swath baseline with a spacing of " << this_float_map["swath_point_spacing"] << endl;

      evenly_spaced_points_along_polyline(UTME, UTMN, spacing, spaced_eastings, spaced_northings, spaced_distances);
      vector< pair<float,float> > points = evenly_spaced_points_along_polyline(UTME, UTMN, spacing);
      cout << "Got the segmented baseline." << endl;  
      
      cout << "Number of points is: " << points.size() << endl;
      int n_points = int(spaced_northings.size());
      //for(int i = 0; i< n_points; i++)
      //{
      // cout << spaced_eastings[i] << "," << spaced_northings[i] <<endl;
      //}

      cout << "Converting pairs to UTM vectors." << endl;
      vector<float> swath_UTME, swath_UTMN;
      convert_pairs_to_x_and_y_vecs(points, swath_UTME, swath_UTMN);
      cout << "Done." << endl;

      ofstream points_out;
      string swath_baseline_fname = OUT_DIR+OUT_ID+"_swath_baseline.csv";
      points_out.open(swath_baseline_fname);
      points_out << "easting,northing" << endl;
      for (int i = 0; i < int(swath_UTME.size()); i++)
      {
        points_out << spaced_eastings[i]<< "," << spaced_northings[i] << endl;
      }
      points_out.close();

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_swath_baseline.geojson";
        LSDSpatialCSVReader thiscsv(swath_baseline_fname);
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    if (this_bool_map["calculate_swath_along_channel"])
    {
      cout << "I am going to calculate a swath along a channel!" << endl;
      cout << "I need to check your channel file to see if it is a full channel record" << endl;
      cout << "Or just the start and end points (i.e., if it has two nodes)" << endl;
      if (n_nodes == 2)
      {
        cout << "Okay, the channel file has two nodes. " << endl;
        cout << "This means I need to find the channel. " << endl;
        cout << "To do this I need to get the flow info object, which will take a little while." << endl;
        cout << "Also I am assuming your first point is the upstream point" << endl;
        cout << "and the second point is the downstream point." << endl;

        LSDRaster filled_topography,carved_topography;
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
          if(this_bool_map["carve_before_fill"])
          {
            carved_topography = topography_raster.Breaching_Lindsay2016();
            filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
          }
          else
          {
            filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
          }
        }
        LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);
        cout << "Finished flow routing for the swath." << endl;

        LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
        vector<int> sources;
        sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);
        LSDRaster FlowDistance = FlowInfo.distance_from_outlet();
        LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

        // Now get the upstream and downstream node
        CSVFile.get_x_and_y_from_latlong(UTME,UTMN); 
        vector<double> latitudes = CSVFile.get_latitude();
        vector<double> longitudes = CSVFile.get_longitude();
        vector<int> NI_of_points =  CSVFile.get_nodeindices_from_lat_long(FlowInfo);

        // flow path in NI
        vector<int> flow_path = FlowInfo.get_flow_path(UTME[0],UTMN[0]);

        // snap the lower point to a channel
        int downstream_NI;
        cout << "The search radius is: " << this_int_map["search_radius_nodes"] << endl;
        downstream_NI = JunctionNetwork.get_nodeindex_of_nearest_channel_for_specified_coordinates(UTME[1],UTMN[1],
                            this_int_map["search_radius_nodes"],  
                            this_int_map["threshold_stream_order_for_snapping"], 
                            FlowInfo);

        // travel down the flow path until you either reach the outlet node or 
        // get to the end
        vector<int> final_channel_list;
        for(int node = 0; node < int(flow_path.size()); node ++)
        {
          final_channel_list.push_back(flow_path[node]);

          // get the flow distance of the point, as well as the x and y locations
          int curr_row,curr_col;
          float curr_E,curr_N;
          FlowInfo.retrieve_current_row_and_col(flow_path[node],curr_row,curr_col);
          spaced_distances.push_back( FlowDistance.get_data_element(curr_row,curr_col));
          FlowInfo.get_x_and_y_from_current_node(flow_path[node],curr_E, curr_N);
          spaced_eastings.push_back( curr_E );
          spaced_northings.push_back( curr_N );


          // If this is the downstream node, jump to the end of the flowpath
          if ( flow_path[node] == downstream_NI)
          {
            node = flow_path.size();
          }
        }
        LSDRaster DA_d8 = FlowInfo.write_DrainageArea_to_LSDRaster();

        string filename = OUT_ID+"_swath_channel";
        FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(final_channel_list,
                          OUT_DIR, filename,filled_topography, FlowDistance,
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
      else
      {
        cout << "The channel file has more than two nodes. " << endl;
        cout << "I am assuming this is a fully functional channel file. " << endl;
        cout << "I need to look at this file for the correct header (flow_distance)" << endl;
        cout << "If it doesn't have it, the file isn't correctly formatted." << endl;

        CSVFile.get_x_and_y_from_latlong(spaced_eastings,spaced_northings); 
        string fd_column_name = "flowdistance(m)";

        if ( CSVFile.is_column_in_csv(fd_column_name) )
        {
          cout << "Super, I've go the distances in the file." << endl;
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
      }
    }


    // Now for the swath!
    // The categorised swath is mainly used with terraces. 
    // See lsdtt-valley-metrics to get the terraces
    if (this_bool_map["make_categorised_swath"])
    {
      // first we need to open the categorised raster
      string category_infile_prefix = DATA_DIR+this_string_map["swath_category_raster_prefix"];
      LSDIndexRaster categories(category_infile_prefix,raster_ext);

      string swath_data_prefix = OUT_DIR+OUT_ID+"_categorised";
      topography_raster.make_swaths_from_categorised(spaced_eastings, spaced_northings,spaced_distances, this_float_map["swath_width"],
                                 this_float_map["swath_bin_spacing"], swath_data_prefix, this_bool_map["print_swath_rasters"],
                                 categories);

    }
    else
    {
      string swath_data_prefix = OUT_DIR+OUT_ID;
      topography_raster.make_swath(spaced_eastings, spaced_northings,spaced_distances, this_float_map["swath_width"],
                                 this_float_map["swath_bin_spacing"], swath_data_prefix, this_bool_map["print_swath_rasters"]);
    }



  }


  //=========================================================================
  // .##..##...####...##......##......######..##..##..........##...##..######..#####...######..##..##.
  // .##..##..##..##..##......##......##.......####...........##...##....##....##..##....##....##..##.
  // .##..##..######..##......##......####......##............##.#.##....##....##..##....##....######.
  // ..####...##..##..##......##......##........##............#######....##....##..##....##....##..##.
  // ...##....##..##..######..######..######....##.............##.##...######..#####.....##....##..##.
  //
  //=========================================================================
  // Below is just a debugging routine
  if(this_bool_map["test_bearing_template"])
  {
    float i = this_float_map["test_bearing"];
    cout << "Bearing is: " << i << endl;
    Array2D<int> test = make_template_for_vector_bearing(float(i), this_int_map["test_scale"]);
    //for (int i = 316; i<=360; i++)
    //{
    //  cout << "Bearing is: " << i << endl;
    //  
    // }
  }
  
  
  if(this_bool_map["channel_and_valley_width_extraction"])
  {
    cout << "I am going to test some code for valley and channel widths" << endl;


    Array2D<int> test = make_template_for_vector_bearing(this_float_map["test_bearing"], this_int_map["test_scale"]);

    // now I am going to use a channel to get the bearings. 
    // I need to load the channel file
    cout << "I am loading a channel file. " << endl;
    LSDRasterInfo RI(topography_raster);
    LSDSpatialCSVReader CSVFile(RI,this_string_map["valley_points_csv"]);

    LSDRaster filled_topography,carved_topography;
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
      if(this_bool_map["carve_before_fill"])
      {
        carved_topography = topography_raster.Breaching_Lindsay2016();
        filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
      }
      else
      {
        filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
      }
    }
    LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);
    vector<int> node_list =  CSVFile.get_nodeindices_from_lat_long(FlowInfo);

    vector<int> bearing_nodes;
    vector<float> bearing_vec;
    string bearing_fname =  OUT_DIR+OUT_ID+"_bearings.csv";
    FlowInfo.calculate_bearings_from_nodelist(node_list, this_int_map["channel_bearing_node_spacing"], 
                                              this_bool_map["print_channel_bearings"], 
                                              bearing_fname,bearing_nodes, bearing_vec);

    if ( this_bool_map["print_channel_bearings"] && this_bool_map["convert_csv_to_geojson"])
    {
      cout << "Printing the geojson of the channel_bearings" << endl;
      string gjson_name = OUT_DIR+OUT_ID+"_bearings.geojson";
      LSDSpatialCSVReader thiscsv(bearing_fname);
      thiscsv.print_data_to_geojson(gjson_name);
    }

  }

  //=================================================================================================
  // Some tools for punching out some nodata using a secondary raster that we can use to make river
  // pathways in noisy DEMs.
  //
  //.#####...##..##..##..##...####...##..##..........##..##...####...#####....####...######...####..
  //.##..##..##..##..###.##..##..##..##..##..........###.##..##..##..##..##..##..##....##....##..##.
  //.#####...##..##..##.###..##......######..........##.###..##..##..##..##..######....##....######.
  //.##......##..##..##..##..##..##..##..##..........##..##..##..##..##..##..##..##....##....##..##.
  //.##.......####...##..##...####...##..##..........##..##...####...#####...##..##....##....##..##. 
  //
  //=================================================================================================
  if(this_bool_map["punch_nodata"])
  {
    cout << "I am going to punch out some nodata. " << endl;
    
    string punch_R_in = DATA_DIR+this_string_map["punch_raster_prefix"];
    cout << "The name of the raster I will use to punch out the nodata is " << punch_R_in+raster_ext << endl;
    cout << "I will change to nodata any pixels in your raster that have a value in the punching raster";
    if (this_bool_map["belowthresholdisnodata"])
    {
      cout << " below ";
    }
    else
    {
      cout << " above ";
    }
    cout << " the threshold " << this_float_map["punch_threshold"] << endl;
    
    LSDRaster punch_raster_in(punch_R_in,raster_ext);

    LSDRaster Punched = topography_raster.mask_to_nodata_using_threshold_using_other_raster(this_float_map["punch_threshold"],this_bool_map["belowthresholdisnodata"], punch_raster_in);

    string punch_name = OUT_DIR+OUT_ID+"_Punched";
    cout << "The punched raster is: " << punch_name << endl;
    Punched.write_raster(punch_name,raster_ext);

  }

  //================================================================================================
  //
  //.##..##...####...##......##......######..##..##...........####....####...#####...##..##..######.
  //.##..##..##..##..##......##......##.......####...........##..##..##..##..##..##..##..##..##.....
  //.##..##..######..##......##......####......##............##......######..#####...##..##..####...
  //..####...##..##..##......##......##........##............##..##..##..##..##..##...####...##.....
  //...##....##..##..######..######..######....##.............####...##..##..##..##....##....######.
  //
  //================================================================================================  
  if(this_bool_map["create_valley_trough"] || this_bool_map["get_valley_centreline"])
  {
    string punch_name;
    if(this_bool_map["punch_nodata"])
    {
      cout << "You punched a raster already so I am using that instead of the" << endl;
      cout << DATA_DIR+this_string_map["Punched_raster_prefix"]+".bil" << endl;
      punch_name = OUT_DIR+OUT_ID+"_Punched";
      cout << "The punch raster I'm using is: " <<  punch_name << endl;
    }
    else
    {

      punch_name = DATA_DIR+this_string_map["Punched_raster_prefix"];
      cout << "I am going to use the raster" << endl;
      cout << punch_name << ".bil as the punched raster." << endl;  
    }
           
    LSDRaster Punched(punch_name,raster_ext);
    
    cout << "Grabbing rasters giving the distance and elevation of the nearest" << endl;
    cout << "pixels in the punched area to pixels with data." << endl;
    vector<LSDRaster> nearest_rs = Punched.get_nearest_distance_and_value_masks(topography_raster);

    string distance_R_name = OUT_DIR+OUT_ID+"_DistToND";
    string value_R_name = OUT_DIR+OUT_ID+"_ValueToND";    

    nearest_rs[0].write_raster(distance_R_name,raster_ext);
    nearest_rs[1].write_raster(value_R_name,raster_ext);

    // now get the lowest value in a window
    cout << "Now get the lowest value in a window." << endl;
    float window_radius = this_float_map["minimum_bank_elevation_window_radius"];
    int window_type = 1;
    bool find_maximum = false;

    cout << "Finding the lowest pixels along the river or valley edge." << endl;
    LSDRaster minimum_river_elev = nearest_rs[1].neighbourhood_statistics_local_min_max(window_radius, window_type, find_maximum);
    string minelev_R_name = OUT_DIR+OUT_ID+"_MinBankElev";
    minimum_river_elev.write_raster(minelev_R_name,raster_ext);

    if(this_bool_map["create_valley_trough"] || this_bool_map["get_valley_centreline"])
    {
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
        if (i < n_loops-2)
        {
          add_trough = add_trough.MapAlgebra_add(trough);
        }
      }

      river_path = river_path.MapAlgebra_add(add_trough);
      river_path.AdjustElevation(this_float_map["river_depth"]);


      // Now merge the rasters
      cout << "Combining rasters." << endl;
      LSDRaster topo_copy = topography_raster;
      topo_copy.OverwriteRaster(river_path);

      string imposedriver_R_name = OUT_DIR+OUT_ID+"_RiverTrough";
      topo_copy.write_raster(imposedriver_R_name,raster_ext);    

      if(this_bool_map["write_hillshade"])
      {
        cout << "Let me print the hillshade for you. " << endl;
        float hs_azimuth = 315;
        float hs_altitude = 45;
        float hs_z_factor = 1;
        LSDRaster hs_raster = topo_copy.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

        string hs_fname = OUT_DIR+OUT_ID+"_RiverTrough_hs";
        hs_raster.write_raster(hs_fname,raster_ext);
      }


      if (this_bool_map["get_valley_centreline"])
      {
        
        LSDRaster carved_topography = topo_copy.Breaching_Lindsay2016();
        topo_copy = carved_topography.fill(this_float_map["min_slope_for_fill"]);

        // Now we need the flow path
        cout << "\t Flow routing for centreline. Note this is memory intensive. If your DEM is very large you may get a segmentation fault here..." << endl;
        // get a flow info object
        LSDFlowInfo FI(boundary_conditions,topo_copy);
        cout << "Finished flow routing." << endl; 

        // get the maximum elevation
        int max_row,max_col;
        float max_elev = river_path.max_elevation(max_row,max_col);

        int start_ni = FI.get_NodeIndex_from_row_col(max_row,max_col);

        vector<int> flow_path = FI.get_flow_path(start_ni);

        cout << "I need to calculate the flow distance now." << endl;
        LSDRaster FD = FI.distance_from_outlet();
        LSDRaster DA = FI.write_DrainageArea_to_LSDRaster();

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



  }


  if(this_bool_map["print_nearest_to_nodata_rasters"])
  {
    cout << "I am going to find the nearest value to any nodata pixels, and the distance to those pixels. " << endl;
    vector<LSDRaster> nearest_rs = topography_raster.get_nearest_distance_and_value_masks();

    string distance_R_name = OUT_DIR+OUT_ID+"_DistToND";
    string value_R_name = OUT_DIR+OUT_ID+"_ValueToND";

    nearest_rs[0].write_raster(distance_R_name,raster_ext);
    nearest_rs[1].write_raster(value_R_name,raster_ext);

    // now get the lowest value in a window
    float window_radius = this_float_map["minimum_bank_elevation_window_radius"];
    int window_type = 1;
    bool find_maximum = false;

    LSDRaster minimum_river_elev = nearest_rs[1].neighbourhood_statistics_local_min_max(window_radius, window_type, find_maximum);
    string minelev_R_name = OUT_DIR+OUT_ID+"_MinBankElev";
    minimum_river_elev.write_raster(minelev_R_name,raster_ext);

    // drop the river by the river depth
    minimum_river_elev.AdjustElevation(this_float_map["river_depth"]);

    // Now merge the rasters
    LSDRaster topo_copy = topography_raster;
    topo_copy.OverwriteRaster(minimum_river_elev);

    string imposedriver_R_name = OUT_DIR+OUT_ID+"_ImposedRiver";
    topo_copy.write_raster(imposedriver_R_name,raster_ext);    

    if(this_bool_map["write_hillshade"])
    {
      cout << "Let me print the hillshade for you. " << endl;
      float hs_azimuth = 315;
      float hs_altitude = 45;
      float hs_z_factor = 1;
      LSDRaster hs_raster = topo_copy.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

      string hs_fname = OUT_DIR+OUT_ID+"_ImposedRiver_hs";
      hs_raster.write_raster(hs_fname,raster_ext);
    }

  }


  if(this_bool_map["print_relief_raster"])
  {
    cout << "Calculating the relief raster bit." << endl;
    LSDRaster relief = topography_raster.calculate_relief(this_float_map["relief_window"], this_int_map["relief_window_kernel_type"]);
    string this_raster_name = OUT_DIR+OUT_ID+"_BASICRELIEF";
    relief.write_raster(this_raster_name,raster_ext);
  }


  //============================================================================
  //
  //.#####....####....####...######..######..#####...........#####...##..##..#####...##..##.
  //.##..##..##..##..##........##....##......##..##..........##..##..##..##..##..##..###.##.
  //.#####...######...####.....##....####....#####...........#####...##..##..#####...##.###.
  //.##..##..##..##......##....##....##......##..##..........##..##..##..##..##..##..##..##.
  //.##..##..##..##...####.....##....######..##..##..........#####....####...##..##..##..##.
  //
  //============================================================================
  // if this gets burned, do it
  if(this_bool_map["burn_raster_to_csv"])
  {
    cout << "You asked me to burn a raster to a csv." << endl;
    cout << "WARNING: This was written in a hurry and has no bug checking." << endl;
    cout << "If you have the wrong filenames it will crash." << endl;

    cout << "First I am going to load the raster." << endl;
    string burn_raster_name;
    bool burn_raster_exists = true;
    if (this_string_map["burn_raster_prefix"] != "NULL")
    {
      burn_raster_name = DATA_DIR+this_string_map["burn_raster_prefix"];
      cout << "I will burn data from the raster: " << burn_raster_name << endl;
    }
    else
    {
      burn_raster_exists = false;
      cout << "You don't have a working burn raster." << endl;
    }

    if(burn_raster_exists)
    {
      LSDRaster BurnRaster(burn_raster_name,raster_ext);

      string header_for_burn_data;
      header_for_burn_data = this_string_map["burn_data_csv_column_header"];

      cout << "I am burning the raster into the column header " << header_for_burn_data << endl;

      string full_csv_name = DATA_DIR+this_string_map["csv_to_burn_name"];
      cout << "I am burning the raster to the csv file: " << full_csv_name << endl;
      LSDSpatialCSVReader CSVFile(RI,full_csv_name);

      cout << "Now let me burn that raster for you." << endl;
      CSVFile.burn_raster_data_to_csv(BurnRaster,header_for_burn_data);

      string full_burned_csv_name = OUT_DIR+OUT_ID+"_burned.csv";
      cout << "Now I'll print the data to a new file, the file is: " << full_burned_csv_name << endl;
      CSVFile.print_data_to_csv(full_burned_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_csv_burned.geojson";
        LSDSpatialCSVReader thiscsv(full_burned_csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
    else
    {
      cout << "Burn raster doesn't exist so I didn't do anything." << endl;
    }
  }
  //============================================================================




  //============================================================================
  //
  //.##..##..######..##......##.......####...##..##...####...#####...######.
  //.##..##....##....##......##......##......##..##..##..##..##..##..##.....
  //.######....##....##......##.......####...######..######..##..##..####...
  //.##..##....##....##......##..........##..##..##..##..##..##..##..##.....
  //.##..##..######..######..######...####...##..##..##..##..#####...######.
  //
  //============================================================================
  // Compute hillshade and print
  //============================================================================
  if (this_bool_map["write_hillshade"])
  {
    cout << "Let me print the hillshade for you. " << endl;
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+OUT_ID+"_hs";
    hs_raster.write_raster(hs_fname,raster_ext);
  }

  //============================================================================
  //
  //..####...##..##..#####...######...####....####...######..........######..######..######..######..######..##..##...####..
  //.##......##..##..##..##..##......##..##..##..##..##..............##........##......##......##......##....###.##..##.....
  //..####...##..##..#####...####....######..##......####............####......##......##......##......##....##.###..##.###.
  //.....##..##..##..##..##..##......##..##..##..##..##..............##........##......##......##......##....##..##..##..##.
  //..####....####...##..##..##......##..##...####...######..........##......######....##......##....######..##..##...####..
  //
  //============================================================================
  // The surface fitting metrics
  //============================================================================
  vector<int> raster_selection(9, 0);  // This controls which surface fitting metrics to compute
  bool doing_polyfit = false;
  if(this_bool_map["print_smoothed_elevation"])
  {
    raster_selection[0] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_slope"])
  {
    raster_selection[1] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_aspect"])
  {
    raster_selection[2] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_curvature"])
  {
    raster_selection[3] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_planform_curvature"])
  {
    raster_selection[4] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_profile_curvature"])
  {
    raster_selection[5] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_tangential_curvature"])
  {
    raster_selection[6] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_point_classification"])
  {
    raster_selection[7] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_directional_gradients"])
  {
    raster_selection[8] = 1;
    doing_polyfit = true;
  }
  if(this_bool_map["print_channel_data_plus_surface_metrics"])
  {
    raster_selection[1] = 1;
    raster_selection[3] = 1;
    doing_polyfit = true;
  }

  // We place the surface fitting vector outside because we will need this information
  // later for the basin statistics.
  vector<LSDRaster> surface_fitting;
  if (doing_polyfit)
  {
    cout << "I am running the polyfit function. This could take some time." << endl;

    surface_fitting = topography_raster.calculate_polyfit_surface_metrics_directional_gradients(this_float_map["surface_fitting_radius"], raster_selection);
    if(this_bool_map["print_smoothed_elevation"])
    {
      cout << "Let me print the smoothed elevation raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_SMOOTH";
      surface_fitting[0].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_slope"])
    {
      cout << "Let me print the slope raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_SLOPE";
      surface_fitting[1].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_aspect"])
    {
      cout << "Let me print the aspect raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_ASPECT";
      surface_fitting[2].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_curvature"])
    {
      cout << "Let me print the curvature raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_CURV";
      surface_fitting[3].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_planform_curvature"])
    {
      cout << "Let me print the planform curvature raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_PLFMCURV";
      surface_fitting[4].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_profile_curvature"])
    {
      cout << "Let me print the profile curvature raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_PROFCURV";
      surface_fitting[5].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_tangential_curvature"])
    {
      cout << "Let me print the tangential curvature raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_TANCURV";
      surface_fitting[6].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_point_classification"])
    {
      cout << "Let me print the point classification curvature raster for you."  << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_CLASS";
      surface_fitting[7].write_raster(this_raster_name,raster_ext);
    }
    if(this_bool_map["print_directional_gradients"])
    {
      cout << "Let me print the directional gradient rasters for you."  << endl;
      string this_raster_name1 = OUT_DIR+OUT_ID+"_DDX";
      string this_raster_name2 = OUT_DIR+OUT_ID+"_DDY";
      surface_fitting[8].write_raster(this_raster_name1,raster_ext);
      surface_fitting[9].write_raster(this_raster_name2,raster_ext);
    }
  }
  else
  {
    cout << "I won't be doing any polyfitting so you have saved yourself some time!" << endl;
  }

  if (this_bool_map["calculate_window_size"])
  {
    cout << "I am going to run the calculations so you can determine the window size." << endl;
    cout << "This is extremely computationally expensive. I suggest finding a good book." << endl;

    vector<int> raster_selection(9, 0);  // This controls which surface fitting metrics to compute
    raster_selection[3] = 1;    // Turns on the curvature

    //array of window sizes to test
    //this is a little arbritrary but reflects the sizes used by Roering et al. 2010
    int WindowSizes[] = {1, 2, 3, 4, 5, 7, 8, 10,15,20,25,50,70,90,100};

    // floats to hold the stats about the fitted surface
    float Curv_mean;
    float Curv_stddev;
    float Curv_median;
    float Curv_UpperQuartile;
    float Curv_LowerQuartile;
    float Curv_MaxValue;
    float Curv_iqr;
    vector<float> Curv_vec;

    string this_window_fitting_name = OUT_DIR+OUT_ID+"_WindowSize.csv";
    ofstream WriteData;
    WriteData.open(this_window_fitting_name.c_str());

    //write headers
    WriteData << "Length_scale,Curv_mean,Curv_stddev,Curv_iqr" << endl;

    for (int w = 0; w < 15; ++w){

      cout << "Processing surface " << w+1 << " of " << "15" << endl;

      vector<LSDRaster> Surfaces = topography_raster.calculate_polyfit_surface_metrics(WindowSizes[w], raster_selection);
      LSDRaster curvature = Surfaces[3];

      //reset values for next run
      Curv_mean = 0;
      Curv_stddev = 0;
      Curv_median = 0;
      Curv_UpperQuartile = 0;
      Curv_LowerQuartile = 0;
      Curv_MaxValue = 0;
      Curv_iqr = 0;
      Curv_vec.clear();

      //go through the landscape and get every curvature value into a 1D vector
      for (int i = 0; i < int(curvature.get_NRows()); ++i){
        for (int j = 0; j < int(curvature.get_NCols()); ++j){
          if (curvature.get_data_element(i,j) != curvature.get_NoDataValue()){
            Curv_vec.push_back(curvature.get_data_element(i,j));
          }
        }
      }

      //calculate the std dev and iqr (this comes from LSDStatsTools)
      get_distribution_stats(Curv_vec, Curv_mean, Curv_median, Curv_UpperQuartile, Curv_LowerQuartile, Curv_MaxValue);
      Curv_stddev = get_standard_deviation(Curv_vec, Curv_mean);
      Curv_iqr = Curv_UpperQuartile - Curv_LowerQuartile;

      //write the values to the output file
      WriteData << WindowSizes[w] << "," << Curv_mean << "," << Curv_stddev << "," << Curv_iqr << endl;

      }

    WriteData.close();

  }

  //============================================================================
  //
  //.#####....####...##..##...####...##..##..##..##..######...####....####..
  //.##..##..##..##..##..##..##......##..##..###.##..##......##......##.....
  //.#####...##..##..##..##..##.###..######..##.###..####.....####....####..
  //.##..##..##..##..##..##..##..##..##..##..##..##..##..........##......##.
  //.##..##...####....####....####...##..##..##..##..######...####....####..
  //
  //============================================================================
  if(this_bool_map["print_REI_raster"])
  {
    cout << "I am caluclating the REI raster." << endl;
    LSDRaster REI_raster = topography_raster.calculate_REI(this_float_map["REI_window_radius"], this_float_map["REI_critical_slope"]);
    string this_raster_name = OUT_DIR+OUT_ID+"_REI";
    REI_raster.write_raster(this_raster_name,raster_ext);
  }

  if(this_bool_map["print_roughness_rasters"])
  {
    cout << "I am printing the roughness rasters S1, S2, and S3. " << endl;
    cout << "See Milodowski et al ESURF 2015 https://doi.org/10.5194/esurf-3-483-2015" << endl;
    // This just ensures all three roughness rasters are printed
    vector<int> file_code;
    file_code.push_back(1);  file_code.push_back(1); file_code.push_back(1);
    string file_prefix = OUT_DIR+OUT_ID;

    // The roughness radius should be around 3m for a 1m DEM: you get more diffuse "rocks with wider radii"
    topography_raster.calculate_roughness_rasters(this_float_map["surface_fitting_radius"], this_float_map["roughness_radius"],file_prefix,file_code);

  }



  //============================================================================
  // Print the wiener filtered raster if that is what you want
  //============================================================================
  if (this_bool_map["print_wiener_filtered_raster"])
  {

    cout << "I am running a filter to print to raster." << endl;
    cout << "This uses spectral analysis and is memory intensive. If you are working on a system with limited memory, " << endl;
    cout << "you may get a segmentation fault here!" << endl;
    LSDRasterSpectral SpectralRaster(topography_raster);
    LSDRaster topo_test_wiener = SpectralRaster.fftw2D_wiener();
    cout << "Finished getting the filtered raster. " << endl;

    string wiener_name = OUT_DIR+OUT_ID+"_Wfilt";
    topo_test_wiener.write_raster(wiener_name,raster_ext);
  }

  //============================================================================================================
  //
  //..####...######..##..##...####...##......######...........####...##..##...####...##..##..##..##..######..##.....
  //.##........##....###.##..##......##......##..............##..##..##..##..##..##..###.##..###.##..##......##.....
  //..####.....##....##.###..##.###..##......####............##......######..######..##.###..##.###..####....##.....
  //.....##....##....##..##..##..##..##......##..............##..##..##..##..##..##..##..##..##..##..##......##.....
  //..####...######..##..##...####...######..######...........####...##..##..##..##..##..##..##..##..######..######.
  //
  //============================================================================================================
  if (this_bool_map["impose_single_channel"])
  {
    cout << "I am going to impose a channel on your DEM." << endl;
    cout << "This will overwrite these pixels. It also fills nodata." << endl;
    LSDRasterInfo RI(topography_raster);
    cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
    LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );
    single_channel_data.print_data_map_keys_to_screen();

    cout << "Subtracting " << this_float_map["single_channel_drop"] << " m elevation from the single channel" << endl;
    single_channel_data.data_column_add_float(this_string_map["single_channel_elev_string"], -this_float_map["single_channel_drop"]);

    if (this_bool_map["force_single_channel_slope"])
    {
      single_channel_data.enforce_slope(this_string_map["single_channel_fd_string"], this_string_map["single_channel_elev_string"], this_float_map["min_slope_for_fill"]);
    }

    LSDRasterMaker RM(topography_raster);
    if(this_bool_map["buffer_single_channel"] and this_bool_map["use_XY_for_buffer"])
    {
      cout << "I am going to impose and buffer the single channel using XY coordinates." << endl;
      RM.impose_channels_with_buffer_use_XY(single_channel_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);
    }
    else if(this_bool_map["buffer_single_channel"] and not this_bool_map["use_XY_for_buffer"])
    {
      cout << "I am going to impose and buffer the single channel using lat-long coordinates." << endl;
      RM.impose_channels_with_buffer(single_channel_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);
    }
    else
    {
      cout << "I am going to impose the single channel using lat-long coordinates, without buffering." << endl;
      RM.impose_channels(single_channel_data, this_string_map["single_channel_elev_string"]);
    }
    LSDRaster imposed_topography = RM.return_as_raster();

    cout << "Let me print the fill raster for you."  << endl;
    string filled_raster_name = OUT_DIR+OUT_ID+"_imposed_channels";
    imposed_topography.write_raster(filled_raster_name,raster_ext);       
  }

  if (this_bool_map["centre_and_interpolate_channel_coordinates"])
  {
    cout << "I am going to centre your channel nodes on the corresponding raster pixel and interpolate across any gaps." << endl;
    cout << "This will overprint your existing XY coordinates in the channel csv file. I am assuming your channel and DEM are in UTM." << endl;
    cout << "WARNING: Make sure your channel data is sorted by flow distance, or your final channel will look like a pretzel." << endl;  

    LSDRasterInfo RI(topography_raster);

    // load the channel nodes csv file, this will also print all the column names to screen
    LSDSpatialCSVReader channel_nodes(RI, (this_string_map["fixed_channel_csv_name"]));

    // assign directory and name of output files
    string centred_csv_name = OUT_DIR+"centred_"+this_string_map["fixed_channel_csv_name"];
    string interpolated_csv_name = OUT_DIR+"interpolated_"+this_string_map["fixed_channel_csv_name"];

    string X_column_name = this_string_map["X_column_name"];
    string Y_column_name = this_string_map["Y_column_name"];
    string fd_column_name = this_string_map["single_channel_fd_string"];

    channel_nodes.centre_XY_by_row_col(topography_raster, X_column_name, Y_column_name);

    channel_nodes.print_data_to_csv(centred_csv_name);   
    cout << "I have saved your file as " << centred_csv_name << endl;

    channel_nodes.interpolate_across_gap(topography_raster, X_column_name, Y_column_name, fd_column_name);

    cout << "I have interpolated, let's centre the coordinates again." << endl;
    channel_nodes.centre_XY_by_row_col(topography_raster, X_column_name, Y_column_name);

    // now print the data with the gaps filled in
    channel_nodes.print_data_to_csv(interpolated_csv_name);
    cout << "I have saved your file as " << interpolated_csv_name << ". Good night." << endl;

  }

  if(this_bool_map["prepare_single_channel_from_xy"])
  {
    cout << "I am preparing a fixed channel based only on xy data and elevation." << endl;
    // In this instance we only have xy data without flow distance or lat long. 
    // This is the kind of data supplied by Andreas from Nagra. 
    // These pixels have elevations, so use this to order them
    // I need to 1) centre the pixels.
    // 2) get the lat-long
    // 3) burn the channel in
    // 4) get the flow distance
    // 5) then print this new file out
    LSDRasterInfo RI(topography_raster);

    // load the channel nodes csv file, this will also print all the column names to screen
    LSDSpatialCSVReader channel_nodes(RI, (this_string_map["fixed_channel_csv_name"]));

    // Add dummy latitude and longitude columns
    channel_nodes.add_placeholder_latitude_and_longitude();

    // assign directory and name of output files
    string centred_csv_name = OUT_DIR+"centred_"+this_string_map["fixed_channel_csv_name"];
    string interpolated_csv_name = OUT_DIR+"interpolated_"+this_string_map["fixed_channel_csv_name"];

    string X_column_name = this_string_map["X_column_name"];
    string Y_column_name = this_string_map["Y_column_name"];
    string fd_column_name = this_string_map["single_channel_fd_string"];

    // this gets the x,y, centres them on the nodes, and then recalculates latitude and longitude
    channel_nodes.centre_XY_by_row_col(topography_raster, X_column_name, Y_column_name);

    // now we find the minimum elevation of the topography raster and drop the channel below this
    //float minimum_elevation = topography_raster.min_elevation();
    
    cout << "Subtracting " << this_float_map["single_channel_drop"] << " m elevation from the single channel" << endl;
    channel_nodes.data_column_add_float(this_string_map["single_channel_elev_string"], -this_float_map["single_channel_drop"]);

    // now impose these pixels
    LSDRasterMaker RM(topography_raster);
    RM.impose_channels(channel_nodes, this_string_map["single_channel_elev_string"]);
    LSDRaster imposed_topography = RM.return_as_raster();

    // Now get flow routing
    LSDRaster carved = imposed_topography.Breaching_Lindsay2016();
    LSDRaster filled_topography = carved.fill(this_float_map["min_slope_for_fill"]);

    // Get the flow routing and and flow distance
    LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);
    LSDRaster FD = FlowInfo.distance_from_outlet();
    string header_for_fd = this_string_map["single_channel_fd_string"];

    // burn the flow distance to csv
    channel_nodes.burn_raster_data_to_csv(FD,header_for_fd);

    // put the channel back up
    channel_nodes.data_column_add_float(this_string_map["single_channel_elev_string"], this_float_map["single_channel_drop"]);

    // print the centred csv that now includes flow distance
    channel_nodes.print_data_to_csv(centred_csv_name);    

    // now interpolate
    channel_nodes.interpolate_across_gap(topography_raster, X_column_name, Y_column_name, fd_column_name);

    cout << "I have interpolated, let's centre the coordinates again." << endl;
    channel_nodes.centre_XY_by_row_col(topography_raster, X_column_name, Y_column_name);

    // now print the data with the gaps filled in
    channel_nodes.print_data_to_csv(interpolated_csv_name);

    cout << "I am all finished prepping some x,y data to serve as a baselevel channel." << endl;
    cout << "I now need to exit because I did some raster manipulation that would mess up any other analysis." << endl;
    cout << "This is a tough computation and I think you should treat yourself to some chocolate." << endl;
    
    exit(0);

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
  //
  // EVERTHING BELOW THIS POINT NEEDS A FILL RASTER AND FLOW ROUTING
  // THIS IS WHERE MEMORY CONSUMPTION BECOMES A PROBLEM
  // The LSDFlowInfo object is ~10-20x as big as the original DEM so
  // if the DEM is really big or you are working on a computer
  // with limited memory you might get segmentation faults after this point
  //
  //============================================================================
  if (this_bool_map["print_dinf_drainage_area_raster"]
        || this_bool_map["print_d8_drainage_area_raster"]
        || this_bool_map["print_QuinnMD_drainage_area_raster"]
        || this_bool_map["print_FreemanMD_drainage_area_raster"]
        || this_bool_map["print_MD_drainage_area_raster"]
        || this_bool_map["print_fill_raster"]
        || this_bool_map["print_stream_order_raster"]
        || this_bool_map["print_channels_to_csv"]
        || this_bool_map["print_junction_index_raster"]
        || this_bool_map["print_junctions_to_csv"]
        || this_bool_map["find_basins"]
        || this_bool_map["print_channel_tips_raster"]
        || this_bool_map["print_chi_data_maps"]
        || this_bool_map["print_junction_angles_to_csv"]
        || this_bool_map["extract_single_channel"]
        || this_bool_map["remove_nodes_influenced_by_edge"]
        || this_bool_map["isolate_pixels_draining_to_fixed_channel"]
        || this_bool_map["calculate_basin_statistics"]
        || this_bool_map["divide_finder"]
        || this_bool_map["tag_nodes"]
        || this_bool_map["tag_nodes_from_categorised_channel"]
        || this_bool_map["extract_ridges"]
        || this_bool_map["calculate_connectivity_index"]
        || this_bool_map["calculate_hypsometric_integral"]
        || this_bool_map["clip_raster_to_basins"])
  {
    cout << "I will need to compute flow information, because you are getting drainage area or channel networks." << endl;
    //==========================================================================
    // Fill the raster
    //==========================================================================

    if(this_bool_map["fill_interior_nodata"])
    {
      cout << "I will try to get rid of pixels in the middle of your raster that are nodata." << endl;
      LSDRaster new_topography = topography_raster.alternating_direction_nodata_fill_irregular_raster(this_int_map["window_for_filling_interior_nodata"]);
      topography_raster = new_topography;
    }

    LSDRaster filled_topography,carved_topography;
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
      if(this_bool_map["carve_before_fill"])
      {
        carved_topography = topography_raster.Breaching_Lindsay2016();
        filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
      }
      else
      {
        filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
      }
    }

    if (this_bool_map["print_fill_raster"])
    {
      cout << "Let me print the fill raster for you."  << endl;
      string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
      filled_topography.write_raster(filled_raster_name,raster_ext);
    }
    //==========================================================================

    // This routine makes sure that the outer edges of the DEM
    // are trimmed to a rectangle and then the edge nodes are replaced
    // with nodata
    if (this_bool_map["remove_nodes_influenced_by_edge"])
    {
      if (this_bool_map["remove_outer_nodes_for_edge_influence"])
      {
        if( this_bool_map["use_spiral_trimmer_for_edge_influence"])
        {
          LSDRaster SpiralTrim = filled_topography.RasterTrimmerSpiral();
          filled_topography = SpiralTrim;
        }
        cout << "Warning, I have replaced the edges of your DEM with nodata" << endl;
        cout << "This usually means I am trying to remove nodes influence by the edge" << endl;
        filled_topography.replace_edges_with_nodata();
      }
    }

    cout << "\t Flow routing. Note this is memory intensive. If your DEM is very large you may get a segmentation fault here..." << endl;
    // get a flow info object
    LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);
    cout << "Finished flow routing." << endl;

    if (this_bool_map["remove_nodes_influenced_by_edge"])
    {
      cout << "I am going to print a raster that has nodata for all nodes influenced by the edge or nodata" << endl;
      LSDRaster NodesRemovedRaster = FlowInfo.remove_nodes_influneced_by_edge(filled_topography);
      string remove_raster_name = OUT_DIR+OUT_ID+"_NoEdge";
      NodesRemovedRaster.write_raster(remove_raster_name,raster_ext);
    }

    if(this_bool_map["isolate_pixels_draining_to_fixed_channel"])
    {
      // first read the
      // Get the latitude and longitude
      cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
      LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );
      vector<float> X_coords; 
      vector<float> Y_coords; 
      vector<int> nodes_from_channel;

      if(this_bool_map["use_xy_for_node_index"])
      {
        cout << "I am going to get node indices from X-Y coordinates." << endl;
        
        source_points_data.get_nodeindices_from_x_and_y_coords(FlowInfo, X_coords, Y_coords, nodes_from_channel);        
      } else{
        cout << "I am going to get node indices from lat-long coordinates." << endl;
        nodes_from_channel = source_points_data.get_nodeindices_from_lat_long(FlowInfo); 
      }
      

      //vector<int> new_nodes_from_channel = source_points_data.get_nodeindex_vector();
      //cout << "Old nodes: " <<  new_nodes_from_channel.size() << " and new: " << nodes_from_channel.size() << endl;
      //for (int i = 0; i< int(nodes_from_channel.size()); i++)
      //{
      //  if (nodes_from_channel[i] < 0)
      //  {
      //    cout << "Invalid node index: " << nodes_from_channel[i] << endl;
      //  }
      //}

      // Now run the flowinfo routine
      LSDRaster NodesRemovedRaster = FlowInfo.find_nodes_not_influenced_by_edge_draining_to_nodelist(nodes_from_channel,filled_topography);
      string remove_raster_name = OUT_DIR+OUT_ID+"_IsolateFixedChannel";
      NodesRemovedRaster.write_raster(remove_raster_name,raster_ext);

    }

    if (this_bool_map["tag_nodes_from_categorised_channel"])
    {
      cout << "Hello, I am going to take a channel file, categorise some quantity, and then tag nodes on that basis" << endl;
      cout << "This is probably most useful for divide finding." << endl;

      // now load the channel to be used
      LSDRasterInfo RI(topography_raster);
      cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
      LSDSpatialCSVReader channel_to_categorise( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

      // now categorise the something
      cout << "The category boundaries string is: " << this_string_map["category_boundaries"] << endl;
      vector<double> new_boundary_vec = LSDPP.parse_double_vector( "category_boundaries" );
      for(int i = 0; i<int(new_boundary_vec.size()); i++)
      {
        cout << "bound val: " << new_boundary_vec[i] << endl;
      }

      string new_column_name = this_string_map["column_to_categorise"] + "_categories";
      channel_to_categorise.data_column_add_categorised(this_string_map["column_to_categorise"], 
                                                        new_column_name, new_boundary_vec);
      

      // Now we get the nodeindex column
      channel_to_categorise.add_nodeindex_vector_to_data_map_using_lat_long(FlowInfo);

      // Now we make a map of the values
      vector<int> ni_vec = channel_to_categorise.data_column_to_int("nodeindex");
      vector<int> cat_vec = channel_to_categorise.data_column_to_int(this_string_map["column_to_categorise"] + "_categories");

      LSDIndexRaster tagged_raster = FlowInfo.tag_donor_pixels_exclude_nodelist(ni_vec,cat_vec);

      // print the tagged raster
      string tag_raster_name = OUT_DIR+OUT_ID+"_channel_tagged_pixels";
      tagged_raster.write_raster(tag_raster_name,raster_ext);

      // Print the channel for bug checking
      channel_to_categorise.print_data_to_csv(OUT_DIR+"categorised_"+this_string_map["fixed_channel_csv_name"]);

      cout << "At present the tag_nodes_from_categorised_channel is still inder testing so I am exiting now" << endl;
      exit(0);
    }

    // Now some logic for looking for upslope influenced pixels, used in the divide finder
    if (this_bool_map["tag_nodes"])
    {
      cout << "I am going to do some pixel tagging!" << endl;
      cout << "My upslope tagging distance is " << this_float_map["upslope_tagging_distance"] << " metres." << endl;
      LSDIndexRaster tagged_raster(DATA_DIR+this_string_map["tagged_raster_input_name"],raster_ext);

      // Now get the uplsope tags
      LSDIndexRaster tr = FlowInfo.tag_upstream_nodes(tagged_raster,this_float_map["upslope_tagging_distance"]);
      string tag_raster_name = OUT_DIR+OUT_ID+"_upslope_tagged_pixels";
      tr.write_raster(tag_raster_name,raster_ext);

      if (this_bool_map["tag_downslope_nodes"])
      {
        cout << "I am tagging downslope pixels." << endl;
        LSDIndexRaster tr_down = FlowInfo.tag_downstream_nodes(tagged_raster,this_float_map["downslope_tagging_distance"]);

        string tag_raster_name = OUT_DIR+OUT_ID+"_downslope_tagged_pixels";
        tr_down.write_raster(tag_raster_name,raster_ext);

      }
    }





    // Now some logic for looking for upslope influenced pixels, used in the divide finder
    if (this_bool_map["divide_finder"])
    {
      // First set up the strips
      LSDRasterMaker RM(filled_topography);

      int NCols = filled_topography.get_NCols();
      int NRows = filled_topography.get_NRows();

      float NDV = filled_topography.get_NoDataValue();
      RM.set_to_constant_value(NDV);

      // E or N strip
      int start_row_or_col, end_row_or_col;
      if(this_bool_map["horizontal_strips"])
      {
        start_row_or_col = 0;
        end_row_or_col = this_int_map["n_row_or_col_for_strips"];

        RM.add_strip(start_row_or_col, end_row_or_col, this_bool_map["horizontal_strips"], 0);
        cout << "Got N strip" << endl;

        start_row_or_col = NRows - 1 -this_int_map["n_row_or_col_for_strips"];
        end_row_or_col = NRows - 1;    

        RM.add_strip(start_row_or_col, end_row_or_col, this_bool_map["horizontal_strips"], 1);
        cout << "Got S strip" << endl;
      }
      else
      {
        start_row_or_col = 0;
        end_row_or_col = this_int_map["n_row_or_col_for_strips"];

        RM.add_strip(start_row_or_col, end_row_or_col, this_bool_map["horizontal_strips"], 0);
        cout << "Got E strip" << endl;

        start_row_or_col = NCols - 1-this_int_map["n_row_or_col_for_strips"];
        end_row_or_col =  NCols - 1;    

        RM.add_strip(start_row_or_col, end_row_or_col, this_bool_map["horizontal_strips"], 1);        
        cout << "Got W strip" << endl;
      }    

      LSDRaster strip = RM.return_as_raster();
      LSDIndexRaster tagged_raster(strip);

      cout << "Printing raster" << endl;
      string tag_raster_name = OUT_DIR+OUT_ID+"_tagged";
      tagged_raster.write_raster(tag_raster_name,raster_ext);

      // Now get the uplsope tags
      LSDIndexRaster tr = FlowInfo.tag_nodes_upstream_of_baselevel(tagged_raster);
      tag_raster_name = OUT_DIR+OUT_ID+"_tagged_basins";
      tr.write_raster(tag_raster_name,raster_ext);

    }


    //=======================================================================================================
    //
    //.#####...#####....####...######..##..##...####....####...######...........####...#####...######...####..
    //.##..##..##..##..##..##....##....###.##..##..##..##......##..............##..##..##..##..##......##..##.
    //.##..##..#####...######....##....##.###..######..##.###..####............######..#####...####....######.
    //.##..##..##..##..##..##....##....##..##..##..##..##..##..##..............##..##..##..##..##......##..##.
    //.#####...##..##..##..##..######..##..##..##..##...####...######..........##..##..##..##..######..##..##.
    //
    //========================================================================================================
    // Now, if you want, calculate drainage areas
    //=================================================================
    if (this_bool_map["print_dinf_drainage_area_raster"])
    {
      cout << "I am writing dinf drainage area to raster." << endl;
      string DA_raster_name = OUT_DIR+OUT_ID+"_dinf_area";
      LSDRaster DA1 = filled_topography.D_inf();
      LSDRaster DA2 = DA1.D_inf_ConvertFlowToArea();
      DA2.write_raster(DA_raster_name,raster_ext);
    }

    if (this_bool_map["print_d8_drainage_area_raster"] ||
        this_bool_map["find_basins"])
    {
      LSDRaster DA_d8 = FlowInfo.write_DrainageArea_to_LSDRaster();

      if (this_bool_map["print_d8_drainage_area_raster"])
      {
        cout << "I am writing d8 drainage area to raster." << endl;
        string DA_raster_name = OUT_DIR+OUT_ID+"_d8_area";
        DA_d8.write_raster(DA_raster_name,raster_ext);
      }
    }

    if (this_bool_map["print_QuinnMD_drainage_area_raster"])
    {
      cout << "I am writing Quinn drainage area to raster." << endl;
      string DA_raster_name = OUT_DIR+OUT_ID+"_QMD_area";
      LSDRaster DA3 = filled_topography.QuinnMDFlow();
      DA3.write_raster(DA_raster_name,raster_ext);
    }

    if (this_bool_map["print_FreemanMD_drainage_area_raster"])
    {
      cout << "I am writing Freeman drainage area to raster." << endl;
      string DA_raster_name = OUT_DIR+OUT_ID+"_FMD_area";
      LSDRaster DA4 = filled_topography.FreemanMDFlow();
      DA4.write_raster(DA_raster_name,raster_ext);
    }
    if (this_bool_map["print_MD_drainage_area_raster"])
    {
      cout << "I am writing mulitdirection drainage area to raster." << endl;
      string DA_raster_name = OUT_DIR+OUT_ID+"_MD_area";
      LSDRaster DA5 = filled_topography.M2DFlow();
      DA5.write_raster(DA_raster_name,raster_ext);
    }

    // Get the distance from outet raster if you want it.
    LSDRaster FD;
    if(this_bool_map["print_distance_from_outlet"] ||
       this_bool_map["find_basins"] ||
       this_bool_map["print_chi_data_maps"] ||
       this_bool_map["extract_single_channel"])
    {
      cout << "I need to calculate the flow distance now." << endl;
      FD = FlowInfo.distance_from_outlet();

      if(this_bool_map["print_distance_from_outlet"])
      {
        cout << "I am writing a distance from outlet raster." << endl;
        string FD_raster_name = OUT_DIR+OUT_ID+"_FDIST";
        FD.write_raster(FD_raster_name,raster_ext);
      }
    }

    // Logic for printing single channel
    if (this_bool_map["extract_single_channel"])
    {
      // Get the latitude and longitude

      LSDRaster DA;
      if (this_bool_map["use_dinf_for_single_channel"])
      {
        LSDRaster DA1 = filled_topography.D_inf();
        DA = DA1.D_inf_ConvertFlowToArea();
      }
      else
      {
        DA = FlowInfo.write_DrainageArea_to_LSDRaster();
      }

      cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
      LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["channel_source_fname"]) );

      // Get the local coordinates
      vector<float> fUTM_easting,fUTM_northing;
      source_points_data.get_x_and_y_from_latlong(fUTM_easting,fUTM_northing);
      float X,Y;

      if ( int(fUTM_easting.size()) != 0)
      {


        X = fUTM_easting[0];
        Y = fUTM_northing[0];

        cout << "Looking for the single source. " << endl;
        cout << "Easting is is: " << X << " and northing is: " << Y << endl;


        vector<int> node_list = FlowInfo.get_flow_path(X,Y);

        cout << "Your node list has " << node_list.size() << "nodes" << endl;

        string fname = "single_channel";
        FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(node_list,DATA_DIR, fname,
                                                      filled_topography, FD, DA);

      }
      else
      {
        cout << "You are trying to make a channel but I cannot find the source." << endl;
      }
    }


    //==============================================================
    //
    //..####...######..#####...######...####...##...##..........##..##..######..######.
    //.##........##....##..##..##......##..##..###.###..........###.##..##........##...
    //..####.....##....#####...####....######..##.#.##..........##.###..####......##...
    //.....##....##....##..##..##......##..##..##...##..........##..##..##........##...
    //..####.....##....##..##..######..##..##..##...##..........##..##..######....##...
    //
    //==============================================================    
    // This is the logic for a simple stream network
    if (this_bool_map["print_channels_to_csv"]
        || this_bool_map["print_junctions_to_csv"]
        || this_bool_map["print_sources_to_csv"]
        || this_bool_map["find_basins"]
        || this_bool_map["print_channel_tips_raster"]
        || this_bool_map["print_chi_data_maps"]
        || this_bool_map["print_junction_angles_to_csv"]
        || this_bool_map["extract_ridges"]
        || this_bool_map["calculate_hypsometric_integral"]
        || this_bool_map["clip_raster_to_basins"])
    {
      // calculate the flow accumulation
      cout << "\t Calculating flow accumulation (in pixels)..." << endl;
      LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      //==============================================================
      //
      //..####....####...##..##..#####....####...######...####..
      //.##......##..##..##..##..##..##..##..##..##......##.....
      //..####...##..##..##..##..#####...##......####.....####..
      //.....##..##..##..##..##..##..##..##..##..##..........##.
      //..####....####....####...##..##...####...######...####..
      //
      //==============================================================
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
      }

      // now get the junction network
      LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);


      // Extract all the ridges
      if(this_bool_map["extract_ridges"])
      {
        string RN_csv_name = OUT_DIR+OUT_ID+"_Ridges.csv";
        map<int, int > RidgeNetwork_nodes = JunctionNetwork.ExtractAllRidges(FlowInfo,RN_csv_name);

        // convert to geojson if that is what the user wants
        // It is read more easily by GIS software but has bigger file size
        if ( this_bool_map["convert_csv_to_geojson"])
        {
          cout << "Let me convert that data to json." << endl;
          string gjson_name = OUT_DIR+OUT_ID+"_Ridges.geojson";
          LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_Ridges.csv");
          thiscsv.print_data_to_geojson(gjson_name);
        }
      }


      if( this_bool_map["print_channel_tips_raster"])
      {
        LSDIndexRaster channel_tips = FlowInfo.get_channel_tip_raster(sources, int_default_map["channel_tips_pixel_distance"]);
        string ct_name = OUT_DIR+OUT_ID+"_chan_tips";
        channel_tips.write_raster(ct_name,"bil");
      }

      // Print channels and junctions if you want them.
      if( this_bool_map["print_channels_to_csv"])
      {
        cout << "I am going to print the channel network." << endl;
        string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
        if ( this_bool_map["use_extended_channel_data"])
        {
          cout << "I am going to use the extended channel network data outputs." << endl;
          JunctionNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, filled_topography);
        }
        if ( this_bool_map["print_channel_data_plus_surface_metrics"])
        {
          cout << "I am going to print the channel network with some surrounding surface metrics (slope, local relief, and curvature." << endl;
          int kernel_type = 1;
          channel_csv_name = OUT_DIR+OUT_ID+"_CN_surface_metrics";
          LSDRaster Relief = filled_topography.calculate_relief(this_float_map["surface_fitting_radius"], kernel_type);
          JunctionNetwork.PrintChannelNetworkToCSV_WithSurfaceMetrics(FlowInfo, channel_csv_name, filled_topography, surface_fitting[1], Relief, surface_fitting[3]);
        }
        else
        {
          JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
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

      // print junctions
      if( this_bool_map["print_junctions_to_csv"])
      {
        cout << "I am writing the junctions to csv." << endl;
        string channel_csv_name = OUT_DIR+OUT_ID+"_JN.csv";
        JunctionNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

        if ( this_bool_map["convert_csv_to_geojson"])
        {
          string gjson_name = OUT_DIR+OUT_ID+"_JN.geojson";
          LSDSpatialCSVReader thiscsv(channel_csv_name);
          thiscsv.print_data_to_geojson(gjson_name);
        }
      }

      // print the junction angles
      if( this_bool_map["print_junction_angles_to_csv"])
      {
        // Calculate flow distance
        LSDRaster FlowDistance = FlowInfo.distance_from_outlet();

        cout << "I am testing the junction angle code with elevations." << endl;
        string JAngles_csv_name = OUT_DIR+OUT_ID+"_FULL_JAngles.csv";
        vector<int> JunctionList;
        JunctionNetwork.print_complete_junction_angles_to_csv(JunctionList, FlowInfo, filled_topography, FlowDistance, this_float_map["SA_vertical_interval"], JAngles_csv_name);

        if ( this_bool_map["convert_csv_to_geojson"])
        {
          string gjson_name = OUT_DIR+OUT_ID+"_FULL_JAngles.geojson";
          LSDSpatialCSVReader thiscsv(JAngles_csv_name);
          thiscsv.print_data_to_geojson(gjson_name);
        }
      }


      // Print sources
      if( this_bool_map["print_sources_to_csv"])
      {
        string sources_csv_name = OUT_DIR+OUT_ID+"_ATsources.csv";

        //write channel_heads to a csv file
        FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(sources, sources_csv_name);
        string sources_csv_name_2 = OUT_DIR+OUT_ID+"_ATsources_rowcol.csv";
        FlowInfo.print_vector_of_nodeindices_to_csv_file(sources, sources_csv_name_2);

        if ( this_bool_map["convert_csv_to_geojson"])
        {
          string gjson_name = OUT_DIR+OUT_ID+"_ATsources.geojson";
          LSDSpatialCSVReader thiscsv(sources_csv_name);
          thiscsv.print_data_to_geojson(gjson_name);
        }
      }   // End print sources logic


      //=================================================================================
      //.#####....####....####...######..##..##...####..
      //.##..##..##..##..##........##....###.##..##.....
      //.#####...######...####.....##....##.###...####..
      //.##..##..##..##......##....##....##..##......##.
      //.#####...##..##...####...######..##..##...####..
      //................................................ 
      //=================================================================================     
      // Now we check if we are going to deal with basins
      if(this_bool_map["find_basins"] ||
         this_bool_map["print_chi_data_maps"] ||
         this_bool_map["calculate_basin_statistics"] ||
         this_bool_map["clip_raster_to_basins"])
      {
        cout << "I am now going to extract some basins for you." << endl;
        vector<int> BaseLevelJunctions;
        vector<int> BaseLevelJunctions_Initial;

        // deal with the baselevel junctions file
        string BaselevelJunctions_file;
        string test_BaselevelJunctions_file = LSDPP.get_BaselevelJunctions_file();
        if(this_string_map["BaselevelJunctions_file"] == "NULL" && test_BaselevelJunctions_file == "NULL")
        {
          cout << "No baselevel junctions file found. I am going to use algorithms to extract basins." << endl;
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

        // Now we try to get the basins using
        //Check to see if a list of junctions for extraction exists
        if (BaselevelJunctions_file == "NULL" || BaselevelJunctions_file == "Null" || BaselevelJunctions_file == "null" || BaselevelJunctions_file.empty() == true)
        {
          if(this_bool_map["get_basins_from_outlets"])
          {
            cout << "I am going to get basins lat-long coordinates" << endl;
            string full_BL_LL_name = DATA_DIR+this_string_map["basin_outlet_csv"];
            cout << "The file is: " << full_BL_LL_name << endl;
            int search_radius_nodes = this_int_map["search_radius_nodes"];
            int threshold_stream_order = 3;
            BaseLevelJunctions = JunctionNetwork.snap_point_locations_to_upstream_junctions_from_latlong_csv(full_BL_LL_name,
                                                               search_radius_nodes, threshold_stream_order,FlowInfo, RI);
          }
          else
          {
            if(this_bool_map["prune_edge_draining_basins"])
            {
              cout << "I am going to select basins for you using an algorithm. " << endl;
              cout << "I am going to look for basins in a pixel window that are not influended by nodata." << endl;
              cout << "I am also going to remove any nested basins." << endl;
              cout << "The pixel limits are: lower: " << this_int_map["minimum_basin_size_pixels"] << " and upper: " << this_int_map["maximum_basin_size_pixels"] << endl;
              BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(FlowInfo, filled_topography, FlowAcc,
                                                      this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
            }
            else
            {
              cout << endl << endl;
              cout << "======================================" << endl;
              cout << "I am going to select basins for you using an algorithm. " << endl;
              cout << "I am going to look for basins in a pixel window." << endl;
              cout << "I am also going to remove any nested basins." << endl;
              cout << "!!!!WARNING!!!! you are selecting basins" << endl;
              cout << "but you have chosen not to remove basins with pixels at the edge." << endl;
              cout << "This could mean your drainage area or chi coordinate will be inaccurate." << endl;
              cout << "You should only do this if you are very confident you either" << endl;
              cout << " don't need drainge area or chi, or if you know the " << endl;
              cout << "effect will be very small." << endl;
              cout << "The pixel limits are: lower: " << this_int_map["minimum_basin_size_pixels"] << " and upper: " << this_int_map["maximum_basin_size_pixels"] << endl;
              cout << "======================================" << endl << endl << endl;   

              BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_Keep_Edge(FlowInfo, filled_topography, FlowAcc,
                                                      this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
            }
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

          if (this_bool_map["prune_edge_draining_basins"])
          {
            // Now make sure none of the basins drain to the edge
            cout << "I am pruning junctions that are influenced by the edge of the DEM!" << endl;
            cout << "This is necessary because basins draining to the edge will have incorrect chi values." << endl;
            BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions_Initial, FlowInfo, filled_topography);
          }
          else
          {
            cout << endl << endl;
            cout << "======================================" << endl;
            cout << "WARNING you are selecting basins" << endl;
            cout << "but you have chosen not to remove basins with pixels at the edge." << endl;
            cout << "This could mean your drainage area or chi coordinate will be inaccurate." << endl;
            cout << "You should only do this if you are very confident you either" << endl;
            cout << " don't need drainge area or chi, or if you know the " << endl;
            cout << "effect will be very small." << endl;
            cout << "======================================" << endl << endl << endl;            
          }

          cout << "The remaining baselevel junctions are: " << endl;
          for (int i = 0; i<int(BaseLevelJunctions.size()); i++)
          {
            cout << BaseLevelJunctions[i] << endl;
          }

        }    // end logic for reading from junctions list

        // Now check for largest basin, if that is what you want.
        if (this_bool_map["only_take_largest_basin"])
        {
          cout << endl << endl << "==========================" << endl;
          cout << "I am only going to take the largest basin." << endl;
          cout << "WARNING: the basins have already been selected within the pixel window" << endl;
          cout << "If you want to get largest basin in the DEM you need to " << endl;
          cout << "Set minimum_basin_size_pixels to a very large number." << endl;
          cout << "============================" << endl << endl;
          BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Largest(BaseLevelJunctions, FlowInfo, FlowAcc);
        }

        // Correct number of base level junctions
        int N_BaseLevelJuncs = BaseLevelJunctions.size();
        cout << "The number of basins I will analyse is: " << N_BaseLevelJuncs << endl;
        if (N_BaseLevelJuncs == 0)
        {
          cout << "I am stopping here since I don't have any basins to analyse." << endl;
          exit(EXIT_FAILURE);
        }
        else
        {
          cout << "The baselevel junction numbers are: " << endl;
          for(int i = 0; i< N_BaseLevelJuncs; i++)
          {
            cout << BaseLevelJunctions[i] << ", ";
          }
          cout << endl;
        }


        // Now for junction angles only within the selected basins
        if(this_bool_map["print_junction_angles_to_csv_in_basins"])
        {
          vector<int> basin_junctions;
          for(int i = 0; i< N_BaseLevelJuncs; i++)
          {
            // get the upslope junctions
            vector<int> these_junctions = JunctionNetwork.get_upslope_junctions(BaseLevelJunctions[i]);

            // now append these to the master list
            basin_junctions.insert( basin_junctions.end(), these_junctions.begin(), these_junctions.end() );
          }

          int N_basin_junctions = int(basin_junctions.size());
          if (N_basin_junctions == 0)
          {
            cout << "I am stopping here since I don't have any junctions over which to measure the junction angle." << endl;
            exit(EXIT_FAILURE);
          }

          // now get all the junction angles
          // Calculate flow distance
          LSDRaster FlowDistance = FlowInfo.distance_from_outlet();

          cout << "I am getting junction angles in your basins." << endl;
          string JAngles_csv_name = OUT_DIR+OUT_ID+"_FULL_Basins_JAngles.csv";
          JunctionNetwork.print_complete_junction_angles_to_csv(basin_junctions, FlowInfo, filled_topography, 
                                                                FlowDistance, this_float_map["SA_vertical_interval"], 
                                                                JAngles_csv_name);

          if ( this_bool_map["convert_csv_to_geojson"])
          {
            string gjson_name = OUT_DIR+OUT_ID+"_FULL_Basins_JAngles.geojson";
            LSDSpatialCSVReader thiscsv(JAngles_csv_name);
            thiscsv.print_data_to_geojson(gjson_name);
          }          
        }



        // Now to clip basins
        if(this_bool_map["clip_raster_to_basins"])
        {
          cout << "I will now clip the raster to your basins" << endl;
          vector<int> outlet_nodes = JunctionNetwork.get_node_list_from_junction_list(BaseLevelJunctions);
          float masking_threshold = 0.5;
          bool belowthresholdisnodata = true;

          LSDRaster BasinMask = FlowInfo.get_upslope_node_mask(outlet_nodes);
          LSDRaster Masked_raster = filled_topography.mask_to_nodata_using_threshold_using_other_raster_expunge_nodata_in_mask(masking_threshold,belowthresholdisnodata, BasinMask);

          cout << "I've got only the masked raster, now I am going to reduce the size of the DEM " << endl;
          cout << "So it only includes the basin you selected." << endl;
          LSDRaster Masked_and_clipped_raster = Masked_raster.RasterTrimmer();   

          string this_raster_name = OUT_DIR+OUT_ID+"_BASINCLIP";
          Masked_and_clipped_raster.write_raster(this_raster_name,raster_ext);         
        }


        // Now for the logic for basin statistics.
        if(this_bool_map["calculate_basin_statistics"])
        {
          cout << "Hello there, I am getting some basin statistics." << endl;

          string basin_details_fname = OUT_DIR+OUT_ID+"_basin_details.csv";
          ofstream bd_out;
          bd_out.open(basin_details_fname.c_str());
          bd_out << "id,junction_number";

          if(this_bool_map["print_slope"])
          {
            bd_out << ",slope_max,slope_84,slope_median,slope_16,slope_min";
          }
          if(this_bool_map["print_curvature"])
          {
            bd_out << ",curv_max,curv_84,curv_median,curv_16,curv_min";
          }
          bd_out << endl;


          int N_BL = int(BaseLevelJunctions.size());
          for(int bsn = 0; bsn < N_BL; bsn++)
          {
            cout << "Processing basin " << bsn << " of " << N_BL << endl;
            LSDBasin ThisBasin(BaseLevelJunctions[bsn], FlowInfo, JunctionNetwork);

            bd_out << bsn << "," << BaseLevelJunctions[bsn];

            if(this_bool_map["print_slope"])
            {
              float slope_max,slope_84,slope_median,slope_16,slope_min;
              slope_max = ThisBasin.CalculateBasinMax(FlowInfo, surface_fitting[1]);
              slope_84 = ThisBasin.CalculateBasinPercentile(FlowInfo, surface_fitting[1], 84);
              slope_median =  ThisBasin.CalculateBasinMedian(FlowInfo, surface_fitting[1]);
              slope_16 = ThisBasin.CalculateBasinPercentile(FlowInfo, surface_fitting[1], 16);
              slope_min = ThisBasin.CalculateBasinMin(FlowInfo, surface_fitting[1]);
              bd_out << "," << slope_max << "," << slope_84 << "," << slope_median << "," << slope_16 << "," << slope_min;
            }

            if(this_bool_map["print_curvature"])
            {
              float curv_max,curv_84,curv_median,curv_16,curv_min;
              curv_max = ThisBasin.CalculateBasinMax(FlowInfo, surface_fitting[3]);
              curv_84 = ThisBasin.CalculateBasinPercentile(FlowInfo, surface_fitting[3], 84);
              curv_median =  ThisBasin.CalculateBasinMedian(FlowInfo, surface_fitting[3]);
              curv_16 = ThisBasin.CalculateBasinPercentile(FlowInfo, surface_fitting[3], 16);
              curv_min = ThisBasin.CalculateBasinMin(FlowInfo, surface_fitting[3]);
              bd_out << "," << curv_max << "," << curv_84 << "," << curv_median << "," << curv_16 << "," << curv_min;
            }
            bd_out << endl;

          }
        }

        //====================================================================================
        //
        //..####...##..##..######..........######..##..##..######..#####....####....####...######..######...####...##..##.
        //.##..##..##..##....##............##.......####.....##....##..##..##..##..##..##....##......##....##..##..###.##.
        //.##......######....##............####......##......##....#####...######..##........##......##....##..##..##.###.
        //.##..##..##..##....##............##.......####.....##....##..##..##..##..##..##....##......##....##..##..##..##.
        //..####...##..##..######..........######..##..##....##....##..##..##..##...####.....##....######...####...##..##.
        //
        //====================================================================================
        if ( this_bool_map["print_chi_data_maps"])
        {
          LSDRaster chi_coordinate;       
          cout << "I am calculating the chi coordinate." << endl;
          chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(this_float_map["m_over_n"],this_float_map["A_0"],this_int_map["threshold_contributing_pixels"]);

          // Now we get the channel segments. This information is used for plotting chi stuff
          vector<int> source_nodes;
          vector<int> outlet_nodes;
          vector<int> baselevel_node_of_each_basin;
          int n_nodes_to_visit = 10;
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
            cout << "Note if you have nested basins the basin numbers only take the larger basin." << endl;
            cout << "=====================================================" << endl << endl;

            JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, FD,
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
            cout << "Note if you have nested basins the basin numbers only take the larger basin." << endl;
            cout << "=====================================================" << endl << endl;

            JunctionNetwork.get_overlapping_channels(FlowInfo, BaseLevelJunctions, FD,
                                          source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);
          }

          cout << "I've got the overlapping channels. The baselevel junctions are now: " << endl;
          for (int i = 0; i<int(BaseLevelJunctions.size()); i++)
          {
            cout << BaseLevelJunctions[i] << endl;
          }
          cout << "And the baselevel nodes of each basin are: " << endl;
          for (int i = 0; i<int(BaseLevelJunctions.size()); i++)
          {
            cout << baselevel_node_of_each_basin[i] << endl;
          }


          //======================================================================
          // Print a basin raster if you want it.
          if(this_bool_map["print_basin_raster"])
          {
            cout << "I am going to print the basins for you. " << endl;
            cout << "This is using the full chi data" << endl;
            LSDChiTools ChiTool_basins(FlowInfo);
            LSDRaster DA_for_chi = FlowInfo.write_DrainageArea_to_LSDRaster();
            ChiTool_basins.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                                    filled_topography, FD,
                                    DA_for_chi , chi_coordinate);
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
            LSDChiTools ChiTool_chi_checker(FlowInfo);
            LSDRaster DA_for_chi = FlowInfo.write_DrainageArea_to_LSDRaster();

            ChiTool_chi_checker.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                                    filled_topography, FD,
                                    DA_for_chi, chi_coordinate);


            string chi_data_maps_string = OUT_DIR+OUT_ID+"_chi_data_map.csv";

            if(this_bool_map["use_extended_channel_data"])
            {
              ChiTool_chi_checker.print_chi_data_map_to_csv_with_junction_information(FlowInfo, JunctionNetwork,chi_data_maps_string);  
            }
            else
            {
              ChiTool_chi_checker.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);  
            }


            if ( this_bool_map["convert_csv_to_geojson"])
            {
              string gjson_name = OUT_DIR+OUT_ID+"_chi_data_map.geojson";
              LSDSpatialCSVReader thiscsv(chi_data_maps_string);
              thiscsv.print_data_to_geojson(gjson_name);
            }
          }
        }
        else
        {
          if (this_bool_map["print_basin_raster"])
          {
            cout << "I am going to print the basins for you. " << endl;
            cout << "I am not printing the chi data maps." << endl;
            LSDChiTools ChiTool_basins(FlowInfo);
            string basin_raster_prefix = OUT_DIR+OUT_ID;
            ChiTool_basins.print_basins(FlowInfo, JunctionNetwork, BaseLevelJunctions, basin_raster_prefix);          
          }
        }
      }   // end logic for basin finding
      cout << "Finished with basins" << endl;


      // Calculate hypsometric integral for the channel network
      if( this_bool_map["calculate_hypsometric_integral"])
      {
          cout << "Calculating the hypsometric integral..." << endl;
          string HI_csv_name = OUT_DIR+OUT_ID+"_HI";
          //JunctionNetwork.calculate_hypsometric_integral(HI_csv_name, FlowInfo, filled_topography);
          JunctionNetwork.calculate_upstream_elevations_network(HI_csv_name, FlowInfo, filled_topography, this_float_map["HI_bin_width"], this_float_map["HI_lower_limit"], this_float_map["HI_upper_limit"]);
      }


      //================================================================================================
      //..####....####...##..##..##..##..######...####...######..######..##..##..######..######..##..##.
      //.##..##..##..##..###.##..###.##..##......##..##....##......##....##..##....##......##.....####..
      //.##......##..##..##.###..##.###..####....##........##......##....##..##....##......##......##...
      //.##..##..##..##..##..##..##..##..##......##..##....##......##.....####.....##......##......##...
      //..####....####...##..##..##..##..######...####.....##....######....##....######....##......##...
      //................................................................................................
      //================================================================================================
      if (this_bool_map["calculate_connectivty_index"] )
      {
        cout << "Hello there, I am going to calculate connectivity" << endl;
        cout << "I need some rasters and various data elements" << endl;

        // Okay so here is what you need to do
        // 1) Check if the dinf raster exists
        // 2) Check if a channel file exists
        // 3) Then get the "roughness" raster (which is just the standard deviation)
        // 4) Get the average upslope roughness and the average upslope gradient
        // 5) Then for each pixel you need to route to the nearest channel accumulating roughness and slope and flow distance
        // 6) You will then have everything you need to calculate CI, although you won't have eliminated areas upstream of sinks

        // Part 1: get the dinfinity raster. 
        // Note that at the moment this is always recalculated
        string Dinf_fname = OUT_DIR+OUT_ID+"_dinf_area";  
        LSDRaster Dinf_raster(Dinf_fname,"bil");
 

      }   // end logic for connectivity



    }     // end logic for tasks related to channel network extraction
    cout << "Done with channel extraction" << endl;

  }       // end logic for tasks requiring flow info and filling
  cout << "I'm all finished! Have a nice day." << endl;
}
