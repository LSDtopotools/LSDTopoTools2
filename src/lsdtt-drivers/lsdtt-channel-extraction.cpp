//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// lsdtt-channel-extraction 
//
// This programme extractes channel metrics using a variety of methods
//
// The different methods are described in this paper:
// Clubb, F.J., Mudd, S.M., Milodowski, D.T., Hurst, M.D., Slater, L.J., 2014. 
// Objective extraction of channel heads from high-resolution topographic data. 
// Water Resources Research 50, 4283–4304. https://doi.org/10.1002/2013WR015167
//
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate
// different kinds of chi analysis
//
// Call with -h to generate a help file
//
// Copyright (C) 2022 Fiona J. Clubb Simon M. Mudd 2022
//
// Developed by:
//  Fiona Clubb
//  Simon M. Mudd
//  David T. Milodowski
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona J. Clubb, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>
#include "../LSDStatsTools.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDParameterParser.hpp"

int main (int nNumberofArgs,char *argv[])
{
  //start the clock
  clock_t begin = clock();

  string version_number = "0.9";
  string citation = "http://doi.org/10.5281/zenodo.4577879";

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the LSDTopoTools channel extraction tool!||" << endl;
  cout << "|| This program has a number of options to extract     ||" << endl;
  cout << "|| channel networks.                                   ||" << endl;
  cout << "|| This program was developed by Fiona J. Clubb        ||" << endl;
  cout << "||  and Simon M. Mudd                                  ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;
  cout << "=========================================================" << endl;  
  cout << "|| Citation for this code is:                          ||" << endl;
  cout << "|| " << citation << endl;
  cout << "|| If you use these routines please cite:              ||" << endl;   
  cout << "|| https://www.doi.org/10.1002/2013WR015167            ||" << endl;
  cout << "|| If you use the Wiener routine please cite:          ||" << endl;   
  cout << "|| https://www.doi.org/10.5194/esurf-4-627-2016        ||" << endl;  
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
    ofs.open("./lsdtt-channel-extraction-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    cout << endl << endl << endl << "==============================================" << endl;    
    cout << "This is lsdtt-channel-extraction version number " << version_number << endl;
    cout << "If the version contains a 'd' then you are using a development version." << endl;
    cout << "=========================================================" << endl;
    ofstream ofs;
    ofs.open("./lsdtt-channel-extraction-version.txt");
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
 
  string_default_map["CHeads_file"] = "NULL";
  help_map["CHeads_file"] = {  "string","NULL","The name of a channel heads file.","You can output this csv file with the channel extraction algorithms. It contains latitude and longitude values of the channel heads."};
 

  bool_default_map["use_spiral_trimmer_for_edge_influence"] = false;
  help_map["use_spiral_trimmer_for_edge_influence"] = {  "bool","false","Makes sure raster is rectangular before running drainage extraction.","Use this if you want to avoid edge effects of your flow routing that can cause seg faults."};



  // set default float parameters
  int_default_map["threshold_contributing_pixels"] = 1000;
  help_map["threshold_contributing_pixels"] = {  "int","1000","The number of contributing pixels needed to start a channel using the threshold method.","This is in pixels not drainage area. More options are in the lsdtt-channel-extraction tool."};

  int_default_map["connected_components_threshold"] = 100;
  help_map["connected_components_threshold"] = {  "int","100","Number of connected pixels to classify as a channel.","Used in Pelletier and Wiener methods."};

  float_default_map["surface_fitting_radius"] = 6;
  help_map["surface_fitting_radius"] = {  "float","6","Our surface fitting routines fit a polynomial over the points with in a radius defined by surface_fitting_radius and then differentiate this surface to get the surface metrics like gradient and curvature","If not bigger than the pixel_size*sqrt(2) then will increase to that number. Unlike lsdtt-basic-metrics this default is tuned for lidar data"};

  float_default_map["pruning_drainage_area"] = 1000;
  help_map["pruning_drainage_area"] = {  "float","1000","In both driech and Wiener channels with less than this drainage area are removed during the pruning process.","In m^2."};

  float_default_map["curvature_threshold"] = 0.1;
  help_map["curvature_threshold"] = {  "float","0.1","For the Pelletier method this is the threshold of planform curvature for a channel head.","In 1/m"};

  float_default_map["minimum_drainage_area"] = 400;
  help_map["minimum_drainage_area"] = {  "float","400","For the Pelletier method this is the minimum drainage area for a channel head.","In m^2"};

  int_default_map["number_of_junctions_dreich"] = 1;
  help_map["number_of_junctions_dreich"] = {  "int","1","For driech the number of junctions to pass to construct the chi profile at the upstream end of the channel.","Applies only to dreich"};

    

  float_default_map["A_0"] = 1.0;
  help_map["A_0"] = {  "float","1.0","The A_0 parameter for chi computation. See https://doi.org/10.1002/esp.3302","Usually set to 1 so that the slope in chi-elevation space is the same as k_sn"};
   
  float_default_map["m_over_n"] = 0.5;
  help_map["m_over_n"] = {  "float","0.5","The concavity index for chi calculations. Usually denoted as the Greek symbol theta.","Default is 0.5 but possibly 0.45 is better as Kwang and Parker suggest 0.5 leads to unrealistic behaviour in landscape evolution models."};


  // The actual channel extraction routines
  bool_default_map["print_area_threshold_channels"] = true;
  help_map["print_area_threshold_channels"] = {  "bool","true","Prints the channel network determined by an area threshold.","Output is a csv file."};

  bool_default_map["print_dreich_channels"] = false;
  help_map["print_dreich_channels"] = {  "bool","false","Prints the channel network determined by the driech method.","Output is a csv file."};

  bool_default_map["print_pelletier_channels"] = false;
  help_map["print_pelletier_channels"] = {  "bool","false","Prints the channel network determined by the Pelletier method.","Output is a csv file."};

  bool_default_map["print_wiener_channels"] = false;
  help_map["print_wiener_channels"] = {  "bool","false","Prints the channel network determined by the Wiener method which is a mashup of Passalacqua and Pelletier methods first reported by grieve et al 2016.","Output is a csv file."};

  bool_default_map["print_perona_malik_channels"] = false;
  help_map["print_perona_malik_channels"] = {  "bool","false","Prints the channel network using the GeoNet method proposed by Passalacqua et al. 2010.","Output is a csv file."};

  bool_default_map["print_stream_order_raster"] = false;
  help_map["print_stream_order_raster"] = {  "bool","false","Prints a raster with _SO in filename with stream orders of channel in the appropriate pixel.","Generates a big file so we suggest printing the network to csv."};

  bool_default_map["print_sources_to_raster"] = false;
  help_map["print_sources_to_raster"] = {  "bool","false","If true the sources to a raster.","Inefficient so print_sources_to_csv is preferred."};
  
  bool_default_map["print_wiener_filtered_raster"] = false;
  help_map["print_wiener_filtered_raster"] = {  "bool","false","Runs the raster through a Wiener filter and prints the result.","Output file has _Wfilt in filename."};


  // some simple raster printing
  bool_default_map["write_hillshade"] = false;
  help_map["write_hillshade"] = {  "bool","false","Write the hillshade raster.","You need this for a lot of our plotting routines. Filename includes _HS"};

  bool_default_map["print_fill_raster"] = false;
  help_map["print_fill_raster"] = {  "bool","false","Prints the fill raster.","Filename includes _FILL"};

  // This converts all csv files to geojson (for easier loading in a GIS)
  bool_default_map["convert_csv_to_geojson"] = false;
  help_map["convert_csv_to_geojson"] = {  "bool","false","Converts csv files to geojson files","Makes csv output easier to read with a GIS. Warning: these files are much bigger than csv files."};

  bool_default_map["print_curvature"]= false;
  help_map["print_curvature"] = {  "bool","false","Prints curvature raster after polynomial fitting.","Part of surface fitting metrics."};

  bool_default_map["print_sources_to_csv"] = false;
  help_map["print_sources_to_csv"] = {  "bool","false","Prints the sources to a csv file.","Each source on its own row with latitude and longitude columns."};

  bool_default_map["print_channels_to_csv"] = true;
  help_map["print_channels_to_csv"] = {  "bool","false","Prints the channel network to a csv file.","This version produces smaller files than the raster version."};

  bool_default_map["use_extended_channel_data"] = false;
  help_map["use_extended_channel_data"] = {  "bool","false","If this is true you get more data columns in your channel network csv.","I will tell you what these columns are one day."};

  bool_default_map["print_channels_to_line_vector"] = true;
  help_map["print_channels_to_line_vector"] = { "bool","true","If this is true print a geojson with the channel network as a polyline.","Used for more efficient visualisation of the channel network."};

  bool_default_map["print_junctions_to_csv"] = false;
  help_map["print_junctions_to_csv"] = {  "bool","false","Prints a csv with the locations and numbers of the junctions.","This is better to use than the raster version."};
 
  // drainage area printing
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
  //LSDPP.print_parameters();

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

  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./lsdtt-channel-extraction-README.csv" << endl;
    string help_prefix = "lsdtt-channel-extraction-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
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
      cout << "I'm replacing your raster with a raster without seas." << endl;
      string this_raster_name = OUT_DIR+OUT_ID;
      topography_raster.write_raster(this_raster_name,raster_ext);
    }
  }
  else
  {
    LSDRaster start_raster((DATA_DIR+DEM_ID), raster_ext);
    topography_raster = start_raster;
  }
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;



  //============================================================================
  // Start gathering necessary rasters
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

  // This routine makes sure that the outer edges of the DEM
  // are trimmed to a rectangle and then the edge nodes are replaced
  // with nodata
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

  // check to see if you need hillshade
  if (this_bool_map["write_hillshade"])
  {
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+OUT_ID+"_hs";
    hs_raster.write_raster(hs_fname,raster_ext);
  }

  // get the flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);


  //=======================================================================================================
  //
  //.#####...#####....####...######..##..##...####....####...######...........####...#####...######...####..
  //.##..##..##..##..##..##....##....###.##..##..##..##......##..............##..##..##..##..##......##..##.
  //.##..##..#####...######....##....##.###..######..##.###..####............######..#####...####....######.
  //.##..##..##..##..##..##....##....##..##..##..##..##..##..##..............##..##..##..##..##......##..##.
  //.#####...##..##..##..##..######..##..##..##..##...####...######..........##..##..##..##..######..##..##.
  //
  //========================================================================================================
  if (this_bool_map["print_dinf_drainage_area_raster"])
  {
    cout << "I am writing dinf drainage area to raster." << endl;
    string DA_raster_name = OUT_DIR+OUT_ID+"_dinf_area";
    LSDRaster DA1 = filled_topography.D_inf();
    LSDRaster DA2 = DA1.D_inf_ConvertFlowToArea();
    DA2.write_raster(DA_raster_name,raster_ext);
  }

  if (this_bool_map["print_d8_drainage_area_raster"])
  {
    string DA_raster_name = OUT_DIR+OUT_ID+"_d8_area";
    LSDRaster DA2 = FlowInfo.write_DrainageArea_to_LSDRaster();
    //cout << "d8:" << endl <<  DA2.get_data_element(452,364) << " " << DA2.get_data_element(1452,762) << endl;
    DA2.write_raster(DA_raster_name,raster_ext);
  }

  if (this_bool_map["print_QuinnMD_drainage_area_raster"])
  {
    string DA_raster_name = OUT_DIR+OUT_ID+"_QMD_area";
    LSDRaster DA3 = filled_topography.QuinnMDFlow();
    //cout << "Qiunn:" << endl <<  DA3.get_data_element(452,364) << " " << DA3.get_data_element(1452,762) << endl;
    DA3.write_raster(DA_raster_name,raster_ext);
  }

  if (this_bool_map["print_FreemanMD_drainage_area_raster"])
  {
    string DA_raster_name = OUT_DIR+OUT_ID+"_FMD_area";
    LSDRaster DA4 = filled_topography.FreemanMDFlow();
    //cout << "Freeman:" << endl <<  DA4.get_data_element(452,364) << " " << DA4.get_data_element(1452,762) << endl;
    DA4.write_raster(DA_raster_name,raster_ext);
  }

  if (this_bool_map["print_MD_drainage_area_raster"])
  {
    string DA_raster_name = OUT_DIR+OUT_ID+"_MD_area";
    LSDRaster DA5 = filled_topography.M2DFlow();
    DA5.write_raster(DA_raster_name,raster_ext);
  }

  //=================================================================
  //==============================================================
  //
  //..####....####...##..##..#####....####...######...####..
  //.##......##..##..##..##..##..##..##..##..##......##.....
  //..####...##..##..##..##..#####...##......####.....####..
  //.....##..##..##..##..##..##..##..##..##..##..........##.
  //..####....####....####...##..##...####...######...####..
  //
  //==============================================================
  //=================================================================
  cout << endl << endl << "The channel heads file is " << CHeads_file << endl;
  if (CHeads_file != "NULL" && CHeads_file != "Null" && CHeads_file != "null")
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    cout << "The channel heads file *MUST* have a csv extension or this will crash!!" << endl;
    LSDRasterInfo ThisRI(filled_topography);
    string csv_filename = DATA_DIR+CHeads_file;
    LSDSpatialCSVReader CHeadCSV(ThisRI,csv_filename);    
    vector<int> sources = CHeadCSV.get_nodeindices_from_lat_long(FlowInfo);
    cout << "\t Got sources!" << endl;

    // now get the junction network
    LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_FromCHF_SO";

      //write stream order array to a raster
      LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }

    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_FromCHF_CN";
      if ( this_bool_map["use_extended_channel_data"])
      {
        cout << "I am going to use the extended channel network data outputs." << endl;
        ChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, filled_topography);
      }
      else
      {
        ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
      }

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_FromCHF_CN.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_FromCHF_CN.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }
  }



  //===============================================================
  // AREA THRESHOLD
  //===============================================================
  if (this_bool_map["print_area_threshold_channels"])
  {
    cout << "I am calculating channels using an area threshold." << endl;
    cout << "Only use this if you aren't that bothered about where the channel heads actually are!" << endl;

    // get some relevant rasters
    LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
    LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    //get the sources: note: this is only to select basins!
    vector<int> sources;
    sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["threshold_contributing_pixels"]);

    // now get the junction network
    LSDJunctionNetwork ChanNetwork(sources, FlowInfo);


    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
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

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_ATsources";

      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(sources);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_AT_SO";

      //write stream order array to a raster
      LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }

    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_AT_CN";
      if ( this_bool_map["use_extended_channel_data"])
      {
        cout << "I am going to use the extended channel network data outputs." << endl;
        ChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, filled_topography);
      }
      else
      {
        ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
      }

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_AT_CN.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_AT_CN.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    // print junctions
    if( this_bool_map["print_junctions_to_csv"])
    {
      cout << "I am writing the junctions to csv." << endl;
      string channel_csv_name = OUT_DIR+OUT_ID+"_AT_JN.csv";
      ChanNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_AT_JN.geojson";
        LSDSpatialCSVReader thiscsv(channel_csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    if (this_bool_map["print_channels_to_line_vector"])
    {
      string channel_geojson_name = OUT_DIR+OUT_ID+"_AT_CN_line";
      cout << "I am going to print your channel network to a line geojson" << endl;
      LSDRaster DA = FlowInfo.write_DrainageArea_to_LSDRaster();
      ChanNetwork.PrintChannelNetworkToLineVector(FlowInfo, channel_geojson_name, filled_topography, DA, DistanceFromOutlet);
    }
  }

  //===============================================================
  // DREICH
  //===============================================================
  if (this_bool_map["print_dreich_channels"])
  {
    cout << "I am calculating channels using the dreich algorighm (DOI: 10.1002/2013WR015167)." << endl;

    // get some relevant rasters
    LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

    LSDRasterSpectral raster(filled_topography);
    string QQ_fname = OUT_DIR+OUT_ID+"__qq.txt";

    cout << "I am am getting the connected components using a weiner QQ filter." << endl;
    LSDIndexRaster connected_components = raster.IsolateChannelsWienerQQ(this_float_map["pruning_drainage_area"],
                                                       this_float_map["surface_fitting_radius"], QQ_fname);

    cout << "Filtering by connected components." << endl;
    LSDIndexRaster connected_components_filtered = connected_components.filter_by_connected_components(this_int_map["connected_components_threshold"]);
    LSDIndexRaster CC_raster = connected_components_filtered.ConnectedComponents();

    cout << "thin network to skeleton" << endl;
    LSDIndexRaster skeleton_raster = connected_components_filtered.thin_to_skeleton();
    cout << "finding end points" << endl;
    LSDIndexRaster Ends = skeleton_raster.find_end_points();
    Ends.remove_downstream_endpoints(CC_raster, raster);


    cout << "Starting channel head processing" << endl;

    //this processes the end points to only keep the upper extent of the channel network
    cout << "getting channel heads" << endl;
    vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);
    cout << "processed all end points" << endl;

    // we need a temp junction network to search for single pixel channels
    LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
    LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

    cout << "removing single px channels" << endl;
    vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);

    // using these sources as the input to run the DrEICH algorithm  - FJC

    //Generate a channel netowrk from the sources
    LSDJunctionNetwork JunctionNetwork(FinalSources, FlowInfo);
    LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();

    LSDIndexRaster StreamNetwork = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

    // Calculate the channel head nodes
    int MinSegLength = 10;
    vector<int> ChannelHeadNodes_temp = JunctionNetwork.GetChannelHeadsChiMethodFromSources(FinalSources,
                    MinSegLength, this_float_map["A_0"], this_float_map["m_over_n"],
                    FlowInfo, DistanceFromOutlet, filled_topography, this_int_map["number_of_junctions_dreich"]);

    LSDIndexRaster Channel_heads_raster_temp = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(ChannelHeadNodes_temp);

    //create a channel network based on these channel heads
    LSDJunctionNetwork NewChanNetwork(ChannelHeadNodes_temp, FlowInfo);


    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
      string sources_csv_name = OUT_DIR+OUT_ID+"_Dsources.csv";

      //write channel_heads to a csv file
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(ChannelHeadNodes_temp, sources_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_Dsources.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_Dsources.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_Dsources";

      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(ChannelHeadNodes_temp);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_D_SO";

      //write stream order array to a raster
      LSDIndexRaster SOArray = NewChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }

    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_D_CN";
      if ( this_bool_map["use_extended_channel_data"])
      {
        cout << "I am going to use the extended channel network data outputs." << endl;
        NewChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, filled_topography);
      }
      else
      {
        NewChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
      }

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_D_CN.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_D_CN.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    // print junctions
    if( this_bool_map["print_junctions_to_csv"])
    {
      cout << "I am writing the junctions to csv." << endl;
      string channel_csv_name = OUT_DIR+OUT_ID+"_D_JN.csv";
      NewChanNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_D_JN.geojson";
        LSDSpatialCSVReader thiscsv(channel_csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    if (this_bool_map["print_channels_to_line_vector"])
    {
      string channel_geojson_name = OUT_DIR+OUT_ID+"_D_CN_line";
      cout << "I am going to print your channel network to a line geojson" << endl;
      LSDRaster DA = FlowInfo.write_DrainageArea_to_LSDRaster();
      NewChanNetwork.PrintChannelNetworkToLineVector(FlowInfo, channel_geojson_name, filled_topography, DA, DistanceFromOutlet);
    }

  }

  //===============================================================
  // PELLETIER
  //===============================================================
  if (this_bool_map["print_pelletier_channels"])
  {
    cout << "I am calculating channels using the Pelletier algorighm (doi:10.1029/2012WR012452)." << endl;
    cout << endl << "!!!!WARNING!!!" << endl;
    cout << "This routine is memory intensive! You will need ~20-30 times as much memory as the size of your DEM!" << endl;
    cout << "On a 3 Gb vagrant box a 100 Mb DEM is likeley to crash the machine." << endl;
    cout << "If you have this problem try reducing DEM resolution (3-5 m is okay, see Grieve et al 2016 ESURF)" << endl;
    cout << "or tile your DEM and run this multiple times. Or get a linux workstation." << endl << endl;
    LSDRasterSpectral SpectralRaster(filled_topography);

    cout << "I am running a Wiener filter" << endl;
    LSDRaster topo_test_wiener = SpectralRaster.fftw2D_wiener();
    int border_width = 100;
    topo_test_wiener = topo_test_wiener.border_with_nodata(border_width);

    cout << "I am going to fill and calculate flow info based on the filtered DEM." << endl;
    LSDRaster filter_fill = topo_test_wiener.fill(this_float_map["min_slope_for_fill"]);
    LSDFlowInfo FilterFlowInfo(boundary_conditions,filter_fill);

    // get some relevant rasters
    LSDRaster DistanceFromOutlet = FilterFlowInfo.distance_from_outlet();
    LSDIndexRaster ContributingPixels = FilterFlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    // get an initial sources network
    vector<int> sources;
    int pelletier_threshold = 250;
    sources = FilterFlowInfo.get_sources_index_threshold(ContributingPixels, pelletier_threshold);

    // now get an initial junction network. This will be refined in later steps.
    LSDJunctionNetwork ChanNetwork(sources, FilterFlowInfo);

    float surface_fitting_window_radius = this_float_map["surface_fitting_radius"];
    float surface_fitting_window_radius_LW = 25;
    vector<LSDRaster> surface_fitting, surface_fitting_LW;
    LSDRaster tan_curvature;
    LSDRaster tan_curvature_LW;
    string curv_name = "_tan_curv";
    vector<int> raster_selection(8, 0);
    raster_selection[6] = 1;

    // Two surface fittings: one for the short wacelength and one long wavelength
    surface_fitting = topo_test_wiener.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
    surface_fitting_LW = topo_test_wiener.calculate_polyfit_surface_metrics(surface_fitting_window_radius_LW, raster_selection);

    for(int i = 0; i<int(raster_selection.size()); ++i)
    {
      if(raster_selection[i]==1)
      {
        tan_curvature = surface_fitting[i];
        //tan_curvature.write_raster((path_name+DEM_name+curv_name), DEM_flt_extension);
        tan_curvature_LW = surface_fitting[i];
        //tan_curvature_LW.write_raster((path_name+DEM_name+curv_name+"_LW"), DEM_flt_extension);
      }
    }


    Array2D<float> topography = filled_topography.get_RasterData();
    Array2D<float> curvature = tan_curvature.get_RasterData();
    Array2D<float> curvature_LW = tan_curvature_LW.get_RasterData();
    cout << "\tLocating channel heads..." << endl;
    vector<int> ChannelHeadNodes = ChanNetwork.calculate_pelletier_channel_heads_DTM(FilterFlowInfo, topography, this_float_map["curvature_threshold"], curvature,curvature_LW);

    // Now filter out false positives along channel according to a threshold
    // catchment area
    cout << "\tFiltering out false positives..." << endl;
    LSDJunctionNetwork ChanNetworkNew(ChannelHeadNodes, FilterFlowInfo);
    vector<int> ChannelHeadNodesFilt;
    int count = 0;
    for(int i = 0; i<int(ChannelHeadNodes.size()); ++i)
    {
      int upstream_junc = ChanNetworkNew.get_Junction_of_Node(ChannelHeadNodes[i], FilterFlowInfo);
      int test_node = ChanNetworkNew.get_penultimate_node_from_stream_link(upstream_junc, FilterFlowInfo);
      float catchment_area = float(FilterFlowInfo.retrieve_contributing_pixels_of_node(test_node)) * FilterFlowInfo.get_DataResolution() * FilterFlowInfo.get_DataResolution();
      if (catchment_area >= this_float_map["minimum_drainage_area"])
      {
        ChannelHeadNodesFilt.push_back(ChannelHeadNodes[i]);
      }
      else
      {
        ++count;
      }
    }
    cout << "\t...removed " << count << " nodes out of " << ChannelHeadNodes.size() << endl;

    vector<int> FinalSources = ChannelHeadNodesFilt;

    //create a channel network based on these channel heads
    cout << "Making a channel network from the filtered channel heads." << endl;
    LSDJunctionNetwork NewChanNetwork(ChannelHeadNodesFilt, FilterFlowInfo);
    cout << "Got the network!" << endl;

    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
      string sources_csv_name = OUT_DIR+OUT_ID+"_Psources.csv";

      //write channel_heads to a csv file
      FilterFlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(FinalSources, sources_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_Psources.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_Psources.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_Psources";

      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FilterFlowInfo.write_NodeIndexVector_to_LSDIndexRaster(FinalSources);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_P_SO";

      //write stream order array to a raster
      LSDIndexRaster SOArray = NewChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }

    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_P_CN";
      if ( this_bool_map["use_extended_channel_data"])
      {
        cout << "I am going to use the extended channel network data outputs." << endl;
        NewChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, filled_topography);
      }
      else
      {
        NewChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
      }

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_P_CN.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_P_CN.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    // print junctions
    if( this_bool_map["print_junctions_to_csv"])
    {
      cout << "I am writing the junctions to csv." << endl;
      string channel_csv_name = OUT_DIR+OUT_ID+"_P_JN.csv";
      NewChanNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_P_JN.geojson";
        LSDSpatialCSVReader thiscsv(channel_csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    if (this_bool_map["print_channels_to_line_vector"])
    {
      string channel_geojson_name = OUT_DIR+OUT_ID+"_P_CN_line";
      cout << "I am going to print your channel network to a line geojson" << endl;
      LSDRaster DA = FlowInfo.write_DrainageArea_to_LSDRaster();
      NewChanNetwork.PrintChannelNetworkToLineVector(FlowInfo, channel_geojson_name, filled_topography, DA, DistanceFromOutlet);
    }

  }

  //===============================================================
  // WIENER
  //===============================================================
  if (this_bool_map["print_wiener_channels"])
  {
    cout << "I am calculating channels using the weiner algorithm (doi:10.1029/2012WR012452)." << endl;
    cout << "This algorithm was used by Clubb et al. (2016, DOI: 10.1002/2015JF003747) " << endl;
    cout << "and Grieve et al. (2016, doi:10.5194/esurf-4-627-2016) " << endl;
    cout << "and combines elements of the Pelletier and Passalacqua et al  methods: " << endl;
    cout << " doi:10.1029/2012WR012452 and doi:10.1029/2009JF001254" << endl;

    // initiate the spectral raster
    LSDRasterSpectral Spec_raster(filled_topography);

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
    vector<int> tmpsources = FlowInfo.ProcessEndPointsToChannelHeads(Ends);

    cout << "got all the end points" << endl;

    // we need a temp junction network to search for single pixel channels
    LSDJunctionNetwork tmpJunctionNetwork(tmpsources, FlowInfo);
    LSDIndexRaster tmpStreamNetwork = tmpJunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

    cout << "removing single px channels" << endl;
    vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);

    //Now we have the final channel heads, so we can generate a channel network from them
    LSDJunctionNetwork ChanNetwork(FinalSources, FlowInfo);

    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
      string sources_csv_name = OUT_DIR+OUT_ID+"_Wsources.csv";

      //write channel_heads to a csv file
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(FinalSources, sources_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_Wsources.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_Wsources.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_Wsources";

      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(FinalSources);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_W_SO";

      //write stream order array to a raster
      LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }

    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_W_CN";
      if ( this_bool_map["use_extended_channel_data"])
      {
        cout << "I am going to use the extended channel network data outputs." << endl;
        ChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, topography_raster);
      }
      else
      {
        ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
      }


      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_W_CN.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_W_CN.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    // print junctions
    if( this_bool_map["print_junctions_to_csv"])
    {
      cout << "I am writing the junctions to csv." << endl;
      string channel_csv_name = OUT_DIR+OUT_ID+"_W_JN.csv";
      ChanNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_W_JN.geojson";
        LSDSpatialCSVReader thiscsv(channel_csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    if (this_bool_map["print_channels_to_line_vector"])
    {
      string channel_geojson_name = OUT_DIR+OUT_ID+"_W_CN_line";
      cout << "I am going to print your channel network to a line geojson" << endl;
      // distance from outlet
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
      LSDRaster DA = FlowInfo.write_DrainageArea_to_LSDRaster();
      ChanNetwork.PrintChannelNetworkToLineVector(FlowInfo, channel_geojson_name, filled_topography, DA, DistanceFromOutlet);
    }

  }

  //===============================================================
  // PERONA-MALIK
  //===============================================================
  if (this_bool_map["print_perona_malik_channels"])
  {
    cout << "I am calculating channels using a QQ threshold algorithm (doi:10.1029/2012WR012452)." << endl;
    cout << "and a Perona-Malik filter similar to GeoNet (Passalacqua et al, 2010): " << endl;
    cout << "doi:10.1029/2009JF001254" << endl;

    // initiate the spectral raster
    LSDRasterSpectral Spec_raster(filled_topography);

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
    vector<int> FinalSources = FlowInfo.RemoveSinglePxChannels(tmpStreamNetwork, tmpsources);

    //Now we have the final channel heads, so we can generate a channel network from them
    LSDJunctionNetwork ChanNetwork(FinalSources, FlowInfo);

    // Print sources
    if( this_bool_map["print_sources_to_csv"])
    {
      string sources_csv_name = OUT_DIR+OUT_ID+"_PMsources.csv";

      //write channel_heads to a csv file
      FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(FinalSources, sources_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_PMsources.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_PMsources.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }

    }

    if( this_bool_map["print_sources_to_raster"])
    {
      string sources_raster_name = OUT_DIR+OUT_ID+"_PMsources";

      //write channel heads to a raster
      LSDIndexRaster Channel_heads_raster = FlowInfo.write_NodeIndexVector_to_LSDIndexRaster(FinalSources);
      Channel_heads_raster.write_raster(sources_raster_name,raster_ext);
    }

    if( this_bool_map["print_stream_order_raster"])
    {
      string SO_raster_name = OUT_DIR+OUT_ID+"_PM_SO";

      //write stream order array to a raster
      LSDIndexRaster SOArray = ChanNetwork.StreamOrderArray_to_LSDIndexRaster();
      SOArray.write_raster(SO_raster_name,raster_ext);
    }

    if( this_bool_map["print_channels_to_csv"])
    {
      string channel_csv_name = OUT_DIR+OUT_ID+"_PM_CN";
      if ( this_bool_map["use_extended_channel_data"])
      {
        cout << "I am going to use the extended channel network data outputs." << endl;
        ChanNetwork.PrintChannelNetworkToCSV_WithElevation_WithDonorJunction(FlowInfo, channel_csv_name, topography_raster);
      }
      else
      {
        ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
      }


      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_PM_CN.geojson";
        LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_PM_CN.csv");
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    // print junctions
    if( this_bool_map["print_junctions_to_csv"])
    {
      cout << "I am writing the junctions to csv." << endl;
      string channel_csv_name = OUT_DIR+OUT_ID+"_PM_JN.csv";
      ChanNetwork.print_junctions_to_csv(FlowInfo, channel_csv_name);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_PM_JN.geojson";
        LSDSpatialCSVReader thiscsv(channel_csv_name);
        thiscsv.print_data_to_geojson(gjson_name);
      }
    }

    if (this_bool_map["print_channels_to_line_vector"])
    {
      string channel_geojson_name = OUT_DIR+OUT_ID+"_PM_CN_line";
      cout << "I am going to print your channel network to a line geojson" << endl;
      // distance from outlet
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
      LSDRaster DA = FlowInfo.write_DrainageArea_to_LSDRaster();
      ChanNetwork.PrintChannelNetworkToLineVector(FlowInfo, channel_geojson_name, filled_topography, DA, DistanceFromOutlet);
    }

  }

  //============================================================================
  // Write some rasters. THis is inefficient because it could duplicate
  // Creation of rasters above, but means the different extraction
  // switches don't need to be turned on for these rasters to be printed
  //============================================================================
  if (this_bool_map["print_wiener_filtered_raster"])
  {

    cout << "I am running a filter to print to raster" << endl;
    LSDRasterSpectral SpectralRaster(topography_raster);
    LSDRaster topo_test_wiener = SpectralRaster.fftw2D_wiener();

    string wiener_name = OUT_DIR+OUT_ID+"_Wfilt";
    topo_test_wiener.write_raster(wiener_name,raster_ext);
  }

  // Prints the curvature raster if you want it
  if (this_bool_map["print_curvature_raster"])
  {
    vector<int> raster_selection(8, 0);
    raster_selection[6] = 1;

    float surface_fitting_window_radius = this_float_map["surface_fitting_radius"];
    float surface_fitting_window_radius_LW = 25;
    vector<LSDRaster> surface_fitting, surface_fitting_LW;

    cout << "I am printing curvature rasters for you. These are not filtered!" << endl;

    // Two surface fittings: one for the short wavelength and one long wavelength
    surface_fitting = topography_raster.calculate_polyfit_surface_metrics(surface_fitting_window_radius, raster_selection);
    surface_fitting_LW = topography_raster.calculate_polyfit_surface_metrics(surface_fitting_window_radius_LW, raster_selection);

    string curv_name = OUT_DIR+OUT_ID+"_tan_curv";
    string curv_name_LW = OUT_DIR+OUT_ID+"_tan_curv_LW";

    surface_fitting[6].write_raster(curv_name,raster_ext);
    surface_fitting_LW[6].write_raster(curv_name_LW,raster_ext);
  }





  // Done, check how long it took
  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  cout << "DONE! Have fun with your channels! Time taken (secs): " << elapsed_secs << endl;

}
