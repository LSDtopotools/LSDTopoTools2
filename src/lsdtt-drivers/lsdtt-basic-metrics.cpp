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
// https://lsdtopotools.github.io/LSDTopoTools_ChiMudd2014/
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2019 Simon M. Mudd 2019
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

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the LSDTopoTools basic metrics tool!     ||" << endl;
  cout << "|| This program has a number of options for calculating||" << endl;
  cout << "|| simple landscape metrics.                           ||" << endl;
  cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| If you use these routines please cite:              ||" << endl;
  cout << "|| https://www.doi.org/10.5281/zenodo.2560223          ||" << endl;
  cout << "|| If you use the roughness routine please cite:       ||" << endl;
  cout << "|| https://www.doi.org/10.5194/esurf-3-483-2015        ||" << endl;
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
  bool_default_map["raster_is_filled"] = false; // assume base raster is already filled
  bool_default_map["carve_before_fill"] = false; // Implements a carving algorithm
  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  string_default_map["CHeads_file"] = "NULL";
  bool_default_map["only_check_parameters"] = false;

  // Various ways of trimming your raster
  bool_default_map["remove_nodes_influenced_by_edge"] = false;
  bool_default_map["isolate_pixels_draining_to_fixed_channel"] = false;
  string_default_map["fixed_channel_csv_name"] = "single_channel_nodes";
  
  // Parameters for swath mapping
  bool_default_map["calculate_swath_profile"] = false;
  bool_default_map["calculate_swath_along_line"] = false;
  bool_default_map["calculate_swath_along_channel"] = false;
  bool_default_map["print_swath_rasters"] = true;
  string_default_map["swath_points_csv"] = "swath.csv";
  float_default_map["swath_point_spacing"] = 500;
  float_default_map["swath_bin_spacing"] = 1000;
  float_default_map["swath_width"] = 1000;

  // some tools for punching out polygons and then getting distance to nearest nodata
  // used for creating rivers in noisy DEMs
  bool_default_map["punch_nodata"] = false;
  bool_default_map["print_nearest_to_nodata_rasters"] = false; 
  string_default_map["punch_raster_prefix"] = "DEM_punch";
  bool_default_map["belowthresholdisnodata"] = true;
  float_default_map["punch_threshold"] = 1000;
  float_default_map["minimum_bank_elevation_window_radius"] = 20;
  float_default_map["river_depth"] = 1.0;

  // More river processing tools
  bool_default_map["channel_and_valley_width_extraction"] = false;
  string_default_map["channel_or_valley_raster_prefix"] = "channel";
  string_default_map["channel_or_valley_skeleton_prefix"] = "skeleton";
  int_default_map["test_scale"] = 4;
  float_default_map["test_bearing"] = 90;

  // raster trimming, to take care of rasters that have a bunch of nodata at the edges
  bool_default_map["print_trimmed_raster"] = false;
  int_default_map["trimming_buffer_pixels"] = 0;



  // Calculate the basic relief within a window in a raster
  bool_default_map["print_relief_raster"] = false;
  float_default_map["relief_window"] = 200;
  int_default_map["relief_window_kernel_type"] = 1;


  // the most basic raster printing
  bool_default_map["write_hillshade"] = false;
  bool_default_map["print_raster_without_seas"] = false;
  bool_default_map["print_distance_from_outlet"] = false;
  bool_default_map["print_fill_raster"] = false;

  // This converts all csv files to geojson (for easier loading in a GIS)
  bool_default_map["convert_csv_to_geojson"] = false;

  // Slope calculations
  float_default_map["surface_fitting_radius"] = 30;
  bool_default_map["print_smoothed_elevation"]= false;
  bool_default_map["print_slope"] = false;
  bool_default_map["print_aspect"]= false;
  bool_default_map["print_curvature"]= false;
  bool_default_map["print_planform_curvature"]= false;
  bool_default_map["print_profile_curvature"]= false;
  bool_default_map["print_tangential_curvature"]= false;
  bool_default_map["print_point_classification"]= false;
  bool_default_map["print_directional_gradients"] = false;
  bool_default_map["calculate_basin_statistics"] = false;

  // Window size estimation
  bool_default_map["calculate_window_size"] = false;

  // Roughness calculations
  float_default_map["REI_critical_slope"] = 1.0;
  float_default_map["REI_window_radius"] = 10;
  bool_default_map["print_REI_raster"] = false;

  bool_default_map["print_roughness_rasters"] = false;
  float_default_map["roughness_radius"] = 3;



  // drainage area
  bool_default_map["print_dinf_drainage_area_raster"] = false;
  bool_default_map["print_d8_drainage_area_raster"] = false;
  bool_default_map["print_QuinnMD_drainage_area_raster"] = false;
  bool_default_map["print_FreemanMD_drainage_area_raster"] = false;
  bool_default_map["print_MD_drainage_area_raster"] = false;


  // Extracting a single channel
  bool_default_map["extract_single_channel"] = false;
  bool_default_map["use_dinf_for_single_channel"] = false;
  string_default_map["channel_source_fname"] = "single_channel_source";

  // imposing a single channel
  bool_default_map["impose_single_channel"] = false;
  string_default_map["fixed_channel_csv_name"] = "NULL";
  bool_default_map["buffer_single_channel"] = true;
  bool_default_map["use_XY_for_buffer"] = false;
  bool_default_map["force_single_channel_slope"] = false;
  bool_default_map["fixed_channel_dig_and_slope"] = false;
  float_default_map["single_channel_drop"] = 0;
  float_default_map["single_channel_dig"] = 0.01;
  string_default_map["single_channel_fd_string"] = "flow distance(m)";
  string_default_map["single_channel_elev_string"] = "elevation(m)";
  string_default_map["test_single_channel_name"] = "test_single_channel";

  // for reading in a channel csv file
  bool_default_map["use_xy_for_node_index"] = false;
  


  // Basic channel network
  int_default_map["threshold_contributing_pixels"] = 1000;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_channels_to_csv"] = false;
  bool_default_map["use_extended_channel_data"] = false;
  bool_default_map["print_junction_index_raster"] = false;
  bool_default_map["print_junctions_to_csv"] = false;

  // Basin-based channel extraction
  bool_default_map["find_basins"] = false;
  int_default_map["minimum_basin_size_pixels"] = 50000;
  int_default_map["maximum_basin_size_pixels"] = 1000000;
  bool_default_map["only_take_largest_basin"] = false;
  string_default_map["BaselevelJunctions_file"] = "NULL";
  bool_default_map["get_basins_from_outlets"] = false;
  int_default_map["search_radius_nodes"] = 8;
  int_default_map["threshold_stream_order_for_snapping"] = 2;
  string_default_map["basin_outlet_csv"] = "NULL";
  bool_default_map["extend_channel_to_node_before_receiver_junction"] = true;
  bool_default_map["print_basin_raster"] = false;

  // finding major drainage divides
  bool_default_map["divide_finder"] = false;
  bool_default_map["horizontal_strips"] = true;
  int_default_map["n_row_or_col_for_strips"] = 10;

  // Getting all the ridges
  bool_default_map["extract_ridges"] = false;


  // Tagging pixels
  bool_default_map["tag_nodes"] = false;
  string_default_map["tagged_raster_input_name"] = "NULL";
  bool_default_map["tag_downslope_nodes"] = false;
  float_default_map["downslope_tagging_distance"] = 100;
  float_default_map["upslope_tagging_distance"] = 100;

  // Some chi coordinate settings
  float_default_map["A_0"] = 1.0;
  float_default_map["m_over_n"] = 0.5;
  bool_default_map["print_chi_data_maps"] = false;

  if (bool_default_map["print_chi_data_maps"] == true)
  {
    cout << "You want the chi data maps, so I am setting the basin finding to true." << endl;
    bool_default_map["print_chi_data_maps"] = true;
  }

  // The wiener filter
  bool_default_map["print_wiener_filtered_raster"] = false;

  // This burns a raster value to any csv output of chi data
  // Useful for appending geology data to chi profiles
  bool_default_map["burn_raster_to_csv"] = false;
  string_default_map["burn_raster_prefix"] = "NULL";
  string_default_map["burn_data_csv_column_header"] = "burned_data";
  string_default_map["csv_to_burn_name"] = "NULL";

  // This is for junction angles
  bool_default_map["print_junction_angles_to_csv"] = false;
  float_default_map["SA_vertical_interval"] = 10;  // The vertical interval
                                                   // for slope measurements in the
                                                   // junction angle printing

  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

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

    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);

  // A little switch to find basins if you supply an outlet file
  if( this_bool_map["get_basins_from_outlets"] )
  {
    this_bool_map["find_basins"] = true;
  }


  // load the  DEM
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
    string swath_data_prefix = OUT_DIR+OUT_ID;
    topography_raster.make_swath(spaced_eastings, spaced_northings,spaced_distances, this_float_map["swath_width"],
                                 this_float_map["swath_bin_spacing"], swath_data_prefix, this_bool_map["print_swath_rasters"]);



  }


  //=========================================================================
  // .##..##...####...##......##......######..##..##..........##...##..######..#####...######..##..##.
  // .##..##..##..##..##......##......##.......####...........##...##....##....##..##....##....##..##.
  // .##..##..######..##......##......####......##............##.#.##....##....##..##....##....######.
  // ..####...##..##..##......##......##........##............#######....##....##..##....##....##..##.
  // ...##....##..##..######..######..######....##.............##.##...######..#####.....##....##..##.
  //
  //=========================================================================
  if(this_bool_map["channel_and_valley_width_extraction"])
  {
    cout << "I am going to test some code for valley and channel widths" << endl;


    Array2D<int> test = make_template_for_vector_bearing(this_float_map["test_bearing"], this_int_map["test_scale"]);
  }

  // Some tools for punching out some nodata using a secondary raster that we can use to make river
  // pathways in noisy DEMs. 
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
    Punched.write_raster(punch_name,raster_ext);

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
  // Raster burning

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
  // Print the roughness rasters
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
      RM.impose_channels_with_buffer_use_XY(single_channel_data, this_float_map["min_slope_for_fill"]*2);
    }
    else if(this_bool_map["buffer_single_channel"] and not this_bool_map["use_XY_for_buffer"])
    {
      cout << "I am going to impose and buffer the single channel using lat-long coordinates." << endl;
      RM.impose_channels_with_buffer(single_channel_data, this_float_map["min_slope_for_fill"]*2);
    }
    else
    {
      cout << "I am going to impose the single channel using lat-long coordinates, without buffering." << endl;
      RM.impose_channels(single_channel_data);
    }
    LSDRaster imposed_topography = RM.return_as_raster();

    cout << "Let me print the fill raster for you."  << endl;
    string filled_raster_name = OUT_DIR+OUT_ID+"_imposed_channels";
    imposed_topography.write_raster(filled_raster_name,raster_ext);       
  }

  //============================================================================
  //
  //.######..##.......####...##...##..........#####....####...##..##..######..######..##..##...####..
  //.##......##......##..##..##...##..........##..##..##..##..##..##....##......##....###.##..##.....
  //.####....##......##..##..##.#.##..........#####...##..##..##..##....##......##....##.###..##.###.
  //.##......##......##..##..#######..........##..##..##..##..##..##....##......##....##..##..##..##.
  //.##......######...####....##.##...........##..##...####....####.....##....######..##..##...####..
  //.................................................................................................
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
        || this_bool_map["print_chi_data_maps"]
        || this_bool_map["print_junction_angles_to_csv"]
        || this_bool_map["extract_single_channel"]
        || this_bool_map["remove_nodes_influenced_by_edge"]
        || this_bool_map["isolate_pixels_draining_to_fixed_channel"]
        || this_bool_map["calculate_basin_statistics"]
        || this_bool_map["divide_finder"]
        || this_bool_map["tag_nodes"]
        || this_bool_map["extract_ridges"])
  {
    cout << "I will need to compute flow information, because you are getting drainage area or channel networks." << endl;
    //==========================================================================
    // Fill the raster
    //==========================================================================
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




    //=================================================================
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


    // This is the logic for a simple stream network
    if (this_bool_map["print_channels_to_csv"]
        || this_bool_map["print_junctions_to_csv"]
        || this_bool_map["print_sources_to_csv"]
        || this_bool_map["find_basins"]
        || this_bool_map["print_chi_data_maps"]
        || this_bool_map["print_junction_angles_to_csv"]
        || this_bool_map["extract_ridges"])
    {
      // calculate the flow accumulation
      cout << "\t Calculating flow accumulation (in pixels)..." << endl;
      LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

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
          string gjson_name = OUT_DIR+OUT_ID+"_JAngles.geojson";
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

      // Now we check if we are going to deal with basins
      if(this_bool_map["find_basins"] ||
         this_bool_map["print_chi_data_maps"] ||
         this_bool_map["calculate_basin_statistics"] )
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
            cout << "I am going to select basins for you using an algorithm. " << endl;
            cout << "I am going to look for basins in a pixel window that are not influended by nodata." << endl;
            cout << "I am also going to remove any nested basins." << endl;
            cout << "The pixel limits are: lower: " << this_int_map["minimum_basin_size_pixels"] << " and upper: " << this_int_map["maximum_basin_size_pixels"] << endl;
            BaseLevelJunctions = JunctionNetwork.Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(FlowInfo, filled_topography, FlowAcc,
                                                      this_int_map["minimum_basin_size_pixels"],this_int_map["maximum_basin_size_pixels"]);
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

          // Now make sure none of the basins drain to the edge
          cout << "I am pruning junctions that are influenced by the edge of the DEM!" << endl;
          cout << "This is necessary because basins draining to the edge will have incorrect chi values." << endl;
          BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge_Ignore_Outlet_Reach(BaseLevelJunctions_Initial, FlowInfo, filled_topography);

          cout << "The remaining baselevel junctions are: " << endl;
          for (int i = 0; i<int(BaseLevelJunctions.size()); i++)
          {
            cout << BaseLevelJunctions[i] << endl;
          }

        }    // end logic for reading from junctions list

        // Now check for largest basin, if that is what you want.
        if (this_bool_map["only_take_largest_basin"])
        {
          cout << "I am only going to take the largest basin." << endl;
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


        // Now we get the channel segments. This information is used for plotting
        vector<int> source_nodes;
        vector<int> outlet_nodes;
        vector<int> baselevel_node_of_each_basin;
        int n_nodes_to_visit = 10;
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

          JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, FD,
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


        // Get the chi coordinate if needed
        LSDRaster chi_coordinate;
        if ( this_bool_map["print_basin_raster"] ||
             this_bool_map["print_chi_data_maps"])
        {
          cout << "I am calculating the chi coordinate." << endl;
          chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(this_float_map["m_over_n"],this_float_map["A_0"],this_int_map["threshold_contributing_pixels"]);

        }

        //======================================================================
        // Print a basin raster if you want it.
        if(this_bool_map["print_basin_raster"])
        {
          cout << "I am going to print the basins for you. " << endl;
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
        //======================================================================

      }   // end logic for basin finding
      cout << "Finished with basins" << endl;
    }     // end logic for tasks related to channel network extraction
    cout << "Done with channel extraction" << endl;
  }       // end logic for tasks requiring flow info and filling
  cout << "I'm all finished! Have a nice day." << endl;
}
