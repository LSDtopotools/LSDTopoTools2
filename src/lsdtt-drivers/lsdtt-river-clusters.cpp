//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// get profiles for river clustering.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2018 Fiona Clubb
//
// Developer can be contacted by clubb _at_ uni-potsdam.de
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "../LSDParameterParser.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDBasin.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDSpatialCSVReader.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the river cluster tool, developed by     ||" << endl;
    cout << "|| Fiona Clubb at the University of Potsdam            ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "lsdtt-river-clusters /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_spaghetti.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);

  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the basin tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // parameters for DEM preprocessing
  float_default_map["min_slope_for_fill"] = 0.0001;
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_filled_raster"] = false;
  bool_default_map["write hillshade"] = false;
  bool_default_map["write_slope_raster"] = false;
  int_default_map["surface_fitting_window_radius"] = 6;

  // parameters for channel extraction
  int_default_map["threshold_contributing_pixels"] = 2000;
  string_default_map["channel heads fname"] = "NULL";
  bool_default_map["print_channels_to_csv"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_junction_index_raster"] = false;

  // parameters for basin selection
  bool_default_map["select_basins_by_order"] = false;
  bool_default_map["all_basins"] = false;
  int_default_map["basin_order"] = 3;
  bool_default_map["print_basin_raster"] = false;
  bool_default_map["print_junctions_to_csv"] = false;
  int_default_map["junction_number"] = 0;

  // parameters for getting the profile plots
  bool_default_map["print_spaghetti_profiles_to_csv"] = false;
  bool_default_map["print_all_tributaries_to_csv"] = false;
  bool_default_map["print_channels_all_sources_to_csv"] = false;
  bool_default_map["get_catchment_info_clusters"] = false;  //only use if you have already done the clustering
  int_default_map["stream_order"] = 1; // stream order that you clustered by
  int_default_map["n_clusters"] = 2; // the number of clusters. needs to be specified for csv reading
  int_default_map["threshold_level"] = 0; // the threshold level from the clustering that you want to select. Either 1 or 0
  bool_default_map["vegetation"] = false; // turn to true if you have a veg height raster and want to include it in the analysis
  float_default_map["channel_length"] = 1000;
  int_default_map["slope_window_size"] = 1;

  // Use the parameter parser to get the maps of the parameters required for the
  // analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();

  cout << "Read filename is:" <<  DATA_DIR+DEM_ID << endl;

    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
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


  if (this_bool_map["print_filled_raster"])
  {
    string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }

  // check to see if you need hillshade
  if (this_bool_map["write hillshade"])
  {
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+OUT_ID+"_hs";
    hs_raster.write_raster(hs_fname,raster_ext);
  }

  if (this_bool_map["write_slope_raster"])
  {
    vector<LSDRaster> surface_fitting;
    vector<int> raster_selection(8, 0);
    raster_selection[1] = 1;
    surface_fitting = filled_topography.calculate_polyfit_surface_metrics(this_int_map["surface_fitting_window_radius"], raster_selection);
    LSDRaster Slope = surface_fitting[1];
    string slope_name = OUT_DIR+OUT_ID+"_slope";
    Slope.write_raster(slope_name, raster_ext);
  }

  //============================================================================
  // Get the flow info
  //============================================================================
  cout << "\t Flow routing..." << endl;
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

  // now deal with the channel network
  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
  {
    // calculate the flow accumulation
    cout << "\t Calculating flow accumulation (in pixels)..." << endl;
    LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    // Now get the sources from flow accumulation
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

    cout << "The number of sources is: " << sources.size() << endl;

  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t Got sources!" << endl;
  }

  // Now create the network
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

  // get the distance from outlet
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  if( this_bool_map["print_channels_to_csv"])
  {
    cout << "I am writing the channels to csv." << endl;
    string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
    JunctionNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);

    if ( this_bool_map["convert_csv_to_geojson"])
    {
      string gjson_name = OUT_DIR+OUT_ID+"_CN.geojson";
      LSDSpatialCSVReader thiscsv(OUT_DIR+OUT_ID+"_CN.csv");
      thiscsv.print_data_to_geojson(gjson_name);
    }

  }

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

  if (this_bool_map["print_stream_order_raster"])
  {
    LSDIndexRaster SOArray = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
    string SO_raster_name = OUT_DIR+OUT_ID+"_SO";
    SOArray.write_raster(SO_raster_name,raster_ext);
  }

  if (this_bool_map["print_junction_index_raster"])
  {
    LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
    string JI_raster_name = OUT_DIR+OUT_ID+"_JI";
    JIArray.write_raster(JI_raster_name,raster_ext);
  }


  if ( this_int_map["junction_number"] != 0 )
  {
    // Pad the DEMs by 1 pixel
    //int padding_pixels = 2;

    cout << "You supplied a junction number, so I'll only analyse channels in that basin." << endl;
    cout << "Getting basin for baselevel junction " << this_int_map["junction_number"] << endl;
    LSDBasin thisBasin(this_int_map["junction_number"],FlowInfo, JunctionNetwork);


    // name the new DEM file
    string INT = string(itoa(this_int_map["junction_number"]));
    string basin_fname = DEM_ID+"_basin_" + INT;

    // write the channel network
    string cn_name = (OUT_DIR+basin_fname+"_CN.csv");
    thisBasin.write_channel_network(cn_name, FlowInfo, JunctionNetwork);

    if (this_bool_map["print_basin_raster"])
    {
        LSDRaster BasinElev = thisBasin.write_raster_data_to_LSDRaster(filled_topography, FlowInfo);
        BasinElev.write_raster(OUT_DIR+basin_fname,raster_ext);
    }

    if ( this_bool_map["print_all_tributaries_to_csv"])
    {
      string csv_filename = OUT_DIR+DEM_ID+"_all_tribs";
      vector<int> BasinJunctions;
      BasinJunctions.push_back(this_int_map["junction_number"]);
      JunctionNetwork.write_river_profiles_to_csv_all_tributaries(BasinJunctions, FlowInfo, DistanceFromOutlet, filled_topography, csv_filename, this_int_map["slope_window_size"]);
    }
  }
  if ( this_bool_map["select_basins_by_order"] )
  {
    // Now get all the basins of the selected order
    vector<int> BasinJunctions = JunctionNetwork.extract_basins_order_outlet_junctions(this_int_map["basin_order"], FlowInfo);

    if( this_bool_map["print_spaghetti_profiles_to_csv"])
    {
      string csv_filename = OUT_DIR+OUT_ID+"_spaghetti_profiles.csv";
      JunctionNetwork.write_river_profiles_to_csv(BasinJunctions, FlowInfo, DistanceFromOutlet, filled_topography, csv_filename, this_int_map["slope_window_size"]);
    }

    // print profiles for all tributaries for each third order basin
    if( this_bool_map["print_all_tributaries_to_csv"])
    {
      string csv_filename = OUT_DIR+OUT_ID+"_all_tribs";
      JunctionNetwork.write_river_profiles_to_csv_all_tributaries(BasinJunctions, FlowInfo, DistanceFromOutlet, filled_topography, csv_filename, this_int_map["slope_window_size"]);
    }

    // print profiles for all tributaries for each third order basin
    if( this_bool_map["print_channels_all_sources_to_csv"])
    {
      string csv_filename = OUT_DIR+OUT_ID+"_all_sources"+itoa(this_float_map["channel_length"]);
      JunctionNetwork.write_river_profiles_to_csv_all_sources(this_float_map["channel_length"], this_int_map["slope_window_size"], FlowInfo, filled_topography, csv_filename);
    }

    if( this_bool_map["print_basin_raster"])
    {
      LSDIndexRaster BasinRaster = JunctionNetwork.extract_basins_from_junction_vector(BasinJunctions, FlowInfo);
      string basins_name = OUT_DIR+OUT_ID+"_basins";
      BasinRaster.write_raster(basins_name, raster_ext);
    }
  }

  if ( this_bool_map["all_basins"] )
  {
    cout << "I am going to get all of the basins (every base level junction will be used). This might take a while!!" << endl;
    // get all the base level junctions
    vector<int> BasinJunctions = JunctionNetwork.get_BaseLevelJunctions();
    // print profiles for all tributaries for each third order basin
    if( this_bool_map["print_all_tributaries_to_csv"])
    {
      string csv_filename = OUT_DIR+OUT_ID+"_all_tribs";
      JunctionNetwork.write_river_profiles_to_csv_all_tributaries(BasinJunctions, FlowInfo, DistanceFromOutlet, filled_topography, csv_filename, this_int_map["slope_window_size"]);
    }
  }

  // get the upstream catchments of all the profiles in each cluster
  if ( this_bool_map["get_catchment_info_clusters"])
  {
    // get the directory for this clustering level
    string this_dir = DATA_DIR+"threshold_"+string(itoa(this_int_map["threshold_level"]))+"/";

    LSDCoordinateConverterLLandUTM Converter;
    float NoDataValue = filled_topography.get_NoDataValue();
    // read csvs up until n clusters
    string string_SO = string(itoa(this_int_map["stream_order"]));
    string cl = string(itoa(this_int_map["n_clusters"]));

    //surface fitting
    vector<int> raster_selection;

    raster_selection.push_back(0);
    raster_selection.push_back(1); //slope
    raster_selection.push_back(1); //aspect
    raster_selection.push_back(1); //curvature
    raster_selection.push_back(1); //plan curvature
    raster_selection.push_back(1); //profile curvature
    raster_selection.push_back(0);
    raster_selection.push_back(0);

    vector<LSDRaster> Surfaces = filled_topography.calculate_polyfit_surface_metrics(this_int_map["surface_fitting_window_radius"], raster_selection);
    LSDRaster TotalCurv = Surfaces[3];
    LSDRaster ProfileCurv = Surfaces[5];
    LSDRaster PlanCurv = Surfaces[4];
    LSDRaster Aspect = Surfaces[2];
    LSDRaster Slope = Surfaces[1];

    // get the stream network as an index raster
    LSDIndexRaster StreamNetwork = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();

    // LSDSpatialCSVReader Veg(filled_topography, DATA_DIR+DEM_ID+"_veg_height.csv");

    // calculate the relief
    LSDRaster Relief = filled_topography.calculate_relief(20, 1);

    for (int i = 1; i < this_int_map["n_clusters"]+1; i++)
    {
      // read the csv
      string csv_filename = this_dir+DEM_ID+"_junctions_SO"+string_SO+"_CL"+string(itoa(i))+".csv";
      LSDSpatialCSVReader Clusters(csv_filename);

      cout << "Reading the file: " << csv_filename << endl;

      // get the downstream_nodes
      vector<double> lat = Clusters.get_latitude();
      vector<double> lon =  Clusters.get_longitude();
      vector<int> nodes = Clusters.data_column_to_int("node");
      vector<int> ids = Clusters.data_column_to_int("id");
      vector<int> basin_junctions;

      // set up the csv with the mean data from each catchment
      string out_filename = this_dir+DEM_ID+"_catchment_info_SO"+string_SO+"_CL"+string(itoa(i))+".csv";
      // open the csv
      ofstream csv_out;
      csv_out.open(out_filename.c_str());
      if ( this_bool_map["vegetation"] )
      {
        csv_out << "id,basin_junc,roughness,mean_slope,total_curv,plan_curv,drainage_density,veg_height,latitude,longitude" << endl;
      }
      else
      {
        csv_out << "id,basin_junc,roughness,mean_slope,total_curv,plan_curv,drainage_density,latitude,longitude" << endl;
      }
      // loop through these nodes and get the basins
      for (int j = 0; j < int(lat.size()); j++)
      {
        // get the nearest junction to this node and initialise the basin object
        int this_junc = JunctionNetwork.find_upstream_junction_from_channel_nodeindex(nodes[j], FlowInfo);
        if (this_junc != NoDataValue)
        {
          cout << "This junction is: " << this_junc << endl;
          basin_junctions.push_back(this_junc);
          LSDBasin Basin(this_junc, FlowInfo, JunctionNetwork);

          // now - what analysis do we want to do?
          // we want the catchment relief
          //float MaxElev = Basin.CalculateBasinMax(FlowInfo, filled_topography);
          //float MinElev = Basin.CalculateBasinMin(FlowInfo, filled_topography);
          //float Relief = MaxElev - MinElev;
          float MeanRoughness = Basin.CalculateBasinMean(FlowInfo, Relief);

          // slopes
          Basin.set_SlopeMean(FlowInfo, Slope);
          float MeanSlope = Basin.get_SlopeMean();

          //curvature
          Basin.set_TotalCurvMean(FlowInfo, TotalCurv);
          float MeanTotalCurv = Basin.get_TotalCurvMean();

          Basin.set_PlanCurvMean(FlowInfo, PlanCurv);
          float MeanPlanCurv = Basin.get_PlanCurvMean();

          // drainage density
          Basin.set_FlowLength(StreamNetwork, FlowInfo);
          Basin.set_DrainageDensity();
          float DrainageDensity = Basin.get_DrainageDensity();

          // if ( this_bool_map["vegetation"] )
          // {
          //   float MeanVeg = Basin.get_basin_mean_from_csv(FlowInfo, Veg, "value");
          //   csv_out << ids[j] << "," << this_junc << "," << MeanRoughness << "," << MeanSlope << "," << MeanTotalCurv  << "," << MeanPlanCurv << "," << DrainageDensity << "," << MeanVeg << ",";
          //   csv_out.precision(9);
          //   csv_out << lat[j] << "," << lon[j] << endl;
          // }
          // else
          // {
          csv_out << ids[j] << "," << this_junc << "," << MeanRoughness << "," << MeanSlope << "," << MeanTotalCurv  << "," << MeanPlanCurv << "," << DrainageDensity << ",";
          csv_out.precision(9);
          csv_out << lat[j] << "," << lon[j] << endl;
          // }

        }
      }
      csv_out.close();

      // write the basins to a raster for visualisation
      LSDIndexRaster Basins = JunctionNetwork.extract_basins_from_junction_vector(basin_junctions, FlowInfo);
      string basins_name = "_basins_CL";
      Basins.write_raster(this_dir+DEM_ID+"_basins_SO"+string_SO+"_CL"+string(itoa(i)), raster_ext);
    }
  }


}
