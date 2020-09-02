//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// lsdtt-cosmo-tool.cpp
//
// This is a program for dealing with cosmogenic data.
// It has routines for both calculating erosion rates on the 
// basis of measured CRN concentrations (10Be, 26Al, etc)
// and also forward modelling of different erosion scenarios
// and predicting cosmogenic concentrations 
//
// This program takes two arguments, the path name and the driver name
// It is an updated version of the CAIRN software package:
// https://www.earth-surf-dynam.net/4/655/2016/
//
// Mudd, S. M., Harel, M.-A., Hurst, M. D., Grieve, S. W. D., and Marrero, S. M.: 
// The CAIRN method: automated, reproducible calculation of catchment-averaged denudation 
// rates from cosmogenic nuclide concentrations, Earth Surf. Dynam., 
// 4, 655–674, https://doi.org/10.5194/esurf-4-655-2016, 2016.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2020 Simon M. Mudd 2020
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
#include <chrono>  // for high_resolution_clock
#include <string>
#include <vector>
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
#include "../LSDSoilHydroRaster.hpp"
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDRasterSpectral.hpp"
#include "../LSDChiTools.hpp"
#include "../LSDStrahlerLinks.hpp"
#include "../LSDBasin.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDCosmoData.hpp"
#include "../LSDCosmoRaster.hpp"
#include "../LSDRasterMaker.hpp"
#include "../LSDRasterAggregator.hpp"

int main (int nNumberofArgs,char *argv[])
{

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the LSDTopoTools cosmogenic tool!        ||" << endl;
  cout << "|| This program has a number of options for calculating||" << endl;
  cout << "|| cosmogenic nuclide concentrations.                  ||" << endl;
  cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;
  cout << "=========================================================" << endl;   
  cout << "|| If you use these routines please cite:              ||" << endl;   
  cout << "|| https://doi.org/10.5194/esurf-4-655-2016            ||" << endl;
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
  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  string_default_map["CHeads_file"] = "NULL";
  bool_default_map["only_check_parameters"] = false;
  
  // the most basic raster printing
  bool_default_map["write_hillshade"] = false;
  bool_default_map["print_raster_without_seas"] = false;
  bool_default_map["print_fill_raster"] = false;
  bool_default_map["print_channels_to_csv"] = false;
  
  // This converts all csv files to geojson (for easier loading in a GIS)
  bool_default_map["convert_csv_to_geojson"] = false;  

  // These are parameters for the cosmogenic data
  string_default_map["cosmo_parameter_prefix"] = "NULL";  
  bool_default_map["check_cosmo_basins"] = false;
  bool_default_map["spawn_cosmo_basins"] = false;
  int_default_map["spawn_padding_pixels"] = 20;

  // Here are some options for shielding
  bool_default_map["make_shielding_rasters"] = false;
  bool_default_map["calculate_shielding"] = false;

  // This is all for snow shielding
  bool_default_map["calculate_snow_shielding"] = false;
  string_default_map["snowpack_method"] = "Bilinear";

  // parameters for bilinear snowpack
  float_default_map["snow_SlopeAscend"] = 0.035;
  float_default_map["snow_SlopeDescend"] = -0.03;
  float_default_map["snow_PeakElevation"] = 1500;
  float_default_map["snow_PeakSnowpack"] = 30;

  // parameters for richards snowpack
  float_default_map["snow_v"] = 0.5;
  float_default_map["snow_lambda"] = 0.5;
  float_default_map["snow_MaximumSlope"] = 0.05;

  // These are for the erosion rate calculations
  bool_default_map["calculate_erosion_rates"] = false;
  bool_default_map["calculate_soil_erosion_rates"] = false;
  bool_default_map["calculate_nested_erosion_rates"] = false;
  
  bool_default_map["use_spawned_basins"] = false;
  bool_default_map["calculate_using_muons_and_errors"] = true;   

  bool_default_map["print_production_raster"] = false; 
  bool_default_map["print_pressure_raster"] = false; 

  bool_default_map["print_scaling_and_shielding_rasters"] = false;
  
  
  // The stuff below here is all for creating rasters for forward prediction of 
  // erosion rates.
  // Note the effective erosion rate is in g/cm^2/yr
  // To convert to mm/kyr you multiply by 10^7/(density in kg/m^3)
  // 0.01 g/cm^2/yr is ~ 37 mm/kyr
  // 0.1 g/cm^2/yr is ~ 370 mm/kyr or 0.37 mm/yr
  //bool_default_map["make_all_shielding_and_scaling_rasters"] = false;
  bool_default_map["single_erosion_rate"] = false;
  float_default_map["effective_erosion_rate"] = 0.01;
  float_default_map["self_shield_eff_thickness"] = 0;
  float_default_map["snow_shield_eff_thickness"] = 0;
  bool_default_map["calculate_CRN_concentration_raster"] = false;
  bool_default_map["load_production_raster"] = false;
  string_default_map["production_raster_suffix"] = "_PROD";
  string_default_map["concentration_column_name"] = "Be10_CONC";
  bool_default_map["calculate_erosion_rates_new"] = true;

  bool_default_map["calculate_accumulated_CRN_concentration"] = false;
  bool_default_map["calculate_accumulated_CRN_concentration_from_points"] = false;

  // some parameters for getting points in the landscape
  // You need a river network for this. 
  // You probably should keep the contributing pixels at the default unless you have either
  // very small or very large basins. Increase the number 
  bool_default_map["read_points_csv"] = false;
  string_default_map["points_filename"] ="CRN_points.csv";
  int_default_map["threshold_contributing_pixels"] = 1000;


  // These parameters are for raster aggregation that is used for sediment routine and 
  // CRN concentration prediction
  bool_default_map["route_cosmo_concentrations"] = false;
  string_default_map["raster_fnames_prefix"] = "NULL";
  bool_default_map["accumulate_cosmo"] = false;

  bool_default_map["check_sediment_routing_rasters"] = false;



  // This burns a raster value to any csv output of chi data
  // Useful for appending geology data to chi profiles
  bool_default_map["burn_raster_to_csv"] = false;
  string_default_map["burn_raster_prefix"] = "NULL";
  string_default_map["burn_data_csv_column_header"] = "burned_data";
  string_default_map["csv_to_burn_name"] = "NULL";

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

  cout << "Read filename is: " <<  DATA_DIR+DEM_ID << endl;
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;

  // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);



  //============================================================================
  // load the  DEM
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
  // Snow shielding
  // This makes a snow shielding raster
  //============================================================================
  if (this_bool_map["calculate_snow_shielding"])
  {
    cout << "Let me make a snow shielding raster for you." << endl;
    
    // now make a snow raster
    LSDSoilHydroRaster SnowRaster(topography_raster);
    
    string snow_ext = "_SnowBL";
    string snow_out_name = OUT_DIR+OUT_ID+snow_ext;
    
    // update it with a bilinear snow function
    if (this_string_map["snowpack_method"] == "Bilinear")
    {
      SnowRaster.SetSnowEffDepthBilinear(this_float_map["snow_SlopeAscend"], this_float_map["snow_SlopeDescend"], 
                                         this_float_map["snow_PeakElevation"], 
                                         this_float_map["snow_PeakSnowpack"], topography_raster);
                                  
      SnowRaster.write_raster(snow_out_name,raster_ext);
    }   
    else if (this_string_map["snowpack_method"] == "Richards")
    {
      SnowRaster.SetSnowEffDepthRichards(this_float_map["snow_PeakSnowpack"], this_float_map["snow_MaximumSlope"], 
                                         this_float_map["snow_v"], this_float_map["snow_lambda"],
                                         topography_raster);
      SnowRaster.write_raster(snow_out_name,raster_ext);
    }
    else
    {
      cout << "You did not give me a valid snow shielding method. " << endl;
      cout << "Options are Richards or Bilinear" << endl;
    } 
  }  

  //============================================================================
  // Raster burning
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
  // This begins the cosmogenic routines
  //============================================================================  
  // load the CRNCosmoData object
  cout << "===============================" << endl;
  if( this_bool_map["check_cosmo_basins"] ||
      this_bool_map["spawn_cosmo_basins"] ||
      this_bool_map["make_shielding_rasters"] ||
      this_bool_map["calculate_nested_erosion_rates"] ||
      this_bool_map["calculate_soil_erosion_rates"] ||
      this_bool_map["calculate_erosion_rates"])
  {
    cout << "You have asked for an analysis that requires cosmogenic data input." << endl;
    cout << "I am now going to try to load cosmogenic data." << endl;
    if (this_string_map["cosmo_parameter_prefix"] == "NULL")
    {
      cout << "You have not specified a cosmogenic parameter prefix. " << endl;
      cout << "I'm afraid I can't do any cosmogenic computations without this information." << endl;
      cout << "Make sure you designate the cosmo_parameter_prefix in the parameter file." << endl;
    }
    else
    {
      LSDCosmoData CosmoData(DATA_DIR,this_string_map["cosmo_parameter_prefix"]);
      cout << "Got the CosmoData, I am going to write the basin data and stream order rasters now." << endl; 
  
      // This checks the cosmo data by printing out the basins based on cosmo data locations
      if(this_bool_map["check_cosmo_basins"])
      {
        // check the basins
        CosmoData.print_basins_to_for_checking();
      }

      // This starts the spawning routine. It finds the basins from the cosmo data and 
      // then creates seperate rasters for computing each of these basins
      // This is only really required when you are using topographic shielding, 
      // which is highly computationally expensive. 
      // A recent ESURF paper by Roman DiBiase suggests these computations are not 
      // necessary so if you believe that paper you don't really need to do this spawning
      if(this_bool_map["spawn_cosmo_basins"])
      {
        cout << "I am spawning rasters. " << endl;
        cout << "Please note that if you want to use these rasters you need to run this program again" << endl;
        cout << "Using the updated parameter file produced by this routine that has _spawned in its name." << endl;
        CosmoData.BasinSpawnerMaster(DATA_DIR,this_string_map["cosmo_parameter_prefix"],this_int_map["spawn_padding_pixels"]);
      }

      // This is for making shielding rasters
      if(this_bool_map["make_shielding_rasters"])
      {
        cout << "I am going to make some shielding rasters for you." << endl;
        // Again note that if you belive DiBiase 2018 ESURF shielding calculations are not necessary.
        if(this_bool_map["calculate_shielding"])
        {
          cout << "I will calculate shielding for you." << endl;
          cout << "Note that if you want spawned basins you will need to run the spawning routine first." << endl;
          cout << "use  spawn_cosmo_basins: true" << endl;
          cout << "to activate spawning." << endl;
          CosmoData.RunShielding(DATA_DIR,this_string_map["cosmo_parameter_prefix"]);
        }      
        else
        {
          cout << "I am going to set all topographic shielding to 1." << endl;
          cout << "This is based on a recent Roman DiBiase ESURF paper that argued shielding doesn't matter." << endl;
          cout << "Note that if you want spawned basins you will need to run the spawning routine first." << endl;
          cout << "use  spawn_cosmo_basins: true" << endl;
          cout << "to activate spawning." << endl;
          CosmoData.RunShielding_Unshielded(DATA_DIR,this_string_map["cosmo_parameter_prefix"]);
        }
      }
       

      if(this_bool_map["calculate_nested_erosion_rates"])
      {
        cout << "I am now going to calculate nested erosion rates for you." << endl;
        cout << "This starts with the smallest nested basin and then" << endl;
        cout << "progressively calculates the erosion rates from the proportions" << endl;
        cout << "of the basin remaining. " << endl;
        CosmoData.calculate_nested_erosion_rates();
        CosmoData.print_results();
      } 

      // This calculates erosion rates based on discrete locations in the landscape.
      // That is, it does not attempt to integrate over contributing pixels, but just
      // calculates erosion rates based on local shielding and scaling. 
      if(this_bool_map["calculate_soil_erosion_rates"])
      {
        cout << "I am now going to calculate erosion rates from soil samples." << endl;
        cout << "This is a purely local calculation." << endl;
        int method_flag = 3;
        CosmoData.calculate_erosion_rates(method_flag);
        CosmoData.print_results();
      } 

      // This is the main erosion rate calculator. It takes locations in channels
      // and calculates, with uncertainties, the inferred erosion rates based on
      // 10Be or other isotope data. 
      if(this_bool_map["calculate_erosion_rates"]) 
      {
        // the method flag determines what sort of analysis you want to do
        int method_flag; 

        if(not this_bool_map["calculate_using_muons_and_errors"])
        {
          cout << "I am going to calculate basic erosion rates without" << endl;
          cout << "muogenic production or errors. This is just a very" << endl;
          cout << "fast way to get ballpark erosion rates." << endl;
          method_flag = 0;            
          CosmoData.calculate_erosion_rates(method_flag);
        }
        else
        {
          if(this_bool_map["use_spawned_basins"])
          {
            cout << "I am calculating erosion rates using spawning. " << endl;
            cout << "If you haven't run the spawning routine yet this won't work!" << endl;
            method_flag = 2;
            CosmoData.calculate_erosion_rates(method_flag);
          }
          else
          {
            cout << "I am calculating erosion rates without spawning. " << endl;
            cout << "If you have a big DEM and haven't done shielding this will take ages!" << endl;
            method_flag = 1;
            CosmoData.calculate_erosion_rates(method_flag); 
          }
        }
        CosmoData.print_results();

      }        // End logic for the main erosion rate 

      // This prints the full scaling and shielding rasters
      if(this_bool_map["print_scaling_and_shielding_rasters"])
      {
        cout << "Let me print the scaling and shielding rasters for you." << endl;
        CosmoData.print_scaling_and_shielding_complete_rasters();
      }
    }          // End logic for routines where the cosmogenic data prefix has been found
  }            // End logic for routines that require checking of the cosmogenic data


  //-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // CRN PREDICTION
  //
  // The following routines are for making predictions about CRN concentrations
  // under a number of different scenarios
  //
  //-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if(this_bool_map["print_production_raster"])
  {
    cout << "I am going to print the production raster for you." << endl;

    // The atmospheric data needs to be in the same directory as the directory
    // from which the programis called.
    // This should be modified later!!
    string path_to_atmospheric_data = "./";

    // We need this to make the production raster
    LSDCosmoRaster ThisCosmoRaster(topography_raster);  

    LSDRaster production = ThisCosmoRaster.calculate_production_raster(topography_raster, path_to_atmospheric_data);

    string this_raster_name = OUT_DIR+OUT_ID+"_PROD";                                      
    production.write_raster(this_raster_name,raster_ext);

  }

  if(this_bool_map["print_pressure_raster"])
  {
    cout << "I am going to print the pressure raster for you." << endl;
    
    // The atmospheric data needs to be in the same directory as the directory
    // from which the programis called.
    // This should be modified later!!
    string path_to_atmospheric_data = "./";

    // We need this to make the production raster
    LSDCosmoRaster ThisCosmoRaster(topography_raster);  

    LSDRaster pressure = ThisCosmoRaster.calculate_atm_pressure_raster(topography_raster, path_to_atmospheric_data);

    string this_raster_name = OUT_DIR+OUT_ID+"_ATMPRESSURE";                                      
    pressure.write_raster(this_raster_name,raster_ext);

  }
 
  // This makes all the appropriate rasters for constant values of the 
  // shielding and erosion. It prints the concentration raster,
  // and the erosion raster
  // The concentration is the concentration at individual pixels,
  // not accumulated concentration. 
  // For areas with landslides (e.g., nonzero pixel values in the self shielding raster)
  // the concentration is the average concentration of the evacuated material. 
  if(this_bool_map["calculate_CRN_concentration_raster"])
  {
    cout << "I will now calculate the concentration of your nuclide at each pixel in your raster." << endl;
    cout << "For pixels with landslides, this is the mean concentration of the eroded material. " << endl;

    // The atmospheric data needs to be in the same directory as the directory
    // from which the program is called.
    // This should be modified later!!
    string path_to_atmospheric_data = "./";
    
    // This will be used to make some secondary rasters
    LSDRasterMaker MakeItYeah(topography_raster);
    
    // We need this to make the production raster
    LSDCosmoRaster ThisCosmoRaster(topography_raster);
    
    cout << "Let me make the production raster." << endl;
    LSDRaster ProductionRaster;
    if(this_bool_map["load_production_raster"])
    {
      cout << "Let me load a production raster. " << endl;
      cout << "At the moment there is no error checking so if you don't have this raster the program will fail. " << endl;
      string Prod_fname = DATA_DIR+DEM_ID+this_string_map["production_raster_suffix"];
      cout << "The production raster filename is: " << Prod_fname << endl;
      LSDRasterInfo PRI(Prod_fname,raster_ext);
      LSDRaster LoadProductionRaster(DATA_DIR+DEM_ID+this_string_map["production_raster_suffix"], raster_ext);
      ProductionRaster = LoadProductionRaster;
    }
    else
    {
      LSDRaster LoadProductionRaster = ThisCosmoRaster.calculate_production_raster(topography_raster,
                                                        path_to_atmospheric_data);
      ProductionRaster = LoadProductionRaster;                                           
    }
    cout << "I've got the production raster." << endl;
    
    // Now make the shielding and other rasters
    // TODO
    // This needs to be make more flexible in the future to allow
    // loading of rasters
    MakeItYeah.set_to_constant_value(this_float_map["effective_erosion_rate"]);
    LSDRaster ErosionRaster =  MakeItYeah.return_as_raster();
    
    MakeItYeah.set_to_constant_value(this_float_map["self_shield_eff_thickness"]);
    LSDRaster SelfShield =  MakeItYeah.return_as_raster();    
                                                        
    MakeItYeah.set_to_constant_value(this_float_map["snow_shield_eff_thickness"]);
    LSDRaster SnowShield =  MakeItYeah.return_as_raster();   
    
    MakeItYeah.set_to_constant_value(1.0);
    LSDRaster TopoShield =  MakeItYeah.return_as_raster();      
    
    // TODO
    // Later this will need to be updated to get the nuclide on other factors from input to the program
    string Nuclide = "Be10";
    string Muon_scaling = "newCRONUS";
    
    bool is_production_uncertainty_plus_on = false;
    bool is_production_uncertainty_minus_on = false;
    
    LSDRaster CRNConc = ThisCosmoRaster.calculate_CRN_concentration_raster(Nuclide, Muon_scaling, ErosionRaster,
                                          ProductionRaster, TopoShield, 
                                          SelfShield, SnowShield, 
                                          is_production_uncertainty_plus_on,
                                          is_production_uncertainty_minus_on);

    // This writes the concentration raster 
    string this_raster_name = OUT_DIR+OUT_ID+"_"+Nuclide+"_Conc";                                      
    CRNConc.write_raster(this_raster_name,raster_ext);
    
    // This writes the erosion raster
    this_raster_name = OUT_DIR+OUT_ID+"_EffEros";                                      
    ErosionRaster.write_raster(this_raster_name,raster_ext);
  }




  //===================================================================================================
  //===================================================================================================
  // This accumulates CRN. It wraps the concentration calculation
  if( this_bool_map["calculate_accumulated_CRN_concentration"] ||
      this_bool_map["calculate_accumulated_CRN_concentration_from_points"] ) 
  {
    // First some housekeeping: we need the flow info object for these routines
    // Get the filled topography
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

    // I am going to route some sediment concentrations for you so I need to do some 
    // flow routing
    cout << "\t Flow routing..." << endl;
    // get a flow info object
    LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

    // We need this to make the production raster
    LSDCosmoRaster ThisCosmoRaster(topography_raster);
    

    cout << "Let me make or get the production raster." << endl;
    LSDRaster ProductionRaster;
    
    // The atmospheric data needs to be in the same directory as the directory
    // from which the programis called.
    // This should be modified later!!
    string path_to_atmospheric_data = "./";
    if(this_bool_map["load_production_raster"])
    {
      cout << "Let me load a production raster. " << endl;
      cout << "If you don't have this raster the program will fail. " << endl;
      string Prod_fname = DATA_DIR+DEM_ID+this_string_map["production_raster_suffix"];
      cout << "The production raster filename is: " << Prod_fname;
      LSDRasterInfo PRI(Prod_fname,raster_ext);
      LSDRaster LoadProductionRaster(DATA_DIR+DEM_ID+this_string_map["production_raster_suffix"], raster_ext);
      ProductionRaster = LoadProductionRaster;
    }
    else
    {
      LSDRaster LoadProductionRaster = ThisCosmoRaster.calculate_production_raster(topography_raster,
                                                        path_to_atmospheric_data);
      ProductionRaster = LoadProductionRaster;                                           
    }
    cout << "I've got the production raster." << endl;

    // This will be used to make some secondary rasters
    LSDRasterMaker MakeItYeah(topography_raster);
    

    
    // Now make the shielding and other rasters
    // TODO
    // This will need to be more flexible and allow for reading of rasters
    MakeItYeah.set_to_constant_value(this_float_map["effective_erosion_rate"]);
    LSDRaster ErosionRaster =  MakeItYeah.return_as_raster();    
    MakeItYeah.set_to_constant_value(this_float_map["self_shield_eff_thickness"]);
    LSDRaster SelfShield =  MakeItYeah.return_as_raster();                                                            
    MakeItYeah.set_to_constant_value(this_float_map["snow_shield_eff_thickness"]);
    LSDRaster SnowShield =  MakeItYeah.return_as_raster();       
    MakeItYeah.set_to_constant_value(1.0);
    LSDRaster TopoShield =  MakeItYeah.return_as_raster();  
    MakeItYeah.set_to_constant_value(TopoShield.get_NoDataValue());
    LSDRaster NoDataRaster =  MakeItYeah.return_as_raster();  

    // These are used to calculate uncertainties so are only switched on when you 
    // are looking for errors. 
    bool is_production_uncertainty_plus_on = false;
    bool is_production_uncertainty_minus_on = false;

    // TODO
    // Later you need to read this from file
    string Nuclide = "Be10";
    string Muon_scaling = "newCRONUS";

    if(this_bool_map["calculate_accumulated_CRN_concentration"])
    {
      cout << "I will now calculate the accumulated concentration of your nuclide." << endl;
      
      auto t1 = std::chrono::high_resolution_clock::now();   
      LSDRaster CRN_conc = ThisCosmoRaster.calculate_CRN_concentration_raster(Nuclide, Muon_scaling, ErosionRaster,
                                            ProductionRaster, TopoShield, 
                                            SelfShield, SnowShield, 
                                            is_production_uncertainty_plus_on,
                                            is_production_uncertainty_minus_on);
      auto t2 = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
      cout << "Full concentration calculation took: " << duration << endl;

      // This writes the concentration raster 
      string this_raster_name = OUT_DIR+OUT_ID+"_"+Nuclide+"_Conc";                                   
      CRN_conc.write_raster(this_raster_name,raster_ext);
      
      // This writes the erosion raster
      this_raster_name = OUT_DIR+OUT_ID+"_EffEros";                                      
      ErosionRaster.write_raster(this_raster_name,raster_ext);

      // Some testing
      // This will be used to make some secondary rasters
      //MakeItYeah.set_to_constant_value(0.5);
      //LSDRaster AreaPixelRaster =  MakeItYeah.return_as_raster();

      // Need to test this function
      //LSDIndexRaster IntContrib = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
      //LSDRaster Contrib(IntContrib);

      // PLACEHOLDER
      //float minimum_val = 0.001;
      //float maximum_val = 0.005;
      //MakeItYeah.random_values(minimum_val,maximum_val);
      //LSDRaster Random_erosion =  MakeItYeah.return_as_raster();  

      cout << endl << endl << "=================================" << endl;
      cout << "This is currently in testing phase. I am trying to accumulate the cosogenic concentration." << endl;
      LSDRaster AccConc = ThisCosmoRaster.calculate_accumulated_CRN_concentration(CRN_conc, ErosionRaster,FlowInfo);   
      cout << "Done with accumulating the concentration." << endl;

      // This writes the accumulation raster
      this_raster_name = OUT_DIR+OUT_ID+"_AccConc";                                      
      AccConc.write_raster(this_raster_name,raster_ext);
    }


    if(this_bool_map["calculate_accumulated_CRN_concentration_from_points"] ||
       this_bool_map["calculate_erosion_rates"])
    {
      cout << "I am reading points from the file: "+ this_string_map["points_filename"] << endl;
      LSDSpatialCSVReader CRN_points_data( RI, (DATA_DIR+this_string_map["points_filename"]) );

      // You need to get the node indices of the points first
      // get some relevant rasters
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
      LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      //get the sources: note: this is only to select basins!
      vector<int> sources;
      sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["threshold_contributing_pixels"]);

      // now get the junction network
      LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

      // Print channels and junctions if you want them.
      if( this_bool_map["print_channels_to_csv"])
      {
        cout << "I am going to print the channel network." << endl;
        string channel_csv_name = OUT_DIR+OUT_ID+"_CN";
        ChanNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);
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

      // Now get the lat long points
      // First get the lat-long
      int UTM_zone;
      bool is_North;
      CRN_points_data.get_UTM_information(UTM_zone, is_North);
      cout << "The UTM zone is: " << UTM_zone << " and it ";
      if(is_North)
      {
        cout << "is north." << endl;
      }
      else
      {
        cout << "is south." << endl;
      }

      // Get the local coordinates
      vector<float> fUTM_easting,fUTM_northing;
      CRN_points_data.get_x_and_y_from_latlong(fUTM_easting,fUTM_northing);
      cout << "The x and y data are: " << endl;
      for (int i = 0; i< int(fUTM_easting.size()); i++)
      {
        cout << fUTM_easting[i] << "," << fUTM_northing[i] << endl;
      }
      
      int search_radius_nodes = 8;
      int threshold_stream_order = 3;
      vector<int> valid_cosmo_points;
      vector<int> snapped_node_indices;
      vector<int> snapped_junction_indices;

      // Now get the snapped points 
      // The snapped points are in a vector called snapped_node_indices  
      cout << "I'm snapping some points" << endl; 
      ChanNetwork.snap_point_locations_to_nearest_channel_node_index(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices);
      cout << "The number of snapped points is: " <<   valid_cosmo_points.size() << endl;   

      // Now get accumulated cosmo for testing
      auto t1 = std::chrono::high_resolution_clock::now();   
      LSDRaster CRN_conc = ThisCosmoRaster.calculate_CRN_concentration_raster(Nuclide, Muon_scaling, ErosionRaster,
                                            ProductionRaster, TopoShield, 
                                            SelfShield, SnowShield, 
                                            is_production_uncertainty_plus_on,
                                            is_production_uncertainty_minus_on);
      auto t2 = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
      cout << "Full concentration calculation took: " << duration << endl;

      t1 = std::chrono::high_resolution_clock::now();   
      LSDRaster AccConc = ThisCosmoRaster.calculate_accumulated_CRN_concentration(CRN_conc, ErosionRaster,FlowInfo);   
      t2 = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
      cout << "Accumulation calculation took: " << duration << endl;       


      if(this_bool_map["calculate_accumulated_CRN_concentration_from_points"])
      {

        // now loop through the valid node indices:
        int n_nodes = int(snapped_node_indices.size());
        for(int node = 0; node<n_nodes; node++)
        {
          cout << "Found the node " << snapped_node_indices[node] << endl;

          t1 = std::chrono::high_resolution_clock::now();   
          LSDRaster CRNConc_node = ThisCosmoRaster.calculate_CRN_concentration_raster(Nuclide, Muon_scaling, ErosionRaster,
                                            ProductionRaster, TopoShield, 
                                            SelfShield, SnowShield, 
                                            FlowInfo,snapped_node_indices[node],
                                            is_production_uncertainty_plus_on,
                                            is_production_uncertainty_minus_on);

          t2 = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
          cout << "Basin concentration calculation took: " << duration << endl;

          // This writes the concentration raster 
          string this_raster_name = OUT_DIR+OUT_ID+"_"+Nuclide+"_Conc_"+itoa(node);                                   
          CRNConc_node.write_raster(this_raster_name,raster_ext);

          // Now accumulate from this node
          float This_accumulated_cosmo = ThisCosmoRaster.calculate_accumulated_CRN_concentration(CRNConc_node, 
                                            ErosionRaster, FlowInfo, 
                                            snapped_node_indices[node]);

          cout << endl << endl << "===========================" << endl;
          cout << "Accumulated concentration is: " << This_accumulated_cosmo << " atoms/g" <<endl;

          int row,col;
          FlowInfo.retrieve_current_row_and_col(snapped_node_indices[node],row,col);
          cout << "And same accumulation from the raster I previously calculated: " <<  AccConc.get_data_element(row,col);
          cout << "========================" << endl << endl;

        }
      }

      if( this_bool_map["calculate_erosion_rates_new"])
      {
        cout << "Let me calculate the erosion rate for you. "  << endl;

        // This looks for concentrations data
        cout << "The concentration column name is: " << this_string_map["concentration_column_name"] << endl;
        vector<float> concentrations = CRN_points_data.data_column_to_float(this_string_map["concentration_column_name"]);

        // Now loop through concentrations, getting the erosion rates
        int n_valid_samples = int(valid_cosmo_points.size());
        for(int sample = 0; sample< n_valid_samples; sample++)
        {
          cout << "The concentration is: " << concentrations[ valid_cosmo_points[sample]  ] << endl;

          float This_erate = ThisCosmoRaster.calculate_eff_erate_from_conc(concentrations[ valid_cosmo_points[sample] ],
                                            Nuclide, Muon_scaling, NoDataRaster,
                                            ProductionRaster, TopoShield, 
                                            SelfShield, SnowShield, 
                                            NoDataRaster,FlowInfo,
                                            snapped_node_indices[sample]);

          cout << "The erosion rate is: " << This_erate << endl;
        }
      }

    } // end logic for getting points
  } // end logic for getting flow info


  //================================================================================================
  // This checks the sediment routing components
  // it is essentially used to both debug the raster aggregator as well
  // as test of the raster file is working properly
  // AT THE MOMENT THIS DOES NOTHING
  if(this_bool_map["check_sediment_routing_rasters"])
  {
    cout << "Let me check some sediment routing for you!" << endl;
    // Make a raster aggregator object
    LSDRasterAggregator LSDRA(DATA_DIR,this_string_map["raster_fnames_prefix"]);
    
    // Print the rasters and keys to screen.
    LSDRA.print_raster_names_and_types_to_screen();

    // See if all the keys are represented
    vector<string> required_rasters;
    required_rasters.push_back("DEM");
    required_rasters.push_back("10BeCONC");
    required_rasters.push_back("Funky CHICKEN");
    LSDRA.check_raster_types(required_rasters);
  }




  // This is the logic for sediment routing
  if(this_bool_map["route_cosmo_concentrations"])
  {

    // Get the filled topography
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
      cout << "Let me print the fill raster for you."  << endl;
      string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
      filled_topography.write_raster(filled_raster_name,raster_ext);
    }

    // I am going to route some sediment concentrations for you so I need to do some 
    // flow routing
    cout << "\t Flow routing..." << endl;
    // get a flow info object
    LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

    // Now create a raster aggregator
    LSDRasterAggregator LSDRA();
  }
  
  bool_default_map["route_cosmo_concentrations"] = false;
  string_default_map["raster_fnames_prefix"] = "NULL";
  bool_default_map["accumulate_cosmo"] = false;  
  

  cout << "I'm all finished! Have a nice day." << endl;
}
