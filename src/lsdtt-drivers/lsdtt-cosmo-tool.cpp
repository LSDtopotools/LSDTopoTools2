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
// Call with -h to generate a help file
//
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

  string version_number = "0.7";
  string citation = "http://doi.org/10.5281/zenodo.4577879";

  cout << "=========================================================" << endl;
  cout << "|| Welcome to the LSDTopoTools cosmogenic tool!        ||" << endl;
  cout << "|| This program has a number of options for calculating||" << endl;
  cout << "|| cosmogenic nuclide concentrations.                  ||" << endl;
  cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
  cout << "||  at the University of Edinburgh                     ||" << endl;
  cout << "=========================================================" << endl;  
  cout << "|| Citation for this code is:                          ||" << endl;
  cout << "|| " << citation << endl; 
  cout << "|| If you use these routines please cite:              ||" << endl;   
  cout << "|| https://doi.org/10.5194/esurf-4-655-2016            ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| Documentation can be found at:                      ||" << endl;
  cout << "|| https://lsdtopotools.github.io/LSDTT_documentation/ ||" << endl;
  cout << "=========================================================" << endl;
  cout << "|| This is LSDTopoTools2 version                       ||" << endl;
  cout << "|| " << version_number << endl;
  cout << "|| If the version number has a d at the end it is a    ||" << endl;
  cout << "||  development version.                               ||" << endl;
  cout << "=========================================================" << endl;

  // this will contain the help file
  map< string, vector<string> > help_map;

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
    ofs.open("./lsdtt-cosmo-tool-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    cout << endl << endl << endl << "==============================================" << endl;    
    cout << "This is lsdtt-cosmo-tool version number " << version_number << endl;
    cout << "If the version contains a 'd' then you are using a development version." << endl;
    cout << "=========================================================" << endl;
    ofstream ofs;
    ofs.open("./lsdtt-cosmo-tool-version.txt");
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
 
  bool_default_map["only_check_parameters"] = false;
  help_map["only_check_parameters"] = {  "bool","false","This just checks parameters without running an analysis.","For bug checking."};


  // the most basic raster printing
  bool_default_map["write_hillshade"] = false;
  help_map["write_hillshade"] = {  "bool","false","Write the hillshade raster.","You need this for a lot of our plotting routines. Filename includes _HS"};

  bool_default_map["print_raster_without_seas"] = false;
  help_map["print_raster_without_seas"] = {  "bool","false","Overwrites the raster without seas.","DANGER this will replace your existing raster with the high and low points replaced by nodata. See the remove_seas flag"};

  bool_default_map["print_fill_raster"] = false;
  help_map["print_fill_raster"] = {  "bool","false","Prints the fill raster.","Filename includes _FILL"};

  // This converts all csv files to geojson (for easier loading in a GIS)
  bool_default_map["convert_csv_to_geojson"] = false;
  help_map["convert_csv_to_geojson"] = {  "bool","false","Converts csv files to geojson files","Makes csv output easier to read with a GIS. Warning: these files are much bigger than csv files."};

  bool_default_map["print_channels_to_csv"] = false;
  help_map["print_channels_to_csv"] = {  "bool","false","Prints the channel network to a csv file.","This version produces smaller files than the raster version."};


  // These are parameters for the cosmogenic data
  string_default_map["cosmo_parameter_prefix"] = "NULL";  
  help_map["cosmo_parameter_prefix"] = {  "string","NULL","The prefix of the cosmogenic parameter files.","You will need three files with this prefix and extensions .CRNParam ._CRNData.csv and ._CRNRasters.csv."};

  bool_default_map["check_cosmo_basins"] = false;
  help_map["check_cosmo_basins"] = {  "bool","false","This checks if the cosmo data are in the dem and prints a basin raster.","Used to check input data locations."};

  bool_default_map["spawn_cosmo_basins"] = false;
  help_map["spawn_cosmo_basins"] = {  "bool","false","If true will create a little raster for each basin.","Used to speed up computation."};

  int_default_map["spawn_padding_pixels"] = 20;
  help_map["spawn_padding_pixels"] = {  "int","20","How many padded ext5ra pixels around each basin after spawning.","Use with spawn_cosmo_basins"};


  // Here are some options for shielding
  bool_default_map["make_shielding_rasters"] = false;
  help_map["make_shielding_rasters"] = {  "bool","false","If true calculate all the shielding rasters.","Makes shielding in a separate step from erosion rate analysis"};

  bool_default_map["calculate_shielding"] = false;
  help_map["calculate_shielding"] = {  "bool","false","You need to make all the shielding rasters with make_shielding_rasters but if this is true it will calculate full topographic shielding.","DiBiase 2018 ESURF suggests topographic shielding calculations are not needed so keep this false to save computational time"};


  // This is all for snow shielding
  bool_default_map["calculate_snow_shielding"] = false;
  help_map["calculate_snow_shielding"] = {  "bool","false","If true this calculates the snow shielding using some snowpack information.","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  string_default_map["snowpack_method"] = "Bilinear";
  help_map["snowpack_method"] = {  "string","Bilinear","Method for calculating snowpack. Options are Bilinear and Richards. Case sensitive","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};


  // parameters for bilinear snowpack
  float_default_map["snow_SlopeAscend"] = 0.035;
  help_map["snow_SlopeAscend"] = {  "float","0.035","A parameter for the bilinear snowpack method","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  float_default_map["snow_SlopeDescend"] = -0.03;
  help_map["snow_SlopeDescend"] = {  "float","-0.03","A parameter for the bilinear snowpack method","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  float_default_map["snow_PeakElevation"] = 1500;
  help_map["snow_PeakElevation"] = {  "float","1500","A parameter for the bilinear snowpack method. The elevation of the peak snowpack","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  float_default_map["snow_PeakSnowpack"] = 30;
  help_map["snow_PeakSnowpack"] = {  "float","30","A parameter for the bilinear snowpack method. The peak snowpack in g/cm^2","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  // parameters for richards snowpack
  float_default_map["snow_v"] = 0.5;
  help_map["snow_v"] = {  "float","0.5", "A parameter for the Richards snowpack method","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  float_default_map["snow_lambda"] = 0.5;
  help_map["snow_lambda"] = {  "float","0.5", "A parameter for the Richards snowpack method","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};

  float_default_map["snow_MaximumSlope"] = 0.05;
  help_map["snow_MaximumSlope"] = {  "float","0.05", "A parameter for the Richards snowpack method","You need to fit all the snowpack parameters using some python code and then adjust the parameters for the snowpack in the parameter file"};


  //..####....####...######..#####...##..##...........####...##...##..######..######...####...##..##..######...####..
  //.##..##..##..##....##....##..##..###.##..........##......##...##....##......##....##..##..##..##..##......##.....
  //.##......######....##....#####...##.###...........####...##.#.##....##......##....##......######..####.....####..
  //.##..##..##..##....##....##..##..##..##..............##..#######....##......##....##..##..##..##..##..........##.
  //..####...##..##..######..##..##..##..##...........####....##.##...######....##.....####...##..##..######...####..
  //
  // These are for the erosion rate calculations
  bool_default_map["calculate_erosion_rates"] = false;
  help_map["calculate_erosion_rates"] = {  "bool","false", "The basic CAIRN routine that just gets the erosion rate for each CRN point","See Mudd et al 2016 ESURF for details"};

  bool_default_map["calculate_soil_erosion_rates"] = false;
  help_map["calculate_soil_erosion_rates"] = {  "bool","false", "This calculates soil erosion rates so does not accumulate erosion from upslope","See Mudd et al 2016 ESURF for details"};

  bool_default_map["calculate_nested_erosion_rates"] = false;
  help_map["calculate_nested_erosion_rates"] = {  "bool","false", "This looks for nested basins. It calculates the smallest basin first and then sets this as the erosion rate for pixels in that basin before it goes on to calculate the next largest basin.","See Mudd et al 2016 ESURF for details"};

  bool_default_map["use_spawned_basins"] = false;
  help_map["use_spawned_basins"] = {  "bool","false", "You need to have previously set spawn_cosmo_basins to true and this will use those for erosion rate calculations. Spawning saves computational expense but is not necessary if you have turned of the topographic shielding calculation as advised by DiBiase ESURF 2018.","See Mudd et al 2016 ESURF for details"};

  bool_default_map["calculate_using_muons_and_errors"] = true;   
  help_map["calculate_using_muons_and_errors"] = {  "bool","true", "If you want a very simple calculation just using spallation set this to false.","See Mudd et al 2016 ESURF for details"};

  bool_default_map["print_production_raster"] = false; 
  help_map["print_production_raster"] = {  "bool","false", "Prints the production raster duh.","You can recycle the production raster in other analyses such as forward modelling of 10Be concentrations under different erosion scenarios"};

  bool_default_map["print_pressure_raster"] = false; 
  help_map["print_pressure_raster"] = {  "bool","false", "Prints the pressure raster for bug checking.","Data comes from NCEP see balco's cronus calculator and this can be used to check if the numbers make sense."};

  bool_default_map["print_scaling_and_shielding_rasters"] = false;
  help_map["print_scaling_and_shielding_rasters"] = {  "bool","false", "Prints scaling and shielding rasters.","You can recycle the production raster in other analyses such as forward modelling of 10Be concentrations under different erosion scenarios"};


  // ..####....####...##......##..##..##...##..##..##..........######..######...####...######.
  // .##..##..##..##..##......##..##..###.###..###.##............##....##......##........##...
  // .##......##..##..##......##..##..##.#.##..##.###............##....####.....####.....##...
  // .##..##..##..##..##......##..##..##...##..##..##............##....##..........##....##...
  // ..####....####...######...####...##...##..##..##............##....######...####.....##...
  //     
  // The stuff below here is for testing columns of particles
  bool_default_map["transient_column_calculator"] = false;
  help_map["transient_column_calculator"] = {  "bool","false", "A tool that allows calculation of multiple nuclides in a column of rock.","Used for testing some of the transient concentration scenarios."};

  string_default_map["column_mode"] = "void";
  help_map["column_mode"] = {  "string","void", "Sets the method of the transient column.","Options are step_change:full_transient:steady_state"};

  float_default_map["effective_erosion_rate_new"] = 0.01;
  help_map["effective_erosion_rate_new"] = {  "float","0.01", "Used for forward modelling of cosmogenic concentrations in the transient mode. This is the rate after the step change. Units are g/cm^2/yr.","To convert to mm/kyr you multiply by 10^7/(density in kg/m^3)  0.01 g/cm^2/yr is ~ 37 mm/kyr  0.1 g/cm^2/yr is ~ 370 mm/kyr or 0.37 mm/yr"};
   
  float_default_map["time_since_step"] = 0;
  help_map["time_since_step"] = {  "float","0", "Used for testing of the particle column. Time since step change.","Can also be used to make a step change raster"};
 
  float_default_map["total_shielding_test"] = 1;
  help_map["total_shielding_test"] = {  "float","1", "The total shielding value for testing the CRN concentrations of particles.","This is used so you don't have to test shielding."};

  float_default_map["transient_effective_depth_test_top"] = 0;
  help_map["transient_effective_depth_test_top"] = {  "float","0", "For testing particle concentrations. This is the top depth for depth integrated calculations.","In g/cm^2."};

  float_default_map["transient_effective_depth_test_bottom"] = 0;
  help_map["transient_effective_depth_test_bottom"] = {  "float","0", "For testing particle concentrations. This is the bottom depth for depth integrated calculations.","In g/cm^2."};

  float_default_map["transient_effective_depth_test"] = 0;
  help_map["transient_effective_depth_test"] = {  "float","0", "For testing the particle concentrations. This is the effective depth for the test.","In g/cm^2."};


  // Some names for shielding and production rasters
  string_default_map["self_shielding_raster_prefix"] = "NULL";
  help_map["self_shielding_raster_prefix"] = {  "string","NULL", "Self shielding effective depth raster prefix. Used for forward modelling of cosmogenic concentrations and the new erosion calculations. Units are g/cm^2/yr.","If this is set to NULL then the value will be taken from self_shield_eff_thickness"};
 
  string_default_map["snow_shielding_raster_prefix"] = "NULL";
  help_map["snow_shielding_raster_prefix"] = {  "string","NULL", "Snow shielding effective depth raster prefix. Used for forward modelling of cosmogenic concentrations and the new erosion calculations. Units are g/cm^2/yr.","If this is set to NULL then the value will be taken from snow_shield_eff_thickness"};
 
  string_default_map["topographic_shielding_raster_prefix"] = "NULL";
  help_map["topographic_shielding_raster_prefix"] = {  "string","NULL", "Topographic shielding effective depth raster prefix. Used for forward modelling of cosmogenic concentrations and the new erosion calculations. Varies from 0 to 1.","If value is NULL the topographic shielding factor will be 1 everywhere based on a paper by DiBiase"};

  string_default_map["background_erosion_rate_raster_prefix"] = "NULL";
  help_map["background_erosion_rate_raster_prefix"] = {  "string","NULL", "Prefix of raster that has a background erosion rate. Used for forward modelling of cosmogenic concentrations and the new erosion calculations. Varies from 0 to 1.","If value is NULL the topographic shielding factor will be 1 everywhere based on a paper by DiBiase"};
 
  string_default_map["quartz_content_raster_prefix"] = "NULL";
  help_map["quartz_content_raster_prefix"] = {  "string","NULL", "Quartz content (as mass fraction). Used for forward modelling of cosmogenic concentrations and the new erosion calculations. Units are g/cm^2/yr.","If value is NULL then the whole raster gets the value of effective_erosion_rate. To convert to mm/kyr you multiply by 10^7/(density in kg/m^3)  0.01 g/cm^2/yr is ~ 37 mm/kyr  0.1 g/cm^2/yr is ~ 370 mm/kyr or 0.37 mm/yr"};
 
  // Nuclide settings for prediction
  string_default_map["muon_scheme_for_prediction"] = "BraucherBorchers";
  help_map["muon_scheme_for_prediction"] = {  "string","BraucherBorchers", "The muon scaling scheme for predictive modelling of CRN concentrations","Options are newCRONUS:Schller:Granger:Borchers:BraucherBorchers"};
   
  string_default_map["nuclide_for_prediction"] = "Be10";
  help_map["nuclide_for_prediction"] = {  "string","Be10", "The nuclide to be used for prediction","Options are Be10:C14:Al26:Cl36"};
   


  //.######...####...#####...##...##...####...#####...#####...........#####....####...#####....####...##...##...####..
  //.##......##..##..##..##..##...##..##..##..##..##..##..##..........##..##..##..##..##..##..##..##..###.###..##.....
  //.####....##..##..#####...##.#.##..######..#####...##..##..........#####...######..#####...######..##.#.##...####..
  //.##......##..##..##..##..#######..##..##..##..##..##..##..........##......##..##..##..##..##..##..##...##......##.
  //.##.......####...##..##...##.##...##..##..##..##..#####...........##......##..##..##..##..##..##..##...##...####..
  //
  // The stuff below here is all for creating rasters for forward prediction of 
  // erosion rates.
  // Note the effective erosion rate is in g/cm^2/yr
  // To convert to mm/kyr you multiply by 10^7/(density in kg/m^3)
  // 0.01 g/cm^2/yr is ~ 37 mm/kyr
  // 0.1 g/cm^2/yr is ~ 370 mm/kyr or 0.37 mm/yr
  //bool_default_map["make_all_shielding_and_scaling_rasters"] = false;
  bool_default_map["calculate_CRN_concentration_raster"] = false;
  help_map["calculate_CRN_concentration_raster"] = {  "bool","false", "Switch to true to forward model CRN concentrations.","This can be used to test the 10Be outcome of different erosion scenarios"};

  float_default_map["effective_erosion_rate"] = 0.01;
  help_map["effective_erosion_rate"] = {  "float","0.01", "Used for forward modelling of cosmogenic concentrations. Units are g/cm^2/yr.","To convert to mm/kyr you multiply by 10^7/(density in kg/m^3)  0.01 g/cm^2/yr is ~ 37 mm/kyr  0.1 g/cm^2/yr is ~ 370 mm/kyr or 0.37 mm/yr"};
 
  float_default_map["self_shield_eff_thickness"] = 0;
  help_map["self_shield_eff_thickness"] = {  "float","0", "Used for forward modelling of cosmogenic concentrations. Units are g/cm^2/yr.","This is for removing chunks of surface material"};
 
  float_default_map["snow_shield_eff_thickness"] = 0;
  help_map["snow_shield_eff_thickness"] = {  "float","0", "Single snow shielding used for forward modelling of cosmogenic concentrations. Units are g/cm^2/yr.","This is for removing chunks of surface material"};

  bool_default_map["load_production_raster"] = false;
  help_map["load_production_raster"] = {  "bool","false", "For forward modelling of CRN concentrations.","You can precalculate a production raster that then gets used for multiple 10Be concentration scenarios"};

  string_default_map["production_raster_suffix"] = "_PROD";
  help_map["production_raster_suffix"] = {  "string","_PROD", "If you are loading a production raster using load_production_raster this is the extension of the production raster.","The suffix goes after the DEM_ID"};

  bool_default_map["calculate_erosion_rates_new"] = false;
  help_map["calculate_erosion_rates_new"] = {  "bool","false", "For cosmogenic calculations this turns on the new more computationally efficient method of calculating denudation rates.","You still need to set calculate_erosion_rates to true"};

  string_default_map["concentration_column_name"] = "Be10_CONC";
  help_map["concentration_column_name"] = {  "string","Be10_CONC", "If you use the calculate_erosion_rates_new this is the column in the points_filename file that has the 10Be concentration.","Only used with calculate_erosion_rates_new"};

  bool_default_map["calculate_accumulated_CRN_concentration"] = false;
  help_map["calculate_accumulated_CRN_concentration"] = {  "bool","false", "Calculates accumulated CRN concentrations across a DEM.","For forward modelling of 10Be concentrations"};

  bool_default_map["calculate_accumulated_CRN_concentration_from_points"] = false;
  help_map["calculate_accumulated_CRN_concentration_from_points"] = {  "bool","false", "Calculates accumulated CRN concentrations from a list of points.","For forward modelling of 10Be concentrations"};

  string_default_map["cosmo_accumulated_concentration_points_fname"] = "NULL";
  help_map["cosmo_accumulated_concentration_points_fname"] = {  "string","accumulated_conc.csv", "The filename of the accumulated cosmo points.","You will need to be a bit careful with the filenames. You probably should include the nuclide in the filename."};

  string_default_map["cosmo_concentration_raster_prefix"] = "NULL";
  help_map["cosmo_concentration_raster_prefix"] = {  "string","NULL", "The prefix of the cosmogenic concentration raster.","You will need to be a bit careful with the filenames"};

  string_default_map["cosmo_accumulated_concentration_raster_prefix"] = "NULL";
  help_map["cosmo_accumulated_concentration_raster_prefix"] = {  "string","NULL", "The prefix of the cosmogenic accumulated concentration raster.","You will need to be a bit careful with the filenames"};

  bool_default_map["print_outlet_accumulated_concentration"] = false;
  help_map["print_outlet_accumulated_concentration"] = {  "bool","false", "This prints the accumulated concentration at the outlet.","This allows back calculation of the apparent erosion rate"};

  bool_default_map["sample_accumulated_conc_in_channels"] = false;
  help_map["sample_accumulated_conc_in_channels"] = {  "bool","false", "This prints the accumulated concentration in samples across the network","This allows back calculation of the apparent erosion rate"};

  bool_default_map["calculate_average_erosion_from_erosion_raster"] = false;
  help_map["calculate_average_erosion_from_erosion_raster"] = {  "bool","false", "Prints an raster with the average erosion rate from all upstream nodes","Used to look at the characteristics of upstream erosion in comparison to cosmogenic rates"};
  

  int_default_map["thinning_factor_for_sampling"] = 100;
  help_map["thinning_factor_for_sampling"] = {  "int","100","When you sample accumulated cosmo this is the nodes to skip in the channel network between samples.","Used to get a representative number of sampling points but not at every point along the channel."};

  int_default_map["threshold_pixels_for_sampling"] = 5000;
  help_map["threshold_pixels_for_sampling"] = {  "int","5000","The number of contributing pixels needed to start a channel used for sampling accumulated CRN. This will be larger than normal to avoid lots of sampling above the knickpoint.","This is in pixels not drainage area."};

  bool_default_map["erate_use_accumulated_points"] = false;
  help_map["erate_use_accumulated_points"] = {  "bool","false", "If you use the erate calculator, this grabs the accumulated points csv.","Allows checking of apparent erosion rate against accumulated rate in one step"};




  // some parameters for getting points in the landscape
  // You need a river network for this. 
  // You probably should keep the contributing pixels at the default unless you have either
  // very small or very large basins. Increase the number 
  bool_default_map["read_points_csv"] = false;
  help_map["read_points_csv"] = {  "bool","false", "For calculating concentrations of 10Be for scenarios this reads locations from a csv.","Assign the file with points_filename"};

  string_default_map["points_filename"] ="CRN_points.csv";
  help_map["points_filename"] = {  "string","CRN_points.csv", "Name of the file where you pick points to calculate accumulated 10Be","You need to set read_points_csv to true for this to work"};

  int_default_map["threshold_contributing_pixels"] = 1000;
  help_map["threshold_contributing_pixels"] = {  "int","1000","The number of contributing pixels needed to start a channel using the threshold method.","This is in pixels not drainage area. More options are in the lsdtt-channel-extraction tool."};


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


  // These are for the step change raster and associated tools
  bool_default_map["make_step_change_rasters"] = false;
  help_map["make_step_change_rasters"] = {  "bool","false","Switches on the step change raster routine.","To get at elevations and rates following a step change in uplift"};

  bool_default_map["carve_before_fill"] = false; // Implements a carving algorithm
  help_map["carve_before_fill"] = {  "bool","false","This implements a breaching algorithm before filling.","Good for landscapes with DEM obstructions (like roads) across the channels."};

  // step change getting the outlets
  bool_default_map["get_basins_from_outlets"] = false;
  help_map["get_basins_from_outlets"] = {  "bool","false","Switches on the outlet based basin finding.","See BaselevelJunctions_file for format of outlets csv."};
  
  int_default_map["search_radius_nodes"] = 8;
  help_map["search_radius_nodes"] = {  "int","8","A parameter for snapping to the nearest channel. It will search for the largest channel (by stream order) within the pixel window.","You will want smaller pixel numbers if you have a dense channel network."};
 
  int_default_map["threshold_stream_order_for_snapping"] = 2;
  help_map["threshold_stream_order_for_snapping"] = {  "int","2","If you are snapping to a channel the routine it will ignore channel with lower stream order than this number.","Set this to a higher number to avoid snapping to small channels."};
  
  string_default_map["basin_outlet_csv"] = "NULL";
  help_map["basin_outlet_csv"] = {  "string","NULL","A csv file with the lat long of basin outlets.","csv should have latitude and longitude columns and rows with basin outlets."};

  // step change upflift settings
  float_default_map["step_change_uplift_old"] = 0.0001;
  help_map["step_change_uplift_old"] = {  "float","0.0001", "Uplift rate in mm/yr of the old rate in the step change scenario","Used to calculate relief upstream of knickpoint"};

  float_default_map["step_change_uplift_new"] = 0.001;
  help_map["step_change_uplift_new"] = {  "float","0.001", "Uplift rate in mm/yr of the new rate in the step change scenario","Used to calculate relief downstream of knickpoint"};

  float_default_map["step_change_old_relief"] = 1000;
  help_map["step_change_old_relief"] = {  "float","1000", "The desired old relief of the landscape in metres","Used to back calculate K"};

  float_default_map["step_change_uplift_age"] = 100000;
  help_map["step_change_uplift_age"] = {  "float","100000", "The number of years ago the step change at the outlet happened","We assume all changes propigate upstream"};

  float_default_map["step_change_Sc_old"] = 0.25;
  help_map["step_change_Sc_old"] = {  "float","0.25", "Critical slope value in step change simulations. This is in the old erosion rate area","Dimensionless gradient on hillslopes"};

  float_default_map["step_change_Sc_new"] = 0.4;
  help_map["step_change_Sc_new"] = {  "float","0.4", "Critical slope value in step change simulations. This is in the new erosion rate area","Dimensionless gradient on hillslopes"};

  int_default_map["step_change_contributing_pixels"] = 500;
  help_map["step_change_contributing_pixels"] = {  "int","500", "For step change simulations this is the minimum conributing pixels for a channel","Anything with fewer than these contributing pixels will be calculated as a hillslope"};

  // step change parameters for chi and knickpoint propigation
  float_default_map["min_slope_for_fill"] = 0.0001;
  help_map["min_slope_for_fill"] = {  "float","0.0001","Minimum slope between pixels for the filling algorithm.","Best not to change the default."};

  float_default_map["A_0"] = 1.0;
  help_map["A_0"] = {  "float","1.0","The A_0 parameter for chi computation. See https://doi.org/10.1002/esp.3302","Usually set to 1 so that the slope in chi-elevation space is the same as k_sn"};
   
  float_default_map["m_over_n"] = 0.5;
  help_map["m_over_n"] = {  "float","0.5","The concavity index for chi calculations. Usually denoted as the Greek symbol theta.","Default is 0.5 but possibly 0.45 is better as Kwang and Parker suggest 0.5 leads to unrealistic behaviour in landscape evolution models."};

  float_default_map["m"] = 0.5;
  help_map["m"] = {  "float","0.5","The area exponent in the stream power law.","Defaults lead to m/n = 0.5"};

  float_default_map["n"] = 1;
  help_map["n"] = {  "float","1","The slope exponent in the stream power law.","Model slows down a lot if n does not equal 1 but there is much evidence that n is greater than 1 for example Harel et al 2016 Geomorphology."};

  // step change printing and reading of raster options
  string_default_map["time_since_change_raster_prefix"] = "NULL";
  help_map["time_since_change_raster_prefix"] = {  "string","NULL","This prefix of the raster with the time since the stepchange","Used for stepchange raster prediction."};

  // step change cosmo concentration prediction
  bool_default_map["calculate_step_change_CRN_concentration_raster"] = false;
  help_map["calculate_step_change_CRN_concentration_raster"] = {  "bool","false","Switches on the routine for calculating the step change concentration raster.","You will need to have the step change rasters calculated first"};
  
  float_default_map["rock_density_kg_m3"] = 2650;
  help_map["rock_density_kg_m3"] = {  "float","2650","The rock density in kg per m^3","Used to convert between L/T and M/L^2/T in effecti e erosion rates."};

  bool_default_map["convert_units_when_input_erosion_in_m_yr"] = false;
  help_map["convert_units_when_input_erosion_in_m_yr"] = {  "bool","false", "If the input erosion raster is in m/yr use this to convert it to effective erosion.","You will need to supply the density if this is turned on"};


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

  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./lsdtt-cosmo-tool-README.csv" << endl;
    string help_prefix = "lsdtt-cosmo-tool-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }


  cout << "Read filename is: " <<  DATA_DIR+DEM_ID << endl;
  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;

  // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);



  //============================================================================
  //
  //.##.......####....####...#####...........#####...######..##...##.
  //.##......##..##..##..##..##..##..........##..##..##......###.###.
  //.##......##..##..######..##..##..........##..##..####....##.#.##.
  //.##......##..##..##..##..##..##..........##..##..##......##...##.
  //.######...####...##..##..#####...........#####...######..##...##.
  //
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
  //
  //..####...##..##...####...##...##..........#####...#####...######..#####...#####....####....####...######...####....####..
  //.##......###.##..##..##..##...##..........##..##..##..##..##......##..##..##..##..##..##..##..##..##......##......##.....
  //..####...##.###..##..##..##.#.##..........#####...#####...####....#####...#####...##..##..##......####.....####....####..
  //.....##..##..##..##..##..#######..........##......##..##..##......##......##..##..##..##..##..##..##..........##......##.
  //..####...##..##...####....##.##...........##......##..##..######..##......##..##...####....####...######...####....####..
  //
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
  //..####...#####...##..##..........######..........#####....####...######..######.
  //.##..##..##..##..###.##..........##..............##..##..##..##....##....##.....
  //.##......#####...##.###..........####............#####...######....##....####...
  //.##..##..##..##..##..##..........##..............##..##..##..##....##....##.....
  //..####...##..##..##..##..........######..........##..##..##..##....##....######.
  //
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

  //============================================================================  
  //
  //..####...#####...##..##...........####....####...##......##..##..##...##..##..##...####..
  //.##..##..##..##..###.##..........##..##..##..##..##......##..##..###.###..###.##..##.....
  //.##......#####...##.###..........##......##..##..##......##..##..##.#.##..##.###...####..
  //.##..##..##..##..##..##..........##..##..##..##..##......##..##..##...##..##..##......##.
  //..####...##..##..##..##...........####....####...######...####...##...##..##..##...####..
  //
  //============================================================================  
  if(this_bool_map["transient_column_calculator"])
  {
    cout << "Let me make a column for you to test some concentrations." << endl;
    double this_conc  = -99;

    // Later this will need to be updated to get the nuclide on other factors from input to the program
    string Nuclide = this_string_map["nuclide_for_prediction"];
    string Muon_scaling = this_string_map["muon_scheme_for_prediction"];

    // the total shielding. A product of snow, topographic and production scaling
    double total_shielding = this_float_map["total_shielding_test"];

    double test_top_effective_depth = double(this_float_map["transient_effective_depth_test_top"]);
    double test_bottom_effective_depth = double(this_float_map["transient_effective_depth_test_bottom"]);

    // initiate a particle. We'll just repeatedly call this particle
    // for the sample.
    int startType = 0;
    double Xloc = 0;
    double Yloc = 0;
    double  startdLoc = 0.0;
    double  start_effdloc = double(this_float_map["transient_effective_depth_test"]);
    double startzLoc = 0.0;

    // create a particle at zero depth
    LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                                startdLoc, start_effdloc, startzLoc);

    // now create the CRN parameters object
    LSDCRNParameters LSDCRNP;

    // set up the scaling
    if (Muon_scaling == "Schaller" )
    {
      LSDCRNP.set_Schaller_parameters();
    }
    else if (Muon_scaling == "Braucher" )
    {
      LSDCRNP.set_Braucher_parameters();
    }
    else if (Muon_scaling == "BraucherBorchers" )
    {
      LSDCRNP.set_BraucherBorchers_parameters();
    }
    else if (Muon_scaling == "Granger" )
    {
      LSDCRNP.set_Granger_parameters();
    }
    else if (Muon_scaling == "newCRONUS" )
    {
      LSDCRNP.set_newCRONUS_parameters();
    }
    else
    {
      cout << "You didn't set the muon scaling." << endl
          << "Options are Schaller, Braucher, newCRONUS, BraucherBorchers, and Granger." << endl
          << "You chose: " << Muon_scaling << endl
          << "Defaulting to Braucher et al (2009) scaling" << endl;
      LSDCRNP.set_Braucher_parameters();
    }

    // set the scaling vector
    vector<bool> nuclide_scaling_switches(4,false);
    if (Nuclide == "Be10")
    {
      nuclide_scaling_switches[0] = true;
    }
    else if (Nuclide == "Al26")
    {
      nuclide_scaling_switches[1] = true;
    }
    else if (Nuclide == "C14")
    {
      nuclide_scaling_switches[3] = true;
    }
    else
    {
      cout << "You didn't choose a valid nuclide for your test column. Defaulting"
          << " to 10Be." << endl;
      Nuclide = "Be10";
      nuclide_scaling_switches[0] = true;
    }

    // Scale the F values
    LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

    cout << endl << endl << endl << "===============================================================" << endl;
    cout << "I am using the " << Muon_scaling << " scaling with a total shielding of: " << total_shielding << endl;
    cout << "The nuclide parameters are: " << endl;
    LSDCRNP.print_parameters_to_screen(nuclide_scaling_switches);
    cout << "===============================================================" << endl << endl;


    if(this_string_map["column_mode"] == "step_change")
    { 
      cout << "I am now going to simulate a step change in erosion rate" << endl;
      // get the nuclide concentration from this node
      if (Nuclide == "Be10")
      {
        //cout << "LInE 2271, 10Be" << endl;
        eroded_particle.update_10Be_step_change(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP);
        this_conc=eroded_particle.getConc_10Be();
      }
      else if (Nuclide == "Al26")
      {
        //cout << "LINE 2278, 26Al" << endl;
        eroded_particle.update_26Al_step_change(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP);
        this_conc=eroded_particle.getConc_26Al();
      }
      else if (Nuclide == "C14")
      {
        eroded_particle.update_14C_step_change(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP);
        this_conc=eroded_particle.getConc_14C();
      }
      else
      {
        cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
        cout << "Choices are Be10, Al26 or C14.  Note these case sensitive and cannot" << endl;
        cout << "contain spaces or control characters. Defaulting to Be10." << endl;
        eroded_particle.update_10Be_step_change(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP);
        this_conc=eroded_particle.getConc_10Be();
      }  
      cout << "Okay, I got particle with a step change." << endl;
      cout << "The nuclide is "<< Nuclide << endl;
      cout << "The depth of the particle, in effective depth units (g/cm^2) is: " << this_float_map["transient_effective_depth_test"] << endl;
      cout << "The initial erosion rate is: " << this_float_map["effective_erosion_rate"] << " and the new erosion rate is: " << this_float_map["effective_erosion_rate_new"] << endl;
      cout << "The time since the step is: " << this_float_map["time_since_step"] << endl;
      cout << " The concentration is: " << this_conc << " atoms/g"<< endl;  
    }
    else if(this_string_map["column_mode"] == "step_change_depth_integrated")
    { 
      cout << "I am now going to simulate a step change in erosion rate with depth integration" << endl;

      // get the nuclide concentration from this node
      if (Nuclide == "Be10")
      {
        //cout << "LInE 2271, 10Be" << endl;
        eroded_particle.update_10Be_step_change_depth_integrated(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP, test_top_effective_depth, test_bottom_effective_depth);
        this_conc=eroded_particle.getConc_10Be();
      }
      else if (Nuclide == "Al26")
      {
        eroded_particle.update_26Al_step_change_depth_integrated(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP, test_top_effective_depth, test_bottom_effective_depth);
        this_conc=eroded_particle.getConc_26Al();
      }
      else if (Nuclide == "C14")
      {
        
        eroded_particle.update_14C_step_change_depth_integrated(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP, test_top_effective_depth, test_bottom_effective_depth);
        this_conc=eroded_particle.getConc_14C();
      }
      else
      {
        cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
        cout << "Choices are Be10, Al26 or C14.  Note these case sensitive and cannot" << endl;
        cout << "contain spaces or control characters. Defaulting to Be10." << endl;
        eroded_particle.update_10Be_step_change_depth_integrated(this_float_map["effective_erosion_rate"], 
                                               this_float_map["effective_erosion_rate_new"], 
                                               this_float_map["time_since_step"],
                                               LSDCRNP, test_top_effective_depth, test_bottom_effective_depth);
        this_conc=eroded_particle.getConc_10Be();
      }  
      cout << "Okay, I got particle with a step change. This is DEPTH INTEGRATED" << endl;
      cout << "The nuclide is "<< Nuclide << endl;
      cout << "The top and  bottom effective depths are: " << test_top_effective_depth << ", " <<  test_bottom_effective_depth << endl;
      cout << "The initial erosion rate is: " << this_float_map["effective_erosion_rate"] << " and the new erosion rate is: " << this_float_map["effective_erosion_rate_new"] << endl;
      cout << "The time since the step is: " << this_float_map["time_since_step"] << endl;
      cout << " The concentration is: " << this_conc << " atoms/g"<< endl;  
    }
    else if(this_string_map["column_mode"] == "full_transient")
    {
      cout << "I haven't implemented the full transient run yet, exiting" << endl;
      exit(0);
    }
    else if(this_string_map["column_mode"] == "steady_state")
    {
      double this_top_eff_depth = 0;
      double this_bottom_eff_depth = 0;
      // get the nuclide concentration from this node
      if (Nuclide == "Be10")
      {
        //cout << "LInE 2271, 10Be" << endl;
        eroded_particle.update_10Be_SSfull_depth_integrated(this_float_map["effective_erosion_rate"],LSDCRNP,
                                            this_top_eff_depth, this_bottom_eff_depth);
        this_conc=eroded_particle.getConc_10Be();
      }
      else if (Nuclide == "Al26")
      {
        //cout << "LINE 2278, 26Al" << endl;
        eroded_particle.update_26Al_SSfull_depth_integrated(this_float_map["effective_erosion_rate"],LSDCRNP,
                                            this_top_eff_depth, this_bottom_eff_depth);
        this_conc=eroded_particle.getConc_26Al();
      }
      else if (Nuclide == "C14")
      {
        //cout << "LINE 2278, 26Al" << endl;
        eroded_particle.update_14C_SSfull_depth_integrated(this_float_map["effective_erosion_rate"],LSDCRNP,
                                            this_top_eff_depth, this_bottom_eff_depth);
        this_conc=eroded_particle.getConc_14C();
      }
      else
      {
        cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
        cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
        cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
        eroded_particle.update_10Be_SSfull_depth_integrated(this_float_map["effective_erosion_rate"],LSDCRNP,
                                            this_top_eff_depth, this_bottom_eff_depth);
        this_conc=eroded_particle.getConc_10Be();
      }  
      cout << "Okay, I got a steady state particle, the nuclide is "<< Nuclide << " and the concentration is: " << this_conc << " atoms/g"<< endl;   
    }
    else
    {
      cout << "You said you wanted a transient column calculator but did not give me an appropriate mode" << endl;
      cout << "Currently the functioning modes are step_change, full_transient, and steady_state" << endl;
      cout << "You chose: " << this_string_map["column_mode"] << endl;
      cout << "I am exiting the program." << endl;
      exit(0);
    }

    if(this_conc != -99)
    {
      cout << "Let me check what the apparent erosion rate would be in COSMOCALC" << endl;
      double this_top_eff_depth = 0;
      double this_bottom_eff_depth = 0;
      double rho = 2650;
      vector<double> erate_vec;
      if (Nuclide == "Be10")
      {
        erate_vec = eroded_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, total_shielding, Muon_scaling,
                                                        this_top_eff_depth, this_bottom_eff_depth);       
      }
      else if (Nuclide == "Al26")
      {
        erate_vec = eroded_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, total_shielding, Muon_scaling,
                                                        this_top_eff_depth, this_bottom_eff_depth);  
      }
      else if (Nuclide == "C14")
      {
        erate_vec = eroded_particle.apparent_erosion_14C_COSMOCALC(rho, LSDCRNP, total_shielding, Muon_scaling,
                                                        this_top_eff_depth, this_bottom_eff_depth);  
      }
      cout << "Apparent erosion rates are:" << endl;
      cout << erate_vec[0] << " g/cm^3/yr" << endl;
      cout << erate_vec[1] << " m/yr" << endl;
    }

  }



  //-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  //..####...######..######..#####...........#####....####....####...######..######..#####..
  //.##........##....##......##..##..........##..##..##..##..##........##....##......##..##.
  //..####.....##....####....#####...........#####...######...####.....##....####....#####..
  //.....##....##....##......##..............##..##..##..##......##....##....##......##..##.
  //..####.....##....######..##..............##..##..##..##...####.....##....######..##..##.
  //
  //-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (this_bool_map["make_step_change_rasters"])
  {
    cout << "Hi there. I am going to make some rasters to model a step change in erosion rate"<< endl;
    cout << "I'll be using some theory from Royden and Perron 2013" << endl;
    cout << "Which was nicely abbreviated by Mitchell and Yanites and Mitchell 2019" << endl;

    // First, get either the largest basin or a basin determined by a source point
    cout << endl << endl << endl << endl << "Before I start, I need to make sure the raster is filled. " << endl;
    LSDRaster carved_topography;
    LSDRaster filled_topography;
    LSDRasterInfo RI(topography_raster);

    if(this_bool_map["carve_before_fill"])
    {
      carved_topography = topography_raster.Breaching_Lindsay2016();
      filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
    }
    else
    {
      filled_topography = topography_raster.fill(this_float_map["min_slope_for_fill"]);
    }

    cout << "Flow routing in the initial DEM." << endl;
    LSDFlowInfo FI(filled_topography);

    cout << "Getting the junction network" << endl;
    LSDIndexRaster FlowAcc = FI.write_NContributingNodes_to_LSDIndexRaster();
    vector<int> sources;
    sources = FI.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);
    LSDRaster FlowDistance = FI.distance_from_outlet();
    LSDJunctionNetwork JunctionNetwork(sources, FI);    

    vector<int> BaseLevelJunctions;
    if(this_bool_map["get_basins_from_outlets"])
    {
      cout << "I am getting a specific basin to simulate" << endl;
      cout << "I am going to get basins lat-long coordinates" << endl;
      cout << "I need a csv with a latitude and a longitude location of the basin outlet." << endl;
      string full_BL_LL_name = DATA_DIR+this_string_map["basin_outlet_csv"];
      cout << "The file is: " << full_BL_LL_name << endl;
      int search_radius_nodes = this_int_map["search_radius_nodes"];
      int threshold_stream_order = 3;
      BaseLevelJunctions = JunctionNetwork.snap_point_locations_to_upstream_junctions_from_latlong_csv(full_BL_LL_name,
                                                               search_radius_nodes, threshold_stream_order,FI, RI);
    }
    else
    {
      cout << "You did not supply a basin outlet file so I will assume you want the largest complete catchment in this DEM." << endl;
      cout << "Let me get that for you now." << endl;
      BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Largest(BaseLevelJunctions, FI, FlowAcc);
    } 

    int single_junc = BaseLevelJunctions[0];
    cout << "The junction number is: " << single_junc << endl;
    int outlet_node = JunctionNetwork.get_Node_of_Junction(single_junc);
    vector<int> outlet_nodes;
    outlet_nodes.push_back(outlet_node);

    
    // get the elevation of the outlet
    int this_row,this_col;
    FI.retrieve_current_row_and_col(outlet_node,this_row,this_col);
    float outlet_elevation = filled_topography.get_data_element(this_row,this_col);
    cout << "The outlet elevation is: " << outlet_elevation << endl;

    // Now get the chi coordinate of every pixel in the basin
    LSDRaster chi_coordinate;       
    cout << "I am calculating the chi coordinate." << endl;
    float area_threshold_all_pixels_for_chi = 0.5;
    vector<int> starting_nodes;

    // get the chi coordinate
    chi_coordinate = FI.get_upslope_chi_from_multiple_starting_nodes(outlet_nodes,
                              this_float_map["m_over_n"], this_float_map["A_0"], 
                              area_threshold_all_pixels_for_chi);

    // print the chi coordinate to file
    string chi_fname = OUT_DIR+OUT_ID+"_SC_CHI";
    chi_coordinate.write_raster(chi_fname,raster_ext);

    // now get the upslope nodes
    vector<int> upslope_nodes = FI.get_upslope_nodes(outlet_node);
    int n_upslope = int(upslope_nodes.size());
    float max_chi = 0;
    for(int i = 0; i< int(upslope_nodes.size()); i++)
    {
      FI.retrieve_current_row_and_col(upslope_nodes[i],this_row,this_col);
      if ( chi_coordinate.get_data_element(this_row,this_col) > max_chi )
      {
        max_chi = chi_coordinate.get_data_element(this_row,this_col);
      }
    }
    FI.retrieve_current_row_and_col(upslope_nodes[n_upslope-1],this_row,this_col);
    float test_max_chi = chi_coordinate.get_data_element(this_row,this_col);
    cout << "max_chi: " << max_chi << " and test max chi = " << test_max_chi << endl;

    // IMPORTANT ALL BELOW ASSUMES A0 == 1
    // YOU WILL NEED TO USE DIFFERENT EQUATIONS IF THIS IS NOT THE CASE
    // MAKE MORE FLEXIBLE LATER

    // Now back calculate K using the max chi coordinate
    float calculated_K = this_float_map["step_change_uplift_old"]*pow(this_float_map["step_change_old_relief"]/max_chi,-this_float_map["n"]);
    cout << "The calculated K value is: " << calculated_K << endl;

    // calculate the various slopes and locations of the knickpoints analytically
    // first, the slope 
    float ks_upstream,ks_downstream;
    if (this_float_map["n"] != 1)
    {
      ks_upstream = pow(this_float_map["step_change_uplift_old"]/calculated_K,(1/this_float_map["n"]));
      ks_downstream = pow(this_float_map["step_change_uplift_new"]/calculated_K,(1/this_float_map["n"]));
    } 
    else 
    {
      ks_upstream = this_float_map["step_change_uplift_old"]/calculated_K;
      ks_downstream = this_float_map["step_change_uplift_new"]/calculated_K;
    }
    cout << "The kn values are: " << ks_upstream << " upstream and " << ks_downstream << " downstream " << endl; 
  
    // now calculate the location and elevation of the knickpoint
    // this comes from the mitchell and yanites paper
    float chi_knickpoint;
    float z_knickpoint;
    
    if (this_float_map["n"] == 1)
    {
      chi_knickpoint = calculated_K*this_float_map["step_change_uplift_age"];

      z_knickpoint = outlet_elevation+this_float_map["step_change_uplift_new"]*this_float_map["step_change_uplift_age"];

    }
    else if ( this_float_map["n"] > 1)
    {
      float numerator = (this_float_map["step_change_uplift_old"]-this_float_map["step_change_uplift_old"])*this_float_map["step_change_uplift_age"];
      float denominator1 = pow(this_float_map["step_change_uplift_old"]/calculated_K,1/this_float_map["n"]);
      float denominator2 = pow(this_float_map["step_change_uplift_new"]/calculated_K,1/this_float_map["n"]);
      chi_knickpoint = numerator/(denominator2-denominator1);


      float n1 = pow(this_float_map["step_change_uplift_new"]/calculated_K,1/this_float_map["n"]);
      float d1 = denominator2-denominator1;
      z_knickpoint = outlet_elevation+(this_float_map["step_change_uplift_new"]-this_float_map["step_change_uplift_old"])*this_float_map["step_change_uplift_age"]*(n1/d1);

    }
    else
    {
      float numerator = this_float_map["n"]*this_float_map["step_change_uplift_old"]*this_float_map["step_change_uplift_age"];
      float denominator = pow(this_float_map["step_change_uplift_old"]/calculated_K,1/this_float_map["n"]);
      chi_knickpoint = numerator/denominator;

      z_knickpoint = outlet_elevation + this_float_map["step_change_uplift_age"]*(this_float_map["step_change_uplift_new"]+(1-this_float_map["n"])*this_float_map["step_change_uplift_old"]);
    }

    cout << "The chi of the knickpoint is: " << chi_knickpoint << " and its elevation is: " << z_knickpoint << endl;

    // Now we calculate the elevations an erosion rates in the basin
    LSDCosmoRaster SC(chi_coordinate);
    LSDRaster erate;
    LSDRaster set_topography;
    LSDRaster time_since_step;



    SC.calculate_step_change_rasters(FI, upslope_nodes, chi_coordinate, erate, set_topography, 
                                        time_since_step, this_float_map["step_change_uplift_age"], 
                                        this_float_map["step_change_uplift_old"], this_float_map["step_change_uplift_new"],
                                        chi_knickpoint, z_knickpoint, this_float_map["n"], ks_upstream, ks_downstream, 
                                        outlet_elevation, calculated_K);  


    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;


    // Now do the hillslopes
    LSDRaster raster_tags;
    LSDRaster hillslope_topography;
    SC.calculate_step_change_hillslopes(set_topography, raster_tags, hillslope_topography, 
                                        erate, time_since_step, boundary_conditions, 
                                        this_int_map["step_change_contributing_pixels"], this_float_map["step_change_uplift_old"], 
                                        this_float_map["step_change_Sc_old"], this_float_map["step_change_Sc_new"]);  

    string ER_fname;
    if (this_string_map["background_erosion_rate_raster_prefix"] == "NULL")
    {
      ER_fname = OUT_ID+"_SC_ERATE";
    }
    else
    {
      ER_fname = this_string_map["background_erosion_rate_raster_prefix"];
    }
    string erate_fname = OUT_DIR+ER_fname;
    erate.write_raster(erate_fname,raster_ext);

    string z_fname = OUT_DIR+OUT_ID+"_SC_DEM";
    hillslope_topography.write_raster(z_fname,raster_ext);

    string T_fname;
    if (this_string_map["time_since_change_raster_prefix"] == "NULL")
    {
      T_fname = OUT_ID+"_SC_STEPTIME";
    }
    else
    {
      T_fname = this_string_map["time_since_change_raster_prefix"];
    }
    string time_fname = OUT_DIR+T_fname;
    time_since_step.write_raster(time_fname,raster_ext);

    string tag_fname = OUT_DIR+OUT_ID+"_SC_TAG";
    raster_tags.write_raster(tag_fname,raster_ext);

    LSDRaster hs_hs_raster = hillslope_topography.hillshade(hs_altitude,hs_azimuth,hs_z_factor);
    string hs_fname = OUT_DIR+OUT_ID+"_SC_DEM_hs";
    hs_hs_raster.write_raster(hs_fname,raster_ext);
  }



  //-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // CRN PREDICTION
  //
  //..####...#####...##..##..........#####...#####...######..#####...######...####...######
  //.##..##..##..##..###.##..........##..##..##..##..##......##..##....##....##..##....##..
  //.##......#####...##.###..........#####...#####...####....##..##....##....##........##..
  //.##..##..##..##..##..##..........##......##..##..##......##..##....##....##..##....##..
  //..####...##..##..##..##..........##......##..##..######..#####...######...####.....##..
  //
  //
  //-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // This makes all the appropriate rasters for constant values of the 
  // shielding and erosion. It prints the concentration raster,
  // and the erosion raster
  // The concentration is the concentration at individual pixels,
  // not accumulated concentration. 
  // For areas with landslides (e.g., nonzero pixel values in the self shielding raster)
  // the concentration is the average concentration of the evacuated material. 
  if(this_bool_map["calculate_CRN_concentration_raster"] || this_bool_map["calculate_step_change_CRN_concentration_raster"])
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
    
    // Now make or load the shield and other rasters
    LSDRaster SelfShield;
    if (this_string_map["self_shielding_raster_prefix"] != "NULL")
    {
      LSDRaster temp(DATA_DIR+this_string_map["self_shielding_raster_prefix"],raster_ext);
      SelfShield= temp;
    }
    else
    {
      MakeItYeah.set_to_constant_value(this_float_map["self_shield_eff_thickness"]);
      SelfShield =  MakeItYeah.return_as_raster(); 
    }
    
    LSDRaster SnowShield;
    if (this_string_map["snow_shielding_raster_prefix"] != "NULL")
    {
      LSDRaster temp(DATA_DIR+this_string_map["snow_shielding_raster_prefix"],raster_ext);
      SnowShield= temp;
    }
    else
    {
      MakeItYeah.set_to_constant_value(this_float_map["snow_shield_eff_thickness"]);
      SnowShield =  MakeItYeah.return_as_raster(); 
    }
    
    LSDRaster TopoShield;
    if (this_string_map["topographic_shielding_raster_prefix"] != "NULL")
    {
      LSDRaster temp(DATA_DIR+this_string_map["topographic_shielding_raster_prefix"],raster_ext);
      TopoShield= temp;
    }
    else
    {
      MakeItYeah.set_to_constant_value(1.0);
      TopoShield =  MakeItYeah.return_as_raster();  
    }

    LSDRaster ErosionRaster;
    if (this_string_map["background_erosion_rate_raster_prefix"] != "NULL")
    {
      cout << "I am loading an erosion raster. This is typically in g/cm^2/yr but for the step change it is in m/yr" << endl;
      LSDRaster temp(DATA_DIR+this_string_map["background_erosion_rate_raster_prefix"],raster_ext);
      ErosionRaster= temp;
    }
    else
    {
      MakeItYeah.set_to_constant_value(this_float_map["effective_erosion_rate"]);
      ErosionRaster =  MakeItYeah.return_as_raster(); 
    }

    LSDRaster time_since_change;
    if (this_string_map["time_since_change_raster_prefix"] != "NULL")
    {
      LSDRaster temp(DATA_DIR+this_string_map["time_since_change_raster_prefix"],raster_ext);
      time_since_change= temp;
    }
    else
    {
      MakeItYeah.set_to_constant_value(0.0);
      time_since_change =  MakeItYeah.return_as_raster();  
    }

       
    // Later this will need to be updated to get the nuclide on other factors from input to the program
    string Nuclide = this_string_map["nuclide_for_prediction"];
    string Muon_scaling = this_string_map["muon_scheme_for_prediction"];

    cout << "The nuclide you are using is: " << this_string_map["nuclide_for_prediction"] << endl;
    
    bool is_production_uncertainty_plus_on = false;
    bool is_production_uncertainty_minus_on = false;

    LSDRaster CRNConc;
    if(this_bool_map["calculate_step_change_CRN_concentration_raster"])    
    {
      CRNConc = ThisCosmoRaster.calculate_CRN_concentration_raster_step_change(Nuclide, Muon_scaling, 
                                           ErosionRaster, time_since_change,
                                           ProductionRaster, TopoShield,   
                                           this_float_map["step_change_uplift_old"],
                                           this_float_map["step_change_uplift_new"],
                                           this_float_map["rock_density_kg_m3"],
                                           is_production_uncertainty_plus_on,
                                           is_production_uncertainty_minus_on);
    }
    else
    {
      CRNConc = ThisCosmoRaster.calculate_CRN_concentration_raster(Nuclide, Muon_scaling, ErosionRaster,
                                          ProductionRaster, TopoShield, 
                                          SelfShield, SnowShield, 
                                          is_production_uncertainty_plus_on,
                                          is_production_uncertainty_minus_on);
    }

    // This writes the concentration raster 
    string conc_fname;
    if (this_string_map["cosmo_concentration_raster_prefix"] == "NULL")
    {
      conc_fname = OUT_ID+"_CONC";
    }
    else
    {
      conc_fname = this_string_map["cosmo_concentration_raster_prefix"];
    }
    string this_raster_name = OUT_DIR+conc_fname;                                  
    CRNConc.write_raster(this_raster_name,raster_ext);
    
  }


  //===================================================================================================
  //..####....####....####...##..##..##...##...........####....####...##..##...####..
  //.##..##..##..##..##..##..##..##..###.###..........##..##..##..##..###.##..##..##.
  //.######..##......##......##..##..##.#.##..........##......##..##..##.###..##.....
  //.##..##..##..##..##..##..##..##..##...##..........##..##..##..##..##..##..##..##.
  //.##..##...####....####....####...##...##...........####....####...##..##...####..
  //
  // This accumulates CRN. It wraps the concentration calculation
  //===================================================================================================
  if( this_bool_map["calculate_accumulated_CRN_concentration"] ||
      this_bool_map["calculate_accumulated_CRN_concentration_from_points"] ||
      this_bool_map["calculate_erosion_rates_new"] ||
      this_bool_map["calculate_average_erosion_from_erosion_raster"]) 
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

    // This will be used to make some secondary rasters
    LSDRasterMaker MakeItYeah(topography_raster);

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
    
    // Now make or load the shield and other rasters
    LSDRaster ErosionRaster;
    if (this_string_map["background_erosion_rate_raster_prefix"] != "NULL")
    {
      cout << "I am loading an erosion raster. This is typically in g/cm^2/yr but for the step change it is in m/yr" << endl;
      LSDRaster temp(DATA_DIR+this_string_map["background_erosion_rate_raster_prefix"],raster_ext);

      if(this_bool_map["convert_units_when_input_erosion_in_m_yr"])
      {
        float multiplier = this_float_map["rock_density_kg_m3"]*0.1;
        temp.raster_multiplier(multiplier);
      }
      ErosionRaster= temp;
    }
    else
    {
      cout << "I am accumulating based on a steady raster" << endl;
      cout << "The erosion rate is based on the parameter effective_erosion_rate" << endl;
      MakeItYeah.set_to_constant_value(this_float_map["effective_erosion_rate"]);
      ErosionRaster =  MakeItYeah.return_as_raster(); 
    }

    if(this_bool_map["calculate_average_erosion_from_erosion_raster"])
    {
      cout << "I am getting the averaged erosion rate. It will be in the same units as the erosion raster" << endl;
      LSDRaster AvgErosion =  FlowInfo.upslope_average(ErosionRaster);
      string this_fname = OUT_DIR+OUT_ID+"_AVGEROS";
      AvgErosion.write_raster(this_fname,raster_ext);

    }
    
    if(this_bool_map["calculate_accumulated_CRN_concentration"])
    {

      LSDRaster CRN_conc;
      if (this_string_map["cosmo_concentration_raster_prefix"] != "NULL")
      {
        LSDRaster temp(DATA_DIR+this_string_map["cosmo_concentration_raster_prefix"],raster_ext);
        CRN_conc= temp;
      }
      else
      {
        cout << "You have asked me to accumulate cosmo but you didn't give me a concentration raster" << endl;
        exit(0);
      }

      cout << endl << endl << "=================================" << endl;
      cout << "I will now calculate the accumulated concentration of your nuclide." << endl;
      LSDRaster AccConc = ThisCosmoRaster.calculate_accumulated_CRN_concentration(CRN_conc, ErosionRaster,FlowInfo);   
      cout << "Done with accumulating the concentration." << endl;

      // This writes the accumulation raster
      string conc_fname;
      if (this_string_map["cosmo_accumulated_concentration_raster_prefix"] == "NULL")
      {
        conc_fname = OUT_ID+"_ACCCONC";
      }
      else
      {
        conc_fname = this_string_map["cosmo_accumulated_concentration_raster_prefix"];
      }
      string this_raster_name = OUT_DIR+conc_fname;                                     
      AccConc.write_raster(this_raster_name,raster_ext);


      if(this_bool_map["sample_accumulated_conc_in_channels"])
      {
        cout << "Let me sample the accumulated concentration in the channels." << endl;
        string value_column_name = this_string_map["concentration_column_name"];
        int thinning_factor = this_int_map["thinning_factor_for_sampling"];
        int thresh_pixels_for_accumulated_conc = this_int_map["threshold_pixels_for_sampling"];

        // You need to get the node indices of the points first
        // get some relevant rasters
        LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();
        LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();     

        vector<int> sources;
        sources = FlowInfo.get_sources_index_threshold(ContributingPixels, thresh_pixels_for_accumulated_conc);

        // now get the junction network
        LSDJunctionNetwork ChanNetwork(sources, FlowInfo);

        // Print channels and junctions if you want them.

        cout << "I am going to print the thinned channel network now." << endl;
        string channel_csv_name = OUT_DIR+OUT_ID+"_SAMPLED_CN";
        ChanNetwork.PrintChannelNetworkToCSV_WithValuesThinned(FlowInfo, channel_csv_name, 
                                                                AccConc, thinning_factor,
                                                                value_column_name);
        cout << "I've printed the channel network. " << endl;        
      }

      if(this_bool_map["print_outlet_accumulated_concentration"])
      {
        int outlet = 0;     // this just gest the first value in the stack
        vector<int> cosmo_nodes;
        cosmo_nodes.push_back(outlet);

        string conc_pts_fname;
        if (this_string_map["cosmo_accumulated_concentration_points_fname"] == "NULL")
        {
          conc_pts_fname = OUT_ID+"_ACCCONC_pts.csv";
        }
        else
        {
          conc_pts_fname = this_string_map["cosmo_accumulated_concentration_points_fname"];
        }
        string acc_points_out_fname = OUT_DIR+conc_pts_fname;  
        string add_column_name = this_string_map["concentration_column_name"];
        FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(cosmo_nodes, acc_points_out_fname, AccConc, add_column_name);
      }
    }   // end logic for accumulation
  

    //=======================================================================================
    //
    //.##..##..######..##...##..........######..#####....####...######..######.
    //.###.##..##......##...##..........##......##..##..##..##....##....##.....
    //.##.###..####....##.#.##..........####....#####...######....##....####...
    //.##..##..##......#######..........##......##..##..##..##....##....##.....
    //.##..##..######...##.##...........######..##..##..##..##....##....######.
    //
    //=======================================================================================    
    if(this_bool_map["calculate_erosion_rates_new"])
    {
      cout << endl << endl << "===================================" << endl;
      cout << "I am now entering the erosion rate calculator." << endl;
      
      string Nuclide = this_string_map["nuclide_for_prediction"];
      string Muon_scaling = this_string_map["muon_scheme_for_prediction"];

      LSDSpatialCSVReader CRN_points_data;
      if (this_bool_map["erate_use_accumulated_points"])
      {
        cout << "I'm going to use accumulated points for the point file." << endl;
        cout << "This is typicall calculated alongside accumulation rasters to check" << endl;
        cout << "If the erosion calculator is working or for step change simulations" << endl;
        string acc_points_out_fname = OUT_DIR+OUT_ID+"_SAMPLED_CN.csv";
        cout << "The points file is " << acc_points_out_fname << endl;
        LSDSpatialCSVReader CRNPD( RI, acc_points_out_fname );
        CRN_points_data = CRNPD;
      }
      else
      {
        cout << "I am reading points from the file: "+ this_string_map["points_filename"] << endl;
        LSDSpatialCSVReader CRNPD( RI, (DATA_DIR+this_string_map["points_filename"]) );
        CRN_points_data = CRNPD;
      } 

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
      
      int search_radius_nodes = 4;
      vector<int> valid_cosmo_points;
      vector<int> snapped_node_indices;
      vector<int> snapped_junction_indices;

      // Now get the snapped points 
      // The snapped points are in a vector called snapped_node_indices  
      cout << "I'm snapping some points" << endl;
      cout << "I am going for the closest channel of first order." << endl; 
      int snap_threshold_stream_order = 1;
      ChanNetwork.snap_point_locations_to_nearest_channel_node_index(fUTM_easting, fUTM_northing, 
            search_radius_nodes, snap_threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices);
      cout << "The number of snapped points is: " <<   valid_cosmo_points.size() << endl;   


      // Now make or load the shield and other rasters
      LSDRaster SelfShield;
      if (this_string_map["self_shielding_raster_prefix"] != "NULL")
      {
        LSDRaster temp(DATA_DIR+this_string_map["self_shielding_raster_prefix"],raster_ext);
        SelfShield= temp;
      }
      else
      {
        MakeItYeah.set_to_constant_value(this_float_map["self_shield_eff_thickness"]);
        SelfShield =  MakeItYeah.return_as_raster(); 
      }
      
      LSDRaster SnowShield;
      if (this_string_map["snow_shielding_raster_prefix"] != "NULL")
      {
        LSDRaster temp(DATA_DIR+this_string_map["snow_shielding_raster_prefix"],raster_ext);
        SnowShield= temp;
      }
      else
      {
        MakeItYeah.set_to_constant_value(this_float_map["snow_shield_eff_thickness"]);
        SnowShield =  MakeItYeah.return_as_raster(); 
      }
      
      LSDRaster TopoShield;
      if (this_string_map["topographic_shielding_raster_prefix"] != "NULL")
      {
        LSDRaster temp(DATA_DIR+this_string_map["topographic_shielding_raster_prefix"],raster_ext);
        TopoShield= temp;
      }
      else
      {
        MakeItYeah.set_to_constant_value(1.0);
        TopoShield =  MakeItYeah.return_as_raster();  
      }

      LSDRaster NoDataRaster;
      MakeItYeah.set_to_constant_value(ProductionRaster.get_NoDataValue());
      NoDataRaster = MakeItYeah.return_as_raster();

      // The erosion rate calculation
      // This version will work from a specific node rather than a channel junction
      // It uses a hybrid newton-raphson and bisection method (the latter if the former fails)
      // to converge on the correct erosion rate
      // There is no error analysis in this version, you just get the closest number. 
      // For a representative error analysis we would need to know local variation in
      // quartz content, erosion rate, storage, etc. So instead I suggest a blanket 20% uncertainty.
      // You only need a concentration for this routine
      // you tell the routine which nuclide to use in the driver file
      // In this case we do not use nested basins so we do not use an erosion raster
      // there is also an option for including a quartz concentration but that has yet to be implemented
      // SMM, 04 June 2022
      if( this_bool_map["calculate_erosion_rates_new"])
      {
        cout << "Let me calculate the erosion rate for you. "  << endl;

        // Set up the output file
        string erate_calc_fname = OUT_DIR+OUT_ID+"_CN_ERATE.csv";
        string erate_calc_column = "apparent_eff_e";
        double this_lat, this_long;
        float predicted_conc;
        LSDCoordinateConverterLLandUTM Converter;


        // open the outfile
        ofstream csv_out;
        csv_out.open(erate_calc_fname.c_str());
        csv_out.precision(8);
        csv_out << "latitude,longitude,concentration,predicted_conc,apparent_eff_e,apparent_e_mm_yr" << endl;

        // This looks for concentrations data
        cout << "The concentration column name is: " << this_string_map["concentration_column_name"] << endl;
        vector<float> concentrations = CRN_points_data.data_column_to_float(this_string_map["concentration_column_name"]);

        // Now loop through concentrations, getting the erosion rates
        int n_valid_samples = int(valid_cosmo_points.size());
        for(int sample = 0; sample< n_valid_samples; sample++)
        {
          cout << endl << endl << "Processing sample " << sample+1 << " of " << n_valid_samples << ". ";
          cout << "The concentration is: " << concentrations[ valid_cosmo_points[sample]  ] << endl;

          float This_erate = ThisCosmoRaster.calculate_eff_erate_from_conc(concentrations[ valid_cosmo_points[sample] ],
                                            Nuclide, Muon_scaling, NoDataRaster,
                                            ProductionRaster, TopoShield, 
                                            SelfShield, SnowShield, 
                                            NoDataRaster,FlowInfo,
                                            snapped_node_indices[sample],
                                            predicted_conc);

          FlowInfo.get_lat_and_long_from_current_node(snapped_node_indices[sample], this_lat, this_long, Converter);
          cout << "The erosion rate is: " << This_erate << endl;
          csv_out << this_lat << "," << this_long << "," << concentrations[ valid_cosmo_points[sample] ]
                  << "," << predicted_conc << "," << This_erate << "," 
                  << 10000*This_erate/this_float_map["rock_density_kg_m3"] << endl;
        }
        csv_out.close();
      }
    }  // end logic for new erate

  } // end logic for accumulation and erate, which both involve flowinfo and production



  cout << "I'm all finished! Have a nice day." << endl;
}
