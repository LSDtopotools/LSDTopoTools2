//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDCosmoData.hpp
// Land Surface Dynamics CosmoData
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for keeping track of cosmogenic data
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <fstream>
#include <math.h>
#include <iostream>
#include <map>
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
using namespace std;

#ifndef LSDCosmoData_HPP
#define LSDCosmoData_HPP

class LSDCosmoData
{
  public:
  
    /// @brief the default constructor. This doesn't do anything. 
    LSDCosmoData()                                { create(); }
    
    /// @brief This constructor requires a filename
    /// @detail The file format is: 
    /// sample_name nuclide latitude longitude concentration uncertainty standard
    ///  the sample name cannot have a space 
    ///  the nuclide name must be Be10 or Al26
    ///  lat and long are in decimal degrees
    ///  concetration is in atoms per gram
    ///  uncertainty is the AMS error in atoms per gram
    ///  standard is the laboratory standardisation. The options can be found
    ///  in the code of this function; the numbers come from Balco et al's (2009)
    ///  CRONUS calculator. 
    /// @param filename the string of the file within which the cosmo data
    ///  is stored. This INCLUDES PATH AND EXTENSION.
    /// @param ftype the type of file, options are csv and txt. If not one of these
    ///  options the code assumes txt
    /// @author SMM
    /// @date 06/02/2015
    LSDCosmoData( string path, string file_prefix)  { create(path,file_prefix); }
    
    //// @brief This initiates some default parameters
    /// @detail These parameters will be overwritten by any parameters
    ///  supplied by the parameter file
    /// @author SMM
    /// @date 14/02/2016
    void initiate_default_parameters();
    
    
    /// @brief  THis function loads the crn data rom either a csv or text file
    /// @param filename the name of the file
    /// @param filetype this is either csv or txt
    /// @author SMM
    /// @date 26/02/2015
    void load_cosmogenic_data(string filename, string filetype);
    
    /// @brief This function load a csv file containing names of DEMs and 
    ///  (possibly) sheilding rasters or shielding parameters
    /// @detail The file should have one DEM name per row. Each row
    ///  can have between 1 and 4 elements. 
    /// The elements are:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    /// If elements are missing, they are considered null arguments. 
    /// @param filename the name of the csv file
    /// @author SMM
    /// @date 26/02/2015
    void load_DEM_and_shielding_filenames_csv(string filename);
    
    /// @brief This loads parameters for the comsogenic calculations
    ///  it uses the parse_files frinction in LSDStatsTools, so the format
    ///  is the paramterer name, followed by ":", followed by a space and then
    ///  the parameter value
    /// @param filename a string of the full filename
    /// @author SMM
    /// @date 02/02/2016
    void load_parameters(string filename);
   
    /// @brief This loads information about any soil samples
    /// The file contains several columns:
    /// sampleID,sample_top_depth,sample_bottom_depth,density
    /// The file prefix needs to be the same as the other files, and should
    /// have the extension _CRNSoilInfo.csv
    /// @detail The top depths and sample thicknesses need to be in cm
    ///   The densities need to be in kg/m^3
    /// @param filename a string of the full filename
    /// @author SMM
    /// @date 02/03/2015
    void load_soil_info(string filename);

    /// @brief This checks on the filenames to see if the files exist.
    /// @param crn_fname a string with the name (including path) of the CRN file
    /// @param Rasters_fname a string with the name (including path) of the Rasters file
    /// @param parameters_fname a string with the name (including path) of the parameters file
    /// @param soil_fname a string with the name (including path) of the soil file
    /// @author SMM
    /// @date 14/02/2016 
    void check_files(string crn_fname,string Rasters_fname,string parameters_fname,
                               string soil_fname);

    /// @brief this gets the names of the DEMs to be used in the analysis
    /// @detail only returns the DEM, not snow shielding, topo shielding, etc
    ///  rasters
    /// @return DEM_fnames a vector of fname strings
    /// @author SMM 
    /// @date 19/03/2015
    vector<string> get_DEM_fnames();
    
    /// @brief this gets the names of the snow shielding rasters
    //  to be used in the analysis
    /// @detail returns only snow shielding names. If name does not exist,  
    ///  returns NULL.
    /// @return Snow_fnames a vector of fname strings
    /// @author SMM 
    /// @date 07/07/2015
    vector<string> get_Snow_fnames();

    /// @brief this gets the names of the self shielding rasters
    //  to be used in the analysis
    /// @detail returns only snow shielding names. If name does not exist,  
    ///  returns NULL.
    /// @return Self_fnames a vector of fname strings
    /// @author SMM 
    /// @date 07/07/2015
    vector<string> get_Self_fnames();

    /// @brief This function checks to make sure parameter values are
    ///  valid for the cosmo data
    /// @author SMM
    /// @date 03/03/2015
    void check_parameter_values();
    
    /// @brief this function checks the existence and georeferencing of 
    ///  the rasters outlined in the file list
    /// @author SMM
    /// @date 03/03/2015
    void check_rasters();

    /// @brief This sets the path to atmospheric data
    /// @author SMM
    /// @param the new path to atmostperic data
    /// @date 22/04/2020
    void set_path_to_atmospheric_data(string atmos_path);
    
    /// @detail This function takes a DEM and then spawns dems that are clipped
    ///  to the basin boundaries, with padding to incorporate surrounding high
    ///  topography and CRN samples that might lie outside the basin due to 
    ///  GPS errors
    /// @param DEM_fname the name of the DEM with FULL PATH but excluding extension
    /// @param basin_padding_px the number of pixels that should be padded from
    ///  the basin
    /// @author SMM
    /// @date 18/03/2015
    vector<string> spawn_clipped_basins(string DEM_fname, int basin_padding_px);

    /// @brief This dirves the spawning of basins
    /// @param path This is a string containing the path to the data files (needs / at the end)
    /// @param prefix the prefix of the data files
    /// @param padding_pixels the number of pixels with which to pad the basins
    /// @author SMM
    /// @date 10/07/2015
    void BasinSpawnerMaster(string path, string prefix, int padding_pixels);

    /// @brief This calculates topographic shielding for basins listed in the 
    ///   _CRNRasters.csv file
    /// @detail Shielding rasters are printed to the same folder as the DEM
    /// @param path This is a string containing the path to the data files (needs / at the end)
    /// @param prefix the prefix of the data files
    /// @author SMM
    /// @date 15/07/2015
    void RunShielding(string path, string prefix);

    /// @brief This sets topographic shielding for basins listed in the 
    ///   _CRNRasters.csv file to 1
    /// @detail Shielding rasters (of 1) are printed to the same folder as the DEM
    /// @param path This is a string containing the path to the data files (needs / at the end)
    /// @param prefix the prefix of the data files
    /// @author SMM
    /// @date 22/04/2020
    void RunShielding_Unshielded(string path, string prefix);

    /// @brief This function calculates and then returns a production raster
    /// @param Elevation_data a raster holding the elevations
    /// @param path_to_atmospheric_data a string that holds the path of the atmospheric data
    /// @author SMM
    /// @date 28/01/2016
    LSDRaster calculate_production_raster(LSDRaster& Elevation_Data,
                                          string path_to_atmospheric_data);

    /// @brief this function calculates the UTM coordinates of all the sample
    ///  points for a given UTM zone. 
    /// @param UTM_zone the UTM zone
    /// @author SMM
    /// @date 06/02/2015
    void convert_to_UTM(int UTM_zone);
    
    /// @brief this function calculates the UTM coordinates of all the sample
    ///  points for a given UTM zone. It determines the UTM zone from a raster
    /// @param Raster the LSDRaster from which the UTM zone is determined 
    /// @author SMM
    /// @date 06/02/2015
    void convert_to_UTM(LSDRaster& Raster);    

    /// @brief This function calculates the CRN erosion rate from the data
    ///  stored within the data elements of the object
    /// @detail This analysis does not calculate snow shielding or landsliding
    /// @param DEM_prefix a string holding the prefix of the DEM (without the bil)
    ///  NOTE the DEM must be in bil format.
    /// @author SMM
    /// @date 09/02/2015 
    void basic_cosmogenic_analysis(string DEM_prefix);

    /// @brief This function computes erosion rates and uncertainties for 
    ///  a given DEM. It is wrapped by a function that goes through
    ///  the list of DEM, providing this function with the raster names
    ///  and the parameters for the model run
    /// @param Raster_names a vector of strings with 4 elements:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    ///  It there is no DEM then this is set to "NULL"
    ///  @param CRN_params this contains the single shielding depths for snow
    ///   and self shielding if the rasters are not supplied. 
    /// @author SMM
    /// @date 28/02/2015
    void full_shielding_cosmogenic_analysis(vector<string> Raster_names,
                            vector<double> CRN_params);

    /// @brief This function computes erosion rates and uncertainties for 
    ///  a given DEM. It is wrapped by a function that goes through
    ///  the list of DEM, providing this function with the raster names
    ///  and the parameters for the model run . This version can take a raster
    ///  with known erosion rates that can be used for nesting of basins
    /// @param Raster_names a vector of strings with 4 elements:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    ///  It there is no DEM then this is set to "NULL"
    ///  @param CRN_params this contains the single shielding depths for snow
    ///   and self shielding if the rasters are not supplied. 
    /// @param known_eff_erosion an LSDRaster with known erosion rates in g/cm^2/yr
    /// @author SMM
    /// @date 12/02/2016
    void full_shielding_cosmogenic_analysis_nested(vector<string> Raster_names,
                            vector<double> CRN_params, 
                            LSDRaster& known_eff_erosion);

    /// @brief This function computes erosion rates and uncertainties for 
    ///  a given DEM. It is wrapped by a function that goes through
    ///  the list of DEM, providing this function with the raster names
    ///  and the parameters for the model run. This version is for spawned rasters
    /// @detail Because this is spawned, each raster only corresponds to one 
    ///   sample, which the function finds from the sample name corresponding
    ///   to the last string in the raster name after the `_` character.
    /// @param Raster_names a vector of strings with 4 elements:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    ///  It there is no DEM then this is set to "NULL"
    ///  @param CRN_params this contains the single shielding depths for snow
    ///   and self shielding if the rasters are not supplied. 
    /// @author SMM
    /// @date 30/01/2016
    void full_shielding_cosmogenic_analysis_for_spawned(vector<string> Raster_names,
                            vector<double> CRN_params);

    /// @brief This function computes erosion rates and uncertainties for 
    ///  a given det of soil data, 
    /// @param Raster_names a vector of strings with 4 elements:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    ///  It there is no DEM then this is set to "NULL"
    ///  @param CRN_params this contains the single shielding depths for snow
    ///   and self shielding if the rasters are not supplied. 
    /// @author SMM
    /// @date 22/01/2016
    void Soil_sample_calculator(vector<string> Raster_names,
                            vector<double> CRN_params);

    /// @brief This function wraps the erosion rate calculator, and returns 
    ///  both the erosion rate as well as the uncertainties
    /// @param Nuclide_conc Concetration of the nuclide
    /// @param Nuclide a string denoting the name of the nuclide (at the moment
    ///  options are 10Be and 26Al)
    /// @param Nuclide_conc_err The instrument error in the nuclide concentration
    /// @param prod_uncert_fracton This is a fraction of the total uncertainty
    ///  for the production rates. It is a lumped parameter that can be used
    ///  for just production, or for snow, topo and porduction uncertainty
    /// @param Muon_scaling string that gives the muon scaling scheme
    ///  options are Schaller, Granger and Braucher
    /// @param snow_eff_depth the effective depth (in g/cm^2) of snow
    /// @param self_eff_depth the effective depth (in g/cm^2) of self shielding
    /// @param topo_shield the topographic shielding at this point
    /// @param production the production rate of the nuclude at this point
    /// @return  a vector of both the erosion rates and the uncertainties of the sample
    /// @author SMM
    /// @date 24/01/2016
    vector<double> full_CRN_erosion_analysis_point(double Nuclide_conc, string Nuclide, 
                            double Nuclide_conc_err, double prod_uncert_factor,
                            string Muon_scaling, double snow_eff_depth,
                            double self_eff_depth, double topo_shield,
                            double production);

    /// @brief this uses Newton Raphson iteration to retrieve the erosion rate
    ///  from a basin given a nuclide concentration
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production 
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme. 
    ///  options are Schaller, Braucher and Granger
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function. 
    /// @param average_production This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms   
    /// @param is_production_uncertainty_plus_on a boolean that is true if the 
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the 
    ///  production rate uncertainty (-) is switched on. If the + switch is 
    ///  true this parameter defauts to false.         
    /// @return The effective erosion rate in g/cm^-2/yr
    /// @author SMM
    /// @date 24/01/2016
    double predict_CRN_erosion_point(double Nuclide_conc, string Nuclide, 
                               double prod_uncert_factor,string Muon_scaling,
                               double& production_uncertainty,
                               double& average_production,
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on,
                               double snow_eff_depth,
                               double self_eff_depth, double topo_shield,
                               double production);


    /// @brief this predicts the mean concentration of a nuclide for a point location.
    ///  It does a full analyitical solution to account for
    ///  snow and self sheilding
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production 
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme. 
    ///  options are Schaller, Braucher and Granger
    /// @param data_from_outlet_only boolean that is true of you want 
    ///  concentration calculated from the outlet only.
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param is_production_uncertainty_plus_on a boolean that is true if the 
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the 
    ///  production rate uncertainty (-) is switched on. If the + switch is 
    ///  true this parameter defauts to false. 
    /// @param snow_eff_depth the effective depth (in g/cm^2) of snow
    /// @param self_eff_depth the effective depth (in g/cm^2) of self shielding
    /// @param topo_shield the topographic shielding at this point
    /// @param production the production rate of the nuclude at this point
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 24/01/2016
    double predict_mean_CRN_conc_point(double eff_erosion_rate, string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            double& production_uncertainty, 
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on,
                                            double snow_eff_depth,
                                            double self_eff_depth,
                                            double topo_shield,
                                            double production);



    /// @brief This function wraps the cosmogenic rate calculators.
    /// @detail Looks throught the vecvecs listing file locations and then
    ///   finds valid CRN data, and runs the erosion rate routine for these
    ///   data.
    /// @param method_flag: an int that sets the method
    ///  default: use snow and self shielding
    ///  method_flag == 0, basic analysis (snow and self set to 0, calculate toposhield)
    ///  method_flag == 1, snow and self sheilding, can load toposhield if provided
    /// @author SMM
    /// @date 28/02/2015
    void calculate_erosion_rates(int method_flag);


    /// @brief This function wraps the cosmogenic rate calculators. THis one is used with
    /// nested basins.
    /// @detail Looks throught the vecvecs listing file locations and then
    ///   finds valid CRN data, and runs the erosion rate routine for these
    ///   data.
    ///  It assumes that any known erosion rates will have the same name as the
    ///   base DEM with the extension _ERKnown
    /// @author SMM
    /// @date 23/02/2016    
    void calculate_nested_erosion_rates();

    /// @brief This function prints sevear rasters to file:
    ///  1) Pixel-by-pixel production scaling
    ///  2) Pixel-by-pixel combined scaling (production plus combined shielding)
    ///  3) Pixel-by-pixel combined shielding
    ///  4) The cosmo concetration coming from each pixel at the proscribed erosion rate
    /// @param Raster_names a vector containing the names of DEM, snow, self and toposhield rasters
    /// @param CRN_params a vector of parameter values
    /// @author SMM
    /// @date 23/03/2015
    void full_shielding_raster_printer(vector<string> Raster_names,
                                        vector<double> CRN_params);

    /// @brief This function uses the COSMOCALC functions of LSDParticles to
    ///  calculate the erosion rates from point measurements using pre-calculated
    ///  vectors of shielding and scaling data
    /// @param valid_samples a vector continaing the indices into the samples
    ///  you want to use
    /// @param snow_thickness a vector of snow thickness in g/cm^2 (or any surface shielding)
    ///  that is the same size as valid_samples
    /// @param self_thickness the self thickness in g/cm^2. For soil samples
    ///  this can be the thickness of the sample. Same size as valid_samples
    /// @param toposhield a vecor with the local topographic shielding
    /// @param production_scaling a vector with the production scaling
    /// @param muon_scaling a string with a valid scaling scheme name
    /// @author SMM
    /// @date 29/01/2016
    void point_measurements(vector<int> valid_samples,vector<double> snow_thickness, 
                                      vector<double> self_thickness,
                                      vector<double> toposhield,
                                      vector<double> production_scaling, 
                                      string muon_scaling);



    /// @brief This function calculates the concentration of a CRN
    ///  (either 10Be or 26Al) given an erosion rate supplied by and 
    ///  erosion rate raster (in g/cm^2/yr)
    /// @param Raster_names a vector of strings with 4 elements:
    ///  [0] = DEM_filename
    ///  [1] = Snow_shield_raster_name OR const_snow_shield in g/cm^2
    ///  [2] = Self_shield_raster_name OR const_self_shield in g/cm^2
    ///  [3] = Toposhield_raster_name 
    ///  It there is no DEM then this is set to "NULL"
    ///  @param CRN_params this contains the single shielding depths for snow
    ///   and self shielding if the rasters are not supplied. 
    /// @param known_eff_erosion an LSDRaster with known erosion rates in g/cm^2/yr
    /// @author SMM
    /// @date 14/06/2016
    void full_shielding_CRN_concentration_predictor(vector<string> Raster_names,
                            vector<double> CRN_params, 
                            LSDRaster& known_eff_erosion);



    /// @brief this function prints the data held in the the data members
    ///  to screen. Is used for bug checking. 
    /// @detail Note the function does not print the standardised values, only raw values.
    /// @author SMM
    /// @date 09/02/2015
    void print_data_to_screen(); 

    /// @brief This prints the cosmo data to a new file 
    /// @param path is the pathname. Needs a "/"  at the end
    /// @param prefix: the name of the file before _CRNData.csv
    /// @author SMM
    /// @date 14/07/2015
    void print_renamed_cosmo_data(string path, string prefix);

    /// @brief This prints the parameter data to a new file 
    /// @param path is the pathname. Needs a "/"  at the end
    /// @param prefix: the name of the file before _CRNData.csv
    /// @author SMM
    /// @date 14/07/2015
    void print_renamed_parameter_data(string path, string prefix);

    /// @brief this prints the file structure data to screen
    ///  it is a list of the DEMs, snow shielding rasters, self shelding rasters
    ///  and topo shielding rasters used in the analysis
    ///  the later three rasters can be NULL values, and the snow and self shielding
    ///  rasters can be replaced by single values
    /// @author SMM
    /// @date 26/02/2015
    void print_file_structures_to_screen();

    /// @brief this prints all the data, parameters and file structures to screen
    ///  It is used to keep a record of CRN erosion rate computations so
    ///  analyses can be reproduced
    /// @param outfilename the name of the outfile including full path
    /// @author SMM
    /// @date 02/03/2015
    void print_all_data_parameters_and_filestructures(string outfilename);

    /// @brief Prints the simple results to the screen
    /// @detail The 'simple' is because it only looks at external, 
    ///  muon, and production uncertainties. 
    /// @param rho the rock denisty in kg/m^3
    /// @author SMM
    /// @date 10/02/2015
    void print_simple_results_to_screen(double rho);

    /// @brief This prints the results to a csv file and to a file for
    ///  passing to CRONUS
    /// @detail The columns in the CSV file are:
    ///  sample_name
    ///  nuclide
    ///  latitude
    ///  longitude
    ///  concentration (atms/g)
    ///  concentration_uncert (atoms/g)
    ///  erosion rate g_percm2_peryr
    ///  erosion rate AMS_uncert g_percm2_peryr
    ///  muon_uncert g_percm2_peryr
    ///  production_uncert g_percm2_peryr
    ///  total_uncert g_percm2_peryr
    ///  AvgProdScaling dimensionless
    ///  AverageTopoShielding dimensionless
    ///  AverageSelfShielding dimensionless
    ///  AverageSnowShielding dimensionless
    ///  AverageCombinedScaling dimensionless (this is averaged production scaling times toposhielding)
    ///  outlet_latitude
    ///  OutletPressure hPa
    ///  OutletEffPressure hPa (pressure needed to get basin averaged production scaling)
    ///  centroid_latitude
    ///  CentroidPressure hPa
    ///  CentroidEffPressure (pressure needed to get basin averaged production scaling)
    ///  ErosionRate_in_cmperkyr (to check against cosmocalc, assumes 2650 kg/m^2)
    ///  ErosionRate_COSMOCALC_in_g_percm2_peryr (assumes 2650 kg/m^2): The erosion
    ///   rate you would get if you took production weighted scaling and used
    ///   cosmocalc. 
    ///  ErosionRate_COSMOCALC_cmperkyr (assumes 2650 kg/m^2): The erosion
    ///   rate you would get if you took production weighted scaling and used
    ///   cosmocalc. 
    /// @author SMM
    /// @date 12/03/2015
    void print_results();

    /// @brief This function prints several rasters to file:
    ///  1) Pixel-by-pixel production scaling
    ///  2) Pixel-by-pixel combined scaling (production plus combined shielding)
    ///  3) Pixel-by-pixel combined shielding
    ///  4) The cosmo concetration coming from each pixel at the proscribed erosion rate
    /// @author SMM
    /// @date 23/03/2015
    void print_rasters();
    
    /// @brief This prints scaling, shielding production and other rasters
    ///  for entire DEMs.  It only prints topo, snow and self shielding if those
    ///  rasters exists. That is, it does not calculate toposhielding automatically.
    /// @detail The printed rasters are:
    ///  extension _PRES the atmospheric pressure in hPa, determined from the
    ///  NCEP atmospheric pressure compilation following Balco et al 2008
    ///  extension _PROD the production scaling, based on stone/lal
    ///  extension _CSHIELD the combined sheilding, that is the combination of
    ///  snow, self and topographic shielding. Snow and self shielding are
    ///  approximated with spallation only production
    ///  extension _CSCALE the product of the combined shielding and production
    ///  for each pixel.
    /// @author SMM
    /// @date 15/02/2016
    void print_scaling_and_shielding_complete_rasters();

    /// @brief This prints the basins to a raster, along with the stream order raster
    ///  so that users can check if their basins are in the correct place
    /// @detail THis prints out a raster for each raster supplied in the CRNRaster
    ///  data file that contains the basins, a file with a key to the basins, 
    ///  and a stream order raster to help in refining the location of the data
    /// @author SMM
    /// @date 11/04/2016
    void print_basins_to_for_checking();

  protected:
    
    /// the number of samples
    int N_samples;
    
    /// the path to the cosmo data. Also used to print results
    string path;
    
    /// the prefix of the parameter files
    string param_name;
    
    /// A vector of the sample names
    vector<string> sample_name;
    
    /// A vector with the indices into the valid samples that have soil data
    vector<int> soil_sample_index;
    
    /// this is an index that points from a sample to an entry in the soil
    /// vectors
    vector<int> has_soil_data_index;
    
    /// The top depth (in soil, in g/cm^2) of a soil sample
    /// This vector is indexed into the other vectors with valid_soil_samples
    vector<double> soil_top_effective_depth;
    
    /// The thickness (in g/cm^2) of a soil sample
    /// This vector is indexed into the other vectors with valid_soil_samples
    vector<double> soil_effective_thickness;
    
    /// a vector holding the latitude of the samples
    vector<double> latitude;
    
    /// a vector holding the longitude of the samples
    vector<double> longitude;
    
    /// a vector holding the UTM of the samples
    vector<double> UTM_easting;
    
    /// a vector holding the UTM northing of the samples
    vector<double> UTM_northing;
    
    /// the nuclide. Only options are Be10 and 26Al. 
    vector<string> nuclide;
    
    /// a string holding the standardisation
    vector<string> standardisation;

    /// a vector holding the concetration, in atoms per gram, of the sample
    vector<double> Concentration;
    
    /// a vector holding the AMS uncertainty of the sample, in atoms per gram
    vector<double> Concentration_uncertainty; 
    
    /// a vector holding the concetration, in atoms per gram, of the sample
    /// this is for data before standardisation
    vector<double> Concentration_unstandardised;
    
    /// a vector holding the AMS uncertainty of the sample, in atoms per gram
    /// this is for data before standardisation
    vector<double> Concentration_uncertainty_unstandardised; 
    
    /// a vector of vectors holding the results of the cosmogenic analysis
    vector< vector<double> > erosion_rate_results;
    
    /// a standardisation map for Be10
    map<string,double> standards_Be10;
    
    /// a standardisation map for Al26
    map<string,double> standards_Al26;
    
    /// a string vecvec for holding the DEM names involved in the analysis
    ///  elements are
    ///  [0] = DEM name
    ///  [1] = snow_shield raster
    ///  [2] = self shield raster
    ///  [3] = topo shield raster
    ///  When there is no raster, this is "NULL"
    vector< vector<string> > DEM_names_vecvec;
    
    /// A vector holding the parameters for snow and self shielding if it is
    ///  a single value.
    vector< vector<double> > snow_self_topo_shielding_params;
    
    /// The minimum slope for the fill function
    float min_slope;
  
    /// The number of pixels for a channel
    int source_threshold;
    
    /// The number of pixels over which to search for a channel
    int search_radius_nodes;
  
    /// The minimum stream order to be considered a channel for cosmo sampling
    int threshold_stream_order;

    /// The azimuth step for topographic shielding calculations
    int theta_step;
  
    /// The inclination step for topographic sheilding calculations
    int phi_step;

    /// an uncertainty parameter which was superceded by the new
    /// error analyses but I have been too lazy to remove it. Does nothing 
    double prod_uncert_factor;
    
    /// The muon production scaling. Options are "Braucher", "Granger" and "Schaller"
    string Muon_scaling;       

    /// the atmospheric data is in the folder with the driver_functions, 
    /// but can be changed if necessary.
    string path_to_atmospheric_data;
  
    /// The boundary conditions for the flow info object
    vector<string> boundary_conditions;
    
    //---------------Flags for writing files---------------------
    /// Write toposheild rasters if they don't exist
    bool write_TopoShield_raster;
    
    /// Write a LSDIndexRaster with the basins
    bool write_basin_index_raster;
    
    /// Write the shielding and scaling rasters
    bool write_full_scaling_rasters;
    
    //-----------------Information used in cosmogenic calculators---------------
    /// This contains data with all sorts of scaling parameters
    /// for calculation of erosion rates using other calculators
    map<string, map<int,double> > MapOfProdAndScaling;

    //---------------basin metrics---------------------
    /// The Basin Relief
    vector<double> MBS;

  private:
  
    /// @brief the empty create function
    void create();
    
    /// @brief the create function used when calling a file
    /// @param the path to the file. Must have the "/" at the end
    /// @param  file_prefice the prefix (without extension) of the parameter files
    /// @author SMM
    /// @date 02/02/2015
    void create(string path, string file_prefix);
    
    /// @brief This loads data from a text file
    /// @detail The data columns are:
    ///  column[0]: sample_name (NO SPACES OR COMMAS!!)
    ///  column[1]: latitude (decimal degrees)
    ///  column[2]: longitude (decimal degrees)
    ///  column[3]: Nuclide (Be10 or Al26)
    ///  column[4]: Nuclide concentration (atoms per gram)
    ///  column[6]: Nuclide uncertainty (atoms per gram)
    ///  column[7]: standardisation
    /// @param filename the name of the file WITH PATH AND EXTENSION
    /// @author SMM
    /// @date 09/02/2015
    void load_txt_cosmo_data(string filename);
    
    /// @brief This loads data from a csv file
    /// @detail The data columns are:
    ///  column[0]: sample_name (NO SPACES OR COMMAS!!)
    ///  column[1]: latitude (decimal degrees)
    ///  column[2]: longitude (decimal degrees)
    ///  column[3]: Nuclide (Be10 or Al26)
    ///  column[4]: Nuclide concentration (atoms per gram)
    ///  column[6]: Nuclide uncertainty (atoms per gram)
    ///  column[7]: standardisation
    /// @param filename the name of the file WITH PATH AND EXTENSION
    /// @author SMM
    /// @date 09/02/2015
    void load_csv_cosmo_data(string filename);
};

#endif
