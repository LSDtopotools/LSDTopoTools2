//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDCosmoData.cpp
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
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <ctype.h>
#include <sstream>
#include <algorithm> 
#include <vector>
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDCosmoData.hpp"
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDBasin.hpp"
#include "LSDRasterInfo.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDCosmoData_CPP
#define LSDCosmoData_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The default create function. Doesn't do anything
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::create()
{
  Muon_scaling = "Braucher";
}

void LSDCosmoData::create(string path_name, string param_name_prefix)
{
  // Fix the path if it doesn't have a trailing slash
  path_name = FixPath(path_name);
  
  cout << "I am loading an LSDCosmoData object" << endl;
  cout << "Your pathname is: " << path_name << endl;
  cout << "You parameter name prefix is: " << param_name_prefix << endl;
  
  // remove control characters from these strings
  path_name.erase(remove_if(path_name.begin(), path_name.end(), ::iscntrl), path_name.end());
  param_name_prefix.erase(remove_if(param_name_prefix.begin(), param_name_prefix.end(), ::iscntrl), param_name_prefix.end());
  
  string combined_filename = path_name+param_name_prefix;
  cout << "The combined filename is: " << combined_filename << endl;

  // set the names of the data members
  path = path_name;
  param_name = param_name_prefix;

  // set up extensions to the files
  string crn_ext = "_CRNData.csv";
  string files_ext = "_CRNRasters.csv";
  string csv_ext = "csv";
  string params_ext = ".CRNParam";
  string outfile_ext = ".CRNParamReport";
  string soil_ext = "_CRNSoilInfo.csv";
  
  // get the filenames to open
  string crn_fname = path_name+param_name_prefix+crn_ext;
  string Rasters_fname =  path_name+param_name_prefix+files_ext;
  string parameters_fname = path_name+param_name_prefix+params_ext;
  string soil_fname = path_name+param_name_prefix+soil_ext;
  
  cout << "The files I am checking are: " << endl;
  cout << crn_fname << endl;
  cout << Rasters_fname << endl;
  cout << parameters_fname << endl;

  // check the parameter files
  check_files(crn_fname,Rasters_fname,parameters_fname,soil_fname);
  
  // Initiate default parameters
  initiate_default_parameters();

  // load the data. 
  cout << "Loading CRN data" << endl;
  load_cosmogenic_data(crn_fname, csv_ext);
  
  cout << "Number of samples: " << N_samples << endl;
  
  // get the correct number of samples in the results vec
  vector<double> for_calculator_empty(N_samples);
  
  // set the map of scaling parameters with an empty map
  map< string, map<int,double> > empty_map;
  MapOfProdAndScaling = empty_map;
  
  cout << "Loading file structures" << endl;
  load_DEM_and_shielding_filenames_csv(Rasters_fname);
  
  cout << "Loading parameters" << endl;
  load_parameters(parameters_fname);
  
  // check the input parameters
  check_parameter_values();
  
  // check the rasters
  cout << "Hold on while I check the rasters. I don't want to eat dirty rasters!" << endl;
  check_rasters();
  cout << "The rasters seem okay. Mmmm, rasters." << endl;
  
  // see if there is additional soil information
  cout << "Checking soil information from the file:" << soil_fname <<  endl;
  load_soil_info(soil_fname);
  
  cout << "I'm all done loading the CRN data!" << endl;
  cout << "=================================================" << endl << endl << endl;
  
  // now print all the data into a report
  string outfile_name = path_name+param_name_prefix+outfile_ext;
  print_all_data_parameters_and_filestructures(outfile_name);
  
}
//==============================================================================


//==============================================================================
// This initiates default parameters in case they are not in the parameter file
//==============================================================================
void LSDCosmoData::initiate_default_parameters()
{
  cout << "I am initiating default parameters, these will be overwittern by the paramfile...";
  // flags for writing files
  write_TopoShield_raster = true;
  write_basin_index_raster = true;
  
  // Set the parameters
  // The default slope parameter for filling. Do not change. 
  min_slope = 0.0001;
  
  //cout << "1" << endl;
  
  // parameters for making stream networks and looking for channels
  source_threshold = 12;
  search_radius_nodes = 1;
  threshold_stream_order = 1;

  // the values of theta and phi step are based on testing by S. Grieve 
  // Note that Codilian recommends 5,5 but 10,15 leads to minimal errors
  theta_step = 30;
  phi_step = 30;

  // some environment variables
  prod_uncert_factor = 1;          // this is a legacy parameter.
  

  //cout << "Default theta and phi steps: " << theta_step << " " << phi_step << endl;

  string test_scaling = "Braucher";
  
  //cout << "Test scaling is: " << test_scaling << endl;
  //cout << "default production_uncer_factor: " <<   prod_uncert_factor << endl;
  
  Muon_scaling = test_scaling;       // default muon scaling

  //cout << "default muon scaling: " << Muon_scaling << endl;

  // the atmospheric data is in the folder with the driver_functions
  string path_to_atmospheric_data = "./";
  
  //cout << "default path_to_atmospheric_data: " << path_to_atmospheric_data  << endl;
  
  // a boundary condition for the flow info object
  vector<string> boundary_conditionst(4);
  boundary_conditionst[0] = "n";
  boundary_conditionst[1] = "n";
  boundary_conditionst[2] = "n";
  boundary_conditionst[3] = "n";
  boundary_conditions = boundary_conditionst;
  
  cout << "done." << endl;
  
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function that actually loads the file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_cosmogenic_data(string filename, string filetype)
{
  if(filetype != "csv" && filetype != "txt")
  {
    cout << "LSDCosmoData line 111 You have not selected a valid filetype!" << endl;
    cout << "Options are csv and txt" << endl;
    cout << "Defaulting to txt" << endl;
    filetype = "txt";
  }
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: The file" << filename
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  // now populate the standardisation maps
  standards_Be10["07KNSTD"] = 1.0;
  standards_Be10["KNSTD"] = 0.9042;
  standards_Be10["NIST_Certified"] = 1.0425;
  standards_Be10["LLNL31000"] = 0.8761;
  standards_Be10["LLNL10000"] = 0.9042;
  standards_Be10["LLNL3000"] = 0.8644;
  standards_Be10["LLNL1000"] = 0.9313;
  standards_Be10["LLNL300"] = 0.8562;
  standards_Be10["NIST_30000"] = 0.9313;
  standards_Be10["NIST_30200"] = 0.9251;
  standards_Be10["NIST_30300"] = 0.9221;
  standards_Be10["NIST_30600"] = 0.9130;
  standards_Be10["NIST_27900"] = 1.0;
  standards_Be10["S555"] = 0.9124;
  standards_Be10["S2007"] = 0.9124;
  standards_Be10["BEST433"] = 0.9124;
  standards_Be10["BEST433N"] = 1.0;
  standards_Be10["S555N"] = 1.0;
  standards_Be10["S2007N"] = 1.0;

  standards_Al26["KNSTD"] = 1.0;
  standards_Al26["ZAL94"] = 0.9134;
  standards_Al26["SMAL11"] = 1.021;
  standards_Al26["0"] = 1.0;
  standards_Al26["ZAL94N"] = 1.0;
  standards_Al26["ASTER"] = 1.021;
  standards_Al26["Z92-0222"] = 1.0;

  if(filetype == "csv")
  {
    load_csv_cosmo_data(filename);
  }
  else
  {
    load_txt_cosmo_data(filename);
  }
  
  // now loop through the data, getting the standardised concentrations
  N_samples = int(sample_name.size());
  
  // create the vec vec for holding sample results
  vector<double> empty_vec;
  vector< vector<double> > result_vecvec;
  
  for(int i = 0; i<N_samples; i++)
  {
    // populate the results vector
    result_vecvec.push_back(empty_vec);
  
    // read the nuclide
    if(nuclide[i] == "Be10")
    {
      // check to see if standard exists
      if(standards_Be10.find(standardisation[i]) == standards_Be10.end())
      {
        cout << "Standardisation not found, assuming 07KNSTD" << endl;
        standardisation[i] = "07KNSTD";
      }
      
      // adjust the 10Be concentration
      Concentration.push_back(Concentration_unstandardised[i]*
                                standards_Be10[ standardisation[i] ]);
      Concentration_uncertainty.push_back(Concentration_uncertainty_unstandardised[i]*
                                standards_Be10[ standardisation[i] ]); 
    }
    else if (nuclide[i] == "Al26")
    {
      // check to see if standard exists
      if(standards_Al26.find(standardisation[i]) == standards_Al26.end())
      {
        cout << "Standardisation not found, assuming KNSTD" << endl;
        standardisation[i] = "KNSTD";
      }
      
      // adjust the 26Al concentration
      Concentration.push_back(Concentration_unstandardised[i]*
                                standards_Al26[ standardisation[i] ]);
      Concentration_uncertainty.push_back(Concentration_uncertainty_unstandardised[i]*
                                standards_Al26[ standardisation[i] ]);   
    }
    else
    {
      cout << "You haven't selected a valid nuclide, choices are Be10 and Al26" << endl;
      cout << "You chose: " << nuclide[i] << ", defaulting to Be10 and 07KNSTD" << endl;
      nuclide[i] = "Be10";
      Concentration.push_back(Concentration_unstandardised[i]);
      Concentration_uncertainty.push_back(Concentration_uncertainty_unstandardised[i]);
    }

  }
  erosion_rate_results = result_vecvec;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a csv file of cosmogenic data
// The first row is a header
// the following rows contain the data. 
// The data columns are:
//  column[0]: sample_name (NO SPACES OR COMMAS!!)
//  column[1]: latitude (decimal degrees)
//  column[2]: longitude (decimal degrees)
//  column[3]: Nuclide (Be10 or Al26)
//  column[4]: Nuclide concentration (atoms per gram)
//  column[6]: Nuclide uncertainty (atoms per gram)
//  column[7]: standardisation
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_csv_cosmo_data(string filename)
{
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv cosmo data file, but the file" << filename
         << "doesn't exist; LINE 245 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  // initiate temporary vectors
  vector<string> temp_sample_name;
  vector<double> temp_latitude;
  vector<double> temp_longitude;
  vector<string> temp_nuclide;
  vector<double> temp_Concentration_unstandardised;
  vector<double> temp_Concentration_uncertainty_unstandardised;
  vector<string> temp_standardisation;

  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  
  // discard the first line
  getline(ifs, line_from_file);
  
  // now loop through the rest of the lines, getting the data. 
  while( getline(ifs, line_from_file))
  {
    // reset the string vec
    this_string_vec = empty_string_vec;
    
    // create a stringstream
    stringstream ss(line_from_file);
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
    if ( int(this_string_vec.size()) < 7)
    {
      cout << "Hey there, I am trying to load your cosmo data but you seem not to have" << endl;
      cout << "enough columns in your file. I am ignoring a line" << endl;
    }
    else
    {
      // now convert the data
      //cout << "Getting sample name: " <<  this_string_vec[0] << endl;
      
      // let the user know about offending underscores, and replace them
      string s = this_string_vec[0];
      string uscore = "_";
      size_t found = s.find(uscore);
      if (found!=std::string::npos)
      {
        cout << "I found an underscore in the sample name. Replacing with a dash." <<endl;
        replace( s.begin(), s.end(), '_', '-');
        cout << "New sample name is: " << s << endl;
      }

      temp_sample_name.push_back( s );
      temp_latitude.push_back( atof( this_string_vec[1].c_str() ) );
      temp_longitude.push_back( atof(this_string_vec[2].c_str() ) );
      temp_nuclide.push_back( this_string_vec[3] );
      temp_Concentration_unstandardised.push_back( atof(this_string_vec[4].c_str() ) );
      temp_Concentration_uncertainty_unstandardised.push_back( atof(this_string_vec[5].c_str() ) );
      temp_standardisation.push_back( this_string_vec[6] );
    }
  
  }
  
  //for(int i = 0; i<temp_sample_name.size(); i++)
  //{
  //  cout << endl << endl << "==============Sample "<< i << "=======" << endl;
  //  cout << temp_sample_name[i] << endl;
  //  cout << temp_latitude[i] << endl;
  //  cout << temp_longitude[i] << endl;
  //  cout << temp_nuclide[i] << endl;  
  //  cout << temp_Concentration_unstandardised[i] << endl;
  //  cout << temp_Concentration_uncertainty_unstandardised[i] << endl;
  //  cout << temp_standardisation[i] << endl;
  //}
  
  
  // now update the data members
  sample_name = temp_sample_name;
  latitude = temp_latitude;
  longitude = temp_longitude;
  nuclide = temp_nuclide;
  Concentration_unstandardised = temp_Concentration_unstandardised;
  Concentration_uncertainty_unstandardised = 
                         temp_Concentration_uncertainty_unstandardised;
  standardisation = temp_standardisation;
  
  //cout << "Done!" << endl;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This loads information from a soil file
// The file contains several columns:
// sampleID,sample_top_depth,sample_bottom_depth,density
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_soil_info(string filename)
{

  // initiate temporary vectors
  vector<string> temp_sample_name;
  vector<double> temp_top_depth;
  vector<double> temp_sample_thickness;
  vector<double> temp_density;
  vector<int> temp_has_soil_data_index(N_samples,-1);

  vector<int> valid_soil_samples;
  vector<double> valid_top_eff_depth;
  vector<double> valid_sample_eff_thickness;
  double this_top_eff_depth;
  double this_eff_thick;
  
  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nThere isn't any soil information. If you think there should be, check your filename." << endl;
    cout << "You looked for: " << filename << endl;
  }
  else
  {
    // discard the first line
    getline(ifs, line_from_file);
  
    // now loop through the rest of the lines, getting the data. 
    while( getline(ifs, line_from_file))
    {
      // reset the string vec
      this_string_vec = empty_string_vec;
    
      // create a stringstream
      stringstream ss(line_from_file);
    
      while( ss.good() )
      {
        string substr;
        getline( ss, substr, ',' );
      
        // remove the spaces
        substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
        // remove control characters
        substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
        // add the string to the string vec
        this_string_vec.push_back( substr );
      }
    
      //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
      if ( int(this_string_vec.size()) < 4)
      {
        cout << "Hey there, I am trying to load your soil data but you seem not to have" << endl;
        cout << "enough columns in your file. I am ignoring a line" << endl;
      }
      else
      {
        // now convert the data
        //cout << "Getting sample name: " <<  this_string_vec[0] << endl;
      
        // let the user know about offending underscores, and replace them
        string s = this_string_vec[0];
        string uscore = "_";
        size_t found = s.find(uscore);
        if (found!=std::string::npos)
        {
          cout << "I found an underscore in the sample name. Replacing with a dash." <<endl;
          replace( s.begin(), s.end(), '_', '-');
          cout << "New sample name is: " << s << endl;
        }

        temp_sample_name.push_back( s );
        temp_top_depth.push_back( atof( this_string_vec[1].c_str() ) );
        temp_sample_thickness.push_back( atof(this_string_vec[2].c_str() ) );
        temp_density.push_back( atof( this_string_vec[3].c_str() ) );
      }  
    }
    

    int n_soil_samples = int(temp_sample_name.size());
    // now you need to check if the soil samples match up with and recorded samples
    int soil_sample_counter = 0;
    for (int i = 0; i< n_soil_samples; i++)
    {
      //cout << "I have a soil sample named " << temp_sample_name[i] 
      //     << " and I'm trying to find it in the sample list" << endl;
    
      for (int samp = 0; samp < N_samples; samp++)
      {
        
        if(sample_name[samp] == temp_sample_name[i])
        {
          cout << "You found a valid sample: " << sample_name[samp] << endl;
          valid_soil_samples.push_back(samp);
          
          // calculate the effective depths and thicknesses
          // effective depths are in g/cm^2
          // these assume that the density is in kg/m^3 and the
          // thicknesses are in cm
          this_top_eff_depth = temp_density[i]*temp_top_depth[i]*0.001;
          valid_top_eff_depth.push_back( this_top_eff_depth );
          
          this_eff_thick = temp_density[i]*temp_sample_thickness[i]*0.001;
          valid_sample_eff_thickness.push_back(this_eff_thick);
          
          
          cout << "This effective depth is: " << this_top_eff_depth << endl;
          cout << "The eff thickness is: " <<  this_eff_thick << endl;
          
          // have and index that points from the CRN sample to the vectors
          // of soil data
          temp_has_soil_data_index[samp] = soil_sample_counter;
          soil_sample_counter++;
        }
        else
        {
          //cout << "No, this sample is called: " << sample_name[samp] << endl;
        }
      }
    }

    // now update the data members
    soil_sample_index = valid_soil_samples;
    soil_top_effective_depth = valid_top_eff_depth;
    soil_effective_thickness = valid_sample_eff_thickness;
    
  }
  
  has_soil_data_index = temp_has_soil_data_index; 
  
  if (soil_sample_index.size() != 0)
  {
    cout << "The number of soil samples are: " << soil_sample_index.size() << " "
         << soil_top_effective_depth.size() << " " <<  soil_effective_thickness.size() << endl;
  }
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads parameters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_parameters(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  string parameter, value, lower;

  if( infile.fail() )
  {
    cout << "Parameter file: " << filename << " does not exist." << endl;
    cout << "Using default parameters." << endl;
  }

  // now ingest parameters
  while (infile.good())
  {
    parse_line(infile, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
    {
      lower[i] = tolower(parameter[i]);
    }

    cout << "parameter is: " << lower << " and value is: " << value << endl;

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);
    
    // now set the parameters
    if (lower == "min_slope")
    {
      min_slope = atof(value.c_str());
    }
    else if (lower == "source_threshold")
    {
      source_threshold = atoi(value.c_str());
    }
    else if (lower == "search_radius_nodes")
    {
      search_radius_nodes = atoi(value.c_str());
    }
    else if (lower == "threshold_stream_order")
    {
      threshold_stream_order = atoi(value.c_str());
    }
    else if (lower == "theta_step")
    {
      theta_step = atoi(value.c_str());
    }
    else if (lower == "phi_step")
    {
      phi_step = atoi(value.c_str());
    }
    else if (lower == "path_to_atmospheric_data")
    {
      path_to_atmospheric_data = value;
    }
    else if (lower == "muon_scaling")
    {
      if(value.find("braucher") == 0 || value.find("Braucher") == 0)
      {
        Muon_scaling = "Braucher";
        cout << "You have selected Braucher scaling" << endl;
      }
      else if(value.find("granger") == 0 || value.find("Granger") == 0)
      {
        Muon_scaling = "Granger";
        cout << "You have selected Granger scaling" << endl;
      }
      else if(value.find("Schaller") == 0 || value.find("schaller") == 0)
      {
        Muon_scaling = "Schaller";
        cout << "You have selected Schaller scaling" << endl;
      }
      else if(value.find("newCRONUS") == 0 || value.find("newCRONUS") == 0)
      {
        Muon_scaling = "newCRONUS";
        cout << "You have selected the new CRONUS scaling" << endl;
      }
      else
      {
        Muon_scaling = "Braucher";
        cout << "You have not selected a valid scaling, defaulting to Braucher" << endl;
      }
    }
    else if (lower == "write_toposhield_raster")
    {
      if(value.find("true") == 0 || value.find("True") == 0)
      {
        write_TopoShield_raster = true;
      }
      else if (value.find("false") == 0 || value.find("False") == 0)
      {
        write_TopoShield_raster = false;
      }
      else
      {
        write_TopoShield_raster = true;
        cout << "You have not selected a valid toposhield write. Defaulting to true." << endl;
      }
    }
    else if (lower == "write_basin_index_raster")
    {
      if(value.find("true") == 0 || value.find("True") == 0)
      {
        write_basin_index_raster = true;
      }
      else if (value.find("false") == 0 || value.find("False") == 0)
      {
        write_basin_index_raster = false;
      }
      else
      {
        write_basin_index_raster = true;
        cout << "You have not selected a valid toposhield write. Defaulting to true." << endl;
      }
    }
    else if (lower == "write_full_scaling_rasters")
    {
      if(value.find("true") == 0 || value.find("True") == 0)
      {
        write_full_scaling_rasters = true;
      }
      else if (value.find("false") == 0 || value.find("False") == 0)
      {
        write_full_scaling_rasters = false;
      }
      else
      {
        write_full_scaling_rasters = true;
        cout << "You have not selected a valid toposhield write. Defaulting to true." << endl;
      }
    }
    else
    {
      cout << "Line " << __LINE__ << ": No parameter '"
           << parameter << "' expected.\n\t> Check spelling." << endl;
    }
  }
  infile.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads the filenames 
// for the DEMs or single value parameters.
// It reads the DEM, either the snow shield raster or a single value
// the self shield raster name or a single value, and the topo shield raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_DEM_and_shielding_filenames_csv(string filename)
{
  //cout << "Getting filenames" << endl;
  
  // this vecvec holds data for determining the dem, the snow shielding
  // the self shielding and the topographic sheilding
  vector< vector<string> > temp_DEM_names_vecvec;
  
  // a string for null values
  string null_str = "NULL";

  // this vecvec holds information about snow, self and toposhielding. 
  vector< vector<double> > temp_snow_self_topo_shielding_params; 

  // a vector of strings for holding the DEM names.
  // Elements without a DEM get a null value
  vector<string> null_string_vec(4,null_str);
  vector<string> this_snow_self_shield_names;
  
  // a vector for holding parameter values. 
  vector<double> empty_snow_self(2,0);
  vector<double> this_snow_self;
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv filenames file, but the file" << filename
         << "doesn't exist; LINE 348 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  
  // now loop through the rest of the lines, getting the data. 
  while( getline(ifs, line_from_file))
  {
    // reset the string vec
    this_string_vec = empty_string_vec;
    
    // create a stringstream
    stringstream ss(line_from_file);
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    // reset the vectors for this line
    this_snow_self_shield_names = null_string_vec;
    this_snow_self = empty_snow_self;
    
    // now we need to see how many data elements we have
    int n_strings_in_line =  int(this_string_vec.size());
    
    // the can hold 1, 2, 3 or 4 elements. 
    // If it holds 1, it only has the name of the DEM
    // If it has 2, it has name of DEM and snow shielding
    // If it has 3, it has name of DEM and self shielding
    // If it has 4, the last element must be the name of the toposhield raster
    //  which has been precalculated
    if(n_strings_in_line == 1)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
    }
    else if(n_strings_in_line == 2)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
      
      // test if the second element is a number
      string second = this_string_vec[1];
      if(isdigit(second[0]))
      {
        this_snow_self[0] = atof(this_string_vec[1].c_str());
      }
      else
      {
        this_snow_self_shield_names[1] = this_string_vec[1];
      }
      
    }
    else if(n_strings_in_line == 3)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
      
      // test if the second element is a number
      string second = this_string_vec[1];
      if(isdigit(second[0]))
      {
        this_snow_self[0] = atof(this_string_vec[1].c_str());
      }
      else
      {
        this_snow_self_shield_names[1] = this_string_vec[1];
      }

      // test if the third element is a number
      string third = this_string_vec[2];
      if(isdigit(third[0]))
      {
        this_snow_self[1] = atof(this_string_vec[2].c_str());
      }
      else
      {
        this_snow_self_shield_names[2] = this_string_vec[2];
      }            
    }    
    else if(n_strings_in_line == 4)
    {
      this_snow_self_shield_names[0] = this_string_vec[0];
      
      // test if the second element is a number
      string second = this_string_vec[1];
      if(isdigit(second[0]))
      {
        this_snow_self[0] = atof(this_string_vec[1].c_str());
      }
      else
      {
        this_snow_self_shield_names[1] = this_string_vec[1];
      }

      // test if the third element is a number
      string third = this_string_vec[2];
      if(isdigit(third[0]))
      {
        this_snow_self[1] = atof(this_string_vec[2].c_str());
      }
      else
      {
        this_snow_self_shield_names[2] = this_string_vec[2];
      }                  
      
      this_snow_self_shield_names[3] = this_string_vec[3];
    }
    else 
    {
      cout << "LSDCosmoData line 383, your cosmo data file names has the" << endl;
      cout << "wrong number of elements on this line. " << endl;
      cout << "I am only going to keep the DEM name" << endl;
      this_snow_self_shield_names[0] = this_string_vec[0];
    } 

    // push the parameters back into the vecvecs
    temp_DEM_names_vecvec.push_back(this_snow_self_shield_names);
    temp_snow_self_topo_shielding_params.push_back(this_snow_self);
  } 
  
  // just check if there are no empty entries (thanks to hidden newline characters)
  bool back_is_not_okay = 1;
  //cout << "HEY FAT ALBERT" << endl;
  while(back_is_not_okay)
  {
    int NDEMS = int(temp_DEM_names_vecvec.size());
    string last_DEM_name = temp_DEM_names_vecvec[NDEMS-1][0];
    cout << "Last_DEM_name:" << last_DEM_name << endl;
    if(last_DEM_name == "")
    {
      cout << "Oh, BUGGER, who put a bunch of empty lines at the back of this file!" << endl;
      cout << "Microsoft is probably to blame if you opened in excel." << endl;
      cout << "I'm getting rid of this one." <<endl;
      temp_DEM_names_vecvec.pop_back();
    }
    else
    {
      back_is_not_okay = 0; 
    }
    
    
  }
  

  DEM_names_vecvec = temp_DEM_names_vecvec;
  snow_self_topo_shielding_params = temp_snow_self_topo_shielding_params;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a text file of cosmogenic data
// The first row is a header
// the following rows contain the data. 
// The data columns are:
//  column[0]: sample_name (NO SPACES OR COMMAS!!)
//  column[1]: latitude (decimal degrees)
//  column[2]: longitude (decimal degrees)
//  column[3]: Nuclide (Be10 or Al26)
//  column[4]: Nuclide concentration (atoms per gram)
//  column[6]: Nuclide uncertainty (atoms per gram)
//  column[7]: standardisation
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::load_txt_cosmo_data(string filename)
{
  cout << "Opening text file: " << filename << endl;
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: The file" << filename
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  // initiate temporary vectors
  vector<string> temp_sample_name;
  vector<double> temp_latitude;
  vector<double> temp_longitude;
  vector<string> temp_nuclide;
  vector<double> temp_Concentration_unstandardised;
  vector<double> temp_Concentration_uncertainty_unstandardised;
  vector<string> temp_standardisation;

  // initiate the string to hold the file
  string line_from_file;

  // now load the file. 
  double this_lat;
  double this_long;
  double this_conc;
  double this_uncert;
  
  string this_sample_name;
  string this_standard;
  string this_nuclide;
  
  // discard the first line
  ifs >> this_sample_name >> this_sample_name  >> this_sample_name >> this_sample_name 
             >> this_sample_name >> this_sample_name >> this_sample_name;
  
  while( ifs >> this_sample_name >> this_lat >> this_long >> this_nuclide 
             >> this_conc >> this_uncert >> this_standard)
  {
    temp_sample_name.push_back(this_sample_name);
    temp_latitude.push_back(this_lat);
    temp_longitude.push_back(this_long);
    temp_nuclide.push_back(this_nuclide);
    temp_Concentration_unstandardised.push_back(this_conc);
    temp_Concentration_uncertainty_unstandardised.push_back(this_uncert);
    temp_standardisation.push_back(this_standard);
  }

  // now update the data members
  sample_name = temp_sample_name;
  latitude = temp_latitude;
  longitude = temp_longitude;
  nuclide = temp_nuclide;
  Concentration_unstandardised = temp_Concentration_unstandardised;
  Concentration_uncertainty_unstandardised = 
                         temp_Concentration_uncertainty_unstandardised;
  standardisation = temp_standardisation;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks to see if all the files exist
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::check_files(string crn_fname,string Rasters_fname,string parameters_fname,
                               string soil_fname)
{
  cout << "=======================================================" << endl;
  cout << "I am checking on your files to see if they exist." << endl;
  // make sure the filename works
  ifstream ifsc(crn_fname.c_str());
  if( ifsc.fail() )
  {
    cout << "FATAL ERROR: The crn file: " << crn_fname
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "The CRN file exists." << endl;
  }

  // make sure the filename works
  ifstream ifsr(Rasters_fname.c_str());
  if( ifsr.fail() )
  {
    cout << "FATAL ERROR: The rasters file: " << Rasters_fname
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "The Rasters file exists." << endl;
  }
  
  // make sure the filename works
  ifstream ifsp(parameters_fname.c_str());
  if( ifsp.fail() )
  {
    cout << "FATAL ERROR: The parameter file: " << parameters_fname
         << "doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "The parameter file exists." << endl;
  }
  
  // make sure the filename works
  ifstream ifss(soil_fname.c_str());
  if( ifss.fail() )
  {
    cout << "WARNING: The soil file: " << soil_fname
         << "doesn't exist." << endl
         << "This is fine as long as you are not doing soil analysis." << endl;
  }
  else
  {
    cout << "There is a soil file." << endl;
  }
  cout << "=====================================================" << endl;    
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the DEM names
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDCosmoData::get_DEM_fnames()
{
  int N_DEMS = int(DEM_names_vecvec.size());
  vector<string> DEM_fnames;

  // loop through the filename vecvec, extracting the DEM names
  for (int i = 0; i<N_DEMS; i++)
  {
    vector<string> vnames = DEM_names_vecvec[i];
    DEM_fnames.push_back(vnames[0]);
  }
  return DEM_fnames;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the names of the snow and self shielding rasters
// If none exists, the name is NULL
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDCosmoData::get_Snow_fnames()
{

  vector<string> Snow_fnames;

  // now print the results to screen
  int n_files = int(DEM_names_vecvec.size());
  
  if(n_files != 0)
  {
    for(int row = 0; row<n_files; row++)
    {
      vector<string> temp_stringvec = DEM_names_vecvec[row];
      Snow_fnames.push_back(temp_stringvec[1]);
    }      
  }
  else
  {
    cout << "LSDCosmoData::print_file_structures_to_screen; files have not been loaded!" << endl;
  }
  return Snow_fnames;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the names of the snow and self shielding rasters
// If none exists, the name is NULL
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDCosmoData::get_Self_fnames()
{

  vector<string> Self_fnames;

  // now print the results to screen
  int n_files = int(DEM_names_vecvec.size());
  
  if(n_files != 0)
  {
    for(int row = 0; row<n_files; row++)
    {
      vector<string> temp_stringvec = DEM_names_vecvec[row];
      Self_fnames.push_back(temp_stringvec[2]);
    }      
  }
  else
  {
    cout << "LSDCosmoData::print_file_structures_to_screen; files have not been loaded!" << endl;
  }
  return Self_fnames;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function checks the parameters of the object and sends warnings if they
// seem incorrect
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::check_parameter_values()
{
  if (min_slope <=0)
  {
    cout << "Your minimum slope is negative! Changing to default 0.0001" << endl;
    min_slope = 0.0001;
  }
  else if (min_slope > 0.05)
  {
    cout << "You minimum slope is a bit high at: " << min_slope << endl;
    cout << "Are you sure about this? Check your parameter file. "<< endl; 
  }

  if (source_threshold  <= 0)
  {
    cout << "Your source threshold is too small! Changing to default 10" << endl;
    source_threshold = 10;
  }
  
  if (search_radius_nodes  <= 0)
  {
    cout << "Your earch_radius_nodes is too small! Changing to default 1" << endl;
    search_radius_nodes= 1;
  }
  
  if (threshold_stream_order  <= 0)
  {
    cout << "Your threshold_stream_order is too small! Changing to default 1" << endl;
    threshold_stream_order = 1;
  }
  
  if (Muon_scaling != "Braucher" && Muon_scaling != "Granger" && 
      Muon_scaling != "Schaller" && Muon_scaling != "newCRONUS")
  {
    cout << "You have not seleceted a valid scaling. Defaulting to Braucher" << endl;
    Muon_scaling = "Braucher";
  }
  
  if (prod_uncert_factor != 1.0)
  {
    cout << "Production uncertainty factor must be 1." << endl;
    prod_uncert_factor = 1;
  }
  
  // Check the atmospheric data files
  string filename = "NCEP2.bin";
  filename = path_to_atmospheric_data+filename;
  //cout << "Loading mean sea level, file is: " << endl << filename << endl;

  ifstream ifs_data(filename.c_str(), ios::in | ios::binary);
  if( ifs_data.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist. You need to put the atmospheric data" << endl
         << "In the correct path" << endl;
    exit(EXIT_FAILURE);
  }  
  
  // now the data with levels
  filename = "NCEP_hgt.bin";
  filename = path_to_atmospheric_data+filename;
  //cout << "Loading hgt, file is: " << endl << filename << endl;

  ifstream ifs_data2(filename.c_str(), ios::in | ios::binary);
  if( ifs_data2.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist. You need to put the atmospheric data" << endl
         << "In the correct path" << endl;
    exit(EXIT_FAILURE);
  }
  
  // now check the phi and theta values. These must be a factor of 360 and 90, respectively
  int temp_step;
  if (360%theta_step != 0)
  {
    temp_step = theta_step-1;
    while(360%temp_step != 0)
    {
      temp_step--;
    }
    cout << "Theta was not a factor of 360, changing from: " << theta_step << endl;
    theta_step = temp_step;
    cout << " to: " << theta_step << endl;
  }
  if (90%phi_step != 0)
  {
    temp_step = phi_step-1;
    while(90%temp_step != 0)
    {
      temp_step--;
    }
    cout << "Phi was not a factor of 90, changing from: " << phi_step << endl;
    phi_step = temp_step;
    cout << " to: " << phi_step << endl;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis checks the rasters for georeferencing and scaling
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::check_rasters()
{
  // loop through the lines in the files, checking to see if the 
  // georeferencing is equivalent
  string bil_ext = "bil";
  string null_str = "NULL";
  int N_DEMS = int(DEM_names_vecvec.size());
  for(int iDEM = 0; iDEM<N_DEMS; iDEM++)
  {
    // get the names from this DEM
    vector<string>  DEM_names_vec = DEM_names_vecvec[iDEM];
    cout << "Checking rasters, this raster is: " << DEM_names_vec[0] << endl;
    
    // get the info from the DEM
    LSDRasterInfo DEM_info(DEM_names_vec[0],bil_ext);
    
    // now compare with the other DEMs
    // first snow shielding
    if(DEM_names_vec[1] != null_str)
    {
      LSDRasterInfo SnowShield_info(DEM_names_vec[1],bil_ext);
      if( SnowShield_info!=DEM_info)
      {
        cout << "Snow shielding raster not the same shape and size as the DEM!" << endl;
        cout << "Setting snow shielding to NULL" << endl;
        DEM_names_vec[1] = null_str;
      }
    }
    // now self shielding
    if(DEM_names_vec[2] != null_str)
    {
      LSDRasterInfo SelfShield_info(DEM_names_vec[2],bil_ext);
      if( SelfShield_info!=DEM_info)
      {
        cout << "Self shielding raster not the same shape and size as the DEM!" << endl;
        cout << "Setting self shielding to NULL" << endl;
        DEM_names_vec[2] = null_str;
      }
    }
    // now toposhielding
    if(DEM_names_vec[3] != null_str)
    {
      LSDRasterInfo TopoShield_info(DEM_names_vec[3],bil_ext);
      if( TopoShield_info!=DEM_info)
      {
        cout << "Topo shielding raster not the same shape and size as the DEM!" << endl;
        cout << "Setting topo shielding to NULL" << endl;
        DEM_names_vec[3] = null_str;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Spawn clipped basins
// This function attempts to speed up shielding calculations by
// finding basins and then clipping each basin to a DEM. It sucks up a whole bunch
// of disk space but is faster than doing entire DEMs
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDCosmoData::spawn_clipped_basins(string DEM_fname, int padding_pixels)
{

  vector<string> new_dem_names;
  
  cout << "\n\n\n====================================================\n";
  cout << "Spawning basins" << endl;

  // now find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // Load the DEM
  string DEM_bil_extension = "bil";
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
  //cout << "Filled raster" << endl;
  
  // remove the seas
  //cout << endl << endl << endl << "--------------------" << endl;
  //cout << "Removing seas" << endl;
  //filled_raster.remove_seas();
  //string JB = DEM_fname+ "_Updated";
  //filled_raster.write_raster(JB,"bil");
  
  
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  //cout << "Got flow info" << endl;

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  //cout << "Got contributing pixels" << endl;
  
  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);
  //cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  //cout << "Got junction network" << endl;
  
  // Now convert the data into this UTM zone
  convert_to_UTM(filled_raster);

  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);

  int n_valid_points = int(valid_cosmo_points.size());

  //========================
  // LOOPING THROUGH BASINS
  //========================
  // now loop through the valid points, getting the cosmo data 
  cout << "I am trying to spawn basins, found " << n_valid_points << " valid points." << endl;
  for(int samp = 0; samp<n_valid_points; samp++)
  {
     
    cout << "Valid point is: " << valid_cosmo_points[samp] << " X: " 
         << UTM_easting[valid_cosmo_points[samp]] << " Y: "
         << UTM_northing[valid_cosmo_points[samp]] << endl;
    cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
         << snapped_junction_indices[samp] << " and sample name is: " << sample_name[valid_cosmo_points[samp]] << endl;
         
    cout << "Getting basin" << endl;
    LSDBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork);
    
    cout << "Now getting the basin raster" << endl;  
    LSDRaster BasinRaster = thisBasin.TrimPaddedRasterToBasin(padding_pixels, 
                                         FlowInfo,filled_raster);
    cout << "Formatting the path" << endl;
    
    
    string DEM_path = DEM_fname+ "/";
    
    //string DEM_newpath = ReformatPath(DEMpath);
    string DEMnewname = DEM_path+ "SpawnedBasin_"+sample_name[valid_cosmo_points[samp]];
    
    cout << "Writing a new basin to: " << DEMnewname << endl;
        
    BasinRaster.write_raster(DEMnewname,"bil");
    new_dem_names.push_back(DEMnewname);
    
    
    //string JB = DEM_path+ "JustBasin_"+itoa(valid_cosmo_points[samp]);
    //LSDRaster JustBasinRaster = thisBasin.write_raster_data_to_LSDRaster(filled_raster, FlowInfo);
    //JustBasinRaster.write_raster(JB,"bil");
    
  }
  cout << "Finished with this DEM!" << endl 
       << "I found " << new_dem_names.size() << " new  basins!" << endl << "----------------" << endl; 
  return new_dem_names;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::BasinSpawnerMaster(string path, string prefix, int padding_pixels)
{
  // first get the names of the DEMs:
  vector<string> dfnames = get_DEM_fnames();
  vector<string> snow_fnames = get_Snow_fnames();
  vector<string> self_fnames = get_Self_fnames();
  string DEM_bil_extension = "bil";
  
  // holds shielding parameters
  vector<double> this_ss_shielding;
  
  // this holds the information from the spawned basins
  // open the file for this dem
  ofstream new_CRNRasters_csv;
  string new_csv_name =  path+prefix+"_Spawned_CRNRasters.csv";
  cout << endl << endl << "=======================" << endl << "SPAWNING" << endl;
  cout << "New csv name is: " << new_csv_name << endl;
  new_CRNRasters_csv.open(new_csv_name.c_str());
  
  int n_DEMs = dfnames.size();
  for(int i = 0; i<n_DEMs; i++)
  {
    
    cout << "Looking at DEM: " <<  dfnames[i] << endl;
    vector<string> new_dem_names = spawn_clipped_basins(dfnames[i], padding_pixels);
    
    if (new_dem_names.size() != 0)
    {
      // start buidling the string that will go into the new file
      vector<string> new_DEMs_vec;
      vector<string> new_snow_vec;
      vector<string> new_self_vec;
      
      // get the snow_self_shielding_params
      this_ss_shielding = snow_self_topo_shielding_params[i];

      // loop though the new rasters
      int NSpawn = new_dem_names.size();
      cout << "Number of new DEMs is: " <<  NSpawn << endl;
      for(int ns = 0; ns<NSpawn; ns++)
      {
        cout << "======================" << endl;
        cout << "On spawn number " << ns << endl;
        
        new_DEMs_vec.push_back(new_dem_names[ns]);

        // now check to see if there is a shielding raster
        if(snow_fnames[i] == "NULL" && self_fnames[i] == "NULL")
        {
          cout << "No snow or self shielding rasters, adding constant values" << endl;
          new_snow_vec.push_back(dtoa(this_ss_shielding[0]));
          new_self_vec.push_back(dtoa(this_ss_shielding[1]));
        }
        else
        {
          // at least one of these rasters must be clipped. So load the basin
          cout << "I need to clip a snow or shielding raster!" << endl;
          cout << "I am clipping from the raster: " << new_dem_names[ns] << endl;
          LSDRaster clipped_raster(new_dem_names[ns],DEM_bil_extension);
        
          cout << "First lets look for the snow raster. It is: " << snow_fnames[i] << endl;
        
          if(snow_fnames[i] == "NULL")
          {
            new_snow_vec.push_back(dtoa(this_ss_shielding[0]));
          }
          else
          {
            cout << "There is a snow raster here, it is the file: " << snow_fnames[i] << endl;
          
            // load the shielding dem
            LSDRaster Snow_Raster(snow_fnames[i],DEM_bil_extension);
            
            cout << "I loaded the snow raster" << endl;

            //string DEMnewname = new_dem_names[ns]+"_hummahumma";

            //clipped_raster.write_raster(DEMnewname,"bil");
            
            // clip this raster
            LSDRaster Snow_clip = Snow_Raster.clip_to_smaller_raster(clipped_raster);
            
            cout << "I clipped that raster for you, buddy" << endl;
            
            // save this raster
            Snow_clip.write_raster(new_dem_names[ns]+"_snowclip",DEM_bil_extension);
            
            cout << "I wrote that snow raster" << endl;
            
            // add the name to the list
            new_snow_vec.push_back(new_dem_names[ns]+"_snowclip");
            
            cout << "The snow name has been pushed back" << endl;
            
          }
          
          cout << "Now lets look for the self raster. It is: " << self_fnames[i] << endl;
          
          if(self_fnames[i] == "NULL")
          {
            cout << "There isn't a self shielding raster, I'm pushing back a constant value" << endl;
            new_self_vec.push_back(dtoa(this_ss_shielding[1]));
          }
          else
          {
            // load the shielding dem
            LSDRaster Self_Raster(self_fnames[i],DEM_bil_extension);
            
            // clip this raster
            LSDRaster Self_clip = Self_Raster.clip_to_smaller_raster(clipped_raster);
            
            // save this raster
            Self_clip.write_raster(new_dem_names[ns]+"_selfclip",DEM_bil_extension);
            
            // add the name to the list
            new_self_vec.push_back(new_dem_names[ns]+"_selfclip");
          }
        }
        cout << "I've finished clipping those rasters" << endl << endl;
        
      }      // end if statement to see if there are shielding rasters
    
      int n_new_dems = new_DEMs_vec.size();
      cout << "The number of names is: " << n_new_dems << endl;
      cout << "New snow vecs: " << new_snow_vec.size() << endl;
      cout << "New self vecs: " << new_self_vec.size() << endl;
      
      
      for(int ndm = 0; ndm < n_new_dems; ndm++)
      {
        new_CRNRasters_csv << new_dem_names[ndm] + "," + new_snow_vec[ndm] + "," +
                              new_self_vec[ndm] << endl;
      }
      
    }
  }
  new_CRNRasters_csv.close();
  
  // now spawn a new cosmodata file with the correct prefix
  string new_prefix = prefix+"_Spawned";
  print_renamed_cosmo_data(path, new_prefix);
  print_renamed_parameter_data(path, new_prefix);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Shielding calculations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::RunShielding(string path, string prefix)
{
  vector<string> dfnames = get_DEM_fnames();
  vector<string> shield_names;
  string DEM_Format = "bil";

  //loop over each DEM and generate the shielding raster.
  int n_DEMs = dfnames.size();
  for(int i = 0; i<n_DEMs; i++)
  {
    
    cout << "Processing " << i+1 << " of " << n_DEMs << endl;
    
    //load the raster
    LSDRaster ForShield(dfnames[i], DEM_Format);

    // run shielding
    LSDRaster Shielded = ForShield.TopographicShielding(theta_step, phi_step);
    
    //write the shielding raster to the working directory
    Shielded.write_raster((dfnames[i]+"_SH"),DEM_Format);
    shield_names.push_back(dfnames[i]+"_SH");
  }
  
  // now rewrite the parameter file
  string rasters_name = path+prefix+"_CRNRasters.csv";
  cout << endl << endl << endl << "---------------------------------" << endl;
  cout << "Rasters name is: " << rasters_name << endl;
  
  
  // make sure the filename works
  ifstream ifs(rasters_name.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv filenames file, but the file" << rasters_name
         << "doesn't exist; LINE 348 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  
  string line_from_file;
  vector<string> lines;
  vector<string> this_string_vec;
  vector<string> empty_string_vec;
  
  string zero = "0";
  string comma = ",";
  
  // now loop through the rest of the lines, getting the data. 
  int i = 0;
  while( getline(ifs, line_from_file))
  {
    // reset the string vec
    this_string_vec = empty_string_vec;
    
    // create a stringstream
    stringstream ss(line_from_file);
    
    // read in the elements in the line
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      this_string_vec.push_back( substr );
    }
    
    // a temporary string for holding the line
    string temp_string;
    
    // now check on the size of the string
    int n_params = this_string_vec.size();
    
    // this is logic to update the file depending on how many arguments are already in the line
    if(n_params ==1 )
    {
      temp_string = this_string_vec[0]+comma+zero+comma+zero+comma+shield_names[i];
    }
    else if (n_params == 2)
    {
      temp_string = this_string_vec[0]+comma+this_string_vec[1]+comma+zero+comma+shield_names[i];
    }
    else if (n_params >= 3)
    {
      temp_string = this_string_vec[0]+comma+this_string_vec[1]+comma+this_string_vec[2]+comma+shield_names[i];
    }
    else
    {
      cout << "\nFATAL ERROR: Trying to load csv filenames file, but the file" << rasters_name
         << "doesn't have the correct number of elements in the line; LINE 1375 LSDCosmoData" << endl;
      exit(EXIT_FAILURE);    
    }

    // add the line to the list of lines
    cout << "line is: " << temp_string << " and shield name is: " << shield_names[i] << endl;
    lines.push_back(temp_string);
    i++;
  }
  ifs.close();
  
  // now write the data
  string rasters_out_name = path+prefix+"_Shield_CRNRasters.csv";
  ofstream rasters_out;
  rasters_out.open(rasters_out_name.c_str());
  
  int n_files = int(lines.size());
  cout << "N_files is: " << n_files << endl;
  for(int i = 0; i<n_files; i++)
  {
    rasters_out << lines[i] << endl;
  }
  rasters_out.close();
  
  // now spawn a new cosmodata file with the correct prefix
  string new_prefix = prefix+"_Shield";
  print_renamed_cosmo_data(path, new_prefix);
  print_renamed_parameter_data(path, new_prefix);
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints a production raster for a basin
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDCosmoData::calculate_production_raster(LSDRaster& Elevation_Data,
                                               string path_to_atmospheric_data)
{
  int NRows = Elevation_Data.get_NRows();
  int NCols = Elevation_Data.get_NCols();
  float NDV =  Elevation_Data.get_NoDataValue();
  
  Array2D<float> Production(NRows,NCols,NDV);

  // variables for converting location and elevation
  double this_elevation, this_pressure;

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;
  
  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();
  
  // a function for scaling stone production, defaults to 1
  double Fsp = 1.0;
  
  // the latitude and longitude
  double lat,longitude;
  
  // decalre converter object
  LSDCoordinateConverterLLandUTM Converter;
  
  for (int row = 0; row < NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      
      //exclude NDV from average
      if (Elevation_Data.get_data_element(row,col) != NDV)
      {
        // To get pressure, first get the lat and long
        Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
        //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);
      
        // now the elevation
        this_elevation = Elevation_Data.get_data_element(row,col);
      
        // now the pressure
        this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude), 
                                        double(this_elevation));

        // get the production
        Production[row][col] = LSDCRNP.stone2000sp(lat,this_pressure, Fsp);
      }
    }
  }

  float XMinimum = Elevation_Data.get_XMinimum();
  float YMinimum = Elevation_Data.get_YMinimum();
  float DataResolution = Elevation_Data.get_DataResolution();
  map<string,string> GeoReferencingStrings = Elevation_Data.get_GeoReferencingStrings();

  LSDRaster Production_raster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NDV, Production,GeoReferencingStrings);
  return Production_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function prints the data to screen
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_data_to_screen()
{
  cout << endl;
  cout << "==========================================================" << endl;
  cout << "PRINTING COSMO DATA HELD IN LSDCOSMODATA OBJECT" << endl;
  cout << "Sample_name\tLatitude\tLongitude\tNuclide\tConcentration\tUncertainty\tscaling\n";
  for(int i = 0; i<N_samples; i++)
  {
    cout << sample_name[i] << "\t"  << latitude[i] << "\t"  << longitude[i] << "\t"
         << nuclide[i] << "\t"  << Concentration_unstandardised[i] << "\t"
         << Concentration_uncertainty_unstandardised[i] << "\t"
         << standardisation[i] << "\n";
  }
  cout << "==========================================================" << endl;
  cout << endl << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function prints the data to screen
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_renamed_cosmo_data(string path, string prefix)
{
  ofstream new_CRN_data;
  
  // make sure you don't lose any information
  new_CRN_data.precision(8);
  

  
  string new_CRN_data_name = path+prefix+"_CRNData.csv";
  new_CRN_data.open(new_CRN_data_name.c_str());

  cout << "I am printing a renamed cosmo data file for you." << endl;
  cout << "The filename is: " <<  new_CRN_data_name << endl;

  new_CRN_data << "Sample_name,Latitude,Longitude,Nuclide,Concentration,Uncertainty,Standardisation" << endl;
  
  for (int i= 0; i<N_samples; i++)
  {
    new_CRN_data << sample_name[i] << "," << latitude[i] << "," << longitude[i] << ","
                 << nuclide[i] << "," << Concentration_unstandardised[i] 
                 << "," << Concentration_uncertainty_unstandardised[i] << ","
                 << standardisation[i] << endl;
  }
  new_CRN_data.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_renamed_parameter_data(string path, string prefix)
{
  ofstream new_param_data;
  string new_param_data_name = path+prefix+".CRNParam";
  
  cout << "I am printing a renamed cosmo data file for you." << endl;
  cout << "The filename is: " <<  new_param_data_name << endl;
  
  new_param_data.open(new_param_data_name.c_str());
  new_param_data << "min_slope: " << min_slope << endl;
  new_param_data << "source_threshold: " << source_threshold << endl;
  new_param_data << "search_radius_nodes: " << search_radius_nodes << endl;
  new_param_data << "threshold_stream_order: " << threshold_stream_order << endl;
  new_param_data << "theta_step: " << theta_step << endl;
  new_param_data << "phi_step: " << phi_step << endl; 
  new_param_data << "Muon_scaling: " << Muon_scaling << endl;
  if (write_basin_index_raster)
  {
    new_param_data << "write_basin_index_raster: True" << endl;
  }
  else
  {
    new_param_data << "write_basin_index_raster: False" << endl;
  }
  if (write_TopoShield_raster)
  {
    new_param_data << "write_toposhield_raster: True" << endl;
  }
  else
  {
    new_param_data << "write_toposhield_raster: False" << endl;
  }    
  if (write_full_scaling_rasters)
  {
    new_param_data << "write_full_scaling_rasters: True" << endl;
  }
  else
  {
    new_param_data << "write_full_scaling_rasters: False" << endl;
  }  
  new_param_data.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the file structures to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_file_structures_to_screen()
{

  // now print the results to screen
  int n_files = int(DEM_names_vecvec.size());
  
  if(n_files != 0)
  {
    for(int row = 0; row<n_files; row++)
    {
      vector<string> temp_stringvec = DEM_names_vecvec[row];
      vector<double> temp_params = snow_self_topo_shielding_params[row];

      cout << "------------------------------" << endl << "sample" << row << endl;
      cout << "DEM: " << temp_stringvec[0] << endl;
      if (temp_stringvec[1] == "NULL")
      {
        cout << "Snow shielding constant depth (g/cm^2): " << temp_params[0] << endl;
      }
      else
      {
        cout << "Snow shield raster: " << temp_stringvec[1] << endl;
      }
      if (temp_stringvec[2] == "NULL")
      {
        cout << "Self shielding constant depth (g/cm^2): " << temp_params[1] << endl;
      }
      else
      {
        cout << "Self shield raster: " << temp_stringvec[2] << endl;
      }      
      cout << "Topo shield raster: " << temp_stringvec[3] << endl;

    }      
  }
  else
  {
    cout << "LSDCosmoData::print_file_structures_to_screen; files have not been loaded!" << endl;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis function prints all the parameters, file structures and cosmo
// data to one file that can be used for reconstructing calculations later
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_all_data_parameters_and_filestructures(string outfilename)
{
  ofstream outfile;
  outfile.open(outfilename.c_str());
  
  // first the cosmo data, comma seperated:
  for(int i = 0; i<N_samples; i++)
  {
    outfile << sample_name[i] << ","  << latitude[i] << ","  << longitude[i] << ","
         << nuclide[i] << ","  << Concentration_unstandardised[i] << ","
         << Concentration_uncertainty_unstandardised[i] << ","
         << standardisation[i] << "\n";
  }
  
  // now the parameters
  outfile << "----------------------------------------------" << endl << endl;
  outfile << "min_slope: " << min_slope << endl;
  outfile << "source_threshold: " << source_threshold << endl;
  outfile << "search_radius_nodes: " << search_radius_nodes << endl;
  outfile << "threshold_stream_order: " << threshold_stream_order << endl;
  outfile << "theta_step: " << theta_step << endl;
  outfile << "phi_step: " << phi_step << endl; 
  outfile << "Muon_scaling: " << Muon_scaling << endl;
  outfile << "----------------------------------------------" << endl << endl;
  
  // now the file structures
  // now print the results to screen
  int n_files = int(DEM_names_vecvec.size());
  
  if(n_files != 0)
  {
    for(int row = 0; row<n_files; row++)
    {
      vector<string> temp_stringvec = DEM_names_vecvec[row];
      vector<double> temp_params = snow_self_topo_shielding_params[row];

      outfile << "------------------------------" << endl << "sample" << row << endl;
      outfile << "DEM: " << temp_stringvec[0] << endl;
      if (temp_stringvec[1] == "NULL")
      {
        outfile << "Snow shielding constant depth (g/cm^2): " << temp_params[0] << endl;
      }
      else
      {
        outfile << "Snow shield raster: " << temp_stringvec[1] << endl;
      }
      if (temp_stringvec[2] == "NULL")
      {
        outfile << "Self shielding constant depth (g/cm^2): " << temp_params[1] << endl;
      }
      else
      {
        outfile << "Self shield raster: " << temp_stringvec[2] << endl;
      }      
      outfile << "Topo shield raster: " << temp_stringvec[3] << endl;

    }      
  }
  else
  {
    outfile << "LSDCosmoData::print_file_structures_to_screen; files have not been loaded!" << endl;
  }
  
  outfile.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the results to screen
// simple is for external, muon and production uncertainty only
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_simple_results_to_screen(double rho)
{
    
    cout << "=======================================================" << endl;
    cout << "Printing results of the analysis" << endl;
    for (int i = 0; i<N_samples; i++)
    {
      
      
      // don't print the results unless they exist
      if (int(erosion_rate_results[i].size()) > 0)
      {
        // get the results from this sample
        vector<double> erate_analysis = erosion_rate_results[i];
        
        cout << "-----------------------------------------------------" << endl;
        cout << "Sample " << sample_name[i] << " , a " << nuclide[i] << " sample" << endl;
        cout << "latitude:\t" << latitude[i] << "\tlongitude:\t" << longitude[i] << endl;
        cout << "Concentration: " << Concentration[i] << " +/- " 
             << Concentration_uncertainty[i] << " atoms/g" << endl;
        cout << "Erate is: " 
             <<  erate_analysis[0] << " g/cm^2/yr" << endl;
        cout << "The erosion rate for rho = "<< rho << " is: " << endl
             << erate_analysis[0]*10/rho << " m/yr and " 
             << erate_analysis[0]*1e6/rho << " cm/kyr" << endl;
        cout << "Ext uncert: " << erate_analysis[1] << " muon uncert: " << erate_analysis[2]
             << " production uncert: " <<  erate_analysis[3] << " g/cm^2/yr"  << endl;
        cout << "Total uncertainty in g/cm^2/yr: " << erate_analysis[4] << endl;
      }
    }
    cout << "=======================================================" << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function converts data to UTM coordinates. You have to tell it 
// the UTM zone
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::convert_to_UTM(int UTM_zone)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;
  
  // set up some temporary vectors
  vector<double> this_UTMN(N_samples,0);
  vector<double> this_UTME(N_samples,0);
  
  double this_Northing;
  double this_Easting;
  
  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    cout << "Converting point " << i << " to UTM." << endl;
    Converter.LLtoUTM_ForceZone(eId, latitude[i], longitude[i], 
                      this_Northing, this_Easting, UTM_zone);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
  }
  
  UTM_easting = this_UTME;
  UTM_northing = this_UTMN;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function converts data to UTM coordinates. It determines the UTM from
// a raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::convert_to_UTM(LSDRaster& Raster)
{
  // get the UTM zone from the raster
  int UTM_zone;
  bool is_North;
  Raster.get_UTM_information(UTM_zone, is_North);
  
  // convert the data
  convert_to_UTM(UTM_zone);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the determination of cosmogenic erosion rates
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::basic_cosmogenic_analysis(string DEM_fname)
{
  // Load the DEM
  string DEM_bil_extension = "bil";
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
    
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);

  // get the topographic shielding
  LSDRaster TopoShield = filled_raster.TopographicShielding(theta_step, phi_step);

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  
  // first, convert the data into this UTM zone
  convert_to_UTM(filled_raster);
  
  // now find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);

  // after this operation the three int vectors valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertainties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }
  
  // some environment variables
  double prod_uncert_factor = 1;          // this is a legacy parameter.
  string Muon_scaling = "Braucher";       // defaul muon scaling
  //bool data_from_outlet_only = true;      // this needs to be turned off later!!!
  
  
  // some temporary doubles to hold the nuclide concentrations
  double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
  double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin
  double test_N, test_dN;       // concentration and uncertainty of the nuclide in basin
  
  // now loop through the valid points, getting the cosmo data 
  for(int samp = 0; samp<n_valid_points; samp++)
  {
    if( valid_nuclide_names[samp] == "Be10")
    {
      test_N10 = valid_concentrations[samp];
      test_dN10 = valid_concentration_uncertainties[samp];
      test_N26 = 1e9;
      test_dN26 = 0;
      test_N = test_N10;
      test_dN = test_dN10;
    }
    else if( valid_nuclide_names[samp] == "Al26")
    {
      test_N10 = 1e9;
      test_dN10 = 0;
      test_N26 = valid_concentrations[samp];
      test_dN26 = valid_concentration_uncertainties[samp];
      test_N = test_N26;
      test_dN = test_dN26;
    }
    else
    {
      cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
      cout << "Defaulting to Be10" << endl;
      valid_nuclide_names[samp] = "Be10";
      test_N10 = valid_concentrations[samp];
      test_dN10 = valid_concentration_uncertainties[samp];
      test_N26 = 1e9;
      test_dN26 = 0;
      test_N = test_N10;
      test_dN = test_dN10;
    }
    
    cout << "Valid point is: " << valid_cosmo_points[samp] << " X: " 
         << UTM_easting[valid_cosmo_points[samp]] << " Y: "
         << UTM_northing[valid_cosmo_points[samp]] << endl;
    cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
         << snapped_junction_indices[samp] << endl;
    LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                            test_N10,test_dN10, test_N26,test_dN26);
    
    // populate the scaling vectors
    thisBasin.populate_scaling_vectors(FlowInfo, filled_raster, TopoShield,
                                       path_to_atmospheric_data);
    
    // now do the analysis
    vector<double> erate_analysis = thisBasin.full_CRN_erosion_analysis(test_N, 
                                        valid_nuclide_names[samp], test_dN, 
                                        prod_uncert_factor, Muon_scaling);
    
    erosion_rate_results[ valid_cosmo_points[samp] ] = erate_analysis;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the determination of cosmogenic erosion rates 
// It allows nesting by accounting for erosion rates known from other 
//  samples (for example and upstream cosmogenic data point)
// THIS IS FOR A SINGLE RASTER
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::full_shielding_cosmogenic_analysis_nested(vector<string> Raster_names,
                            vector<double> CRN_params, 
                            LSDRaster& known_eff_erosion)
{

  cout << endl << endl << "================================================" << endl;
  cout << "Looking for basins in raster: " << Raster_names[0] << endl << endl;

  // some parameters for printing the basins, if that is called for
  int basin_number;
  int basin_pixel_area;
  map<int,int> basin_area_map;

  //cout << "\n\n\n====================================================\n";
  //cout << "Performing analysis with snow and self shielding" << endl;

  // now find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // Load the DEM
  string DEM_bil_extension = "bil";
  //string fill_ext = "_fill";
  string DEM_fname = Raster_names[0];
  //cout << "Loading raster: " << DEM_fname << endl;
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  topo_test.remove_seas();
  
  // check to see if the rasters are the same
  LSDRasterInfo DEM_info(topo_test);
  LSDRasterInfo Erosion_info(known_eff_erosion);
  if (DEM_info != Erosion_info)
  {
    cout << "LSDCosmoData::full_shielding_cosmogenic_analysis_nested ERROR!" << endl;
    cout << "The erosion raster and DEM are not the same dimesions." <<endl;
    exit(EXIT_SUCCESS);
  }
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
  //cout << "Filled raster" << endl;
  
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  //cout << "Got flow info" << endl;

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  //cout << "Got contributing pixels" << endl;
  
  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);
  //cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  //cout << "Got junction network" << endl;
  
  // Now convert the data into this UTM zone
  convert_to_UTM(filled_raster);

  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  cout << "Getting snapped basins" << endl;
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);
  cout << "Snapped points" << endl;
  
  
  // after this operation the three int vectors: valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertanties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }
  cout << "Got valid points, there are " << n_valid_points << " of them." << endl;
  
  
  // Initiate pointers to the rasters
  LSDRaster Topographic_shielding;
  LSDRaster Snow_shielding;
  LSDRaster Self_shielding;
  
  
  // a flag for writing the initial basin index
  bool written_inital_basin_index = false;
  
  // now, IF there are valid points, go on to the rest of the analysis
  if (valid_nuclide_names.size() != 0)
  {
    // if you need to write the inital raster, do so 
    LSDIndexRaster BasinIndex;   // initiate and empty raster

    
    // first check if topographic shielding raster exists
    cout << "Looking for toposhield raster, name is: " << Raster_names[3] << endl;
    if( Raster_names[3] != "NULL")
    {
      LSDRaster T_shield(Raster_names[3], DEM_bil_extension);
      Topographic_shielding = T_shield;
    }
    else
    {
      // get the topographic shielding
      cout << "Starting topographic shielding" << endl;
      LSDRaster T_shield = filled_raster.TopographicShielding(theta_step, phi_step);
      Topographic_shielding = T_shield;
      
      if(write_TopoShield_raster)
      {
        string TShield_name = Raster_names[0]+"_SH";
        T_shield.write_raster(TShield_name,DEM_bil_extension);
      }
      
    }
    
    // Now check if snow and self shielding rasters exist
    bool have_snow_raster;
    bool have_self_raster;
    double constant_snow_depth = 0;
    double constant_self_depth = 0;
    if (Raster_names[1] != "NULL")
    {
      cout << "LSDCosmoData, line 2058: Loading the snow shielding raster, " 
           << Raster_names[1] << ".bil" <<  endl;
      LSDRaster Snow_shield(Raster_names[1], DEM_bil_extension);
      Snow_shielding = Snow_shield;
      have_snow_raster = true;
    }
    else
    {
      // snow shielding with a single parameter
      have_snow_raster = false;
      constant_snow_depth = CRN_params[0];
    }
    if (Raster_names[2] != "NULL")
     {
      cout << "LSDCosmoData, line 977: Loading the self shielding raster, " 
           << Raster_names[2] << ".bil" <<  endl;
      LSDRaster Self_shield(Raster_names[2], DEM_bil_extension);
      Self_shielding = Self_shield;
      have_self_raster = true;
    }
    else
    {
      // self shielding with a single parameter
      have_self_raster = false;
      constant_self_depth = CRN_params[1];
    }

    // some temporary doubles to hold the nuclide concentrations
    double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
    double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin
    double test_N, test_dN;       // concentration and uncertainty of the nuclide in basin  

    //========================
    // LOOPING THROUGH BASINS
    //========================
    // now loop through the valid points, getting the cosmo data 
    cout << "-----------------------------------------------------------" << endl;
    cout << "I found " << n_valid_points << " valid CRN basins in this raster! " << endl;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
      if( valid_nuclide_names[samp] == "Be10")
      {
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      else if( valid_nuclide_names[samp] == "Al26")
      {
        test_N10 = 1e9;
        test_dN10 = 0;
        test_N26 = valid_concentrations[samp];
        test_dN26 = valid_concentration_uncertainties[samp];
        test_N = test_N26;
        test_dN = test_dN26;
      }
      else
      {
        cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
        cout << "Defaulting to Be10" << endl;
        valid_nuclide_names[samp] = "Be10";
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      

      cout << endl << "Valid point is: " << valid_cosmo_points[samp]
           << " Sample name: " << sample_name[ valid_cosmo_points[samp] ] << " Easting: " 
           << UTM_easting[valid_cosmo_points[samp]] << " Northing: "
           << UTM_northing[valid_cosmo_points[samp]] << endl;
      cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
           << snapped_junction_indices[samp] << endl;
      cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -" << endl;
      LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                              test_N10,test_dN10, test_N26,test_dN26);

      // write the index basin if flag is set to true
      if(write_basin_index_raster)
      {
        cout << "I'm writing a basin index number for you" << endl;
      
        if (not written_inital_basin_index)
        {
        
          basin_number = valid_cosmo_points[samp];
          basin_pixel_area = thisBasin.get_NumberOfCells();
          basin_area_map[basin_number] = basin_pixel_area;
          LSDIndexRaster NewBasinIndex = 
             thisBasin.write_integer_data_to_LSDIndexRaster(basin_number, FlowInfo);
          BasinIndex = NewBasinIndex;
          written_inital_basin_index = true;
        }
        else
        {
          basin_number = valid_cosmo_points[samp];
          thisBasin.add_basin_to_LSDIndexRaster(BasinIndex, FlowInfo,
                                                basin_area_map,basin_number);
        }
      }

      // we need to scale the shielding parameters
      // now do the snow and self shielding
      if (have_snow_raster)
      {
        cout << "I've got the snow raster" << endl;
      
        if(have_self_raster)
        {
          cout << "I've also got the self raster." << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, Self_shielding);
          cout << "Done with effective depths" << endl;                           
        }
        else
        {
          cout << "No self raster, but I'm getting the effective depths" << endl;
          cout << "The constant self depth is: " <<  constant_self_depth << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, constant_self_depth);
          cout << "Done with effective depths" << endl;        
        }
      }
      else
      {
        if(have_self_raster)
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                constant_snow_depth, Self_shielding);
          cout << "Done with effective depths" << endl;                                
        }
        else
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(constant_snow_depth, 
                                             constant_self_depth);
          cout << "Done with effective depths" << endl;                                 
        }
      }


      cout << "Now I will populate the scaling vectors." << endl;
      // Now topographic shielding and production scaling
      thisBasin.populate_scaling_vectors(FlowInfo, filled_raster, 
                                         Topographic_shielding,
                                         path_to_atmospheric_data);
      cout << "The scaling vectors are populated. I am moving on to the analysis" << endl;

      // now do the analysis
      cout << "Line 2571, doing analysis" << endl;
      vector<double> erate_analysis = thisBasin.full_CRN_erosion_analysis_nested(known_eff_erosion, FlowInfo, test_N, 
                                          valid_nuclide_names[samp], test_dN, 
                                          prod_uncert_factor, Muon_scaling);
       cout << "erate: " << erate_analysis[0] << endl;
      
    
    
      // now get parameters for cosmogenic calculators
      vector<double> param_for_calc = 
          thisBasin.calculate_effective_pressures_for_calculators_nested(filled_raster,
                                            FlowInfo, path_to_atmospheric_data, 
                                            known_eff_erosion);

      cout << "Paramforcalc size: " << param_for_calc.size() << endl;              
      cout << "Getting pressures" << endl;


      // get the relief of the basin
      float R = thisBasin.CalculateBasinRange(FlowInfo, filled_raster);
      double relief = double(R);

      MapOfProdAndScaling["BasinRelief"][ valid_cosmo_points[samp] ] = relief;
      MapOfProdAndScaling["AverageProdScaling"][ valid_cosmo_points[samp] ] = param_for_calc[0];
      MapOfProdAndScaling["AverageTopoShielding"][ valid_cosmo_points[samp] ] = param_for_calc[1];
      MapOfProdAndScaling["AverageSelfShielding"][ valid_cosmo_points[samp] ] = param_for_calc[2];
      MapOfProdAndScaling["AverageSnowShielding"][ valid_cosmo_points[samp] ] = param_for_calc[3];
      MapOfProdAndScaling["AverageShielding"][ valid_cosmo_points[samp] ] =  param_for_calc[11];
      MapOfProdAndScaling["AverageCombinedScaling"][ valid_cosmo_points[samp] ] = param_for_calc[4];
      MapOfProdAndScaling["outlet_lat"][ valid_cosmo_points[samp] ] = param_for_calc[5];
      MapOfProdAndScaling["OutletPressure"][ valid_cosmo_points[samp] ] = param_for_calc[6];
      MapOfProdAndScaling["OutletEffectivePressure"][ valid_cosmo_points[samp] ] = param_for_calc[7];
      MapOfProdAndScaling["centroid_lat"][ valid_cosmo_points[samp] ] = param_for_calc[8];
      MapOfProdAndScaling["CentroidPressure"][ valid_cosmo_points[samp] ] = param_for_calc[9];
      MapOfProdAndScaling["CentroidEffectivePressure"][ valid_cosmo_points[samp] ] = param_for_calc[10];
    
      // add the erosion rate results to the holding data member
      erosion_rate_results[ valid_cosmo_points[samp] ] = erate_analysis;

      //cout << "finished adding data" << endl;

    }  // finished looping thorough basins
    
    // now print the basin LSDIndexRaster
    if(write_basin_index_raster)
    {
      string basin_ext = "_BASINS";
      string basin_fname =  DEM_fname+basin_ext;
      BasinIndex.write_raster(basin_fname, DEM_bil_extension);
    }
     
  }    // finsiehd logic for a DEM with valid points
  else
  {
    cout << "There are no valid CRN points in this raster" << endl;
  }
  cout << "==========================================" << endl << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the determination of cosmogenic erosion rates
// THIS IS FOR A SINGLE RASTER
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::full_shielding_cosmogenic_analysis(vector<string> Raster_names,
                            vector<double> CRN_params)
{

  cout << endl << endl << "================================================" << endl;
  cout << "Looking for basins in raster: " << Raster_names[0] << endl << endl;

  // some parameters for printing the basins, if that is called for
  int basin_number;
  int basin_pixel_area;
  map<int,int> basin_area_map;

  //cout << "\n\n\n====================================================\n";
  //cout << "Performing analysis with snow and self shielding" << endl;

  // now find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // Load the DEM
  string DEM_bil_extension = "bil";
  //string fill_ext = "_fill";
  string DEM_fname = Raster_names[0];
  //cout << "Loading raster: " << DEM_fname << endl;
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  topo_test.remove_seas();
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
  //cout << "Filled raster" << endl;
  
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  //cout << "Got flow info" << endl;

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  //cout << "Got contributing pixels" << endl;
  
  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);
  //cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  //cout << "Got junction network" << endl;
  
  // Now convert the data into this UTM zone
  convert_to_UTM(filled_raster);

  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  cout << "Getting snapped basins" << endl;
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);
  cout << "Snapped points" << endl;
  
  
  // after this operation the three int vectors valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertainties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }
  cout << "Got valid points, there are " << n_valid_points << " of them." << endl;
  
  
  // Initiate pointers to the rasters
  LSDRaster Topographic_shielding;
  LSDRaster Snow_shielding;
  LSDRaster Self_shielding;
  
  
  // a flag for writing the initial basin index
  bool written_inital_basin_index = false;
  
  // now, IF there are valid points, go on to the rest of the analysis
  if (valid_nuclide_names.size() != 0)
  {
    // if you need to write the inital raster, do so 
    LSDIndexRaster BasinIndex;   // initiate and empty raster

    
    // first check if topographic shielding raster exists
    cout << "Looking for toposhield raster, name is: " << Raster_names[3] << endl;
    if( Raster_names[3] != "NULL")
    {
      LSDRaster T_shield(Raster_names[3], DEM_bil_extension);
      Topographic_shielding = T_shield;
    }
    else
    {
      // get the topographic shielding
      cout << "Starting topographic shielding" << endl;
      LSDRaster T_shield = filled_raster.TopographicShielding(theta_step, phi_step);
      Topographic_shielding = T_shield;
      
      if(write_TopoShield_raster)
      {
        string TShield_name = Raster_names[0]+"_SH";
        T_shield.write_raster(TShield_name,DEM_bil_extension);
      }
      
    }
    
    // Now check if snow and self shielding rasters exist
    bool have_snow_raster;
    bool have_self_raster;
    double constant_snow_depth = 0;
    double constant_self_depth = 0;
    if (Raster_names[1] != "NULL")
    {
      cout << "LSDCosmoData, line 2058: Loading the snow shielding raster, " 
           << Raster_names[1] << ".bil" <<  endl;
      LSDRaster Snow_shield(Raster_names[1], DEM_bil_extension);
      Snow_shielding = Snow_shield;
      have_snow_raster = true;
    }
    else
    {
      // snow shielding with a single parameter
      have_snow_raster = false;
      constant_snow_depth = CRN_params[0];
    }
    if (Raster_names[2] != "NULL")
     {
      cout << "LSDCosmoData, line 977: Loading the self shielding raster, " 
           << Raster_names[2] << ".bil" <<  endl;
      LSDRaster Self_shield(Raster_names[2], DEM_bil_extension);
      Self_shielding = Self_shield;
      have_self_raster = true;
    }
    else
    {
      // self shielding with a single parameter
      have_self_raster = false;
      constant_self_depth = CRN_params[1];
    }

    // some temporary doubles to hold the nuclide concentrations
    double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
    double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin
    double test_N, test_dN;       // concentration and uncertainty of the nuclide in basin  

    //========================
    // LOOPING THROUGH BASINS
    //========================
    // now loop through the valid points, getting the cosmo data 
    cout << "-----------------------------------------------------------" << endl;
    cout << "I found " << n_valid_points << " valid CRN basins in this raster! " << endl;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
      if( valid_nuclide_names[samp] == "Be10")
      {
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      else if( valid_nuclide_names[samp] == "Al26")
      {
        test_N10 = 1e9;
        test_dN10 = 0;
        test_N26 = valid_concentrations[samp];
        test_dN26 = valid_concentration_uncertainties[samp];
        test_N = test_N26;
        test_dN = test_dN26;
      }
      else
      {
        cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
        cout << "Defaulting to Be10" << endl;
        valid_nuclide_names[samp] = "Be10";
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      

      cout << endl << "Valid point is: " << valid_cosmo_points[samp]
           << " Sample name: " << sample_name[ valid_cosmo_points[samp] ] << " Easting: " 
           << UTM_easting[valid_cosmo_points[samp]] << " Northing: "
           << UTM_northing[valid_cosmo_points[samp]] << endl;
      cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
           << snapped_junction_indices[samp] << endl;
      cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -" << endl;
      LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                              test_N10,test_dN10, test_N26,test_dN26);

      // write the index basin if flag is set to true
      if(write_basin_index_raster)
      {
        cout << "I'm writing a basin index number for you" << endl;
      
        if (not written_inital_basin_index)
        {
        
          basin_number = valid_cosmo_points[samp];
          basin_pixel_area = thisBasin.get_NumberOfCells();
          basin_area_map[basin_number] = basin_pixel_area;
          LSDIndexRaster NewBasinIndex = 
             thisBasin.write_integer_data_to_LSDIndexRaster(basin_number, FlowInfo);
          BasinIndex = NewBasinIndex;
          written_inital_basin_index = true;
        }
        else
        {
          basin_number = valid_cosmo_points[samp];
          thisBasin.add_basin_to_LSDIndexRaster(BasinIndex, FlowInfo,
                                                basin_area_map,basin_number);
        }
      }

      // we need to scale the shielding parameters
      // now do the snow and self shielding
      if (have_snow_raster)
      {
        cout << "I've got the snow raster" << endl;
      
        if(have_self_raster)
        {
          cout << "I've also got the self raster." << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, Self_shielding);
          cout << "Done with effective depths" << endl;                           
        }
        else
        {
          cout << "No self raster, but I'm getting the effective depths" << endl;
          cout << "The constant self depth is: " <<  constant_self_depth << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, constant_self_depth);
          cout << "Done with effective depths" << endl;        
        }
      }
      else
      {
        if(have_self_raster)
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                constant_snow_depth, Self_shielding);
          cout << "Done with effective depths" << endl;                                
        }
        else
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(constant_snow_depth, 
                                             constant_self_depth);
          cout << "Done with effective depths" << endl;                                 
        }
      }


      cout << "Now I will populate the scaling vectors." << endl;
      // Now topographic shielding and production scaling
      thisBasin.populate_scaling_vectors(FlowInfo, filled_raster, 
                                         Topographic_shielding,
                                         path_to_atmospheric_data);
      cout << "Done populating the scaling vectors. " << endl;

      // now do the analysis
      vector<double> erate_analysis = thisBasin.full_CRN_erosion_analysis(test_N, 
                                          valid_nuclide_names[samp], test_dN, 
                                          prod_uncert_factor, Muon_scaling);
    
      cout << "Line 2205, doing analysis" << endl;
    
    
      // now get parameters for cosmogenic calculators
      vector<double> param_for_calc = 
          thisBasin.calculate_effective_pressures_for_calculators(filled_raster,
                                            FlowInfo, path_to_atmospheric_data);
        
      cout << "Paramforcalc size: " << param_for_calc.size() << endl;              
      cout << "Getting pressures" << endl;


      // get the relief of the basin
      float R = thisBasin.CalculateBasinRange(FlowInfo, filled_raster);
      double relief = double(R);

      MapOfProdAndScaling["BasinRelief"][ valid_cosmo_points[samp] ] = relief;
      MapOfProdAndScaling["AverageProdScaling"][ valid_cosmo_points[samp] ] = param_for_calc[0];
      MapOfProdAndScaling["AverageTopoShielding"][ valid_cosmo_points[samp] ] = param_for_calc[1];
      MapOfProdAndScaling["AverageSelfShielding"][ valid_cosmo_points[samp] ] = param_for_calc[2];
      MapOfProdAndScaling["AverageSnowShielding"][ valid_cosmo_points[samp] ] = param_for_calc[3];
      MapOfProdAndScaling["AverageShielding"][ valid_cosmo_points[samp] ] =  param_for_calc[11];
      MapOfProdAndScaling["AverageCombinedScaling"][ valid_cosmo_points[samp] ] = param_for_calc[4];
      MapOfProdAndScaling["outlet_lat"][ valid_cosmo_points[samp] ] = param_for_calc[5];
      MapOfProdAndScaling["OutletPressure"][ valid_cosmo_points[samp] ] = param_for_calc[6];
      MapOfProdAndScaling["OutletEffectivePressure"][ valid_cosmo_points[samp] ] = param_for_calc[7];
      MapOfProdAndScaling["centroid_lat"][ valid_cosmo_points[samp] ] = param_for_calc[8];
      MapOfProdAndScaling["CentroidPressure"][ valid_cosmo_points[samp] ] = param_for_calc[9];
      MapOfProdAndScaling["CentroidEffectivePressure"][ valid_cosmo_points[samp] ] = param_for_calc[10];
    
      // add the erosion rate results to the holding data member
      erosion_rate_results[ valid_cosmo_points[samp] ] = erate_analysis;

      //cout << "finished adding data" << endl;

    }  // finished looping thorough basins
    
    // now print the basin LSDIndexRaster
    if(write_basin_index_raster)
    {
      string basin_ext = "_BASINS";
      string basin_fname =  DEM_fname+basin_ext;
      BasinIndex.write_raster(basin_fname, DEM_bil_extension);
    }
     
  }    // finsiehd logic for a DEM with valid points
  else
  {
    cout << "There are no valid CRN points in this raster" << endl;
  }
  cout << "==========================================" << endl << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the determination of cosmogenic erosion rates
// It is intended for spawned data where the rasters are for one sample only
// It only works with ONE raster
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::full_shielding_cosmogenic_analysis_for_spawned(vector<string> Raster_names,
                            vector<double> CRN_params)
{


  // some parameters for printing the basins, if that is called for
  int basin_number;
  int basin_pixel_area;
  map<int,int> basin_area_map;

  // We need to determine which sample number corresponds to the current raster
  // The spawning routine names the rasters so that the raster name contains 
  // the sample name after the final underscore character. So we need to find
  // the sample name
  
  int valid_samp = 0;
  // reset the string vec
  vector<string> string_vec;

  // create a stringstream
  stringstream ss(Raster_names[0]);
    
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, '_' );
      
    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
    // remove constrol characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
    // add the string to the string vec
    string_vec.push_back( substr );
  }
  
  // this extracts the sample name
  int n_elements = int(string_vec.size());
  string this_sample_name = string_vec[n_elements-1];
  
  // now find the index of the correct sample
  bool found_sample = 0;
  for (int i = 0; i< N_samples; i++)
  {
    if ( sample_name[i] == this_sample_name)
    {
      valid_samp = i;
      found_sample = 1;
    }
  }
  if(found_sample == 0)
  {
    cout << "LINE 2287 FATAL ERROR: Could not find sample!" << endl;
    exit(EXIT_SUCCESS);
  }

  cout << "Correct sample name: " << this_sample_name 
       << " the found sample name:" << sample_name[valid_samp] << endl;

  
  // Load the DEM
  string DEM_bil_extension = "bil";
  //string fill_ext = "_fill";
  string DEM_fname = Raster_names[0];
  cout << "Loading raster: " << DEM_fname << endl;
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  topo_test.remove_seas();
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
  cout << "Filled raster" << endl;
  
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  cout << "Got flow info" << endl;

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  cout << "Got contributing pixels" << endl;
  
  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);
  cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  cout << "Got junction network" << endl;
  
  // Now convert the data into this UTM zone
  convert_to_UTM(filled_raster);
  cout << "Converted to UTM" << endl;



  // Some vectors to hold the valid node
  // We need vectors rather than ints because the functions that are called
  // during this routine expect vectors
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // convert UTM vectors to float
  // Again, we only pass one value to these vectors, since that is needed by the 
  // snap to points fuction
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  fUTM_easting.push_back( float(UTM_easting[valid_samp]));
  fUTM_northing.push_back( float(UTM_northing[valid_samp]));
  cout << "Got point locations" << endl;
  
  // This snaps the junction network to the valid point
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);
  cout << "Snapped points" << endl;
  
  
  // Some parameters for debugging
  // these are the values for the single snapped point
  int single_snapped_node_index = snapped_node_indices[0];
  int single_snapped_junction_index = snapped_junction_indices[0];
  string single_snapped_nuclide_name =  nuclide[valid_samp];
  double single_snapped_concentration = Concentration[valid_samp];
  double single_snapped_uncertainty = Concentration_uncertainty[valid_samp];
  cout << endl << endl <<"=============================================" << endl;
  cout << "I've snapped the point, here are the vitalstatistix, chief:"  << endl;
  cout << "snapped NI: " << single_snapped_node_index << " JI: " << single_snapped_junction_index
       << " nuclide: " << single_snapped_nuclide_name << endl;
  cout << "Conc: " << single_snapped_concentration << " uncert: " <<  single_snapped_uncertainty << endl;
  cout << "=============================================" << endl;

  // after this operation the three int vectors valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertainties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  valid_nuclide_names.push_back(nuclide[valid_samp] );
  valid_concentrations.push_back( Concentration[valid_samp] );
  valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[valid_samp] );
  //cout << "Got valid point" << endl;
  
  // Initiate pointers to the rasters
  LSDRaster Topographic_shielding;
  LSDRaster Snow_shielding;
  LSDRaster Self_shielding;
  
  
  // a flag for writing the initial basin index
  bool written_inital_basin_index = false;
  
  int YoYoMa = int(valid_nuclide_names.size());
  
  // now, IF there are valid points, go on to the rest of the analysis
  if ( YoYoMa != 0)
  {
    // if you need to write the inital raster, do so 
    LSDIndexRaster BasinIndex;   // initiate and empty raster

    
    // first check if topographic shielding raster exists
    cout << "Looking for toposhield raster, name is: " << Raster_names[3] << endl;
    if( Raster_names[3] != "NULL")
    {
      cout << "You have given me a toposheld raster" << endl;
      LSDRaster T_shield(Raster_names[3], DEM_bil_extension);
      Topographic_shielding = T_shield;
    }
    else
    {
      // get the topographic shielding
      cout << "No toposheild raster. Starting topographic shielding." << endl;
      LSDRaster T_shield = filled_raster.TopographicShielding(theta_step, phi_step);
      Topographic_shielding = T_shield;
      
      if(write_TopoShield_raster)
      {
        string TShield_name = Raster_names[0]+"_SH";
        T_shield.write_raster(TShield_name,DEM_bil_extension);
      }
      
    }
    
    // Now check if snow and self shielding rasters exist
    bool have_snow_raster;
    bool have_self_raster;
    double constant_snow_depth = 0;
    double constant_self_depth = 0;
    if (Raster_names[1] != "NULL")
    {
      cout << "LSDCosmoData, line 2367: Loading the snow shielding raster, " 
           << Raster_names[1] << ".bil" <<  endl;
      LSDRaster Snow_shield(Raster_names[1], DEM_bil_extension);
      Snow_shielding = Snow_shield;
      have_snow_raster = true;
    }
    else
    {
      // snow shielding with a single parameter
      have_snow_raster = false;
      constant_snow_depth = CRN_params[0];
      cout << "There is no snow shielding raster. I am assuming a constant depth of " << constant_snow_depth << endl;
    }
    if (Raster_names[2] != "NULL")
     {
      cout << "LSDCosmoData, line 2381: Loading the self shielding raster, " 
           << Raster_names[2] << ".bil" <<  endl;
      LSDRaster Self_shield(Raster_names[2], DEM_bil_extension);
      Self_shielding = Self_shield;
      have_self_raster = true;
    }
    else
    {
      // self shielding with a single parameter
      have_self_raster = false;
      constant_self_depth = CRN_params[1];
      cout << "There is no self shielding raster. I am assuming a constant depth of " << constant_self_depth << endl;
    }

    // some temporary doubles to hold the nuclide concentrations
    double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
    double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin
    double test_N, test_dN;       // concentration and uncertainty of the nuclide in basin  

    //========================
    // Erosion rate in BASIN
    //========================
    // now calculate the erosion from the valid sample 
    cout << "-----------------------------------------------------------" << endl;
    if (found_sample == 0)
    {
      cout << "WARNING, FATAL ERROR: I didn't find the sample" << endl;
      exit(EXIT_SUCCESS);
    }
    else
    {
      cout << "The valid sample is: " << valid_samp << endl;
      
      // The samp index is the index into the valid_concentrations, uncertainties, 
      // etc vectors. These should all have only one element, since the routine
      // should have only found one valid sample. 
      int samp = 0;
      if( valid_nuclide_names[samp] == "Be10")
      {
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      else if( valid_nuclide_names[samp] == "Al26")
      {
        test_N10 = 1e9;
        test_dN10 = 0;
        test_N26 = valid_concentrations[samp];
        test_dN26 = valid_concentration_uncertainties[samp];
        test_N = test_N26;
        test_dN = test_dN26;
      }
      else
      {
        cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
        cout << "Defaulting to Be10" << endl;
        valid_nuclide_names[samp] = "Be10";
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      

      cout << endl << "Valid point is: " << valid_cosmo_points[samp]
           << " Sample name: " << sample_name[ valid_cosmo_points[samp] ] << " Easting: " 
           << UTM_easting[valid_cosmo_points[samp]] << " Northing: "
           << UTM_northing[valid_cosmo_points[samp]] << endl;
      cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
           << snapped_junction_indices[samp] << endl;
      cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -" << endl;
      
      // extract the basin
      LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                              test_N10,test_dN10, test_N26,test_dN26);

      // write the index basin if flag is set to true
      if(write_basin_index_raster)
      {
        if (not written_inital_basin_index)
        {
        
          basin_number = valid_cosmo_points[samp];
          basin_pixel_area = thisBasin.get_NumberOfCells();
          basin_area_map[basin_number] = basin_pixel_area;
          LSDIndexRaster NewBasinIndex = 
             thisBasin.write_integer_data_to_LSDIndexRaster(basin_number, FlowInfo);
          BasinIndex = NewBasinIndex;
          written_inital_basin_index = true;
        }
        else
        {
          basin_number = valid_cosmo_points[samp];
          thisBasin.add_basin_to_LSDIndexRaster(BasinIndex, FlowInfo,
                                                basin_area_map,basin_number);
        }
      }

      // we need to scale the shielding parameters
      // now do the snow and self shielding
      if (have_snow_raster)
      {
        cout << "I have a snow raster" << endl;
        if(have_self_raster)
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, Self_shielding);
          cout << "Done with effective depths" << endl;                       
        }
        else
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, constant_self_depth);
          cout << "Done with effective depths" << endl;                               
        }
      }
      else
      {
        cout << "I don't have a snow raster and am using a constant depth." << endl;
        if(have_self_raster)
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                constant_snow_depth, Self_shielding);
          cout << "Done with effective depths" << endl;                             
        }
        else
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(constant_snow_depth, 
                                             constant_self_depth);
          cout << "Done with effective depths" << endl;                                    
        }
      }

      // Now topographic sheidling and production scaling
      cout << "I am now populating the scaling vectors." << endl;
      thisBasin.populate_scaling_vectors(FlowInfo, filled_raster, 
                                         Topographic_shielding,
                                         path_to_atmospheric_data);

      // GET THE EROSION RATES
      cout << "I have finished with the scaling vectors and am now doing the erosion rate analysis." << endl;
      vector<double> erate_analysis = thisBasin.full_CRN_erosion_analysis(test_N, 
                                          valid_nuclide_names[samp], test_dN, 
                                          prod_uncert_factor, Muon_scaling);
      cout << "Done with the erosion rate analysis" << endl;
    
      //cout << "Line 1493, doing analysis" << endl;
    
    
      // now get parameters for cosmogenic calculators
      cout << "Let me calculate some effective pressures to test against calculators like CRONUS." << endl;
      vector<double> param_for_calc = 
          thisBasin.calculate_effective_pressures_for_calculators(filled_raster,
                                            FlowInfo, path_to_atmospheric_data);
        
      //cout << "Paramforcalc size: " << param_for_calc.size() << endl;              
      cout << "Got pressures for calculators." << endl;


      // get the relief of the basin
      float R = thisBasin.CalculateBasinRange(FlowInfo, filled_raster);
      double relief = double(R);

      MapOfProdAndScaling["BasinRelief"][ valid_samp ] = relief;
      MapOfProdAndScaling["AverageProdScaling"][ valid_samp ] = param_for_calc[0];
      MapOfProdAndScaling["AverageTopoShielding"][ valid_samp ] = param_for_calc[1];
      MapOfProdAndScaling["AverageSelfShielding"][ valid_samp ] = param_for_calc[2];
      MapOfProdAndScaling["AverageSnowShielding"][ valid_samp ] = param_for_calc[3];
      MapOfProdAndScaling["AverageShielding"][ valid_samp ] =  param_for_calc[11];
      MapOfProdAndScaling["AverageCombinedScaling"][ valid_samp ] = param_for_calc[4];
      MapOfProdAndScaling["outlet_lat"][ valid_samp ] = param_for_calc[5];
      MapOfProdAndScaling["OutletPressure"][ valid_samp ] = param_for_calc[6];
      MapOfProdAndScaling["OutletEffectivePressure"][ valid_samp ] = param_for_calc[7];
      MapOfProdAndScaling["centroid_lat"][ valid_samp ] = param_for_calc[8];
      MapOfProdAndScaling["CentroidPressure"][ valid_samp ] = param_for_calc[9];
      MapOfProdAndScaling["CentroidEffectivePressure"][ valid_samp ] = param_for_calc[10];
    
      // add the erosion rate results to the holding data member
      erosion_rate_results[ valid_samp ] = erate_analysis;

      //cout << "Added the result to the " << valid_samp << " sample." << endl;
      //cout << "finished adding data" << endl;
      
      //cout << endl << endl << "LINE 2606=-=-=-=-=-=-=" << endl;
      //cout << "Valid samp: " << valid_samp << " and fmor the vector: " 
      //     << valid_cosmo_points[samp] << endl << endl << endl;

    }  // finished looping thorough basins
    
    // now print the basin LSDIndexRaster
    if(write_basin_index_raster)
    {
      string basin_ext = "_BASINS";
      string basin_fname =  DEM_fname+basin_ext;
      BasinIndex.write_raster(basin_fname, DEM_bil_extension);
    }
     
  }    // finished logic for a DEM with valid points
  else
  {
    cout << "There are no valid CRN points in this raster" << endl;
  }
  cout << "==========================================" << endl << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function gets the erosion rate from a soil sample.
// IMPORTANT: it assumes that there is no sediment transport from upslope!
// This is most appropriate for ridgetops or from samples at the soil-saprolite bondary
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::Soil_sample_calculator(vector<string> Raster_names,
                            vector<double> CRN_params)
{

  // Load the DEM
  string DEM_bil_extension = "bil";
  //string fill_ext = "_fill";
  string DEM_fname = Raster_names[0];
  cout << "Loading raster: " << DEM_fname << endl;
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  cout << "Loaded the raster" << endl;
  
  // Some vectors to hold the valid node
  // We need vectors rather than ints because the functions that are called
  // during this routine expect vectors
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_row;       // a vector to hold the valid rows
  vector<int> snapped_col;   // a vector to hold the valid columns
  vector<double> snapped_lat;
  vector<double> snapped_long;
  
  // convert UTM vectors to float
  // Again, we only pass one value to these vectors, since that is needed by the 
  // snap to points fuction
  convert_to_UTM(topo_test);
  
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i= 0; i<N_samples; i++)
  {
    //cout << "East:" << UTM_easting[i]<< " and north: " << UTM_northing[i] << endl;
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }

  cout << "Got point locations" << endl;
  // now see if the point is in the raster
  float X_coordinate,Y_coordinate;
  for (int i = 0; i< N_samples; i++)
  {
    X_coordinate = fUTM_easting[i];
    Y_coordinate = fUTM_northing[i];
    
    // if it is in the raster, record it
    if(topo_test.check_if_point_is_in_raster(X_coordinate,Y_coordinate))
    {
      valid_cosmo_points.push_back(i);
      
      // get the row and column
      int this_row,this_col;
      topo_test.get_row_and_col_of_a_point(X_coordinate,Y_coordinate,this_row, this_col);
      snapped_row.push_back(this_row);
      snapped_col.push_back(this_col);
      snapped_lat.push_back( latitude[i] );
      snapped_long.push_back( longitude[i] );
      
      cout << "X: " << X_coordinate << " Y: " << Y_coordinate << " R:" << this_row << " C: " << this_col << endl;
    
    }
  }

  // after this operation the three int vectors valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertainties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }
  cout << "Got valid points, there are " << n_valid_points << " of them." << endl;
  
  // now the data elements that will be populated with the scaling and shielding 
  // information
  vector<double> snow_eff_depth_vec;
  vector<double> self_eff_depth_vec;
  vector<double> topo_shield_vec;
  vector<double> production_vec;
  vector<double> pressure_vec;

  // Initiate pointers to the rasters
  LSDRaster Topographic_shielding;
  LSDRaster Snow_shielding;
  LSDRaster Self_shielding;
    

  // now see if there are self and snow rasters, and extract values from them if there are
  if (snapped_row.size()!=0)
  {

    // first check if topographic shielding raster exists
    cout << "Looking for toposhield raster, name is: " << Raster_names[3] << endl;
    if( Raster_names[3] != "NULL")
    {
      LSDRaster T_shield(Raster_names[3], DEM_bil_extension);
      Topographic_shielding = T_shield;
    }
    else
    {
      // get the topographic shielding
      cout << "Starting topographic shielding" << endl;
      topo_test.remove_seas();
      
      
      LSDRaster T_shield = topo_test.TopographicShielding(theta_step, phi_step);
      Topographic_shielding = T_shield;
        
      if(write_TopoShield_raster)
      {
        string TShield_name = Raster_names[0]+"_SH";
        T_shield.write_raster(TShield_name,DEM_bil_extension);
      }
      
    }
    
    // Now check if snow and self shielding rasters exist
    bool have_snow_raster;
    bool have_self_raster;
    double constant_snow_depth = 0;
    double constant_self_depth = 0;
    if (Raster_names[1] != "NULL")
    {
      cout << "LSDCosmoData, line 2367: Loading the snow shielding raster, " 
           << Raster_names[1] << ".bil" <<  endl;
      LSDRaster Snow_shield(Raster_names[1], DEM_bil_extension);
      Snow_shielding = Snow_shield;
      have_snow_raster = true;
    }
    else
    {
      // snow shielding with a single parameter
      have_snow_raster = false;
      constant_snow_depth = CRN_params[0];
    }
    if (Raster_names[2] != "NULL")
    {
      cout << "LSDCosmoData, line 2381: Loading the self shielding raster, " 
           << Raster_names[2] << ".bil" <<  endl;
      LSDRaster Self_shield(Raster_names[2], DEM_bil_extension);
      Self_shielding = Self_shield;
      have_self_raster = true;
    }
    else
    {
      // self shielding with a single parameter
      have_self_raster = false;
      constant_self_depth = CRN_params[1];
    }
    
    // now get the snow depth and topo and self shielding from the point
    // we need to scale the sheilding parameters
    // now do the snow and self sheilding
    for(int vs = 0; vs<n_valid_points; vs++)
    {
      // first the topo shielding
      int row = snapped_row[vs];
      int col = snapped_col[vs];
      
      topo_shield_vec.push_back(Topographic_shielding.get_data_element(row,col));
      if(have_self_raster)
      {
        // get the shielding for this point
        self_eff_depth_vec.push_back(Self_shielding.get_data_element(row,col));
      }
      else
      {
        self_eff_depth_vec.push_back(constant_self_depth);
      }
      if(have_snow_raster)
      {
        snow_eff_depth_vec.push_back(Snow_shielding.get_data_element(row,col));
      }
      else
      {
        snow_eff_depth_vec.push_back(constant_snow_depth);
      }
      
      // now get the scalings
      // now create the CRN parameters object
      LSDCRNParameters LSDCRNP;
      double this_elevation, this_pressure;
  
      // get the atmospheric parameters
      LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
      LSDCRNP.set_CRONUS_data_maps();
  
      // a function for scaling stone production, defaults to 1
      double Fsp = 1.0;
      
      // now the elevation
      this_elevation = double(topo_test.get_data_element(row,col));
      cout << "This elevation " << this_elevation << endl;
      
      // now the pressure
      this_pressure = LSDCRNP.NCEPatm_2(snapped_lat[vs], snapped_long[vs], 
                                        double(this_elevation));
    
      pressure_vec.push_back(this_pressure);
    
      // now get the scaling
      production_vec.push_back(LSDCRNP.stone2000sp(snapped_lat[vs],this_pressure, Fsp));
    }

  }
  
  // now you need to check to see if there is additional information about the soil
  // samples. This will have the file data_name_CRNSoilInfo.csv
  

  // so now you have all the shielding and production in the four vecs:
  //vector<double> snow_eff_depth_vec;
  //vector<double> self_eff_depth_vec;
  //vector<double> topo_shield_vec;
  //vector<double> production_vec;
  // you then use these to get at the eroion rate. 
  for (int i = 0; i<n_valid_points; i++)
  {
    cout << "Sample: " <<  sample_name[ valid_cosmo_points[i]] << " snow: " << snow_eff_depth_vec[i] << " self: " << self_eff_depth_vec[i] 
         << endl << "    topo_shield: " << topo_shield_vec[i] << " production: " << production_vec[i] << endl;
    
    double test_N,test_dN;
    test_N = valid_concentrations[i];
    test_dN =  valid_concentration_uncertainties[i];
    
    cout << "    Nuclude: " << valid_nuclide_names[i] << " N: " << test_N << " sigma: " << test_dN << endl;
    
    
    // find any additional soil data, if it is there:
    double this_top_eff_depth = -99;
    double this_eff_thickness = -99;
    if (has_soil_data_index[i] != -1)
    {
      // check if the vectors reference back to the correct point
      int this_index = has_soil_data_index[i];
      cout << "This index is: " << this_index;
      
      if (this_index >= int(soil_sample_index.size()))
      {
        cout << "WARNING: your soil sample index is: " << this_index 
             << " which is bigger than the soil vectors." << endl;
      }
      else
      {
        this_top_eff_depth = soil_top_effective_depth[this_index];
        this_eff_thickness = soil_effective_thickness[this_index];
        cout << "The top eff depth of the soil sample is: " << this_top_eff_depth << endl;
        cout << "The effective thickness is: " << this_eff_thickness << endl;
      }
    }
    
    // if there is a surface thickness, add it to the "snow" shielding
    if(this_top_eff_depth>0)
    {
      snow_eff_depth_vec[i] = snow_eff_depth_vec[i]+this_top_eff_depth;
      cout << "The updated snow is: "  << snow_eff_depth_vec[i] << endl;
    }
    
    // if there is a sample thickness, add it to the "self" shielding
    if(this_eff_thickness>0)
    {
      self_eff_depth_vec[i] = self_eff_depth_vec[i]+this_eff_thickness;
      cout << "The updated self is: "  << self_eff_depth_vec[i] << endl;
    }
         
    vector<double> erate_analysis = full_CRN_erosion_analysis_point(test_N, 
                                          valid_nuclide_names[i], test_dN, 
                                          prod_uncert_factor, Muon_scaling,
                                          snow_eff_depth_vec[i],self_eff_depth_vec[i],
                                          topo_shield_vec[i],production_vec[i]);
                                          
    double gamma_spallation = 160;
                                          
    // get the snow shielding. Assumes spallation only
    double this_snow_shield = 1;
    if (snow_eff_depth_vec[i] != 0)
    {
      this_snow_shield = exp(-snow_eff_depth_vec[i]/gamma_spallation);
    }
      
    // now the self shiedling. Again assume spallation only
    double this_self_shield = 1.0;
    if (self_eff_depth_vec[i] != 0)
    {
      this_self_shield = gamma_spallation/self_eff_depth_vec[i]*
                             (1-exp(-self_eff_depth_vec[i]/gamma_spallation));
    }

    double this_shielding = topo_shield_vec[i]*this_self_shield*this_snow_shield;
    double this_comb_scaling = this_shielding*production_vec[i];

    MapOfProdAndScaling["BasinRelief"][ valid_cosmo_points[i] ] = 0.0;
    MapOfProdAndScaling["AverageProdScaling"][ valid_cosmo_points[i] ] = production_vec[i];
    MapOfProdAndScaling["AverageTopoShielding"][ valid_cosmo_points[i] ] = topo_shield_vec[i];
    MapOfProdAndScaling["AverageSelfShielding"][ valid_cosmo_points[i] ] = this_self_shield;
    MapOfProdAndScaling["AverageSnowShielding"][ valid_cosmo_points[i] ] = this_snow_shield;
    MapOfProdAndScaling["AverageShielding"][ valid_cosmo_points[i] ] =  this_shielding;
    MapOfProdAndScaling["AverageCombinedScaling"][ valid_cosmo_points[i] ] = this_comb_scaling;
    MapOfProdAndScaling["outlet_lat"][ valid_cosmo_points[i] ] = latitude[i];
    MapOfProdAndScaling["OutletPressure"][ valid_cosmo_points[i] ] = pressure_vec[i];
    MapOfProdAndScaling["OutletEffectivePressure"][ valid_cosmo_points[i] ] = pressure_vec[i];
    MapOfProdAndScaling["centroid_lat"][ valid_cosmo_points[i] ] = latitude[i];
    MapOfProdAndScaling["CentroidPressure"][ valid_cosmo_points[i] ] = pressure_vec[i];
    MapOfProdAndScaling["CentroidEffectivePressure"][ valid_cosmo_points[i] ] = pressure_vec[i];
    
    // add the erosion rate results to the holding data member
    erosion_rate_results[ valid_cosmo_points[i] ] = erate_analysis;

  }


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This runs the erosion rate analysis for a point in the landcape
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCosmoData::full_CRN_erosion_analysis_point(double Nuclide_conc, string Nuclide, 
                              double Nuclide_conc_err, double prod_uncert_factor,
                              string Muon_scaling, double snow_eff_depth,
                              double self_eff_depth, double topo_shield,
                              double production)
{
  // the vector for holding the erosion rates and uncertainties
  vector<double> erate_uncert_vec;
  
  double erate;                   // effective erosion rate g/cm^2/yr
  double erate_external_plus;     // effective erosion rate g/cm^2/yr for AMS uncertainty +
  double erate_external_minus;    // effective erosion rate g/cm^2/yr for AMS uncertainty -
  double dEdExternal;             // change in erosion rate for change in AMS atoms/g
  double External_uncert;         // uncertainty of effective erosion rate g/cm^2/yr for AMS
  
  
  double production_uncertainty;  // a lumped production uncertainty value. 
                                  // not generally used but needs to be passed
                                  // to the erosion finding routines as a parameter
  
  // variables for the muon uncertainty
  double average_production_rate; // The average production rate, used in uncertainty
                                  // calculations
  double erate_muon_scheme_schaller;  // erosion rate using schaller scheme
  double erate_muon_scheme_braucher;  // erosion rate using braucher scheme
  
  double dEdMuonScheme;           // change in erosion rate for change in Muon Scheme
  double Muon_uncert;             // uncertainty of effective erosion rate 
                                  // in g/cm^2/yr for different muon schemes
  
  // variable for the production uncertainty
  double erate_prod_plus;   // erosion rate for positive production uncertainty
  double erate_prod_minus;  // erosion rate for negative production uncertainty
  
  double dEdProduction;        // change in erosion rate for change in production
  double Prod_uncert;          // uncertainty of effective erosion rate 
                               // in g/cm^2/yr for production uncertainty
  
  double this_prod_difference; // the difference in production for production uncertainty
  
  // initially we do not modify production rates
  bool is_production_uncertainty_plus_on = false;
  bool is_production_uncertainty_minus_on = false;
  
  // first get the prediction of the erosion rate
  cout << "LINE 2993 I'm predicting the erosion rate of this sample" << endl;
  erate = predict_CRN_erosion_point(Nuclide_conc, Nuclide, prod_uncert_factor, 
                              Muon_scaling, production_uncertainty,
                              average_production_rate,
                              is_production_uncertainty_plus_on,
                              is_production_uncertainty_minus_on, snow_eff_depth,
                              self_eff_depth, topo_shield,production);
  cout << "The erosion rate is: " << erate << endl;
  
  double no_prod_uncert = 1.0;    // set the scheme to no production uncertainty
                                  // for the external uncertainty
  // now get the external uncertainty                                        
  erate_external_plus = predict_CRN_erosion_point(Nuclide_conc+Nuclide_conc_err, Nuclide, 
                                             no_prod_uncert, Muon_scaling,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on, snow_eff_depth,
                                             self_eff_depth, topo_shield,production);
                                             
  erate_external_minus = predict_CRN_erosion_point(Nuclide_conc-Nuclide_conc_err, Nuclide, 
                                             no_prod_uncert, Muon_scaling,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on, snow_eff_depth,
                                             self_eff_depth, topo_shield,production);
                                             
  dEdExternal = (erate_external_plus-erate_external_minus)/(2*Nuclide_conc_err);
  External_uncert = fabs(dEdExternal*Nuclide_conc_err);
  
  //cout << "LSDCosmoBasin, line 1160, erate: " << erate << " and uncertainty: " 
  //     << External_uncert << endl;

  // now calculate uncertainty from different muon scaling schemes. 
  // The end members are Braucher and Schaller
  string braucher_string = "Braucher";
  string schaller_string = "Schaller";
  
  // get the difference in the pair
  LSDCRNParameters LSDCRNP;
  int pair_key = 0;       // this is for braucher-schaller
  vector<double> muon_uncert_diff = LSDCRNP.get_uncertainty_scaling_pair(pair_key);
  
  double this_muon_uncert_dif;
  if(Nuclide == "Be10")
  {
    this_muon_uncert_dif = muon_uncert_diff[0];
  }
  else if (Nuclide == "Al26")
  {
    this_muon_uncert_dif = muon_uncert_diff[1];
  }
  else
  {
    cout << "LINE 1295 LSDBasin you did not supply a valid nuclide, defaulting to 10Be" << endl;
    Nuclide = "Be10";
    this_muon_uncert_dif = muon_uncert_diff[0];
  }
  
  // now get the muon uncertainty
  erate_muon_scheme_schaller = predict_CRN_erosion_point(Nuclide_conc, Nuclide, 
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on, snow_eff_depth,
                                             self_eff_depth, topo_shield,production);
                                             
  erate_muon_scheme_braucher = predict_CRN_erosion_point(Nuclide_conc, Nuclide, 
                                             no_prod_uncert, braucher_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on, snow_eff_depth,
                                             self_eff_depth, topo_shield,production);
                                             
  dEdMuonScheme = (erate_muon_scheme_schaller-erate_muon_scheme_braucher)/
                  this_muon_uncert_dif;
  Muon_uncert = fabs(dEdMuonScheme*this_muon_uncert_dif);
  
  //cout << "LSDCosmoBasin, Line 1292, change in scaling production rate: " 
  //     << this_muon_uncert_dif << " erate Schal: "
  //     << erate_muon_scheme_schaller << " erate Braucher: " 
  //     << erate_muon_scheme_braucher << " and erate uncert: " << Muon_uncert << endl;
  
  // now get the production uncertainty
  // first set the scaling
  // reset scaling parameters. This is necessary since the F values are
  // reset for local scaling
  if (Muon_scaling == "Schaller" )
  {
    LSDCRNP.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    LSDCRNP.set_Braucher_parameters();
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
         << "Options are Schaller, Braucher, newCRONUS, and Granger." << endl
         << "You chose: " << Muon_scaling << endl
         << "Defaulting to Braucher et al (2009) scaling" << endl;
    LSDCRNP.set_Braucher_parameters();     
  }
  
  // now get the uncertainty parameters
  vector<double> prod;
  double prod_plus,prod_minus;
  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[0];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[0];
  }
  else if (Nuclide == "Al26")
  {
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[1];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[1];
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[0];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[0];
  }
  //cout << "Prod plus: " << prod_plus << " prod minus: " << prod_minus << endl;
  this_prod_difference = prod_plus+prod_minus;     
  
  
  is_production_uncertainty_plus_on = true;
  is_production_uncertainty_minus_on = false;
  erate_prod_plus = predict_CRN_erosion_point(Nuclide_conc, Nuclide, 
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on, snow_eff_depth,
                                             self_eff_depth, topo_shield,production);

  is_production_uncertainty_plus_on = false;
  is_production_uncertainty_minus_on = true;
  erate_prod_minus = predict_CRN_erosion_point(Nuclide_conc, Nuclide, 
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on, snow_eff_depth,
                                             self_eff_depth, topo_shield,production);
  
  dEdProduction = (erate_prod_plus-erate_prod_minus)/
                   this_prod_difference;
  Prod_uncert = fabs(dEdProduction*this_prod_difference);

  //cout << "LSDCosmoBasin, Line 1368, change in production rate for production uncertainty: " 
  //     << this_prod_difference << " erate plus: "
  //     << erate_prod_plus << " erate minus: " 
  //     << erate_prod_minus << " and erate uncert: " << Prod_uncert << endl;



  // now calculate the total uncertainty
  double total_uncert = sqrt( External_uncert*External_uncert +
                              Muon_uncert*Muon_uncert +
                              Prod_uncert*Prod_uncert);

  erate_uncert_vec.push_back(erate);
  erate_uncert_vec.push_back(External_uncert);
  erate_uncert_vec.push_back(Muon_uncert);
  erate_uncert_vec.push_back(Prod_uncert);
  erate_uncert_vec.push_back(total_uncert);
  
  return erate_uncert_vec;
} 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function calculates the erosion from a point location using Newton-Raphson
//  iteration
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoData::predict_CRN_erosion_point(double Nuclide_conc, string Nuclide, 
                                          double prod_uncert_factor,
                                          string Muon_scaling,
                                          double& production_uncertainty,
                                          double& average_production,
                                          bool is_production_uncertainty_plus_on,
                                          bool is_production_uncertainty_minus_on,
                                          double snow_eff_depth,
                                          double self_eff_depth,
                                          double topo_shield,
                                          double production)
{
  // effective erosion rates (in g/cm^2/yr) for running the Newton Raphson
  // iterations
  double erate_guess;
  double eff_erate_guess;
  //double this_eff_erosion_rate;
  //double d_eff_erosion_rate;
  
  double rho = 2650;  // this is the rock density but in fact it doesn't 
                      // really play a role since it is factored into the
                      // apparent erosion to get erosion in mm/yr but the divided
                      // out again. 
                      // The value 2650 is used because this is the default in cosmocalc

  // production uncertainty factor is a multiplier that sets the production 
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

  // initiate a particle. We'll just repeatedly call this particle
  // for the sample. 
  int startType = 0; 
  double Xloc = 0;
  double Yloc = 0;
  double  startdLoc = 0.0;
  double  start_effdloc = 0.0;
  double startzLoc = 0.0;
  // create a particle at zero depth
  LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                               startdLoc, start_effdloc, startzLoc);

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }

  
  // now set the scaling parameters
  if (Muon_scaling == "Schaller" )
  {
    LSDCRNP.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    LSDCRNP.set_Braucher_parameters();
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
         << "Options are Schaller, Braucher, newCRONUS and Granger." << endl
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
  else
  {
    cout << "LSDBasin line 1583, You didn't choos a valid nuclide. Defaulting"
         << " to 10Be." << endl;
    Nuclide = "Be10";
    nuclide_scaling_switches[0] = true; 
  }

  // now get the guess from the particle
  // the initial guess just takes scaling from the outlet, and then 
  // uses that for the entire basin. This guess will probably be quite
  // far off, but provides a useful starting point
  
  // the elevation, snow shielding, topographic shielding
  // and production scaling are all independent of the erosion rate
  // and are calculated seperately. 
  double total_shielding =  production*topo_shield;
                        
  //cout << "LSDBasin line 1128 Prod scaling is: " << production_scaling[0] << endl;
                        
  //cout << "LSDBasin line 1129; total scaling is: " << total_shielding << endl;
                      
  // now recalculate F values to match the total shielding
  LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);
  double initial_snow_shield = 1.0;
  LSDCRNP.set_neutron_scaling(production,topo_shield,
                             initial_snow_shield);
  
  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    //cout << "LINE 1134 LSDBasin Nuclide conc is: " << Nuclide_conc << endl;
    eroded_particle.setConc_10Be(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_10Be_neutron_only(rho, LSDCRNP);
    //cout << "Be10, initial erate guess in m/yr with density " << rho << ": " << erate_guess << endl;
  }
  else if (Nuclide == "Al26")
  {
    eroded_particle.setConc_26Al(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_26Al_neutron_only(rho, LSDCRNP);
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    eroded_particle.setConc_10Be(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_10Be_neutron_only(rho, LSDCRNP);
  }
  
  // convert to  g/cm^2/yr
  eff_erate_guess = 0.1*erate_guess*rho;
  
  // now using this as the initial guess, use Newton-Raphson to zero in on the
  // correct erosion rate
  double eff_e_new = eff_erate_guess; // the erosion rate upon which we iterate
  double eff_e_change;                // the change in erosion rate between iterations
  double tolerance = 1e-10;           // tolerance for a change in the erosion rate
                                      // between Newton-Raphson iterations
  double eff_e_displace = 1e-6;       // A small displacment in the erosion rate used
                                      // to calculate the derivative
  double N_this_step;                 // the concentration of the nuclide reported this step
  double N_displace;                  // the concentration at the displaced erosion rate
  double N_derivative;                // dN/de derivative for Newton-Raphson
  double f_x;                         // the function being tested by newton raphson
  double f_x_displace;                // the displaced function (for calculating the derivative)
  
  double this_step_prod_uncert;       // the uncertainty in the production rate
                                      // from this step
  double displace_uncertainty;        // the uncertainty from the displaced calculations
                                      // is not used so a dummy variable is used here
  do
  {
    // get the new values
    //cout << "LSDBasin line 1649 You are doing this wih the effective depth driven shielding" << endl;
    N_this_step = predict_mean_CRN_conc_point(eff_e_new, Nuclide,
                                        prod_uncert_factor,
                                        Muon_scaling,
                                        this_step_prod_uncert,
                                        is_production_uncertainty_plus_on,
                                        is_production_uncertainty_minus_on,
                                        snow_eff_depth, self_eff_depth,
                                        topo_shield, production);
      //cout << " Conc: " << N_this_step << endl;

    // now get the derivative
    N_displace = predict_mean_CRN_conc_point(eff_e_new+eff_e_displace,Nuclide,
                                       prod_uncert_factor,Muon_scaling, 
                                       displace_uncertainty,
                                       is_production_uncertainty_plus_on,
                                       is_production_uncertainty_minus_on,
                                        snow_eff_depth, self_eff_depth,
                                        topo_shield, production);

    f_x =  N_this_step-Nuclide_conc;
    f_x_displace =  N_displace-Nuclide_conc;
    
    N_derivative = (f_x_displace-f_x)/eff_e_displace;
      
    if(N_derivative != 0)
    {
      eff_e_new = eff_e_new-f_x/N_derivative;
      
      // check to see if the difference in erosion rates meet a tolerance
      eff_e_change = f_x/N_derivative;
      //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
    }
    else
    {
      eff_e_change = 0;
    }
  
  } while(fabs(eff_e_change) > tolerance);

  // replace the production uncertainty
  production_uncertainty = this_step_prod_uncert;
  
  return eff_e_new;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function calculates the Nuclide concentration from a point sample,
//  given an erosion rate (in g/cm^2/yr)
// it can cope with snow and self shielding
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoData::predict_mean_CRN_conc_point(double eff_erosion_rate, string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            double& production_uncertainty, 
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on,
                                            double snow_eff_depth,
                                            double self_eff_depth,
                                            double topo_shield,
                                            double production)
{
  double N;       // The number of atoms
  
  // production uncertainty factor is a multiplier that sets the production 
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

  // the total shielding. A product of snow, topographic and production scaling
  double total_shielding;
  
  // parameters for the shielding
  double this_top_eff_depth;
  double this_bottom_eff_depth;
  
  // initiate a particle. We'll just repeatedly call this particle
  // for the sample. 
  int startType = 0; 
  double Xloc = 0;
  double Yloc = 0;
  double  startdLoc = 0.0;
  double  start_effdloc = 0.0;
  double startzLoc = 0.0;
  
  // create a particle at zero depth
  LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                               startdLoc, start_effdloc, startzLoc);

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

  // at this stage we will try to replicate the basin averaging that goes on in 
  // most paper
  // set scaling parameters. This is necessary since the F values are
  // reset for local scaling
  if (Muon_scaling == "Schaller" )
  {
    LSDCRNP.set_Schaller_parameters();
  }
    else if (Muon_scaling == "Braucher" )
  {
    LSDCRNP.set_Braucher_parameters();
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
         << "Options are Schaller, Braucher, newCRONUS, and Granger." << endl
         << "You chose: " << Muon_scaling << endl
         << "Defaulting to Braucher et al (2009) scaling" << endl;
    LSDCRNP.set_Braucher_parameters();     
  }

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }

  // set the scaling to the correct production uncertainty
  vector<double> test_uncert;
  if(is_production_uncertainty_plus_on)
  {
    test_uncert = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
  }
  else if(is_production_uncertainty_minus_on)
  {
    test_uncert = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
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
  else
  {
    cout << "LSDBasin line 1583, You didn't choos a valid nuclide. Defaulting"
         << " to 10Be." << endl;
    Nuclide = "Be10";
    nuclide_scaling_switches[0] = true; 
  }

  // now you need logic to test if you are accounting for self shielding
  double total_shielding_no_uncert = production*topo_shield;
  total_shielding = prod_uncert_factor*total_shielding_no_uncert;

  // scale the F values
  LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);


  // get the snow shielding
  this_top_eff_depth = snow_eff_depth;

  // get the self_shielding
  this_bottom_eff_depth = this_top_eff_depth+self_eff_depth;
      
      
  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    //cout << "LInE 2271, 10Be" << endl;
    eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                       this_top_eff_depth, this_bottom_eff_depth);
    N=eroded_particle.getConc_10Be();
  }
  else if (Nuclide == "Al26")
  {
    //cout << "LINE 2278, 26Al" << endl;
    eroded_particle.update_26Al_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                       this_top_eff_depth, this_bottom_eff_depth);
    N=eroded_particle.getConc_26Al();
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                       this_top_eff_depth, this_bottom_eff_depth);
    N=eroded_particle.getConc_10Be();         
  }

  // return the nuclide concentration
  return N;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This takes a raster list name and a CRN params name and spits out a number of
// derivative rasters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::full_shielding_raster_printer(vector<string> Raster_names,
                                        vector<double> CRN_params)
{
  // find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // Load the DEM
  string DEM_bil_extension = "bil";
  //string fill_ext = "_fill";
  string DEM_fname = Raster_names[0];
  //cout << "Loading raster: " << DEM_fname << endl;
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
  //cout << "Filled raster" << endl;
  
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  //cout << "Got flow info" << endl;

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  //cout << "Got contributing pixels" << endl;
  
  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);
  //cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  //cout << "Got junction network" << endl;
  
  // Now convert the data into this UTM zone
  convert_to_UTM(filled_raster);

  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);
  //cout << "Snapped points" << endl;
  
  
  // after this operation the three int vectors valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertainties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }

  // Initiate pointers to the rasters
  LSDRaster Topographic_shielding;
  LSDRaster Snow_shielding;
  LSDRaster Self_shielding;

  // now, IF there are valid points, go on to the rest of the analysis
  if (valid_nuclide_names.size() != 0)
  {
    // if you need to write the inital raster, do so 
    LSDIndexRaster BasinIndex;   // initiate and empty raster

    // first check if topographic shielding raster exists
    cout << "Looking for toposhield raster, name is: " << Raster_names[3] << endl;
    if( Raster_names[3] != "NULL")
    {
      LSDRaster T_shield(Raster_names[3], DEM_bil_extension);
      Topographic_shielding = T_shield;
    }
    else
    {
      // get the topographic shielding
      cout << "Starting topographic shielding" << endl;
      LSDRaster T_shield = filled_raster.TopographicShielding(theta_step, phi_step);
      Topographic_shielding = T_shield;
      
      if(write_TopoShield_raster)
      {
        string TShield_name = Raster_names[0]+"_SH";
        T_shield.write_raster(TShield_name,DEM_bil_extension);
      }
      
    }
    
    // Now check if snow and self shielding rasters exist
    bool have_snow_raster;
    bool have_self_raster;
    double constant_snow_depth = 0;
    double constant_self_depth = 0;
    if (Raster_names[1] != "NULL")
    {
      cout << "LSDCosmoData, line 971: Loading the snow shielding raster, " 
           << Raster_names[1] << ".bil" <<  endl;
      LSDRaster Snow_shield(Raster_names[1], DEM_bil_extension);
      Snow_shielding = Snow_shield;
      have_snow_raster = true;
    }
    else
    {
      // snow shielding with a single parameter
      have_snow_raster = false;
      constant_snow_depth = CRN_params[0];
    }
    if (Raster_names[2] != "NULL")
     {
      cout << "LSDCosmoData, line 977: Loading the self shielding raster, " 
           << Raster_names[2] << ".bil" <<  endl;
      LSDRaster Self_shield(Raster_names[2], DEM_bil_extension);
      Self_shielding = Self_shield;
      have_self_raster = true;
    }
    else
    {
      // self shielding with a single parameter
      have_self_raster = false;
      constant_self_depth = CRN_params[1];
    }

    // some temporary doubles to hold the nuclide concentrations
    double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
    double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin

    //========================
    // LOOPING THROUGH BASINS
    //========================
    // now loop through the valid points, getting the cosmo data 
    cout << "-----------------------------------------------------------" << endl;
    cout << "I found " << n_valid_points << " valid CRN basins in this raster! " << endl;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
      if( valid_nuclide_names[samp] == "Be10")
      {
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
      }
      else if( valid_nuclide_names[samp] == "Al26")
      {
        test_N10 = 1e9;
        test_dN10 = 0;
        test_N26 = valid_concentrations[samp];
        test_dN26 = valid_concentration_uncertainties[samp];
      }
      else
      {
        cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
        cout << "Defaulting to Be10" << endl;
        valid_nuclide_names[samp] = "Be10";
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
      }
      

      cout << endl << "Valid point is: " << valid_cosmo_points[samp]
           << " Sample name: " << sample_name[ valid_cosmo_points[samp] ] << " Easting: " 
           << UTM_easting[valid_cosmo_points[samp]] << " Northing: "
           << UTM_northing[valid_cosmo_points[samp]] << endl;
      cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
           << snapped_junction_indices[samp] << endl;
      cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -" << endl;
      
      cout << "Making the basin" << endl;
      LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                              test_N10,test_dN10, test_N26,test_dN26);

      // we need to scale the sheilding parameters
      // now do the snow and self sheilding
      if (have_snow_raster)
      {
        if(have_self_raster)
        {
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, Self_shielding);
        }
        else
        {
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, constant_self_depth);        
        }
      }
      else
      {
        if(have_self_raster)
        {
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                constant_snow_depth, Self_shielding);      
        }
        else
        {
          thisBasin.populate_snow_and_self_eff_depth_vectors(constant_snow_depth, 
                                             constant_self_depth);
        }
      }

      // Now topographic shielding and production scaling
      cout << "Populating shielding vectors" << endl;
      thisBasin.populate_scaling_vectors(FlowInfo, filled_raster, 
                                         Topographic_shielding,
                                         path_to_atmospheric_data);

      // get the prefix of the basin vectors
      string basin_ID =  itoa(valid_cosmo_points[samp]);
      string basin_ext = "_BASIN"+basin_ID;
      string thisfilename = Raster_names[0]+basin_ext;

      // get the erosion rate from this sample
      vector<double> erate_analysis = erosion_rate_results[ valid_cosmo_points[samp] ];
      cout << "The erosion rate for this basin is: " << erate_analysis[0] << " g/m^2/yr" << endl;
      

      // now print to file
      cout << "Printing the concentration raster" << endl;
      thisBasin.print_CRN_conc_raster(thisfilename, erate_analysis[0], 
                                      nuclide[ valid_cosmo_points[samp] ],
                                      Muon_scaling, FlowInfo);
      cout << "Printing the scaling rasters"  << endl;
      thisBasin.print_scaling_and_shielding_rasters(thisfilename,FlowInfo);

      cout << "Finished printing rasters from this basin." << endl;

    }   // finished looping through basins
    

     
  }    // finished logic for a DEM with valid points
  else
  {
    cout << "There are no valid CRN points in this raster" << endl;
  }

}                                        


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function loops though the file structures calculating 
// cosmogenic-derived erosion rates and uncertainties
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::calculate_erosion_rates(int method_flag)
{

  if (method_flag != 0 && method_flag != 1 && method_flag != 2 && method_flag != 3)
  {
    cout << "Warning, you did not give a valid method flag. Assuming you've used spawned basins" << endl;
    method_flag = 2;
  }

  // find out how many DEMs there are:
  int n_DEMS = int(DEM_names_vecvec.size());

  vector<string> this_Raster_names;
  vector<double> this_Param_names;
  
  // now loop through the DEMs
  for (int iDEM = 0; iDEM< n_DEMS; iDEM++)
  {
    this_Raster_names = DEM_names_vecvec[iDEM];
    this_Param_names = snow_self_topo_shielding_params[iDEM];
    
    if (method_flag == 0)
    {
      basic_cosmogenic_analysis(this_Raster_names[0]);
    }
    else if (method_flag == 1)
    {
      full_shielding_cosmogenic_analysis(this_Raster_names,this_Param_names);
    }
    else if (method_flag == 2)
    {
      full_shielding_cosmogenic_analysis_for_spawned(this_Raster_names,this_Param_names);
    }
    else if (method_flag == 3)
    {
      cout << "You are doing soil samples on DEM: " <<  this_Raster_names[0] << endl;
      Soil_sample_calculator(this_Raster_names,this_Param_names);
    }
    else
    {
      full_shielding_cosmogenic_analysis(this_Raster_names,this_Param_names);
    }
    
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function loops though the file structures calculating 
// cosmogenic-derived denudation rates and uncertainties
// This version uses nesting: it points to rasters with known erosion rates
// in order to calculate the nested erosion rate. 
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::calculate_nested_erosion_rates()
{

  // find out how many DEMs there are:
  int n_DEMS = int(DEM_names_vecvec.size());

  vector<string> this_Raster_names;
  vector<double> this_Param_names;
  
  // now loop through the DEMs
  for (int iDEM = 0; iDEM< n_DEMS; iDEM++)
  {
    this_Raster_names = DEM_names_vecvec[iDEM];
    this_Param_names = snow_self_topo_shielding_params[iDEM];
    
    // check if the raster has a known erosion rate raster
    string DEM_name = this_Raster_names[0];
    string known_erate_name = DEM_name+"_ERKnown";
    string known_erate_header = known_erate_name+".hdr";
    
    // see if the known erate file exists
    // make sure the filename works
    ifstream ifs(known_erate_header.c_str());
    if( ifs.fail() )
    {
      cout << "\nThere is no known erosion rate raster for this DEM." << endl;
    }
    else
    {
      // check to make sure the dimensions of this raster match
      string bil_ext = "bil";
      LSDRasterInfo RI_ER(known_erate_name,bil_ext);
      LSDRasterInfo RI_DEM(this_Raster_names[0],bil_ext);
      if(RI_ER == RI_DEM)
      {
        LSDRaster known_rate_raster(known_erate_name,bil_ext);
        full_shielding_cosmogenic_analysis_nested(this_Raster_names,this_Param_names, 
                            known_rate_raster);
      }
      else
      {
        cout << "Your known erosion rate raster does not have the same dimensions" << endl;
        cout << "as your DEM." << endl;
        cout << "The rasters are:" << endl;
        cout << this_Raster_names[0] << endl;
        cout << known_erate_name << endl;
      }
    }
    
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::point_measurements(vector<int> valid_samples,vector<double> snow_thickness, 
                                      vector<double> self_thickness,
                                      vector<double> toposhield,
                                      vector<double> production_scaling, 
                                      string muon_string)
{
  int N_valid_samples =  int(valid_samples.size());

  // now loop through the samples, printing the ones with data
  for (int i = 0; i<N_valid_samples; i++)
  {
    // initiate the vector for holding the erate
    vector<double> erate_info;
    vector<double> erate_info_AMS_err;
    vector<double> erate_info_scaling_err;
    vector<double> erate_info_prod_err;  
    
    // get this sample
    int this_sample = valid_samples[i];
    double this_C = Concentration[this_sample];
    double this_dC = Concentration_uncertainty[this_sample];
    string this_nuclide = nuclide[this_sample];
    double top_eff_depth = snow_thickness[i];     // even if there is shielding these
    double bottom_eff_depth = snow_thickness[i]+self_thickness[i];  // get subsumed into the combined scaling 
    double combined_scaling = toposhield[i]* production_scaling[i];
    
    double rho = 2650;
    
    // now get apparent erosion rates that cosmocalc would have produced
    int startType = 0; 
    double Xloc = 0;
    double Yloc = 0;
    double  startdLoc = 0.0;
    double  start_effdloc = 0.0;
    double startzLoc = 0.0;
      
    LSDCRNParticle test_particle(startType, Xloc, Yloc,
                              startdLoc, start_effdloc, startzLoc);
      
    LSDCRNParameters LSDCRNP;
    // at this stage we will try to replicate the basin averaging that goes on in 
    // most paper
    // set scaling parameters. This is necessary since the F values are
    // reset for local scaling
    if (Muon_scaling == "Schaller" )
    {
      LSDCRNP.set_Schaller_parameters();
    }
      else if (Muon_scaling == "Braucher" )
    {
      LSDCRNP.set_Braucher_parameters();
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
           << "Options are Schaller, Braucher, newCRONUS, and Granger." << endl
           << "You chose: " << Muon_scaling << endl
           << "Defaulting to Braucher et al (2009) scaling" << endl;
      LSDCRNP.set_Braucher_parameters();     
    }
    
    
    if (nuclide[i] == "Be10")
    {
      test_particle.setConc_10Be(this_C);
      erate_info=test_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling,
                                 muon_string, top_eff_depth, bottom_eff_depth);
    }
    else if(nuclide[i] == "Al26")    
    {
      test_particle.setConc_26Al(this_C);
      erate_info=test_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling,
                                 muon_string, top_eff_depth, bottom_eff_depth);
    }
    
    // now do the error from the AMS measurement
    if (nuclide[i] == "Be10")
    {
      test_particle.setConc_10Be(this_C+this_dC);
      erate_info_AMS_err=test_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling,
                                 muon_string, top_eff_depth, bottom_eff_depth);
    }
    else if(nuclide[i] == "Al26")    
    {
      test_particle.setConc_26Al(this_C+this_dC);
      erate_info_AMS_err=test_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling,
                                 muon_string, top_eff_depth, bottom_eff_depth);
    }
    
    // now do the error from the scaling
    string new_muon_string = "Schaller";
    if (nuclide[i] == "Be10")
    {
      test_particle.setConc_10Be(this_C);
      erate_info_scaling_err=test_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling,
                                 new_muon_string, top_eff_depth, bottom_eff_depth);
    }
    else if(nuclide[i] == "Al26")    
    {
      test_particle.setConc_26Al(this_C);
      erate_info_scaling_err=test_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling,
                                 new_muon_string, top_eff_depth, bottom_eff_depth);
    }
    
    // now do the error from the production, assumes 10% production error
    if (nuclide[i] == "Be10")
    {
      test_particle.setConc_10Be(this_C);
      erate_info_scaling_err=test_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling*(1.1),
                                 muon_string, top_eff_depth, bottom_eff_depth);
    }
    else if(nuclide[i] == "Al26")    
    {
      test_particle.setConc_26Al(this_C);
      erate_info_scaling_err=test_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, 
                                 combined_scaling*(1.1),
                                 muon_string, top_eff_depth, bottom_eff_depth);
    }
  }
      
}









//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the determination of cosmogenic erosion rates 
// It allows nesting by accounting for erosion rates known from other 
//  samples (for example and upstream cosmogenic data point)
// THIS IS FOR A SINGLE RASTER
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::full_shielding_CRN_concentration_predictor(vector<string> Raster_names,
                            vector<double> CRN_params, 
                            LSDRaster& known_eff_erosion)
{

  cout << endl << endl << "================================================" << endl;
  cout << "Looking for basins in raster: " << Raster_names[0] << endl << endl;

  // some parameters for printing the basins, if that is called for
  int basin_number;
  int basin_pixel_area;
  map<int,int> basin_area_map;

  // now find valid points
  vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
  vector<int> snapped_node_indices;       // a vector to hold the valid node indices
  vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
  // Load the DEM
  string DEM_bil_extension = "bil";
  //string fill_ext = "_fill";
  string DEM_fname = Raster_names[0];
  //cout << "Loading raster: " << DEM_fname << endl;
  LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  topo_test.remove_seas();
  
  // check to see if the rasters are the same
  LSDRasterInfo DEM_info(topo_test);
  LSDRasterInfo Erosion_info(known_eff_erosion);
  if (DEM_info != Erosion_info)
  {
    cout << "LSDCosmoData::full_shielding_cosmogenic_analysis_nested ERROR!" << endl;
    cout << "The erosion raster and DEM are not the same dimesions." <<endl;
    exit(EXIT_SUCCESS);
  }
  
  // Fill this raster
  LSDRaster filled_raster = topo_test.fill(min_slope);
  //cout << "Filled raster" << endl;
  
  // get the flow info
  LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);
  //cout << "Got flow info" << endl;

  // get contributing pixels (needed for junction network)
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  //cout << "Got contributing pixels" << endl;
  
  // get the sources
  vector<int> sources;
  sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);
  //cout << "Got sources" << endl;

  // now get the junction network
  LSDJunctionNetwork JNetwork(sources, FlowInfo);
  //cout << "Got junction network" << endl;
  
  // Now convert the data into this UTM zone
  convert_to_UTM(filled_raster);

  // convert UTM vectors to float
  vector<float> fUTM_easting;
  vector<float> fUTM_northing;
  for (int i = 0; i< int(UTM_easting.size()); i++)
  {
    fUTM_easting.push_back( float(UTM_easting[i]));
    fUTM_northing.push_back( float(UTM_northing[i]));
  }
  
  cout << "Getting snapped basins" << endl;
  JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
            search_radius_nodes, threshold_stream_order, FlowInfo, 
            valid_cosmo_points, snapped_node_indices, snapped_junction_indices);
  cout << "Snapped points" << endl;
  
  
  // after this operation the three int vectors: valid_cosmo_points, snapped_node_indices,
  // and snapped_junction_indices should be populated with valid points
  // you now need to get the concentrations and uncertanties from these
  // points
  
  // first, get the data elements from the object
  vector<string> valid_nuclide_names;
  vector<double> valid_concentrations;
  vector<double> valid_concentration_uncertainties;
  
  int n_valid_points = int(valid_cosmo_points.size());
  for (int i = 0; i< int(n_valid_points); i++)
  {
    valid_nuclide_names.push_back(nuclide[ valid_cosmo_points[i] ] );
    valid_concentrations.push_back( Concentration[ valid_cosmo_points[i] ] );
    valid_concentration_uncertainties.push_back( 
                          Concentration_uncertainty[ valid_cosmo_points[i] ] );
  }
  cout << "Got valid points, there are " << n_valid_points << " of them." << endl;
  
  
  // Initiate pointers to the rasters
  LSDRaster Topographic_shielding;
  LSDRaster Snow_shielding;
  LSDRaster Self_shielding;
  
  
  // a flag for writing the initial basin index
  bool written_inital_basin_index = false;
  
  // now, IF there are valid points, go on to the rest of the analysis
  if (valid_nuclide_names.size() != 0)
  {
    // if you need to write the inital raster, do so 
    LSDIndexRaster BasinIndex;   // initiate and empty raster

    
    // first check if topographic shielding raster exists
    cout << "Looking for toposhield raster, name is: " << Raster_names[3] << endl;
    if( Raster_names[3] != "NULL")
    {
      LSDRaster T_shield(Raster_names[3], DEM_bil_extension);
      Topographic_shielding = T_shield;
    }
    else
    {
      // get the topographic shielding
      cout << "Starting topographic shielding" << endl;
      LSDRaster T_shield = filled_raster.TopographicShielding(theta_step, phi_step);
      Topographic_shielding = T_shield;
      
      if(write_TopoShield_raster)
      {
        string TShield_name = Raster_names[0]+"_SH";
        T_shield.write_raster(TShield_name,DEM_bil_extension);
      }
      
    }
    
    // Now check if snow and self shielding rasters exist
    bool have_snow_raster;
    bool have_self_raster;
    double constant_snow_depth = 0;
    double constant_self_depth = 0;
    if (Raster_names[1] != "NULL")
    {
      cout << "LSDCosmoData, line 2058: Loading the snow shielding raster, " 
           << Raster_names[1] << ".bil" <<  endl;
      LSDRaster Snow_shield(Raster_names[1], DEM_bil_extension);
      Snow_shielding = Snow_shield;
      have_snow_raster = true;
    }
    else
    {
      // snow shielding with a single parameter
      have_snow_raster = false;
      constant_snow_depth = CRN_params[0];
    }
    if (Raster_names[2] != "NULL")
     {
      cout << "LSDCosmoData, line 977: Loading the self shielding raster, " 
           << Raster_names[2] << ".bil" <<  endl;
      LSDRaster Self_shield(Raster_names[2], DEM_bil_extension);
      Self_shielding = Self_shield;
      have_self_raster = true;
    }
    else
    {
      // self shielding with a single parameter
      have_self_raster = false;
      constant_self_depth = CRN_params[1];
    }

    // some temporary doubles to hold the nuclide concentrations
    double test_N10, test_dN10;   // concetration and uncertainty of 10Be in basin
    double test_N26, test_dN26;   // concetration and uncertainty of 26Al in basin
    double test_N, test_dN;       // concentration and uncertainty of the nuclide in basin  

    //========================
    // LOOPING THROUGH BASINS
    //========================
    // now loop through the valid points, getting the cosmo data 
    cout << "-----------------------------------------------------------" << endl;
    cout << "I found " << n_valid_points << " valid CRN basins in this raster! " << endl;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
      if( valid_nuclide_names[samp] == "Be10")
      {
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      else if( valid_nuclide_names[samp] == "Al26")
      {
        test_N10 = 1e9;
        test_dN10 = 0;
        test_N26 = valid_concentrations[samp];
        test_dN26 = valid_concentration_uncertainties[samp];
        test_N = test_N26;
        test_dN = test_dN26;
      }
      else
      {
        cout << "You did not select a valid nuclide name, options are Be10 and Al26" << endl;
        cout << "Defaulting to Be10" << endl;
        valid_nuclide_names[samp] = "Be10";
        test_N10 = valid_concentrations[samp];
        test_dN10 = valid_concentration_uncertainties[samp];
        test_N26 = 1e9;
        test_dN26 = 0;
        test_N = test_N10;
        test_dN = test_dN10;
      }
      

      cout << endl << "Valid point is: " << valid_cosmo_points[samp]
           << " Sample name: " << sample_name[ valid_cosmo_points[samp] ] << " Easting: " 
           << UTM_easting[valid_cosmo_points[samp]] << " Northing: "
           << UTM_northing[valid_cosmo_points[samp]] << endl;
      cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
           << snapped_junction_indices[samp] << endl;
      cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -" << endl;
      LSDCosmoBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork,
                              test_N10,test_dN10, test_N26,test_dN26);

      // write the index basin if flag is set to true
      if(write_basin_index_raster)
      {
        cout << "I'm writing a basin index number for you" << endl;
      
        if (not written_inital_basin_index)
        {
        
          basin_number = valid_cosmo_points[samp];
          basin_pixel_area = thisBasin.get_NumberOfCells();
          basin_area_map[basin_number] = basin_pixel_area;
          LSDIndexRaster NewBasinIndex = 
             thisBasin.write_integer_data_to_LSDIndexRaster(basin_number, FlowInfo);
          BasinIndex = NewBasinIndex;
          written_inital_basin_index = true;
        }
        else
        {
          basin_number = valid_cosmo_points[samp];
          thisBasin.add_basin_to_LSDIndexRaster(BasinIndex, FlowInfo,
                                                basin_area_map,basin_number);
        }
      }

      // we need to scale the shielding parameters
      // now do the snow and self shielding
      if (have_snow_raster)
      {
        cout << "I've got the snow raster" << endl;
      
        if(have_self_raster)
        {
          cout << "I've also got the self raster." << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, Self_shielding);
          cout << "Done with effective depths" << endl;                           
        }
        else
        {
          cout << "No self raster, but I'm getting the effective depths" << endl;
          cout << "The constant self depth is: " <<  constant_self_depth << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                Snow_shielding, constant_self_depth);
          cout << "Done with effective depths" << endl;        
        }
      }
      else
      {
        if(have_self_raster)
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(FlowInfo, 
                                constant_snow_depth, Self_shielding);
          cout << "Done with effective depths" << endl;                                
        }
        else
        {
          cout << "Getting effective depths" << endl;
          thisBasin.populate_snow_and_self_eff_depth_vectors(constant_snow_depth, 
                                             constant_self_depth);
          cout << "Done with effective depths" << endl;                                 
        }
      }


      cout << "Now I will populate the scaling vectors." << endl;
      // Now topographic shielding and production scaling
      thisBasin.populate_scaling_vectors(FlowInfo, filled_raster, 
                                         Topographic_shielding,
                                         path_to_atmospheric_data);
      cout << "The scaling vectors are populated. I am moving on to the analysis" << endl;

      // now do the analysis
      cout << "Line 2571, doing analysis" << endl;
      vector<double> erate_analysis = thisBasin.full_CRN_erosion_analysis_nested(known_eff_erosion, FlowInfo, test_N, 
                                          valid_nuclide_names[samp], test_dN, 
                                          prod_uncert_factor, Muon_scaling);
       cout << "erate: " << erate_analysis[0] << endl;
      
    
    
      // now get parameters for cosmogenic calculators
      vector<double> param_for_calc = 
          thisBasin.calculate_effective_pressures_for_calculators_nested(filled_raster,
                                            FlowInfo, path_to_atmospheric_data, 
                                            known_eff_erosion);

      cout << "Paramforcalc size: " << param_for_calc.size() << endl;              
      cout << "Getting pressures" << endl;


      // get the relief of the basin
      float R = thisBasin.CalculateBasinRange(FlowInfo, filled_raster);
      double relief = double(R);

      MapOfProdAndScaling["BasinRelief"][ valid_cosmo_points[samp] ] = relief;
      MapOfProdAndScaling["AverageProdScaling"][ valid_cosmo_points[samp] ] = param_for_calc[0];
      MapOfProdAndScaling["AverageTopoShielding"][ valid_cosmo_points[samp] ] = param_for_calc[1];
      MapOfProdAndScaling["AverageSelfShielding"][ valid_cosmo_points[samp] ] = param_for_calc[2];
      MapOfProdAndScaling["AverageSnowShielding"][ valid_cosmo_points[samp] ] = param_for_calc[3];
      MapOfProdAndScaling["AverageShielding"][ valid_cosmo_points[samp] ] =  param_for_calc[11];
      MapOfProdAndScaling["AverageCombinedScaling"][ valid_cosmo_points[samp] ] = param_for_calc[4];
      MapOfProdAndScaling["outlet_lat"][ valid_cosmo_points[samp] ] = param_for_calc[5];
      MapOfProdAndScaling["OutletPressure"][ valid_cosmo_points[samp] ] = param_for_calc[6];
      MapOfProdAndScaling["OutletEffectivePressure"][ valid_cosmo_points[samp] ] = param_for_calc[7];
      MapOfProdAndScaling["centroid_lat"][ valid_cosmo_points[samp] ] = param_for_calc[8];
      MapOfProdAndScaling["CentroidPressure"][ valid_cosmo_points[samp] ] = param_for_calc[9];
      MapOfProdAndScaling["CentroidEffectivePressure"][ valid_cosmo_points[samp] ] = param_for_calc[10];
    
      // add the erosion rate results to the holding data member
      erosion_rate_results[ valid_cosmo_points[samp] ] = erate_analysis;

      //cout << "finished adding data" << endl;

    }  // finished looping thorough basins
    
    // now print the basin LSDIndexRaster
    if(write_basin_index_raster)
    {
      string basin_ext = "_BASINS";
      string basin_fname =  DEM_fname+basin_ext;
      BasinIndex.write_raster(basin_fname, DEM_bil_extension);
    }
     
  }    // finsiehd logic for a DEM with valid points
  else
  {
    cout << "There are no valid CRN points in this raster" << endl;
  }
  cout << "==========================================" << endl << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-















//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function prints derivative rasters from the cosmogenic analysis
// It takes map of bools to determine which rasters to print. 
// The function recalculates a number of the basinwide parameters: 
// efficiency is sacrificed for better code readability
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_rasters()
{

  if(write_full_scaling_rasters)
  {
  
    cout << endl << endl <<"============================================" << endl;
    cout << "I am now going to print some rasters for you. You will find these" << endl;
    cout << "In the folders with the DEMs" << endl;
    cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -" << endl;
  
    // find out how many DEMs there are:
    int n_DEMS = int(DEM_names_vecvec.size());

    vector<string> this_Raster_names;
    vector<double> this_Param_names;
  
    // now loop through the DEMs
    for (int iDEM = 0; iDEM< n_DEMS; iDEM++)
    {
      //cout << "Doing DEM " << iDEM << " of " <<  n_DEMS << endl;
      
      this_Raster_names = DEM_names_vecvec[iDEM];
      this_Param_names = snow_self_topo_shielding_params[iDEM];
      
      full_shielding_raster_printer(this_Raster_names,this_Param_names);
    }
  }
  else
  {
    cout << "LSDCosmoData::print_rasters, raster printing set to false" << endl;
    cout << "If you wish to print out production, shielding, scaling and" << endl;
    cout << "predicted concentration rasters, add this line to the parameter file: " << endl;
    cout << "write_full_scaling_rasters: true" << endl;
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints all valid results to a csv file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_results()
{
  string combined_filename = path+param_name;

  // set up extensions to the files
  string results_ext = "_CRNResults.csv";
  string CRONUS_ext = "_CRONUSInput.txt";

  // get the filenames to open
  string crn_results_fname = path+param_name+results_ext;
  string CRONUS_results_fname = path+param_name+CRONUS_ext;

  // open the file
  ofstream results_out;
  results_out.open(crn_results_fname.c_str());
  results_out.precision(8);
  
  ofstream CRONUS_out;
  CRONUS_out.open(CRONUS_results_fname.c_str());
  CRONUS_out.precision(8);
  
  // a message for cronus file
  CRONUS_out << "->IMPORTANT nuclide concentrations are not original!" << endl
             << "  They are scaled to the 07KNSTD!!" << endl;
  CRONUS_out << "->Scaling is averaged over the basin for snow, self and topographic shielding." << endl;
  CRONUS_out << "->Snow and self shielding are considered by neutron spallation only." << endl;
  CRONUS_out << "->Pressure is an effective pressure that reproduces Stone scaled production" << endl
             << "  that is calculated on a pixel by pixel basis." << endl;
  CRONUS_out << "->Self shielding is embedded in the shielding calculation and so" << endl
             << "  sample thickness is set to 0." << endl;
  CRONUS_out << "->These results generated by the LSDCosmoBasin program," << endl
             << " written by Simon M. Mudd Martin D. Hurst and Stuart W.D. Grieve" << endl;
  CRONUS_out << "=========================================================" << endl;           
    
  // first the header line
  results_out << "All concentrations are standardised to 07KNSTD using Balco et al 2008" << endl;
  results_out << "Results generated by the LSDCosmoBasin program," << endl
              << "written by Simon M. Mudd Martin D. Hurst and Stuart W.D. Grieve" << endl;
  results_out << "basin_ID,sample_name,nuclide,latitude,longitude,concentration,concentration_uncert,"
       << "erate_g_percm2_peryr,AMS_uncert,muon_uncert,production_uncert,total_uncert,"
       << "AvgProdScaling,AverageTopoShielding,AverageSelfShielding,"
       << "AverageSnowShielding,AverageShielding,AvgShield_times_AvgProd,AverageCombinedScaling,outlet_latitude,"
       << "OutletPressure,OutletEffPressure,centroid_latitude,CentroidPressure,"
       << "CentroidEffPressure,eff_erate_COSMOCALC,erate_COSMOCALC_mmperkyr_rho2650,"
       << "eff_erate_COSMOCALC_emulating_CRONUS,erate_COSMOCALC_emulating_CRONUS_mmperkyr_rho2650,"
       << "erate_mmperkyr_rho2650,erate_totalerror_mmperkyr_rho2650,basin_relief"
       << endl;
  
  double rho = 2650;
  
  // now loop through the sampes, printing the ones with data
  for (int i = 0; i<N_samples; i++)
  {
  
    //cout << "On sample number " << i << ", size of results vector: " << erosion_rate_results[i].size() << endl;
    // don't print the results unless they exist
    if (int(erosion_rate_results[i].size()) > 0)
    {
      // get the results from this sample
      vector<double> erate_analysis = erosion_rate_results[i];

      // now get apparent erosion rates that cosmocalc would have produced
      int startType = 0; 
      double Xloc = 0;
      double Yloc = 0;
      double  startdLoc = 0.0;
      double  start_effdloc = 0.0;
      double startzLoc = 0;
      
      LSDCRNParticle test_particle(startType, Xloc, Yloc,
                              startdLoc, start_effdloc, startzLoc);
      LSDCRNParameters LSDCRNP;
      string muon_string = "Braucher";
      double top_eff_depth = 0;     // even if there is shielding these
      double bottom_eff_depth = 0;  // get subsumed into the combined scaling 
      
      vector<double> erate_info;
      vector<double> erate_info_CCCR;   // this gets a cosmocalc erosion rate but
                                        // using a CRONUS-like scaling
      
      if (nuclide[i] == "Be10")
      {
        test_particle.setConc_10Be(Concentration[i]);
        erate_info=test_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, 
                                 MapOfProdAndScaling["AverageCombinedScaling"][i],
                                 muon_string, top_eff_depth, bottom_eff_depth);
      }
      else if(nuclide[i] == "Al26")    
      {
        test_particle.setConc_26Al(Concentration[i]);
        erate_info=test_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, 
                                 MapOfProdAndScaling["AverageCombinedScaling"][i],
                                 muon_string, top_eff_depth, bottom_eff_depth);
      }
      
      // do a second calculation based on the CRONUS equivalent erosion
      if (nuclide[i] == "Be10")
      {
        test_particle.setConc_10Be(Concentration[i]);
        erate_info_CCCR=test_particle.apparent_erosion_10Be_COSMOCALC(rho, LSDCRNP, 
                                 MapOfProdAndScaling["AverageShielding"][i]*
                                 MapOfProdAndScaling["AverageProdScaling"][i],  
                                 muon_string, top_eff_depth, bottom_eff_depth);
      }
      else if(nuclide[i] == "Al26")    
      {
        test_particle.setConc_26Al(Concentration[i]);
        erate_info_CCCR=test_particle.apparent_erosion_26Al_COSMOCALC(rho, LSDCRNP, 
                                 MapOfProdAndScaling["AverageShielding"][i]*
                                 MapOfProdAndScaling["AverageProdScaling"][i],
                                 muon_string, top_eff_depth, bottom_eff_depth);
      }

      // print to the results file
      results_out << i << "," << sample_name[i] << "," << nuclide[i] << "," << latitude[i] 
                  << "," << longitude[i] << "," << Concentration[i] 
                  << "," << Concentration_uncertainty[i] << "," 
                  << erate_analysis[0] << "," << erate_analysis[1] << "," 
                  << erate_analysis[2] << "," << erate_analysis[3] << "," 
                  << erate_analysis[4] << "," 
                  << MapOfProdAndScaling["AverageProdScaling"][i] << ","
                  << MapOfProdAndScaling["AverageTopoShielding"][i] << "," 
                  << MapOfProdAndScaling["AverageSelfShielding"][i] << "," 
                  << MapOfProdAndScaling["AverageSnowShielding"][i] << "," 
                  << MapOfProdAndScaling["AverageShielding"][i] << ","
                  << MapOfProdAndScaling["AverageShielding"][i]*
                     MapOfProdAndScaling["AverageProdScaling"][i] << "," 
                  << MapOfProdAndScaling["AverageCombinedScaling"][i] << ","
                  << MapOfProdAndScaling["outlet_lat"][i] << "," 
                  << MapOfProdAndScaling["OutletPressure"][i] << "," 
                  << MapOfProdAndScaling["OutletEffectivePressure"][i] << ","  
                  << MapOfProdAndScaling["centroid_lat"][i] << "," 
                  << MapOfProdAndScaling["CentroidPressure"][i] << "," 
                  << MapOfProdAndScaling["CentroidEffectivePressure"][i] << ","
                  << erate_info[0] << "," << erate_info[0]*1e7/rho << ","
                  << erate_info_CCCR[0] << "," << erate_info_CCCR[0]*1e7/rho << ","
                  << erate_analysis[0]*1e7/rho <<","<< erate_analysis[4]*1e7/rho 
                  << "," << MapOfProdAndScaling["BasinRelief"][i] << endl;


      
      // now print to the CRONUS file
      // Note this subsumes self shielding into the shielding factor so 
      // the self shielding thickness is set to 0!
      CRONUS_out << sample_name[i] << " " << latitude[i] << " " << longitude[i]
                 << " " << MapOfProdAndScaling["CentroidEffectivePressure"][i] 
                 << " pre 0 2.65 " << MapOfProdAndScaling["AverageShielding"][i] << " ";
      if(nuclide[i] == "Be10")
      {
        CRONUS_out << Concentration[i] << " " << Concentration_uncertainty[i] 
                   << " 07KNSTD 0 0 KNSTD" << endl;
      }
      else if(nuclide[i] == "Al26")
      {
        CRONUS_out << "0 0 07KNSTD " << Concentration[i] << " " 
                   << Concentration_uncertainty[i] << " KNSTD" << endl;      
      }           
        
    }
  }
  results_out.close();
  CRONUS_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints a production raster for any of the rasters in the list
// that has a valid cosmo point
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_scaling_and_shielding_complete_rasters()
{
  // vector that will contain the snow and self shielding values, if they are constant
  vector<double> this_Params;
  
  // loop through the lines in the files, checking to see if the 
  // georeferencing is equivalent
  string bil_ext = "bil";
  string null_str = "NULL";
  int N_DEMS = int(DEM_names_vecvec.size());
  for(int iDEM = 0; iDEM<N_DEMS; iDEM++)
  {
    // get snow and self shielding parameters (this will only be used if the now and
    // self shielding rasters are empty)
    this_Params = snow_self_topo_shielding_params[iDEM];
            
    // get the names from this DEM
    vector<string>  DEM_names_vec = DEM_names_vecvec[iDEM];
    cout << "Checking rasters, this raster is: " << DEM_names_vec[0] << endl;
    string DEM_fname = DEM_names_vec[0];
    
    // load the DEM
    LSDRaster topo_test(DEM_fname, bil_ext);
    
    // Fill the DEM, since the production rates for the calculators are derived 
    // the filled rasters
    LSDRaster filled_raster = topo_test.fill(min_slope);
    
    // get the rows and columns
    int NRows, NCols;
    NRows = filled_raster.get_NRows();
    NCols = filled_raster.get_NCols();
    float NoDataValue =  filled_raster.get_NoDataValue();
    float XMinimum =  filled_raster.get_XMinimum();
    float YMinimum =  filled_raster.get_YMinimum();
    float DataResolution = filled_raster.get_DataResolution();
    map<string,string> GeoReferencingStrings = filled_raster.get_GeoReferencingStrings();
    
    // initiate the arrays that will contain the data
    Array2D<float> NewPressure(NRows,NCols,NoDataValue);
    Array2D<float> NewScaling(NRows,NCols,NoDataValue);
    Array2D<float> NewCombinedShielding(NRows,NCols,NoDataValue);
    Array2D<float> NewCombinedScaling(NRows,NCols,NoDataValue);

    // see if there are raster for the snow, self and toposhield rasters
    // now toposhielding. Initially, we assume there are none. 
    bool there_is_toposhield = false;
    bool there_is_snowshield = false;
    bool there_is_selfshield = false;
    
    // First initiate pointers to the rasters
    LSDRaster Topographic_shielding;
    LSDRaster Snow_shielding;
    LSDRaster Self_shielding;
    


    if(DEM_names_vec[1] != null_str)
    {
      there_is_snowshield = true;
      LSDRaster Sn_shield(DEM_names_vec[1], bil_ext);
      Snow_shielding = Sn_shield;
    }
    // now self shielding
    if(DEM_names_vec[2] != null_str)
    {
      there_is_selfshield = true;
      LSDRaster Slf_shield(DEM_names_vec[2], bil_ext);
      Self_shielding = Slf_shield;
    }
    // now toposhielding
    if(DEM_names_vec[3] != null_str)
    {
      there_is_toposhield = true;
      LSDRaster T_shield(DEM_names_vec[3], bil_ext);
      Topographic_shielding = T_shield;
    }

    // now create the CRN parameters object
    LSDCRNParameters LSDCRNP;
    double gamma_spallation = 160;      // in g/cm^2: spallation attentuation depth

    // get the atmospheric parameters
    LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
    LSDCRNP.set_CRONUS_data_maps();
  
    // a function for scaling stone production, defaults to 1
    double Fsp = 1.0;
  
    // the latitude and longitude
    double lat,longitude;
  
    // decalre converter object
    LSDCoordinateConverterLLandUTM Converter;
    
    // Get some temporary variables for holding the data. These will get printed
    // to and array.
    float this_elevation;
    float this_pressure;
    float this_scaling;
    float this_shielding;
    float this_combined_scaling;
    
    // teporary values for the shielding
    float this_snowshield;
    float this_selfshield;
    float this_toposhield;

    // now loop through the rows and columns, populating the scaling values
    for(int row = 0; row<NRows; row++)
    {
      for(int col = 0; col<NCols; col++)
      {
        // now the elevation
        this_elevation = filled_raster.get_data_element(row,col);

        // Only do something if there is data
        if(this_elevation != NoDataValue)
        {
          // To get pressure, first get the lat and long
          filled_raster.get_lat_and_long_locations(row, col, lat, longitude, Converter);

          // now the pressure
          this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude), 
                                        double(this_elevation));

          // now get the scaling
          this_scaling = LSDCRNP.stone2000sp(lat,this_pressure, Fsp);
        
          // update the scaling array
          NewScaling[row][col] = this_scaling;
          NewPressure[row][col] = this_pressure;
          
          // now get combined shielding and scaling
          // first topographic shielding
          if (there_is_toposhield)
          {
            this_toposhield = Topographic_shielding.get_data_element(row,col);
            //cout << "Toposhield is: " << this_toposhield << endl;
          }
          else
          {
            this_toposhield = 1;
          }
          
          // now snow shielding
          if (there_is_snowshield)
          {

            this_snowshield = exp(-Snow_shielding.get_data_element(row,col)/gamma_spallation);
            //cout <<  "Snowshield is: " << this_snowshield << endl;
          }
          else
          {
            this_snowshield = exp(-this_Params[0]/gamma_spallation);
          }
          if (there_is_selfshield)
          {
            this_selfshield = gamma_spallation/Self_shielding.get_data_element(row,col)*
                             (1-exp(-Self_shielding.get_data_element(row,col)/gamma_spallation));
          }
          else
          {
            if(this_Params[1] == 0)
            {
              this_selfshield = 1;
            }
            else
            {
              this_selfshield = gamma_spallation/this_Params[1]*
                             (1-exp(-this_Params[1]/gamma_spallation));
            }
          }
          
          // now get the products
          
          this_shielding = this_toposhield*this_snowshield*this_selfshield;
          this_combined_scaling = this_shielding*this_scaling;
          //cout << "shield: " << this_shielding << " scale " << this_combined_scaling << endl;
          
          NewCombinedScaling[row][col]= this_combined_scaling;
          NewCombinedShielding[row][col] = this_shielding;
        }           // END test of no data
      }             // END columns loop
    }               // END rows loop
    
    // now get the rasters
    LSDRaster CombinedScalingRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, NewCombinedScaling, GeoReferencingStrings);
    LSDRaster CombinedShieldingRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, NewCombinedShielding, GeoReferencingStrings);
    LSDRaster ScalingRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, NewScaling, GeoReferencingStrings);
    LSDRaster PressureRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, NewPressure, GeoReferencingStrings);
                               
    // Format the filenames
    string CShield_ext = "_CSHIELD";
    string CScale_ext = "_CSCALE";
    string P_ext = "_PRES";
    string Prod_ext = "_PROD";
    
    // print the rasters
    PressureRaster.write_raster(DEM_fname+P_ext,bil_ext);
    ScalingRaster.write_raster(DEM_fname+Prod_ext,bil_ext);
    CombinedScalingRaster.write_raster(DEM_fname+CScale_ext,bil_ext);
    CombinedShieldingRaster.write_raster(DEM_fname+CShield_ext,bil_ext); 
  }                 // END DEM loop
  
  
  cout << endl << endl << "=====================================================" << endl;
  cout << "++WARNING++"  << endl;
  cout << "You have printed the scaling, sheilding and production rasters." << endl;
  cout << " The snow and self shielding are determined by an approximation, " << endl;
  cout << " where shielding is due to attenuation of spallation only." << endl;
  cout << " This differs from the full analysis which uses an analytical solution" << endl;
  cout << " of all the production mechanisms (incl. muons)." << endl;
  cout << "See Mudd et al (2016) ESURF for details." << endl;
  cout << "=====================================================" << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints all the basins. Smaller basins should nest within bigger
// basins
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoData::print_basins_to_for_checking()
{
  // first get the names of the DEMs:
  vector<string> dfnames = get_DEM_fnames();
  string DEM_bil_extension = "bil";
  
  int n_DEMs = dfnames.size();
  for(int i = 0; i<n_DEMs; i++)
  {
    string DEM_fname = dfnames[i];
    string basin_raster_name = DEM_fname+"_AllBasins";
    
    // now find valid points
    vector<int> valid_cosmo_points;         // a vector to hold the valid nodes
    vector<int> snapped_node_indices;       // a vector to hold the valid node indices
    vector<int> snapped_junction_indices;   // a vector to hold the valid junction indices
  
    // Load the DEM
    string DEM_bil_extension = "bil";
    LSDRaster topo_test(DEM_fname, DEM_bil_extension);
  
    // Fill this raster
    LSDRaster filled_raster = topo_test.fill(min_slope);
  
    // get the flow info
    LSDFlowInfo FlowInfo(boundary_conditions, filled_raster);

    // get contributing pixels (needed for junction network)
    LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  
    // get the sources
    vector<int> sources;
    sources = FlowInfo.get_sources_index_threshold(ContributingPixels, source_threshold);

    // now get the junction network
    LSDJunctionNetwork JNetwork(sources, FlowInfo);
    
    // print the stream order raster (to check against points)
    LSDIndexRaster SO_raster = JNetwork.StreamOrderArray_to_LSDIndexRaster();
    string SO_filename = DEM_fname+"_SO";
    SO_raster.write_raster(SO_filename,DEM_bil_extension);

    // Also print a csv of the channel nodes
    string channel_csv_name = DEM_fname+"_CN";
    JNetwork.PrintChannelNetworkToCSV(FlowInfo, channel_csv_name);


    // Now convert the data into this UTM zone
    convert_to_UTM(filled_raster);

    // convert UTM vectors to float
    vector<float> fUTM_easting;
    vector<float> fUTM_northing;
    for (int i = 0; i< int(UTM_easting.size()); i++)
    {
      fUTM_easting.push_back( float(UTM_easting[i]));
      fUTM_northing.push_back( float(UTM_northing[i]));
    }
  
    JNetwork.snap_point_locations_to_channels(fUTM_easting, fUTM_northing, 
              search_radius_nodes, threshold_stream_order, FlowInfo, 
              valid_cosmo_points, snapped_node_indices, snapped_junction_indices);

    int n_valid_points = int(valid_cosmo_points.size());
    
    string sample_string_name = DEM_fname+"_BasinKey.csv";
    
    // a file for holding a key to the basins. 
    ofstream basin_key_out;
    basin_key_out.open(sample_string_name.c_str());
    basin_key_out << "basin_ID,sample_name,basin_area" << endl;
    
    LSDIndexRaster BasinMasterRaster;   // Bow before the BasinMasterRaster
    
    // this holds all the basins. We need this because we need to know the area of the other basins
    // for the nesting. 
    vector<LSDBasin> AllTheBasins;
    
    map<int,int> drainage_of_other_basins;
    map<int,float> drainage_of_other_basins_area;

    //========================
    // LOOPING THROUGH BASINS
    //========================
    // now loop through the valid points, getting the cosmo data 
    cout << "I am trying to print basins, found " << n_valid_points << " valid points." << endl;
    for(int samp = 0; samp<n_valid_points; samp++)
    {
       
      cout << "Valid point is: " << valid_cosmo_points[samp] << " X: " 
           << UTM_easting[valid_cosmo_points[samp]] << " Y: "
           << UTM_northing[valid_cosmo_points[samp]] << endl;
      cout << "Node index is: " <<  snapped_node_indices[samp] << " and junction is: " 
           << snapped_junction_indices[samp] << " and sample name is: " << sample_name[valid_cosmo_points[samp]] << endl;
           
      cout << "Getting basin" << endl;
      LSDBasin thisBasin(snapped_junction_indices[samp],FlowInfo, JNetwork);
      AllTheBasins.push_back(thisBasin);
      
      drainage_of_other_basins[snapped_node_indices[samp]] = thisBasin.get_NumberOfCells();
      drainage_of_other_basins_area[snapped_node_indices[samp]] = thisBasin.get_Area();
      
      basin_key_out << valid_cosmo_points[samp] << "," << sample_name[valid_cosmo_points[samp]] 
                    << "," << thisBasin.get_Area() << endl;
    }
    basin_key_out.close();
    
    // now loop through everything again getting the raster
    if (n_valid_points > 0)     // this gets the first raster
    {
      BasinMasterRaster = AllTheBasins[0].write_integer_data_to_LSDIndexRaster(valid_cosmo_points[0], FlowInfo);
    }
    
    // now add on the subsequent basins
    for(int samp = 1; samp<n_valid_points; samp++)
    {
      AllTheBasins[samp].add_basin_to_LSDIndexRaster(BasinMasterRaster, FlowInfo,
                              drainage_of_other_basins, valid_cosmo_points[samp]);
    }
    
    BasinMasterRaster.write_raster(basin_raster_name, DEM_bil_extension);
    cout << "Finished with this DEM!" << endl; 
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#endif
