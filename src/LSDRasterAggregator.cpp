//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterAggregator.cpp
// Land Surface Dynamics RasterAggregator
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
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDRasterAggregator.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDRasterAggregator_CPP
#define LSDRasterAggregator_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The default create function. Doesn't do anything
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::create()
{}

void LSDRasterAggregator::create(string path_name, string param_name_prefix)
{
  // set the names of the data members
  path = path_name;
  param_name = param_name_prefix;

  // Set the parameters
  // The default slope parameter for filling. Do not change. 
  min_slope = 0.0001;

  // a boundary condition for the flow info object
  vector<string> boundary_conditionst(4);
  boundary_conditionst[0] = "n";
  boundary_conditionst[1] = "n";
  boundary_conditionst[2] = "n";
  boundary_conditionst[3] = "n";
  boundary_conditions = boundary_conditionst;

  // remove control characters from these strings
  path_name.erase(remove_if(path_name.begin(), path_name.end(), ::iscntrl), path_name.end());
  param_name_prefix.erase(remove_if(param_name_prefix.begin(), param_name_prefix.end(), ::iscntrl), param_name_prefix.end());
  
  string combined_filename = path_name+param_name_prefix;
  
  // load the file that contains the path to both the cosmo data and the 
  // DEMs
  
  // set up extensions to the files
  string rasters_ext = ".rasters";
  string params_ext = ".param";
  
  // get the filenames to open
  string rasters_fname =  path_name+param_name_prefix+rasters_ext;
  string parameters_fname = path_name+param_name_prefix+params_ext;

  // load the data. 
  cout << "Loading the names of the rasters." << endl;
  load_raster_filenames(rasters_fname);


  // check the rasters
  cout << "Hold on while I check the rasters. I don't want to eat dirty rasters!" << endl;
  check_rasters();
  cout << "The rasters seem okay. Mmmm, rasters." << endl << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads the filenames 
// of the rasters. The input file is a csv file with the first line being
// the column names and the data thereafter
// The file should have two columns, comma seperated, with the raster
// type first and the raster name second. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::load_raster_filenames(string filename)
{
  // a string for null values
  string null_str = "NULL";
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv filenames file, but the file" << filename
         << "doesn't exist"<< endl;
    exit(EXIT_FAILURE);
  }
  
  string line_from_file;
  // now loop through the rest of the lines, getting the data. 
  while( getline(ifs, line_from_file))
  {
    vector<string> this_string_vec;
  
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
    
    // The first element in the raster filenames file is the raster type, 
    // the second is the filename (with full path)
    raster_filenames[this_string_vec[0]] =  this_string_vec[1];
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks the rasters for georeferencing and scaling
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::check_rasters()
{
  // check to see if the read extension is set
  
  string bil_ext = "bil";
  string null_str = "NULL";
  string raster_ext = bil_ext;
  
  // check to make sure the read extension make sense
  if (raster_ext != "asc" && raster_ext != "flt" && raster_ext != "bil")
  {
    cout << "You didn't choose a valid raster format." << endl;
    cout << "Options are asc, bil and flt." << endl;
    cout << "You chose: " << raster_ext << endl;
    cout << "defaulting to bil, but your DEMs might not load" << endl;
    
    raster_ext = bil_ext;
  }

  // loop through the lines in the files, checking to see if the 
  // georeferencing is equivalent
  int N_rasters = int(raster_filenames.size());
  
  if(N_rasters < 1)
  {
    cout << "Sorry, I didn't find any rasters." << endl;
  }
  else if (N_rasters == 1)
  {
    cout << "You've only got one raster. Are you sure your raster file is correct?" << endl;
  }
  else
  {
    // get the information from the DEM
    // get the info from the DEM
    map<string, string >::iterator it = raster_filenames.begin();
    string this_raster_name = path+(it->second);
    cout << "The first raster is: " << this_raster_name << endl;
    cout << "All other rasters need to match the size of this raster" << endl; 
    LSDRasterInfo DEM_info(this_raster_name,raster_ext);

    // now loop through the rasters checking if they are the same size as the original raster
    for( it = raster_filenames.begin(); it != raster_filenames.end(); ++it)
    {
      this_raster_name = path+(it->second);
      LSDRasterInfo ThisRaster_info(it->second,bil_ext);
      if( ThisRaster_info!=DEM_info)
      {
        cout << "A raster " << it->second << " is not same shape and size as the DEM!" << endl;
        cout << "THIS WILL NOT WORK!!!" <<endl;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints out the various rasters and their type to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterAggregator::print_raster_names_and_types_to_screen()
{
  // Loop through the raster list printing out the types and filenames.
  for( map<string, string >::iterator it = raster_filenames.begin(); it != raster_filenames.end(); ++it)
  {
    string key = it->first;
    cout << "Raster type is: " <<it->first << " and name is: " << it->second << endl;
  }
  cout << endl;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is provided a vector of strings that contains raster types
// These are arbitrary, but we do have some naming conventions (which I'll write later,
// but some are DEM, SLOPE, CURV, HS, PROD, SHIELD, etc). 
// The function then looks for these keys in the data map and if they are all there, 
// then it returns true, if not false. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterAggregator::check_raster_types(vector<string> raster_types)
{
  cout << "Let me check the raster types for you. " << endl;
  int n_types = int(raster_types.size());
  string this_type;
  bool all_there = true;
  for (int i = 0; i< n_types; i++)
  {
    this_type = raster_types[i];
    cout << "I'm checking the type " << this_type << "; ";
    if(raster_filenames.find(this_type) == raster_filenames.end())
    {
      cout << "I am sad because I didn't find that type." << endl;
      all_there = false;
    }
    else
    {
      cout << "Good news, I found this raster type!" << endl;
    }
  }

  return all_there;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis function takes a list of raster types and crease
// a flux raster that represents the flux
// (usually either mass per time for solids or atoms per time for cosmo)
// This is done by mulitplying relevant rasters
//
// For example, for cosmo applications, you would want a mass flux of quartz
// and a concentration flux of 10Be
// 
// So the mass flux of quartz is the erosion rate (in g/m^2/time)
// times the quartz fraction
// 
// The flux of atoms is the erosion rate times the quartz fraction
// times the concentration 
// 
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDRasterAggregator::create_flux_raster(vector<string> raster_types)
{
  cout << "Let me check the raster types for you. " << endl;
  bool can_proceed = check_raster_types(raster_types);

  string raster_ext = "bil";
  float multiplier;

  vector<LSDRaster> RasterVec;
  if(can_proceed)
  {
    cout << "The rasters exist, I will continue." << endl;
    int n_types = int(raster_types.size());
    string this_type;

    // This loads the rasters
    for (int i = 0; i< n_types; i++)
    {
      string full_raster_name = path+raster_filenames[raster_types[i]];
      LSDRaster this_raster(full_raster_name,raster_ext); 
      RasterVec.push_back(this_raster);  
    }

    // now loop through the rasters, multiplying them
    int NRows = RasterVec[0].get_NRows();
    int NCols = RasterVec[0].get_NCols();
    for(int row = 0; row<NRows; row++)
    {
      for(int col = 0; col<NCols; col++)
      {
        cout << "WORKING HERE LINE 334" << endl;

      }

    }
  }
  else
  {
    cout << "You can't perform this operation since I don't have the rasters." << endl;

  }

  

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates and LSDSedimentRouting function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSedimentRouting::create(string path_name, string param_name_prefix)
{
  cout << "Hello friends. I am creating a LSDSedimentRouting object" << endl;

  // set the names of the data members
  path = path_name;
  param_name = param_name_prefix;

  // Set the parameters
  // The default slope parameter for filling. Do not change. 
  min_slope = 0.0001;

  // a boundary condition for the flow info object
  vector<string> boundary_conditionst(4);
  boundary_conditionst[0] = "n";
  boundary_conditionst[1] = "n";
  boundary_conditionst[2] = "n";
  boundary_conditionst[3] = "n";
  boundary_conditions = boundary_conditionst;

  // remove control characters from these strings
  path_name.erase(remove_if(path_name.begin(), path_name.end(), ::iscntrl), path_name.end());
  param_name_prefix.erase(remove_if(param_name_prefix.begin(), param_name_prefix.end(), ::iscntrl), param_name_prefix.end());
  
  string combined_filename = path_name+param_name_prefix;
  
  // load the file that contains the path to both the cosmo data and the 
  // DEMs
  
  // set up extensions to the files
  string rasters_ext = ".rasters";
  string params_ext = ".param";
  
  // get the filenames to open
  string rasters_fname =  path_name+param_name_prefix+rasters_ext;
  string parameters_fname = path_name+param_name_prefix+params_ext;

  // load the data. 
  cout << "Loading the names of the rasters." << endl;
  load_raster_filenames(rasters_fname);

  // check the rasters
  cout << "Hold on while I check the rasters. I don't want to eat dirty rasters!" << endl;
  check_rasters();
  cout << "The rasters seem okay. Mmmm, rasters." << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// TODO
// Currently the create function is loading the parameters for fertility etc
// into vectoers but these need to be changed to maps so we can have different
// lithologies where there is an integer index into each lithology



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks the parameter values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSedimentRouting::check_parameter_values()
{
  vector<int> lithology_indices;
  
  string this_string; 
  
  if(parameter_map.find("n_lithologies") == parameter_map.end())
  {
    cout << "You are missing the number of lithologies!!" << endl;
  }
  else
  {
    N_lithologies = atoi(parameter_map["n_lithologies"].c_str());
    
    // this sets the default lithology indices. If no lithology indicies are given
    // Lithologies are assumed to start at index 0 and the fertility, erodibility
    // etc vectors have lithology indices that correspond to their vector indices. 
    for(int i = 0; i< N_lithologies; i++)
    {
      lithology_indices.push_back(i);
    }
  }
  
  if(parameter_map.find("lithology_indices") == parameter_map.end())
  {
    cout << "You are missing the lithology indices!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["lithology_indices"]);
    
    // a temporary vector for holding the data
    vector<int> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }
    lithology_indices = temp_vec;
  }
  
  // check on the erodibility_coefficients
  if(parameter_map.find("erodibility_coefficients") == parameter_map.end())
  {
    cout << "You are missing the erodibility_coefficients!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["erodibility_coefficients"]);
    
    // a temporary vector for holding the data
    vector<float> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }
    
    if (int(temp_vec.size())!=N_lithologies)
    {
      cout << "Fatal error, you do not have the same number of erodibility_coefficients" << endl;
      cout << "as you have lithologies" << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      for(int i = 0; i< N_lithologies; i++)
      {
        int lith_index = lithology_indices[i];
        erodibility_coefficients[lith_index] = temp_vec[i];
      }
    }
  }
  
  
  // check on the fertility_coefficients
  if(parameter_map.find("fertility_coefficients") == parameter_map.end())
  {
    cout << "You are missing the fertility_coefficients!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["fertility_coefficients"]);
    
    // a temporary vector for holding the data
    vector<float> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }

    if (int(temp_vec.size())!=N_lithologies)
    {
      cout << "Fatal error, you do not have the same number of fertility_coefficients" << endl;
      cout << "as you have lithologies" << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      for(int i = 0; i< N_lithologies; i++)
      {
        int lith_index = lithology_indices[i];
        fertility_coefficients[lith_index] = temp_vec[i];
      }
    }
  }

  // check on the source_1mm
  if(parameter_map.find("source_1mm") == parameter_map.end())
  {
    cout << "You are missing the source_1mm!!" << endl;
  }
  else
  {
    // create a stringstream
    stringstream ss(parameter_map["source_1mm"]);
    
    // a temporary vector for holding the data
    vector<float> temp_vec;
    
    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );
      
      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());
      
      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());
      
      // add the string to the string vec
      temp_vec.push_back( atof(substr.c_str()) );
    }

    if (int(temp_vec.size())!=N_lithologies)
    {
      cout << "Fatal error, you do not have the same number of source_1mm" << endl;
      cout << "as you have lithologies" << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      for(int i = 0; i< N_lithologies; i++)
      {
        int lith_index = lithology_indices[i];
        source_1mm[lith_index] = temp_vec[i];
      }
    }
  }
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks the parameter values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSedimentRouting::print_parameter_values_to_screen()
{
  cout << "===================================================" << endl;
  cout << "Here are your parameter values"  << endl;
  cout << "===================================================" << endl;
  cout << "N_lithologies: " << N_lithologies << endl;
  cout << "erodibility_coefficients:" << endl;
  for (int i = 0; i< int(erodibility_coefficients.size()); i++)
  {
    cout << " " <<erodibility_coefficients[i];
  }
  cout << endl;
  cout << "fertility_coefficients:" << endl;
  for (int i = 0; i< int(erodibility_coefficients.size()); i++)
  {
    cout << " " <<fertility_coefficients[i];
  }
  cout << endl;
  cout << "source_1mm:" << endl;
  for (int i = 0; i< int(erodibility_coefficients.size()); i++)
  {
    cout << " " <<source_1mm[i];
  }
  cout << endl;
  
  cout << "===================================================" << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-








/*



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the required rasters
// The required rasters are: 
// 1) the DEM
// 2) The lithology raster
// 3) The distance to outlet raster 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<LSDRaster> LSDSedimentRouting::get_required_rasters(LSDFlowInfo& FlowInfo)
{
  int N_rasters = int(raster_filenames.size());
  cout << endl << endl << endl << "====================================" << endl;
  cout << "I am getting rasters for the sediment routing routine." << endl;
  cout << "You need at least the DEM and lithology rasters. " << endl;
  cout << "The rasters MUST be denoted by the appropriate raster type in the raster filenames file" << endl;
  cout << "DEM must be called DEM" << endl;
  cout << "Lithology -> Lithology" << endl;
  cout << "Erosion -> Erosion" << endl;
  cout << "Distance from Outlet -> DistFromOutlet" << endl;
  cout << "==================================" << endl;

  // check the rasters
  check_rasters();

  string bil_ext = "bil";
  string null_str = "NULL";
  string raster_ext;
  if(parameter_map.find("dem read extension") == parameter_map.end())
  {
    cout << "You did not set the read extension. Defaulting to bil" << endl;
    raster_ext = bil_ext;
  }
  else
  {
    raster_ext = parameter_map["dem read extension"];
  }

  map<string,LSDRaster> RasterMap;
  if (N_rasters <2)
  {
    cout << "Fatal error!" << endl;
    cout << "You need to supply a DEM and a lithology raster." << endl;
    cout << "Go back and check your raster filenames." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    
    if(raster_filenames.find("DEM") == raster_filenames.end())
    {
      cout << "Fatal error, you have not designated a raster of type DEM" << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      LSDRaster ThisRaster(raster_filenames["DEM"],bil_ext);
      RasterMap["DEM"] = ThisRaster;
    }

    if(raster_filenames.find("Lithology") == raster_filenames.end())
    {
      cout << "Fatal error, you have not designated a raster of type Lithology" << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      LSDRaster ThisRaster(raster_filenames["Lithology"],bil_ext);
      RasterMap["Lithology"] = ThisRaster;
    }

    if(raster_filenames.find("Erosion") == raster_filenames.end())
    {
      cout << "I did not find a raster of type erosion, so I will assume a constant erosion rate." << endl;
    }
    else
    {
      LSDRaster ThisRaster(raster_filenames["Erosion"],bil_ext);
      RasterMap["Erosion"] = ThisRaster;
    }

    if(raster_filenames.find("Flow_distance") == raster_filenames.end())
    {
      cout << "I did not find a raster of typeFlow_distance" << endl;
      cout << "I will compute that from the DEM now." << endl;
      
      // create a fill raster
      LSDRaster filled_raster = RasterMap["DEM"].fill(min_slope);

      LSDFlowInfo FI(boundary_conditions,filled_raster);
    
      LSDRaster DistFromOutlet = FI.distance_from_outlet();
      RasterMap["DistFromOutlet"] = DistFromOutlet;
      
    }
    else
    {
      LSDRaster ThisRaster(raster_filenames["DistFromOutlet"],bil_ext);
      RasterMap["DistFromOutlet"] = ThisRaster;
    }
  }

  
  return RasterMap;
}

//-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates the bedload from a given node
//-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDSedimentRouting::calculate_suspended_and_bedload(int node, LSDFlowInfo& FI,
                                                    vector<LSDRaster> RasterVec)
{
  // See if there is an erosion rate raster
  bool is_there_an_erosion_raster = false;
  if(raster_map.find("erosion") != raster_map.end())
  {
    cout << "I found an erosion raster!" << endl;
    is_there_an_erosion_raster = true;
  }

  // get the vectors for each lithology
  map<int,vector<float>> suspended;
  map<int,vector<float>> bedload;
  map<int,int> count;
  
  // Get all the nodes upslope of the outlet node
  vector<int> upslope_nodes = FI.get_upslope_nodes(node);
  
  // now loop through these nodes
  int n_nodes  = upslope_nodes.size();
  int outlet_row,outlet_col;
  int this_row,this_col;
  int this_lithology;
  float this_bedload;
  float this_suspended;
  FI.retrieve_current_row_and_col(node,outlet_row,outlet_col);
  
  // the  RasterVec[2] is the flow distance raster
  float distance_at_outlet_node = RasterVec[2].get_data_element(outlet_row,outlet_col);
  float distance_at_this_node;
  float distance_from_node;
  
  // Now loop through all the upslope nodes
  for (upslope_node = 0; upslope_node<n_nodes; upslope_node++)
  {
    // get the node number from the upslope_nodes vector
    this_node = upslope_nodes[upslope_node];
    
    // get the row and col of this upslope node
    FI.retrieve_current_row_and_col(this_node,this_row,this_col);
    
    // get the distance rfom outlet of this node
    distance_at_this_node = RasterVec[2].get_data_element(this_row,this_col);
    
    // get the distance between the two nodes
    distance_from_node =  distance_at_this_node-distance_at_outlet_node;
    
    // get the lithology
    this_lithology = int(RasterVec[1].get_data_element(outlet_row,outlet_col));
    
    // now calculate suspended and bedload
    // check if the lithology has been recorded before:
    if(suspended.find(this_lithology) == suspended.end())
    {
      cout "I have found the first instance of lithology " << this_lithology << endl;
      
      // if it hasn't been recorded before, set the suspended to zero
      suspended[lithology] = 0;
      bedload[lithology] = 0;
    }
    
    
    
    // now add the bedload and suspended load for this lithology to the map 
    // of bedload and suspended loads
    float this_erosion
    if(is_there_an_erosion_raster)
    {
      this_erosion = RasterVec[3].get_data_element(this_row,this_col);
    }
    else
    {
      this_erosion = erosion;
    }
    
    // first the bedoad
    bedload[lithology]=bedload[lithology] + (1-source_1mm[lithology])*this_erosion*pixel_size*pixel_size
                    *exp(-distance_from_node/1000*erodibility_coefficients[lithology]/100);
    
    // Now suspended sediment
    suspended[lithology] = suspended[lithology] + source_1mm[lithology]*this_erosion*pixel_size*pixel_size 
                           + (1-source_1mm[lithology])*this_erosion*pixel_size*pixel_size
                           *(1-exp(-distance_from_node/1000*erodibility_coefficients[lithology]/100));
    // count the number of pixels of lithology
    pixel_count[lithology] = pixel_count[lithology] + 1;

  }
}
*/

#endif