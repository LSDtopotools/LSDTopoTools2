//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDModelDriver.cpp
// Land Surface Dynamics Model(s) Driver
//
// This object parses parameter files and drives the models:
//
// LSDRasterModel
// LSDCatchmentModel (coming soon to an svn-repo near you!)
//
// Its purpose is to stop having to write a bunch of .cpp driver functions and instead
// be able to write a parameter files without compiling
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2014 Simon M. Mudd 2014
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
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
/// @author David T. Milodowski, University of Edinburgh
/// @author Martin D. Hurst, British Geological Survey
/// @author Fiona Clubb, University of Edinburgh
/// @author Stuart Grieve, University of Edinburgh
/// @author James Jenkinson, University of Edinburgh
/// @author Declan Valters, University of Manchester
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1		2015-01-15
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include "LSDModelDriver.hpp"
#include "LSDRasterModel.hpp"
#include "LSDCatchmentModel.hpp"

#ifndef LSDModelDriver_CPP
#define LSDModelDriver_CPP



/// @brief This is a class to manage running LSDTopoTools. It parses a parameter
/// file and then manages running of analyses. 
/// @details The intention of this object is to run analyses via parameter
/// files and not through numerous compiled driver functions. We eventually
/// will want some kind of 'recorder' so that any time this object
/// runs an analysis it gives a full report of what analyses were run so that
/// results are reproducable


// constructor and destructor functions
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// The default constructor. This asks the user for a pathname and
// param filename
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDModelDriver::create()
{
  std::cout << "I need a model parameter file to run. Please enter the path: " << std::endl;
  std::cin >> pathname;
  check_pathname_for_slash();

  std::cout << "Now I need a model parameter filename: " << std::endl;
  std::cin >> param_fname;
  //got_flowinfo = false;
  //got_polyfit = false;
  //got_JunctionNetwork = false;

  ingest_data(pathname, param_fname);
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This constructor runs with the path and filename
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDModelDriver::create(string pname, string fname)
{
  pathname = pname;
  check_pathname_for_slash();
  param_fname = fname;
  //got_flowinfo = false;
  //got_polyfit = false;
  //got_JunctionNetwork = false;

  ingest_data(pathname, param_fname);
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// ingest data 

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function gets all the data from a parameter file
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDModelDriver::ingest_data(string pname, string p_fname)
{
  // the full name of the file
  string full_name = pname+p_fname;

  ifstream infile;
  infile.open(full_name.c_str());
  string parameter, value, lower, lower_val;
  string bc;

  cout << "Parameter filename is: " << full_name << endl;

  // now ingest parameters
  while (infile.good())
  {
    parse_line(infile, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
    {
      lower[i] = std::tolower(parameter[i]);  // converts to lowercase
    }

    cout << "parameter is: " << lower << " and value is: " << value << endl;

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);

    if (lower == "dem read extension")
    {
      dem_read_extension = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      dem_read_extension = RemoveControlCharactersFromEndOfString(dem_read_extension);
    }
    else if (lower == "dem write extension")
    {
      dem_write_extension = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      dem_write_extension = RemoveControlCharactersFromEndOfString(dem_write_extension);
    }
    else if (lower == "write path")
    {
      write_path = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      write_path = RemoveControlCharactersFromEndOfString(write_path);
    }
    else if (lower == "write fname")
    {
      write_fname = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      write_fname = RemoveControlCharactersFromEndOfString(write_fname);
      //cout << "Got the write name, it is: "  << write_fname << endl;
    }
    else if (lower == "read path")
    {
      read_path = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      read_path = RemoveControlCharactersFromEndOfString(read_path);
      //cout << "Got the write name, it is: "  << write_fname << endl;
    }
    else if (lower == "read fname")
    {
      read_fname = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      read_fname = RemoveControlCharactersFromEndOfString(read_fname);
      //cout << "Got the read name, it is: " << read_fname << endl;
    }
    
    //=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // parameters for landscape numerical modelling
    // (LSDCatchmentModel: DAV 2015-01-14)
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    // The File Information is read at the start of this method
    // no need to duplicate it here
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Supplementary Input Files
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "hydroindex_file")
    {
      CM_support_file_names["hydroindex_file"] = atof(value.c_str());
    }
    else if (lower == "rainfall_data_file")
    {
      CM_support_file_names["rainfall_data_file"] = atof(value.c_str());
    }
    else if (lower == "grain_data_file")
    {
      CM_support_file_names["grain_data_file"] = atof(value.c_str());
    }
    else if (lower == "bedrock_data_file")
    {
      CM_support_file_names["bedrock_data_file"] = atof(value.c_str());
    }
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Numerical
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "no_of_iterations")
    {
      CM_float_parameters["no_of_iterations"] = atof(value.c_str());
    }
    else if (lower == "min_time_step")
    {
      CM_float_parameters["min_time_step"] = atof(value.c_str());
    }
    else if (lower == "max_time_step")
    {
      CM_float_parameters["max_time_step"] = atof(value.c_str());
    }
    else if (lower == "run_time_start")
    {
      CM_float_parameters["run_time_start"] = atof(value.c_str());
    }
    else if (lower == "max_run_duration")
    {
      CM_float_parameters["max_run_duration"] = atof(value.c_str());
    }
    else if (lower == "memory_limit")
    {
      CM_float_parameters["memory_limit"] = atof(value.c_str());
    }
    else if (lower == "max_time_step")
    {
      CM_float_parameters["max_time_step"] = atof(value.c_str());
    }
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Output and Save Options
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "save_interval")
    {
      CM_int_parameters["save_interval"] = atof(value.c_str());
    }
    else if (lower == "time_series_interval")
    {
      CM_int_parameters["time_series_interval"] = atof(value.c_str());
    }
    else if (lower == "elevation_file")
    {
      //bool temp_bool = (value == "true") ? true : false;
      CM_model_switches["write_elevation_file"] = atof(value.c_str()); //temp_bool;
    }
    else if (lower == "grainsize_file")
    {
      CM_model_switches["write_grainsize_file"] = atof(value.c_str());
    }
    else if (lower == "parameters_file")
    {
      CM_model_switches["write_parameters_file"] = atof(value.c_str());
    }    
    else if (lower == "flowvelocity_file")
    {
      CM_model_switches["write_flowvelocity_file"] = atof(value.c_str());
    }   
    else if (lower == "waterdepth_file")
    {
      CM_model_switches["write_waterdepth_file"] = atof(value.c_str());
    }   
    else if (lower == "timeseries_file")
    {
      CM_model_switches["write_timeseries_file"] = atof(value.c_str());
    }
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Sediment
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "transport_law")
    {
      CM_method_map["transport_law"] = atof(value.c_str());
    }
    else if (lower == "max_tau_velocity")
    {
      CM_float_parameters["max_tau_velocity"] = atof(value.c_str());
    } 
    else if (lower == "active_layer_thickness")
    {
      CM_float_parameters["active_layer_thickness"] = atof(value.c_str());
    }     
    else if (lower == "recirculate_proportion")
    {
      CM_float_parameters["recirculate_proportion"] = atof(value.c_str());
    }     
    else if (lower == "lateral_erosion_on")
    {
      CM_model_switches["lateral_erosion_on"] = atof(value.c_str());
    }     
    else if (lower == "lateral_ero_rate")
    {
      CM_float_parameters["lateral_ero_rate"] = atof(value.c_str());
    }     
    else if (lower == "edge_filter_passes")
    {
      CM_float_parameters["edge_filter_passes"] = atof(value.c_str());
    }     
    else if (lower == "cells_shift_lat")
    {
      CM_float_parameters["cells_shift_lat"] = atof(value.c_str());
    }     
    else if (lower == "max_diff_cross_chann")
    {
      CM_float_parameters["max_diff_cross_chan"] = atof(value.c_str());
    }    
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Grain Size Options
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "suspended_sediment_on")
    {
      CM_model_switches["suspended_sediment_on"] = atof(value.c_str());
    }     
    else if (lower == "fall_velocity")
    {
      CM_float_parameters["fall_velocity"] = atof(value.c_str());
    }    
    else if (lower == "grain_size_frac_file")
    {
      CM_support_file_names["grain_size_frac_file"] = atof(value.c_str());
    }      
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Hydrology and Flow
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "TOPMODEL_m_value")
    {
      CM_float_parameters["TOPMODEL_m_value"] = atof(value.c_str());
    } 
    else if (lower == "rain_data_time_step")
    {
      CM_int_parameters["rain_data_time_step"] = atof(value.c_str());
    } 
    else if (lower == "spatial_var_rain")
    {
      CM_model_switches["spatial_var_rain"] = atof(value.c_str());
    } 
    else if (lower == "input_output_diff")
    {
      CM_float_parameters["input_output_diff"] = atof(value.c_str());
    }     
    else if (lower == "min_Q_for_depth_calc")
    {
      CM_float_parameters["min_Q_for_depth_calc"] = atof(value.c_str());
    }     
    else if (lower == "max_Q_for_depth_calc")
    {
      CM_float_parameters["max_Q_for_depth_calc"] = atof(value.c_str());
    }     
    else if (lower == "depth_ero_threshold")
    {
      CM_float_parameters["depth_ero_threshold"] = atof(value.c_str());
    } 
    else if (lower == "slope_on_edge_cell")
    {
      CM_float_parameters["slope_on_edge_cell"] = atof(value.c_str());
    } 
    else if (lower == "evaporation_rate")
    {
      CM_float_parameters["evaporation_rate"] = atof(value.c_str());
    } 
    else if (lower == "courant_number")
    {
      CM_float_parameters["courant_number"] = atof(value.c_str());
    } 
    else if (lower == "hflow_threshold")
    {
      CM_float_parameters["hflow_threshold"] = atof(value.c_str());
    } 
    else if (lower == "froude_num_limit")
    {
      CM_float_parameters["froude_num_limit"] = atof(value.c_str());
    } 
    else if (lower == "mannings_n")
    {
      CM_int_parameters["mannings_n"] = atof(value.c_str());
    } 
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Vegetation
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "vegetation_on")
    {
      CM_model_switches["vegetation_on"] = atof(value.c_str());
    }    
    else if (lower == "grass_grow_rate")
    {
      CM_float_parameters["grass_grow_rate"] = atof(value.c_str());
    }   
    else if (lower == "vegetation_crit_shear")
    {
      CM_float_parameters["vegetation_crit_shear"] = atof(value.c_str());
    }   
    else if (lower == "veg_erosion_prop")
    {
      CM_float_parameters["veg_erosion_prop"] = atof(value.c_str());
    } 
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Hillslope
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "creep_rate")
    {
      CM_float_parameters["creep_rate"] = atof(value.c_str());
    } 
    else if (lower == "slope_failure_thresh")
    {
      CM_float_parameters["slope_failure_thresh"] = atof(value.c_str());
    } 
    else if (lower == "soil_erosion_rate")
    {
      CM_float_parameters["soil_erosion_rate"] = atof(value.c_str());
    } 
    else if (lower == "soil_j_mean_depends_on")
    {
      CM_model_switches["soil_j_mean_depends_on"] = atof(value.c_str());
    } 
    else if (lower == "call_muddpile_model")
    {
      CM_model_switches["call_muddpile_model"] = atof(value.c_str());
    } 
    //=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // parameters for LSD RASTER MODEL
    // (LSDRasterModel: DAV 2015-01-14)
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    // INITIAL AND BOUNDARY CONDITIONS
    else if (lower == "NRows")
    {
      RM_int_parameters["NRows"] = atof(value.c_str());
    } 
    else if (lower == "NCols")
    {
      RM_int_parameters["NCols"] = atof(value.c_str());
    } 
    else if (lower == "Resolution")
    {
      RM_int_parameters["Resolution"] = atof(value.c_str());
    } 
    else if (lower == "init_parabolic_surface")
    {
      RM_model_switches["init_parabolic_surface"] = atof(value.c_str());
    } 
    else if (lower == "init_random_fractal_surface")
    {
      RM_model_switches["init_random_fractal_surface"] = atof(value.c_str());
    } 
    else if (lower == "parabola_max_elev")
    {
      RM_float_parameters["parabola_max_elev"] = atof(value.c_str());
    } 
    else if (lower == "surface_noise_value")
    {
      RM_float_parameters["surface_noise_value"] = atof(value.c_str());
    } 
    else if (lower == "boundary_code")
    {
      RM_float_parameters["boundary_code"] = atof(value.c_str());
    } 
    
    // NUMERICAL OPTIONS
    else if (lower == "time_step")
    {
      RM_int_parameters["time_step"] = atof(value.c_str());
    } 
    else if (lower == "end_time")
    {
      RM_int_parameters["end_time"] = atof(value.c_str());
    } 
    else if (lower == "end_time_mode")
    {
      RM_method_map["end_time_mode"] = atof(value.c_str());
    } 
    
    // UPLIFT OPTIONS
    else if (lower == "uplift_mode")
    {
      RM_method_map["uplift_mode"] = atof(value.c_str());
    } 
    else if (lower == "max_uplift")
    {
      RM_float_parameters["max_uplift"] = atof(value.c_str());
    } 
    else if (lower == "tolerance")
    {
      RM_float_parameters["tolerance"] = atof(value.c_str());
    } 
    else if (lower == "print_interval")
    {
      RM_int_parameters["print_interval"] = atof(value.c_str());
    } 
    else if (lower == "periodicity")
    {
      RM_int_parameters["periodicity"] = atof(value.c_str());
    } 
    
    // FLUVIAL
    else if (lower == "fluvial_on")
    {
      RM_model_switches["fluvial_on"] = atof(value.c_str());
    } 
    else if (lower == "K")
    {
      RM_float_parameters["K"] = atof(value.c_str());
    } 
    else if (lower == "m")
    {
      RM_float_parameters["m"] = atof(value.c_str());
    } 
    else if (lower == "n")
    {
      RM_float_parameters["n"] = atof(value.c_str());
    } 
    else if (lower == "K_mode")
    {
      RM_method_map["K_mode"] = atof(value.c_str());
    } 
    
    // HILLSLOPE
    else if (lower == "hillslope_on")
    {
      RM_model_switches["hillslope_on"] = atof(value.c_str());
    } 
    else if (lower == "non-linear")
    {
      RM_model_switches["non-linear"] = atof(value.c_str());
    } 
    else if (lower == "threshold_drainage")
    {
      RM_float_parameters["threshold_drainage"] = atof(value.c_str());
    } 
    else if (lower == "D")
    {
      RM_float_parameters["D"] = atof(value.c_str());
    } 
    else if (lower == "S_c")
    {
      RM_float_parameters["S_c"] = atof(value.c_str());
    } 
    else if (lower == "D_mode")
    {
      RM_method_map["D_mode"] = atof(value.c_str());
    } 
    else if (lower == "D_amplitude")
    {
      RM_float_parameters["D_amplitude"] = atof(value.c_str());
    } 
    
    // TECTONICS
	else if (lower == "isostacy")
    {
      RM_model_switches["isostacy"] = atof(value.c_str());
    } 
	else if (lower == "flexure")
    {
      RM_model_switches["flexure"] = atof(value.c_str());
    } 
	else if (lower == "rigidity")
    {
      RM_float_parameters["rigidity"] = atof(value.c_str());
    } 
    
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Select the model type
// I.e. Catchment modelling or raster landscape (square)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=	
void LSDModelDriver::select_model_from_model_map()	
{
  if(model_map.find("LSDRasterModel") != model_map.end())
  {
    run_LSDRasterModel_components();
  }
  else if(model_map.find("LSDCatchmentModel") != model_map.end())
  {
	run_LSDCatchmentModel_components();
  }
  else
  {
	  std::cout << "You did not choose a valid model in the parameter file!" << std::endl;
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Run the appropriate model's components
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDModelDriver::run_LSDRasterModel_components()
{
	// to do
	
	
}

void LSDModelDriver::run_LSDCatchmentModel_components()
{
	// to do
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a private function that makes sure the path has a slash at the end
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDModelDriver::check_pathname_for_slash()
{
  string lchar = pathname.substr(pathname.length()-1,1);
  string slash = "/";
  //cout << "lchar is " << lchar << " and slash is " << slash << endl;

  if (lchar != slash)
  {
    std::cout << "You forgot the frontslash at the end of the path. Appending." << std::endl;
    pathname = pathname+slash;
  }
  std::cout << "The pathname is: " << pathname << std::endl;
}  
	
	
#endif	
	
	
	
	
	
	
	
	
