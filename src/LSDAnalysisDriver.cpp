//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDAnalysisDriver.cpp
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object parses parameter files and drives analysis. Its purpose is
// to stop having to write a bunch of .cpp driver functions and instead
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
#include <vector>
#include <string>
#include <map>
#include "LSDStatsTools.hpp"
#include "LSDRaster.hpp"
#include "LSDAnalysisDriver.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDChannel.hpp"
using namespace std;

#ifndef LSDAnalysisDriver_CPP
#define LSDAnalysisDriver_CPP

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// The default constructor. This asks the user for a pathname and
// param filename
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::create()
{
  cout << "I need a parameter file to run. Please enter the path: " << endl;
  cin >> pathname;
  check_pathname_for_slash();

  cout << "Now I need a parameter filename: " << endl;
  cin >> param_fname;
  got_flowinfo = false;
  got_polyfit = false;
  got_JunctionNetwork = false;

  ingest_data(pathname, param_fname);
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This constructor runs with the path and filename
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::create(string pname, string fname)
{
  pathname = pname;
  check_pathname_for_slash();
  param_fname = fname;
  got_flowinfo = false;
  got_polyfit = false;
  got_JunctionNetwork = false;

  ingest_data(pathname, param_fname);
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function gets all the data from a parameter file
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::ingest_data(string pname, string p_fname)
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
      lower[i] = tolower(parameter[i]);
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

    //=-=-=-=-=-=--=-=-=-=-
    // paramters for fill
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "min_slope_for_fill")
    {
      float_parameters["min_slope_for_fill"] = atof(value.c_str());
    }
    else if (lower == "fill_method")
    {
      method_map["fill_method"] = value;
      method_map["fill_method"] = RemoveControlCharactersFromEndOfString(method_map["fill_method"]);
    }

    //=-=-=-=-=-=--=-=-=-=-
    // paramters for hillshade
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "hs_altitude")
    {
      float_parameters["hs_altitude"] = atof(value.c_str());
    }
    else if (lower == "hs_azimuth")
    {
      float_parameters["hs_azimuth"] = atof(value.c_str());
    }
    else if (lower == "hs_z_factor")
    {
      float_parameters["hs_z_factor"] = atof(value.c_str());
    }
    else if (lower == "hs_use_fill")
    {
      //cout << "Use hs bool: " << value << endl;
      bool temp_bool = (value == "true") ? true : false;
      //cout << "Temp bool: " << temp_bool << endl;
      //bool tbool = true;
      //bool fbool = false;
      //cout << "True is " << tbool << " and false is: " << fbool << endl;
      analyses_switches["hs_use_fill"] = temp_bool;
      cout << "You have set use of the fill raster for the hillshade to " 
           << analyses_switches["hs_use_fill"] << endl;
    }


    //=-=-=-=-=-=--=-=-=-=-
    // parameters for flow info
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "boundary conditions")
    {
      // get the boundary value in lowercase
      lower_val = value;
      for (unsigned int i=0; i<value.length(); ++i)
      {
        lower_val[i] = tolower(value[i]);
      }
      vector<string> temp_bc(4);
      bc = lower_val;

      // now loop through collecting boundary conditions
      for (int i = 0; i<4; i++)
      {
        string this_bc = bc.substr(i,1);
        //cout << "Component " << i << " of the bc string: " << this_bc << endl;
        if (this_bc.find("p") != 0 && this_bc.find("b") != 0 && this_bc.find("n") != 0)
        {
          cout << "boundary condition not periodic, baselevel or noflux!" << endl;
          cout << "defaulting to no flux" << endl;
          temp_bc[i] = "n";
        }
        else
        {
          temp_bc[i] = this_bc;
        }
      }
      boundary_conditions = temp_bc;
    }

    //=-=-=-=-=-=--=-=-=-=-
    // parameters for chi
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "nodeindex fname for chi map")
    {
      support_file_names["nodeindex_fname_for_chi_map"] = atof(value.c_str());
    }
    else if (lower == "a_0")
    {
      float_parameters["A_0"] = atof(value.c_str());
    }
    else if (lower == "m_over_n")
    {
      float_parameters["m_over_n"] = atof(value.c_str());
    }
    else if (lower == "threshold_area_for_chi")
    {
      float_parameters["threshold_area_for_chi"] = atof(value.c_str());
    }



    //=-=-=-=-=-=--=-=-=-=-
    // parameters for polyfit
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "polyfit_window_radius")
    {
      float_parameters["polyfit_window_radius"] = atof(value.c_str());
      cout << "Your polyfit window radius is: "  <<  float_parameters["polyfit_window_radius"] << endl;
    }
    else if (lower == "slope_method")
    {
      method_map["slope_method"] = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      method_map["slope_method"] = RemoveControlCharactersFromEndOfString(method_map["slope_method"]);
      cout << "Your slope method is: "  <<  method_map["slope_method"] << endl;
    }


    //=-=-=-=-=-=-=-=-=-=-=-=-
    // parameters for drainage area extraction
    //=-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "drainage_area_method")
    {
      method_map["drainage_area_method"] = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      method_map["drainage_area_method"] = RemoveControlCharactersFromEndOfString(method_map["drainage_area_method"]);
    }

    //=-=-=-=-=-=-=-=-=-=-=-=-
    // parameters for single thread channel extraction
    //=-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "single_thread_channel_method")
    {
      method_map["single_thread_channel_method"] = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      method_map["single_thread_channel_method"] =
      RemoveControlCharactersFromEndOfString(method_map["single_thread_channel_method"]);
    }

    //=-=-=-=-=-=--=-=-=-=-
    // parameters for area threshold channel network
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "pixel_threshold_for_channel_net")
    {
      float_parameters["pixel_threshold_for_channel_net"] = atof(value.c_str());
    }

    //=-=-=-=-=-=--=-=-=-=-
    // parameters for landscape properties
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "root_cohesion")
    {
      float_parameters["root_cohesion"] = atof(value.c_str());
      cout << "Your root_cohesion is: "  <<  float_parameters["root_cohesion"] << endl;
    }
    else if (lower == "soil_density")
    {
      float_parameters["soil_density"] = atof(value.c_str());
      cout << "Your soil_density is: "  <<  float_parameters["soil_density"] << endl;
    }
    else if (lower == "hydraulic_conductivity")
    {
      float_parameters["hydraulic_conductivity"] = atof(value.c_str());
      cout << "Your hydraulic_conductivity is: "  <<  float_parameters["hydraulic_conductivity"] << endl;
    }
    else if (lower == "soil_thickness")
    {
      float_parameters["soil_thickness"] = atof(value.c_str());
      cout << "Your soil_thickness is: "  <<  float_parameters["soil_thickness"] << endl;
    }
    else if (lower == "tan_phi")
    {
      float_parameters["tan_phi"] = atof(value.c_str());
      cout << "Your tan_phi is: "  <<  float_parameters["tan_phi"] << endl;
    }

    //=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
    // parameters for nodata hole filling
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "nodata_hole_filling_window_width")
    {
      float_parameters["nodata_hole_filling_window_width"] = atof(value.c_str());
      cout << "Your hole_filling_window is: "  
           <<  float_parameters["nodata_hole_filling_window_width"] <<  " pixels" << endl;
    }
    
    //=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
    // parameters for masking
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "curvature_mask_threshold")
    {
      float_parameters["curvature_mask_threshold"] = atof(value.c_str());
      cout << "Your curvature_mask_threshold is: "  
           <<  float_parameters["curvature_mask_threshold"] <<  " 1/m" << endl;
    }
    else if (lower == "curvature_mask_nodataisbelowthreshold")
    {
      int_parameters["curvature_mask_nodataisbelowthreshold"] = atoi(value.c_str());
      cout << "Your curvature_mask_nodataisbelowthreshold is: "  
           <<  int_parameters["curvature_mask_nodataisbelowthreshold"] <<  "; anything other than 0 means true." << endl;
    }
    else if (lower == "mask_threshold")
    {
      float_parameters["mask_threshold"] = atof(value.c_str());
      cout << "Mask_threshold is: "  
           <<  float_parameters["mask_threshold"] <<  " (dimensions depend on raster)" << endl;
    }
    else if (lower == "mask_nodataisbelowthreshold")
    {
      int_parameters["mask_nodataisbelowthreshold"] = atoi(value.c_str());
      cout << "Your mask_nodataisbelowthreshold is: "  
           <<  int_parameters["mask_nodataisbelowthreshold"] <<  "; anything other than 0 means true." << endl;
    }


    //=-=-=-=-=-=--=-=-=-=-
    // what to write
    //-=-=-=-=-=-=-=-=-=-=-=-
    else if (lower == "write fill")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_fill"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
    }
    else if (lower == "write trimmed and nodata filled")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_trim_ndfill"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_trimmed_hole_filled"] = temp_bool;
    }
    else if (lower == "write hillshade")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_hillshade"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_hillshade"] = temp_bool;
    }
    else if (lower == "write mask threshold")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_mask_threshold"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_mask_threshold"] = temp_bool;
    }
    else if (lower == "write slope")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_slope"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_slope"] = temp_bool;
    }
    else if (lower == "write curvature")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_curvature"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_curvature"] = temp_bool;
    }
    else if (lower == "write curvature mask threshold")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_curvature_mask_threshold"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_curvature"] = temp_bool;
      raster_switches["need_curvature_mask_threshold"] = temp_bool;
    }
    else if (lower == "write planform curvature")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_planform_curvature"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_planform_curvature"] = temp_bool;
    }
    else if (lower == "write tangential curvature")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_tangential_curvature"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_tangential_curvature"] = temp_bool;
    }
    else if (lower == "write profile curvature")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_profile_curvature"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_profile_curvature"] = temp_bool;
    }
    else if (lower == "write aspect")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_aspect"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_aspect"] = temp_bool;
    }
    else if (lower == "write topographic classification")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_topographic_classification"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_topographic_classification"] = temp_bool;
    }
    else if (lower == "write drainage area")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_drainage_area"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_drainage_area"] = temp_bool;
    }
    else if (lower == "write channel net")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_channel_net"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_flowinfo"] = temp_bool;
      raster_switches["need_ContributingPixels"] = temp_bool;
      raster_switches["need_JunctionNetwork"] = temp_bool;
      raster_switches["need_sources"] = temp_bool;
      raster_switches["need_SOArray"] = temp_bool;
      raster_switches["need_JunctionIndex"] = temp_bool;
    }
    else if (lower == "write nodeindex")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_nodeindex"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_flowinfo"] = temp_bool;
      raster_switches["need_nodeindex"] = temp_bool;
    }
    else if (lower == "write single thread channel")
    {
      bool temp_bool =  (value == "true") ? true : false;
      analyses_switches["write_single_thread_channel"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_flowinfo"] = temp_bool;
      raster_switches["need_flow_distance"] = temp_bool;
      raster_switches["need_drainage_area"] = temp_bool;
    }
    else if (lower == "write chi map")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_chi_map"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_flowinfo"] = temp_bool;
      raster_switches["need_chi_map"] = temp_bool;
    }
    else if (lower == "write factor of safety at saturation")
    {
      bool temp_bool = (value == "true") ? true : false;
      analyses_switches["write_FS_sat"] = temp_bool;
      raster_switches["need_base_raster"] = temp_bool;
      raster_switches["need_fill"] = temp_bool;
      raster_switches["need_slope"] = temp_bool;
      raster_switches["need_slope_angle"] = temp_bool;
      raster_switches["need_FS_sat"] = temp_bool;
    }  
    else
    {
      cout << "Line " << __LINE__ << ": No parameter '"
           << parameter << "' expected.\n\t> Check spelling." << endl;
    }

    //cout << "Got " << lower << " and value is: " << value << endl;

  }
  infile.close();

  cout << "I'm checking to make sure the filenames are compatible now." << endl;
  check_file_extensions_and_paths();
  cout << "Ingestion of parameter file complete, and pathnames checked.\n"
       << "I am now moving on to computing the rasters. \n\n";

  compute_rasters_from_raster_switches();
  cout << "I've finished computing the rasters.\n"
       << "I am now moving on to writing the data. \n\n";
  write_rasters_from_analysis_switches();

  cout << "Well I guess I am all finished now. I sure hope you got what you wanted! Have a nice day." << endl;


}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a wrapper function that loops through the maps of raster switches
// gets the desired rasters, and then prints where necessary
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::compute_rasters_from_raster_switches()
{
  cout << "LINE 250, computing rasters" << endl;

  // read the base raster
  if(raster_switches.find("need_base_raster") != raster_switches.end())
  {
    //cout<<"LINE 255 I am loading the base raster now"<<endl;
    read_base_raster();
  }

  // get the fill raster
  if(raster_switches.find("need_fill") != raster_switches.end())
  {
    //cout<<"LINE 262 I need to compute fill!!!"<<endl;

    // check to see if the base raster is loaded
    if(map_of_LSDRasters.find("base_raster") == map_of_LSDRasters.end())
    {
      //cout << "Base raster hasn't been loaded. Loading it now." << endl;
      read_base_raster();

      // now the base raster is in the map
      fill_raster();
    }
    else
    {
      fill_raster();
    }
  }

  // get trimmed and nodata filled raster
  if(raster_switches.find("need_trimmed_hole_filled") != raster_switches.end())
  {

    if(map_of_LSDRasters.find("trimmed_hole_filled") == map_of_LSDRasters.end())
    {
      calculate_trimmed_and_nodata_filled();
    }
  }

  // get the hillshade
  if(raster_switches.find("need_hillshade") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("hillshade") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_hillshade();
    }
  }

  // get a thresholded map
  if(raster_switches.find("need_mask_threshold") != raster_switches.end())
  {
    // check to see if curvature maskhas already been calculated
    if(map_of_LSDRasters.find("mask_threshold") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_mask_threshold();
    }
  }

  // get the slope
  if(raster_switches.find("need_slope") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("slope") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_slope();
    }
  }

  // get the angle
  if(raster_switches.find("need_slope_angle") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("slope_angle") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_slope_angle();
    }
  }

  // get the aspect
  if(raster_switches.find("need_aspect") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("aspect") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }
  }

  // get the curvature
  if(raster_switches.find("need_curvature") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("curvature") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }
  }

  // get the curvature  threshold map
  if(raster_switches.find("need_curvature_mask_threshold") != raster_switches.end())
  {
    // check to see if curvature has already been calculated
    if(map_of_LSDRasters.find("curvature") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }

    // check to see if curvature maskhas already been calculated
    if(map_of_LSDIndexRasters.find("curvature_mask_threshold") == map_of_LSDIndexRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_curvature_mask_threshold();
    }
  }


  // get the planform curvature
  if(raster_switches.find("need_planform_curvature") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("planform_curvature") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }
  }

  // get the profile curvature
  if(raster_switches.find("need_profile_curvature") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("profile_curvature") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }
  }

  // get the tangential curvature
  if(raster_switches.find("need_tangential_curvature") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("tangential_curvature") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }
  }

  // get the classification
  if(raster_switches.find("need_polyfit_classification") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("polyfit_classification") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_polyfit();
    }
  }

  // get the drainage area
  if(raster_switches.find("need_drainage_area") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("drainage_area") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_drainage_area();
    }

  }

  // get the angle
  if(raster_switches.find("need_FS_sat") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("FS_sat") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_FS_sat();
    }
  }


  // get the flow info
  if(raster_switches.find("need_flowinfo") != raster_switches.end())
  {
    cout << "Hey buddy, I need to get the FlowInfo object" << endl;
    // only calculate flow info if it has not already been calculated
    if (not got_flowinfo)
    {
      calculate_flowinfo();
    }
  }

  LSDJunctionNetwork JN;

  // get the junction network
  if(raster_switches.find("need_JunctionNetwork") != raster_switches.end())
  {
    cout << "Hey buddy, I need to get the JunctionNetwork object" << endl;
    // only calculate flow info if it has not already been calculated
    if (not got_JunctionNetwork)
    {
      JN = calculate_JunctionNetwork();
    }
  }


  // get the sources
  if(raster_switches.find("need_sources") != raster_switches.end())
  {
    // check to see if it has already been calculated
    if(integer_vector_map.find("sources") ==integer_vector_map.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_sources();
    }
  }

  // get the Stream order array
  if(raster_switches.find("need_SOArray") != raster_switches.end())
  {
    //cout << "Hey buddy, I need to get the FlowInfo object" << endl;
    // only calculate flow info if it has not already been calculated
    if(map_of_LSDIndexRasters.find("SOArray") == map_of_LSDIndexRasters.end())
    {
      calculate_SOArray(JN);
    }
  }

  // get the JunctionIndex aray
  if(raster_switches.find("need_JunctionIndex") != raster_switches.end())
  {
    //cout << "Hey buddy, I need to get the FlowInfo object" << endl;
    // only calculate flow info if it has not already been calculated
    if(map_of_LSDIndexRasters.find("JunctionIndex") == map_of_LSDIndexRasters.end())
    {
      calculate_JunctionIndex(JN);
    }
  }

  // get the contributing pixels
  if(raster_switches.find("need_ContributingPixel") != raster_switches.end())
  {
    //cout << "Hey buddy, I need to get the FlowInfo object" << endl;
    // only calculate flow info if it has not already been calculated
    if(map_of_LSDIndexRasters.find("ContributingPixels") == map_of_LSDIndexRasters.end())
    {
      calculate_ContributingPixels();
    }
  }

  // check to see if you need the nodeindex raster
  if(raster_switches.find("need_nodeindex") != raster_switches.end())
  {
    //cout << "Hey buddy, I need to get the nodeindex object" << endl;

    // check to see if it has already been calculated
    if(map_of_LSDIndexRasters.find("nodeindex") == map_of_LSDIndexRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_nodeindex();
    }
  }

  // check to see if you need the chi raster
  if(raster_switches.find("need_chi_map") != raster_switches.end())
  {
    //cout << "Hey buddy, I sure can compute the chi map for ya. " << endl;
    // check to see if it has already been calculated
    if(map_of_LSDRasters.find("chi_map") == map_of_LSDRasters.end())
    {
      // it hasn't been calculated. Calculate it now.
      calculate_chi_map();
    }
  }

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Reads the base raster for analysis
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::write_rasters_from_analysis_switches()
{
  if(analyses_switches.find("write_fill") != analyses_switches.end())
  {
    //cout << "LINE 323, so you want me to write fill? Okay." << endl;

    if(map_of_LSDRasters.find("fill") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string fill_seperator = "_fill";
    string fill_fname = write_path+write_fname+fill_seperator;
    map_of_LSDRasters["fill"].write_raster(fill_fname,dem_write_extension);
  }

  // write the trimmed and no data filled
  if(analyses_switches.find("write_trim_ndfill") != analyses_switches.end())
  {
    if(map_of_LSDRasters.find("trimmed_hole_filled") == map_of_LSDRasters.end())
    {
      calculate_trimmed_and_nodata_filled();
    }
    
    string ndh_fill_seperator = "_trim_ndf";
    string ndh_fill_fname = write_path+write_fname+ndh_fill_seperator;
    map_of_LSDRasters["trimmed_hole_filled"].write_raster(ndh_fill_fname,dem_write_extension);  
  }

  // write the hillshade
  if(analyses_switches.find("write_hillshade") != analyses_switches.end())
  {
    //cout << "LINE 323, so you want me to write fill? Okay." << endl;

    if(map_of_LSDRasters.find("hillshade") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string r_seperator = "_hs";
    string r_fname = write_path+write_fname+r_seperator;
    map_of_LSDRasters["hillshade"].write_raster(r_fname,dem_write_extension);
  }

  // write curvature
  if(analyses_switches.find("write_mask_threshold") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDIndexRasters.find("mask_threshold") == map_of_LSDIndexRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_mask_threshold();
    }

    string mask_seperator = "_THMASK";
    string mask_fname = write_path+write_fname+mask_seperator;
    map_of_LSDRasters["mask_threshold"].write_raster(mask_fname,dem_write_extension);
  }


  // write slope
  if(analyses_switches.find("write_slope") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("slope") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_slope();
    }

    string slope_seperator = "_slope";
    string slope_fname = write_path+write_fname+slope_seperator;
    map_of_LSDRasters["slope"].write_raster(slope_fname,dem_write_extension);
  }

  // write slope
  if(analyses_switches.find("write_FS_sat") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("FS_sat") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_FS_sat();
    }

    string slope_seperator = "_FSsat";
    string slope_fname = write_path+write_fname+slope_seperator;
    map_of_LSDRasters["FS_sat"].write_raster(slope_fname,dem_write_extension);
  }


  // write aspect
  if(analyses_switches.find("write_aspect") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("aspect") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_polyfit();
    }

    string aspect_seperator = "_aspect";
    string aspect_fname = write_path+write_fname+aspect_seperator;
    map_of_LSDRasters["aspect"].write_raster(aspect_fname,dem_write_extension);
  }

  // write curvature
  if(analyses_switches.find("write_curvature") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("curvature") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_polyfit();
    }

    string curvature_seperator = "_curvature";
    string curvature_fname = write_path+write_fname+curvature_seperator;
    map_of_LSDRasters["curvature"].write_raster(curvature_fname,dem_write_extension);
  }

  // write curvature
  if(analyses_switches.find("write_curvature_mask_threshold") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDIndexRasters.find("curvature_mask_threshold") == map_of_LSDIndexRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_curvature_mask_threshold();
    }

    string curvature_seperator = "_curvature_mask";
    string curvature_fname = write_path+write_fname+curvature_seperator;
    map_of_LSDIndexRasters["curvature_mask_threshold"].write_raster(curvature_fname,dem_write_extension);
  }

  // write profile curvature
  if(analyses_switches.find("write_profile_curvature") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("profile_curvature") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_polyfit();
    }

    string profile_curvature_seperator = "_profile_curvature";
    string profile_curvature_fname = write_path+write_fname+profile_curvature_seperator;
    map_of_LSDRasters["profile_curvature"].write_raster(profile_curvature_fname,dem_write_extension);
  }

  // write planform curvature
  if(analyses_switches.find("write_planform_curvature") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("planform_curvature") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_polyfit();
    }

    string planform_curvature_seperator = "_planform_curvature";
    string planform_curvature_fname = write_path+write_fname+planform_curvature_seperator;
    map_of_LSDRasters["planform_curvature"].write_raster(planform_curvature_fname,dem_write_extension);
  }

  // write tangential curvature
  if(analyses_switches.find("write_tangential_curvature") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("tangential_curvature") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_polyfit();
    }

    string tangential_curvature_seperator = "_tangential_curvature";
    string tangential_curvature_fname = write_path+write_fname+tangential_curvature_seperator;
    map_of_LSDRasters["tangential_curvature"].write_raster(tangential_curvature_fname,dem_write_extension);
  }

  // write classification
  if(analyses_switches.find("write_polyfit_classification") != analyses_switches.end())
  {
    // check to see if the slope map exists
    if(map_of_LSDRasters.find("polyfit_classification") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      calculate_polyfit();
    }

    string polyfit_classification_seperator = "_polyfit_classification";
    string polyfit_classification_fname = write_path+write_fname+polyfit_classification_seperator;
    map_of_LSDRasters["polyfit_classification"].write_raster(polyfit_classification_fname,dem_write_extension);
  }


  // write nodeindex
  if(analyses_switches.find("write_nodeindex") != analyses_switches.end())
  {
    //cout << "LINE 323, so ou want me to write nodeindex? okay." << endl;
    if(map_of_LSDIndexRasters.find("nodeindex") == map_of_LSDIndexRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string NI_seperator = "_NI";
    string NI_fname = write_path+write_fname+NI_seperator;
    map_of_LSDIndexRasters["nodeindex"].write_raster(NI_fname,dem_write_extension);
  }

  // write the channel net
  // at the moment this is limited to pixel accumulation
  if(analyses_switches.find("write_channel_net") != analyses_switches.end())
  {
    if(map_of_LSDIndexRasters.find("SOArray") == map_of_LSDIndexRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }
    if(map_of_LSDIndexRasters.find("JunctionNetwork") == map_of_LSDIndexRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string SO_seperator = "_SO";
    string SO_fname = write_path+write_fname+SO_seperator;
    map_of_LSDIndexRasters["SOArray"].write_raster(SO_fname,dem_write_extension);

    string JI_seperator = "_JI";
    string JI_fname = write_path+write_fname+JI_seperator;
    map_of_LSDIndexRasters["JunctionIndex"].write_raster(JI_fname,dem_write_extension);

  }



  // write the chi map
  if(analyses_switches.find("write_chi_map") != analyses_switches.end())
  {
    if(map_of_LSDRasters.find("chi_map") == map_of_LSDRasters.end())
    {
      //cout << "You've not run the get raster routine. Running now. " << endl;
      compute_rasters_from_raster_switches();
    }

    string chi_seperator = "_chiMap";
    string chi_fname = write_path+write_fname+chi_seperator;
    map_of_LSDRasters["chi_map"].write_raster(chi_fname,dem_write_extension);
  }


}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Writes other data elements from the switches
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::write_shapes_from_switches()
{

  // write a single thread channel
  if(analyses_switches.find("write_single_thread_channel") != analyses_switches.end())
  {
    LSDChannel this_channel = get_single_thread_channel();

    // make sure flow distance exists
    if( map_of_LSDRasters.find("flow_distance") == map_of_LSDRasters.end())
    {
      calculate_flow_distance();
    }

    // write the channel
    this_channel.write_channel_to_csv(write_path,write_fname,map_of_LSDRasters["flow_distance"]);
  } 
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Reads the base raster for analysis
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::read_base_raster()
{
  //cout << "Getting base raster " << endl;
  // check to see if you've already got the base raster
  if( map_of_LSDRasters.find("base_raster") ==  map_of_LSDRasters.end())
  {
    //cout << "Base raster doesn't exist, loading" << endl;

    // you don't have it. Calculate it here.
    string full_raster_name = read_path+read_fname;
    cout <<"Reading the raster: " << endl;
    cout << full_raster_name + "." + dem_read_extension << endl;
    LSDRaster BaseRaster(full_raster_name,dem_read_extension);
    map_of_LSDRasters["base_raster"] =  BaseRaster;
  }
  else
  {
    //cout << "Hey dude, I've already got the base raster" << endl;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the fill DEM
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::fill_raster()
{
  // first check to make sure the base raster exists
  if(map_of_LSDRasters.find("base_raster") == map_of_LSDRasters.end())
  {
    //cout << "LINE 416 Base raster doesn't exist! Reading it now." << endl;
    read_base_raster();
  }


  // see if the min_slope_for_fill has been initialised
  // If not,  set to default.
  if(float_parameters.find("min_slope_for_fill") == float_parameters.end())
  {
    float_parameters["min_slope_for_fill"] = 0.0001;
  }

  // now run the fill
  // first check to see if it has already been calculated
  if(map_of_LSDRasters.find("fill") == map_of_LSDRasters.end())
  {
      cout << "Filling raster"  << endl;
  
    // check to see if a method has been designated
    if(method_map["fill_method"] == "old_fill")
    {
      cout << "Using the old fill function" << endl;
      LSDRaster temp_fill = map_of_LSDRasters["base_raster"].fill();
      map_of_LSDRasters["fill"] =  temp_fill;
    }
    if(method_map["fill_method"] == "remove_seas")
    {
      map_of_LSDRasters["base_raster"].remove_seas();
      LSDRaster temp_fill = map_of_LSDRasters["base_raster"].fill( float_parameters["min_slope_for_fill"] );
      map_of_LSDRasters["fill"] =  temp_fill;
    }
    else
    {
      //cout << "LINE 432, the base raster index is: " << base_raster_index << endl;
      LSDRaster temp_fill = map_of_LSDRasters["base_raster"].fill( float_parameters["min_slope_for_fill"] );
      map_of_LSDRasters["fill"] =  temp_fill;
    }
  }


}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


void LSDAnalysisDriver::calculate_trimmed_and_nodata_filled()
{
  // first check to make sure the base raster exists
  if(map_of_LSDRasters.find("base_raster") == map_of_LSDRasters.end())
  {
    //cout << "LINE 416 Base raster doesn't exist! Reading it now." << endl;
    read_base_raster();
  }

  // see if the min_slope_for_fill has been initialised
  // If not,  set to default.
  if(float_parameters.find("nodata_hole_filling_window_width") == float_parameters.end())
  {
    float_parameters["nodata_hole_filling_window_width"] = 1;
  }
  
  if(map_of_LSDRasters.find("trimmed_hole_filled") == map_of_LSDRasters.end())
  {
    LSDRaster temp_hole_filled = 
        map_of_LSDRasters["base_raster"].alternating_direction_nodata_fill_with_trimmer( 
                    int(float_parameters["nodata_hole_filling_window_width"]));
    map_of_LSDRasters["trimmed_hole_filled"] = temp_hole_filled;
  }
  
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the slope raster
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_slope()
{

  // first check the method
  if(method_map.find("slope_method") == method_map.end())
  {
    cout << "You did not choose a slope method. Defaulting to d8" << endl;
    method_map["slope_method"] = "d8";
  }
  if(method_map["slope_method"] != "d8" && method_map["slope_method"] != "polyfit")
  {
    cout << "You have not selected a valid slope method. Options are d8 and polyfit" << endl
         << "You chose: " << method_map["slope_method"] << endl
         << "note slope options are case sensitive. Defualting to d8" << endl;
    method_map["slope_method"] = "d8";
  }
  cout << "Calculating slope, your slope method is: " <<  method_map["slope_method"] << endl;

  if(method_map["slope_method"] == "d8")
  {
    // d8 method
    // first check to make sure the fill raster exists
    if(map_of_LSDRasters.find("fill") == map_of_LSDRasters.end())
    {
      fill_raster();
    }
    if (not got_flowinfo)
    {
      calculate_flowinfo();
    }

    // now calculate d8 slope
    map_of_LSDRasters["slope"] = FlowInfo.calculate_d8_slope(map_of_LSDRasters["fill"]);

  }
  else if(method_map["slope_method"] == "polyfit")
  {

    // first check to see if the raster has already been calculated
    if(map_of_LSDRasters.find("slope") == map_of_LSDRasters.end())
    {
      if(integer_vector_map.find("polyfit") == integer_vector_map.end())
      {
        check_polyfit();
      }

      // run the polyfit functions
      calculate_polyfit();
    }
  }

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the slope raster
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_slope_angle()
{

  // first check to see if the raster has already been calculated
  if(map_of_LSDRasters.find("slope_angle") == map_of_LSDRasters.end())
  {
    // now check to see if slope raster exists
    if(map_of_LSDRasters.find("slope") == map_of_LSDRasters.end())
    {
      // it doesn't exists. Calculate it.
      calculate_slope();
    }

    // now calculate the slope angles
    map_of_LSDRasters["slope_angle"] = map_of_LSDRasters["slope"].calculate_slope_angles();
  }

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This calculates the drainage area
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_drainage_area()
{

  // check to make sure area doesnt already exist
  if(map_of_LSDRasters.find("drainage_area") == map_of_LSDRasters.end())
  {
    // first see if a method has been assigned
    if(method_map.find("drainage_area_method") != method_map.end())
    {
      cout << "You haven't assigned a valid drainage area method" << endl
           << " These are case sensitive. Options are: " << endl
           << "d8" << endl << "dinf" << endl << "QuinnMD" << endl
           << "FreemanMD" << endl << "M2D" <<endl;
      cout << "Defaulting to dinf" << endl;
      method_map["drainage_area_method"] = "dinf";
    }

    // check to see if you need the flow info object
    if(method_map["drainage_area_method"] == "d8")
    {
      // only calculate flow info if it has not already been calculated
      if (not got_flowinfo)
      {
        calculate_flowinfo();
      }
      // now get flow area from flowinfo
      map_of_LSDRasters["drainage_area"] = FlowInfo.write_DrainageArea_to_LSDRaster();
    }
    else
    {
      // make sure you've got a filled raster
      if(map_of_LSDRasters.find("fill") == map_of_LSDRasters.end())
      {
        // now the base raster is in the map
        fill_raster();
      }

      // now calculate the various area rasters based on the method
      if(method_map["drainage_area_method"] == "M2D")
      {
        map_of_LSDRasters["drainage_area"] = map_of_LSDRasters["fill"].M2DFlow();
      }
      else if(method_map["drainage_area_method"] == "QuinnMD")
      {
        map_of_LSDRasters["drainage_area"] = map_of_LSDRasters["fill"].QuinnMDFlow();
      }
      else if(method_map["drainage_area_method"] == "FreemanMD")
      {
        map_of_LSDRasters["drainage_area"] = map_of_LSDRasters["fill"].FreemanMDFlow();
      }
      else if(method_map["drainage_area_method"] == "dinf")
      {
        map_of_LSDRasters["drainage_area"] = map_of_LSDRasters["fill"].D_inf_units();
      }
      else
      {
        map_of_LSDRasters["drainage_area"] = map_of_LSDRasters["fill"].D_inf_units();
      }
    }
  }
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This calculates the polyfit rasters
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_polyfit()
{
  // check to see if the polyfit vector has been calculated
  if(integer_vector_map.find("polyfit") == integer_vector_map.end())
  {
    check_polyfit();
  }

  if (not got_polyfit)
  {
    // now get the polyfit rasters
    if (map_of_LSDRasters.find("base_raster") == map_of_LSDRasters.end())
    {
      read_base_raster();
    }
    else
    {
      vector<LSDRaster> pfit_rasters;
      vector<int> ivs = integer_vector_map["polyfit"];

      cout << "Calculating polyfit surfaces." << endl
           << "Your designated window radius is: "  << float_parameters["polyfit_window_radius"]
           << " dx is: " <<  map_of_LSDRasters["base_raster"].get_DataResolution()
           << " and the polyfit vector is: " << endl;
      for (int pf = 0; pf<8; pf++)
      {
        cout << ivs[pf] << " ";
      }
      cout << endl;

      float dx = map_of_LSDRasters["base_raster"].get_DataResolution();
      float WR = float_parameters["polyfit_window_radius"];
      // check to make sure the window radius is not too small
      if (WR <= 2*dx)
      {
        cout << "Warning, window radius less than twice the data resolution, defaulting 2* window resolution" << endl;
        WR = 2*sqrt(2)*dx+0.001;
      }


      pfit_rasters = map_of_LSDRasters["base_raster"].calculate_polyfit_surface_metrics(WR,ivs);

      if(ivs[1] == 1)
      {
        map_of_LSDRasters["slope"] = pfit_rasters[1];
      }
      if(ivs[2] == 1)
      {
        map_of_LSDRasters["aspect"] = pfit_rasters[2];
      }
      if(ivs[3] == 1)
      {
        map_of_LSDRasters["curvature"] = pfit_rasters[3];
      }
      if(ivs[4] == 1)
      {
        map_of_LSDRasters["planform_curvature"] = pfit_rasters[4];
      }
      if(ivs[5] == 1)
      {
        map_of_LSDRasters["profile_curvature"] = pfit_rasters[5];
      }
      if(ivs[6] == 1)
      {
        map_of_LSDRasters["tangential_curvature"] = pfit_rasters[6];
      }
      if(ivs[7] == 1)
      {
        map_of_LSDRasters["polyfit_classification"] = pfit_rasters[7];
      }
    }
    got_polyfit = true;
  }
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This creates a mask for the curvature
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_curvature_mask_threshold()
{
  bool nodataisbelowthreshold = true;
  float threshold;
  
  if(float_parameters.find("curvature_mask_threshold") == float_parameters.end())
  {
    float_parameters["curvature_mask_threshold"] = 0;
  }
  threshold = float_parameters["curvature_mask_threshold"];
  
  if(int_parameters.find("curvature_mask_nodataisbelowthreshold") == int_parameters.end())
  {
    int_parameters["curvature_mask_nodataisbelowthreshold"] = 1;
    
    // this just sets the flag to 1 is it is not 0
    if(int_parameters["curvature_mask_nodataisbelowthreshold"] != 0)
    {
      int_parameters["curvature_mask_nodataisbelowthreshold"] = 1;
      nodataisbelowthreshold = true;
    }
    else
    {
      nodataisbelowthreshold = false;
    }
  }
  
  LSDIndexRaster curv_thresh = map_of_LSDRasters["curvature"].mask_to_indexraster_using_threshold(threshold,nodataisbelowthreshold);
  
  map_of_LSDIndexRasters["curvature_mask_threshold"] = curv_thresh;
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This makes a mask retaining the original data but masking above or below
// a threshold
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_mask_threshold()
{
  bool nodataisbelowthreshold = true;
  float threshold;
  
  if(float_parameters.find("mask_threshold") == float_parameters.end())
  {
    float_parameters["mask_threshold"] = 0;
  }
  threshold = float_parameters["mask_threshold"];
  
  if(int_parameters.find("mask_nodataisbelowthreshold") == int_parameters.end())
  {
    int_parameters["mask_nodataisbelowthreshold"] = 1;
    
    // this just sets the flag to 1 is it is not 0
    if(int_parameters["curvature_mask_nodataisbelowthreshold"] != 0)
    {
      int_parameters["mask_nodataisbelowthreshold"] = 1;
      nodataisbelowthreshold = true;
    }
    else
    {
      nodataisbelowthreshold = false;
    }
  }
  
  LSDRaster mask_thresh = map_of_LSDRasters["base_raster"].mask_to_nodata_using_threshold(threshold,nodataisbelowthreshold);
  
  map_of_LSDRasters["mask_threshold"] = mask_thresh;
}



//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the hillshade raster
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_hillshade()
{
  // see if the parameters for hillshade hace been set. If not use default values
  if(float_parameters.find("hs_altitude") == float_parameters.end())
  {
    float_parameters["hs_altitude"] = 45;
  }
  if(float_parameters.find("hs_azimuth") == float_parameters.end())
  {
    float_parameters["hs_azimuth"] = 315;
  }
  if(float_parameters.find("hs_z_factor") == float_parameters.end())
  {
    float_parameters["hs_z_factor"] = 1;
  }

  // see if you need the fill
  if(analyses_switches["hs_use_fill"] == true)
  {
    // first check to make sure the fill raster exists
    if(map_of_LSDRasters.find("fill") == map_of_LSDRasters.end())
    {
      cout << "The fill raster doesn't exist. Calculating." << endl;
      fill_raster();
    }
    
    // now run the hillshade
    // first check to see if it has already been calculated
    if(map_of_LSDRasters.find("hillshade") == map_of_LSDRasters.end())
    {
      cout << "Hillshading from fill raster " << endl;
      LSDRaster temp_hs = map_of_LSDRasters["fill"].hillshade(
                              float_parameters["hs_altitude"],
                              float_parameters["hs_azimuth"],
                              float_parameters["hs_z_factor"] );

      map_of_LSDRasters["hillshade"] =  temp_hs;
    }
  }
  else
  {
    if(map_of_LSDRasters.find("hillshade") == map_of_LSDRasters.end())
    {
      cout << "Hillshading from base raster " << endl;
      LSDRaster temp_hs = map_of_LSDRasters["base_raster"].hillshade(
                              float_parameters["hs_altitude"],
                              float_parameters["hs_azimuth"],
                              float_parameters["hs_z_factor"] );

      map_of_LSDRasters["hillshade"] =  temp_hs;
    }    
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the factor of safety when the hillslopes are saturated
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_FS_sat()
{

  // check to see if this has already been calculated
  if(map_of_LSDRasters.find("FS_sat") == map_of_LSDRasters.end())
  {
    // it hasn't. calculate it

    // check to see if the slope_angle raster is present
    if(map_of_LSDRasters.find("slope_angle") == map_of_LSDRasters.end())
    {
      calculate_slope_angle();
    }


    // see if the parameters for FS calculation have been set. If not use default values
    if(float_parameters.find("root_cohesion") == float_parameters.end())
    {
      cout << "You didn't define root_cohesion. Defaulting to 10000 N/m^2" << endl;
      float_parameters["root_cohesion"] = 10000;
    }
    if(float_parameters.find("soil_density") == float_parameters.end())
    {
      cout << "You didn't define soil_density, degfaulting to 1300 kg/m^3" << endl;
      float_parameters["soil_density"] = 1300;
    }
    if(float_parameters.find("soil_thickness") == float_parameters.end())
    {
      cout << "You didn't define soil_thickness, defaulting to 1m" << endl;
      float_parameters["soil_thickness"] = 1;
    }
    if(float_parameters.find("tan_phi") == float_parameters.end())
    {
      cout << "You didn't define tan_phi, defaulting to 0.8" << endl;
      float_parameters["tan_phi"] = 0.8;
    }

    map_of_LSDRasters["FS_sat"]
     = map_of_LSDRasters["SlopeAngle"].calculate_factor_of_safety_at_saturation(
            float_parameters["root_cohesion"], float_parameters["soil_density"],
            float_parameters["soil_thickness"], float_parameters["tan_phi"],
            map_of_LSDRasters["SlopeAngle"]);
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the flow info object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_flowinfo()
{
  // first check to see if this already exists
  if (not got_flowinfo)
  {
    // it doens't exist. Calculate it here.
    //cout << "LINE 457 Flow info doesn't exist. Getting it from the fill raster" << endl;
    // this requires the fill raster. See if it exists
    if(map_of_LSDRasters.find("fill") == map_of_LSDRasters.end())
    {
      //cout << "LINE 461, fill hasn't been computed yet, getting it." << endl;
      // it doesn't exit. Calculate it.
      fill_raster();
    }

    // now the fill exists. Check to see if boundary conditions exist
    check_boundary_conditions();

    // okay, everything should be ready for flow info calculation
    LSDFlowInfo temp_FI(boundary_conditions, map_of_LSDRasters["fill"]);

    // set the data members
    FlowInfo = temp_FI;
    got_flowinfo = true;
    //cout << "LINE 715, got flowinfo" << endl;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the flow distance
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_flow_distance()
{
  // first check to see if this already exists
  if (map_of_LSDRasters.find("flow_distance") ==  map_of_LSDRasters.end())
  {
    // it doens't exist. Calculate it here.
    if(not got_flowinfo)
    {
      calculate_flowinfo();
    }

    map_of_LSDRasters["flow_distance"] = FlowInfo.distance_from_outlet();
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the junction network object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDJunctionNetwork LSDAnalysisDriver::calculate_JunctionNetwork()
{
  // first check to see if this already exists
  //if (not got_JunctionNetwork)
  //{
    // this needs flowInfo
    if(not got_flowinfo)
    {
      //cout << "LINE 733, fill hasn't been computed yet, getting it." << endl;
      // it doesn't exit. Calculate it.
      calculate_flowinfo();
    }

    // it also needs sources
    if(map_of_LSDIndexRasters.find("sources") == map_of_LSDIndexRasters.end())
    {
      //cout << "LINE 741, you haven't got the sources yet, getting them." << endl;
      // you don't have the sources. Calculate them
      calculate_sources();
    }

    // okay, everything should be ready for flow info calculation
    LSDJunctionNetwork temp_JN(integer_vector_map["sources"], FlowInfo);



    got_JunctionNetwork = true;
    return temp_JN;
    //}
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the nodeindex object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_nodeindex()
{
  // see if you've already got this raster
  if(map_of_LSDIndexRasters.find("nodeindex") == map_of_LSDIndexRasters.end())
  {
    // you don't have it. Calculate it here.
    // this requires the flow info object. See if it exists
    if(not got_flowinfo)
    {
      // it doesn't exit. Calculate it.
      calculate_flowinfo();
    }

    // now you have the LSDFlowInfo object. Spit out the nodeindex raster
    LSDIndexRaster temp_NI = FlowInfo.write_NodeIndex_to_LSDIndexRaster();
    map_of_LSDIndexRasters["nodeindex"] = temp_NI;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the contributing pixels object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_ContributingPixels()
{
  // see if you've already got this raster
  if(map_of_LSDIndexRasters.find("ContributingPixels") == map_of_LSDIndexRasters.end())
  {
    // you don't have it. Calculate it here.
    // this requires the flow info object. See if it exists
    if(not got_flowinfo)
    {
      // it doesn't exit. Calculate it.
      calculate_flowinfo();
    }

    // now you have the LSDFlowInfo object. Spit out the nodeindex raster
    //cout << "Line 799, getting CP" << endl;
    LSDIndexRaster temp_CP = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

    //cout << "Rows: " << temp_CP.get_NRows() << endl;
    map_of_LSDIndexRasters["ContributingPixels"] =  temp_CP;
    //cout << "Line 802, got CP" << endl;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets a single thread channel
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDChannel LSDAnalysisDriver::get_single_thread_channel()
{

  LSDChannel return_channel;

  // check to see if the method has been set
  if(method_map.find("single_thread_channel_method") == method_map.end())
  {
    method_map["single_thread_channel_method"] = "start_and_end_node";
  }

  // now check to see if rasters and flow info exist
  if( map_of_LSDRasters.find("base_raster") ==  map_of_LSDRasters.end())
  {
    read_base_raster();
  }
  if( map_of_LSDRasters.find("fill") ==  map_of_LSDRasters.end())
  {
    fill_raster();
  }
  if( map_of_LSDRasters.find("drainage_area") ==  map_of_LSDRasters.end())
  {
    calculate_drainage_area();
  }
  if (not got_flowinfo)
  {
    calculate_flowinfo();
  }

  // now calculate channel based on these inputs.
  if (method_map["single_thread_channel_method"] == "start_and_end_node")
  {
    // first check to see if nodes exist
    if(float_parameters.find("starting_channel_node") == float_parameters.end())
    {
      cout << "You are trying to extract a channel but haven't assigned a starting" << endl
           << "channel node."  << endl
           << "Note that the float paramter keyword is starting_channel_node." << endl
           << "Enter the starting node now: " << endl;
      cin >> float_parameters["starting_channel_node"];
      cout << endl;
    }
    // first check to see if nodes exist
    if(float_parameters.find("ending_channel_node") == float_parameters.end())
    {
      cout << "You are trying to extract a channel but haven't assigned an ending" << endl
           << "channel node."  << endl
           << "Note that the float paramter keyword is ending_channel_node." << endl
           << "Enter the ending node now: " << endl;
      cin >> float_parameters["ending_channel_node"];
      cout << endl;
    }
    // now make sure the nodes are not out of bounds
    int NNodes = FlowInfo.get_NDataNodes();
    int SN = int(float_parameters["starting_channel_node"]);
    int EN = int(float_parameters["ending_channel_node"]);
    if(SN > NNodes-1)
    {
      cout << "You are extracting a channel profile but your starting node" << endl
           << "is not on the DEM!"  << endl;
      SN = 0;
    }
    if(EN > NNodes-1)
    {
      cout << "You are extracting a channel profile but your ending node" << endl
           << "is not on the DEM!"  << endl;
      EN = 0;
    }

    // now check to see if the A_0 and m/n paramters are assigned
    // now see if you've got the parameters
    if(float_parameters.find("A_0") == float_parameters.end())
    {
      // you don't have A_0. Replace with default
      float_parameters["A_0"] = 1000;
    }
    if(float_parameters.find("m_over_n") == float_parameters.end())
    {
      // you don't have m_over_n. Replace with default
      cout << "You haven't assigned m/n. Defaulting to 0.45" << endl;
      float_parameters["m_over_n"] = 0.45;
    }

    float downslope_chi = 0;
    LSDChannel this_channel(SN, EN, downslope_chi,float_parameters["m_over_n"],
                      float_parameters["m_over_nA_0"], FlowInfo,
                      map_of_LSDRasters["fill"], map_of_LSDRasters["drainage_area"]);

    return_channel = this_channel;
  }

  return return_channel;
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates sources for channel network
// THIS NEEDS WORK TO ACCOMODATE PELLETIER AND DREICH METHODS
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_sources()
{
  // see if you've already got this raster
  if(map_of_LSDIndexRasters.find("sources") == map_of_LSDIndexRasters.end())
  {
    // you don't have it. Calculate it here.
    // this requires the flow info object. See if it exists
    if(not got_flowinfo)
    {
      // it doesn't exit. Calculate it.
      calculate_flowinfo();
    }

    // it also requires contributing pixels. See if that exists
    if(map_of_LSDIndexRasters.find("ContributingPixels") == map_of_LSDIndexRasters.end())
    {
      // it doesn't exist. Calculate it
      cout << "LINE 826, you've not got the contributing pixels yet. Getting them." << endl;
      calculate_ContributingPixels();
    }

    // set to default if no parameters
    if(float_parameters.find("pixel_threshold_for_channel_net") == float_parameters.end())
    {
      float_parameters["pixel_threshold_for_channel_net"] = 10;
    }

    int thres =   int(float_parameters["pixel_threshold_for_channel_net"]);
    vector<int> temp_sources = FlowInfo.get_sources_index_threshold(
        map_of_LSDIndexRasters["ContributingPixels"],thres);

    integer_vector_map["sources"] = temp_sources;
    //cout << "LINE 841 Got the sources" << endl;

  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the stream order array object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_SOArray(LSDJunctionNetwork& JunctionNetwork)
{
  // see if you've already got this raster
  if(map_of_LSDIndexRasters.find("SOArray") == map_of_LSDIndexRasters.end())
  {
    // you don't have it. Calculate it here.
    // this requires the flow info object. See if it exists
    //if(not got_JunctionNetwork)
    //{
    //  // it doesn't exit. Calculate it.
    //  calculate_JunctionNetwork();
    //}

    // now you have the LSDFlowInfo object. Spit out the nodeindex raster
    LSDIndexRaster temp_SO = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
    map_of_LSDIndexRasters["SOArray"] = temp_SO;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates the JunctionIndex array object
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_JunctionIndex(LSDJunctionNetwork& JunctionNetwork)
{
  // see if you've already got this raster
  if(map_of_LSDIndexRasters.find("JunctionIndex") == map_of_LSDIndexRasters.end())
  {
    // you don't have it. Calculate it here.
    // this requires the flow info object. See if it exists
    //if(not got_JunctionNetwork)
    //{
      // it doesn't exit. Calculate it.
    //  calculate_JunctionNetwork();
    //}

    // now you have the LSDFlowInfo object. Spit out the nodeindex raster
    LSDIndexRaster temp_JI = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
    map_of_LSDIndexRasters["JunctionIndex"] =  temp_JI;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculates chi map
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::calculate_chi_map()
{
  //cout << "I'm calculating chi now LINE 503" << endl;
  LSDRaster temp_chi;

  // see if you've already got this raster
  if(map_of_LSDRasters.find("chi_map") == map_of_LSDRasters.end())
  {
    //cout << "I don't have the chi_map LINE 509" << endl;
    // you don't have it. Calculate it here.
    // this requires the flow info object. See if it exists
    if(not got_flowinfo)
    {
      //cout << "I don't have the flow info line 514" << endl;
      // it doesn't exit. Calculate it.
      calculate_flowinfo();
    }

    // now see if you've got the parameters
    if(float_parameters.find("A_0") == float_parameters.end())
    {
      // you don't have A_0. Replace with default
      float_parameters["A_0"] = 1000;
    }
    if(float_parameters.find("m_over_n") == float_parameters.end())
    {
      // you don't have m_over_n. Replace with default
      cout << "You haven't assigned m/n. Defaulting to 0.45" << endl;
      float_parameters["m_over_n"] = 0.45;
    }
    if(float_parameters.find("threshold_area_for_chi_map") == float_parameters.end())
    {
      cout << "You haven't assigned the threshold_area_for_chi_map."
           << "Defaulting to 0" << endl;
      float_parameters["threshold_area_for_chi"] = 0;
    }

    // now check to see if the user has given a nodeindex file
    if(support_file_names.find("threshold_area_for_chi_map") == support_file_names.end())
    {
      // no nodindex file, calculate chi assuming all baselevel nodes have
      // a chi = 0
      temp_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(
                                 float_parameters["m_over_n"],
                                 float_parameters["A_0"],
                                 float_parameters["threshold_area_for_chi"]);
    }
    else
    {
      string sf_name =  support_file_names["threshold_area_for_chi_map"];
      // add the path
      string full_sf_name = read_path+ sf_name;
      cout << "reading nodeindices for chi map from file: "
           << full_sf_name << endl;
      ifstream NI_in;
      NI_in.open(full_sf_name.c_str());
      int this_NI;
      vector<int> NI_for_chi_map;

      while(NI_in >> this_NI)
      {
        NI_for_chi_map.push_back(this_NI);
      }
      NI_in.close();

      if (NI_for_chi_map.size() == 0)
      {
        cout << "Calculating chi map. Your node index file appears to be empty\n"
             << "Running chi from all baselevel nodes" << endl;
        temp_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(
                                 float_parameters["m_over_n"],
                                 float_parameters["A_0"],
                                 float_parameters["threshold_area_for_chi"]);
      }
      else
      {
        temp_chi = FlowInfo.get_upslope_chi_from_multiple_starting_nodes(NI_for_chi_map,
                                 float_parameters["m_over_n"],
                                 float_parameters["A_0"],
                                 float_parameters["threshold_area_for_chi"]);
      }
    }

    // Copy the chi raster to the vector of rasters
    map_of_LSDRasters["chi_map"] =  temp_chi;
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a private function that makes sure the path has a slash at the end
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_pathname_for_slash()
{
  string lchar = pathname.substr(pathname.length()-1,1);
  string slash = "/";
  //cout << "lchar is " << lchar << " and slash is " << slash << endl;

  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;
    pathname = pathname+slash;
  }
  cout << "The pathname is: " << pathname << endl;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a private function that makes sure the path has a slash at the end
// overloaded to take an arbitrary string
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string LSDAnalysisDriver::check_pathname_for_slash(string this_pathname)
{
  string lchar = this_pathname.substr(this_pathname.length()-1,1);
  string slash = "/";
  //cout << "Checking pathname, pathname is: " << this_pathname << endl;
  //cout << "lchar is " << lchar << " and slash is " << slash << endl;

  if (lchar != slash)
  {
    //cout << "You forgot the frontslash at the end of the path. Appending." << endl;
    this_pathname = this_pathname+slash;
  }
  //cout << "The pathname is: " << pathname << endl;
  return this_pathname;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this is a private function that makes sure the path has a slash at the end
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_boundary_conditions()
{
  if( int(boundary_conditions.size()) != 4)
  {
    cout << "Boundary conditions not assigned! Defaulting to no flux."  << endl;
    vector<string> temp_bc(4);
    for (int i = 0; i< 4; i++)
    {
      temp_bc[i] = "n";
    }
    boundary_conditions = temp_bc;
  }

  for (int i =0; i< 4; i++)
  {
    cout << "Boundary["<<i<<"]: "<<boundary_conditions[i]<< endl;
  }


}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function checks the file extensions for reading and writing DEMs
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_file_extensions_and_paths()
{

  // first check the extensions
  if (dem_read_extension != "asc"  && dem_read_extension != "flt" && dem_read_extension != "bil" &&
      dem_write_extension != "asc"  && dem_write_extension != "flt" && dem_write_extension != "bil")
  {
    cout << "Line 1006 Raster file extension not assigned! Defaulting to bil format." << endl;
    cout << "You entered: " << dem_read_extension << "!" <<endl;
    dem_read_extension = "bil";
    dem_write_extension = "bil";
  }
  else
  {
    if (dem_read_extension != "asc"  && dem_read_extension != "flt" && dem_read_extension != "bil")
    {
      cout << "DEM read extension not assigned, defaulting to write extension." << endl;
      dem_read_extension = dem_write_extension;
    }
    else
    {
      //cout << "DEM write extension not assigned, defaulting to read extension." << endl;
      dem_write_extension = dem_read_extension;
    }
  }

  // now check the paths
  //cout << "Write path length is: " << write_path.length() << endl;
  if (write_path.length() == 0)
  {
    write_path = pathname;
    if (read_path.length() != 0)
    {
      write_path = read_path;
    }
  }

  //cout << "CHECKING NAMES, Write fname is: " << write_fname << endl;
  //cout << "The write fname length is " << write_fname.length() << endl;
  if (write_fname.length() == 0)
  {
    if (read_fname.length() != 0)
    {
      write_fname = read_fname;
    }
    write_fname = get_string_before_dot(param_fname);
    //cout << "Write fname not assigned, defaulting to name of parameter file." << endl;
    //cout << "The write fname is: " << write_fname << endl;
  }

  // now check the path
  //cout << "Read path length is: " << read_path.length() << endl;
  if (read_path.length() == 0)
  {
    read_path = write_path;
  }
  if (read_fname.length() == 0)
  {
    read_fname = get_string_before_dot(param_fname);
    //cout << "Read fname not assigned, defaulting to name of parameter file." << endl;
    //cout << "The read fname is: " << read_fname << endl;
  }

  // make sure the read and write paths have the slash at the end
  write_path = check_pathname_for_slash(write_path);
  read_path = check_pathname_for_slash(read_path);

  cout << "The full read fname is:\n " << read_path+read_fname << endl;
  cout << "The full write fname is:\n " << write_path+write_fname << endl;
  cout << "The read and write extensions are: " << dem_read_extension
       << " " << dem_write_extension << endl;

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function checks to see which polyfit methods are needed and then
// creates the correct polyfit vector
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDAnalysisDriver::check_polyfit()
{
  // the vector used for the polyfit switches
  vector<int> polyfit_switches(8,0);

  // need to check for polyfit window radius here
  if(float_parameters.find("polyfit_window_radius") == float_parameters.end())
  {
      // you don't have the polyfit radius. Replace with default
      float_parameters["polyfit_window_radius"] = 7;
  }
  if (raster_switches["need_slope"] == true)
  {
    // first check the method
    if(method_map.find("slope_method") == method_map.end())
    {
      cout << "You did not choose a slope method. Defaulting to d8" << endl;
      method_map["slope_method"] = "d8";
    }
    if(method_map["slope_method"] != "d8" && method_map["slope_method"] != "polyfit")
    {
      cout << "You have not selected a valid slope method. Options are d8 and polyfit" << endl
           << "note these options are case sensitive. Defualting to d8" << endl;
      method_map["slope_method"] = "d8";
    }
    if(method_map["slope_method"] == "polyfit")
    {
      polyfit_switches[1] = 1;
    }
  }
  if (raster_switches["need_aspect"] == true)
  {
    polyfit_switches[2] = 1;
  }
  if (raster_switches["need_curvature"] == true)
  {
    polyfit_switches[3] = 1;
  }
  if (raster_switches["need_planform_curvature"] == true)
  {
    polyfit_switches[4] = 1;
  }
  if (raster_switches["need_profile_curvature"] == true)
  {
    polyfit_switches[5] = 1;
  }
  if (raster_switches["need_tangential_curvature"] == true)
  {
    polyfit_switches[6] = 1;
  }
  if (raster_switches["need_classification"] == true)
  {
    polyfit_switches[7] = 1;
  }

  integer_vector_map["polyfit"] = polyfit_switches;
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function strips the text after the final dot in a string
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string LSDAnalysisDriver::get_string_before_dot(string this_string)
{
  string cut_string;
  unsigned found = this_string.find_last_of(".");
  cut_string = this_string.substr(0,found);
  return cut_string;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#endif
