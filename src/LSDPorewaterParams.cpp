//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDPorewaterParams
// Land Surface Dynamics PorewterParams object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object interfaces with teh porewater column object
//  In a landsacpe each pixel will have its own pore ressures but 
//  the parameters will be constant (or have similar statistical properties)
//  across a landscape. This structure tries to minimize memory requirements
//
// Developed by:
//  Simon M. Mudd
//  Stuart W.D. Grieve
//
// Copyright (C) 2016 Simon M. Mudd 2013 6
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


/** @file LSDPorewaterColumn.hpp
@author Simon M. Mudd, University of Edinburgh
@author Stuart W. D. Grieve, University of Edinburgh

**/



#ifndef LSDPorewaterParams_CPP
#define LSDPorewaterParams_CPP

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include "LSDPorewaterParams.hpp"
#include "LSDParameterParser.hpp"
#include "LSDRaster.hpp" // addign that to get the slope raster
#include "LSDStatsTools.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Empty create function
// Starts with some defaults. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::create()
{
  D_0 = 0.00001;
  d = 2;
  alpha = 0.1;
  Iz_over_K_steady = 0.2;       // Iz is infiltration so this is the steady (long term) infiltration rate divided by the K 
  K_sat = 0.0000001;
  friction_angle = 0.66322511;
  cohesion = 500;

  // This is graviy times water and soil density
  weight_of_soil = 19000;
  weight_of_water = 9800;

  vector<float> this_z;
  for(int i = 0; i<31; i++)
  {
    this_z.push_back(float(i)*0.1);
  }
  
  // YOu get divide by zero if you have a zero at first depth, so give it finite depth
  this_z[0]= 0.05;
  
  Depths = this_z;

  calculate_beta();
  calculate_D_hat();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function reads from a parameter file. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::create(string paramfile_path, string paramfile_name)
{
  LSDParameterParser LSDPP(paramfile_path, paramfile_name);

  saving_prefix = LSDPP.get_write_fname();
  string load_path = LSDPP.get_read_path();
  
  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;
  
  // set default float parameters
  float_default_map["D_0"] = 0.00001;
  float_default_map["K_sat"] = 0.0000001;
  float_default_map["d"] = 2;
  float_default_map["Iz_over_K_steady"] = 0.2;
  float_default_map["depth_spacing"] = 0.1;
  float_default_map["alpha"] = 0.1;
  float_default_map["friction_angle"] = 0.66322511;
  float_default_map["cohesion"] = 500;
  float_default_map["weight_of_soil"] = 19000;
  float_default_map["weight_of_water"] = 9800;
  
  // set default in parameter
  int_default_map["n_depths"] = 31;
  
  bool_default_map["use_depth_vector"] = false;

  bool_default_map["full_1D_output"] = false;
  
  string depth_vector_key = "depth_vector";
  string_default_map["depth_vector"] = "";
  string_default_map["rainfall_csv"] = "";


  // parameter for 2D columning
  bool_default_map["spatial_analysis"] = false; // activate the slope analysis
  bool_default_map["reload_alpha"] = false; // load a slope raster rather than calculating it
  bool_default_map["resample_slope"] = false; // load a slope raster rather than calculating it
  float_default_map["resample_slope_res"] = false; // load a slope raster rather than calculating it
  string_default_map["topo_raster"] = ""; // name of the raster without the .bil extension
  string_default_map["alpha_to_reload"] = ""; // name of the raster without the .bil extension
  float_default_map["polyfit_window_radius"] = 6; // radius to calculate the slope
  int_default_map["n_threads"] = 4; //  number of threads for multiprocessing
  float_default_map["time_of_spatial_analysis"] = 100000; // Time in second (same than in the preprocessed input csv file) of the spatial analysis

  
  // Get the parameters
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();
  
  cout << "Yo! Iz_over_K_steady is: " << this_float_map["Iz_over_K_steady"] << endl;
    
  D_0 = this_float_map["D_0"];
  d = this_float_map["d"];
  alpha = this_float_map["alpha"];
  Iz_over_K_steady = this_float_map["Iz_over_K_steady"];
  K_sat = this_float_map["K_sat"];
  friction_angle = this_float_map["friction_angle"];
  cohesion = this_float_map["cohesion"];
  weight_of_soil = this_float_map["weight_of_soil"];
  weight_of_water = this_float_map["weight_of_water"];
  rainfall_csv_name = this_string_map["rainfall_csv"];
  rainfall_csv_path = paramfile_path; // I assume here this is the same path for your csv file. It should be anyway.
  full_1D_output = this_bool_map["full_1D_output"];
  output_2D = this_bool_map["spatial_analysis"];
  n_threads = this_int_map["n_threads"];
  time_of_spatial_analysis = this_float_map["time_of_spatial_analysis"];

  cout << "In the params, Iz_over_K_steady: " << Iz_over_K_steady << endl;

  // parameters for getting the depths
  float depth_spacing = this_float_map["depth_spacing"];
  int n_depths = this_int_map["n_depths"];
  
  vector<float> these_extracted_depths = LSDPP.parse_float_vector(depth_vector_key);
  vector<float> these_calculated_depths;
  for(int i = 0; i<n_depths; i++)
  {
    these_calculated_depths.push_back(float(i)*depth_spacing);
  }
  
  if (these_extracted_depths.size() == 0)
  {
    these_extracted_depths = these_calculated_depths;
  }
  
  // If we are going to use the depth vector, read it from the 
  // comma separated list
  if( this_bool_map["use_depth_vector"]  )
  {
    Depths = these_extracted_depths;
  }
  else
  {
    Depths = these_calculated_depths;
  }
  
  // Make sure the first depth is not 0
  if(Depths.size() == 0)
  {
    Depths.push_back(0.1);
  }
  else if(Depths.size() == 1)
  {
    Depths[0] = 0.1;
  }
  else
  {
    if(Depths[0] == 0)
    {
      Depths[0] = 0.5*Depths[1];
    }
  }
  
  calculate_beta();
  calculate_D_hat();
  
  if(this_bool_map["spatial_analysis"])
  {
    if(this_bool_map["reload_alpha"] == false)
    {
      cout << "I am now preprocessing your data  for spatial analysis of soil columns" << endl;
      
      LSDRaster this_topo(load_path + this_string_map["topo_raster"], "bil");

      vector<int> relevant_variable_name_14 = {0,1,0,0,0,0,0,0}; // selection of raster to calculate with the polyfit fitting function, with a relevant variable name LOLz
      vector<LSDRaster> tempRast;
      cout << "Let me calculate the slope by fitting a polyfit raster to it." << endl;
      tempRast = this_topo.calculate_polyfit_surface_metrics(this_float_map["polyfit_window_radius"], relevant_variable_name_14);
      alpha_raster = tempRast[1]; // getting the right slpe and pushing it to the parameter
      // Need to veonvert that into radian

      cout << "Alright I got the slope, now I am converting it to radian (IMPORANT)" << endl;
      float res_of_rast = alpha_raster.get_DataResolution();
      size_t nrows = alpha_raster.get_NRows();
      size_t ncols = alpha_raster.get_NCols();
      for(size_t i=0; i<nrows; i++)
      {
        for(size_t j=0; j<ncols;j++)
        {
          alpha_raster.set_data_element(i,j, abs(atan(alpha_raster.get_data_element(i,j))));
        }
      }
      cout << "Done, saving it ..." << endl;
      alpha_raster.write_raster(load_path+saving_prefix + "slope_in_radian", "bil");
      cout << "Saved." << endl;

    }
    else
    {
      // Reloading the slope from previous analysis. Needs to be in radians
      cout << "Loading a slope raster for you. THIS IS IMPORTANT I HOPE YOU CHECKED IT WAS IN RADIANS OTHERWISE YOU WILL GET WRONG RESULTS. SORRY NOT SORRY" << endl;
      LSDRaster this_alpha(load_path + this_string_map["alpha_raster"], "bil");
      alpha_raster = this_alpha;
    }

    // Dealing with raster potential resampling
    if(this_bool_map["resample_slope"])
    {
      cout << "I am resampling the slope raster, because I am calculating Factor of safety for any value!!!" << endl;
      cout << "This is a quick resampling: calculate center pixel of old resolution and apply it to the whole pixel of the new one." << endl;
      cout << "Alternative would be to use GDAL to resample with more complex methods such as nearest neighboors or cubic and directly load it without resampling" << endl;
      LSDRaster relevant_variable_name_3 = alpha_raster.Resample(this_float_map["resample_slope_res"]);
      alpha_raster = relevant_variable_name_3;
      alpha_raster.write_raster(load_path+saving_prefix + "slope_in_radian_resampled", "bil");


    }
    else
    {
      cout << "You chose not to resample the slope raster. Just keep in mind that I will model all the pixels of the raster, it can take a loooot of time then" << endl;
    }

  }

  // Commenting some cout
  // cout << "In the params v2, Iz_over_K_steady: " << Iz_over_K_steady << endl;
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates beta
// Comes from Iverson's eq 27 (it is actually near the top of the first column)
// on page 1902
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::calculate_beta()
{

  beta = cos(alpha)*cos(alpha)-Iz_over_K_steady;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This calculates D_hat
// Comes from equation 26c
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::calculate_D_hat()
{

  D_hat = 4*D_0*(cos(alpha)*cos(alpha));
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This uses the parameters to get a steady state pressure profile
// See the first column in Iverson 2000 page 1902 for explanation
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterParams::calculate_steady_psi()
{
  vector<float> Psi;
  for(int i = 0 ; i < int(Depths.size()); i++ )
  {
    Psi.push_back(beta*(Depths[i]-d));
  }
  return Psi;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// These are a bunch of time manipulating fuctions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterParams::weeks_to_seconds(vector<float> weeks)
{
  vector<float> seconds;
  for (int i = 0; i< int(weeks.size()); i++)
  {
    seconds.push_back(weeks_to_seconds(weeks[i]));
  }
  return seconds;
}

vector<float> LSDPorewaterParams::days_to_seconds(vector<float> days)
{
  vector<float> seconds;
  for (int i = 0; i< int(days.size()); i++)
  {
    seconds.push_back(days_to_seconds(days[i]));
  }
  return seconds;
}

float LSDPorewaterParams::weeks_to_seconds(float weeks)
{
  float seconds = weeks*7.0*24.0*3600.0;
  return seconds;
}

float LSDPorewaterParams::days_to_seconds(float days)
{
  float seconds = days*24.0*3600.0;
  return seconds;
}


vector<float> LSDPorewaterParams::seconds_to_weeks(vector<float> seconds)
{
  vector<float> weeks;
  for (int i = 0; i< int(seconds.size()); i++)
  {
    weeks.push_back(seconds_to_weeks(seconds[i]));
  }
  return weeks;
}

vector<float> LSDPorewaterParams::seconds_to_days(vector<float> seconds)
{
  vector<float> days;
  for (int i = 0; i< int(seconds.size()); i++)
  {
    days.push_back(seconds_to_days(seconds[i]));
  }
  return days;
}

float LSDPorewaterParams::seconds_to_weeks(float seconds)
{
  float weeks = seconds/(7.0*24.0*3600.0);
  return weeks;
}

float LSDPorewaterParams::seconds_to_days(float seconds)
{
  float days = seconds/(24.0*3600.0);
  return days;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a csv file into a data map that contains a map into string vectors
// You will need to convert and clean data later.
// It is inteded to be flexible so you will end up with as many vectors 
// as there were columns in the csv file and these vectors can be accessed by
// referening the key, which is the column header.
// All of the vectors contain strings so you will need to know what data is in
// specific files and convert to those data types at a later stage.  
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<string, vector<string> > LSDPorewaterParams::load_csv_data_into_data_map(string filename)
{
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv data file, but the file" << filename
         << " doesn't exist;  LSDSpatialCSVReader::load_csv_data" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "I have opened the csv file." << endl;
  }

  // Initiate the data map
  map<string, int > temp_vec_vec_key;
  vector< vector<string> > temp_vec_vec;
  map<string, vector<string> > temp_data_map;

  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  // get the headers from the first line
  getline(ifs, line_from_file);

  // reset the string vec
  this_string_vec = empty_string_vec;

  // create a stringstream
  stringstream ss(line_from_file);
  ss.precision(9);

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
  // now check the data map
  int n_headers = int(this_string_vec.size());
  vector<string> header_vector = this_string_vec;
  for (int i = 0; i<n_headers; i++)
  {
    temp_data_map[header_vector[i]] = empty_string_vec;
  }


  // now loop through the rest of the lines, getting the data.
  while( getline(ifs, line_from_file))
  {
    //cout << "Getting line, it is: " << line_from_file << endl;
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
    if ( int(this_string_vec.size()) <= 0)
    {
      cout << "Hey there, I am trying to load your csv data but you seem not to have" << endl;
      cout << "enough columns in your file. I am ignoring a line" << endl;
    }
    else
    {
      int n_cols = int(this_string_vec.size());
      //cout << "N cols is: " << n_cols << endl;
      for (int i = 0; i<n_cols; i++)
      {
        temp_data_map[header_vector[i]].push_back(this_string_vec[i]);
      }
      //cout << "Done with this line." << endl;
    }

  }


  return temp_data_map;
  //cout << "Done reading your file." << endl;
}
//==============================================================================

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Converts a string vec to a float vec
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterParams::convert_string_vec_to_float(vector<string>& string_vec)
{
  float NoData = -9999;
  vector<float> number_vec;
  string this_string;
  string NaN = "NaN";
  string nan = "nan";
  int n_nodes = int(string_vec.size());
  
  // Loop through data looking for NaNs, and if not record the number. 
  for(int i = 0; i< n_nodes; i++)
  {
    if (string_vec[i].find(NaN) != string::npos  || string_vec[i].find(nan) != string::npos)
    {
      cout << "I found a NaN" << endl;
      number_vec.push_back(NoData);
    }
    else
    {
      // Warning!!! I don't check if this is actually a number or not. To be added. 
      number_vec.push_back(atof(string_vec[i].c_str()));
    } 
  }

  return number_vec;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Converts a string vec to an int vec
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDPorewaterParams::convert_string_vec_to_int(vector<string>& string_vec)
{
  int NoData = -9999;
  vector<int> number_vec;
  string this_string;
  string NaN = "NaN";
  string nan = "nan";
  int n_nodes = int(string_vec.size());
  
  // Loop through data looking for NaNs, and if not record the number. 
  for(int i = 0; i< n_nodes; i++)
  {
    if (string_vec[i].find(NaN) != string::npos  || string_vec[i].find(nan) != string::npos)
    {
      cout << "I found a NaN" << endl;
      number_vec.push_back(NoData);
    }
    else
    {
      // Warning!!! I don't check if this is actually a number or not. To be added. 
      number_vec.push_back(atoi(string_vec[i].c_str()));
    } 
  }

  return number_vec;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses a rainfall file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_rainfall_file(string path, string filename, vector<double>& durations, vector<double>& intensities)
{
  string full_fname;
  string PathName = FixPath(path);

  full_fname = PathName+filename;
  
  map<string, vector<string> > csv_data_map = load_csv_data_into_data_map(full_fname);

  // BORIS WORK FROM HERE

  // Alright, let's get the data
  size_t size_of_csv = csv_data_map["timestamp_utc"].size(); // size of time column and of all the csv
  // First, I need to conver tthe data from the files
  vector<double> raw_time_from_file(size_of_csv), raw_intensity_from_file(size_of_csv);

  for(size_t i =0; i<size_of_csv; i++)
  {
    // Getting the duration data and converting it to double
    string str_of_time = csv_data_map["timestamp_utc"][i];
    stringstream ss;
    ss << setprecision(12);
    ss << str_of_time ;
    double double_of_time = atof(ss.str().c_str()); 
    // Getting the rainfall intensity now
    string str_of_prec = csv_data_map["precip_MLP_mm"][i]; double double_of_prec = atof(str_of_prec.c_str()); // WARNING the MLP code site will probably change depending on the site
    // Feeding my vectors
    raw_intensity_from_file.push_back(double_of_prec); raw_time_from_file.push_back(double_of_time);
    cout << setprecision(9);
    cout << str_of_time << "||" << double_of_time << "||" << atof(ss.str().c_str()) << endl;

  }

  // I now have my timing, let's combine it:
  combine_time_to_durations_with_intensities(raw_time_from_file,raw_intensity_from_file,durations,intensities);
  // should be done

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This takes a vector of time and a vector of corresponding intensity and combine it into vector of duration and vector of cumulated intensity for that duration.
// The out vectors have to be empty when feeded
// B.G.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::combine_time_to_durations_with_intensities( vector<double>& raw_time, vector<double>& raw_prec, vector<double>& durations, vector<double>& intensities)
{
  size_t size_of_vecs = raw_time.size();
  double last_prec = raw_prec[0], cumulative_time = 0, cumulative_prec = 0;
  for(size_t i=0;i<size_of_vecs;i++)
  {
    double this_prec = raw_prec[i], this_time = raw_time[i];
    if(this_prec == last_prec)
    {
      if(i!=0)
        cumulative_time = cumulative_time + (this_time - raw_time[i-1]);
      cumulative_prec = cumulative_prec+this_prec;
    }
    else
    {
      durations.push_back(cumulative_time);
      intensities.push_back(cumulative_prec/cumulative_time);
      cumulative_prec =0;
      cumulative_prec = cumulative_prec+this_prec;
      cumulative_time = (this_time - raw_time[i-1]);

    }
    last_prec = this_prec;
  }
  // The last data
  durations.push_back(cumulative_time);
  intensities.push_back(last_prec);
  // Done


}

void LSDPorewaterParams::get_duration_intensity_from_preprocessed_input(vector<float>& duration_s, vector<float>& this_intensity)
{

  // First I need to load the csv file
  map<string, vector<string> > csv_data_map = load_csv_data_into_data_map(rainfall_csv_path + rainfall_csv_name);
  // Simply reading and converting data from it (thanks to python preprocessing)
  duration_s = convert_string_vec_to_float(csv_data_map["duration_s"]);
  this_intensity = convert_string_vec_to_float(csv_data_map["intensity_mm_sec"]);


  // Converting to the right units
  // first taking care of dimensionnalizing the intensity 
  for(size_t i=0; i<duration_s.size(); i++)
  {
    this_intensity[i] = (this_intensity[i]*0.001)/(K_sat);
    if(this_intensity[i]>1) // Cannot be >1
      this_intensity[i]=1;
    else if(this_intensity[i]>1) // Cannot be >1
      this_intensity[i]=0;

  }
  cout << "здесь" << endl;


}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses a rainfall file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_rainfall_file(string path, string filename, vector<float>& intensities)
{
  
  string fname = FixPath(path)+ filename;
  

  
  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;
  float this_rain; 
  vector<float> rain_vec;
  
  vector<string> HeaderInfo = ReadCSVHeader(path, filename);
  
  // now find the rainfall column
  string rain_string = "rainfall_rate";
  string this_string;
  int rain_column = 0; 
  for(int i = 0; i< int(HeaderInfo.size()); i++)
  {
    cout << "Header["<<i<<"]: " << HeaderInfo[i] << endl;
    this_string = HeaderInfo[i];
    if (this_string.compare(rain_string) == 0)
    {
      cout << "I found the rain rate, it is column " << i << endl;
      rain_column = i;
    }
  }
  
  // now we work through the file. 
  // make sure the filename works
  ifstream ifs(fname.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv rainfall file, but the file" << filename
         << "doesn't exist; LINE 357 LSDPorewaterParams" << endl;
    exit(EXIT_FAILURE);
  }
  
  // get the first line  and discard
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
    
    // Now extract the rain rate
    this_rain =  atof(this_string_vec[rain_column].c_str());
    rain_vec.push_back(this_rain);
    
  }
  intensities = rain_vec;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses a rainfall file, specifically a rainfall file derived from MIDAS data

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_MIDAS_rainfall_file(string path, string filename,vector<int>& days, vector<float>& intensities )
{
  
  
  cout << "Path: " << path << endl;
  cout << "filename: " << filename << endl;
  
  string fname = FixPath(path)+filename;
  cout << "Fname is:"<<fname << endl;

  ifstream ifs;
  ifs.open(fname.c_str());
  
  if( ifs.fail() )
  {
    cout << "\nERROR: The parameter \"" << fname
         << "\" doesn't exist. I am not doing anything." << endl;
  }
  else
  {
    // These are 
    vector<int> day_since_1900;
  
    // initiate the string to hold the file
    string line_from_file;
    vector<string> empty_string_vec;
    vector<string> this_string_vec;
    string temp_string;
    float this_rain; 
    int this_date;
    vector<float> rain_vec;
  
    vector<string> HeaderInfo = ReadCSVHeader(path, filename);
  
    // now find the data columns column
    string rain_string = "prcp_amt";
    string this_string;
    int rain_column = 0; 
    string date_string = "days_since_1900";
    int date_column = 0; 
    for(int i = 0; i< int(HeaderInfo.size()); i++)
    {
      cout << "Header["<<i<<"]: " << HeaderInfo[i] << endl;
      this_string = HeaderInfo[i];
      if (this_string.compare(rain_string) == 0)
      {
        cout << "I found the rain rate, it is column " << i << endl;
        rain_column = i;
      }
      if (this_string.compare(date_string) == 0)
      {
        cout << "I found the date, it is column " << i << endl;
        date_column = i;
      }
    }

    // now we work through the file. 
  
    // get the first line  and discard
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
      
      // Now extract the rain rate
      this_rain =  atof(this_string_vec[rain_column].c_str());
      this_date = atoi(this_string_vec[date_column].c_str());
      rain_vec.push_back(this_rain);
      day_since_1900.push_back(this_date);
    }
    intensities = rain_vec;
    days = day_since_1900;
    
    vector<float> new_intensities;
    vector<float> durations;
    parse_MIDAS_duration_intensities(days, intensities, durations);
  
    // bug checking
    for(int i = 0; i<int(durations.size()); i++)
    {
      cout << "duration["<<i<<"]: " << durations[i] << " intensity: "  << intensities[i] << endl;
    }
  
    // now change the units so they are in m/s
    // MIDAS data is mm/day
    // These are divided by K_sat to give the Iz/Kz numbers
    for(int i = 0; i<int(durations.size()); i++)
    {
      durations[i] = durations[i]*3600.0*24.0;
      intensities[i] = (intensities[i]*0.001)/(3600.0*24.0*K_sat);
      cout << "duration["<<i<<"]: " << durations[i] << " intensity: "  << intensities[i] << endl;
      
      // the maximum intensity is 1
      if(intensities[i] > 1)
      {
        intensities[i] = 1;
      }
    }
  
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This parses MIDAS intensity and day data to give multiday durations if
// there are empty data slots or repeated rainfall amounts
// The end goal is a vector with durations and intensities. 
// Some of these can be 0 intensity!
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterParams::parse_MIDAS_duration_intensities(vector<int>& days, vector<float>& intensities, vector<float>& durations)
{
  vector<float> new_intensities;
  vector<float> new_durations;

  // we first do a sweep getting intermediate days
  vector<float> intermediate_intensities;
  intermediate_intensities.push_back(intensities[0]);
  int n_records = int(intensities.size());
  for(int i = 1; i<n_records; i++)
  {
    // If we have incremented by a single day, just push back the data
    if(days[i] == days[i-1]+1)
    {
      intermediate_intensities.push_back(intensities[i]);
    }
    else
    {
      // here we need to push back zeros for the missing days
      int day_dif = days[i]-days[i-1]-1;
      for(int j = 0; j<day_dif; j++)
      {
        intermediate_intensities.push_back(0.0);
      }
      // now push back the next intensity
      intermediate_intensities.push_back(intensities[i]);
    }
  }
  
  // Now go back over the intensities and set the durations and intensities
  int n_int_intensities = int(intermediate_intensities.size());
  int this_consecutive = 0;
  float last_intensity = intermediate_intensities[0];
  for(int i = 1; i<n_int_intensities; i++)
  {
    // add a consecutive day
    this_consecutive++;
    
    // check to see if the intensity has changed
    // if the intensity has not changed, you don't do anything.
    if (intermediate_intensities[i] != last_intensity)
    {
      // the intensity has changed. Add the last one to the record and update 
      // the last changed intensity
      new_intensities.push_back(last_intensity);
      new_durations.push_back(float(this_consecutive));
      last_intensity = intermediate_intensities[i];
      this_consecutive = 0;
    }
  }
  // now you need to deal with the last point. There is always data there
  if(this_consecutive == 0)
  {
    // This means that the last data point was a new intensity, so there is only
    // one day left
    new_intensities.push_back(last_intensity);
    new_durations.push_back(1.0);
  }
  else
  {
    // this means that the last data point was repeated, so we need to use
    // different logic
    new_intensities.push_back(last_intensity);
    new_durations.push_back(float(this_consecutive));
  }
  
  
  // check to see if it worked
  //int n_intensities = int(new_intensities.size());
  //for(int i = 0; i<n_intensities; i++)
  //{
  //  cout << "Duration: " << new_durations[i] << " i: " << new_intensities[i] << endl;
  //}
  
  // now reset the vectors
  intensities = new_intensities;
  durations = new_durations; 
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis prints the parameters to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterParams::print_parameters_to_screen()
{
  cout << "Your parameters are: " << endl;
  cout << "D_0: " << D_0 << endl;  
  cout << "K_sat: " << K_sat << endl; 
  cout << "D_hat: " << D_hat << endl;
  cout << "alpha: " << alpha << endl;
  cout << "d: " << d << endl;
  cout << "beta: " << beta << endl;
  cout << "Iz_over_K_steady: " << Iz_over_K_steady << endl;
  cout << "And the depths are: " << endl;
  for (int i = 0; i< int(Depths.size()); i++)
  {
    cout << Depths[i] << ",";
  }
  cout << endl;
}

#endif