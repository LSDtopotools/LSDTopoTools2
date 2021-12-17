//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDSpatialCSVReader.hpp
// Land Surface Dynamics SpatialCSVReader
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for reading csv data. The data needs to have latitude and longitude
//  in WGS84 coordinates.
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2017 Simon M. Mudd 2017
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
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <ctype.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <utility>
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDCosmoData.hpp"
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDBasin.hpp"
#include "LSDSpatialCSVReader.hpp"
#include "LSDRasterInfo.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDSpatialCSVReader_CPP
#define LSDSpatialCSVReader_CPP



// empty create function
void LSDSpatialCSVReader::create()
{
  cout << "Size data map: " << data_map.size() << endl;
  map<string, vector<string> > empty_map;
  cout << "Size empty map: " << empty_map.size() << endl;
  data_map = empty_map;
}

//==============================================================================
// Basic create function
//==============================================================================
void LSDSpatialCSVReader::create(string csv_fname)
{
  //cout << "I am creating something for you" << endl;
  NRows = -9999;
  NCols = -9999;
  XMinimum = -9999;
  YMinimum = -9999;
  DataResolution = -9999;
  NoDataValue = -9999;

  ///A map of strings for holding georeferencing information
  map<string,string> EmptyString;
  GeoReferencingStrings = EmptyString;

  //cout << "Size data map: " << data_map.size() << endl;
  //cout << "Size test map: " << test_map.size() << endl;

  //map<string, vector<string> > empty_map;
  //cout << "Size empty map: " << empty_map.size() << endl;
  //data_map = empty_map;


  load_csv_data(csv_fname);

}



//==============================================================================
// Basic create function
//==============================================================================
void LSDSpatialCSVReader::create(LSDRasterInfo& ThisRasterInfo, string csv_fname)
{
  cout << "I am creating a csv object from a raster info object and a csv name." << endl;
  NRows = ThisRasterInfo.get_NRows();
  NCols = ThisRasterInfo.get_NCols();
  XMinimum = ThisRasterInfo.get_XMinimum();
  YMinimum = ThisRasterInfo.get_YMinimum();
  DataResolution = ThisRasterInfo.get_DataResolution();
  NoDataValue = ThisRasterInfo.get_NoDataValue();
  GeoReferencingStrings = ThisRasterInfo.get_GeoReferencingStrings();

  //cout << "Size data map: " << data_map.size() << endl;
  //cout << "Size test map: " << test_map.size() << endl;

  //map<string, vector<string> > empty_map;
  //cout << "Size empty map: " << empty_map.size() << endl;
  //data_map = empty_map;


  load_csv_data(csv_fname);

}

//==============================================================================
// Basic create function
//==============================================================================
void LSDSpatialCSVReader::create(LSDRaster& ThisRaster, string csv_fname)
{
  cout << "I am creating a csv object from a raster info object and a csv name." << endl;
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();

  load_csv_data(csv_fname);
}



// A create function for getting all the elements to copy or duplicate a csv object
void LSDSpatialCSVReader::create(int nrows, int ncols, float xmin, float ymin,
           float cellsize, float ndv, map<string,string> temp_GRS,
           vector<double>& this_latitude, vector<double>& this_longitude,
           vector<bool>& this_is_point_in_raster, map<string, vector<string> >& this_data_map)
{


  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;

  //cout << "Now for the georeferencing strings" << endl;
  GeoReferencingStrings = temp_GRS;

  latitude = this_latitude;
  longitude = this_longitude;
  is_point_in_raster = this_is_point_in_raster;
  data_map = this_data_map;
  }


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a csv file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::load_csv_data(string filename)
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

  data_map = temp_data_map;

  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  vector<double> temp_longitude;
  vector<double> temp_latitude;

  // get the headers from the first line
  getline(ifs, line_from_file);

  // reset the string vec
  this_string_vec = empty_string_vec;

  // create a stringstream
  stringstream ss(line_from_file);
  ss.precision(14);

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
  int latitude_index = -9999;
  int longitude_index = -9999;
  for (int i = 0; i<n_headers; i++)
  {
    //cout << "This header is: " << this_string_vec[i] << endl;
    if (this_string_vec[i]== "latitude" || this_string_vec[i] == "Latitude" || this_string_vec[i] == "lat" || this_string_vec[i] == "Lat")
    {
      latitude_index = i;
      cout << "The latitude index is: " << latitude_index << endl;

    }
    else if (this_string_vec[i] == "longitude" || this_string_vec[i] == "Longitude" || this_string_vec[i] == "long" || this_string_vec[i] == "Lon")
    {
      longitude_index = i;
      cout << "The longitude index is: " << longitude_index << endl;
    }
    else
    {
      temp_data_map[header_vector[i]] = empty_string_vec;
    }
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
        if (i == latitude_index)
        {
          temp_latitude.push_back( atof(this_string_vec[i].c_str() ) );
        }
        else if (i == longitude_index)
        {
          float this_longitude = atof(this_string_vec[i].c_str() );

          /*
          if (this_longitude < -180)
          {
            this_longitude = 360+this_longitude;
          }
          if (this_longitude > 180)
          {
            this_longitude = this_longitude-360;
          }
          */
          temp_longitude.push_back( this_longitude );

        }
        else
        {
          temp_data_map[header_vector[i]].push_back(this_string_vec[i]);
        }

      }
      //cout << "Done with this line." << endl;
    }

  }



  //cout << "Assigning the vectors." << endl;
  latitude = temp_latitude;
  longitude = temp_longitude;
  data_map = temp_data_map;
  //cout << "Done reading your file." << endl;
}
//==============================================================================




//==============================================================================
// Some functions to check the data members
//==============================================================================
bool LSDSpatialCSVReader::check_if_latitude_and_longitude_exist()
{
  int n_lat, n_long;
  n_lat = int(latitude.size());
  n_long = int(longitude.size());
  bool lat_and_long_exist = false;
  if (n_lat == n_long && n_lat > 0)
  {
    lat_and_long_exist = true;
  }
  return lat_and_long_exist;

}
//==============================================================================


//==============================================================================
// Some functions to check the data members
//==============================================================================
bool LSDSpatialCSVReader::check_if_all_data_columns_same_length()
{

  bool all_data_columns_same_legth = true;
  int n_lat;

  n_lat = int(latitude.size());
  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
  {
    int n_this_column;
    n_this_column = int((it->second).size());

    cout << "The size of this data column is: " <<n_this_column << "\n";

    // if the columns being teh same length is still true, check if the next
    // column is the same length
    if(all_data_columns_same_legth)
    {
      if (n_this_column != n_lat)
      {
        all_data_columns_same_legth = false;
      }

    }
  }

  return all_data_columns_same_legth;

}
//==============================================================================


//==============================================================================
// Some functions to check the data members
//==============================================================================
bool LSDSpatialCSVReader::add_data_column(string column_name, vector<string> column_data)
{

  bool added_column_works = false;
  int n_lat,n_col;

  n_lat = int(latitude.size());

  n_col = int(column_data.size());

  if(n_lat == n_col)
  {
    data_map[column_name] = column_data;
  }
  else
  {
    cout << "The data column is not the same size as the other columns. The addition of this column has failed" << endl;
  }

  return added_column_works;

}
//==============================================================================



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Checks if a column is in the csv
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDSpatialCSVReader::is_column_in_csv(string column_name)
{
  bool is_in_csv = false;
  if ( data_map.find(column_name) == data_map.end() )
  {
    // not found
    cout << "I'm afraid the column "<< column_name << " is not in this dataset" << endl;
  }
  else
  {
    is_in_csv = true;
  }
  return  is_in_csv;
}

string LSDSpatialCSVReader::find_string_in_column_name(string column_name_fragment)
{
  
  string column_name = "NULL";
  string this_key;
  int column_count = 0;
  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
  {
    this_key = it->first;
    if (this_key.find(column_name_fragment) != std::string::npos) 
    {
      column_name = this_key;
      column_count++;
    }
  }

  if(column_count == 0)
  {
    cout << "I didn't find the string fragment. Returning a null column name." << endl;
    cout << "Check if the fragment has a matching case." << endl;
  }
  if(column_count > 1)
  {
    cout << "I found more than one column names with that string fragment!" << endl;
  }
  return  column_name;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns the string vector of data from a given column name
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDSpatialCSVReader::get_data_column(string column_name)
{
  vector<string> data_vector;
  if ( data_map.find(column_name) == data_map.end() )
  {
    // not found
    cout << "I'm afraid the column "<< column_name << " is not in this dataset" << endl;
  }
  else
  {
    data_vector = data_map[column_name];
  }
  return data_vector;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns data that is in the raster for snapping
// It focuses on the ID vector, and converts to Easting-Northing coordinates
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_data_in_raster_for_snapping(string column_name,
                                    vector<float>& UTMEasting,
                                    vector<float>& UTMNorthing,
                                    vector<string>& data_vector)
{
  vector<string> this_data_vector;
  vector<string> thinned_data_vector;
  vector<float> easting;
  vector<float> northing;
  if ( data_map.find(column_name) == data_map.end() )
  {
    // not found
    cout << "I am afraid you tried to access data in the csv file that isn't there." << endl;
    cout << "The missing column name is: " << column_name << endl;
    exit(EXIT_SUCCESS);
  }
  else
  {
    this_data_vector = data_map[column_name];

  }

  // make sure the vector of booleans stating if the point is in the DEM exists
  check_if_points_are_in_raster();

  // convert to UTM
  vector<float> UTME;
  vector<float> UTMN;
  get_x_and_y_from_latlong(UTME,UTMN);

  // now loop through the samples
  int N_samples = int(longitude.size());
  for (int i = 0; i<N_samples; i++)
  {
    if (is_point_in_raster[i])
    {
      easting.push_back(UTME[i]);
      northing.push_back(UTMN[i]);
      thinned_data_vector.push_back(this_data_vector[i]);

      cout << "I am pushing E: " << UTME[i] << " N: " << UTMN[i] << " and data: " << this_data_vector[i] << endl;

    }
  }

  UTMEasting = easting;
  UTMNorthing = northing;
  data_vector = thinned_data_vector;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Converts a data column to a float vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDSpatialCSVReader::data_column_to_float(string column_name)
{
  vector<string> string_vec = get_data_column(column_name);
  vector<float> float_vec;
  int N_data_elements = string_vec.size();
  for(int i = 0; i<N_data_elements; i++)
  {
    float_vec.push_back( atof(string_vec[i].c_str()));
  }
  return float_vec;
}

// Converts a data column to a double vector
vector<int> LSDSpatialCSVReader::data_column_to_int(string column_name)
{
  vector<string> string_vec = get_data_column(column_name);
  vector<int> int_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {
    int_vec.push_back( atoi(string_vec[i].c_str()));
  }
  return int_vec;
}

// Converts a data column to a double vector
vector<double> LSDSpatialCSVReader::data_column_to_double(string column_name)
{
  vector<string> string_vec = get_data_column(column_name);
  vector<double> double_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {
    double_vec.push_back( atof(string_vec[i].c_str()));
  }
  return double_vec;
}

// Adds a float to a data column
void LSDSpatialCSVReader::data_column_add_float(string column_name, float add_value)
{
  //cout << "Adding " << add_value << " to the column " << column_name << endl;
  vector<string> string_vec = get_data_column(column_name);
  float this_value;
  vector<string> new_string_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {

    this_value = atof(string_vec[i].c_str())+add_value;

    //if (i< 10)
    //{
    // cout << "Orig: " << atof(string_vec[i].c_str()) << " final " << this_value << endl;
    //}

    new_string_vec.push_back(to_string(this_value));
  }
  data_map[column_name] = new_string_vec;
}


// Adds a float to a data column
void LSDSpatialCSVReader::data_column_add_float(string column_name, vector<float> add_value)
{
  //cout << "Adding " << add_value << " to the column " << column_name << endl;
  vector<string> string_vec = get_data_column(column_name);
  float this_value;
  vector<string> new_string_vec;
  int N_data_elements = string_vec.size();

  if (N_data_elements != int( add_value.size()))
  {
    cout << "Can't add vectors of different size to csv object" << endl;
    exit(0);
  }

  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {

    this_value = atof(string_vec[i].c_str())+add_value[i];

    //if (i< 10)
    //{
    // cout << "Orig: " << atof(string_vec[i].c_str()) << " final " << this_value << endl;
    //}

    new_string_vec.push_back(to_string(this_value));
  }
  data_map[column_name] = new_string_vec;
}




// Adds a float to a data column
void LSDSpatialCSVReader::data_column_replace(string column_name, vector<float> new_column)
{
  // cout << "Adding " << add_value << " to the column " << column_name << endl;
  vector<string> string_vec = get_data_column(column_name);
  float this_value;
  vector<string> new_string_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }

  if (N_data_elements != int(new_column.size()))
  {
    cout << "You are trying to replace a data column using a column that is a different size. " << endl;
    cout << "I'm afraid I can't do that." << endl;
    cout << "Aborting replace." << endl;
  }
  else
  {
    for(int i = 0; i<N_data_elements; i++)
    {
      new_string_vec.push_back(to_string(new_column[i]));
    }
    data_map[column_name] = new_string_vec;
  }
  

}


// Multiplies values in a data column
void LSDSpatialCSVReader::data_column_multiply_float(string column_name, float multiply_value)
{
  vector<string> string_vec = get_data_column(column_name);
  float this_value;
  vector<string> new_string_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {
    this_value = atof(string_vec[i].c_str())*multiply_value;
    new_string_vec.push_back(to_string(this_value));
  }
  data_map[column_name] = new_string_vec;
}

// This is very specific: it enforces a slope for a particular data colum
void LSDSpatialCSVReader::enforce_slope(string fd_column_name, string elevation_column_name, float minimum_slope)
{
  double min_slope = double(minimum_slope);
  cout << "I am going to enforce a minimum slope." << endl;
  vector<double> flow_distance = data_column_to_double(fd_column_name);
  vector<double> elevation = data_column_to_double(elevation_column_name);

  // the single channel starts from the top, so the last node is the base level and doesn't get modified.
  int N_data_elements = int(flow_distance.size());
  vector<double> new_elevation = elevation;
  float dist, min_elev;

  cout.precision(9);

  //cout << "Minimum slope is: " << min_slope << endl;
  for (int i = N_data_elements-2; i>=0; i--)
  {
    dist = flow_distance[i]-flow_distance[i+1];

    min_elev = dist*min_slope+new_elevation[i+1];
    cout << "dist: " << dist << " z[i+1]: " << new_elevation[i+1] << " z[i]: " << elevation[i] << " min_elev: " << min_elev << endl;
    if (dist < 0 || fabs(dist) >128)
    {
      cout << "There seems to be a big changes in flow distance (greater than a diagonal pixel at 90m resolution" << endl;
      cout << "I am considering this a new channel and resetting the elevation values to this pixel" << endl;
      new_elevation[i] = elevation[i];
    }
    else
    {
      if (elevation[i] < min_elev)
      {
        cout << "Found something where I need to increase slope!" << endl;
        new_elevation[i] = min_elev;
      }
      else
      {
        new_elevation[i] = elevation[i];
      }
    }

  }

  vector<string> new_elev_string;
  for(int i = 0; i< N_data_elements; i++)
  {
    new_elev_string.push_back(to_string(new_elevation[i]));
  }

  data_map[elevation_column_name] = new_elev_string;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_UTM_information(int& UTM_zone, bool& is_North)
{

  // set up strings and iterators
  map<string,string>::iterator iter;

  //check to see if there is already a map info string
  string mi_key = "ENVI_map_info";
  iter = GeoReferencingStrings.find(mi_key);

  if (iter != GeoReferencingStrings.end() )
  {
    string info_str = GeoReferencingStrings[mi_key];

    //cout << "info str is: " << info_str << endl;

    // now parse the string
    vector<string> mapinfo_strings;
    istringstream iss(info_str);
    while( iss.good() )
    {
      string substr;
      getline( iss, substr, ',' );
      mapinfo_strings.push_back( substr );
    }
    UTM_zone = atoi(mapinfo_strings[7].c_str());
    //cout << "Line 1041, UTM zone: " << UTM_zone << endl;
    //cout << "LINE 1042 LSDRaster, N or S: " << mapinfo_strings[7] << endl;

    // find if the zone is in the north
    string n_str = "n";
    string N_str = "N";
    is_North = false;
    size_t found = mapinfo_strings[8].find(N_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    //cout << "is_North is: " << is_North << endl;

  }
  else
  {
    UTM_zone = NoDataValue;
    is_North = false;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets some coordinate system strings for UTM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::set_UTM_information(int UTM_zone, bool is_North)
{
  string cs_key = "coordinate system string";

  string fpart = "{PROJCS[\"WGS_1984_UTM_Zone_";
  string spart = itoa(UTM_zone);

  string tpart;
  if(is_North)
  {
    tpart = "N";
  }
  else
  {
    tpart = "S";
  }
  string fopart = ",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",";

  int start_cm = -183;
  int cm = start_cm+UTM_zone*6;
  string fipart = itoa(cm);
  string sipart = "],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]}";

  string cs_string = fpart+spart+tpart+fopart+fipart+sipart;
  GeoReferencingStrings[cs_key] = cs_string;
  cout << "The coordinate system string is: " << cs_string << endl;

  string mi_key = "ENVI_map_info";
  string mi_p1 = "{UTM, 1, 1, -9999, -9999, -9999, -9999,";
  string mi_p2 = itoa(UTM_zone);

  string mi_p3;
  if(is_North)
  {
    mi_p3 = ", North,WGS-84}";
  }
  else
  {
    mi_p3 = ", South,WGS-84}";
  }

  string mi = mi_p1+mi_p2+mi_p3;
  GeoReferencingStrings[mi_key] = mi;
  cout << "the map info string is:" << mi << endl;

}


//==============================================================================
// This gets the x and y locations from the latitude and longitude
//==============================================================================
void LSDSpatialCSVReader::get_x_and_y_from_latlong(vector<float>& UTME,vector<float>& UTMN)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  int N_samples =  int(latitude.size());

  // set up some temporary vectors
  vector<float> this_UTMN(N_samples,0);
  vector<float> this_UTME(N_samples,0);

  double this_Northing;
  double this_Easting;

  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);


  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    //cout << "Converting point " << i << " to UTM." << endl;
    Converter.LLtoUTM_ForceZone(eId, latitude[i], longitude[i],
                      this_Northing, this_Easting, UTM_zone);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
    //cout << "Easting: " << this_Easting << " and northing: " << this_Northing << endl;
  }

  UTME = this_UTME;
  UTMN = this_UTMN;
}

//==============================================================================
// This gets the x and y locations from specified columns
//==============================================================================
void LSDSpatialCSVReader::get_x_and_y_from_latlong_specify_columns(string lat_column_name,
  string long_column_name, vector<float>& UTME,vector<float>& UTMN)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  // get the data to vectors
  vector<float> lat_data = data_column_to_float(lat_column_name);
  vector<float> long_data = data_column_to_float(long_column_name);

  int N_samples =  int(long_data.size());

  // set up some temporary vectors
  vector<float> this_UTMN(N_samples,0);
  vector<float> this_UTME(N_samples,0);

  double this_Northing;
  double this_Easting;

  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);


  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    //cout << "Converting point " << i << " to UTM." << endl;
    Converter.LLtoUTM_ForceZone(eId, lat_data[i], long_data[i],
                      this_Northing, this_Easting, UTM_zone);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
    //cout << "Easting: " << this_Easting << " and northing: " << this_Northing << endl;
  }

  UTME = this_UTME;
  UTMN = this_UTMN;
}


//==============================================================================
// This gets the latitude and longitude from x and y columns
//==============================================================================
void LSDSpatialCSVReader::get_latlong_from_x_and_y(string X_column_name, string Y_column_name)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  vector<double> new_lat;
  vector<double> new_long;

  int UTM_zone;
  bool is_North;
  int eId = 22;
  get_UTM_information(UTM_zone,is_North);
  cout << "Getting lat and long from UTM Easting and Northing. " << endl;
  cout << "Zone: " << UTM_zone << " and is north? ";
  if (is_North)
  {
    cout <<  "youbetcha!" << endl;
  }
  else
  {
    cout << " no, it is south." << endl;
  }

  vector<float> X_data = data_column_to_float(X_column_name);
  vector<float> Y_data = data_column_to_float(Y_column_name);

  int N_x = int(X_data.size());
  int N_y = int(Y_data.size());
  if (N_x != N_y)
  {
    cout << "Your X and Y columns don't have the same lengths, something has gone wrong." << endl;
    cout << "I am not updating the latitude and longitude" << endl;
  }
  else
  {
    for (int i = 0; i<N_x; i++)
    {
      double thisX = X_data[i];
      double thisY = Y_data[i];

      double Lat;
      double Long;
      Converter.UTMtoLL(eId, thisY, thisX, UTM_zone, is_North,Lat, Long);
      new_lat.push_back(Lat);
      new_long.push_back(Long);

      //if (i == 0)
      //{
      //  cout << "X: " << thisX << " Y: " << thisY << " Lat: " << Lat << " Long: " << Long << endl;
      //}

    }

    latitude = new_lat;
    longitude = new_long;
  }

}


//==============================================================================
// Burns data from a raster to the shapefile
//==============================================================================
void LSDSpatialCSVReader::burn_raster_data_to_csv(LSDRaster& ThisRaster,string column_name)
{
  vector<float> UTME;
  vector<float> UTMN;
  float this_UTME, this_UTMN;
  float this_value;

  vector<string> new_column_data;

  // The csv file needs to have lat-long data
  if (check_if_latitude_and_longitude_exist() == false)
  {
    cout << "You must have lat-long data for burning to work. " << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "Let me get the x and y data." << endl;
    get_x_and_y_from_latlong(UTME,UTMN);
    cout << "Got the x and y" << endl;

    int n_nodes = int(UTME.size());
    for(int i = 0; i<n_nodes; i++)
    {
      stringstream s;
      s.precision(9);
      this_UTME = UTME[i];
      this_UTMN = UTMN[i];

      this_value = ThisRaster.get_value_of_point(this_UTME, this_UTMN);
      //cout << "Node is: " << i << " and value is: " << this_value << endl;
      s << this_value;
      new_column_data.push_back(s.str());
    }
    data_map[column_name] = new_column_data;
  }

}

//==============================================================================
// Burns data from a raster to the shapefile
//==============================================================================
void LSDSpatialCSVReader::burn_raster_data_to_csv(LSDIndexRaster& ThisRaster,string column_name)
{
  vector<float> UTME;
  vector<float> UTMN;
  float this_UTME, this_UTMN;
  int this_value;

  vector<string> new_column_data;

  // The csv file needs to have lat-long data
  if (check_if_latitude_and_longitude_exist() == false)
  {
    cout << "You must have lat-long data for burning to work. " << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    get_x_and_y_from_latlong(UTME,UTMN);
    int n_nodes = int(UTME.size());
    for(int i = 0; i<n_nodes; i++)
    {
      this_UTME = UTME[i];
      this_UTMN = UTMN[i];
      this_value = ThisRaster.get_value_of_point(this_UTME, this_UTMN);
      new_column_data.push_back(itoa(this_value));
    }

    data_map[column_name] = new_column_data;

  }
}


//==============================================================================
// This checks if points are in raster
//==============================================================================
void LSDSpatialCSVReader::check_if_points_are_in_raster()
{
  vector<float> UTME;   // easting
  vector<float> UTMN;   // northing
  vector<bool> temp_is_point_in_raster;

  // get the easting and northing
  get_x_and_y_from_latlong(UTME,UTMN);

  bool is_in_raster;

  int N_samples = int(latitude.size());
  for(int i = 0; i<N_samples; i++)
  {

    is_in_raster = true;

    // Shift origin to that of dataset
    float X_coordinate_shifted_origin = UTME[i] - XMinimum;
    float Y_coordinate_shifted_origin = UTMN[i] - YMinimum;

    // Get row and column of point
    int col_point = int(X_coordinate_shifted_origin/DataResolution);
    int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

    if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1)
    {
      is_in_raster = false;
    }
    temp_is_point_in_raster.push_back(is_in_raster);
  }

  is_point_in_raster = temp_is_point_in_raster;
}

vector<bool> LSDSpatialCSVReader::get_if_points_are_in_raster_vector()
{
  return is_point_in_raster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to get vectors of x and y coordinates, and the node indices of these
// points - this DOES NOT use latitude and longitude, instead it assumes that
// your csv file has columns labelled "X" and "Y". Can be used when you want to
// read in points without converting from latitude and longitude.
// FJC 21/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_nodeindices_from_x_and_y_coords(LSDFlowInfo& FlowInfo, vector<float>& X_coords, vector<float>& Y_coords, vector<int> & NodeIndices)
{
  // get vectors from the columns in the csv file
  string X_coords_name = "X";
  string Y_coords_name = "Y";
  vector<float> X_coords_temp = data_column_to_float(X_coords_name);
  vector<float> Y_coords_temp = data_column_to_float(Y_coords_name);
  vector<int> NodeIndices_temp;

  for (int i = 0; i < int(X_coords_temp.size()); ++i)
  {
    int NodeIndex = FlowInfo.get_node_index_of_coordinate_point(X_coords_temp[i], Y_coords_temp[i]);
    if (NodeIndex != NoDataValue) { NodeIndices_temp.push_back(NodeIndex); }
  }

  //copy to output vectors
  X_coords = X_coords_temp;
  Y_coords = Y_coords_temp;
  NodeIndices = NodeIndices_temp;
}

// get the node indices from lat-long coords in the csv file
// If the nodeindex is nodata it will return nodatavalue
vector<int> LSDSpatialCSVReader::get_nodeindices_from_lat_long(LSDFlowInfo& FlowInfo)
{
  vector<int> NIs;
  vector<float> X_coords, Y_coords;
  get_x_and_y_from_latlong(X_coords,Y_coords);
  for (int i = 0; i < int(X_coords.size()); i++)
  {
    int NodeIndex = FlowInfo.get_node_index_of_coordinate_point(X_coords[i], Y_coords[i]);
    if (NodeIndex != NoDataValue) {
      NIs.push_back(NodeIndex);
    }
  }

  return NIs;

}

void LSDSpatialCSVReader::add_nodeindex_vector_to_data_map_using_lat_long(LSDFlowInfo& FlowInfo)
{
  vector<int> NIs;
  vector<string> NI_strings;
  vector<float> X_coords, Y_coords;
  get_x_and_y_from_latlong(X_coords,Y_coords);
  for (int i = 0; i < int(X_coords.size()); i++)
  {
    int NodeIndex = FlowInfo.get_node_index_of_coordinate_point(X_coords[i], Y_coords[i]);
    NIs.push_back(NodeIndex);
    NI_strings.push_back(itoa(NodeIndex));
  }

  string nistr1 = "nodeindex";
  add_data_column( nistr1, NI_strings);
}



vector<int> LSDSpatialCSVReader::get_nodeindex_vector()
{
  vector<int> ni_vec;

  bool is_nodeindex = false;
  bool is_id = false;
  bool is_node = false;

  //bool ni_exists = false;

  string nistr1 = "nodeindex";
  string nistr2 = "node";
  string nistr3 = "id";
  string ni_column = nistr1;

  // Because we have been sloppy, the nodeindex can appear in csv files as "id", "node", or "nodeindex"
  // we check for all three of these

  cout << "The column headers are: " << endl;
  print_data_map_keys_to_screen();

  cout << "Checking columns to find node index" << endl;
  is_nodeindex = is_column_in_csv(nistr1);
  is_node = is_column_in_csv(nistr2);
  is_id = is_column_in_csv(nistr3);
  cout << "Okay, done checking columns." << endl;

  if ( is_nodeindex == false &&  is_node == false &&  is_id == false)
  {
    cout << "I could not find a nodeindex column. Returning and empty map." << endl;
  }
  else
  {
    // This load of switches basically says that the order of preference if there are
    // more than one liklely columns is nodeindex, node, id.
    if (is_nodeindex)
    {
      ni_column = nistr1;
    }
    else
    {
      if (is_node)
      {
        ni_column = nistr2;
      }
      else
      {
        cout << "I found the code 'id'" << endl;
        ni_column = nistr3;
      }
    }

    // now get the data
    cout << "Grabbing the data from columns " << ni_column << endl;
    ni_vec = data_column_to_int(ni_column);
  }

  return ni_vec;
}

// This uses a nodeindex colum and creates a data map with the column data
map<int,float> LSDSpatialCSVReader::get_nodeindex_map_float(string column_name)
{
  map<int,float> nodeindex_map;

  bool is_nodeindex = false;
  bool is_id = false;
  bool is_node = false;

  //bool ni_exists = false;

  string nistr1 = "nodeindex";
  string nistr2 = "node";
  string nistr3 = "id";
  string ni_column = nistr1;

  // Because we have been sloppy, the nodeindex can appear in csv files as "id", "node", or "nodeindex"
  // we check for all three of these

  cout << "The column headers are: " << endl;
  print_data_map_keys_to_screen();

  cout << "Checking columns to find node index" << endl;
  is_nodeindex = is_column_in_csv(nistr1);
  is_node = is_column_in_csv(nistr2);
  is_id = is_column_in_csv(nistr3);
  cout << "Okay, done checking columns." << endl;

  cout << "The column name is: " << column_name << endl;
  bool is_data_column = is_column_in_csv(column_name);

  if (is_data_column == false)
  {
    cout << "I can't find your data column" << endl;
  }
  else
  {
    cout << "I found your data column." << endl;
    if (is_nodeindex == false && is_node == false && is_id == false)
    {
      cout << "I could not find a nodeindex column. Returning and empty map." << endl;
    }
    else
    {
      // This load of switches basically says that the order of preference if there are
      // more than one liklely columns is nodeindex, node, id.
      if (is_nodeindex)
      {
        ni_column = nistr1;
      }
      else
      {
        if (is_node)
        {
          ni_column = nistr2;
        }
        else
        {
          cout << "I found the code 'id'" << endl;
          ni_column = nistr3;
        }
      }

      // now get the data
      cout << "Grabbing the data from columns " << ni_column << " and " << column_name << endl;
      vector<int> ni_vec = data_column_to_int(ni_column);
      vector<float> data_vec = data_column_to_float(column_name);

      // now make the map
      for (int i = 0; i< int(ni_vec.size()); i++)
      {
        nodeindex_map[ ni_vec[i] ] = data_vec[i];
      }

    }
  }


  return nodeindex_map;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This selects data and then returns a new LSDSpatialCSVobject with only
// the selected data
// Note: this is brute force appraoch: there is probably a faster way to do this!
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDSpatialCSVReader LSDSpatialCSVReader::select_data_to_new_csv_object(string selection_column, vector<string> data_for_selection)
{
  // get the data column
  vector<string> select_column = get_data_column(selection_column);
  int n_nodes = int(select_column.size());
  vector<int> selected_indices;
  string this_item;

  for(int i = 0; i<n_nodes; i++)
  {
    //cout << "Looking for: " << select_column[i] << endl;
    if (std::find(data_for_selection.begin(), data_for_selection.end(), select_column[i]) != data_for_selection.end())
    {
      //cout << "Found it!" << endl;
      selected_indices.push_back(i);
    }
  }

  // now we need to go through the data and remove the rows that don't meet selection
  int n_selected_nodes = int(selected_indices.size());
  map<string, vector<string> > new_data_map;
  vector<string> empty_vec;
  vector<double> new_latitude;
  vector<double> new_longitude;

  // set up the new data map
  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
  {
    //cout << "Key is: " <<it->first << "\n";
    new_data_map[it->first] = empty_vec;
  }

  int this_index;
  for(int i = 0; i<n_selected_nodes; i++)
  {
    this_index = selected_indices[i];

    new_latitude.push_back(latitude[this_index]);
    new_longitude.push_back(longitude[this_index]);

    // now loop through the data map
    for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
    {
      //cout << "Key is: " <<it->first << "\n";
      string element = it->second[this_index];
      new_data_map[it->first].push_back(element);
    }
  }

  // now create the new csv object
  vector<bool> new_in_raster_vec;
  LSDSpatialCSVReader new_csv(NRows,NCols,XMinimum,YMinimum,DataResolution,
                             NoDataValue,GeoReferencingStrings,new_latitude,
                             new_longitude,new_in_raster_vec,new_data_map);
  return new_csv;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints the data column keys to screen
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::print_data_map_keys_to_screen()
{
  cout << "These are the keys: " << endl;
  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
  {
    cout << "Key is: " <<it->first << "\n";
  }
}

//==============================================================================
// This prints the lat and long to screen
//==============================================================================
void LSDSpatialCSVReader::print_lat_long_to_screen()
{
  int N_data = int(latitude.size());
  cout << "latitude,longitude"<< endl;
  cout.precision(9);
  for (int i = 0; i< N_data; i++)
  {
    cout << latitude[i] << "," << longitude[i] << endl;
  }
}

//==============================================================================
// This prints the lat and long to screen
//==============================================================================
void LSDSpatialCSVReader::print_lat_long_to_screen(bool only_print_in_raster)
{
  if (only_print_in_raster)
  {
    check_if_points_are_in_raster();
  }


  int N_data = int(latitude.size());
  cout << "latitude,longitude"<< endl;
  cout.precision(9);
  for (int i = 0; i< N_data; i++)
  {
    if (only_print_in_raster)
    {
      if (is_point_in_raster[i])
      {
        cout << latitude[i] << "," << longitude[i] << endl;
      }

    }
    else
    {
      cout << latitude[i] << "," << longitude[i] << endl;
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Gets the row and column of a point
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_row_and_col_of_a_point(float X_coordinate,float Y_coordinate,int& row, int& col)
{
  int this_row = NoDataValue;
  int this_col = NoDataValue;

  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(ceil(Y_coordinate_shifted_origin/DataResolution)-0.5);

  //cout << "Getting row and col, " << row_point << " " << col_point << endl;

  if(col_point >= 0 && col_point <= NCols-1)
  {
    this_col = col_point;
  }
  if(row_point >= 0 && row_point <= NRows -1)
  {
    this_row = row_point;
  }

  row = this_row;
  col = this_col;
}

void LSDSpatialCSVReader::get_row_and_col_of_a_point(double X_coordinate,double Y_coordinate,int& row, int& col)
{
  int this_row = NoDataValue;
  int this_col = NoDataValue;

  // Shift origin to that of dataset
  double X_coordinate_shifted_origin = X_coordinate - XMinimum;
  double Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(ceil(Y_coordinate_shifted_origin/DataResolution)-0.5);

  //cout << "Getting row and col, " << row_point << " " << col_point << endl;

  if(col_point >= 0 && col_point <= NCols-1)
  {
    this_col = col_point;
  }
  if(row_point >= 0 && row_point <= NRows -1)
  {
    this_row = row_point;
  }

  row = this_row;
  col = this_col;
}



//==============================================================================
// This prints the UTM coordinates to csv for checking
// FJC 03/03/17
//==============================================================================
void LSDSpatialCSVReader::print_UTM_coords_to_csv(vector<float> UTME, vector<float> UTMN, string csv_outname)
{
  ofstream outfile;
  outfile.open(csv_outname.c_str());

  outfile << "PointNo,x,y" << endl;
  outfile.precision(9);
  for (int i = 0; i < int(UTME.size()); i++)
  {
    outfile << i+1 << "," << UTME[i] << "," << UTMN[i] << endl;
  }
  outfile.close();
}

//==============================================================================
// This prints the UTM coordinates to csv for checking
// SMM 30/09/2020
//==============================================================================
void LSDSpatialCSVReader::print_row_and_col_to_csv(string csv_outname)
{
  ofstream outfile;
  outfile.open(csv_outname.c_str());

  int rf,rd,cf,cd;


  vector<float> UTME, UTMN;
  get_x_and_y_from_latlong(UTME,UTMN);



  outfile << "Easting, Norting, row_float,col_float,row_double,col_double" << endl;
  outfile.precision(9);
  for (int i = 0; i < int(latitude.size()); i++)
  {

    get_row_and_col_of_a_point(UTME[i],UTMN[i],rf,cf);
    get_row_and_col_of_a_point(double(UTME[i]),double(UTMN[i]),rd,cd);
    outfile << UTME[i] << "," << UTMN[i] << "," << rf << "," << cf << "," << rd << "," << cd << endl;
  }
  outfile.close();
}




//==============================================================================
// This prints a new csv name
//==============================================================================
void LSDSpatialCSVReader::print_data_to_csv(string csv_outname)
{
  ofstream outfile;

  // make sure the file has the csv extension
  string ext_str = ".csv";
  if (csv_outname.find(ext_str) == std::string::npos) 
  {
    csv_outname = csv_outname + ext_str;
  }


  outfile.open(csv_outname.c_str());

  outfile << "latitude,longitude";
  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
  {
    outfile << "," <<it->first;
  }
  outfile << endl;

  int N_nodes = int(latitude.size());
  for (int i = 0; i < N_nodes; i++)
  {
    outfile.precision(9);

    double this_longitude;
    if (longitude[i] > 180)
    {
      this_longitude = longitude[i]-360.0;
    }
    else if (longitude[i] < -180)
    {
      this_longitude = 360.0+longitude[i];
    }
    else
    {
      this_longitude = longitude[i];
    }

    outfile << latitude[i] << "," << this_longitude;
    for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
    {
      outfile << "," <<it->second[i];
    }
    outfile << endl;
  }

  outfile.close();
}


//==============================================================================
// This prints a new geojson
//==============================================================================
void LSDSpatialCSVReader::print_data_to_geojson(string json_outname)
{

  // the file will be projected in WGS84 so you need lat-long coordinates
  if (check_if_latitude_and_longitude_exist())
  {
    ofstream outfile;
    outfile.precision(9);
    outfile.open(json_outname.c_str());

    outfile << "{" << endl;
    outfile << "\"type\": \"FeatureCollection\"," << endl;
    outfile << "\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\" } }," << endl;
    outfile << "\"features\": [" << endl;

    int n_nodes = int(latitude.size());
    for(int i = 0; i< n_nodes; i++)
    {
      double this_longitude;
      if (longitude[i] > 180)
      {
        this_longitude = longitude[i]-360.0;
      }
      else if (longitude[i] < -180)
      {
        this_longitude = 360.0+longitude[i];
      }
      else
      {
        this_longitude = longitude[i];
      }


      string first_bit = "{ \"type\": \"Feature\", \"properties\": { \"latitude\": ";
      //string second_bit = dtoa(latitude[i])+", \"longitude\": "+ dtoa(this_longitude);
      string second_bit = ", \"longitude\": ";

      string third_bit;
      for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
      {
        third_bit += ", \""+it->first+"\": "+(it->second)[i];
      }
      string fourth_bit = " }, \"geometry\": { \"type\": \"Point\", \"coordinates\": [ ";
      //string fifth_bit = dtoa(this_longitude) +","+ dtoa(latitude[i]) +" ] } },";

      string fifth_bit;
      if (i == (n_nodes-1) )
      {
        fifth_bit = " ] } }";
      }
      else
      {
        fifth_bit = " ] } },";
      }

      outfile << first_bit << latitude[i] << second_bit << this_longitude
              << third_bit+fourth_bit << this_longitude << "," << latitude[i] << fifth_bit
              << endl;
    }
    outfile << "]" << endl;
    outfile << "}" << endl;

    outfile.close();
  }
  else
  {
    cout << "LSDSpatialCSVReader::print_data_to_geojson error." << endl;
    cout << "This dataset does not have lat-long information so I cannot print a geojson" << endl;
  }



}

void LSDSpatialCSVReader::print_data_to_geojson_linestring(string json_outname)
{
  // the file will be projected in WGS84 so you need lat-long coordinates
  if (check_if_latitude_and_longitude_exist())
  {
    ofstream outfile;
    outfile.precision(9);
    outfile.open(json_outname.c_str());

    outfile << "{" << endl;
    outfile << "\"type\": \"FeatureCollection\"," << endl;
    outfile << "\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\" } }," << endl;
    outfile << "\"features\": [" << endl;

    int n_nodes = int(latitude.size());
    for(int i = 0; i< n_nodes-1; i++)
    {
      double this_longitude, next_longitude;
      if (longitude[i] > 180)
      {
        this_longitude = longitude[i]-360.0;
        next_longitude = longitude[i+1]-360.0;
      }
      else if (longitude[i] < -180)
      {
        this_longitude = 360.0+longitude[i];
        next_longitude = 360.0+longitude[i+1];
      }
      else
      {
        this_longitude = longitude[i];
        next_longitude = longitude[i+1];
      }


      string first_bit = "{ \"type\": \"Feature\", \"properties\": { \"latitude\": ";
      //string second_bit = dtoa(latitude[i])+", \"longitude\": "+ dtoa(this_longitude);
      string second_bit = ", \"longitude\": ";

      string third_bit;
      for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
      {
        third_bit += ", \""+it->first+"\": "+(it->second)[i];
      }
      string fourth_bit = " }, \"geometry\": { \"type\": \"Linestring\", \"coordinates\": [ [";
      //string fifth_bit = dtoa(this_longitude) +","+ dtoa(latitude[i]) +" ] } },";

      string fifth_bit;
      if (i == (n_nodes-2) )
      {
        fifth_bit = " ] ] } }";
      }
      else
      {
        fifth_bit = " ] ] } },";
      }

      outfile << first_bit << latitude[i] << second_bit << this_longitude
              << third_bit+fourth_bit << this_longitude << "," << latitude[i] << "], [" << next_longitude << "," << latitude[i+1] << fifth_bit
              << endl;
    }
    outfile << "]" << endl;
    outfile << "}" << endl;

    outfile.close();
  }
  else
  {
    cout << "LSDSpatialCSVReader::print_data_to_geojson error." << endl;
    cout << "This dataset does not have lat-long information so I cannot print a geojson" << endl;
  }

}

//==============================================================================
// This updates XY coordinates with the coordinates centred on a corresponding
// raster pixel.
//==============================================================================
void LSDSpatialCSVReader::centre_XY_by_row_col(LSDRaster& ThisRaster, string X_column_name, string Y_column_name)
{
 // Get the local coordinate as well as the elevations (elevation not really necessary)
  vector<float> UTME = data_column_to_float(X_column_name);
  vector<float> UTMN = data_column_to_float(Y_column_name);
  // vector<float> elev = channel_nodes.data_column_to_float(elevation_column_name);

  // Initialise vectors and variables to store XY data from the raster
  vector<string> new_X_vector;
  vector<string> new_Y_vector;
  float new_X;
  float new_Y;

  cout << "I am comparing your channel node coordinates with the centre coordinates of the corresponding raster pixel and collecting the raster coordinates." << endl;
  int row,col;
  for(int i = 0; i< int(UTME.size()); i++)
  {
    get_row_and_col_of_a_point(UTME[i],UTMN[i],row, col);
    ThisRaster.get_x_and_y_locations(row, col, new_X, new_Y);

    // Some shenanigans to convert the coordinates into a string, letting us add them to the vector and then the csv data map
    UTME[i] = new_X;
    UTMN[i] = new_Y;

    stringstream new_X_string;
    stringstream new_Y_string;
    new_X_string.precision(9);
    new_Y_string.precision(9);
    new_X_string << UTME[i];
    new_Y_string << UTMN[i];
    // Let's add the new X and Y strings to our vectors
    new_X_vector.push_back(new_X_string.str());
    new_Y_vector.push_back(new_Y_string.str());

  }
  cout << "New channel X coordinate vector stores " << int(new_X_vector.size()) << " X coordinates." << endl;

  cout << "I will try to add the new coordinates to your data now. \nWARNING: This will fail if you do not have a latitude column in your data set!" << endl;
  cout << "We must all bow to our lord and saviour lat-long. \nThis will also overwrite your old XY coordinates." << endl;
  add_data_column(X_column_name, new_X_vector);
  add_data_column(Y_column_name, new_Y_vector);

  cout << "Updating lat-long from the new XY coordinates. This will overwrite any old lat-long coordinates." << endl;
  get_latlong_from_x_and_y(X_column_name, Y_column_name);
}

//==============================================================================
// This finds gaps in channel node data and fills the gap by interpolating
// across the gap and adding an extra node with a new value for all columns of
// the csv file. This is not a particularly elegant function but it does what it
// it is supposed to do.
// NOTE: This only adds one node per gap, so you would need multiple iterations
// of this function if you need to add more then one node.
//==============================================================================

void LSDSpatialCSVReader::interpolate_across_gap(LSDRaster& ThisRaster, string X_column_name, string Y_column_name, string fd_column_name)
{
  cout << "Let's find some gaps!" << endl;
  cout << "WARNING: This will only work if your data is SORTED BY FLOW DISTANCE. \nOtherwise your interpolated channel will be a pretzel." << endl;

  vector<float> flow_distance = data_column_to_float(fd_column_name);
  vector<float> UTME = data_column_to_float(X_column_name);
  vector<float> UTMN = data_column_to_float(Y_column_name);

  int row, col;
  float nw_x, nw_y, n_x, n_y, ne_x, ne_y, e_x, e_y, w_x, w_y, sw_x, sw_y, s_x, s_y, se_x, se_y;

  int gap_counter = 0;
  vector<float> guilty_X_vector, guilty_Y_vector;
  float next_x, next_y;

  for(int i = 0; i< int(flow_distance.size()); i++)
  {
    // cout << "Flow distance is " << flow_distance[i] << endl;
    if(i!=int(flow_distance.size()-1))
    {
      // cout << "Next flow distance is " << flow_distance[i+1] << endl;
      next_x = UTME[i+1];
      next_y = UTMN[i+1];
      // cout << "Next coordinates are " << next_x << next_y << endl;
      get_row_and_col_of_a_point(UTME[i],UTMN[i],row, col);
      // cout << "Row is " << row << " and col is " << col << endl;
      ThisRaster.get_x_and_y_locations(row+1, col-1, nw_x, nw_y);
      ThisRaster.get_x_and_y_locations(row+1, col, n_x, n_y);
      ThisRaster.get_x_and_y_locations(row+1, col+1, ne_x, ne_y);
      ThisRaster.get_x_and_y_locations(row, col-1, w_x, w_y);
      ThisRaster.get_x_and_y_locations(row, col+1, e_x, e_y);
      ThisRaster.get_x_and_y_locations(row-1, col-1, sw_x, sw_y);
      ThisRaster.get_x_and_y_locations(row-1, col, s_x, s_y);
      ThisRaster.get_x_and_y_locations(row-1, col+1, se_x, se_y);

      if((next_x == nw_x && next_y == nw_y)
      || (next_x == n_x && next_y == n_y)
      || (next_x == ne_x && next_y == ne_y)
      || (next_x == e_x && next_y == e_y)
      || (next_x == w_x && next_y == w_y)
      || (next_x == sw_x && next_y == sw_y)
      || (next_x == s_x && next_y == s_y)
      || (next_x == se_x && next_y == se_y)
      )
      {
        // cout << "We've got the next node, all good" << endl;
      }
      else
      {
        cout.precision(12);

        cout << "There seems to be a gap." << endl;
        gap_counter = gap_counter +1;

        float diff_flowdist;
        cout << "Flow distance is " << flow_distance[i] << endl;
        cout << "Next flow distance is " << flow_distance[i+1] << endl;

        diff_flowdist = flow_distance[i] - flow_distance[i+1];

        cout << "Flow distance jump is: " << diff_flowdist << endl;

        if(abs(diff_flowdist) > 200)
        {
          cout << "I think we are at the head of a channel, let's ignore this point." << endl;
        }

        else
        {

          // Store coordinates of node before gap
          guilty_X_vector.push_back(UTME[i]);
          guilty_Y_vector.push_back(UTMN[i]);

          for (const auto& kv : get_data_map())
          {
            // cout << "I AM NAME OF COLUMN:" << kv.first << endl;
            vector<string> this_string_vector = kv.second;

            // Converting the vectors to float, adapted from https://stackoverflow.com/questions/35419046/converting-from-vectorstring-to-vectordouble-without-stdstod
            vector<float> this_float_vector;
            // iterate over vector and convert each object to float, then add to new float vector
            for (vector<string>::const_iterator iter = this_string_vector.begin(); iter != this_string_vector.end(); ++iter)
            {
              string const& element = *iter;
              // use a stringstream to get a float value:
              // istringstream is(element);
              // float result;
              // is >> result;
              std::ostringstream out;
              out << std::setprecision(12) << std::stof(element);
              float precise = std::stof(out.str());



              // add the float value to the result vector:
              this_float_vector.push_back(precise);
            }
            // cout << "Downstream of gap, " << kv.first << " is " << this_float_vector[i] << endl;
            // cout << "Upstream of gap, " << kv.first << " is " << this_float_vector[i+1] << endl;
            // cout << "Let's do a super basic interpolation!" << endl;

            float interpolated_value;
            interpolated_value = this_float_vector[i] + ((this_float_vector[i+1] - this_float_vector[i])/2);
            // cout << "Interpolated value of " << kv.first << " is " << interpolated_value << endl;
            stringstream interpolated_string;
            interpolated_string.precision(12);
            interpolated_string << interpolated_value;
            append_to_col(kv.first, to_string(interpolated_value) ); // previous, outdated version: append_to_col(kv.first, interpolated_string.str());
          }
        }
      }
    }
    else
    {
      cout << "That's the last node, I don't need to check for gaps." << endl;
    }
  }

  // We should probably also centre all the coordinates again to make sure that the new ones are good
  cout << "There were " << gap_counter << " gaps, captain. I have filled them." << endl;
  // Let's check that each column is actually increasing in size
        for (const auto& kv : get_data_map())
        {
          vector<string> new_string_vector = kv.second;
          cout << "Column " << kv.first << " stores " << int(new_string_vector.size()) << " values." << endl;
        }
  // let's print the coordinates before each gap to csv to check them
  print_UTM_coords_to_csv(guilty_X_vector, guilty_Y_vector, "bad_channel_points.csv");

  cout << "Getting lat-long from the new XY coordinates. This will overwrite any old lat-long coordinates." << endl;
  get_latlong_from_x_and_y(X_column_name, Y_column_name);

}

//==============================================================================
// This returns the data map
//==============================================================================
map<string,vector<string> >& LSDSpatialCSVReader::get_data_map() {return this->data_map;}

//==============================================================================
// This appends a value to a column
//==============================================================================
void LSDSpatialCSVReader::append_to_col(string colname, string val){this->data_map[colname].push_back(val);}






#endif
