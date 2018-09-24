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
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "LSDRaster.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
using namespace std;

#ifndef LSDSpatialCSVReader_HPP
#define LSDSpatialCSVReader_HPP

class LSDSpatialCSVReader
{
  public:


    /// @brief Create an LSDSpatialCSVReader from a raster and csv filenname
    /// @param csv_fname The name of the csv file including extension  and path
    /// @author SMM
    /// @date 22/02/2017
    LSDSpatialCSVReader(string csv_fname)  { create(csv_fname); }

    /// @brief Create an LSDSpatialCSVReader from a raster and csv filenname
    /// @param ThisRaster An LSDRaster object
    /// @param csv_fname The name of the csv file including extension  and path
    /// @author SMM
    /// @date 16/02/2017
    LSDSpatialCSVReader(LSDRaster& ThisRaster, string csv_fname)  { create(ThisRaster,csv_fname); }

    /// @brief Create an LSDSpatialCSVReader from a raster.
    /// @param ThisRaster An LSDIndexRaster object
    /// @param csv_fname The name of the csv file including extension  and path
    /// @author SMM
    /// @date 16/02/2017
    LSDSpatialCSVReader(LSDRasterInfo& ThisRaster, string csv_fname)  { create(ThisRaster,csv_fname); }

    /// @brief Create an LSDSpatialCSVReader from all the data elements
    /// @param nrows An integer of the number of rows.
    /// @param ncols An integer of the number of columns.
    /// @param xmin A float of the minimum X coordinate.
    /// @param ymin A float of the minimum Y coordinate.
    /// @param cellsize A float of the cellsize.
    /// @param ndv An integer of the no data value.
    /// @param data An Array2D of floats in the shape nrows*ncols,
    /// @param temp_GRS a map of strings containing georeferencing information. Used
    /// mainly with ENVI format files
    /// @param this_latitude the latitudes in WGS84
    /// @param this_longitude the longitudes in WGS84
    /// @param this_is_point_in_raster bool vec for testing if point in raster domain
    /// @param this_data_map the data elements
    /// @author SMM
    /// @date 15/03/2017
   LSDSpatialCSVReader(int nrows, int ncols, float xmin, float ymin,
           float cellsize, float ndv, map<string,string> temp_GRS,
           vector<double>& this_latitude, vector<double>& this_longitude,
           vector<bool>& this_is_point_in_raster, map<string, vector<string> >& this_data_map)
             { create(nrows, ncols, xmin, ymin, cellsize, ndv, temp_GRS,
                    this_latitude, this_longitude,this_is_point_in_raster, this_data_map); }


    /// @brief This loads a csv file, grabbing the latitude and longitude,
    ///  and putting the rest of the data into data maps
    /// @param filename The name of the csv file including path and extension
    /// @author SMM
    /// @date 16/02/2017
    void load_csv_data(string filename);

    /// @brief This tests to see if there are latitude and longitude vectors
    /// @return lat_and_long_exist This is true if the lat and long vectors exist
    ///  are the same length.
    /// @author SMM
    /// @date 14/03/2017
    bool check_if_latitude_and_longitude_exist();

    /// @brief This tests to see if  all the data columns are the same length
    ///  as the latitude and longitude data columns
    /// @return all_data_columns_exist This is true if the data columns are the
    ///  same length as the lat and long vectyors
    /// @author SMM
    /// @date 14/03/2017
    bool check_if_all_data_columns_same_length();


    /// @brief this function sets a UTM_ coordinate system string
    /// the map is in the northern hemisphere
    /// @param UTM_zone the UTM zone. Replaced in function.
    /// @param is_North a boolean that is true if the DEM is in the northern hemisphere.
    ///  replaced in function
    /// @author SMM
    /// @date 22/02/2017
    void set_UTM_information(int UTM_zone, bool is_North);

    /// @brief this function gets the UTM_zone and a boolean that is true if
    /// the map is in the northern hemisphere
    /// @param UTM_zone the UTM zone. Replaced in function.
    /// @param is_North a boolean that is true if the DEM is in the northern hemisphere.
    ///  replaced in function
    /// @author SMM
    /// @date 22/12/2014
    void get_UTM_information(int& UTM_zone, bool& is_North);

    /// @brief This takes latitude and longitude (in WGS 84) and converts to vectors
    ///  of easting and northing in UTM
    /// @param UTME The easting coordinate (is overwritten)
    /// @param UTMN The northing coordinate (is overwritten)
    /// @author SMM
    /// @date 17/02/2017
    void get_x_and_y_from_latlong(vector<float>& UTME,vector<float>& UTMN);

    /// @brief This takes latitude and longitude (in WGS 84) and converts to vectors
    ///  of easting and northing in UTM. User must specify the column for lat and long (this can be
    /// used for csv files where the column is not headed simply "latitude" or "longitude".
    /// @param lat_column_name The header of the latitude column (string)
    /// @param long_column_name The header of the longitude column (string)
    /// @param UTME The easting coordinate (is overwritten)
    /// @param UTMN The northing coordinate (is overwritten)
    /// @author FJC
    /// @date 01/08/17
    void get_x_and_y_from_latlong_specify_columns(string lat_column_name,
      string long_column_name, vector<float>& UTME,vector<float>& UTMN);

    /// @brief This takes the X and Y columns in a csv, assumes they are UTM, and converts to
    ///  latitude and loingitude in WGS84
    /// @param X_column_name Name of the column in the csv with X coordinate (easting)
    /// @param Y_column_name Name of the column in the csv with Y coordinate (northing)
    /// @author SMM
    /// @date 22/02/2017
    void get_latlong_from_x_and_y(string X_column_name, string Y_column_name);

    /// @brief Takes raster data and adds it to the data column
    /// @param TheRaster an LSDRaster
    /// @author SMM
    /// @date 14/03/2017
    void burn_raster_data_to_csv(LSDRaster& ThisRaster,string column_name);

    /// @brief Takes raster data and adds it to the data column
    /// @param TheRaster an LSDIndexRaster
    /// @author SMM
    /// @date 14/03/2017
    void burn_raster_data_to_csv(LSDIndexRaster& ThisRaster,string column_name);

    /// @brief This gets UTM coordinates and a data vector (usually IDs)
    ///  for snapping to a channel network
    /// @param column_name a string that holds the column name
    /// @param UTMEasting The easting coordinate (is overwritten)
    /// @param UTMNorthing The northing coordinate (is overwritten)
    /// @param data_vector the data in a strong vector from the column. Usually sample ID for snapping
    /// @author SMM
    /// @date 20/02/2017
    void get_data_in_raster_for_snapping(string column_name,
                                    vector<float>& UTMEasting,
                                    vector<float>& UTMNorthing,
                                    vector<string>& data_vector);




    /// @brief This gets a data column from the csv file
    /// @param column_name a string that holds the column name
    /// @return a vector of strings: this holds the data.
    /// @author SMM
    /// @date 17/02/2017
    vector<string> get_data_column(string column_name);

    /// @brief This gets a data column from the csv file, and converts it to a
    ///   float vector
    /// @param column_name a string that holds the column name
    /// @return a vector of floats: this holds the data.
    /// @author SMM
    /// @date 17/02/2017
    vector<float> data_column_to_float(string column_name);

    /// @brief This gets a data column from the csv file, and converts it to an
    ///   int vector
    /// @param column_name a string that holds the column name
    /// @return a vector of ints: this holds the data.
    /// @author SMM
    /// @date 17/02/2017
    vector<int> data_column_to_int(string column_name);


    /// @brief this check to see if a point is within the raster
    /// @param X_coordinate the x location of the point
    /// @param Y_coordinate the y location of the point
    /// @return is_in_raster a boolean telling if the point is in the raster
    /// @author SMM
    /// @date 13/11/2014
    void check_if_points_are_in_raster();

    /// @brief This function gets vectors of x and y coordinates and node indices
    /// from these points
    /// @details This DOES NOT use latitude and longitude, instead
    /// it assumes that your csv file has columns labelled "X" and "Y".
    /// Can be used when you want to read in points without converting
    /// from latitude and longitude.
    /// @param FlowInfo LSDFLowInfo object
    /// @param X_coords vector to write X coords to
    /// @param Y_coords vector to write Y coords to
    /// @param NodeIndices vector to write node indices
    /// @author FJC
    /// @date 21/02/17
    void get_nodeindices_from_x_and_y_coords(LSDFlowInfo& FlowInfo, vector<float>& X_coords, vector<float>& Y_coords, vector<int>& NodeIndices);

    /// @brief This selects specified data and crease a new csv object with just that data
    /// @param selection_column The name of the column from which the data will be selected
    /// @param data_for_selection a vector holding the values of the data column that
    ///  will be selected. Rows without these values will be removed
    /// @author SMM
    /// @date 15/03/2017
    LSDSpatialCSVReader select_data_to_new_csv_object(string selection_column, vector<string> data_for_selection);

    /// @brief this prints keys of the data map to screen
    /// @author SMM
    /// @date 20/02/2017
    void print_data_map_keys_to_screen();

    /// @brief this prints the latitude and longitude to screen
    /// @author SMM
    /// @date 17/02/2017
    void print_lat_long_to_screen();

    /// @brief this prints the latitude and longitude to screen
    /// @param only_print_in_raster a bool when true only prints points in raster
    /// @author SMM
    /// @date 20/02/2017
    void print_lat_long_to_screen(bool only_print_in_raster);

    /// @brief print the UTM coords to a csv file for checking
    /// @param UTME eastings
    /// @param UTMN northings
    /// @author FJC
    /// @date 03/03/17
    void print_UTM_coords_to_csv(vector<float> UTME, vector<float> UTMN, string csv_outname);

    /// @brief print the data to a csv. Used after updating data
    /// @param csv_outname the name of the new file
    /// @author SMM
    /// @date 13/03/17
    void print_data_to_csv(string csv_outname);


    /// @brief print the data to a geojson. Used after updating data
    /// @param json_outname the name of the new file
    /// @author SMM
    /// @date 14/03/17
    void print_data_to_geojson(string json_outname);

    /// Gets the various data members
    vector<double> get_latitude() const {return latitude;}
    vector<double> get_longitude() const {return longitude;}


  protected:

    ///Number of rows.
    int NRows;
    ///Number of columns.
    int NCols;
    ///Minimum X coordinate.
    float XMinimum;
    ///Minimum Y coordinate.
    float YMinimum;

    ///Data resolution.
    float DataResolution;
    ///No data value.
    int NoDataValue;

    ///A map of strings for holding georeferencing information
    map<string,string> GeoReferencingStrings;

    /// A vector of the latitude (in WGS84)
    vector<double> latitude;

    /// A vector of the longitude (in WGS84)
    vector<double> longitude;

    /// A vector of bools telling if the points are in the raster
    vector<bool> is_point_in_raster;

    /// The map
    map<string, vector<string> > data_map;


  private:

    void create();

    void create(string);

    void create(LSDRaster&,string);

    void create(LSDRasterInfo&,string);

    void create(int nrows, int ncols, float xmin, float ymin,
           float cellsize, float ndv, map<string,string> temp_GRS,
           vector<double>& this_latitude, vector<double>& this_longitude,
           vector<bool>& this_is_point_in_raster, map<string, vector<string> >& this_data_map);


};

#endif
