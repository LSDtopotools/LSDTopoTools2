//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRaster
// Land Surface Dynamics Raster Info
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for examining the georeferencing and existence of rasters
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


#ifndef LSDRasterInfo_H
#define LSDRasterInfo_H

#include <map>
#include <string>
using namespace std;

// declare classes. Implementation is included in cpp file
class LSDRaster;
class LSDTindexRaster;

///@brief Object that stores georeferencing information. This information is
/// also stored with the raster, it is seperated here mainly to compare 
/// different datasets
class LSDRasterInfo
{
  public:

    /// an empty create function
    LSDRasterInfo()              { create(); }
    
    /// This function reads the file
    LSDRasterInfo(string filename, string extension)
                               { create(filename, extension); }
                               
    /// This gets the information from an existing raster
    LSDRasterInfo(LSDRaster& Raster)
                               { create(Raster); }

    /// This gets the information from an existing index raster
    LSDRasterInfo(LSDIndexRaster& IRaster)
                               { create(IRaster); }

    /// The equality operator
    bool operator==(LSDRasterInfo& LSDRI);

    /// The inequality operator
    bool operator!=(LSDRasterInfo& LSDRI);

    /// @brief this function gets the UTM_zone and a boolean that is true if
    /// the map is in the northern hemisphere
    /// @param UTM_zone the UTM zone. Replaced in function.
    /// @param is_North a boolean that is true if the DEM is in the northern hemisphere.
    ///  replaced in function
    /// @author SMM
    /// @date 22/12/2014
    void get_UTM_information(int& UTM_zone, bool& is_North);

    /// @brief this check to see if a point is within the raster
    /// @param X_coordinate the x location of the point
    /// @param Y_coordinate the y location of the point
    /// @return is_in_raster a boolean telling if the point is in the raster
    /// @author SMM
    /// @date 13/11/2014
    bool check_if_point_is_in_raster(float X_coordinate,float Y_coordinate);

    /// @brief Prints the raster information to screen
    /// @author SMM
    /// @date 19/11/2019
    void print_raster_information();

    // Get functions
    /// @return Number of rows as an integer.
    int get_NRows() const        { return NRows; }
    /// @return Number of columns as an integer.
    int get_NCols() const        { return NCols; }
    /// @return Minimum X coordinate as an integer.
    float get_XMinimum() const        { return XMinimum; }
    /// @return Minimum Y coordinate as an integer.
    float get_YMinimum() const        { return YMinimum; }
    /// @return Data resolution as an integer.
    float get_DataResolution() const        { return DataResolution; }
    /// @return No Data Value as an integer.
    int get_NoDataValue() const        { return NoDataValue; }
    /// @return map containing the georeferencing strings
    map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }

  protected:

    /// @brief Read a header 
    /// The supported formats are .asc, .flt and .bil which are
    /// both exported and imported by arcmap.
    /// The filename is the string of characters before the '.' in the extension
    /// and the extension is the characters after the '.'.
    /// If the full filename is my_dem.01.asc then:
    /// filename = "my_dem.01" and extension = "asc".
    /// For float files both a data file and a header are read
    /// the header file must have the same filename, before extention, of
    /// the raster data, and the extension must be .hdr.
    /// @param filename the prefix of the file
    /// @param extension this is either "asc", "flt", or "bil"
    /// @author SMM
    /// @date 01/01/12
    void read_header(string filename, string extension);



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


  private:
    void create();
    void create(LSDRaster& Raster);
    void create(LSDIndexRaster& IRaster);
    void create(string filename, string extension);

};

#endif
