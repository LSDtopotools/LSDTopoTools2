//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexRaster
// Land Surface Dynamics IndexRaster
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing raster data, with a particular focus on topography
//
// The IndexRaster object stores only integer values and is used mostly
//  for storing indices into raster data.
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

/** @file LSDIndexRaster.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 0.0.1
@brief Object to handle integer rasters.

@date 20/08/2012
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDIndexRaster_H
#define LSDIndexRaster_H

#include <string>
#include <vector>
#include <map>
#include "TNT/tnt.h"
#include "LSDShapeTools.hpp"
//#include "LSDRaster.hpp"

using namespace std;
using namespace TNT;

class LSDRaster;        // do not include since it wouldn't compile since there
                        // would be a looped dependency

/// @brief Object to handle integer rasters.
class LSDIndexRaster
{
  public:
  /// @brief  The create function. This is default and throws an error.
  /// @author SMM
  /// @date 01/01/12
  LSDIndexRaster()              { create(); }

  /// @brief Create an LSDIndexRaster from a file.
  /// Uses a filename and file extension
  /// @return LSDIndexRaster
  /// @param filename A String, the file to be loaded.
  /// @param extension A String, the file extension to be loaded.
  /// @author SMM
  /// @date 01/01/12
  LSDIndexRaster(string filename, string extension)  { create(filename, extension); }

  /// @brief Create an LSDIndexRaster from memory.
  /// @return LSDIndexRaster
  /// @param nrows An integer of the number of rows.
  /// @param ncols An integer of the number of columns.
  /// @param xmin A float of the minimum X coordinate.
  /// @param ymin A float of the minimum Y coordinate.
  /// @param cellsize A float of the cellsize.
  /// @param ndv An integer of the no data value.
  /// @param data An Array2D of integers in the shape nrows*ncols,
  ///containing the data to be written.
  /// @author SMM
  /// @date 01/01/12
  LSDIndexRaster(int nrows, int ncols, float xmin, float ymin,
            float cellsize, int ndv, Array2D<int> data)
                  { create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }

  /// @brief Create an LSDIndexRaster from memory.
  /// @return LSDIndexRaster
  /// @param nrows An integer of the number of rows.
  /// @param ncols An integer of the number of columns.
  /// @param xmin A float of the minimum X coordinate.
  /// @param ymin A float of the minimum Y coordinate.
  /// @param cellsize A float of the cellsize.
  /// @param ndv An integer of the no data value.
  /// @param data An Array2D of integers in the shape nrows*ncols,
  /// @param GRS_map a map containing information about the georeferencing
  /// containing the data to be written.
  /// @author SMM
  /// @date 09/09/14
  LSDIndexRaster(int nrows, int ncols, float xmin, float ymin,
            float cellsize, int ndv, Array2D<int> data, map<string,string> GRS_map)
             { create(nrows, ncols, xmin, ymin, cellsize, ndv, data, GRS_map); }

  /// @brief Create an LSDIndexRaster with a constant value.
  /// @return LSDIndexRaster
  /// @param nrows An integer of the number of rows.
  /// @param ncols An integer of the number of columns.
  /// @param xmin A float of the minimum X coordinate.
  /// @param ymin A float of the minimum Y coordinate.
  /// @param cellsize A float of the cellsize.
  /// @param ndv An integer of the no data value.
  /// @param GRS_map a map containing information about the georeferencing
  /// @param ConstValue the value all elements in array have
  /// containing the data to be written.
  /// @author SMM
  /// @date 20/05/16
  LSDIndexRaster(int nrows, int ncols, float xmin, float ymin,
            float cellsize, int ndv, map<string,string> GRS_map, int ConstValue)
             { create(nrows, ncols, xmin, ymin, cellsize, ndv, GRS_map, ConstValue); }

  /// @brief Create an LSDIndexRaster from an LSDRaster object, rounding to nearest int
  /// @return LSDIndexRaster
  /// @param NonIntLSDRaster an LSDRaster object containing flaoting point data
  /// @author MDH
  /// @date 17/02/15
  LSDIndexRaster(LSDRaster& NonIntLSDRaster)   { create(NonIntLSDRaster); }

  /// @brief Create an LSDIndexRaster that is the same size as a raster
  ///  but all values are some constant value
  /// @param NonIntLSDRaster an LSDRaster object
  /// @param CanstValue a value that will be assigned to all data points
  /// @author SMM
  /// @date 19/05/16
  LSDIndexRaster(LSDRaster& ARaster,int ConstValue)   { create(ARaster, ConstValue); }


  // Get functions

  /// @return Number of rows as an integer.
  int get_NRows() const           { return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const           { return NCols; }
  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const           { return XMinimum; }
  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const           { return YMinimum; }
  /// @return Data resolution as an integer.
  float get_DataResolution() const           { return DataResolution; }
  /// @return No Data Value as an integer.
  int get_NoDataValue() const           { return NoDataValue; }
  /// @return Raster values as a 2D Array.
  Array2D<int> get_RasterData() const { return RasterData; }
  /// @return Map of strings containing georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }


  /// Assignment operator.
  LSDIndexRaster& operator=(const LSDIndexRaster& LSDIR);

  /// @brief Read a raster into memory from a file.
  ///
  /// The supported formats are .asc and .flt which are
  /// both exported and imported by arcmap.
  ///
  /// The filename is the string of characters before the '.' in the extension
  /// and the extension is the characters after the '.'.
  ///
  /// If the full filename is my_dem.01.asc then:
  /// filename = "my_dem.01" and extension = "asc".
  ///
  ///
  /// For float files both a data file and a header are read
  /// the header file must have the same filename, before extention, of
  /// the raster data, and the extension must be .hdr.
  /// @author SMM
  /// @date 01/01/12
  void read_raster(string filename, string extension);

  /// @brief Return the list of value but check if it has been initialized
  /// @author BG
  /// @date 19/09/2017
  vector<int> get_list_of_values();

  /// @brief Read a raster from memory to a file.
  ///
  /// The supported formats are .asc and .flt which are
  /// both exported and imported by arcmap.
  ///
  /// The filename is the string of characters before the '.' in the extension
  /// and the extension is the characters after the '.'.
  ///
  /// If the full filename is my_dem.01.asc then:
  /// filename = "my_dem.01" and extension = "asc".
  ///
  /// For float files both a data file and a header are written
  /// the header file must have the same filename, before extention, of
  /// the raster data, and the extension must be .hdr.
  /// @author SMM
  /// @date 01/01/12
  void write_raster(string filename, string extension);

  /// @brief Get the raster data at a specified location.
  /// @param row An integer, the X coordinate of the target cell.
  /// @param column An integer, the Y coordinate of the target cell.
  /// @return The raster value at the position (row, column).
  /// @author SMM
  /// @date 01/01/12
  int get_data_element(int row, int column)  { return RasterData[row][column]; }

  /// @brief Sets the raster data at a specified location.
  /// @param row An integer, the X coordinate of the target cell.
  /// @param column An integer, the Y coordinate of the target cell.
  /// @param value the vaule of the data element at the give row and column
  /// @return The raster value at the position (row, column).
  /// @author SMM
  /// @date 18/03/15
  void set_data_element(int row, int column, int data)
                                   { RasterData[row][column] = data; }

  /// @brief this gets the x and y location of a node at row and column
  /// @param row the row of the node
  /// @param col the column of the node
  /// @param x_loc the x location (Northing) of the node
  /// @param y_loc the y location (Easting) of the node
  /// @author SMM
  /// @date 22/12/2014
  void get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc);

  /// @brief this gets the x and y location of a node at row and column
  /// @param row the row of the node
  /// @param col the column of the node
  /// @param x_loc the x location (Northing) of the node
  /// @param y_loc the y location (Easting) of the node
  /// @author SMM
  /// @date 22/12/2014
  void get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc);

  /// @brief a function to get the lat and long of a node in the raster
  /// @detail Assumes WGS84 ellipsiod
  /// @param row the row of the node
  /// @param col the col of the node
  /// @param lat the latitude of the node (in decimal degrees, replaced by function)
  ///  Note: this is a double, because a float does not have sufficient precision
  ///  relative to a UTM location (which is in metres)
  /// @param long the longitude of the node (in decimal degrees, replaced by function)
  ///  Note: this is a double, because a float does not have sufficient precision
  ///  relative to a UTM location (which is in metres)
  /// @param Converter a converter object (from LSDShapeTools)
  /// @author SMM
  /// @date 24/05/2015
  void get_lat_and_long_locations(int row, int col, double& lat,
                  double& longitude, LSDCoordinateConverterLLandUTM Converter);

  /// @brief This gets the value at a point in UTM coordinates
  /// @param UTME the easting coordinate
  /// @param UTMN the northing coordinate
  /// @return The value at that point
  /// @author SMM
  /// @date 14/03/2017
  int get_value_of_point(float UTME, float UTMN);

  /// @brief this function gets the UTM_zone and a boolean that is true if
  /// the map is in the northern hemisphere
  /// @param UTM_zone the UTM zone. Replaced in function.
  /// @param is_North a boolean that is true if the DEM is in the northern hemisphere.
  ///  replaced in function
  /// @author SMM
  /// @date 22/12/2014
  void get_UTM_information(int& UTM_zone, bool& is_North);

  /// @brief Method to flatten an LSDIndexRaster and place the non NDV values in a csv file.
  /// @detail Each value is placed on its own line, so that it can be read more quickly in python etc.
  ///   It includes the x and y locations so it can be read by GIS software
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @author SMM
  /// @date 29/6/15
  void FlattenToCSV(string FileName);

  /// @brief Method to flatten an LSDIndexRaster and place the non NDV values in a csv file.
  ///
  /// @detail Each value is placed on its own line, so that it can be read more quickly in python etc.
  ///   It includes the lat long coordinates in CSV
  /// @param FileName_prefix The prefix of the file to write, if no path is included it will write to the current directory.
  ///  The csv extension is added automatically.
  /// @author SMM
  /// @date 12/11/16
  void FlattenToWGS84CSV(string FileName);


  /// @brief Checks to see if two rasters have the same dimensions
  /// @detail Does NOT check georeferencing
  /// @param Compare_raster: the raster to compare
  /// @author SMM
  /// @date 04/05/2015
  bool does_raster_have_same_dimensions(LSDRaster& Compare_raster);

  /// @brief Checks to see if two rasters have the same dimensions
  /// @detail Does NOT check georeferencing
  /// @param Compare_raster: the raster to compare
  /// @author SMM
  /// @date 04/05/2015
  bool does_raster_have_same_dimensions(LSDIndexRaster& Compare_raster);

  /// @brief Checks to see if two rasters have the same georeferencing
  /// @param Compare_raster: the raster to compare
  /// @author SMM
  /// @date 02/03/2015
  bool does_raster_have_same_dimensions_and_georeferencing(LSDRaster& Compare_raster);

  /// @brief Checks to see if two rasters have the same georeferencing
  /// @param Compare_raster: the raster to compare
  /// @author SMM
  /// @date 02/03/2015
  bool does_raster_have_same_dimensions_and_georeferencing(LSDIndexRaster& Compare_raster);

  /// @brief This returns a clipped raster that has the same dimensions as the
  ///  smaller raster
  /// @param smaller_raster the raster to which the bigger raster should be
  ///  clipped
  /// @author SMM
  /// @date 20/03/2015
  LSDIndexRaster clip_to_smaller_raster(LSDRaster& smaller_raster);

  /// @brief This returns a clipped raster that has the same dimensions as the
  ///  smaller raster
  /// @param smaller_raster the raster to which the bigger raster should be
  ///  clipped
  /// @author SMM
  /// @date 20/03/2015
  LSDIndexRaster clip_to_smaller_raster(LSDIndexRaster& smaller_raster);


  /// @brief Method which takes a new xmin and ymax value and modifys the GeoReferencingStrings
  /// map_info line to contain these new values.
  ///
  /// @details Intended for use in the rastertrimmer methods and is called from within these methods.
  /// Modifying georeferencing information by hand is messy and should be avoided if
  /// at all possible.
  /// @param NewXmin floating point value of the new minimum x value in the raster.
  /// @param NewYmax floating point value of the new maximum y value in the raster.
  /// @return An updated GeoReferencingStrings object.
  ///
  /// @author SWDG
  /// @date 6/11/14
  map<string, string> Update_GeoReferencingStrings(float NewXmin, float NewYmax);

  /// @brief Method which updates the map info element of the georeferencing strings based on
  /// information within the datamembers of the raster
  ///
  /// @details Intended for use when changing raster dimesions
  ///
  /// @author SMM
  /// @date 6/11/14
  void Update_GeoReferencingStrings();

  /// @brief This method imposes georefereing strings assuming the coordinate
  /// system is UTM
  /// @param zone the UTM zone
  /// @param NorS a string containing characters that start either N (for north)
  /// or S for south. The letter is not case sensitive
  /// @author SMM
  /// @date 6/11/14
  void impose_georeferencing_UTM(int zone, string NorS);

  /// @brief This method looks up the central meridian given a UTM zone
  /// @param UTM_zone the UTM zone
  /// @return central_meridian an integer of the central meridian of this UTM zone
  /// @author SMM
  /// @date 6/11/14
  int Find_UTM_central_meridian(int UTM_zone);

  /// @brief this check to see if a point is within the raster
  /// @param X_coordinate the x location of the point
  /// @param Y_coordinate the y location of the point
  /// @return is_in_raster a boolean telling if the point is in the raster
  /// @author SMM
  /// @date 13/11/2014
  bool check_if_point_is_in_raster(float X_coordinate,float Y_coordinate);

  /// @brief Gets the row and column of a point in the raster
  /// @param X_coordinate the x location of the point
  /// @param Y_coordinate the y location of the point
  /// @param row the row of the point, replaced upon running the routine
  /// @param col the col of the point, replaced upon running the routine
  /// @author SMM
  /// @date 22/01/2016
  void get_row_and_col_of_a_point(float X_coordinate,float Y_coordinate,int& row, int& col);


  /// @brief Calculate the minimum bounding rectangle for an LSDIndexRaster Object and crop out
  /// all the surrounding NoDataValues to reduce the size and load times of output rasters.
  ///
  /// @details Ideal for use with chi analysis tools which output basin and chi m value rasters
  /// which can be predominantly no data. As an example, a 253 Mb file can be reduced to
  /// ~5 Mb with no loss or resampling of data.
  ///
  /// @return A trimmed LSDIndexRaster object.
  /// @author SWDG
  /// @date 22/08/13
  LSDIndexRaster RasterTrimmer();

  /// @brief This function runs the hole finding algorithm but instead of printing
  ///  a raster it returns the points that are in holes. This can be used
  ///  for interpolation.
  /// @param NSteps the number of steps for each cellular automata bot
  /// @param NSweeps the number of times the raster is swept
  /// @param UTME the easting coordinates of the hole points (overwritten)
  /// @param UTMN the northing coordinates of the hole points (overwritten)
  /// @param row_nodes the row numbers of the hole nodes (overwritten)
  /// @param col_nodes the col numbers of the hole nodes (overwritten)
  /// @author SMM
  /// @date 17/03/2017
  void get_points_in_holes_for_interpolation(int NSteps, int NSweeps,
                                          vector<float>& UTME, vector<float>& UTMN,
                                          vector<int>& row_nodes, vector<int>& col_nodes);

  /// @brief This is a brute force method for finding all the nodata regions connected to
  ///  the edge of the raster
  /// @param NSteps the number of steps for each cellular automata bot
  /// @param NSweeps the number of times the raster is swept
  /// @return a raster with 0 for non-visited points and and integer elsewhere
  /// @author SMM
  /// @date 17/3/2017
  LSDIndexRaster find_holes_with_nodata_bots(int NSteps, int NSweeps);

  /// @brief This takes a starting position and releases a random bot that
  ///  moves about in nodata regions, marking its presence
  /// @param Visited and array of numbers for the number of times a pixel is visited
  /// @param startrow the starting row
  /// @param startcol the starting column
  /// @param NSteps the number of steps the bot takes
  /// @author SMM
  /// @date 17/03/2017
  void release_random_bot(Array2D<int>& Visited, int startrow,int startcol, int NSteps);

  /// @brief Make LSDIndexRaster object using a 'template' raster and an Array2D of data.
  /// @param InputData 2DArray of ints to be written to LSDIndexRaster.
  /// @return LSDRaster containing the data passed in.
  /// @author SWDG
  /// @date 02/9/13
  LSDIndexRaster LSDRasterTemplate(Array2D<int> InputData);

  /// @brief Method to resample an LSDIndexRaster to a lower resolution.
  /// @param OutputResolution the resolution in spatial units to be resampled to.
  /// @return An LSDIndexRaster resampled to the OutputResolution.
  /// @author SWDG
  /// @date 17/3/14
  LSDIndexRaster Resample(float OutputResolution);

  /// @brief Method to combine two rasters, ignoring nodata.
  /// @param Network1 The first raster to be combined.
  /// @param Network2 The second raster to be combined.
  /// @return An LSDIndexRaster of the combined inputs.
  /// @author SWDG
  /// @date 17/6/14
  LSDIndexRaster CombineBinaryNetwork(LSDIndexRaster& Network1, LSDIndexRaster& Network2);

  /// @brief Method to merge a floodplain raster with a channel raster.
  ///
  /// @details Creates an output LSDIndexRaster which is coded coded
  /// channel == input channel index, floodplain == 500, NDV == hillslopes.
  /// This allows the preservation of the stream order in addition to the floodplain geometry.
  /// @param FloodPlain an LSDIndexRaster of the floodplains coded with any integer value.
  /// @return An LSDIndexRaster of the merged channels and floodplains.
  /// @author SWDG
  /// @date 05/03/15
  LSDIndexRaster MergeChannelWithFloodplain(LSDIndexRaster FloodPlain);


  /// @brief Method to identify connected components using a two pass method
  ///
  /// @details Takes a binary array, where components parts are identifed by 1
  /// and separates into connected components according to a two-pass algorithm
  /// based on that described in He et al. (2008), "A Run-Based Two-Scan
  /// Labeling Algorithm," Image Processing, IEEE Transactions on , vol.17, no.5,
  /// pp.749,756, doi: 10.1109/TIP.2008.919369
  /// @return an LSDRaster with labelled connected components
  /// @author DTM
  /// @date 13/07/2015
  LSDIndexRaster ConnectedComponents();

  /// @brief Method to filter a binary array according to a connected components threshold
  ///
  /// @author DTM
  /// @date 22/07/2015
  LSDIndexRaster filter_by_connected_components(int connected_components_threshold);


  /// @brief A method to thin a multipixel feature (binary) to a single thread skeleton
  ///
  /// @details Takes a binary array in which features are identified as 1, and
  /// background 0, and thins the features down to a skeleton, using the thinning
  /// algorithm described by Zhang and Suen, 1984, "A fast algorithm for thinning
  /// digital patterns", Communications of the ACM.
  /// @return an LSDIndexRaster with the skeleton
  /// @author DTM
  /// @date 15/07/2015
  LSDIndexRaster thin_to_skeleton();
  void thinningIteration(Array2D<int>& binary, int iter);

  LSDIndexRaster find_end_points();
  void remove_downstream_endpoints(LSDIndexRaster CC, LSDRaster Topo);

  /// @brief Method to convert all values in an LSDIndexRaster to a single value.
  /// @param Value, an integer value that will be assigned to every non NDV cell in the raster.
  /// @param ndv an integer no data value.
  /// @author SWDG
  /// @date 24/07/2015
  LSDIndexRaster ConvertToBinary(int Value, int ndv);

  /// @brief Method to remove patches generated by the connected components analysis that are
  /// smaller than a user defined threshold.
  /// @param minimum_segment_size Size in pixels below which a patch should be removed.
  /// @author SWDG
  /// @date 17/9/15
  LSDIndexRaster RemoveSmallPatches(int minimum_segment_size);

  /// @brief Method to remove small holes in patches from a binary raster.
  /// @param window_radius radius over which to search to remove holes (size of hole)
  /// @author FJC
  /// @date 22/10/15
  LSDIndexRaster remove_holes_in_patches(int window_radius);

  /// @brief Method to remove small holes in patches from a connected components raster. Holes
  /// will only be filled if surrounded by pixels with the same CC value.
  /// @param window_radius radius over which to search to remove holes (size of hole)
  /// @author FJC
  /// @date 20/01/16
  LSDIndexRaster remove_holes_in_patches_connected_components(int window_radius);

  /// @brief Method to fill in checkerboard pattern from a binary raster.
  /// @author FJC
  /// @date 30/10/15
  LSDIndexRaster remove_checkerboard_pattern();

	/// @brief Function to calculate the reliability of floodplain method
  /// @param ActualRaster raster of actual values
  /// @author FJC
  /// @date 26/06/16
  vector<float> AnalysisOfQuality(LSDIndexRaster& ActualRaster);

	/// @brief Function to calculate the percentage area difference between
	/// two binary rasters
  /// @param ActualRaster raster of actual values
	/// @return Percentage difference betwen the area of the two rasters
  /// @author FJC
  /// @date 18/01/17
  float GetAreaDifference(LSDIndexRaster& ActualRaster);

  /// @brief Function to merge data from two LSDIndexRasters WITH SAME EXTENT
  /// together. The data from the raster specified as an argument will be added
  /// (will overwrite the original raster if there is a conflict.)
  /// @param RasterToAdd the raster to merge
	/// @return LSDIndexRaster merged raster
  /// @author FJC
  /// @date 07/04/17
  void MergeIndexRasters(LSDIndexRaster& RasterToAdd);

  /// @brief Function to pad values in an LSDIndexRaster by a certain number of pixels
  /// with values taken from the nearest pixel to be padded.
  /// @param NPixels the number of pixels to pad by
	/// @return LSDIndexRaster padded raster
  /// @author MDH
  /// @date 20/07/17
  void PadRaster(int NPixels);

  /// @brief Function to detect and store all the unique values
  /// into the vector list_unique_values. Detect if this has already been launched to avoid relaunch.
  /// @param No param
	/// @return Nothing, change directly the protected vector attribute
  /// @author BG
  /// @date 17/09/17
  void detect_unique_values();

  /// @brief Function to copy NoData Region from another raster
  /// @param LSDRaster OtherRaster
  /// @return Nothing, change directly the value of the raster
  /// @author BG
  /// @date 20/09/2017 
  void NoData_from_another_raster(LSDRaster& other_raster);

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
  /// list of unique values, for example in case of a lithologic raster
  vector<int> list_unique_values;

  ///A map of strings for holding georeferencing information
  map<string,string> GeoReferencingStrings;

  /// Raster data.
  Array2D<int> RasterData;

  private:
  void create();
  void create(string filename, string extension);
  void create(int ncols, int nrows, float xmin, float ymin,
              float cellsize, int ndv, Array2D<int> data);
  void create(int ncols, int nrows, float xmin, float ymin,
              float cellsize, int ndv, Array2D<int> data,
              map<string,string> GRS_map);
  void create(int nrows, int ncols, float xmin, float ymin,
              float cellsize, int ndv, map<string,string> GRS_map,
              int ConstValue);
  void create(LSDRaster& NonIntLSDRaster);
  void create(LSDRaster& ARaster, int ConstValue);
};

#endif
