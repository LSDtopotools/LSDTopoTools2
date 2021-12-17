///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDRasterMaker.cpp
/// header for the RasterMaker object
/// The raster maker is a series of simple functions to make some rasters
/// with different properties.
/// The initial use is mainly to make rasters for use in the raster model
/// for uplift and K
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1    01/09/2017
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDSpatialCSVReader.hpp"
using namespace std;
using namespace TNT;

// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif

#ifndef LSDRasterMaker_H
#define LSDRasterMaker_H

///@brief Create model objects to use LSDRaster methods on synthetic landscapes.
class LSDRasterMaker: public LSDRaster
{
  public:

    /// @brief Constructor. Create an LSDRasterMaker from an LSDRaster.
    /// @return LSDRasterMaker
    LSDRasterMaker()
    {
      create();
    }


    /// @brief Constructor. Create an LSDRasterMaker from an LSDRaster.
    /// @return LSDRasterMaker
    /// @param An_LSDRaster LSDRaster object.
    LSDRasterMaker(LSDRaster& An_LSDRaster)
    {
      create(An_LSDRaster);
    }



    /// @brief Constructor. Create an LSDRasterMaker from a file.
    /// Uses a filename and file extension
    /// @return LSDRasterMaker
    /// @param filename A String, the file to be loaded.
    /// @param extension A String, the file extension to be loaded.
    LSDRasterMaker(string filename, string extension)
    {
      create(filename, extension);
    }

    /// @brief Constructor. Create a blank raster nodel
    /// @return LSDRasterModel
    /// @param NCols Height of raster
    /// @param NRows Width of raster
    LSDRasterMaker(int NRows, int NCols)
    {
      create(NRows,NCols);
    }

    /// @brief This just returns the raster model object data as a raster
    /// @return A raster with the data from the LSDRasterModel
    /// @author SMM
    /// @date 01/09/2017
    LSDRaster return_as_raster();

    /// @brief This resizes the LSDRasterModel, resetting some flags in the process,
    /// as well as setting many of the Array2D data members to be empty arrays
    /// The raster data in the end is a random surface (determined by the noise
    /// data member)
    /// This overloaded version also changes the data resolution
    /// @param new_rows the new number of rows
    /// @param new_cols the new number of columns
    /// @param new_resolution the new data resolution
    /// @param new_value the new value of the raster
    /// @author SMM
    /// @date 01/09/2017
    void resize_and_reset( int new_rows, int new_cols, float new_resolution, float new_value );

    /// @brief Gets the row and column of a point in the raster
    /// @param X_coordinate the x location of the point
    /// @param Y_coordinate the y location of the point
    /// @param row the row of the point, replaced upon running the routine
    /// @param col the col of the point, replaced upon running the routine
    /// @author SMM
    /// @date 22/01/2016
    void get_row_and_col_of_a_point(float X_coordinate,float Y_coordinate,int& row, int& col);
    void get_row_and_col_of_a_point(double X_coordinate,double Y_coordinate,int& row, int& col);

    /// @brief Gets the minimum and maximum values from a raster
    /// @return vector with the first element is the minimum and second is the maximum
    /// @author SMM
    /// @date 03/09/2017
    vector<float> minimum_and_maximum_value();

    /// @brief This sets all non nodata pixels to a constant value
    /// @param new_value does what it says on the tin.
    /// @author SMM
    /// @date 18/11/2018
    void set_to_constant_value(float new_value);

    /// @brief This just adds elevation to any non nodata pixel
    /// @param value_to_add The value to add. To subtract, this will be negative. 
    /// @author SMM
    /// @date 09/02/2021
    void add_value(float value_to_add);

    /// @brief Adds a strip of values to the raster
    /// @param start_row_or_col the starting row or column
    /// @param end_row_or_col the ending row or column
    /// @param use_rows if true this creates horizontal strips
    /// @param value the new value of the strip
    /// @author SMM
    /// @date 11/11/2020
    void add_strip(int start_row_or_col, int end_row_or_col, bool horizontal, float value);

    /// @brief This linearly scales the raster to new minimum and maximum values
    /// @param new_minimum does what it says on the tin.
    /// @param new_maxuimum does what it says on the tin
    /// @author SMM
    /// @date 03/09/2017
    void scale_to_new_minimum_and_maximum_value(float new_minimum, float new_maximum);

    /// @brief This fixes a channel, derived from source points data
    ///  onto the model DEM
    /// @param source_points_data an LSDSpatialCSVReader object. It needs lat and long and elevation columns
    /// @param column_name the name of the elevation column
    /// @author SMM
    /// @date 04/03/2020
    void impose_channels(LSDSpatialCSVReader& source_points_data, string column_name);

    /// @brief This fixes a channel, derived from source points data. It also makes sure none of these nodes is 
    ///   adjacent to a nodata pixel
    /// @param source_points_data an LSDSpatialCSVReader object. It needs lat and long and elevation columns
    /// @param slope the slope between the data point and the adjacent pixel
    /// @param column_name the name of the elevation column
    /// @author SMM
    /// @date 21/11/2020
    void impose_channels_with_buffer(LSDSpatialCSVReader& source_points_data, float slope, string column_name);

    /// @brief This fixes a channel, derived from source points data. It also makes sure none of these nodes is 
    ///   adjacent to a nodata pixel. Same as the above function, but uses XY data instead of lat-long
    /// @param source_points_data an LSDSpatialCSVReader object. It needs XY and elevation columns
    /// @param slope the slope between the data point and the adjacent pixel
    /// @param column_name the name of the elevation column
    /// @author ELSG
    /// @date 23/02/2021
    void impose_channels_with_buffer_use_XY(LSDSpatialCSVReader& source_points_data, float slope, string column_name);

    /// @brief This function adds elevation to the nodata pixels around the edge of a DEM
    ///  It is used to crinkle up the edge of of a raster to ensure drainage 
    /// @param slope the slope between the data point and the adjacent pixel
    /// @param source_points_data this gives a channel that defines the outlet
    /// @author SMM
    /// @date 19/03/2021
    void buffer_basin_to_single_outlet(LSDSpatialCSVReader& source_points_data, float slope);


    /// @brief This function adds elevation to the nodata pixels around the edge of a DEM
    ///  it is used on a DEM derived from a basin to ensure internal drainage to 
    ///  a single outlet
    /// @param slope the slope between the data point and the adjacent pixel
    /// @author SMM
    /// @date 10/02/2020
    void buffer_basin_to_single_outlet(float slope);

    /// @brief This smooths the raster. At some point in the future I'll
    ///  add more options but at the moment it just uses 4 neighbours and has
    ///  double weighting on the central pixel. It assumes periodic boundaries
    ///  on the E/W
    /// @param boundary type: this just goes to a default at the moment.
    /// @author SMM
    /// @date 03/09/2017
    void smooth(int boundary_type);

    /// @brief Caps elevations using the initial raster
    ///  WARNING no testing if the raster is the correct shape!
    /// @param InitialRaster The initial raster above which the new surface cannot rise.
    /// @author SMM
    /// @date 31/08/2020
    void cap_elevations(LSDRaster& InitialRaster);

    //void random_horizontal_strips(int minimum_strip_size, int maximum_strip_size, float minimum_value, float maximum_value);

    /// @brief This makes square blobs randomly in the DEM
    /// @param minimum_blob_size minimum blob size
    /// @param maximum_blob_size maximum blob size
    /// @param minimum_value minimum value on the raster
    /// @param maximum_value maximum value on the raster
    /// @parameter n_blobs the number of square blobs
    /// @author SMM
    /// @date 01/09/2017
    void random_square_blobs(int minimum_blob_size, int maximum_blob_size, float minimum_value, float maximum_value, int n_blobs);

    /// @brief This function adds some sine waves together to get patterns
    /// @param x_coefficients Coefficients of the sin waves in the x direction.
    /// @param y_coefficients Coefficients of the sin waves in the y direction.
    /// @author SMM
    /// @date 01/09/2017
    void sine_waves(vector<float> x_coefficients, vector<float> y_coefficients);

    /// @brief This function makes some random values
    /// @param new_minimum does what it says on the tin.
    /// @param new_maxuimum does what it says on the tin
    /// @author SMM
    /// @date 30/07/2019
    void random_values(float minimum_value, float maximum_value);

    /// @brief This returns a clipped raster that has the same dimensions as the
    ///  smaller raster
    /// @param smaller_raster the raster to which the bigger raster should be
    ///  clipped
    /// @author SMM
    /// @date 20/03/2015
    LSDRaster clip_to_smaller_raster(LSDRaster& smaller_raster);

    /// @brief This returns a clipped raster that has the same dimensions as the
    ///  smaller raster
    /// @param smaller_raster the raster to which the bigger raster should be
    ///  clipped
    /// @author SMM
    /// @date 20/03/2015
    LSDRaster clip_to_smaller_raster(LSDIndexRaster& smaller_raster);

    ///@brief This function returns the raster data as text file
  ///@return text file with raster data
  ///@author FJC
  ///@date 30/09/16
    void write_RasterData_to_text_file(string filename);

    ///@brief This function returns the raster data as a vector
  ///@return vector<float> with raster data
  ///@author FJC
  ///@date 06/11/15
  vector<float> get_RasterData_vector();

  /// @brief This function returns a vector with the X adn Y minimum and max
  ///   values
  /// @return XYMinMax a vector with four elements
  ///   XYMinMax[0] = XMinimum
  ///   XYMinMax[1] = YMinimum
  ///   XYMinMax[2] = XMaximum
  ///   XYMinMax[3] = XMaximum
  /// @author SMM
  /// @date 3/7/2015
  vector<float> get_XY_MinMax();


  /// @brief This function creates a raster where the north half of the raster
  /// values are increased by a certain number of times
  /// @param increase_amt Raster values in north half will be multiplied by this.
  /// @author FJC
  /// @date 28/06/2018
  void increase_north_half_raster_values(int increase_amt);

  /// @brief This function creates a raster where the south half of the raster
  /// values are increased by a certain number of times
  /// @param increase_amt Raster values in south half will be multiplied by this.
  /// @author FJC
  /// @date 17/07/18
  void increase_south_half_raster_values(int increase_amt);
  void increase_south_quarter_raster_values(int increase_amt);
  void increase_west_half_raster_values(int increase_amt);

  /// @brief This function creates a raster of tilted values (e.g. for making a
  /// tilted uplift field)
  /// @param angle tilt angle
  /// @param tilt_boundary the boundary which will stay fixed (can be N, S, E, or W)
  /// @author FJC
  /// @date 17/07/18
  void tilted_block(float angle, string tilt_boundary);

  protected:



  private:
    /// @brief Make a 100x100 raster
    void create();

    /// @brief Load a raster
    void create(string filename, string extension);

    /// @brief Make a rastermaker from another raster
    void create(LSDRaster& An_LSDRaster);

    /// @brief create a raster from NRows and NCols information
    void create(int NRows, int NCols);

};

#endif
