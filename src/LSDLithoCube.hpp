///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDLithoCube.hpp
/// This object has a stack of LSDRasters that have lithologic information
/// currently this is read from a csv supplied by NAGRA
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1    25/01/2019
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <string>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;

// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif

#ifndef LSDLithoCube_H
#define LSDLithoCube_H

///@brief Create model objects to use LSDRaster methods on synthetic landscapes.
class LSDLithoCube
{
  public:

    /// @brief Constructor. 
    /// @return LSDLithoCube
    LSDLithoCube()
    {
      create();
    }


    /// @brief Constructor. Create an LSDLithoCube from an LSDRaster.
    /// @return LSDLithoCube
    /// @param An_LSDRaster LSDRaster object.
    LSDLithoCube(LSDRaster& An_LSDRaster)
    {
      create(An_LSDRaster);
    }



    /// @brief Constructor. Create an LSDLithoCube from a file.
    /// Uses a filename and file extension
    /// @return LSDLithoCube
    /// @param filename A String, the file to be loaded.
    /// @param extension A String, the file extension to be loaded.
    LSDLithoCube(string filename)
    {
      create(filename);
    }

    /// @brief This takes a litho layer and writes the layer to a raster
    /// @param layer The layer number you want
    /// @param fname the filename (without extension but with path)
    /// @param f_ext the raster extension (almost always bil)
    /// @author SMM
    /// @date 30/01/2020
    void write_litho_raster_from_layer(int layer, string fname, string f_ext);

    /// @brief This takes an elevation raster and writes the layer to a raster
    /// @param ElevationRaster A raster of elevations
    /// @param forbidden_codes The lithocodes that cannot be kept in the raster (e.g., and air code)
    /// @param fname the filename (without extension but with path)
    /// @param f_ext the raster extension (almost always bil)
    /// @author SMM
    /// @date 30/01/2020
    void write_litho_raster_from_elevations(LSDRaster& ElevationRaster, list<int> forbidden_codes, string fname, string f_ext);

    /// @brief Impose UTM georeferencing strings
    /// @param zone The UTM zone
    /// @param NorS a string "N" or "S" for north or south
    /// @author SMM
    /// @date 30/01/2020
    void impose_georeferencing_UTM(int zone, string NorS);

    /// @brief For use with UTM conversion/ Does what it says on the tin
    /// @param zone The UTM zone
    /// @author SMM
    /// @date 30/01/2020
    int Find_UTM_central_meridian(int UTM_zone);

    /// @brief This ingests the VO file (created in autocad)
    /// @author SMM
    /// @param filename The full filename (including path and extension)
    /// @date 25/01/2020
    void ingest_vo_file(string filename);

    /// @brief This ingests the data
    /// @author SMM
    /// @param filename The full filename (including path and extension)
    /// @date 25/01/2020
    void ingest_litho_data(string filename);

    /// @brief This ingests the stratigraphic information
    /// @author SMM
    /// @param filename The full filename (including path and extension)
    /// @date 25/01/2020
    void ingest_stratigraphic_info(string filename);


    /// @brief This gets the row and column of a particular point
    /// @param x_location the x location of a particular point in local coordinates
    /// @param y_location the y location of a particular point in local coordinates
    /// @param row This gets replaced with the correct row
    /// @param col this gets replaced with the correct column
    /// @author SMM
    /// @date 30/01/2020
    void get_row_and_col_of_a_point(float X_coordinate,float Y_coordinate,int& row, int& col);

    /// @brief Checks to see if the point is in the raster
    /// @param x_location the x location of a particular point in local coordinates
    /// @param y_location the y location of a particular point in local coordinates
    /// @return true if in raster
    /// @author SMM
    /// @date 30/01/2020
    bool check_if_point_is_in_raster(float X_coordinate,float Y_coordinate);


    /// @brief Gets the layer number from the elevation
    /// @param elevation duh
    /// @return the layer number
    /// @author SMM
    /// @date 30/01/2020
    int get_layer_from_elevation(float elevation);

    /// @brief Gets the lithocode from a location. 
    /// @detail Includes logic to drill through forbidden lithocode values until
    ///  it finds a valid value
    /// @param LLayer the starting layer
    /// @param LRow the row in the lithocube
    /// @param LCol the col in the lithocube
    /// @param forbidden_codes the list of codes that are not allowed to be returned
    ///  for example if you want to remove "air" or quaternary pixels in the lithocube
    /// @return the lithocode
    /// @author SMM
    /// @date 31/08/2020
    int get_lithocode_from_location(int LLayer, int LRow, int LCol, list<int> forbidden_codes);

    /// @brief Takes an elevation raster and returns an erodibility raster
    /// @param ElevationRaster A raster of elevations
    /// @return a raster of erodbility (K) values
    /// @author ELSG
    /// @date 30/01/2020
    LSDRaster get_K_raster_from_raster(LSDRaster& ElevationRaster);

    /// @brief Reads in a csv file (lookup table for stratigrapy code and K value)
    /// @param lookup_filename
    /// @author ELSG
    /// @date 30/01/2020
    void read_strati_to_K_csv(string lookup_filename);


    /// @brief Reads in a csv file (lookup table for stratigrapy code and Sc value)
    /// @param lookup_filename
    /// @author SMM
    /// @date 31/08/2020
    void read_strati_to_Sc_csv(string lookup_filename);

    /// @brief Prints the K and Sc code to screen
    /// @author SMM
    /// @date 31/08/2020
    void print_K_and_Sc_maps_to_screen();


    /// @brief Fills a map with stratigraphy code (key) and corresponding K (value)
    /// @author ELSG
    /// @date 30/01/2020
    void hardcode_fill_strati_to_K_map();

    /// @brief Fills a map with stratigraphy code (key) and corresponding Sc (value)
    /// @author ELSG
    /// @date 24/02/2020
    void hardcode_fill_strati_to_Sc_map();

    /// @brief This gets a rater from a specific layer 
    /// @param layer the layer number
    /// @author SMM
    /// @date 30/01/2020
    LSDIndexRaster get_raster_from_layer(int layer);

    /// @brief This takes an elevation raster and gets a stratigraphic layer raster
    /// @param ElevationRaster A raster of elevations
    /// @param forbidden_codes A list of forbidden lithocube codes (e.g., for removing air from codes)
    /// @return StratData The raster containing the stratigraphic information
    /// @author SMM
    /// @date 30/01/2020
    LSDIndexRaster get_lithology_codes_from_raster(LSDRaster& ElevationRaster, list<int> forbidden_codes );

    /// @brief This takes an elevation raster and writes a K raster
    /// @param ElevationRaster A raster of elevations
    /// @param fname the filename (without extension but with path)
    /// @param f_ext the raster extension (almost always bil)
    /// @author ELSG
    /// @date 30/01/2020
    void write_K_raster_from_elevations(LSDRaster& ElevationRaster, string fname, string f_ext);

    /// @brief Translates XY coordinates to row, col
    /// @param X_coord (float): the x coordinate
    /// @param Y_coord (float): the y coordinate
    /// @param row (int): will get the value of the row
    /// @param col (int): will get the value of the col
    /// @param snap_to_closest_point_on_raster: if true and if the points is outside the raster, it snaps it t the closest boundary
    /// @author BG
    /// @date Added to LSDLithoCube 30/01/2020, originally from LSDDEM_xtensor
    void XY_to_rowcol_careful(float X_coord, float Y_coord, int& row, int& col, bool snap_to_closest_point_on_raster);

    
    /// @brief Takes an index raster and creates a raster with corresponding erodibility values
    /// @param index_raster An index raster
    /// @param default_K A default K value for replacing nodata
    /// @author BG
    /// @date 30/01/2020
    LSDRaster index_to_K(LSDIndexRaster& index_raster, float default_K);

     /// @brief Takes an index raster and creates a raster with corresponding critical slope values
    /// @param index_raster An index raster
    /// @param default_Sc A default Sc value for replacing nodata
    /// @author ELSG
    /// @date 24/02/2020
    LSDRaster index_to_Sc(LSDIndexRaster& index_raster, float default_Sc);




  protected:

    ///Number of rows.
    int NRows;
    ///Number of columns.
    int NCols;
    /// Number of layers
    int NLayers;
    ///Minimum X coordinate.
    float XMinimum;
    ///Minimum Y coordinate.
    float YMinimum;
    ///Minimin Z coordinate
    float ZMinimum;

    ///Data resolution.
    float DataResolution;

    ///X spacing
    float XSpacing;
    ///Y spacing;
    float YSpacing;
    ///Z Spacing
    float ZSpacing;

    bool loaded_vo_file;

    ///No data value.
    int NoDataValue;


    ///A map of strings for holding georeferencing information
    map<string,string> GeoReferencingStrings;

    /// The actual data, in LSDRaster format
    vector< Array2D<int> > LithoLayers;

    /// A vector containing the bottom elevations of the layers
    vector<float> bottom_elevations;

    /// A vector containting 
    vector<float> layer_thicknesses;

    /// A map containing the erodibility value for each stratigraphy code
    map<int,float> StratiToK;

    /// A map containing the critical slope value for each stratigraphy code
    map<int,float> StratiToSc;





  private:
    /// @brief Make a 100x100 raster
    void create();

    /// @brief Load a raster
    void create(string filename);

    /// @brief Make a rastermaker from another raster
    void create(LSDRaster& An_LSDRaster);


};

#endif
