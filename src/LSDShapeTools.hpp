//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDShapeTools
// Land Surface Dynamics Shapefile tools
//
// A collection of routines for maipulating the binary ESRI shapefile format
// for use within the Edinburgh Land Surface Dynamics group topographic toolbox
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

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

using namespace std;
#include <vector>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <vector>
#include <utility>
#include <fstream>
#include <cmath>

// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif


#ifndef ShapeTools_H
#define ShapeTools_H

// An unsigned char can store 1 Byte (8bits) of data
typedef unsigned char BYTE;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to test the Byte order of the system.
// Returns a boolean value where true is little endian.
//
// SWDG 11/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
bool SystemEndiannessTest();


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Structure to store X and Y point data.
//
// SWDG 12/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
struct PointData
{

  vector<double> X;
  vector<double> Y;

};

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to swap the byte order of a word in memory. Used if the system's byte order
// does not match the data in the shapefile.
//
// SWDG 11/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void ByteSwap(int length, void * ByteData);


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to load an ESRI ShapeFile.
//
// Only works for X,Y point shapefiles at present and it's behavious is totally undefined
// if you pass in any other type of file.
//
// In future this will be rebuilt into a full class that can support shapefiles of
// different types.
//
// Built in part from:
// http://www.dreamincode.net/forums/topic/170054-understanding-and-reading-binary-files-in-c/
//
// SWDG 13/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
PointData LoadShapefile(string Filename);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to load an ESRI polyline Shapefile.
//
// Only works for polyline shapefiles at present and it's behaviour is totally undefined
// if you pass in any other type of file.
//
// In future this will be rebuilt into a full class that can support shapefiles of
// different types.
//
// Returns a vector of points. So that each item in the vector represents a single polyline.
//
// Built in part from:
// http://www.dreamincode.net/forums/topic/170054-understanding-and-reading-binary-files-in-c/
//
// SWDG 17/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<PointData> LoadPolyline(string Filename);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to convert an IndexChannelTree to a PointData object.
//
// Returns a vector of points.
//
// multistem_option -> 0 = mainstem only, 1 = all tributaries, 2 specify tributary number (DAV 09/04/2015)
// DTM 11/07/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
PointData LoadChannelTree(string Filename, int multistem_option = 0, int trib_number = 0);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to convert vectors of coordinates to a PointData object.
//
// Returns the point data
// FJC 17/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

PointData get_point_data_from_coordinates(vector<float>& X_coordinates, vector<float>& Y_coordinates);
PointData get_point_data_from_coordinates(vector<double>& X_coordinates, vector<double>& Y_coordinates);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to get the size of the binary file being loaded.
//
// Taken from http://www.dreamincode.net/forums/topic/170054-understanding-and-reading-binary-files-in-c/
//
// SWDG 10/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
long getFileSize(FILE *file);




/// @brief Ellipsoid class for converting coordinates between UTM and lat long
/// @detail This is the old version
class LSDEllipsoid
{
  public:

  /// @detail the default constructor
  LSDEllipsoid()    {};

  /// @detail assigment constructor for the ellipsiod class
  /// @param id a reference into the ellipsoid
  /// @param name the name of the ellipsoid
  /// @param radius the radius of the equator in km
  /// @param fr This is the inverse flattening see https://en.wikipedia.org/wiki/Earth_ellipsoid
  ///   see also https://en.wikipedia.org/wiki/Geodetic_datum
  ///   the eccentricity is defined by flattening. 
  LSDEllipsoid(int id, char* name, double radius, double fr)
      { Name=name;  EquatorialRadius=radius;  eccSquared=2/fr-1/(fr*fr);}

  /// name of the ellipsoid
  char* Name;

  /// equatorial radius in km
  double EquatorialRadius;

  /// square of the eccentricity
  double eccSquared;
};

/// @brief Ellipsoid class for converting coordinates between UTM and lat long
/// @author SMM
/// @date 21/01/2019
class LSDReferenceEllipsoid
{
  public:

    /// @brief the default constructor
    LSDReferenceEllipsoid()    { create();} 

    /// @brief assigment constructor for the ellipsiod class
    /// @param id a reference into the ellipsoid
    /// @param name the name of the ellipsoid
    /// @param radius the radius of the equator in km
    /// @param fr This is the inverse flattening see https://en.wikipedia.org/wiki/Earth_ellipsoid
    ///   see also https://en.wikipedia.org/wiki/Geodetic_datum
    ///   the eccentricity is defined by flattening. 
    LSDReferenceEllipsoid(string ellipsoid_name)
        { create(ellipsoid_name);}

    /// @brief This is intended to check if the ellipsoid is correct but it just returns true 
    ///  since checking is done during object creation
    /// @author SMM
    /// @date 21/01/2019
    bool check_ellipsoid();

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Getter functions
    string get_Name() { return Name; }
    double get_EquatorialRadius() { return EquatorialRadius; }
    double get_f() { return f; }
    double get_f_inverse() { return f_inverse; }
    double get_eccSquared() { return eccSquared; }
    double get_ecc() { return ecc; }


    /// @brief This prints the parameters of the ellipsoid to screen. 
    /// @author SMM
    /// @date 17/07/2019
    void print_ellipsoid_parameters_to_screen();

    /// @brief Sets the ellipsoid parameters
    /// @author SMM
    /// @date 17/07/2019
    void set_ellipsoid_parameters(string Ellipsoid_Name);


  protected:
    /// name of the ellipsoid
    string Name;

    /// equatorial radius in km
    double EquatorialRadius;

    /// Flattening
    double f;

    /// Inverse flattening
    double f_inverse;

    /// square of the eccentricity
    double eccSquared;

    /// the eccentricity
    double ecc;


  private:
    /// The create function. This one assumes WGS84.
    void create();

    /// Create function that uses a name to set the parameters
    void create (string Ellipsoid_Name);
};


/// @brief Datum class for converting coordinates between UTM and lat long
class LSDDatum
{
  public:
    LSDDatum(){};
    LSDDatum(int id, char* name, int eid, double dx, double dy, double dz)
      { Name=name;  eId=eid;  dX=dx;  dY=dy;  dZ=dz;}

  /// name of the datum
  char* Name;

  /// the ellipsoid id
  int   eId;

  double dX;
  double dY;
  double dZ;
};

/// @brief A class for storing projection information
/// We implement a minamalistic interface: it is not meant as a replacement for GDAL
/// Projections and transformations of rasters should happen in GDAL
/// These tools are simply from getting points from projected coordinates
/// to lat-long in EPSG:4326 (WGS84)
/// Initially the only option was UTM but we have expanded to 
/// a hacky version of the british national grid
/// and reasonable implementations of lambert conical conformal and 
/// albers equal area.
/// The latter projections require a bunch of parameters that are stored in a map
/// @author SMM
/// @date 19/01/2019
class LSDProjectionInfo
{
  public:
    /// @brief Default constructor function. Assumes a default projection with parameters
    ///  taken from Snyder 1987 page 296
    /// @author SMM
    /// @date 19/01/2019
    LSDProjectionInfo()  { create(); }

    /// @brief The constructor that takes and EPSG code
    /// @detail The number of available EPSG codes are very limited
    /// @param EPSG_code: The EPSG code of the desired projection.
    ///  For the code to work someone needs to hard code the parameters in the cpp file!!! 
    /// @author SMM
    /// @date 19/01/2019
    LSDProjectionInfo( int EPSG_code)  { create(EPSG_code);}

    /// @brief This create function just assigns all the variables
    /// @detail This is used with the OGCWKT reader object to quickly
    ///  create a map of Projections within it
    /// @param tEPSG_code The epsg code
    /// @param tEPSG_type The type of projection. Various options are avilable
    /// @param tIntProjParameter a map of integer project parameters. 
    /// @param tProjParameter a map of double project parameters. 
    /// @param tRefEllipsoid The reference ellipsoid. 
    /// @author SMM
    /// @date 29/01/2019
    LSDProjectionInfo(int tEPSG_code, string tEPSG_type, 
                      map<string,int> tIntProjParameters, 
                      map<string,double> tProjParameters, 
                      LSDReferenceEllipsoid tRefEllipsoid)
       { create(tEPSG_code, tEPSG_type, 
                tIntProjParameters,tProjParameters, 
                tRefEllipsoid); }

    /// @brief This is intended to check if the projection is correct but it just returns true 
    ///  since checking is done during object creation
    /// @author SMM
    /// @date 21/01/2019
    bool check_projection();  

    /// @detail Prints the projection parameters to screen
    /// @author SMM
    /// @date 19/01/2019
    void print_projection_parameters();

    /// @detail Sets the projection. 
    /// @param EPSG_code: The EPSG code of the desired projection.
    ///  For the code to work someone needs to hard code the parameters in the cpp file!!! 
    /// @author SMM
    /// @date 19/01/2019    
    void set_projection(int EPSG_code);

    /// Get the EPSG code
    int get_EPSG() { return EPSG; }

    /// Get the parameters for the projection
    map<string,double> get_ProjParameters() { return ProjParameters;}

    /// Get the integer parameters for the projection
    map<string,int> get_IntProjParameters() { return IntProjParameters;}

    /// Get the EPSG type
    string get_EPSG_type() { return EPSG_type; }

    /// Get the referene ellipsoid
    LSDReferenceEllipsoid get_ReferenceEllipsoid() { return ReferenceEllipsoid; }

  protected:

    /// The EPSG code
    int EPSG; 

    /// This is the type of projection
    /// IMPORTANT it needs to be derived from the Well KNown Text (WKT) inputs for coordinate systems
    /// For example:
    ///  Transverse_Mercator (e.g., see http://spatialreference.org/ref/epsg/32624/ogcwkt/) (The British National Grid is included here)
    ///  Lambert_Conformal_Conic_2SP  (e.g., see http://spatialreference.org/ref/esri/north-america-lambert-conformal-conic/ogcwkt/)
    ///  Albers_Conic_Equal_Area (e.g., see http://spatialreference.org/ref/esri/102001/ogcwkt/)
    string EPSG_type;

    /// Values for the parameters of the projection
    map<string,double> ProjParameters;

    /// These are integer values
    map<string,int> IntProjParameters;

    LSDReferenceEllipsoid ReferenceEllipsoid;
 
    
  
  private:
    /// @brief Default create function. Assumes a default projection with parameters
    ///  taken from Snyder 1987 page 296
    /// @author SMM
    /// @date 19/01/2019
    void create();

    /// @brief The constructor that takes and EPSG code
    /// @detail The number of available EPSG codes are very limited
    /// @param EPSG_code: The EPSG code of the desired projection.
    ///  For the code to work someone needs to hard code the parameters in the cpp file!!! 
    /// @author SMM
    /// @date 19/01/2019
    void create(int EPSG_code);

    /// @brief The constructor that builds a projection with all the necessary information. 
    /// @author SMM
    /// @date 28/01/2019
    void create(int tEPSG_code, string tEPSG_type, 
                      map<string,int> tIntProjParameters, 
                      map<string,double> tProjParameters, 
                      LSDReferenceEllipsoid tRefEllipsoid);

};

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// A class for converting datums and coordinates
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
/// @brief This is the old version that can only hold UTM to WGS84
class LSDCoordinateConverterLLandUTM
{
  public:
    // default constructor. This sets up the data elements.
    LSDCoordinateConverterLLandUTM()     { create(); }

    /// @brief converts LatLong to UTM coords
    /// 3/22/95: by ChuckGantz chuck.gantz@globalstar.com, from USGS Bulletin 1532.
    /// @param eID the ellipsoid ID. Options are:
    ///  0 = "Airy1830";
    ///  1 = "AiryModified";
    ///  2 = "AustralianNational";
    ///  3 = "Bessel1841Namibia";
    ///  4 = "Bessel1841";
    ///  5 = "Clarke1866";
    ///  6 = "Clarke1880";
    ///  7 = "EverestIndia1830";
    ///  8 = "EverestSabahSarawak";
    ///  9 = "EverestIndia1956";
    ///  10 = "EverestMalaysia1969";
    ///  11 = "EverestMalay_Sing";
    ///  12 = "EverestPakistan";
    ///  13 = "Fischer1960Modified";
    ///  14 = "Helmert1906";
    ///  15 = "Hough1960";
    ///  16 = "Indonesian1974";
    ///  17 = "International1924";
    ///  18 = "Krassovsky1940";
    ///  19 = "GRS80";
    ///  20 = "SouthAmerican1969";
    ///  21 = "WGS72";
    ///  22 = "WGS84";
    /// @param Lat the latitude in decimal degrees
    /// @param Long the longitude in decimal degrees
    /// @param Northing in metres. This argument is replaced by the function
    /// @param Easting in metres. This argument is replaced by the function
    /// @param Zone the UTM zone. This argument is replaced by the function
    /// @author SMM, modified from Chuck Gantz
    /// @date 07/12/2014
    void LLtoUTM(int eId, double Lat, double Long,
             double& Northing, double& Easting, int& Zone);

    /// @brief converts LatLong to UTM coords, but forces the data to a specific zone
    /// 3/22/95: by ChuckGantz chuck.gantz@globalstar.com, from USGS Bulletin 1532.
    /// Updated by Simon M Mudd
    /// @param eID the ellipsoid ID. Options are:
    ///  0 = "Airy1830";
    ///  1 = "AiryModified";
    ///  2 = "AustralianNational";
    ///  3 = "Bessel1841Namibia";
    ///  4 = "Bessel1841";
    ///  5 = "Clarke1866";
    ///  6 = "Clarke1880";
    ///  7 = "EverestIndia1830";
    ///  8 = "EverestSabahSarawak";
    ///  9 = "EverestIndia1956";
    ///  10 = "EverestMalaysia1969";
    ///  11 = "EverestMalay_Sing";
    ///  12 = "EverestPakistan";
    ///  13 = "Fischer1960Modified";
    ///  14 = "Helmert1906";
    ///  15 = "Hough1960";
    ///  16 = "Indonesian1974";
    ///  17 = "International1924";
    ///  18 = "Krassovsky1940";
    ///  19 = "GRS80";
    ///  20 = "SouthAmerican1969";
    ///  21 = "WGS72";
    ///  22 = "WGS84";
    /// @param Lat the latitude in decimal degrees
    /// @param Long the longitude in decimal degrees
    /// @param Northing in metres. This argument is replaced by the function
    /// @param Easting in metres. This argument is replaced by the function
    /// @param Zone the UTM zone. This argument is replaced by the function
    /// @author SMM, modified from Chuck Gantz
    /// @date 11/04/2016
    void LLtoUTM_ForceZone(int eId, double Lat, double Long,
             double& Northing, double& Easting, int Zone);


    /// @brief converts LatLong to UTM coords
    /// 3/22/95: by ChuckGantz chuck.gantz@globalstar.com, from USGS Bulletin 1532.
    /// @param eID the ellipsoid ID. Options are:
    ///  0 = "Airy1830";
    ///  1 = "AiryModified";
    ///  2 = "AustralianNational";
    ///  3 = "Bessel1841Namibia";
    ///  4 = "Bessel1841";
    ///  5 = "Clarke1866";
    ///  6 = "Clarke1880";
    ///  7 = "EverestIndia1830";
    ///  8 = "EverestSabahSarawak";
    ///  9 = "EverestIndia1956";
    ///  10 = "EverestMalaysia1969";
    ///  11 = "EverestMalay_Sing";
    ///  12 = "EverestPakistan";
    ///  13 = "Fischer1960Modified";
    ///  14 = "Helmert1906";
    ///  15 = "Hough1960";
    ///  16 = "Indonesian1974";
    ///  17 = "International1924";
    ///  18 = "Krassovsky1940";
    ///  19 = "GRS80";
    ///  20 = "SouthAmerican1969";
    ///  21 = "WGS72";
    ///  22 = "WGS84";
    /// @param Northing in metres.
    /// @param Easting in metres.
    /// @param Zone the UTM zone.
    /// @param isNorth is a boolean that states if the map is in the northern hemisphere
    /// @param Lat the latitude in decimal degrees.
    ///  This argument is replaced by the function
    /// @param Long the longitude in decimal degrees
    ///  This argument is replaced by the function
    /// @author SMM, modified from Chuck Gantz
    /// @date 07/12/2014
    void UTMtoLL(int eId, double Northing, double Easting, int Zone, bool isNorth,
             double& Lat, double& Long);


    /// @brief converts British national grid to WGS84 lat-long
    void BNGtoLL(double Northing, double Easting, double& Lat, double& Long);

    /// @brief converts LatLongHt in datum dIn, to LatLongHt in datum dTo;
    /// @detail 2002dec: by Eugene Reimer, from PeterDana equations.
    ///   Lat and Long params are in degrees;
    /// North latitudes and East longitudes are positive;  Height is in meters;
    /// ==This approach to Datum-conversion is a waste of time;
    /// to get acceptable accuracy a large table is needed -- see NADCON, NTv2...
    void DatumConvert(int dIn, double LatIn, double LongIn, double HtIn,
                  int dTo,  double& LatTo, double& LongTo, double& HtTo);

  protected:

    /// @brief a vector holding the ellipsoids
    vector<LSDEllipsoid> Ellipsoids;

    /// @brief a vectro holding the datums
    vector<LSDDatum> Datums;

    double RADIANS_PER_DEGREE;
    double DEGREES_PER_RADIAN;

    /** Useful constants **/
    double TWOPI;
    double HALFPI;

    // Grid granularity for rounding UTM coordinates to generate MapXY.
    double grid_size;    // 100 km grid

    // WGS84 Parameters
    double WGS84_A;                     // major axis
    double WGS84_B;                     // minor axis
    double WGS84_F;                     // ellipsoid flattening
    double WGS84_E;                     // first eccentricity
    double WGS84_EP;                    // second eccentricity

    // UTM Parameters
    double UTM_K0;              // scale factor
    double UTM_FE;              // false easting
    double UTM_FN_N;            // false northing, northern hemisphere
    double UTM_FN_S;            // false northing, southern hemisphere
    double UTM_E2;              // e^2
    double UTM_E4;              // e^4
    double UTM_E6;              // e^6
    double UTM_EP2;             // e'^2

  private:

    /// @brief This create function sets up the data membeers that hold the
    ///  ellipsoid and datum data
    void create();

};





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This will hold an object for reading a coordinate system from a OGC WKT file
// Will need to do this later. Thinking about it makes my head hurt. 
class LSDOGCWKTCRSReader
{
  public:
    /// @brief the default constructor. Does nothing
    LSDOGCWKTCRSReader()  {create();}

    /// @brief This creates the object based on a filename that includes the path
    /// @detail This only works if the OGC WKT entries have no line returns in them. That
    ///   is, each projection only occupies one line. 
    /// @param filename the name of the OGC WKT text
    /// @author SMM
    /// @date 29/01/2019
    LSDOGCWKTCRSReader(string filename)  {create(filename);}

    /// @brief This creates the object based on a filename and a path
    /// @detail This only works if the OGC WKT entries have no line returns in them. That
    ///   is, each projection only occupies one line. 
    /// @param path The full path to the OGC WKT file
    /// @param filename the name of the OGC WKT text without path
    /// @author SMM
    /// @date 29/01/2019
    LSDOGCWKTCRSReader(string path,string filename)  {create(path,filename);}

    /// @brief This reads a file with OGC WKT entries on lines
    /// @param The filename with full path
    /// @author SMM
    /// @date 28/01/2019
    void read_OGCWKT_file(string filename);

    /// @brief This creates the data entries for UTM zones
    /// @author SMM
    /// @date 29/01/2019
    void populate_UTM_data();

    /// @brief This is a wrapper function that parses all the relevant information from an OGC WKT file
    /// @detail See https://www.opengeospatial.org/standards/wkt-crs
    /// @param line A string that contains a single OGC WKT entry. Read from a file
    /// @author SMM
    /// @date 29/01/2019
    void parse_OGCWKT_line(string line);

    /// @brief This parses the projection from an OGC WKT entry. We only allow certain options. See source code for those options. 
    /// @detail See https://www.opengeospatial.org/standards/wkt-crs
    /// @param line A string that contains a single OGC WKT entry. Read from a file
    /// @author SMM
    /// @date 29/01/2019
    string parse_projection(string line);

    /// @brief This parses a reference ellipse from an OGC WKT entry. We only allow certain options. See source code for those options. 
    /// @param line A string that contains a single OGC WKT entry. Read from a file
    /// @author SMM
    /// @date 29/01/2019    
    string parse_ellipsoid(string line);

    /// @brief This parses the EPSG code from an OGC WKT entry. 
    /// @param line A string that contains a single OGC WKT entry. Read from a file
    /// @author SMM
    /// @date 29/01/2019 
    int parse_EPSG(string line);

    /// @brief This parses parameters from a conic coordinate system. 
    /// @detail See Snyder 1987. The parameters are phi_0, lambda_0, phi_1 and phi_2. 
    ///  These prepresent the location of the origin (lat,long for phi and lambda)
    ///   and then the two standard parallels (these are in latitude)
    /// @param line A string that contains a single OGC WKT entry. Read from a file
    /// @author SMM
    /// @date 29/01/2019    
    map<string,double> parse_conic_parameters(string line);

    /// @brief This parses a specific keyword, defined by the user
    /// @param line A string that contains a single OGC WKT entry. Read from a file
    /// @param keyword The keyword you wnat to extract from the string
    /// @author SMM
    /// @date 29/01/2019  
    string parse_keyword(string line, string keyword);

    /// @brief this gets the projection from a specific EPSG code
    ///  WARNING no checking is the code is valid!!!!
    LSDProjectionInfo get_projection(int EPSG_code)  { return projections_map[EPSG_code];}



  protected:

    /// A map holding projections, with the EPSG code as the key
    /// see http://spatialreference.org/ for information about various EPSG codes 
    map<int, LSDProjectionInfo> projections_map;


  private:
    void create();

    void create(string filename);

    void create(string path, string filename);


};




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// A class for converting datums and coordinates
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
/// @brief This is a coordinate converter class that is relatively flexible
/// You can add different ellipsoids and coordinate systems
/// @detail You will need to hard code EPSG codes but since all the maths are now
///  in these functions it will be fairly easy to add new lambert and albers 
///  projections.
/// @author SMM
/// @date 20/01/2019
class LSDCoordinateConverter
{
  public:
    /// @brief default constructor. Assumes both a projection and ellipsoid
    /// @detail Assumed ellisoid is WGS84. Assumed projection is TBD
    /// @author SMM
    /// @date 20/01/2019
    LSDCoordinateConverter()     { create(); }

    /// @brief Constructor that takes a projection and assumes and ellipsoid
    /// @param ProjInfo the LSDProjectionInfo projection object
    /// @author SMM
    /// @date 20/01/2019
    LSDCoordinateConverter(int EPSG_code)     
                 { create(EPSG_code); }

    /// @brief Constructor that takes a projection and assumes and ellipsoid
    /// @param ProjInfo the LSDProjectionInfo projection object
    /// @param RefEllipse The LSDReferenceEllipsoid object
    /// @author SMM
    /// @date 20/01/2019
    LSDCoordinateConverter(int EPSG_code, string Ellipse_name)     
                        { create(EPSG_code, Ellipse_name); }



    /// @brief This wraps computation of the Easting and Northing from lat-long
    /// @detail The function reads the projection to determine what type of conversion to call
    /// @param Lat_Long the latitude and logitude in decimal degrees as a pair (Lat_Long.first is the latitude)
    /// @return A pair containing easting and northing (easting is the first element in the pair)
    /// @author SMM
    /// @date 21/01/2019
    pair<double,double> Convert_LL_to_EN(pair<double,double> Lat_Long);
    
    /// @brief This wraps computation of the Easting and Northing from lat-long
    /// @detail The function reads the projection to determine what type of conversion to call
    /// @param A pair containing easting and northing (easting is the first element in the pair)
    /// @return Lat_Long the latitude and logitude in decimal degrees as a pair (Lat_Long.first is the latitude)
    /// @author SMM
    /// @date 22/01/2019
    pair<double,double> Convert_EN_to_LL(pair<double,double> Easting_Northing);


    /// @brief converts LatLongHt in datum dIn, to LatLongHt in datum dTo;
    /// @detail 2002dec: by Eugene Reimer, from PeterDana equations.
    ///   Lat and Long params are in degrees;
    /// North latitudes and East longitudes are positive;  Height is in meters;
    /// ==This approach to Datum-conversion is a waste of time;
    /// to get acceptable accuracy a large table is needed -- see NADCON, NTv2...
    void DatumConvert(int dIn, double LatIn, double LongIn, double HtIn,
                  int dTo,  double& LatTo, double& LongTo, double& HtTo);

    /// @brief A forward conversion from lat long to lambert conical conformal system
    /// @detail See page 108 in Snyder 1987
    /// @param Lat_Long a pair containing the latitude and longitude
    /// @return a pair containing the false easting and false northing in local coordinates
    /// @author SMM 
    /// @date 21/01/2019
    pair<double,double> LCC_forward(pair<double,double> Lat_Long);

    /// @brief An inverse conversion from easting and northing to lat-long
    /// @detail See page 108 in Snyder 1987
    /// @param a pair containing the false easting and false northing in local coordinates
    /// @return Lat_Long a pair containing the latitude and longitude
    /// @author SMM 
    /// @date 21/01/2019
    pair<double,double> LCC_inverse(pair<double,double> Easting_Northing);

    /// @brief A forward conversion from lat long to albers equal area conic system
    /// @detail See page 101 in Snyder 1987
    /// @param Lat_Long a pair containing the latitude and longitude
    /// @return a pair containing the false easting and false northing in local coordinates
    /// @author SMM 
    /// @date 22/01/2019
    pair<double,double> AEAC_forward(pair<double,double> Lat_Long);

    /// @brief An inverse conversion from easting and northing to lat-long
    /// @detail See page 101 in Snyder 1987
    /// @param a pair containing the false easting and false northing in local coordinates
    /// @return Lat_Long a pair containing the latitude and longitude
    /// @author SMM 
    /// @date 22/01/2019
    pair<double,double> AEAC_inverse(pair<double,double> Easting_Northing);

    /// @brief A forward conversion from lat long to universal transverse mercator system
    /// @param Lat_Long a pair containing the latitude and longitude
    /// @return a pair containing the false easting and false northing in local coordinates
    /// @author SMM 
    /// @date 21/01/2019
    pair<double,double> UTM_forward(pair<double,double> Lat_Long);

    /// @brief An inverse conversion from universal transverse mercator system to lat long
    /// @param a pair containing the false easting and false northing in local coordinates
    /// @return Lat_Long a pair containing the latitude and longitude
    /// @author SMM 
    /// @date 23/01/2019
    pair<double,double> UTM_inverse(pair<double,double> Easting_Northing);


               

  protected:

    /// The reference ellipsoid
    LSDReferenceEllipsoid RefEllipsoid;

    /// The parameters 
    LSDProjectionInfo ProjParams;

    /// Some parameters for setting degree conversions
    double RADIANS_PER_DEGREE;
    double DEGREES_PER_RADIAN;



  private:

    /// @brief This create function sets up the data membeers that hold the
    ///  ellipsoid and datum data
    /// @author SMM
    /// @date 20/01/2019    
    void create();

    /// @brief Create function that takes a projection and assumes and ellipsoid
    /// @param ProjInfo the LSDProjectionInfo projection object
    /// @param RefEllipse The LSDReferenceEllipsoid object
    /// @author SMM
    /// @date 20/01/2019
    void create(int EPSG_code);

    /// @brief Create function that takes a projection and assumes and ellipsoid
    /// @param ProjInfo the LSDProjectionInfo projection object
    /// @param RefEllipse The LSDReferenceEllipsoid object
    /// @author SMM
    /// @date 20/01/2019
    void create(int EPSG_code, string Ellipse_name);

};


#endif
