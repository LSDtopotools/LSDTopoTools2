///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDLithoCube.cpp
/// This object has a stack of LSDRasters that have lithologic information
/// currently this is read from a csv supplied by NAGRA
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1    25/01/2020
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDSpatialCSVReader.hpp"
#include "LSDStatsTools.hpp"
#include "LSDLithoCube.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDLithoCube_CPP
#define LSDLithoCube_CPP



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// operators
// SMM, 2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDLithoCube& LSDLithoCube::operator=(const LSDLithoCube& rhs)
{
  if (&rhs != this)
  {
    create(rhs.get_NRows(),rhs.get_NCols(),rhs.get_NLayers(),
           rhs.get_XMinimum(),rhs.get_YMinimum(),rhs.get_ZMinimum(),          
           rhs.get_DataResolution(),
           rhs.get_XSpacing(),rhs.get_YSpacing(),rhs.get_ZSpacing(),
           rhs.get_loaded_vo_file(), rhs.get_NoDataValue(), 
           rhs.get_GeoReferencingStrings(),
           rhs.get_LithoLayers(),
           rhs.get_bottom_elevations(),
           rhs.get_layer_thicknesses(),
           rhs.get_StratiToK(),
           rhs.get_StratiToSc() );
  }
  return *this;
}

// This is the create function that is used by the copy constructor
void LSDLithoCube::create(int nrows, int ncols, int nlayers, float xmin, float ymin, float zmin,
                          float cellsize, float xspace, float yspace, float zspace, 
                          bool loaded, int ndv, map<string,string> temp_GRS,
                          vector< Array2D<int> > ll, vector<float> be, vector<float> lt, 
                          map<int,float> StK, map<int,float> StS)
{
  NRows = nrows;
  NCols = ncols;
  NLayers = nlayers, 
  XMinimum = xmin;
  YMinimum = ymin;
  ZMinimum = zmin;
  DataResolution = cellsize;
  XSpacing = xspace;
  YSpacing = yspace;
  ZSpacing = zspace;
  loaded_vo_file = loaded;
  NoDataValue = ndv;
  GeoReferencingStrings = temp_GRS;
  LithoLayers = ll;
  bottom_elevations = be;
  layer_thicknesses = lt;
  StratiToK = StK;
  StratiToSc = StS;
}


// A minimal create function. 
void LSDLithoCube::create()
{
  NRows = 100;
  NCols = 100;
  DataResolution = 10;
  NoDataValue = -9999;
  XMinimum = 0;
  YMinimum = 0;

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);

  loaded_vo_file = false;
}

// this creates a raster using an infile
void LSDLithoCube::create(string filename)
{
  cout << "Let me eat a vo file for you." << endl;
  ingest_vo_file(filename);
}


// this creates a LSDRasterModel raster from another LSDRaster
void LSDLithoCube::create(LSDRaster& An_LSDRaster)
{
  cout << "Lets get some info!" << endl;
  NRows = An_LSDRaster.get_NRows();
  NCols = An_LSDRaster.get_NCols();
  XMinimum = An_LSDRaster.get_XMinimum();
  YMinimum = An_LSDRaster.get_YMinimum();
  DataResolution = An_LSDRaster.get_DataResolution();
  NoDataValue = An_LSDRaster.get_NoDataValue();
  GeoReferencingStrings =  An_LSDRaster.get_GeoReferencingStrings();

  loaded_vo_file = false;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function imposes the mapinfo strings. It assumes UTM
// THIS HAS NOT BEEN TESTED!!!!!!!!!!!!
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::impose_georeferencing_UTM(int zone, string NorS)
{
  string str_NorS;
  string str_NorSlong;
  string cs_string_fnbit;
  if (NorS.find("N") == 0 || NorS.find("n") == 0)
  {
    str_NorS = "N";
    str_NorSlong = "North";
    cs_string_fnbit = "0";
  }
  else if (NorS.find("S") == 0 || NorS.find("s") == 0)
  {
    str_NorS = "S";
    str_NorSlong = "South";
    cs_string_fnbit = "10000000";
  }
  else
  {
    cout << "imposing georeferencing, but I didn't understand N or S, defaulting to North." << endl;
    str_NorS = "N";
    str_NorSlong = "North";
    cs_string_fnbit = "0";
  }

  string delim = ", ";
  string str_UTM = "UTM";
  string str_x  = "1";
  string str_y = "1";
  string xmin = dtoa(XMinimum);
  float YMax =  YMinimum + NRows*DataResolution;
  string ymax = dtoa(YMax);

  string DR = dtoa(DataResolution);
  string str_UTMZ = itoa(zone);
  string str_hemis = str_NorSlong;
  string str_spheroid = "WGS-84";

  string new_string = str_UTM+delim+str_x+delim+str_y+delim+xmin+delim
                       +ymax+delim+DR+delim+DR+delim+str_UTMZ+delim+str_hemis
                       +delim+str_spheroid;
  GeoReferencingStrings["ENVI_map_info"]= new_string;

  string cs_string_firstbit = "PROJCS[\"WGS_1984_UTM_Zone_";
  string cs_string_secondbit = str_UTMZ+str_NorS;
  string cs_string_thirdbit =  "\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",";
  string cs_string_fifthbit = "],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",";
  string cs_string_seventhbit = "],UNIT[\"Meter\",1]]";

  int central_meridian = Find_UTM_central_meridian(zone);
  string cs_string_central_merid = itoa(central_meridian);


  string cs_str = cs_string_firstbit+cs_string_secondbit+cs_string_thirdbit
                 +cs_string_central_merid+cs_string_fifthbit+cs_string_fnbit
                 +cs_string_seventhbit;

  GeoReferencingStrings["ENVI_coordinate_system"]= cs_str;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This is a utility function to find the central meridian of a UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDLithoCube::Find_UTM_central_meridian(int UTM_zone)
{
  // initiate the vector
  vector<int> zone(61,NoDataValue);

  // here is the lookuptable
  zone[1]=-177;
  zone[2]=-171;
  zone[3]=-165;
  zone[4]=-159;
  zone[5]=-153;
  zone[6]=-147;
  zone[7]=-141;
  zone[8]=-135;
  zone[9]=-129;
  zone[10]=-123;
  zone[11]=-117;
  zone[12]=-111;
  zone[13]=-105;
  zone[14]=-99;
  zone[15]=-93;
  zone[16]=-87;
  zone[17]=-81;
  zone[18]=-75;
  zone[19]=-69;
  zone[20]=-63;
  zone[21]=-57;
  zone[22]=-51;
  zone[23]=-45;
  zone[24]=-39;
  zone[25]=-33;
  zone[26]=-27;
  zone[27]=-21;
  zone[28]=-15;
  zone[29]=-9;
  zone[30]=-3;
  zone[31]=3;
  zone[32]=9;
  zone[33]=15;
  zone[34]=21;
  zone[35]=27;
  zone[36]=33;
  zone[37]=39;
  zone[38]=45;
  zone[39]=51;
  zone[40]=57;
  zone[41]=63;
  zone[42]=69;
  zone[43]=75;
  zone[44]=81;
  zone[45]=87;
  zone[46]=93;
  zone[47]=99;
  zone[48]=105;
  zone[49]=111;
  zone[50]=117;
  zone[51]=123;
  zone[52]=129;
  zone[53]=135;
  zone[54]=141;
  zone[55]=147;
  zone[56]=153;
  zone[57]=159;
  zone[58]=165;
  zone[59]=171;
  zone[60]=177;

  int central_meridian;

  // now look up the table
  if(UTM_zone <1 || UTM_zone > 60)
  {
    cout << "Trying to assign central meridian but you have chosen an invalid UTM zone" << endl;
    cout << "defaulting to central meridian of 0";
    central_meridian = 0;
  }
  else
  {
    central_meridian = zone[UTM_zone];
  }

  return central_meridian;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



// This ingests the data from the vo file
void LSDLithoCube::ingest_vo_file(string filename)
{
  cout << "Now I will ingest the vo file." << endl;
  // Open the file
  ifstream vo_in;
  vo_in.open(filename.c_str());

  // check if the parameter file exists
  if( vo_in.fail() )
  {
    cout << "\nFATAL ERROR: The parameter file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  

  // Read in the file
  vector<string> lines;
  string str;

  int n_rows = 0;
  int n_cols = 0;
  int n_layers = 0;
  float x_dim = 0;
  float y_dim = 0;
  float z_dim = 0;
  float x_min = 0;
  float y_min = 0;
  float z_min = 0;

  string axis_n = "AXIS_N";
  string xd = "AXIS_U";
  string yd = "AXIS_V";
  string zd = "AXIS_W";
  string ao = "AXIS_O";

  while (std::getline(vo_in, str)) 
  {
    string this_line = strip_and_clean_string(str);
    //cout << endl << "This line's string is: " << this_line << endl;
    
    // print out the first word of the line
    int sz = this_line.size();
    vector<string> stringvec;
    //string c_space = " ";
    char space =' ';
    if (sz != 0)
    {
      split_delimited_string( this_line, space, stringvec);
      //cout << "The string is: " << this_line << endl; 
      if (stringvec.size() > 0)
      {
        if (stringvec[0] == axis_n) 
        {
          //cout << "Found axis_N" << endl;
          n_cols = atoi(stringvec[1].c_str());
          n_rows = atoi(stringvec[2].c_str());
          n_layers = atoi(stringvec[3].c_str());
        }
        if (stringvec[0] == xd) 
        {
          x_dim = atof(stringvec[1].c_str());
        }
        if (stringvec[0] == yd) 
        {
          y_dim = atof(stringvec[2].c_str());
        }
        if (stringvec[0] == zd)
        {
          z_dim = atof(stringvec[3].c_str());
        }
        if (stringvec[0] == ao) 
        {
          x_min = atof(stringvec[1].c_str());
          y_min = atof(stringvec[2].c_str());
          z_min = atof(stringvec[3].c_str());
        }  
      }
    }    
  }

  cout.precision(9);
  cout << "Here are the vitalstatistix, chief!" << endl;
  cout << "NRows: " << n_rows << " NCols: " << n_cols << " NLayers: " << n_layers << endl;
  cout << "xdim: " << x_dim << " ydim: " << y_dim << " z_dim: " << z_dim << endl;
  cout << "x_min: " << x_min << " x_min: " << y_min << " z_min: " << z_min << endl; 

  NRows = n_rows;
  NCols = n_cols;
  NLayers = n_layers;
  XMinimum = x_min;
  YMinimum = y_min;
  ZMinimum = z_min;

  XSpacing = x_dim/(float(n_cols-1));
  YSpacing = y_dim/(float(n_rows-1));
  ZSpacing = z_dim/(float(n_layers-1));  

  cout <<"dx: " << XSpacing << " dy: " << YSpacing << " dz: " << ZSpacing << endl;



  DataResolution = XSpacing;

  NoDataValue = -9999;

  loaded_vo_file = true; 
  cout << "The loaded vo file bool is: " << loaded_vo_file << endl; 

}

// This ingests the csv data
void LSDLithoCube::ingest_litho_data(string filename)
{
  cout << "The loaded vo file bool is: " << loaded_vo_file << endl;
  
  if (loaded_vo_file)
  {
    cout << "Great news, fellow analyist!! I have the vo file!" << endl;
  }
  else
  {
    cout << "Fatal error. The vo file has not yet been loaded" << endl;
    exit(EXIT_FAILURE);
  }
  

  cout << "Let me build the stack of arrays." << endl;
  Array2D<int> temp_raster(NRows,NCols,NoDataValue);

  vector< Array2D<int> > ResetLithoLayers;
  LithoLayers = ResetLithoLayers;
  for (int layer = 0; layer< NLayers; layer++)
  {
    LithoLayers.push_back(temp_raster.copy());
  }

  cout << "I've made the data element for the data. Now I will fill it with the asc file." << endl;
  // Open the file
  ifstream data_in;
  data_in.open(filename.c_str());

  // check if the parameter file exists
  if( data_in.fail() )
  {
    cout << "\nFATAL ERROR: The file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }    

  // get the header. Print and discard.
  string str;
  getline(data_in,str);
  //getline(data_in,str);
  cout << "The header is: " << str << endl;
  getline(data_in,str);
  cout << "Let me get the data for you." << endl;
  cout << "If this is a big file loading it will take a little time." << endl;

  int this_row, this_col, this_layer, this_strat;
  
  //for (int i = 0; i< 10; i++)
  while (getline(data_in, str)) 
  {
    string this_line = strip_and_clean_string(str);    

    int sz = this_line.size();
    vector<string> stringvec;
    char comma =',';
    if (sz != 0)
    {
      split_delimited_string( this_line, comma, stringvec);
      //cout << "The string is: " << this_line << endl; 
      // Note: the csv input file has I,J,K indices, which are the X,Y,Z columns
      // The J index is the node and the I index is the row
      // Need to flip the index because the bottom row (J = 0) is NRows-1. 
      // and the top row is J = NRows-1
      if (stringvec.size() > 0)
      {
        this_row = NRows - atoi(stringvec[1].c_str()) -1;
        this_col = atoi(stringvec[0].c_str());
        this_layer = atoi(stringvec[2].c_str());
        this_strat = atoi(stringvec[3].c_str());

        //cout << "R: " << this_row << " C: " << this_col << " L: " << this_layer << " S: " << this_strat << endl;
        LithoLayers[this_layer][this_row][this_col] = this_strat;
      }
    }
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints a layter as a raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDLithoCube::get_raster_from_layer(int layer)
{
  Array2D<int> this_layer;
  
  if (layer <0 || layer >= NLayers)
  {
    cout << "FATAL ERROR You have seleced a layter that doesn't exist." << endl;
    cout << "Cou asked for " << layer << " and we only have " << NLayers << " layers." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    this_layer = LithoLayers[layer].copy();
  }

  LSDIndexRaster Layer_Raster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, this_layer, GeoReferencingStrings);
  return Layer_Raster;
  
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints a layer as a raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::write_litho_raster_from_layer(int layer, string fname, string f_ext)
{
  Array2D<int> this_layer;
  
  if (layer <0 || layer >= NLayers)
  {
    cout << "FATAL ERROR You have seleced a layter that doesn't exist." << endl;
    cout << "Cou asked for " << layer << " and we only have " << NLayers << " layers." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    this_layer = LithoLayers[layer].copy();
  }

  LSDIndexRaster Layer_Raster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, this_layer, GeoReferencingStrings);
  Layer_Raster.write_raster(fname,f_ext);
  
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints a layer as a raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::write_litho_raster_from_elevations(LSDRaster& ElevationRaster, list<int> forbidden_codes, string fname, string f_ext)
{
 
  LSDIndexRaster Layer_raster = get_lithology_codes_from_raster(ElevationRaster, forbidden_codes);
  Layer_raster.write_raster(fname,f_ext);
  
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Gets the row and column of a point
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::get_row_and_col_of_a_point(float X_coordinate,float Y_coordinate,int& row, int& col)
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Checks to see is a point is in the raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDLithoCube::check_if_point_is_in_raster(float X_coordinate,float Y_coordinate)
{
  bool is_in_raster = true;

  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1)
  {
    is_in_raster = false;
  }

  return is_in_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the layer number from an elevation
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDLithoCube::get_layer_from_elevation(float elevation)
{
  //float ZMaximum = float(NLayers)*ZSpacing+ZMinimum;
  int this_layer = 0;

  this_layer = floor( (elevation-ZMinimum)/ZSpacing );

  if (this_layer < 0)
  {
    this_layer = 0;
  }
  if (this_layer > NLayers-1 ) 
  {
    this_layer = NLayers-1;
  }
  
  return this_layer;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the elevation from a layer number
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDLithoCube::get_elevation_from_layer(int layer_number)
{
  float this_elevation;

  this_elevation = ZMinimum + (layer_number * ZSpacing) + 0.5 * ZSpacing;
  cout << "Layer is: " << layer_number << " and elevation is: " << this_elevation << endl;

  return this_elevation;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets an elevation raster of the lithocube surface.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDLithoCube::get_lithocube_surface_elevation(LSDIndexRaster& SurfaceLayerRaster)
{
  int SLR_NRows = SurfaceLayerRaster.get_NRows();
  int SLR_NCols = SurfaceLayerRaster.get_NCols();
  float SLR_NDV = SurfaceLayerRaster.get_NoDataValue();

  Array2D<float> elevation(SLR_NRows,SLR_NCols,SLR_NDV);


  // now loop through topmost layer raster, getting the information
  for(int row = 0; row < SLR_NRows; row++)
  {
    for(int col = 0; col < SLR_NCols; col++)
    {
      elevation[row][col] = get_elevation_from_layer(SurfaceLayerRaster.get_data_element(row,col));    
    }
  }
  cout << "I have assigned the elevations." << endl;
  // now make the raster of the stratigraphic information
  LSDRaster ThisSurfaceElevation(SLR_NRows, SLR_NCols, SurfaceLayerRaster.get_XMinimum(), SurfaceLayerRaster.get_YMinimum(), 
                           SurfaceLayerRaster.get_DataResolution(), SLR_NDV, elevation, SurfaceLayerRaster.get_GeoReferencingStrings());

  return ThisSurfaceElevation;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the topmost layer in the lithocube that represents stratigraphy 
// (i.e. for each location, gets the layer number that represents the surface)
// It has logic to "drill" down the lithocube so that it returns the 
// topmost allowed values
// the forbidden values are in the list called forbidden_codes
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDLithoCube::get_lithocube_surface_layers(list<int> forbidden_codes)
{
  int starting_layer = NLayers-1;
  Array2D<int> layer_number(NRows,NCols,NoDataValue);

  cout << "Starting at layer " << starting_layer << endl;


  // now loop through topmost layer raster, getting the information
  for(int row = 0; row < NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      layer_number[row][col] = get_layer_from_location(starting_layer, row, col, forbidden_codes);    
    }
  }
  cout << "I have assigned the layer numbers." << endl;
  // now make the raster of the stratigraphic information
  LSDIndexRaster ThisSurfaceLayer(NRows, NCols, XMinimum, YMinimum, XSpacing, NoDataValue, layer_number, GeoReferencingStrings);


  return ThisSurfaceLayer;


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets a lithocode from a location
// It has logic to "drill" down the lithocube so that it returns the 
// topmost allowed values
// the forbidden values are in the list called forbidden_codes
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDLithoCube::get_lithocode_from_location(int LLayer, int LRow, int LCol, list<int> forbidden_codes)
{

  int this_lcode = LithoLayers[LLayer][LRow][LCol];

  //cout << "Here are the forbidden lithocodes" << endl;
  //list<int>::iterator it;
  //for (it = forbidden_codes.begin(); it != forbidden_codes.end(); ++it)
  // {
  //  cout << *it << endl;
  //}


  if (forbidden_codes.size() > 0)
  {
    
    int this_layer = LLayer;
    // cout << "This lithocode is: " << this_lcode << endl;
    while (  (find(forbidden_codes.begin(),forbidden_codes.end(),this_lcode) != forbidden_codes.end()) && this_layer > 0)
    {
      // cout << "Got a forbidden code! Let me dig a layer!" << endl;
      this_layer--;
      this_lcode = LithoLayers[this_layer][LRow][LCol];
      // cout << "New layer is: " << this_layer << " and new code is: " << this_lcode << endl;
    } 
  }
    
  return this_lcode;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the layer number from a location
// It has logic to "drill" down the lithocube so that it returns the 
// topmost layer with allowed lithocode values
// the forbidden values are in the list called forbidden_codes
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDLithoCube::get_layer_from_location(int LLayer, int LRow, int LCol, list<int> forbidden_codes)
{

  int this_lcode = LithoLayers[LLayer][LRow][LCol];
  int this_layer = LLayer;

 
  if (forbidden_codes.size() > 0)
  {
    cout << "This lithocode is: " << this_lcode << endl;
    while (  (find(forbidden_codes.begin(),forbidden_codes.end(),this_lcode) != forbidden_codes.end()) && this_layer > 0)
    {
      cout << "Got a forbidden code! Let me dig a layer!" << endl;
      this_layer--;
      this_lcode = LithoLayers[this_layer][LRow][LCol];
      cout << "New layer is: " << this_layer << " and new code is: " << this_lcode << endl;
    } 
  }
    
  return this_layer;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This takes an elevation raster and gets a stratigraphic layer raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDLithoCube::get_lithology_codes_from_raster(LSDRaster& ElevationRaster, list<int> forbidden_codes)
{
  vector<float> Eastings;
  vector<float> Northings;
  vector<float> New_northings;
  ElevationRaster.get_easting_and_northing_vectors(Eastings, Northings);
  int ER_NRows = ElevationRaster.get_NRows();
  int ER_NCols = ElevationRaster.get_NCols(); 
  float ER_NDV = ElevationRaster.get_NoDataValue();
  int LR_NDV = -9999;
  int LRow, LCol;
  int this_layer = 0;

  //cout << "NoDataValue in the elevation raster is " << ER_NDV << endl;

  Array2D<int> lithocode(ER_NRows,ER_NCols,LR_NDV);

  // // we need to reverse the northing vecor
  // int n_north = int(Northings.size());
  // for (int i = 0; i<n_north; i++)
  // {
  //   New_northings.push_back(Northings[ER_NRows-1-i]);
  // }

  std::map<int,float> furb;
  for(map<int,float>::iterator gurg = StratiToK.begin(); gurg != StratiToK.end(); gurg++ )

  // now loop through elevation raster, getting the information
  for(int row = 0; row < ER_NRows; row++)
  {
    for(int col = 0; col < ER_NCols; col++)
    {
      // Check if there is elevation data
      // cout << "A";
      if (ElevationRaster.get_data_element(row,col) != ER_NDV)
      {
        
        // check if the point is in the litho cube
        // cout << "B";
        
        if ( check_if_point_is_in_raster(Eastings[col],Northings[row]))
        {
          //cout << "I am in row " << row << " and col " << col << ". The elevation is " << ElevationRaster.get_data_element(row,col);
          // now get the row and column of the point in the elevation raster within the litho cube
          // get_row_and_col_of_a_point(Eastings[row],New_northings[col],LRow, LCol);

          // cout << "C";
          XY_to_rowcol_careful(Eastings[col],Northings[row],LRow, LCol, true);
          if (LRow == -9999 || LCol == -9999)
          {
            cout << "I am afraid this pixel does not seem to be within the lithocube.";
          }
          // cout << "D";
          this_layer = get_layer_from_elevation(ElevationRaster.get_data_element(row,col)); // this is where you change the relative uplift in a hacky way yo!!!!! < ----
          // cout << "E";

          // if(LRow < 0 || LRow >= NRows || LCol < 0 || LCol>= NCols )
          //cout << "LRow is " << LRow << " LCol is" << LCol << " row is" << row << " col is " << col << std::endl;
          lithocode[row][col] = get_lithocode_from_location(this_layer, LRow, LCol, forbidden_codes);
          // cout << "F";
        }
        else
        {
          // cout << "I am in row " << row << " and col " << col << ". The elevation is " << ElevationRaster.get_data_element(row,col) << endl;
          // cout << "This point is not in the lithocube." << endl;
        }

      }
    }
  }
  //cout << "I have assigned the lithocodes." << endl;
  // now make the raster of the stratigraphic information
  LSDIndexRaster ThisStrat(ER_NRows, ER_NCols, ElevationRaster.get_XMinimum(), ElevationRaster.get_YMinimum(), 
                           ElevationRaster.get_DataResolution(), LR_NDV, lithocode, ElevationRaster.get_GeoReferencingStrings());


  return ThisStrat;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This takes an elevation raster and gets an erodibility value raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDLithoCube::get_K_raster_from_raster(LSDRaster& ElevationRaster )
{
  cout << "Let's get a K raster!" << endl;
  vector<float> Eastings;
  vector<float> Northings;
  vector<float> New_northings;
  ElevationRaster.get_easting_and_northing_vectors(Eastings, Northings);
  int ER_NRows = int(Northings.size());
  int ER_NCols = int(Eastings.size()); 
  float ER_NDV = ElevationRaster.get_NoDataValue();
  int LR_NDV = -9999;
  int LRow, LCol;
  int this_layer = 2;

   
  // initialises array - should I initialise a K array straight away? 
  Array2D<int> lithocode(ER_NRows,ER_NCols,LR_NDV);

  // now loop through elevation raster, getting the information
  for(int row = 0; row < ER_NRows; row++)
  {
    for(int col = 0; col < ER_NCols; col++)
    {
      // Check if there is elevation data
      if (ElevationRaster.get_data_element(row,col) != ER_NDV)
      {
        // check if the point is in the litho cube
        if ( check_if_point_is_in_raster(Eastings[col],New_northings[row]))
        {
          // now get the row and column of the point in the elevation raster within the litho cube
          get_row_and_col_of_a_point(Eastings[col],New_northings[row],LRow, LCol);

          lithocode[row][col] = LithoLayers[this_layer][LRow][LCol];
        }
      }
    }
  }

  // now make the raster of the stratigraphic information
  LSDIndexRaster ThisStrat(ER_NRows, ER_NCols, ElevationRaster.get_XMinimum(), ElevationRaster.get_YMinimum(), 
                           ElevationRaster.get_DataResolution(), LR_NDV, lithocode, ElevationRaster.get_GeoReferencingStrings());


  return ThisStrat;
}


LSDRaster LSDLithoCube::index_to_K(LSDIndexRaster& index_raster, float default_K)
{
  Array2D<float> data_K(index_raster.get_NRows(),index_raster.get_NCols(), default_K);

  for(int row=0;row<index_raster.get_NRows();row++)
  for(int col=0;col<index_raster.get_NCols();col++)
  {
    if(StratiToK.count(index_raster.get_data_element(row,col))>0)
    {
      //std:: cout << "this is happening and not default :: " << index_raster.get_data_element(row,col);
      data_K[row][col] = StratiToK[index_raster.get_data_element(row,col)];
    }
    else
      data_K[row][col] = default_K;
  }

  LSDRaster K_raster(index_raster.get_NRows(), index_raster.get_NCols(), index_raster.get_XMinimum(), index_raster.get_YMinimum(),
      index_raster.get_DataResolution(), NoDataValue, data_K, index_raster.get_GeoReferencingStrings());

  return K_raster;
}

LSDRaster LSDLithoCube::index_to_Sc(LSDIndexRaster& index_raster, float default_Sc)
{
  Array2D<float> data_Sc(index_raster.get_NRows(),index_raster.get_NCols(), default_Sc);

  for(int row=0;row<index_raster.get_NRows();row++)
  for(int col=0;col<index_raster.get_NCols();col++)
  {
    if(StratiToSc.count(index_raster.get_data_element(row,col))>0)
    {
      //std:: cout << "this is happening and not default :: " << index_raster.get_data_element(row,col);
      data_Sc[row][col] = StratiToSc[index_raster.get_data_element(row,col)];
    }
    else
      data_Sc[row][col] = default_Sc;
  }

  LSDRaster Sc_raster(index_raster.get_NRows(), index_raster.get_NCols(), index_raster.get_XMinimum(), index_raster.get_YMinimum(),
      index_raster.get_DataResolution(), NoDataValue, data_Sc, index_raster.get_GeoReferencingStrings());

  return Sc_raster;
}


// Does what the name suggests
void LSDLithoCube::write_K_raster_from_elevations(LSDRaster& ElevationRaster, string fname, string f_ext)
{

  LSDRaster K_raster = get_K_raster_from_raster(ElevationRaster );
  K_raster.write_raster(fname,f_ext);
  
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Reads in a csv file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::read_strati_to_K_csv(string lookup_filename) 
{

  cout << "I am reading a lithocode to K value file." << endl;
  cout << "This needs to be a csv, and it uses a header." << endl;
  
  // Open the file
  ifstream data_in;
  data_in.open(lookup_filename.c_str());

  // check if the parameter file exists
  if( data_in.fail() )
  {
    cout << "\nFATAL ERROR: The file \"" << lookup_filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }    

  // get the header. Print and discard.
  string str;
  getline(data_in,str);
  cout << "The header is: " << str << endl;
  cout << "Let me get the data for you." << endl;
  
  int this_strati_code;
  float this_K;

  map<int,float> TempStratiToK;
  
  while (getline(data_in, str)) 
  {
    string this_line = strip_and_clean_string(str);    

    int sz = this_line.size();
    vector<string> stringvec;
    char comma =',';
    if (sz != 0)
    {
      split_delimited_string( this_line, comma, stringvec);
      cout << "The string is: " << this_line << endl; 
      if (stringvec.size() > 0)
      {
        setprecision(9);
        this_strati_code = stoi(stringvec[0]);
        this_K = stof(stringvec[1]);
      
      
        cout << "Strati code: " << this_strati_code << " K: " << this_K << endl;
        TempStratiToK[this_strati_code] = this_K;
      }
    }
  }

  StratiToK = TempStratiToK;

}

void LSDLithoCube::read_strati_to_Sc_csv(string lookup_filename) 
{

  cout << "I am reading a lithocode to K value file." << endl;
  cout << "This needs to be a csv, and it uses a header." << endl;

  // Open the file
  ifstream data_in;
  data_in.open(lookup_filename.c_str());

  // check if the parameter file exists
  if( data_in.fail() )
  {
    cout << "\nFATAL ERROR: The file \"" << lookup_filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }    

  // get the header. Print and discard.
  string str;
  getline(data_in,str);
  cout << "The header is: " << str << endl;
  cout << "Let me get the data for you." << endl;
  
  int this_strati_code;
  float this_Sc;

  map<int,float> TempStratiToSc;
  
  while (getline(data_in, str)) 
  {
    string this_line = strip_and_clean_string(str);    

    int sz = this_line.size();
    vector<string> stringvec;
    char comma =',';
    if (sz != 0)
    {
      split_delimited_string( this_line, comma, stringvec);
      cout << "The string is: " << this_line << endl; 
      if (stringvec.size() > 0)
      {
        setprecision(9);
        this_strati_code = stoi(stringvec[0]);
        this_Sc = stof(stringvec[1]);
      
      
        cout << "Strati code: " << this_strati_code << " Sc: " << this_Sc << endl;
        TempStratiToSc[this_strati_code] = this_Sc;
      }
    }
  }
  StratiToSc = TempStratiToSc;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Prints stratigraphic codes to screenb
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::print_K_and_Sc_maps_to_screen() 
{
  cout << "Hello, I cannot wait to tell you the K and Sc values for the stratigraphic codes. " << endl;
  cout << "Here they are!!" << endl;
  cout << "The K values are: " << endl;

  for(auto it = StratiToK.begin(); it != StratiToK.end(); ++it)
  {
    cout << "Strati code: " << it->first << " => K value: " << it->second << endl;
  }

  cout << endl << endl << "and the Sc values are:" << endl;
  for(auto it = StratiToSc.begin(); it != StratiToSc.end(); ++it)
  {
    cout << "Strati code: " << it->first << " => Sc value: " << it->second << endl;
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Fills a map with stratigraphy codes and corresponding K value
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::hardcode_fill_strati_to_K_map() {

  // SMM I have updated this with Andreas' codes on 23/03/2020
  // This is from the lithocode strigraphy table:
  // Strati	Erod_Class
  //0	-1
  //1	0
  //2	2
  //3	2
  //4	2
  //5	3
  //6	2
  //7	2
  //8	3
  //9	2
  //10	2
  //11	1
  //12	1
  //13	3
  //14	2
  //15	3
  //-99999	2

  // And this is from Andreas' email
  // -	1: SC = 0.21
  // -	2: SC = 0.27
  // -	3 and 4: SC = 0.45

  // -	1: K = 2.0 x 10e-6
  // -	2: K = 1.5 x 10e-6
  // -	3 and 4: K = 1.0 x 10e-6


  //StratiToK[0] = -9999;
  StratiToK[0] = 0.000003; //Slightly dodgy way to ensure erosion still happens when we are above the lithocube surface
  StratiToK[1] = 0.000003;
  StratiToK[2] = 0.0000015;
  StratiToK[3] = 0.0000015;
  StratiToK[4] = 0.0000015;
  StratiToK[5] = 0.000001;
  StratiToK[6] = 0.0000015;
  StratiToK[7] = 0.0000015;
  StratiToK[8] = 0.000001;
  StratiToK[9] = 0.0000015;
  StratiToK[10] = 0.0000015;
  StratiToK[11] = 0.000002;
  StratiToK[12] = 0.000002;
  StratiToK[13] = 0.000001;
  StratiToK[14] = 0.0000015;
  StratiToK[15] = 0.000001;
  StratiToK[-9999] = 0.0000015;


  for(auto it = StratiToK.begin(); it != StratiToK.end(); ++it)
  {
    cout << "Strati code: " << it->first << " => K value: " << it->second << endl;
  }

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Fills a map with stratigraphy codes and corresponding K value
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDLithoCube::hardcode_fill_strati_to_Sc_map() {

  // SMM I have updated this with Andreas' codes on 23/03/2020
  // This is from the lithocode strigraphy table:
  // Strati	Erod_Class
  //0	-1
  //1	0
  //2	2
  //3	2
  //4	2
  //5	3
  //6	2
  //7	2
  //8	3
  //9	2
  //10	2
  //11	1
  //12	1
  //13	3
  //14	2
  //15	3
  //-99999	2

  // And this is from Andreas' email
  // -	1: SC = 0.21
  // -	2: SC = 0.27
  // -	3 and 4: SC = 0.45

  // -	1: K = 2.0 x 10e-6
  // -	2: K = 1.5 x 10e-6
  // -	3 and 4: K = 1.0 x 10e-6

  //StratiToSc[0] = -9999;
  StratiToSc[0] = 0.21; //Slightly dodgy way to ensure erosion still happens when we are above the lithocube surface
  StratiToSc[1] = 0.21;
  StratiToSc[2] = 0.27;
  StratiToSc[3] = 0.27;
  StratiToSc[4] = 0.27;
  StratiToSc[5] = 0.45;
  StratiToSc[6] = 0.27;
  StratiToSc[7] = 0.27;
  StratiToSc[8] = 0.45;
  StratiToSc[9] = 0.27;
  StratiToSc[10] = 0.27;
  StratiToSc[11] = 0.21;
  StratiToSc[12] = 0.21;
  StratiToSc[13] = 0.45;
  StratiToSc[14] = 0.27;
  StratiToSc[15] = 0.45;
  StratiToSc[-9999] = 0.27;


  for(auto it = StratiToSc.begin(); it != StratiToSc.end(); ++it)
  {
    cout << "Strati code: " << it->first << " => Sc value: " << it->second << endl;
  }

}


void LSDLithoCube::XY_to_rowcol_careful(float X_coord, float Y_coord, int& row, int& col, bool snap_to_closest_point_on_raster)
{
  row = -9999;
  col = -9999;
  // First, shifting to origin
  bool is_X_recasted = false, is_Y_recasted = false;
  int ncols = NCols;
  int nrows = NRows;
  float cellsize = DataResolution;
  float xmin = XMinimum;
  float ymin = YMinimum;
  // Deailing wih point X
  float corrected_X = X_coord - xmin;
  if(corrected_X<0)
  {
    if(snap_to_closest_point_on_raster)
    {
      is_X_recasted = true;
      //std::cout << "WARNING::X coordinate is West to the X minimum by " << corrected_X << ", I am recasting to column 0" << std::endl;
      col = 0;
    }
    else
    {
      std::cout << "FATALERROR::Point offset from the raster by " << corrected_X << " unit West." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if(corrected_X > (ncols+1) * cellsize)
  {
    if(snap_to_closest_point_on_raster)
    {
      is_X_recasted = true;
      //std::cout << "WARNING::X coordinate is East to the X minimum by " << corrected_X - ((ncols+1) * cellsize) << ", I am recasting to column " << ncols - 1 << std::endl;
      col = ncols - 1;
    }
    else
    {
      std::cout << "FATALERROR::Point offset from the raster by " << corrected_X - ((ncols+1) * cellsize) << " unit East." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Deailing wih point Y
  float corrected_Y = Y_coord - ymin;
  if(corrected_Y < 0)
  {
    if(snap_to_closest_point_on_raster)
    {
      is_Y_recasted = true;
      //std::cout << "WARNING::Y coordinate is South to the Y minimum by " << - corrected_Y << ", I am recasting to row " << nrows - 1 << std::endl;
      row = nrows - 1;
    }
    else
    {
      std::cout << "FATALERROR::Point offset from the raster by " << - corrected_Y << " unit West." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if(corrected_Y > (nrows+1) * cellsize)
  {
    if(snap_to_closest_point_on_raster)
    {
      is_Y_recasted = true;
      //std::cout << "WARNING::X coordinate is North to the Y maximum by " <<  (nrows+1) * cellsize - corrected_Y << ", I am recasting to row " << 0 << std::endl;
      row = 0;
    }
    else
    {
      std::cout << "FATALERROR::Point offset from the raster by " <<   ((nrows+1) * cellsize) - corrected_Y << " unit East." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Done if I recasted my coordinates
  if(is_Y_recasted && is_X_recasted)
    return;

  // Otherwise I need to convert them to row col
  if(is_X_recasted == false)
  {
    col = floor(corrected_X/cellsize);
  }
  if(is_Y_recasted == false)
  {
    row = (nrows - 1) - floor(corrected_Y/cellsize);
  }

  //DOne





}



#endif
