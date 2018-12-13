///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// LSDCosmoRaster.cpp
/// header for the CosmoRaster object
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
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDCosmoRaster.hpp"
#include "LSDStatsTools.hpp"
#include "LSDCosmoRaster.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDParticle.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDCosmoRaster_CPP
#define LSDCosmoRaster_CPP




void LSDCosmoRaster::create()
{
  NRows = 100;
  NCols = 100;
  DataResolution = 10;
  NoDataValue = -9999;
  XMinimum = 0;
  YMinimum = 0;
  RasterData = Array2D <float> (NRows, NCols, 0.0);

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);
}

// this creates a raster using an infile
void LSDCosmoRaster::create(string filename, string extension)
{
  read_raster(filename,extension);
}



void LSDCosmoRaster::create(int NRows, int NCols)
{
  this->NRows = NRows;
  this->NCols = NCols;
  this->DataResolution = 10;
  this->NoDataValue = -9999;
  XMinimum = 0;
  YMinimum = 0;
  RasterData = Array2D <float> (NRows, NCols, 0.0);

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);
}


// this creates a LSDRasterModel raster from another LSDRaster
void LSDCosmoRaster::create(LSDRaster& An_LSDRaster)
{
  NRows = An_LSDRaster.get_NRows();
  NCols = An_LSDRaster.get_NCols();
  XMinimum = An_LSDRaster.get_XMinimum();
  YMinimum = An_LSDRaster.get_YMinimum();
  DataResolution = An_LSDRaster.get_DataResolution();
  NoDataValue = An_LSDRaster.get_NoDataValue();
  GeoReferencingStrings =  An_LSDRaster.get_GeoReferencingStrings();
  RasterData = An_LSDRaster.get_RasterData();
}





// This returns the data in the raster model as a raster
LSDRaster LSDCosmoRaster::return_as_raster()
{
  LSDRaster NewRaster(NRows, NCols, XMinimum, YMinimum,
                      DataResolution, NoDataValue, RasterData,
                      GeoReferencingStrings);
  return NewRaster;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This resizes and resets the model
// This overloaded version also resets the data resolution
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoRaster::resize_and_reset( int new_rows, int new_cols, float new_resolution, float new_value )
{
  // set up some empty arrays
  Array2D<float> empty_array_sized(new_rows,new_cols,new_value);

  // reset the size of the RasterData
  RasterData = empty_array_sized.copy();

  // reset the rows and columns
  NRows = new_rows;
  NCols = new_cols;

  DataResolution = new_resolution;

  int zone = 1;
  string NorS = "N";
  impose_georeferencing_UTM(zone, NorS);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets the minimum and maximum values in the raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDCosmoRaster::minimum_and_maximum_value()
{

  // The vector min max will contain minimum and maximum values. It is initiated
  // with a very high minimum and a very low maximum to guarantee that one will
  // always get sensible if the raster has non-nodata values.
  vector<float> min_max;
  min_max.push_back(1e12);
  min_max.push_back(-9998.0);

  for(int row = 0; row < NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if(RasterData[row][col] != NoDataValue)
      {
        if(RasterData[row][col] > min_max[1])
        {
          min_max[1]= RasterData[row][col];
        }
        if(RasterData[row][col] < min_max[0])
        {
          min_max[0]= RasterData[row][col];
        }
      }
    }
  }
  return min_max;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This clips to a smaller raster. The smaller raster does not need
// to have the same data resolution as the old raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDCosmoRaster::clip_to_smaller_raster(LSDRaster& smaller_raster)
{
  // Get the MinX, MaxX, MinY, MaxY from the rasters
  //float XMaximum = XMinimum + (NCols * DataResolution -1);
  //float YMaximum = YMinimum + (NRows * DataResolution -1);

  float SR_XMinimum = smaller_raster.get_XMinimum();
  float SR_YMinimum = smaller_raster.get_YMinimum();

  float SR_NRows = smaller_raster.get_NRows();
  float SR_NCols = smaller_raster.get_NCols();
  float SR_DataR = smaller_raster.get_DataResolution();

  float SR_XMaximum = SR_XMinimum+(SR_NCols)*SR_DataR;
  float SR_YMaximum = SR_YMinimum+(SR_NRows)*SR_DataR;

  cout << "Small Xmin: " << SR_XMinimum << " YMin: " << SR_YMinimum << " Xmax: "
       << SR_XMaximum << " YMax: " << SR_YMaximum << endl;

  cout << "This data resolution: " << DataResolution << " and smaller raster data resolution: " << SR_DataR << endl;


  // find the col of old raster that has the same Xlocations as the XLL of smaller raster
  // the 0.5*DataResolution is in case of rounding errors
  int XLL_col = int((SR_XMinimum-XMinimum+0.5*DataResolution)/DataResolution);
  int XUL_col = int((SR_XMaximum-XMinimum+0.5*DataResolution)/DataResolution);

  // check these columns
  if (XLL_col < 0)
  {
    XLL_col = 0;
  }
  if (XUL_col >= NCols)
  {
    XUL_col = NCols-1;
  }

  // find the row of old raster that has the same Xlocations as the XLL of smaller raster
  // the 0.5*DataResolution is in case of rounding errors
  // Slightly different logic for y because the DEM starts from the top corner
  int YLL_row = NRows - int((SR_YMinimum-YMinimum+0.5*DataResolution)/DataResolution);
  int YUL_row = NRows - int((SR_YMaximum-YMinimum+0.5*DataResolution)/DataResolution);

  // check on the lower row:
  cout << "Checking lower left row." << endl;
  cout << "integer subtraction: " << int((SR_YMinimum-YMinimum+0.5*DataResolution)/DataResolution) << endl;
  cout << "float subtraction: " <<  (SR_YMinimum-YMinimum+0.5*DataResolution)/DataResolution << endl;

  // this catches a weird rounding error.
  double int_sub =  double(int((SR_YMinimum-YMinimum+0.5*DataResolution)/DataResolution));
  double flt_sub =  (SR_YMinimum-YMinimum+0.5*DataResolution)/DataResolution;

  if ((flt_sub- int_sub) > 0.9975)
  {
    YUL_row = YUL_row+1;
  }


  // check these rows
  if (YLL_row < 0)
  {
    YLL_row = 0;
  }
  if (YUL_row >= NRows)
  {
    YUL_row = NRows-1;
  }

  cout << "Small XLLCol: " << XLL_col << " XLR_col: " << XUL_col << " YLLrow: "
       << YLL_row << " YUL_row: " << YUL_row << endl;


  // get the new number of rows and columns:
  int New_NRows = YLL_row-YUL_row;
  int New_NCols = XUL_col-XLL_col;

  cout << "New NRows: " << New_NRows  << " New_NCols: " << New_NCols << endl;

  // now extract the data for the new raster
  float NewR_XMinimum = XMinimum+float(XLL_col)*DataResolution;
  float NewR_YMinimum = YMinimum + ((NRows - YLL_row ) * DataResolution);


  Array2D<float> NewData(New_NRows,New_NCols, NoDataValue);
  for(int row = 0; row< New_NRows; row++)
  {
    for(int col = 0; col<New_NCols; col++)
    {
       NewData[row][col] = RasterData[row+YUL_row][col+XLL_col];
    }
  }

  LSDRaster TrimmedRaster(New_NRows, New_NCols, NewR_XMinimum,
                          NewR_YMinimum, DataResolution, NoDataValue, NewData,
                          GeoReferencingStrings);

  cout << "Made a new raster" << endl;

  TrimmedRaster.Update_GeoReferencingStrings();

  return TrimmedRaster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the minimum and maximum values and returns them as a vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDCosmoRaster::get_XY_MinMax()
{

  vector<float> XYMaxMin(4,0);

  XYMaxMin[0] = XMinimum;
  XYMaxMin[1] = YMinimum;
  XYMaxMin[2] = XMinimum+(NCols)*DataResolution;
  XYMaxMin[3] = YMinimum+(NRows)*DataResolution;

  return XYMaxMin;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the raster data into a vector
// FJC 06/11/15
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDCosmoRaster::get_RasterData_vector()
{
  vector<float> Raster_vector;
  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      Raster_vector.push_back(RasterData[row][col]);
    }
  }

  return Raster_vector;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function writes the raster data (where != NoDataValue) to a text file
// FJC 30/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDCosmoRaster::write_RasterData_to_text_file(string filename)
{
  string string_filename;
  string dot = ".";
  string extension = "txt";
  string_filename = filename+dot+extension;
  cout << "The filename is " << string_filename << endl;

  // open the data file
  ofstream data_out(string_filename.c_str());

  if( data_out.fail() )
  {
    cout << "\nFATAL ERROR: unable to write to " << string_filename << endl;
    exit(EXIT_FAILURE);
  }

  data_out << "row col raster_data" << endl;
  for (int row = 0; row < NRows; row++)
  {
    for (int col =0; col < NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        data_out << row << " " << col << " " << RasterData[row][col] << endl;
      }
    }
  }

  data_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This clips to a smaller raster. The smaller raster does not need
// to have the same data resolution as the old raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDCosmoRaster::clip_to_smaller_raster(LSDIndexRaster& smaller_raster)
{
  // Get the MinX, MaxX, MinY, MaxY from the rasters
  //float XMaximum = XMinimum + (NCols * DataResolution -1);
  //float YMaximum = YMinimum + (NRows * DataResolution -1);

  float SR_XMinimum = smaller_raster.get_XMinimum();
  float SR_YMinimum = smaller_raster.get_YMinimum();

  float SR_NRows = smaller_raster.get_NRows();
  float SR_NCols = smaller_raster.get_NCols();
  float SR_DataR = smaller_raster.get_DataResolution();

  float SR_XMaximum = SR_XMinimum+(SR_NCols)*SR_DataR;
  float SR_YMaximum = SR_YMinimum+(SR_NRows)*SR_DataR;

  cout << "Small Xmin: " << SR_XMinimum << " YMin: " << SR_YMinimum << " Xmax: "
       << SR_XMaximum << " YMax: " << SR_YMaximum << endl;


  // find the col of old raster that has the same Xlocations as the XLL of smaller raster
  // the 0.5*DataResolution is in case of rounding errors
  int XLL_col = int((SR_XMinimum-XMinimum+0.5*DataResolution)/DataResolution);
  int XUL_col = int((SR_XMaximum-XMinimum+0.5*DataResolution)/DataResolution);

  // check these columns
  if (XLL_col < 0)
  {
    XLL_col = 0;
  }
  if (XUL_col >= NCols)
  {
    XUL_col = NCols-1;
  }

  // find the row of old raster that has the same Xlocations as the XLL of smaller raster
  // the 0.5*DataResolution is in case of rounding errors
  // Slightly different logic for y because the DEM starts from the top corner
  int YLL_row = NRows - int((SR_YMinimum-YMinimum+0.5*DataResolution)/DataResolution);
  int YUL_row = NRows - int((SR_YMaximum-YMinimum+0.5*DataResolution)/DataResolution);

  // check these rows
  if (YLL_row < 0)
  {
    YLL_row = 0;
  }
  if (YUL_row >= NRows)
  {
    YUL_row = NRows-1;
  }

  cout << "Small XLLCol: " << XLL_col << " XLR_col: " << XUL_col << " XLLrow: "
       << YLL_row << " YUL_row: " << YUL_row << endl;


  // get the new number of rows and columns:
  int New_NRows = YLL_row-YUL_row;
  int New_NCols = XUL_col-XLL_col;

  cout << "New NRows: " << New_NRows  << " New_NCols: " << New_NCols << endl;

  // now extract the data for the new raster
  float NewR_XMinimum = XMinimum+float(XLL_col)*DataResolution;
  float NewR_YMinimum = YMinimum + ((NRows - YLL_row ) * DataResolution);

  Array2D<float> NewData(New_NRows,New_NCols, NoDataValue);

  //cout << "Writing the array" << endl;

  for(int row = 0; row< New_NRows; row++)
  {
    for(int col = 0; col<New_NCols; col++)
    {
       NewData[row][col] = RasterData[row+YUL_row][col+XLL_col];
    }
  }

  //cout << "Wrote the array" << endl;

  LSDRaster TrimmedRaster(New_NRows, New_NCols, NewR_XMinimum,
                          NewR_YMinimum, DataResolution, NoDataValue, NewData,
                          GeoReferencingStrings);

  cout << "Making the raster" << endl;

  TrimmedRaster.Update_GeoReferencingStrings();

  return TrimmedRaster;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function creates a production raster. 
// It doens't actually use any underlying data in the cosmo raster.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDCosmoRaster::calculate_production_raster(LSDRaster& Elevation_Data,
                                               string path_to_atmospheric_data)
{
  int NRows = Elevation_Data.get_NRows();
  int NCols = Elevation_Data.get_NCols();
  float NDV =  Elevation_Data.get_NoDataValue();
  
  Array2D<float> Production(NRows,NCols,NDV);

  // variables for converting location and elevation
  double this_elevation, this_pressure;

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;
  
  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();
  
  // a function for scaling stone production, defaults to 1
  double Fsp = 1.0;
  
  // the latitude and longitude
  double lat,longitude;
  
  // declare converter object
  LSDCoordinateConverterLLandUTM Converter;
  
  for (int row = 0; row < NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      
      //exclude NDV from average
      if (Elevation_Data.get_data_element(row,col) != NDV)
      {
        // To get pressure, first get the lat and long
        Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
        //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);
      
        // now the elevation
        this_elevation = Elevation_Data.get_data_element(row,col);
      
        // now the pressure
        this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude), 
                                        double(this_elevation));

        // get the production
        Production[row][col] = LSDCRNP.stone2000sp(lat,this_pressure, Fsp);
      }
    }
  }

  float XMinimum = Elevation_Data.get_XMinimum();
  float YMinimum = Elevation_Data.get_YMinimum();
  float DataResolution = Elevation_Data.get_DataResolution();
  map<string,string> GeoReferencingStrings = Elevation_Data.get_GeoReferencingStrings();

  LSDRaster Production_raster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NDV, Production,GeoReferencingStrings);
  return Production_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

LSDRaster LSDCosmoRaster::calculate_atm_pressure_raster(LSDRaster& Elevation_Data, string path_to_atmospheric_data)
{
  int NRows = Elevation_Data.get_NRows();
  int NCols = Elevation_Data.get_NCols();
  float NDV =  Elevation_Data.get_NoDataValue();
  
  Array2D<float> Pressure(NRows,NCols,NDV);

  // variables for converting location and elevation
  double this_elevation, this_pressure;

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;
  
  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();
  
  // a function for scaling stone production, defaults to 1
  double Fsp = 1.0;
  
  // the latitude and longitude
  double lat,longitude;
  
  // declare converter object
  LSDCoordinateConverterLLandUTM Converter;
  
  for (int row = 0; row < NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      
      //exclude NDV from average
      if (Elevation_Data.get_data_element(row,col) != NDV)
      {
        // To get pressure, first get the lat and long
        Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
        //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);
      
        // now the elevation
        this_elevation = Elevation_Data.get_data_element(row,col);
      
        // now the pressure
        Pressure[row][col] = LSDCRNP.NCEPatm_2(double(lat), double(longitude), 
                                        double(this_elevation));

      }
    }
  }

  float XMinimum = Elevation_Data.get_XMinimum();
  float YMinimum = Elevation_Data.get_YMinimum();
  float DataResolution = Elevation_Data.get_DataResolution();
  map<string,string> GeoReferencingStrings = Elevation_Data.get_GeoReferencingStrings();

  LSDRaster Pressure_raster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NDV, Pressure,GeoReferencingStrings);
  return Pressure_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This resets the scaling of the CRN particle which is required since the 
// snow and self shielding alter the scaling values
// The scaling is updated within the object which is passed by reference
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoRaster::reset_scaling(LSDCRNParameters& LSDCRNP, string Muon_scaling)
{
  if (Muon_scaling == "Schaller" )
  {
    LSDCRNP.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    LSDCRNP.set_Braucher_parameters();
  }
  else if (Muon_scaling == "Granger" )
  {
    LSDCRNP.set_Granger_parameters();
  }
  else if (Muon_scaling == "newCRONUS" )
  {
    LSDCRNP.set_newCRONUS_parameters();
  }
  else
  {
    cout << "You didn't set the muon scaling." << endl
         << "Options are Schaller, Braucher, newCRONUS, and Granger." << endl
         << "You chose: " << Muon_scaling << endl
         << "Defaulting to Braucher et al (2009) scaling" << endl;
    LSDCRNP.set_Braucher_parameters();
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This drives the cosmo concentration calculation on the basis of a single
// erosion rate value
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoRaster::predict_mean_CRN_conc(string Nuclide,
                                           string Muon_scaling, 
                                           LSDRaster& eff_erosion_rate,
                                           LSDRaster& ProductionScale, 
                                           LSDRaster& TopoShield, 
                                           LSDRaster& SelfShield,
                                           LSDRaster& SnowShield, 
                                           bool is_production_uncertainty_plus_on,
                                           bool is_production_uncertainty_minus_on)
{
  // the CRN conc in each pixel
  float this_conc;
  
  float this_top_eff_depth;
  float this_bottom_eff_depth;
  
  // set the scaling vector
  vector<bool> nuclide_scaling_switches(4,false);
  if (Nuclide == "Be10")
  {
    nuclide_scaling_switches[0] = true;
  }
  else if (Nuclide == "Al26")
  {
    nuclide_scaling_switches[1] = true;
  }
  else
  {
    cout << "LSDBasin line 1583, You didn't choose a valid nuclide. Defaulting"
         << " to 10Be." << endl;
    Nuclide = "Be10";
    nuclide_scaling_switches[0] = true;
  }

  // the total shielding. A product of snow, topographic and production scaling
  double total_shielding;

  // initiate a particle. We'll just repeatedly call this particle
  // for the sample.
  int startType = 0;
  double Xloc = 0;
  double Yloc = 0;
  double  startdLoc = 0.0;
  double  start_effdloc = 0.0;
  double startzLoc = 0.0;

  // create a particle at zero depth
  LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                               startdLoc, start_effdloc, startzLoc);

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

  // loop through the raster
  for (int col = 0; col < NCols; col++)
  {  
    for (int row = 0; row < NRows; row++)
    {
      //exclude NDV from average
      if (ProductionScale.get_data_element(row,col) == NoDataValue ||
          TopoShield.get_data_element(row,col) == NoDataValue ||
          SelfShield.get_data_element(row,col) == NoDataValue ||
          SnowShield.get_data_element(row,col) == NoDataValue || 
          eff_erosion_rate.get_data_element(row,col) == NoDataValue )
      {
        RasterData[row][col] == NoDataValue;    
      }
      
      // Only calculate if there is a data value at the pixel. 
      if(RasterData[row][col] != NoDataValue)
      { 
        // reset scaling parameters. This is necessary since the F values are
        // reset for local scaling
        reset_scaling(LSDCRNP, Muon_scaling);
  
        // set the scaling to the correct production uncertainty
        vector<double> test_uncert;
        if(is_production_uncertainty_plus_on)
        {
          test_uncert = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
        }
        else if(is_production_uncertainty_minus_on)
        {
          test_uncert = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
        }

        // Get the scaling and shielding
        total_shielding = ProductionScale.get_data_element(row,col)*TopoShield.get_data_element(row,col);

        // scale the F values
        LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

        // Get the top and bottom effective depths from the snow and self shielding rasters
        this_top_eff_depth = SnowShield.get_data_element(row,col);
        this_bottom_eff_depth = this_top_eff_depth+SnowShield.get_data_element(row,col);
 
        // get the nuclide concentration from this node
        if (Nuclide == "Be10")
        {
          //cout << "LInE 2271, 10Be" << endl;
          eroded_particle.update_10Be_SSfull_depth_integrated(double(eff_erosion_rate.get_data_element(row,col)),LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          this_conc=eroded_particle.getConc_10Be();
        }
        else if (Nuclide == "Al26")
        {
          //cout << "LINE 2278, 26Al" << endl;
          eroded_particle.update_26Al_SSfull_depth_integrated(double(eff_erosion_rate.get_data_element(row,col)),LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          this_conc=eroded_particle.getConc_26Al();
        }
        else
        {
          cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
          cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
          cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
          eroded_particle.update_10Be_SSfull_depth_integrated(double(eff_erosion_rate.get_data_element(row,col)),LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          this_conc=eroded_particle.getConc_10Be();
        }
        
        RasterData[row][col] = this_conc;
        
        
      }           // End nodata logic                 
    }             // End cols loop                 
  }               // End rows loop    
}                        

#endif
