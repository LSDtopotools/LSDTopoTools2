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
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDRasterMaker.hpp"
#include "LSDStatsTools.hpp"
#include "LSDSpatialCSVReader.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterMaker_CPP
#define LSDRasterMaker_CPP




void LSDRasterMaker::create()
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
void LSDRasterMaker::create(string filename, string extension)
{
  read_raster(filename,extension);
}



void LSDRasterMaker::create(int NRows, int NCols)
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
void LSDRasterMaker::create(LSDRaster& An_LSDRaster)
{
  cout << "Lets get some info!" << endl;
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
LSDRaster LSDRasterMaker::return_as_raster()
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
void LSDRasterMaker::resize_and_reset( int new_rows, int new_cols, float new_resolution, float new_value )
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
vector<float> LSDRasterMaker::minimum_and_maximum_value()
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
// This function resets all non nodata nodes to a constant value
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::set_to_constant_value(float new_value)
{
  // now loop through the matrix rescaling the values.
  for (int row = 0; row< NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if(RasterData[row][col] != NoDataValue)
      {
        RasterData[row][col] = new_value;
      }
    }
  }
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes the existing raster data and then linearly scales it
// to new minimum and maximum values.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::scale_to_new_minimum_and_maximum_value(float new_minimum, float new_maximum)
{
  // first we get the existing minimum and maximum
  vector<float> min_max = minimum_and_maximum_value();

  float scaling_fraction;
  float original_range = min_max[1] - min_max[0];
  float new_range = new_maximum-new_minimum;

  // now loop through the matrix rescaling the values.
  for (int row = 0; row< NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if(RasterData[row][col] != NoDataValue)
      {
        // first find where the value is between min and max
        if(original_range == 0)
        {
          scaling_fraction = 0;
        }
        else
        {
          scaling_fraction = (RasterData[row][col] - min_max[0])/original_range;
        }
        RasterData[row][col] = (scaling_fraction*new_range + new_minimum);
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


////------------------------------------------------------------------------------
//// impose_channels: this imposes channels onto the landscape
//// You need to print a channel to csv and then load the data
////------------------------------------------------------------------------------
void LSDRasterMaker::impose_channels(LSDSpatialCSVReader& source_points_data)
{

  string column_name = "elevation(m)";


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // need to fill the raster to ensure there are no internal base level nodes
  cout << "I am going to fill" << endl;
  float slope_for_fill = 0.0001; 
  cout << "Filling." << endl;
  LSDRaster filled_topography = temp.fill(slope_for_fill);

  vector <string> bc(4, "b");          // Initialise boundaries to baselevel


  cout << "Getting the flow info. This might take some time." << endl;
  LSDFlowInfo flow(bc, filled_topography);
  // update the raster
  zeta = filled_topography.get_RasterData();

  // Get the local node index as well as the elevations
  vector<int> ni = source_points_data.get_nodeindices_from_lat_long(flow);
  vector<float> elev = source_points_data.data_column_to_float(column_name);
  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  int row,col;
  for(int i = 0; i< int(ni.size()); i++)
  {
    flow.retrieve_current_row_and_col( ni[i], row, col);
    zeta[row][col] = elev[i];
  }


  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This smooths the raster by taking a weighted average of the given pixel
// and neighboring pixels.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::smooth(int boundary_type)
{
  // at the moment the boundary type can only be 0 and this is a periodic
  // boundary type at the E and W boundaries.

  Array2D<float> new_data(NRows,NCols,NoDataValue);
  float total_weighting;
  float total_sum;
  int rp1, rm1,cp1, cm1;


  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      total_weighting = 0;
      total_sum = 0;

      rp1 = row+1;
      rm1 = row-1;
      cp1 = col+1;
      cm1 = col-1;

      // implement boundary conditions.
      if(boundary_type == 0)
      {
        if (rp1 == NRows)
        {
          rp1 = rm1;
        }
        if (rm1 == -1)
        {
          rm1 = rp1;
        }
        if (cp1 == NCols)
        {
          cp1 = 0;
        }
        if(cm1 == -1)
        {
          cm1 = NCols-1;
        }
      }
      else
      {
        if (rp1 == NRows)
        {
          rp1 = rm1;
        }
        if (rm1 == -1)
        {
          rm1 = rp1;
        }
        if (cp1 == NCols)
        {
          cp1 = 0;
        }
        if(cm1 == -1)
        {
          cm1 = NCols-1;
        }
      }

      if( RasterData[row][col] != NoDataValue)
      {
        total_weighting += 2;
        total_sum += 2*RasterData[row][col];

        // now go through all the other directions.
        if (RasterData[row][cp1] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[row][cp1];
        }
        if (RasterData[row][cm1] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[row][cm1];
        }
        if (RasterData[rp1][col] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[rp1][col];
        }
        if (RasterData[rm1][col] != NoDataValue)
        {
          total_weighting +=1;
          total_sum += RasterData[rm1][col];
        }
      }
      // Now update the array

      new_data[row][col] = total_sum/total_weighting;
    }
  }

  RasterData = new_data.copy();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Some functions for making random values in the rasters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::random_values(float minimum_value, float maximum_value)
{
  long seed = time(NULL);

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        RasterData[row][col] = ran3(&seed);  
      }
    }
  }

  // Now scale to min and max
  scale_to_new_minimum_and_maximum_value(minimum_value, minimum_value);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Some functions for making random values in the rasters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::random_square_blobs(int minimum_blob_size, int maximum_blob_size, float minimum_value, float maximum_value, int n_blobs)
{
  long seed = time(NULL);

  // Lets make some blobs!!!
  for (int blob_n = 0; blob_n < n_blobs; blob_n++)
  {
    // get the centrepoint of the blob
    float row_frac = ran3(&seed);
    float col_frac = ran3(&seed);

    // get the row and column of the centre of the blob
    float frow = row_frac*float(NRows);
    int this_row = floor(frow);

    if (this_row < 0)
    {
      this_row = 0;
    }
    if(this_row >= NRows)
    {
      this_row = NRows-1;
    }

    // get the row and column of the centre of the blob
    float fcol = col_frac*float(NCols);
    int this_col = floor(fcol);

    if (this_col < 0)
    {
      this_col = 0;
    }
    if(this_col >= NCols)
    {
      this_col = NCols-1;
    }

    // Get the size of the blob. This will need to be odd
    int size_range =  maximum_blob_size-minimum_blob_size;
    float this_size;

    if(size_range == 0)
    {
      cout << "Check you prarameters, the size range is zero." << endl;
      this_size = minimum_blob_size;
    }
    else
    {
      this_size = ran3(&seed)*float(size_range)+float(minimum_blob_size);
    }

    int size = int(this_size);

    // get the starting rows and ending rows. Note that I am not being very careful about
    // this being exactly the right dimension
    int start_row = this_row - size/2;
    int end_row = this_row + size/2;
    int start_col = this_col - size/2;
    int end_col = this_col + size/2;

    // get the new value
    float value_range = maximum_value-minimum_value;
    float this_blob_value;
    if(value_range == 0)
    {
      cout << "Check you prarameters, the value range is zero." << endl;
      this_blob_value = minimum_value;
    }
    else
    {
      this_blob_value = ran3(&seed)*value_range+minimum_value;
    }

    //cout << "This blob is: " << blob_n << " with a K of: " << this_blob_value << endl;
    //cout << "The size is: " << this_size << " or " << size <<endl;

    // now update the values
    for(int row = start_row; row<=end_row; row++)
    {
      for(int col = start_col; col<= end_col; col++)
      {

        if( row >= 0 && row<NRows && col >= 0 && col<NCols)
        {
          RasterData[row][col] = this_blob_value;
        }
      }
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Some functions for making random values in the rasters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::sine_waves(vector<float> x_coefficients, vector<float> y_coefficients)
{
  int n_x_coeff = int(x_coefficients.size());
  int n_y_coeff = int(y_coefficients.size());

  float x_factor = M_PI/ float(NCols-1);
  float y_factor = M_PI / float(NRows-1);

  float this_x_value;
  float this_y_value;

  // so the wavelengths of the sin waves depend on the number
  for (int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      this_x_value = 0;
      for(int xv = 0; xv<n_x_coeff; xv++)
      {
        this_x_value+=x_coefficients[xv]*sin(x_factor*(xv+1)*float(col));
      }
      this_y_value = 0;
      for(int yv = 0; yv<n_y_coeff; yv++)
      {
        this_y_value+=y_coefficients[yv]*sin(y_factor*(yv+1)*float(row));
      }
      RasterData[row][col] = this_x_value+this_y_value;
    }
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis clips to a smaller raster. The smaller raster does not need
// to have the same data resolution as the old raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDRasterMaker::clip_to_smaller_raster(LSDRaster& smaller_raster)
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
vector<float> LSDRasterMaker::get_XY_MinMax()
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
vector<float> LSDRasterMaker::get_RasterData_vector()
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
void LSDRasterMaker::write_RasterData_to_text_file(string filename)
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
// THis clips to a smaller raster. The smaller raster does not need
// to have the same data resolution as the old raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDRasterMaker::clip_to_smaller_raster(LSDIndexRaster& smaller_raster)
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
//--------------------------------------------------------------------------//
// Function to uplift the north half of the raster at a defined rate.
// e.g. to uplift the north half an order of magnitude faster than the original
// value, set increase_amt to 10.
// FJC 28/06/18
//--------------------------------------------------------------------------//
void LSDRasterMaker::increase_north_half_raster_values(int increase_amt)
{
  // loop through the north half of the raster and change the values.
  for (int i = 0; i < int(NRows/2); i++)
  {
    for (int j = 0; j < NCols; j++)
    {
      RasterData[i][j] = RasterData[i][j] * increase_amt;
    }
  }
}

//--------------------------------------------------------------------------//
// Function to uplift the south half of the raster at a defined rate.
// e.g. to uplift the north half an order of magnitude faster than the original
// value, set increase_amt to 10.
// FJC 28/06/18
//--------------------------------------------------------------------------//
void LSDRasterMaker::increase_south_half_raster_values(int increase_amt)
{
  // loop through the south half of the raster and change the values.
  for (int i = int(NRows/2); i < NRows; i++)
  {
    for (int j = 0; j < NCols; j++)
    {
      RasterData[i][j] = RasterData[i][j] * increase_amt;
    }
  }
}

void LSDRasterMaker::increase_south_quarter_raster_values(int increase_amt)
{
  // loop through the south half of the raster and change the values.
  for (int i = int(NRows - NRows/4); i < NRows; i++)
  {
    for (int j = 0; j < NCols; j++)
    {
      RasterData[i][j] = RasterData[i][j] * increase_amt;
    }
  }
}

void LSDRasterMaker::increase_west_half_raster_values(int increase_amt)
{
  // loop through the south half of the raster and change the values.
  for (int i = 0; i < NRows; i++)
  {
    for (int j = 0; j < int(NCols/2); j++)
    {
      RasterData[i][j] = RasterData[i][j] * increase_amt;
    }
  }
}

//----------------------------------------------------------------------------//
// Tilted uplift field
// for progressive block tilting over time
//----------------------------------------------------------------------------//
void LSDRasterMaker::tilted_block(float angle, string tilt_boundary)
{
  Array2D<float> old_values = RasterData.copy();
  Array2D<float> new_values(NRows, NCols, NoDataValue);

  // loop through each possible tilt direction and get the new array of values
  // after tilting
  if (tilt_boundary == "N")  // north is tilt boundary - max at the S
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // first row stays at same elevation
        if (i == 0) { new_values[i][j] = old_values[i][j]; }
        // other rows, calculate based on angle
        else
        {
          float this_value = old_values[i][j];
          float length = (i + 1) * DataResolution;
          new_values[i][j] = (length * tan(angle)) + this_value;
          //cout << "old elev: " << this_value << " new elev: " << (length * tan(angle)) + this_value << endl;
        }
      }
    }
  }
  else if (tilt_boundary == "S")  // south is tilt boundary - max elevation at the N
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // last row stays at same elevation
        if (i == NRows-1) { new_values[i][j] = old_values[i][j]; }
        // other rows, calculate based on angle
        else
        {
          float this_value = old_values[i][j];
          float length = (NRows - i) * DataResolution;
          new_values[i][j] = (length * tan(angle)) + this_value;
          //cout << "old elev: " << this_value << " new elev: " << (length * tan(angle)) + this_value << endl;
        }
      }
    }
  }
  else if (tilt_boundary == "E")  // east is tilt boundary - max elevation at the W
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // last col stays at same elevation
        if (j == NCols-1) { new_values[i][j] = old_values[i][j]; }
        // other cols, calculate based on angle
        else
        {
          float this_value = old_values[i][j];
          float length = (NCols - i) * DataResolution;
          new_values[i][j] = (length * tan(angle)) + this_value;
        //  cout << "old value: " << this_value << " new elev: " << (length * tan(angle)) + this_value << endl;
        }
      }
    }
  }
  else if (tilt_boundary == "W")  // west is tilt boundary - max elevation at the E
  {
    for (int i = 0; i < NRows; i++)
    {
      for (int j = 0; j < NCols; j++)
      {
        // first col stays at same elevation
        if (j == 0) { new_values[i][j] = old_values[i][j]; }
        // other cols, calculate based on angle
        else
        {
          float this_value = old_values[i][j];
          float length = (j + 1) * DataResolution;
          new_values[i][j] = (length * tan(angle)) + this_value;
          //cout << "old elev: " << this_value << " new elev: " << (length * tan(angle)) + this_value << endl;
        }
      }
    }
  }
  else
  {
    cout << "Warning - you haven't set your boundary to N, W, E, or S. Returning the original raster" << endl;
    new_values = old_values;
  }

  // set the model to the array of new elevations
  RasterData = new_values;
}

#endif
