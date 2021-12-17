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
//
// Gets the row and column of a point
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::get_row_and_col_of_a_point(float X_coordinate,float Y_coordinate,int& row, int& col)
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

void LSDRasterMaker::get_row_and_col_of_a_point(double X_coordinate,double Y_coordinate,int& row, int& col)
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


// Add a float to all pixels in the raster
void LSDRasterMaker::add_value(float value_to_add)
{
  // now loop through the matrix rescaling the values.
  for (int row = 0; row< NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if(RasterData[row][col] != NoDataValue)
      {
        RasterData[row][col] = RasterData[row][col]+value_to_add;
      }
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function add strips of a given value.
// This happily overwrites NoData
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterMaker::add_strip(int start_row_or_col, int end_row_or_col, bool horizontal, float value)
{
  if (start_row_or_col < 0)
  {
    start_row_or_col = 0;
  }

  if (horizontal)
  {
    if(start_row_or_col >= NRows)
    {
      start_row_or_col = NRows-1;
    }
    if(end_row_or_col > NRows)
    {
      end_row_or_col = NRows-1;
    }
  }
  else
  {
    if(start_row_or_col >= NCols)
    {
      start_row_or_col = NCols-1;
    }
    if(end_row_or_col > NCols)
    {
      end_row_or_col = NCols-1;
    }
  }

  if (horizontal)
  {
    for (int row = start_row_or_col; row <= end_row_or_col; row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        RasterData[row][col] = value;
      }
    }
  }
  else
  {
     for (int row = 0; row < NRows; row++)
    {
      for (int col = start_row_or_col; end_row_or_col<= NCols; col++)
      {
        RasterData[row][col] = value;
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
void LSDRasterMaker::impose_channels(LSDSpatialCSVReader& source_points_data, string column_name)
{

  // string column_name = "elevation(m)";


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // Get the local coordinate as well as the elevations
  vector<float> UTME, UTMN;
  source_points_data.get_x_and_y_from_latlong(UTME,UTMN);
  vector<float> elev = source_points_data.data_column_to_float(column_name);


  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  int row,col;
  for(int i = 0; i< int(elev.size()); i++)
  {

    source_points_data.get_row_and_col_of_a_point(UTME[i],UTMN[i],row, col);
    zeta[row][col] = elev[i];
  }


  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}


////------------------------------------------------------------------------------
//// This takes a raster an increases the elevation of nodata pixels at the edge nodes
//// so that the outside pixels always drain inwards
////------------------------------------------------------------------------------
void LSDRasterMaker::buffer_basin_to_single_outlet(LSDSpatialCSVReader& source_points_data, float slope)
{
  
  string column_name = "elevation(m)";

  float elev_diff = 10;  //DataResolution*sqrt(2)*slope;

  vector<int> min_elev_nd_rows;
  vector<int> min_elev_nd_cols;

  Array2D<float> zeta=RasterData.copy();
  Array2D<float> old_zeta = RasterData.copy();

  // Get the local coordinate as well as the elevations
  vector<float> UTME, UTMN;
  source_points_data.get_x_and_y_from_latlong(UTME,UTMN);

  int min_row,min_col;
  int rp1,rm1,cp1,cm1;

  vector<float> elev;
  if ( UTME.size() == 1)
  {
    cout << "Only one element!" << endl;
    elev.push_back(0);
  }
  else
  {
    elev = source_points_data.data_column_to_float(column_name);
    cout << "There are " << elev.size() << " data points" << endl;
  }

  // find the minimum elevation
  int n_elev = int(elev.size());
  float min_elev = 10000000000;
  int node_of_min_elev = 0;
  cout << "n_elev: " << n_elev << endl;
  for(int i = 0; i< n_elev; i++)
  {
    cout << "elev: " << elev[i] << endl;
    if(elev[i] < min_elev)
    {
      min_elev = elev[i];
      node_of_min_elev = i;
    }
  }
  get_row_and_col_of_a_point(UTME[node_of_min_elev],UTMN[node_of_min_elev],min_row, min_col);
  cout << "The minimum elevation of the source point is: " << min_elev << " at node: " << node_of_min_elev << endl;
  cout << "minimum elev in source points is at " << UTME[node_of_min_elev] << " , " << UTMN[node_of_min_elev] << endl;

  // now logic for finding the nodata around the minimum elevation
  rp1 = min_row+1;
  if(rp1 == NRows)
  {
    rp1 = min_row;
  }
  rm1 = min_row-1;
  if (rm1 == -1)
  {
    rm1 = 0;
  }
  cp1 = min_col+1;
  if (cp1 == NCols)
  {
    cp1 = min_col;
  }
  cm1 = min_col-1;
  if (cm1 == -1)
  {
    cm1 = 0;
  }   
 
  // now search all adjacent nodes
  if ( old_zeta[rp1][cp1] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(rp1);
    min_elev_nd_cols.push_back(cp1);
  }
  if ( old_zeta[rp1][min_col] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(rp1);
    min_elev_nd_cols.push_back(min_col);
  }
  if ( old_zeta[rp1][cm1] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(rp1);
    min_elev_nd_cols.push_back(cm1);
  }
  if ( old_zeta[min_row][cp1] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(min_row);
    min_elev_nd_cols.push_back(cp1);
  }
  if ( old_zeta[min_row][cm1] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(min_row);
    min_elev_nd_cols.push_back(cm1);
  }
  if ( old_zeta[rm1][cp1] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(rm1);
    min_elev_nd_cols.push_back(cm1);
  }
  if ( old_zeta[rm1][min_col] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(rm1);
    min_elev_nd_cols.push_back(min_col);
  }
  if ( old_zeta[rm1][cm1] == NoDataValue)
  {
    cout << "Buffering, found ndv" << endl;
    min_elev_nd_rows.push_back(rm1);
    min_elev_nd_cols.push_back(cm1);
  }     

  // Now crinkle up the side
  for(int row = 0; row< NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      if (row == min_row && col == min_col)
      {
        cout << "Found the minimum elevation" << endl;
      }
      else
      {
        if (old_zeta[row][col] != NoDataValue)
        {
          // now search neighbouring nodes for nodata
          rp1 = row+1;
          if(rp1 == NRows)
          {
            rp1 = row;
          }
          rm1 = row-1;
          if (rm1 == -1)
          {
            rm1 = 0;
          }
          cp1 = col+1;
          if (cp1 == NCols)
          {
            cp1 = col;
          }
          cm1 = col-1;
          if (cm1 == -1)
          {
            cm1 = 0;
          }

          // now search all adjacent nodes
          if ( old_zeta[rp1][cp1] == NoDataValue)
          {
            zeta[rp1][cp1] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[rp1][col] == NoDataValue)
          {
            zeta[rp1][col] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[rp1][cm1] == NoDataValue)
          {
            zeta[rp1][cm1] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[row][cp1] == NoDataValue)
          {
            zeta[row][cp1] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[row][cm1] == NoDataValue)
          {
            zeta[row][cm1] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[rm1][cp1] == NoDataValue)
          {
            zeta[rm1][cp1] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[rm1][col] == NoDataValue)
          {
            zeta[rm1][col] = old_zeta[row][col]+elev_diff;
          }
          if ( old_zeta[rm1][cm1] == NoDataValue)
          {
            zeta[rm1][cm1] = old_zeta[row][col]+elev_diff;
          }
        }
      }
    }
  }

  // now return the data around the minimum value to nodata
  for (int i = 0; i< int(min_elev_nd_rows.size()); i++)
  {
    //cout << "yoyoma" << endl;
    zeta[ min_elev_nd_rows[i] ][ min_elev_nd_cols[i] ] = NoDataValue;
  }


  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}



void LSDRasterMaker::buffer_basin_to_single_outlet(float slope)
{
  Array2D<float> zeta=RasterData.copy();  

  vector<int> edge_rows;
  vector<int> edge_cols;

  float elev_diff = DataResolution*sqrt(2)*slope;

  // we look for any edge nodes. We don't look along the edge of the DEM
  // any node that is not itself nodata but has an edge pixel as nodata
  // gets added as an edge node
  int rp1,rm1,cp1,cm1;
  bool is_edge_node;
  int min_edge_row = 0;
  int min_edge_col = 0;
  float min_edge_elev = 99999999999999;
  for(int row = 1; row < NRows-1; row++)
  {
    for(int col = 1; col < NCols-1; col++)
    {
      if (zeta[row][col] != NoDataValue)
      {
        is_edge_node = false;
        rp1 = row+1;
        rm1 = row-1;
        cp1 = col+1;
        cm1 = col-1;

        // now search all adjacent nodes
        if ( zeta[rp1][cp1] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[rp1][col] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[rp1][cm1] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[row][cp1] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[row][cm1] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[rm1][cp1] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[rm1][col] == NoDataValue)
        {
          is_edge_node = true;
        }
        if ( zeta[rm1][cm1] == NoDataValue)
        {
          is_edge_node = true;
        }

        if (is_edge_node)
        {
          edge_rows.push_back(row);
          edge_cols.push_back(col);

          if (min_edge_elev > zeta[row][col])
          {
            min_edge_elev = zeta[row][col];
            min_edge_row = row;
            min_edge_col = col;
          }
        }
      }
    }
  }


  // Okay, now look through the edge rows and columns, buffering up the nodes.
  int n_edge_nodes = int(edge_rows.size()); 
  int row,col;
  for(int i = 0; i< n_edge_nodes; i++)
  {
    row = edge_rows[i];
    col = edge_cols[i];
   
    rp1 = row+1;
    rm1 = row-1;
    cp1 = col+1;
    cm1 = col-1;
    
    if (row ==  min_edge_row && col ==  min_edge_col)
    {
      cout << "Found the minimum elevation edge node." << endl;
    }
    else
    { 
      // now search all adjacent nodes
      if ( zeta[rp1][cp1] == NoDataValue)
      {
        zeta[rp1][cp1] = zeta[row][col]+elev_diff;
      }
      if ( zeta[rp1][col] == NoDataValue)
      {
        zeta[rp1][col] = zeta[row][col]+elev_diff;
      }
      if ( zeta[rp1][cm1] == NoDataValue)
      {
        zeta[rp1][cm1] = zeta[row][col]+elev_diff;
      }
      if ( zeta[row][cp1] == NoDataValue)
      {
        zeta[row][cp1] = zeta[row][col]+elev_diff;
      }
      if ( zeta[row][cm1] == NoDataValue)
      {
        zeta[row][cm1] = zeta[row][col]+elev_diff;
      }
      if ( zeta[rm1][cp1] == NoDataValue)
      {
        zeta[rm1][cp1] = zeta[row][col]+elev_diff;
      }
      if ( zeta[rm1][col] == NoDataValue)
      {
        zeta[rm1][col] = zeta[row][col]+elev_diff;
      }
      if ( zeta[rm1][cm1] == NoDataValue)
      {
        zeta[rm1][cm1] = zeta[row][col]+elev_diff;
      }
    }
  }

  this->RasterData = zeta.copy();

  RasterData = zeta.copy(); 
}

////------------------------------------------------------------------------------
//// impose_channels: this imposes channels onto the landscape
//// You need to print a channel to csv and then load the data
////------------------------------------------------------------------------------
void LSDRasterMaker::impose_channels_with_buffer(LSDSpatialCSVReader& source_points_data, float slope, string column_name)
{

  // string column_name = "elevation(m)";

  float elev_diff = DataResolution*sqrt(2)*slope;


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // Get the local coordinate as well as the elevations
  vector<float> UTME, UTMN;
  source_points_data.get_x_and_y_from_latlong(UTME,UTMN);
  vector<float> elev = source_points_data.data_column_to_float(column_name);

  // find the minimum elevation
  float min_elev = 10000000000;
  int node_of_min_elev = 0;
  for(int i = 0; i< int(elev.size()); i++)
  {
    if(elev[i] < min_elev)
    {
      min_elev = elev[i];
      node_of_min_elev = i;
    }
  }

  cout << "The minimum elevation is: " << min_elev << endl;


  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  int row,col;
  int rp1,rm1,cp1,cm1;
  for(int i = 0; i< int(elev.size()); i++)
  {

    source_points_data.get_row_and_col_of_a_point(UTME[i],UTMN[i],row, col);
    zeta[row][col] = elev[i];

    if (i != node_of_min_elev)
    {
      // now search neighbouring nodes for nodata
      rp1 = row+1;
      if(rp1 == NRows)
      {
        rp1 = row;
      }
      rm1 = row-1;
      if (rm1 == -1)
      {
        rm1 = 0;
      }
      cp1 = col+1;
      if (cp1 == NCols)
      {
        cp1 = col;
      }
      cm1 = col-1;
      if (cm1 == -1)
      {
        cm1 = 0;
      }

      // now search all adjacent nodes
      if ( zeta[rp1][cp1] == NoDataValue)
      {
        zeta[rp1][cp1] = elev[i]+elev_diff;
      }
      if ( zeta[rp1][col] == NoDataValue)
      {
        zeta[rp1][col] = elev[i]+elev_diff;
      }
      if ( zeta[rp1][cm1] == NoDataValue)
      {
        zeta[rp1][cm1] = elev[i]+elev_diff;
      }
      if ( zeta[row][cp1] == NoDataValue)
      {
        zeta[row][cp1] = elev[i]+elev_diff;
      }
      if ( zeta[row][cm1] == NoDataValue)
      {
        zeta[row][cm1] = elev[i]+elev_diff;
      }
      if ( zeta[rm1][cp1] == NoDataValue)
      {
        zeta[rm1][cp1] = elev[i]+elev_diff;
      }
      if ( zeta[rm1][col] == NoDataValue)
      {
        zeta[rm1][col] = elev[i]+elev_diff;
      }
      if ( zeta[rm1][cm1] == NoDataValue)
      {
        zeta[rm1][cm1] = elev[i]+elev_diff;
      }
    }
  }


  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}

////------------------------------------------------------------------------------
//// impose_channels: this imposes channels onto the landscape using XY data
//// You need to print a channel to csv and then load the data
//// ELSG 23/02/2021
////------------------------------------------------------------------------------
void LSDRasterMaker::impose_channels_with_buffer_use_XY(LSDSpatialCSVReader& source_points_data, float slope, string column_name)
{

  // string column_name = "elevation(m)";
  string x_column_name = "X";
  string y_column_name = "Y";

  float elev_diff = DataResolution*sqrt(2)*slope;


  Array2D<float> zeta=RasterData.copy();

  // Step one, create donor "stack" etc. via FlowInfo
  LSDRaster temp(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, zeta, GeoReferencingStrings);

  // Get the local coordinate as well as the elevations
  vector<float> X = source_points_data.data_column_to_float(x_column_name);
  vector<float> Y = source_points_data.data_column_to_float(y_column_name);
  //source_points_data.get_x_and_y_from_latlong(X,Y);
  vector<float> elev = source_points_data.data_column_to_float(column_name);

  // find the minimum elevation
  float min_elev = 10000000000;
  int node_of_min_elev = 0;
  for(int i = 0; i< int(elev.size()); i++)
  {
    if(elev[i] < min_elev)
    {
      min_elev = elev[i];
      node_of_min_elev = i;
    }
  }

  cout << "The minimum elevation is: " << min_elev << endl;


  // make the map
  cout << "I am making an elevation map. This will not work if points in the raster lie outside of the csv channel points." << endl;
  int row,col;
  int rp1,rm1,cp1,cm1;
  for(int i = 0; i< int(elev.size()); i++)
  {
    cout << "At i " << i << ", X is: " << X[i] << ", and Y is: " << Y[i] << "; elevation is: " << elev[i] << endl;
    source_points_data.get_row_and_col_of_a_point(X[i],Y[i],row, col);
    cout << "Row is: " << row << ", and col is: " << col << endl;
    zeta[row][col] = elev[i];

    if (i != node_of_min_elev)
    {
      // now search neighbouring nodes for nodata
      rp1 = row+1;
      if(rp1 == NRows)
      {
        rp1 = row;
      }
      rm1 = row-1;
      if (rm1 == -1)
      {
        rm1 = 0;
      }
      cp1 = col+1;
      if (cp1 == NCols)
      {
        cp1 = col;
      }
      cm1 = col-1;
      if (cm1 == -1)
      {
        cm1 = 0;
      }

      // now search all adjacent nodes
      if ( zeta[rp1][cp1] == NoDataValue)
      {
        zeta[rp1][cp1] = elev[i]+elev_diff;
      }
      if ( zeta[rp1][col] == NoDataValue)
      {
        zeta[rp1][col] = elev[i]+elev_diff;
      }
      if ( zeta[rp1][cm1] == NoDataValue)
      {
        zeta[rp1][cm1] = elev[i]+elev_diff;
      }
      if ( zeta[row][cp1] == NoDataValue)
      {
        zeta[row][cp1] = elev[i]+elev_diff;
      }
      if ( zeta[row][cm1] == NoDataValue)
      {
        zeta[row][cm1] = elev[i]+elev_diff;
      }
      if ( zeta[rm1][cp1] == NoDataValue)
      {
        zeta[rm1][cp1] = elev[i]+elev_diff;
      }
      if ( zeta[rm1][col] == NoDataValue)
      {
        zeta[rm1][col] = elev[i]+elev_diff;
      }
      if ( zeta[rm1][cm1] == NoDataValue)
      {
        zeta[rm1][cm1] = elev[i]+elev_diff;
      }
    }
  }


  this->RasterData = zeta.copy();

  RasterData = zeta.copy();
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cap elevations using an initial raster
//------------------------------------------------------------------------------
void LSDRasterMaker::cap_elevations(LSDRaster& InitialRaster)
{
  cout << "Capping elevations. WARNING: no checking rasters are same dimension" << endl;
  for(int row=0; row<NRows; ++row)
  {
    for(int col=0; col<NCols; ++col)
    {
      if(RasterData[row][col]!=NoDataValue)
      {
        if (RasterData[row][col] > InitialRaster.get_data_element(row,col))
        {
          RasterData[row][col] = InitialRaster.get_data_element(row,col);
        }
      }
    }
  }  
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
