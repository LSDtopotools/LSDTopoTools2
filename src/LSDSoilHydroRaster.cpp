//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDSoilHydroRaster
// Land Surface Dynamics Raster for manipulating soil and hydroligcal data
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing raster data, with a particular focus on topography
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDSoilHydroRaster.cpp
// cpp file for the LSDSoilHydroRaster object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David T. Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
// Stuart Grieve, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 1.0.0		12/11/2015
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDSoilHydroRaster_CPP
#define LSDSoilHydroRaster_CPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include <cmath>
#include "LSDRaster.hpp"
#include "LSDSoilHydroRaster.hpp"
#include "LSDStatsTools.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDShapeTools.hpp"
#include "LSDFlowInfo.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// empty default create function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::create()
{}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Create function that just copies a raster into the hydro raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::create(LSDRaster& OtherRaster)
{
  NRows = OtherRaster.get_NRows();
  NCols = OtherRaster.get_NCols();
  XMinimum = OtherRaster.get_XMinimum();
  YMinimum = OtherRaster.get_YMinimum();
  DataResolution = OtherRaster.get_DataResolution();
  NoDataValue = OtherRaster.get_NoDataValue();
  GeoReferencingStrings = OtherRaster.get_GeoReferencingStrings();
  RasterData = OtherRaster.get_RasterData();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Create function that takes the dimensions and georeferencing of a raster
// but then sets all data to value, setting the NoDataValues to
// the NoData of the raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::create(LSDRaster& OtherRaster, float value)
{
  NRows = OtherRaster.get_NRows();
  NCols = OtherRaster.get_NCols();
  XMinimum = OtherRaster.get_XMinimum();
  YMinimum = OtherRaster.get_YMinimum();
  DataResolution = OtherRaster.get_DataResolution();
  NoDataValue = OtherRaster.get_NoDataValue();
  GeoReferencingStrings = OtherRaster.get_GeoReferencingStrings();

  // set the raster data to be a certain value
  Array2D<float> data(NRows,NCols,NoDataValue);

  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {

      if (OtherRaster.get_data_element(row,col) != NoDataValue)
      {
        data[row][col] = value;
        //cout << value << endl;
      }
    }
  }

  RasterData = data.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDSoilHydroRaster::create(LSDRaster& DEM, LSDRaster& OtherRaster, int min_max)
{
  NRows = OtherRaster.get_NRows();
  NCols = OtherRaster.get_NCols();
  XMinimum = OtherRaster.get_XMinimum();
  YMinimum = OtherRaster.get_YMinimum();
  DataResolution = OtherRaster.get_DataResolution();
  NoDataValue = OtherRaster.get_NoDataValue();
  GeoReferencingStrings = OtherRaster.get_GeoReferencingStrings();

  // set the raster data to be a certain value
  Array2D<float> data(NRows,NCols,NoDataValue);

  float min_max_val = NoDataValue;

  if (min_max == 0){
    // get the minimum value of OtherRaster
    Array2D<float> tmp = OtherRaster.get_RasterData();
    min_max_val = Get_Minimum(tmp, NoDataValue);
  }
  else if (min_max == 1){
    // get the maximum value of OtherRaster
    Array2D<float> tmp = OtherRaster.get_RasterData();
    min_max_val = Get_Maximum(tmp, NoDataValue);
  }
  else if (min_max == 2){
    // Fill with NDVs
    min_max_val = NoDataValue;
  }

  // for each cell, if there is no paramter data but there is topo, fill in the data with the minimum/maximum value
  // otherwise, just keep the minimum value.
  for (int i = 0; i < NRows; ++i){
    for (int j = 0; j < NCols; ++j){
      if (DEM.get_data_element(i, j) != NoDataValue && OtherRaster.get_data_element(i,j) == NoDataValue){
        data[i][j] = min_max_val;
      }
      else if (DEM.get_data_element(i, j) != NoDataValue && OtherRaster.get_data_element(i, j) != NoDataValue){
        data[i][j] = OtherRaster.get_data_element(i,j);
      }
    }
  }

  RasterData = data.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates a raster from raw data
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::create(int nrows, int ncols, float xmin, float ymin,
                float cellsize, float ndv, Array2D<float> data)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;

  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NCols!" << endl;
    exit(EXIT_FAILURE);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates a raster from raw data, this time with the georeferencing strings
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::create(int nrows, int ncols, float xmin, float ymin,
            float cellsize, float ndv, Array2D<float> data, map<string,string> GRS)
{
  NRows = nrows;
  NCols = ncols;
  XMinimum = xmin;
  YMinimum = ymin;
  DataResolution = cellsize;
  NoDataValue = ndv;
  GeoReferencingStrings = GRS;

  RasterData = data.copy();

  if (RasterData.dim1() != NRows)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NRows!" << endl;
    exit(EXIT_FAILURE);
  }
  if (RasterData.dim2() != NCols)
  {
    cout << "LSDSoilHydroRaster dimension of data is not the same as stated in NCols!" << endl;
    exit(EXIT_FAILURE);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets all non nodata pixels to value
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::SetHomogenousValues(float value)
{
  // set the raster data to be a certain value
  Array2D<float> data(NRows,NCols,NoDataValue);

  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        RasterData[row][col] = value;
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates a snow thickenss (effective, in g cm^-2 for cosmogenic
// applications)  based on a bilinear model such as that of P Kirchner: http://escholarship.org/uc/item/9zn1c1mk#page-8
// The paper is here: http://www.hydrol-earth-syst-sci.net/18/4261/2014/hess-18-4261-2014.html
// This paper also agrees withy this general trend:
// http://www.the-cryosphere.net/8/2381/2014/tc-8-2381-2014.pdf
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::SetSnowEffDepthBilinear(float SlopeAscend, float SlopeDescend,
                          float PeakElevation, float PeakSnowpack, LSDRaster& Elevation)
{
  float LocalElevation;
  float ascendEffDepth;
  float descendEffDepth;
  float thisEffDepth = 0;

  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        LocalElevation = Elevation.get_data_element(row,col);

        if (LocalElevation != NoDataValue)
        {
          // get the effective depth on both the ascending and descending limb
          ascendEffDepth = SlopeAscend*(LocalElevation-PeakElevation)+PeakSnowpack;
          descendEffDepth = SlopeDescend*(LocalElevation-PeakElevation)+PeakSnowpack;

          // the correct depth is the lesser of the two
          if (ascendEffDepth < descendEffDepth)
          {
            thisEffDepth =  ascendEffDepth;
          }
          else
          {
            thisEffDepth = descendEffDepth;
          }

          // if the depth is less than zero, then set to zero
          if(thisEffDepth <0)
          {
            thisEffDepth = 0;
          }

          RasterData[row][col] =  thisEffDepth;

        }
        else        // if there ins't any elevation data, set the snow data to NoData
        {
          RasterData[row][col] = NoDataValue;
        }
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates a snow thickenss (effective, in g cm^-2 for cosmogenic
// applications)  based on a richard's equsion sigmoidal growth model
// It was propoesd to represent peak SWE so we cruedly apply it to average annual SWE
// see
// http://onlinelibrary.wiley.com/doi/10.1002/2015GL063413/epdf
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::SetSnowEffDepthRichards(float MaximumEffDepth, float MaximumSlope, float v,
                          float lambda, LSDRaster& Elevation)
{
  // Don't let V be less than or equal to zero
  if (v <= 0)
  {
    v = 0.001;
  }

  // some variables to speed up compuation
  float exp_term;
  float thisEffDepth = 0;
  float elev_mulitplier = (MaximumSlope/MaximumEffDepth)*pow((1+v),1+(1/v));
  float LocalElevation;

  for (int row = 0; row <NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if (RasterData[row][col] != NoDataValue)
      {
        LocalElevation = Elevation.get_data_element(row,col);

        if (LocalElevation != NoDataValue)
        {
          // get the effective depth using the richards sigmoidal gorth function
          exp_term = 1+v*exp(elev_mulitplier*(lambda-LocalElevation));
          thisEffDepth = MaximumEffDepth*pow(exp_term,-(1/v));

          // if the depth is less than zero, then set to zero
          if(thisEffDepth <0)
          {
            thisEffDepth = 0;
          }

          // update the data
          RasterData[row][col] =  thisEffDepth;

        }
        else        // if there ins't any elevation data, set the snow data to NoData
        {
          RasterData[row][col] = NoDataValue;
        }
      }
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is an incredibly rudimentary function used to modify landslide raster
// It takes a few rasters from the
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSoilHydroRaster::NaiveLandslide(LSDRaster& FilledElevation, int initiationPixels,
                                      int MinPixels, float landslide_thickness)
{
  // Get a flow info object
  // set no flux boundary conditions
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";


  // some values from the rasters
  float local_elev;
  float local_mask;

  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,FilledElevation);

  // get the contributing pixels
  LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  vector<int> sources = FlowInfo.get_sources_index_threshold(ContributingPixels, initiationPixels);

  // get a value vector for the landslides
  vector<float> landslide_thicknesses;
  for (int i = 0; i< int(sources.size()); i++)
  {
    landslide_thicknesses.push_back(landslide_thickness);
  }


  // get the mask
  LSDRaster Mask = FlowInfo.get_upslope_node_mask(sources,landslide_thicknesses);

  // now set all points that have elevation data but not landslide data to
  // the value of the landslide thickness, removing data that is below the minium
  // pixel area
  for (int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      local_elev =  FilledElevation.get_data_element(row,col);
      local_mask =  Mask.get_data_element(row,col);

      RasterData[row][col] = local_mask;

      // Turn nodata points into 0s
      if( local_mask == NoDataValue)
      {
        RasterData[row][col] = 0.0;
      }

      // remove data where there is no topographic information
      if( local_elev == NoDataValue)
      {
        RasterData[row][col] = NoDataValue;
      }
    }
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate h, the soil depth normal to the slope, used in the factor of safety equation.
// Call with the soil thickness raster.
// SWDG 13/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::Calculate_h(LSDRaster& Slope){

  Array2D<float> h(NRows, NCols, NoDataValue);

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      if (RasterData[i][j] != NoDataValue){
        h[i][j] = RasterData[i][j] * cos(Slope.get_data_element(i,j));

      }
    }
  }

  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,h,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate w, a hyrdological index, used in the factor of safety equation.
// Call with the ratio of recharge to transmissivity.
// SWDG 13/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::Calculate_w(LSDRaster& Slope, LSDRaster& DrainageArea){

  Array2D<float> w(NRows, NCols, NoDataValue);

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      if (RasterData[i][j] != NoDataValue){

        float value = RasterData[i][j] * (DrainageArea.get_data_element(i,j)/sin(Slope.get_data_element(i,j)));
        if (value < 1.0){
          w[i][j] = value;
        }
        else{
          w[i][j] = 1.0;
        }

      }
    }
  }

  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,w,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate r, the soil to water density ratio, used in the factor of safety equation.
// Call with the soil density raster.
// SWDG 13/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::Calculate_r(float& rhoW){

  Array2D<float> r(NRows, NCols, NoDataValue);

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      if (RasterData[i][j] != NoDataValue){
        r[i][j] = rhoW / RasterData[i][j];

      }
    }
  }

  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,r,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// @brief Calculate C, a cohesion index, used in the factor of safety equation.
// Call with the root cohesion raster.
// SWDG 13/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::Calculate_C(LSDSoilHydroRaster& Cs, LSDSoilHydroRaster& h, LSDSoilHydroRaster& rhoS, float& g){

  Array2D<float> C(NRows, NCols, NoDataValue);

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      if (RasterData[i][j] != NoDataValue){
        C[i][j] = (RasterData[i][j] + Cs.get_data_element(i,j)) / (h.get_data_element(i,j)*rhoS.get_data_element(i,j)*g);

      }
    }
  }

  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,C,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate the factor of safety using the sinmap definition.
// Call with the dimensionless cohesion (C) raster.
// SWDG 13/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::Calculate_sinmap_Fs(LSDRaster& Slope, LSDSoilHydroRaster& w, LSDSoilHydroRaster& r, LSDSoilHydroRaster& phi){

  Array2D<float> Fs(NRows, NCols, NoDataValue);

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      if (RasterData[i][j] != NoDataValue){
        Fs[i][j] = ( RasterData[i][j] + cos(Slope.get_data_element(i,j)) * (1.0-w.get_data_element(i,j)*r.get_data_element(i,j)) * tan(phi.get_data_element(i,j)) ) / sin(Slope.get_data_element(i,j));
      }
    }
  }

  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,Fs,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate the the sinmap Stability Index (SI).
// This is a wrapper around a lightly modified port of the original sinmap 2.0 implementation.
// call with any SoilHydroRaster, it's values are used for identification of NoDataValues.
//
// SWDG 15/6/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::Calculate_sinmap_SI(LSDRaster Slope, LSDRaster DrainageArea, LSDSoilHydroRaster lo_C, LSDSoilHydroRaster hi_C, LSDSoilHydroRaster lo_phi, LSDSoilHydroRaster hi_phi, LSDSoilHydroRaster lo_RoverT, LSDSoilHydroRaster hi_RoverT, LSDSoilHydroRaster lo_r, LSDSoilHydroRaster hi_r, LSDSoilHydroRaster lo_FS, LSDSoilHydroRaster hi_FS){

  Array2D<float> SI(NRows, NCols, NoDataValue);

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      if (RasterData[i][j] != NoDataValue){
        SI[i][j] = StabilityIndex(Slope.get_data_element(i, j), DrainageArea.get_data_element(i, j), lo_C.get_data_element(i, j), hi_C.get_data_element(i, j), lo_phi.get_data_element(i, j), hi_phi.get_data_element(i, j), lo_RoverT.get_data_element(i, j), hi_RoverT.get_data_element(i, j), lo_r.get_data_element(i, j), hi_r.get_data_element(i, j), lo_FS.get_data_element(i, j), hi_FS.get_data_element(i, j));
      }
    }
  }

  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,SI,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate the average stability index value for a given set of upslope nodes.
// SWDG 25/7/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float LSDSoilHydroRaster::AverageSI(vector<int>& Nodes, LSDFlowInfo& FlowInfo){

  int i;
  int j;
  float total = 0;
  int count = 0;

  for (int q = 0; q < int(Nodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(Nodes[q], i, j);
    if (RasterData[i][j] != NoDataValue){
      total += RasterData[i][j];
      count++;
    }
  }

  cout << total/count << endl;

  return total/count;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate the average stability index value for a vector containing many sets of upslope nodes.
// SWDG 25/7/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDSoilHydroRaster::AverageSIs(vector< vector<int> >& Vector_of_Nodes, LSDFlowInfo& FlowInfo){
  vector<float> output;
  for (int q = 0; q < int(Vector_of_Nodes.size()); ++q){
    output.push_back(AverageSI(Vector_of_Nodes[q], FlowInfo));
  }
  return output;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// @brief Write average SI values to spatial units defined by a vector of vectors of nodes.
// Pass in a vector of average values generated by AverageSIs().
//SWDG 25/7/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDSoilHydroRaster LSDSoilHydroRaster::WriteAvgSIs(vector<float> SIs, vector< vector<int> >& Vector_of_Nodes, LSDFlowInfo& FlowInfo){

  int i;
  int j;

  Array2D<float> output(NRows, NCols, NoDataValue);
  for (int q = 0; q < int(Vector_of_Nodes.size()); ++q){
    for (int w = 0; w < int(Vector_of_Nodes[q].size()); ++w){

      FlowInfo.retrieve_current_row_and_col(Vector_of_Nodes[q][w], i, j);
      output[i][j] = SIs[q];

    }
  }

  LSDSoilHydroRaster AvgSIs(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,output,GeoReferencingStrings);
  return AvgSIs;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Implementing the terrain shape index as used by Bolstad et al (1998)
// Predicting Southern Appalachian overstory vegetation with digital terrain data
// which is based on work by McNab (1989) Terrain Shape Index: Quantifying Effect of
// Minor Landforms on Tree Height
//
// These papers use a lot of magic numbers, which turn out to be fudge factors for
// working with DEMs with resolutions reported in YARDS. This code will cope with any
// stupid spatial unit you choose.
//
//SWDG - 24/1/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDSoilHydroRaster LSDSoilHydroRaster::TerrainShapeIndex(){

  float elevationSum;
  int dX[] = {1, 1, 1, 0, -1, -1, -1, 0};
  int dY[] = {-1, 0, 1, 1, 1, 0, -1, -1};

  Array2D<float> Output(NRows, NCols, NoDataValue);

  //compute the mean intercell distance - Assumes square cells!
  float diagonal = sqrt((2 * (DataResolution * DataResolution)));
  float intercellDist = ((12 * DataResolution) + (8 * diagonal)) / 20.0;

  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){

      //check 9 cells for any ndv
      if ((RasterData[i][j] != NoDataValue) && (RasterData[i-1][j-1] != NoDataValue) && (RasterData[i][j-1] != NoDataValue) && (RasterData[i+1][j-1] != NoDataValue) && (RasterData[i+1][j] != NoDataValue) && (RasterData[i+1][j+1] != NoDataValue) && (RasterData[i][j+1] != NoDataValue) && (RasterData[i-1][j+1] != NoDataValue) && (RasterData[i-1][j] != NoDataValue)){

        elevationSum = 0;

        for (int c=0; c < 8; ++c){
          elevationSum += (RasterData[i][j] - RasterData[i + dY[c]][j + dX[c]]);
        }
        Output[i][j] = (elevationSum / 8.0) / intercellDist;

      }
    }
  }

  Array2D<float> Normed = normalize_terrain_index(Output, 2, 98, 10000, NoDataValue);
  LSDSoilHydroRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,Normed,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Method to create a vegetation map based on Bolstad 1998. Calibrated for the
// Southern Appalachians.
//
// Call with an LSDSoilHydroRaster generated by TerrainShapeIndex()
// SWDG - 26/1/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDSoilHydroRaster::BolstadVegetationMap(LSDRaster& DEM){

  Array2D<int> Output(NRows, NCols, NoDataValue);
  Array2D<float> DEMData = DEM.get_RasterData();

  // Create arrays to store the values of each species type
  float CoveX[] = {0, 3500, 3500, 0};
  float CoveY[] = {1065, 1182, 600, 600};

  float DecidX[] = {3500, 6562.5, 7812.5, 7812.5, 7375, 6562.5, 5250, 3500};
  float DecidY[] = {1182, 1287, 1322, 839, 755, 720, 600, 600};

  float XericX[] = {5250, 6562.5, 7375, 7812.5, 7812.5, 10001, 10001};
  float XericY[] = {600, 720, 755, 839, 1322, 1322, 600};

  float NorthX[] = {0, 3500, 6562.5, 7812.5, 10001, 10001, 0};
  float NorthY[] = {1065, 1182, 1287, 1322, 1322, 1600, 1600};

  //for each point in the Raster, check if the point falls within any of the 4 species bands
  //code the output raster 1 = Cove Hardwood, 2 = Mixed Deciduous, 3 = Xeric Oak-Pine 4 = Northern Hardwood
  for (int i = 1; i < NRows - 1; ++i){
    for (int j = 1; j < NCols - 1; ++j){
      if (RasterData[i][j] != NoDataValue){
        int cove = 0;
        int decid = 0;
        int xeric = 0;
        int north = 0;

        cove = PointInPolygon(sizeof(CoveX)/sizeof(*CoveX), CoveX, CoveY, RasterData[i][j], DEMData[i][j]);

        if (cove % 2 == 1){
          Output[i][j] = 1;
        }

        else{
          decid = PointInPolygon(sizeof(DecidX)/sizeof(*DecidX), DecidX, DecidY, RasterData[i][j], DEMData[i][j]);

          if (decid % 2 == 1){
            Output[i][j] = 2;
          }
          else{
            xeric = PointInPolygon(sizeof(XericX)/sizeof(*XericX), XericX, XericY, RasterData[i][j], DEMData[i][j]);

            if (xeric % 2 == 1){
              Output[i][j] = 3;
            }
            else{
              north = PointInPolygon(sizeof(NorthX)/sizeof(*NorthX), NorthX, NorthY, RasterData[i][j], DEMData[i][j]);

              if (north % 2 == 1){
                Output[i][j] = 4;
              }
              else{
                Output[i][j] = NoDataValue;
              }
            }
          }
        }
      }
    }
  }

  LSDIndexRaster output(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,Output,GeoReferencingStrings);
  return output;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#endif
