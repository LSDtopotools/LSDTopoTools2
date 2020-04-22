//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDBasin
// Land Surface Dynamics Basin
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing basins
//
// Developed by:
//  Stuart W.D. Grieve
//  Simon M. Mudd
//  Fiona Clubb
//
// Copyright (C) 2017 Simon M. Mudd 2017
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

#ifndef LSDBasin_CPP
#define LSDBasin_CPP

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDBasin.hpp"
#include "LSDParticle.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDSpatialCSVReader.hpp"
using namespace std;
using namespace TNT;

void LSDBasin::create()
{


  // now we set up empty variables to store properties of the basin
  // these are populated as they are required using the set methods in LSDBasin

  SlopeMean = NoDataValue;
  ElevationMean = NoDataValue;
  AspectMean = NoDataValue;
  ReliefMean = NoDataValue;
  PlanCurvMean = NoDataValue;
  ProfileCurvMean = NoDataValue;
  TotalCurvMean = NoDataValue;
  PlanCurvMax = NoDataValue;
  ProfileCurvMax = NoDataValue;
  TotalCurvMax = NoDataValue;
  HillslopeLength_HFR = NoDataValue;
  HillslopeLength_Binned = NoDataValue;
  HillslopeLength_Spline = NoDataValue;
  HillslopeLength_Density = NoDataValue;
  FlowLength = NoDataValue;
  DrainageDensity = NoDataValue;
  Perimeter_i = vector<int>(1,NoDataValue);
  Perimeter_j =  vector<int>(1,NoDataValue);
  Perimeter_nodes =  vector<int>(1,NoDataValue);
  CosmoErosionRate = NoDataValue;
  OtherErosionRate = NoDataValue;
  CHTMean = NoDataValue;
  EStar = NoDataValue;
  RStar = NoDataValue;
  HilltopPx = NoDataValue;
  BedrockFraction = NoDataValue;
  Biomass = NoDataValue;
  AlternativeIndex=int(NoDataValue);
  DD_preprocessed = false;

  //finished creating empty variables

}

void LSDBasin::create(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet)
{

  //NO BOUNDS CHECKING ON JunctionNumber

  //setting all of the instance variables for the given junction

  NRows = ChanNet.get_NRows();
  NCols = ChanNet.get_NCols();
  XMinimum = ChanNet.get_XMinimum();
  YMinimum = ChanNet.get_YMinimum();
  DataResolution = ChanNet.get_DataResolution();
  NoDataValue = ChanNet.get_NoDataValue();
  GeoReferencingStrings = ChanNet.get_GeoReferencingStrings();

  Junction = JunctionNumber;

  vector <int> JunctionVector = ChanNet.get_JunctionVector();
  vector <int> ReceiverVector = ChanNet.get_ReceiverVector();

//	int ReceiverJunction = 4688;
//	int ReceiverNode = ChanNet.get_Node_of_Junction(ReceiverJunction);

  LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction,JunctionVector[Junction],
                                                     ReceiverVector[Junction], JunctionVector[ReceiverVector[Junction]], FlowInfo);

//	LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction, JunctionVector[Junction],
//                                                     ReceiverJunction, ReceiverNode, FlowInfo);

  int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  int basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  BasinNodes = FlowInfo.get_upslope_nodes(basin_outlet);

  NumberOfCells = int(BasinNodes.size());
  Area = NumberOfCells * (DataResolution*DataResolution);

  Beheaded = ChanNet.node_tester(FlowInfo, Junction);

  FlowInfo.retrieve_current_row_and_col(ChanNet.get_Node_of_Junction(Junction), Outlet_i, Outlet_j);

  vector<int> StreamOrderVector = ChanNet.get_StreamOrderVector();

  BasinOrder = StreamOrderVector[Junction];


  int i_max = 0;
  int i_min = 9999999; //a very large number
  int j_max = 0;
  int j_min = 9999999; //a very large number

  int i = 0;
  int j = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    if (i > i_max){i_max = i;}
    else if (i < i_min){i_min = i;}
    if (j > j_max){j_max = j;}
    else if (j < j_min){j_min = j;}

  }

  Centroid_i = i_min + ((i_max - i_min)/2);
  Centroid_j = j_min + ((j_max - j_min)/2);   //how do these handle 0.5s ??


  //finished setting all the instance variables


  // now we set up empty variables to store properties of the basin
  // these are populated as they are required using the set methods in LSDBasin

  SlopeMean = NoDataValue;
  ElevationMean = NoDataValue;
  AspectMean = NoDataValue;
  ReliefMean = NoDataValue;
  PlanCurvMean = NoDataValue;
  ProfileCurvMean = NoDataValue;
  TotalCurvMean = NoDataValue;
  PlanCurvMax = NoDataValue;
  ProfileCurvMax = NoDataValue;
  TotalCurvMax = NoDataValue;
  HillslopeLength_HFR = NoDataValue;
  HillslopeLength_Binned = NoDataValue;
  HillslopeLength_Spline = NoDataValue;
  HillslopeLength_Density = NoDataValue;
  FlowLength = NoDataValue;
  DrainageDensity = NoDataValue;
  Perimeter_i = vector<int>(1,NoDataValue);
  Perimeter_j =  vector<int>(1,NoDataValue);
  Perimeter_nodes =  vector<int>(1,NoDataValue);
  CosmoErosionRate = NoDataValue;
  OtherErosionRate = NoDataValue;
  CHTMean = NoDataValue;
  EStar = NoDataValue;
  RStar = NoDataValue;
  BedrockFraction = NoDataValue;
  Biomass = NoDataValue;
  AlternativeIndex=int(NoDataValue);

  DD_preprocessed = false;
  //finished creating empty variables

}

void LSDBasin::create(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet, int extended_baslevel_nodes_below_junction)
{  
//NO BOUNDS CHECKING ON JunctionNumber

  //setting all of the instance variables for the given junction

  NRows = ChanNet.get_NRows();
  NCols = ChanNet.get_NCols();
  XMinimum = ChanNet.get_XMinimum();
  YMinimum = ChanNet.get_YMinimum();
  DataResolution = ChanNet.get_DataResolution();
  NoDataValue = ChanNet.get_NoDataValue();
  GeoReferencingStrings = ChanNet.get_GeoReferencingStrings();

  Junction = JunctionNumber;

  vector <int> JunctionVector = ChanNet.get_JunctionVector();
  vector <int> ReceiverVector = ChanNet.get_ReceiverVector();

//  int ReceiverJunction = 4688;
//  int ReceiverNode = ChanNet.get_Node_of_Junction(ReceiverJunction);

  LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction,JunctionVector[Junction],
                                                     ReceiverVector[Junction], JunctionVector[ReceiverVector[Junction]], FlowInfo);

//  LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction, JunctionVector[Junction],
//                                                     ReceiverJunction, ReceiverNode, FlowInfo);

  int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  int basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  BasinNodes = FlowInfo.get_upslope_nodes(extended_baslevel_nodes_below_junction);

  NumberOfCells = int(BasinNodes.size());
  Area = NumberOfCells * (DataResolution*DataResolution);

  Beheaded = ChanNet.node_tester(FlowInfo, Junction);

  FlowInfo.retrieve_current_row_and_col(extended_baslevel_nodes_below_junction, Outlet_i, Outlet_j);

  vector<int> StreamOrderVector = ChanNet.get_StreamOrderVector();

  BasinOrder = StreamOrderVector[Junction];


  int i_max = 0;
  int i_min = 9999999; //a very large number
  int j_max = 0;
  int j_min = 9999999; //a very large number

  int i = 0;
  int j = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    if (i > i_max){i_max = i;}
    else if (i < i_min){i_min = i;}
    if (j > j_max){j_max = j;}
    else if (j < j_min){j_min = j;}

  }

  Centroid_i = i_min + ((i_max - i_min)/2);
  Centroid_j = j_min + ((j_max - j_min)/2);   //how do these handle 0.5s ??


  //finished setting all the instance variables


  // now we set up empty variables to store properties of the basin
  // these are populated as they are required using the set methods in LSDBasin

  SlopeMean = NoDataValue;
  ElevationMean = NoDataValue;
  AspectMean = NoDataValue;
  ReliefMean = NoDataValue;
  PlanCurvMean = NoDataValue;
  ProfileCurvMean = NoDataValue;
  TotalCurvMean = NoDataValue;
  PlanCurvMax = NoDataValue;
  ProfileCurvMax = NoDataValue;
  TotalCurvMax = NoDataValue;
  HillslopeLength_HFR = NoDataValue;
  HillslopeLength_Binned = NoDataValue;
  HillslopeLength_Spline = NoDataValue;
  HillslopeLength_Density = NoDataValue;
  FlowLength = NoDataValue;
  DrainageDensity = NoDataValue;
  Perimeter_i = vector<int>(1,NoDataValue);
  Perimeter_j =  vector<int>(1,NoDataValue);
  Perimeter_nodes =  vector<int>(1,NoDataValue);
  CosmoErosionRate = NoDataValue;
  OtherErosionRate = NoDataValue;
  CHTMean = NoDataValue;
  EStar = NoDataValue;
  RStar = NoDataValue;
  BedrockFraction = NoDataValue;
  Biomass = NoDataValue;
  AlternativeIndex=int(NoDataValue);

  DD_preprocessed = false;
  //finished creating empty variables this-create(JunctionNumber,FlowInfo, ChanNet);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate mean basin value.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMean(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float TotalData = 0;
  int CountNDV = 0;
  float BasinAverage;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    //exclude NDV from average
    if (Data.get_data_element(i,j) != NoDataValue){
      TotalData += Data.get_data_element(i,j);
    }
    else {
      ++CountNDV;
    }
  }

  BasinAverage = TotalData/(NumberOfCells-CountNDV);

  return BasinAverage;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate mean basin value. overloaded for lsdindexraster - FJC
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMean(LSDFlowInfo& FlowInfo, LSDIndexRaster Data){

  int i;
  int j;
  int TotalData = 0;
  int CountNDV = 0;
  float BasinAverage;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    //exclude NDV from average
    if (Data.get_data_element(i,j) != NoDataValue){
      TotalData += Data.get_data_element(i,j);
    }
    else {
      ++CountNDV;
    }
  }

  BasinAverage = TotalData/(NumberOfCells-CountNDV);

  return BasinAverage;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate max basin value.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMax(LSDFlowInfo& FlowInfo, LSDRaster Data){

  //could use max_element here? how would that cope with NDVs??

  int i;
  int j;
  float MaxData = -10000000;   //a very small number
  float CurrentData;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);

    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData > MaxData){
      MaxData = CurrentData;
    }
  }

  return MaxData;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate min basin value.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMin(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float MinData = 100000000; // a large number
  float CurrentData;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);

    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData < MinData){
      MinData = CurrentData;
    }
  }

  return MinData;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate median basin value.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMedian(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  vector<float> UnsortedData;
  vector<float> SortedData;
  vector<size_t> index_map;
  float Median;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      UnsortedData.push_back(Data.get_data_element(i,j));
    }
  }

  //get size of dataset
  size_t n = UnsortedData.size() / 2;

  //sort all non NDV values
  matlab_float_sort(UnsortedData, SortedData, index_map);

  if (n % 2 != 0){
    Median = (SortedData[n] + SortedData[n+1])/2;
  }
  else{
    Median = SortedData[n];
  }

  return Median;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate percentile basin value.
// MDH 5/2/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinPercentile(LSDFlowInfo& FlowInfo, LSDRaster Data, int Percentile)
{

	int i;
	int j;
	vector<float> UnsortedData;
	vector<float> SortedData;
	vector<size_t> index_map;
	float P, PercentileValue, Residual;

	for (int q = 0; q < int(BasinNodes.size()); ++q)
	{
		FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
		//exclude NDV
		if (Data.get_data_element(i,j) != NoDataValue)
		{
			UnsortedData.push_back(Data.get_data_element(i,j));
		}
	}

	//get size of dataset
	size_t n = UnsortedData.size();

	//sort all non NDV values
	matlab_float_sort(UnsortedData, SortedData, index_map);

	//find index for percentile
	P = Percentile*(n/100.);
	int Pint = round(P);

	//Interpolate to get percentile value
	Residual = P-Pint;
	if (Residual > 0)
	{
		PercentileValue = SortedData[Pint] + Residual*(SortedData[Pint+1]-SortedData[Pint]);
	}
	else
	{
		PercentileValue = SortedData[Pint] + Residual*(SortedData[Pint]-SortedData[Pint-1]);
	}

	return PercentileValue;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate Standard devaition of the basin values.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinStdDev(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  vector<float> DataValues;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      DataValues.push_back(Data.get_data_element(i,j));
    }
  }

  float mean = get_mean(DataValues);
  float StdDev = get_standard_deviation(DataValues, mean);

  return StdDev;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate Standard Error of the basin values.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinStdError(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  vector<float> DataValues;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      DataValues.push_back(Data.get_data_element(i,j));
    }
  }

  float mean = get_mean(DataValues);
  float StdDev = get_standard_deviation(DataValues, mean);
  float StdError = get_standard_error(DataValues, StdDev);

  return StdError;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate basin range.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinRange(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float MinData = 100000000; // a large number
  float MaxData = -100000000; // a small number
  float CurrentData;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);

    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData < MinData){
      MinData = CurrentData;
    }
    if (CurrentData != NoDataValue && CurrentData > MaxData){
      MaxData = CurrentData;
    }
  }
  float Range = MaxData - MinData;
  return Range;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate basin range.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDBasin::CalculateNumDataPoints(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  int count = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      ++count;
    }
  }

  return count;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate E* and R* values for the basin, using hilltop flow routed hillslope
// lengths.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_EStar_RStar(float CriticalSlope){

    EStar = (2 * (abs(CHTMean)) * HillslopeLength_HFR) / CriticalSlope;
    RStar = ReliefMean / (HillslopeLength_HFR * CriticalSlope);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate flow length for the basin using the D8 flow directions.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_FlowLength(LSDIndexRaster& StreamNetwork, LSDFlowInfo& FlowInfo){

  int j;
  int i;
  float LengthSum = 0;
  float two_times_root2 = 2.828427;
  Array2D<int> FlowDir = FlowInfo.get_FlowDirection();


  //Loop over every pixel and record it's stream length and basin ID in two vectors
  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

     if (StreamNetwork.get_data_element(i,j) != NoDataValue){
       if ((FlowDir[i][j] % 2) != 0 && (FlowDir[i][j] != -1 )){ //is odd but not -1
         LengthSum += (DataResolution * two_times_root2); //diagonal
       }
       else if (FlowDir[i][j] % 2 == 0){  //is even
         LengthSum +=  DataResolution; //cardinal
       }
     }
  }

  FlowLength = LengthSum;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate hillslope lengths from boomerang plots.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_HillslopeLengths_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold){

  int j;
  int i;
  Array2D<float> slope(NRows, NCols, NoDataValue);
  Array2D<float> area(NRows, NCols, NoDataValue);

  //create subset arrays for just the basin data - this should be rolled into its own method.
  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

      slope[i][j] = Slope.get_data_element(i,j);
      area[i][j] = DinfArea.get_data_element(i,j);

  }

  //do some log binning
  vector<float> Mean_x_out;
  vector<float> Mean_y_out;
  vector<float> Midpoints_out;
  vector<float> STDDev_x_out;
  vector<float> STDDev_y_out;
  vector<float> STDErr_x_out;
  vector<float> STDErr_y_out;
  vector<int> number_observations;

  log_bin_data(area, slope, log_bin_width, Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, NoDataValue);

  //remove empty bins
  RemoveSmallBins(Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, bin_threshold);

  //index value of max slope
  int slope_max_index = distance(Mean_y_out.begin(), max_element(Mean_y_out.begin(), Mean_y_out.end()));

  //hillslope length from the maximum binned values
  HillslopeLength_Binned = Mean_x_out[slope_max_index]/DataResolution;

  // Fit splines through the binned data to get the LH
  vector<float> Spline_X;
  vector<float> Spline_Y;
  PlotCubicSplines(Mean_x_out, Mean_y_out, SplineResolution, Spline_X, Spline_Y);

  //index value of max spline slope
  int slope_max_index_spline = distance(Spline_Y.begin(), max_element(Spline_Y.begin(), Spline_Y.end()));

  //hillslope length from spline curve
  HillslopeLength_Spline = Spline_X[slope_max_index_spline]/DataResolution;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Generate data to create boomerang plots.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::Plot_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold, string Path){

  int j;
  int i;
  Array2D<float> slope(NRows, NCols, NoDataValue);
  Array2D<float> area(NRows, NCols, NoDataValue);

  //create subset arrays for just the basin data - this should be rolled into its own method.
  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

      slope[i][j] = Slope.get_data_element(i,j);
      area[i][j] = DinfArea.get_data_element(i,j);

  }

  //do some log binning
  vector<float> Mean_x_out;
  vector<float> Mean_y_out;
  vector<float> Midpoints_out;
  vector<float> STDDev_x_out;
  vector<float> STDDev_y_out;
  vector<float> STDErr_x_out;
  vector<float> STDErr_y_out;
  vector<int> number_observations;

  log_bin_data(area, slope, log_bin_width, Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, NoDataValue);

  //remove empty bins
  RemoveSmallBins(Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, bin_threshold);

  // Fit splines through the binned data to get the LH
  vector<float> Spline_X;
  vector<float> Spline_Y;
  PlotCubicSplines(Mean_x_out, Mean_y_out, SplineResolution, Spline_X, Spline_Y);

  //set up a filestream object to write the binned data
  ofstream file;

  stringstream ss_bin;
  ss_bin << Path << Junction << "_boom_binned.txt";
  file.open(ss_bin.str().c_str());   //needs a null terminated character array, not a string. See pg 181 of accelerated c++

  for(int q = 0; q < int(Mean_x_out.size()); q++){
    file << Mean_x_out[q] << " " << Mean_y_out[q] << " " << STDDev_x_out[q] << " " << STDDev_y_out[q] << " " << STDErr_x_out[q] << " " << STDErr_y_out[q] << endl;
  }
  file.close();

  //set up a filestream object to write the spline data
  ofstream SplineFile;

  stringstream ss_spline;
  ss_spline << Path << Junction << "_boom_spline.txt";
  SplineFile.open(ss_spline.str().c_str());   //needs a null terminated character array, not a string. See pg 181 of accelerated c++

  for(int q = 0; q < int(Spline_X.size()); q++){ //fixed bug here where I looped over the wrong vector - SWDG 7/11/13
    SplineFile << Spline_X[q] << " " << Spline_Y[q] << endl;

  }
  SplineFile.close();

  //set up a filestream object to write the data cloud
  ofstream cloud;

  stringstream ss_cloud;
  ss_cloud << Path << Junction << "_boom_cloud.txt";
  cloud.open(ss_cloud.str().c_str());     //needs a null terminated character array, not a string. See pg 181 of accelerated c++

  for (int i = 1; i < NRows-1; ++i){
    for (int j = 1; j < NCols-1; ++j){
      if(area[i][j] != NoDataValue && slope[i][j] != NoDataValue){
        cloud << area[i][j] << " " << slope[i][j] << endl;
      }
    }
  }
  cloud.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the mean basin aspect. Does not use the normal basin mean method as angles
// need to be handled differently.
// Bug fixed in the average calculation when values wrapped around 0
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_AspectMean(LSDFlowInfo& FlowInfo, LSDRaster Aspect){

  int i;
  int j;
  float avg_r;
  float angle_r;
  float x_component = 0.0;
  float y_component = 0.0;
  int ndv_cell_count = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    if (Aspect.get_data_element(i,j) != NoDataValue){

      angle_r = rad(Aspect.get_data_element(i,j));
      x_component += cos(angle_r);
      y_component += sin(angle_r);

    }
    else{
      ++ndv_cell_count;
    }

  }

  avg_r = atan2(y_component, x_component);
  AspectMean = deg(avg_r);

  if (AspectMean < 0){
    AspectMean = 360 + AspectMean;
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the perimeter pixels using a simple edge detection algorithm. This is quite
// messy and will be improved soon.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_Perimeter(LSDFlowInfo& FlowInfo)
{

  int i;
  int j;
  vector<int> I;
  vector<int> J;
  vector<int> B;
  int NDVCount = 0;
  Array2D<float> BasinData(NRows, NCols, NoDataValue);

  //create subset arrays for just the basin data - this should be rolled into its own method.
  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
      BasinData[i][j] = BasinNodes[q];

  }

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    NDVCount = 0;
    // cout << i << "||" << j << endl;

      if (i != 0 && j != 0 && i<NRows-1 && j < NCols-1)
      {
        //count border cells that are NDV

        if (BasinData[i-1][j-1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i][j-1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i+1][j-1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i-1][j] == NoDataValue){ ++NDVCount; }
        if (BasinData[i+1][j] == NoDataValue){ ++NDVCount; }
        if (BasinData[i-1][j+1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i][j+1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i+1][j+1] == NoDataValue){ ++NDVCount; }

        if (NDVCount >= 1 && NDVCount < 8)
        {  //increase the first value to get a simpler polygon  (changed to 1 by FJC 23/03/15 to get only internal hilltops.
          //edge pixel                       // Otherwise not all external ridges were being excluded from the analysis).
          I.push_back(i);
          J.push_back(j);
          B.push_back(BasinNodes[q]);
        }
      }
      else
      {
        ++i;
      }

  }

  //now have 2 vectors of i and j indexes of every point
  Perimeter_i = I;
  Perimeter_j = J;
  Perimeter_nodes = B;


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the perimeter to csv. It can then be ingested to find
// concave hull of the basin (or the basin outline)
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::print_perimeter_to_csv(LSDFlowInfo& FlowInfo, string perimeter_fname)
{
  // make sure we have found the perimeter
  if (int(Perimeter_nodes.size()) == 0)
  {
    set_Perimeter(FlowInfo);
  }

  // open the file
  ofstream perim_out;
  perim_out.open(perimeter_fname.c_str());
  perim_out << "node,x,y,latitude,longitude,row,col" << endl;
  perim_out.precision(9);

  float curr_x,curr_y;
  double curr_lat,curr_long;
  int this_row =0, this_col =0 ;

  LSDCoordinateConverterLLandUTM converter;
  int n_nodes = int(Perimeter_nodes.size());
  for(int i = 0; i< n_nodes; i++)
  {
    FlowInfo.get_x_and_y_from_current_node(Perimeter_nodes[i], curr_x, curr_y);
    FlowInfo.retrieve_current_row_and_col(Perimeter_nodes[i], this_row, this_col);
    FlowInfo.get_lat_and_long_from_current_node(Perimeter_nodes[i], curr_lat, curr_long,converter);
    perim_out << Perimeter_nodes[i] << "," << curr_x << "," << curr_y <<"," << curr_lat << "," << curr_long << "," << this_row << "," << this_col << endl;
  }
  perim_out.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the perimeter to csv. It can then be ingested to find
// concave hull of the basin (or the basin outline)
// Also prints the elevations so can investigate the perimeter hypsometry
// FJC 10/01/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::print_perimeter_hypsometry_to_csv(LSDFlowInfo& FlowInfo, string perimeter_fname, LSDRaster& ElevationRaster)
{

  set_Perimeter(FlowInfo);
  clean_perimeter(FlowInfo);

  // TESTING
  string perim_test = "/home/clubb/Data_for_papers/drainage_capture/Santa_Cruz/FUCK_THIS.csv";
  print_perimeter_to_csv(FlowInfo, perim_test);

  //  open the file
  ofstream perim_out;
  perim_out.open(perimeter_fname.c_str());
  perim_out << "node_key,node,elevation,x,y,latitude,longitude,dist_from_outlet_node" << endl;
  perim_out.precision(9);

  double curr_lat,curr_long;

  LSDCoordinateConverterLLandUTM converter;
  int n_nodes = int(Perimeter_nodes.size());
  cout << "N perimeter nodes: " << n_nodes << endl;

  int outlet_node = get_Outlet_node();

  // sort perimeter nodes
  vector<int> Reordered_nodes = order_perimeter_nodes(FlowInfo);

  // sanity checks

  cout << "n unordered nodes: " << n_nodes << " n sorted nodes: " << Reordered_nodes.size() << endl;

  for(int i = 0; i< int(Reordered_nodes.size()); i++)
  {
    // get this elevation
    int this_row, this_col;
    FlowInfo.retrieve_current_row_and_col(Reordered_nodes[i], this_row, this_col);
    float this_elev = ElevationRaster.get_data_element(this_row, this_col);

    // get the x and y from this node
    float curr_x, curr_y;
    FlowInfo.get_x_and_y_from_current_node(Reordered_nodes[i], curr_x, curr_y);

    // get the coordinates
    FlowInfo.get_lat_and_long_from_current_node(Reordered_nodes[i], curr_lat, curr_long,converter);

    // get the euclidian distance from the outlet junction
    float dist = FlowInfo.get_Euclidian_distance(outlet_node, Reordered_nodes[i]);

    // write to csv
    perim_out << i << "," << Reordered_nodes[i] << "," << this_elev << "," << curr_x << "," << curr_y <<"," << curr_lat << "," << curr_long << "," << dist << endl;
  }

  perim_out.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Order perimeter nodes from the outlet
// Must CLEAN the perimeter first using clean_perimeter function
// FJC 16/01/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDBasin::order_perimeter_nodes(LSDFlowInfo& FlowInfo)
{
  // first clean the perimeter nodes
  //clean_perimeter(FlowInfo);
  cout << "Ordering the perimeter nodes..." << endl;
  vector<int> sorted_nodes;
  Array2D<int> PerimeterNodes(NRows,NCols,0);
  Array2D<int> VisitedBefore(NRows,NCols,0);

  // first get the outlet node
  int outlet_node = get_Outlet_node();
  int outlet_row, outlet_col;
  FlowInfo.retrieve_current_row_and_col(outlet_node, outlet_row, outlet_col);
  Outlet_i = outlet_row;
  Outlet_j = outlet_col;
  cout << "The outlet node is: " << outlet_node << endl;

  // get an array of perimeter nodes
  for (int i = 0; i < int(Perimeter_nodes.size()); i++)
  {
    int this_row, this_col;
    FlowInfo.retrieve_current_row_and_col(Perimeter_nodes[i], this_row, this_col);
    PerimeterNodes[this_row][this_col] = 1;
  }
  // The outlet always needs to be in the perimeter
  PerimeterNodes[Outlet_i][Outlet_j] = 1;
  cout << "Got the array of perimeter nodes" << endl;

  //find the nearest perimeter to the outlet node
  VisitedBefore[Outlet_i][Outlet_j] = 2;
  // push back the outlet node, node 0
  sorted_nodes.push_back(outlet_node);
  int next_i, next_j;
  bool first_node = true; // bool to check if you're at the first node. This makes a difference because if you're at the first node we don't want to go back to the outlet.

  // N, S, E, and W will always be the shortest distances, so do these first
  if (PerimeterNodes[Outlet_i][Outlet_j-1] == 1) // West
  {
    next_i = Outlet_i;
    next_j = Outlet_j-1;
    cout << "The closest node is to the west" << endl;
  }
  else if (PerimeterNodes[Outlet_i-1][Outlet_j] == 1) // North
  {
    next_i = Outlet_i-1;
    next_j = Outlet_j;
    cout << "The closest node is to the north" << endl;
  }
  else if (PerimeterNodes[Outlet_i][Outlet_j+1] == 1)  // east
  {
    next_i = Outlet_i;
    next_j = Outlet_j+1;
    cout << "The closest node is to the east" << endl;
  }
  else if (PerimeterNodes[Outlet_i+1][Outlet_j] == 1)  // south
  {
    next_i = Outlet_i+1;
    next_j = Outlet_j;
    cout << "The closest node is to the south" << endl;
  }
  else if (PerimeterNodes[Outlet_i-1][Outlet_j-1] == 1) // northwest
  {
    next_i = Outlet_i-1;
    next_j = Outlet_j-1;
    cout << "The closest node is to the NW" << endl;
  }
  else if (PerimeterNodes[Outlet_i-1][Outlet_j+1] == 1) // northeast
  {
    next_i = Outlet_i-1;
    next_j = Outlet_j+1;
    cout << "The closest node is to the NE" << endl;
  }
  else if (PerimeterNodes[Outlet_i+1][Outlet_j+1] == 1) // southeast
  {
    next_i = Outlet_i+1;
    next_j = Outlet_j+1;
    cout << "The closest node is to the SE" << endl;
  }
  else if (PerimeterNodes[Outlet_i-1][Outlet_j+1] == 1) // southwest
  {
    next_i = Outlet_i-1;
    next_j = Outlet_j+1;
    cout << "The closest node is to the SW" << endl;
  }
  else
  {
    cout << "None of these were perimeter nodes, oops" << endl;
  }

  // push back the next node to the sorted node vector
  int next_node = FlowInfo.retrieve_node_from_row_and_column(next_i, next_j);
  sorted_nodes.push_back(next_node);

  bool reached_outlet = false;
  int this_i, this_j;
  // now start at the outlet node and find the nearest perimeter node.
  while (reached_outlet == false)
  {
    // start at the next node and find the one with the closest distance that
    // hasn't already been visited
    this_i = next_i;
    this_j = next_j;
    VisitedBefore[this_i][this_j] = 1;

    vector<float> Distances(8, 100); // distances to each node in the order N, NE, E, SE, S, SW, W, NW

    if (PerimeterNodes[this_i-1][this_j] == 1) // north
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i-1][this_j] == 0)
        {
          // get the distance
          Distances[0] = DataResolution;
        }
      }
      else
      {
        if (VisitedBefore[this_i-1][this_j] != 1)
        {
          // get the distance
          Distances[0] = DataResolution;
        }
      }
    }
    if (PerimeterNodes[this_i-1][this_j+1] == 1) // northeast
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i-1][this_j+1] == 0)
        {
          // get the distance
          Distances[1] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
      else
      {
        if (VisitedBefore[this_i-1][this_j+1] != 1)
        {
          // get the distance
          Distances[1] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
    }
    if (PerimeterNodes[this_i][this_j+1] == 1) // east
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i][this_j+1] == 0)
        {
          // get the distance
          Distances[2] = DataResolution;
        }
      }
      else
      {
        if (VisitedBefore[this_i][this_j+1] != 1)
        {
          // get the distance
          Distances[2] = DataResolution;
        }
      }
    }
    if (PerimeterNodes[this_i+1][this_j+1] == 1) // southeast
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i+1][this_j+1] == 0)
        {
          // get the distance
          Distances[3] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
      else
      {
        if (VisitedBefore[this_i+1][this_j+1] != 1)
        {
          // get the distance
          Distances[3] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
    }
    if (PerimeterNodes[this_i+1][this_j] == 1) //south
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i+1][this_j] == 0)
        {
          // get the distance
          Distances[4] = DataResolution;
        }
      }
      else
      {
        if (VisitedBefore[this_i+1][this_j] != 1)
        {
          // get the distance
          Distances[4] = DataResolution;
        }
      }
    }
    if (PerimeterNodes[this_i+1][this_j-1] == 1) // southwest
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i+1][this_j-1] == 0)
        {
          // get the distance
          Distances[5] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
      else
      {
        if (VisitedBefore[this_i+1][this_j-1] != 1)
        {
          // get the distance
          Distances[5] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
    }
    if (PerimeterNodes[this_i][this_j-1] == 1) // west
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i][this_j-1] == 0)
        {
          // get the distance
          Distances[6] = DataResolution;
        }
      }
      else
      {
        if (VisitedBefore[this_i][this_j-1] != 1)
        {
          // get the distance
          Distances[6] = DataResolution;
        }
      }
    }
    if (PerimeterNodes[this_i-1][this_j-1] == 1) // northwest
    {
      if (first_node == true)
      {
        if (VisitedBefore[this_i-1][this_j-1] == 0)
        {
          // get the distance
          Distances[7] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
      else
      {
        if (VisitedBefore[this_i-1][this_j-1] != 1)
        {
          // get the distance
          Distances[7] = sqrt(DataResolution*DataResolution+DataResolution*DataResolution);
        }
      }
    }

    // now search the vector of distances for the one with the smallest distance.
    int min_idx = distance(Distances.begin(), min_element(Distances.begin(), Distances.end()));
    if ( min_idx == 0 )
    {
      next_i = this_i-1;
      next_j = this_j;
    }
    else if ( min_idx == 1 )
    {
      next_i = this_i-1;
      next_j = this_j+1;
    }
    else if ( min_idx == 2 )
    {
      next_i = this_i;
      next_j = this_j+1;
    }
    else if ( min_idx == 3 )
    {
      next_i = this_i+1;
      next_j = this_j+1;
    }
    else if ( min_idx == 4 )
    {
      next_i = this_i+1;
      next_j = this_j;
    }
    else if ( min_idx == 5 )
    {
      next_i = this_i+1;
      next_j = this_j-1;
    }
    else if ( min_idx == 6 )
    {
      next_i = this_i;
      next_j = this_j-1;
    }
    else if ( min_idx == 7 )
    {
      next_i = this_i-1;
      next_j = this_j-1;
    }
    // push back the node to the sorted vector
    next_node = FlowInfo.retrieve_node_from_row_and_column(next_i, next_j);
    if (next_node == outlet_node)
    {
      reached_outlet = true;
      cout << "You've reached the outlet, wooohooo" << endl;
      int this_node = FlowInfo.retrieve_node_from_row_and_column(this_i, this_j);
      sorted_nodes.push_back(this_node);
    }
    else if ( Distances[min_idx] == 100)
    {
      VisitedBefore[next_i][next_j] = 1;
      int last_node = sorted_nodes.back();
      int last_i, last_j;;
      FlowInfo.retrieve_current_row_and_col(last_node, last_i, last_j);
      next_i = last_i;
      next_j = last_j;
      if (PerimeterNodes[next_i][next_j] != 1)
      {
        cout << "This isn't even a perimeter node. WHAT THE FUCK" << endl;
      }
      break;
    }
    else
    {
      int this_node = FlowInfo.retrieve_node_from_row_and_column(this_i, this_j);
      sorted_nodes.push_back(this_node);
    }
    first_node = false;
  }

  return sorted_nodes;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the four different hillslope length measurements for the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_all_HillslopeLengths(LSDFlowInfo& FlowInfo, LSDRaster& HillslopeLengths, LSDRaster& Slope, LSDRaster& DinfArea, float log_bin_width, int SplineResolution, float bin_threshold){

  set_HillslopeLength_HFR(FlowInfo, HillslopeLengths);
  set_HillslopeLengths_Boomerang(Slope, DinfArea, FlowInfo, log_bin_width, SplineResolution, bin_threshold);

  if (DrainageDensity != NoDataValue){
    set_HillslopeLength_Density();
  }
  else{
    cout << "\nDrainage Density has not been set, so the hillslope length cannot be set." << endl;
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set all of the basin parameters with one call.
//
// Runs polyfit to get the elevation derivatives, so can be quite memory intensive. Method
// calls all the setters one by one, to populate all the basin parameters. So a
// basin can be created and all it's properties set with 2 calls. The erosion rates have default
// parameters of -9999 as these are rarely used variables.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_All_Parameters(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& CHT, LSDIndexRaster& StreamNetwork,
                                  LSDRaster& HillslopeLengths, LSDRaster& Relief, float window_radius, float log_bin_width,
                                  int SplineResolution, float bin_threshold, float CriticalSlope, float CosmoErosionRate,
                                  float OtherErosionRate){


  //surface fitting
  vector<int> raster_selection;

  raster_selection.push_back(0);
  raster_selection.push_back(1); //slope
  raster_selection.push_back(1); //aspect
  raster_selection.push_back(1); //curvature
  raster_selection.push_back(1); //plan curvature
  raster_selection.push_back(1); //profile curvature
  raster_selection.push_back(0);
  raster_selection.push_back(0);

  vector<LSDRaster> Surfaces = Elevation.calculate_polyfit_surface_metrics(window_radius, raster_selection);
  LSDRaster TotalCurv = Surfaces[3];
  LSDRaster ProfileCurv = Surfaces[5];
  LSDRaster PlanCurv = Surfaces[4];
  LSDRaster Aspect = Surfaces[2];
  LSDRaster Slope = Surfaces[1];
  LSDRaster DinfArea = Elevation.D_inf_units();

  set_SlopeMean(FlowInfo, Slope);
  set_ElevationMean(FlowInfo, Elevation);
  set_ReliefMean(FlowInfo, Relief);
  set_PlanCurvMean(FlowInfo, PlanCurv);
  set_ProfileCurvMean(FlowInfo, ProfileCurv);
  set_TotalCurvMean(FlowInfo, TotalCurv);
  set_PlanCurvMax(FlowInfo, PlanCurv);
  set_ProfileCurvMax(FlowInfo, ProfileCurv);
  set_TotalCurvMax(FlowInfo, TotalCurv);
  set_CHTMean(FlowInfo, CHT);
  set_AspectMean(FlowInfo, Aspect);
  set_FlowLength(StreamNetwork, FlowInfo);
  set_DrainageDensity();
  set_all_HillslopeLengths(FlowInfo, HillslopeLengths, Slope, DinfArea, log_bin_width, SplineResolution, bin_threshold);
  set_Perimeter(FlowInfo);
  set_EStar_RStar(CriticalSlope);
  set_HilltopPx(FlowInfo, CHT);
  set_CosmoErosionRate(CosmoErosionRate);
  set_OtherErosionRate(OtherErosionRate);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write integer basin parameters into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_HilltopPx(LSDFlowInfo& FlowInfo, LSDRaster Hilltops){

  int i;
  int j;
  int HilltopCount = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    //only count hilltop pixels
    if (Hilltops.get_data_element(i,j) != NoDataValue){
      ++HilltopCount;
    }

  }

  HilltopPx = HilltopCount;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write integer basin parameters into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDBasin::write_integer_data_to_LSDIndexRaster(int Param, LSDFlowInfo FlowInfo)
{

  int i;
  int j;
  Array2D<int> Output(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Param;
  }

  LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output, GeoReferencingStrings);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write real basin parameters into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::write_real_data_to_LSDRaster(float Param, LSDFlowInfo FlowInfo)
{

  int i;
  int j;
  Array2D<float> Output(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Param;
  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output, GeoReferencingStrings);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cookie cut data from an LSDRaster into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::write_raster_data_to_LSDRaster(LSDRaster& Data, LSDFlowInfo& FlowInfo)
{

  int i;
  int j;
  Array2D<float> Output(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Data.get_data_element(i,j);
  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output, GeoReferencingStrings);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cookie cut data from an LSDIndexRaster into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDBasin::write_raster_data_to_LSDIndexRaster(LSDIndexRaster Data, LSDFlowInfo FlowInfo){

  int i;
  int j;
  Array2D<int> Output(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Data.get_data_element(i,j);
  }

  LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output, GeoReferencingStrings);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to check whether a node is in the basin
// FJC 21/02/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDBasin::is_node_in_basin(int test_node)
{
  int node_checker = 0;
  for (int i =0; i < int(BasinNodes.size()) && node_checker == 0; i++)
  {
    if (test_node == BasinNodes[i])
    {
      node_checker = 1;
    }
  }
  return node_checker;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis function test if the basins are from the same base DEM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDBasin::are_basins_from_same_base_DEM(LSDBasin& other)
{
  map<string,string> DRS = other.get_GeoReferencingStrings();
  string mi_str_key = "ENVI_map_info";
  string cs_str_key = "ENVI_coordinate_system";

  bool is_georef_and_dimensions_same;

  if (NRows == other.get_NRows() &&
      NCols == other.get_NCols() &&
      XMinimum == other.get_XMinimum() &&
      YMinimum == other.get_YMinimum() &&
      DataResolution == other.get_DataResolution() &&
      NoDataValue == other.get_NoDataValue() &&
      GeoReferencingStrings[mi_str_key] == DRS[mi_str_key] &&
      GeoReferencingStrings[cs_str_key] == DRS[cs_str_key])
  {
    is_georef_and_dimensions_same = true;
  }
  else
  {
    is_georef_and_dimensions_same = false;
  }

  return is_georef_and_dimensions_same;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to check whether a node is in the basin
// FJC 21/02/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDBasin::is_this_a_subbasin(LSDBasin& DifferentBasin)
{
  bool is_this_a_subbasin = false;

  // first check to see if the basins are from the same DEM
  if (are_basins_from_same_base_DEM(DifferentBasin))
  {
    // if the second basin is a subbasin, the outlet will be within the base
    // DEM
    int outlet_node = DifferentBasin.get_Outlet_node();
    int in_basin = is_node_in_basin(outlet_node);
    if (in_basin == 1)
    {
      is_this_a_subbasin = true;
    }

  }

  return is_this_a_subbasin;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to provide an alternative (integer) index (e.g. lithology) based on
// majority
// DTM 26/08/2015
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_AlternativeIndex(LSDFlowInfo& FlowInfo, LSDIndexRaster& AltIndex)
{
  int max = 0;
  int i,j;
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    if(AltIndex.get_data_element(i,j)>max) max = AltIndex.get_data_element(i,j);
  }
  //  for(int i=0; i<AltIndex.get_NRows(); ++i){
  //    for(int j=0; j<AltIndex.get_NCols(); ++j){
  //      if(AltIndex.get_data_element(i,j)>max) max=AltIndex.get_data_element(i,j);
  //    }
  //  }

  vector<int> indices;
  vector<int> counts;
  for(int index = 0; index < max+1; ++index){
    indices.push_back(index);
    counts.push_back(0);
  }
  //cout << max << " " << indices.back() << endl;
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    ++counts[AltIndex.get_data_element(i,j)];
  }
  vector<size_t> index_map;
  matlab_int_sort(counts, counts, index_map);
  matlab_int_reorder(indices,index_map,indices);
  AlternativeIndex = indices.back();
  //cout << "yoo hoo " << AlternativeIndex << endl;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to remove any hilltop curvature values that are not internal to the
// basin
// FJC 19/03/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::keep_only_internal_hilltop_curvature(LSDRaster hilltop_curvature, LSDFlowInfo FlowInfo)
{
  //Array2D<float> CHT_array(NRows, NCols, NoDataValue);
  //Array2D<int> perimeter(NRows, NCols, 0);
  //Array2D<int> ridge_pixels(NRows, NCols,0);
  Array2D<float> CHT_array = hilltop_curvature.get_RasterData();
  //cout << "Number of perimeter pixels: " << Perimeter_i.size() << endl;
  for (int q =0; q < int(Perimeter_i.size()); q++){
    int row_perim = Perimeter_i[q];
    int col_perim = Perimeter_j[q];
    if (row_perim != NoDataValue && col_perim != NoDataValue){
      //perimeter[row_perim][col_perim] = 1;
      CHT_array[row_perim][col_perim]=NoDataValue;
    }
  }

  //  for (int row = 0; row < NRows; row++)
  //{
  //for (int col = 0; col < NCols; col++)
  //{
  //  float curvature = hilltop_curvature.get_data_element(row, col);
  //  if (curvature != NoDataValue)
  //  {
  //    ridge_pixels[row][col] = 1;
  //  }
  //}
  //}


  //for (int row = 0; row < NRows; row++)
  //{
  //for (int col = 0; col < NCols; col++)
  //{
  //  if (perimeter[row][col] == 0 || ridge_pixels[row][col] == 0)
  //  {
  //    float curvature = hilltop_curvature.get_data_element(row, col);
  //    CHT_array[row][col] = curvature;
  //  }
  //}
  //}
  LSDRaster CHT(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, CHT_array,GeoReferencingStrings);

   return CHT;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function trims a padded raster to the basin
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::TrimPaddedRasterToBasin(int padding_pixels, LSDFlowInfo& FlowInfo,
                                            LSDRaster& Raster_Data)
{
  int max_row = -1;
  int max_col = -1;
  int min_row = 1e8;
  int min_col = 1e8;

  int row,col;

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    // get the maximum and minimum columns.
    if(row > max_row)
    {
      max_row = row;
    }
    if(col > max_col)
    {
      max_col = col;
    }
    if(row < min_row)
    {
      min_row = row;
    }
    if(col < min_col)
    {
      min_col = col;
    }
  }

  // pad the rows and columms
  min_row = min_row-padding_pixels;
  min_col = min_col-padding_pixels;
  max_row = max_row+padding_pixels;
  max_col = max_col+padding_pixels;

  // restrict the rows and columns to the size of the DEM
  if (min_row < 0)
  {
    min_row = 0;
  }
  if (min_col < 0)
  {
    min_col = 0;
  }
  if (max_row > NRows-1)
  {
    max_row = NRows-1;
  }
  if (max_col > NCols -1)
  {
    max_col = NCols-1;
  }

  // now generate the new DEM
  // create new row and col sizes taking account of zero indexing
  int new_row_dimension = (max_row-min_row) + 1;
  int new_col_dimension = (max_col-min_col) + 1;

  Array2D<float>TrimmedData(new_row_dimension, new_col_dimension, NoDataValue);

  //loop over min bounding rectangle and store it in new array of shape new_row_dimension x new_col_dimension
  int TrimmedRow;
  int TrimmedCol;

  for (int row = min_row; row < max_row+1; ++row)
  {
    for(int col = min_col; col < max_col+1; ++col)
    {
      TrimmedRow = row-min_row;
      TrimmedCol = col-min_col;

      TrimmedData[TrimmedRow][TrimmedCol] = Raster_Data.get_data_element(row,col);
    }
  }

  //calculate lower left corner coordinates of new array
  float new_XLL = (min_col * DataResolution) + XMinimum;
  float new_YLL = YMinimum + ((NRows - max_row - 1) * DataResolution);
  float YMax = new_YLL + (new_row_dimension* DataResolution);


  LSDRaster TrimmedRaster(new_row_dimension, new_col_dimension, new_XLL,
                          new_YLL, DataResolution, NoDataValue, TrimmedData,
                          GeoReferencingStrings);

  cout << "Trimming DataResolution: " << DataResolution << endl;


  // update the georeferencing
  TrimmedRaster.Update_GeoReferencingStrings(new_XLL,YMax);

  // need to do this a second time to make sure the float data is passed to the
  // header file
  TrimmedRaster.Update_GeoReferencingStrings();

  return TrimmedRaster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This adds basins to an LSDIndexRaster
// Smaller basins overwrite larger basins
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::add_basin_to_LSDIndexRaster(LSDIndexRaster& basin_raster,
                                           LSDFlowInfo& FlowInfo,
                                           map<int,int>& drainage_of_other_basins,
                                           int this_basin_index)
{
  int row, col;
  int existing_basin_index;
  int this_basin_pixel_area = int(BasinNodes.size());

  // add the drainage area of the basin
  drainage_of_other_basins[this_basin_index] = this_basin_pixel_area;

  // make sure that the index raster is the same size as the basin
  // check_georeferencing
  int RNR = basin_raster.get_NRows();
  int RNC = basin_raster.get_NCols();
  int NDV = basin_raster.get_NoDataValue();
  if (RNR != NRows || RNC != NCols)
  {
    cout << "LSDBasin::Add_basin_to_LSDIndexRaster, failed to add basins" << endl
         << "since LSDIndexRaster does not have same georeferncing as basin." << endl;
  }
  else
  {
    for (int q = 0; q < int(BasinNodes.size()); ++q)
    {
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      // check to see if there is already data there
      if( basin_raster.get_data_element(row,col) == NDV)
      {
        basin_raster.set_data_element(row,col,this_basin_index);
      }
      else
      {
        // if there is already data there, the smaller basin overrwites the larger
        // basin. You will be able to see the larger basins since the smaller
        // basins are nested within the larger basins
        existing_basin_index = basin_raster.get_data_element(row,col);

        if(drainage_of_other_basins.find(existing_basin_index) == drainage_of_other_basins.end())
        {
          cout << "LSDBasin::Add_basin_to_LSDIndexRaster, Something has gone wrong." << endl
               << "The existing basin is not in the basin map. " << endl
               << "Overwriting data." << endl;
          basin_raster.set_data_element(row,col,this_basin_index);
        }
        else if (this_basin_pixel_area < drainage_of_other_basins[existing_basin_index])
        {
          basin_raster.set_data_element(row,col,this_basin_index);
        }
      }
    }
  }

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Method to merge a vector of LSDRaster basins generated using LSDBasin into
// a single LSDRaster for visualisation.
//
// 50% less computationally expesnive than the old method, but still very
// inefficient. Does not test for overlaps in the data, will simply overwrite
// so that the last value to occupy a cell will be written.
//
// SWDG 07/12/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::Merge_Basins(vector<LSDRaster> Basins)
{

  Array2D<float> Output(NRows,NCols,NoDataValue);

  for (int q = 0; q < int(Basins.size()); ++q){

    for (int i = 0; i < NRows; ++i){
      for (int j = 0; j < NCols; ++j){
        if (Basins[q].get_data_element(i,j) != NoDataValue){
          Output[i][j] = Basins[q].get_data_element(i,j);
        }
      }
    }

  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output, GeoReferencingStrings);

  return OutputRaster;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Count the number of each unique lithology value contain in the basin from a topographic raster
//
// take a lithologic raster and a flowinfo raster in argument
// BG 17/09/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,int> LSDBasin::count_unique_values_from_litho_raster(LSDIndexRaster& litho, LSDFlowInfo& topo)
{
  //cout << "I am now proceeding to counting the lithology per basin you gave me." <<endl;
  // First get the unique values in litho
  vector<int> values = litho.get_list_of_values();
  // Initialisation of the map
  map<int,int> lithost;
  for(size_t i = 0; i<values.size(); i++)
  {
    lithost[values[i]] = 0;
    //cout << values[i] << endl;

  }

  // Initializing the NoData as well
  // map initialized with all the values of lithology present on the map

  // Now counting through the nodes

  //int trow,tcol,tvalue; // temporary variables to store the row, col and value - OBSOLETE
  int tvalue,tcbbtf; // temporary variables to store the row, col and value
  float teasting, tnorthing; // temporary variables to store the easting/northing.
  // OPTIMIZITION THOUGH:
  // rather than calculating the easting/northing/row/col each times, maybe cropping the
  // litho raster to the extent of the toporaster to avoid this??
  // END OF OPTIMIZATION THOUGH
  for(size_t j = 0; j<BasinNodes.size();j++)
  {
    // getting the easting_northing for this node
    topo.get_x_and_y_from_current_node(BasinNodes[j],teasting,tnorthing);
    // getting the corresponding row-col for the litho raster
    //litho.get_row_and_col_of_a_point(teasting,tnorthing,trow,tcol);
    // implementing the counting
    tvalue = litho.get_value_of_point(teasting,tnorthing);
    tcbbtf = lithost[tvalue];
    lithost[tvalue] = tcbbtf+1;
    // ALRIGHT TO DO LATER:
    // OPTIMIZE THE STRUCTURE TO LOOP ONCE THROUGH THE NODES AND IMPLEMENT DIRECTLY THE map
    // WRITE THE MAP IN A CSV FILE FROM THE driver
    // ADD A PARAMETER PERCENTAGE
  }

  // //DEBUG STATMENT TO KEEP IN CASE
  // for(map<int,int>::iterator it = lithost.begin(); it!=lithost.end(); ++it)
  // {
  //   cout << "Litho " << it->first << " is represented " << it->second << " times." << endl;
  // }


  //cout << "I am done with counting the lithology per basins" <<endl;
  //return the results
  return lithost;
}


bool LSDBasin::is_adjacent(LSDBasin& DifferentBasin, LSDFlowInfo& flowpy)
{
  // Check if the perimeter has been calculated, and set it if not
  if(Perimeter_nodes.size() == 1){set_Perimeter(flowpy);}
  //
  //cout << "I will check if the two basins are adjacents" << endl;
  // this is a comment
  // Perimeter nodes of the other basin + proceeding to a test if the other basin's perimeter has been calculated
  vector<int> DB_perimeter_nodes = DifferentBasin.get_Perimeter_nodes();
  if(DB_perimeter_nodes.size() == 1){DifferentBasin.set_Perimeter(flowpy); DB_perimeter_nodes =  DifferentBasin.get_Perimeter_nodes();}
  //

  // cout << Perimeter_nodes.size() << "|||||" << DB_perimeter_nodes.size() << endl;
  bool adjacenty = false;
  int this_row = 0;
  int this_col = 0;
  int tested_node = 0;
  //Preprocessing looping
  static const int arr[] = {-1,0,1};
  vector<int> mongoose (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  // Now looping through the neiboor nodes of each perimeter nodes of the first basin, then check if it is in the second one. stop the process if one node is found
  for(size_t flude = 0; flude<Perimeter_nodes.size() && adjacenty == false ;flude++)
  {
    flowpy.retrieve_current_row_and_col(Perimeter_nodes[flude], this_row, this_col);
    for(size_t impala = 0; impala < mongoose.size() && adjacenty == false; impala ++)
    {
      for(size_t peccary = 0; peccary < mongoose.size() && adjacenty == false; peccary ++)
        if((this_row + mongoose[impala]>=0)&&(this_row + mongoose[impala]< flowpy.get_NRows())&&(this_col+ mongoose[peccary] >=0)&&(this_col+ mongoose[peccary]<flowpy.get_NCols()))
        {
          tested_node = flowpy.retrieve_node_from_row_and_column(this_row + mongoose[impala],this_col+ mongoose[peccary]);
          if(find(DB_perimeter_nodes.begin(), DB_perimeter_nodes.end(), tested_node) != DB_perimeter_nodes.end()){adjacenty=true;}
        }
    }

  }
  return adjacenty;
}

/// @brief merge and contour the perimeter from a vector of adjacent basins
/// @detail WARNING There may be 1 pixel-size holes in the perimeter.
/// @param vector of LSDBasin objects
/// @return vector of node indices of the new perimeter
/// @author BG
/// @date 10/10/17
vector<int> LSDBasin::merge_perimeter_nodes_adjacent_basins(vector<LSDBasin> budgerigar,LSDFlowInfo& flowpy)
{
  // getting and checking the perimeter for the the first basin
  vector<int> out_perimeter;
  vector<int> first_perimeter;

  // other variables
  vector<int> temp_perimeter;
  bool adjacenty = false;
  int this_row = 0;
  int this_col = 0;
  int tested_node = 0;


  for(int basilisk = 0; basilisk< int(budgerigar.size());basilisk++)
  {
    // Getting and checking the existence of the perimeter node
    first_perimeter = budgerigar[basilisk].get_Perimeter_nodes();
    if(first_perimeter.size() == 1){budgerigar[basilisk].set_Perimeter(flowpy);first_perimeter = budgerigar[basilisk].get_Perimeter_nodes();}

    // now check the nodes around each perimeter nodes to see if it is adjacent to the other perimiter (diagonal excluded) and select it if not
    for(int itb = 0; itb< int(budgerigar.size());itb++)
    {
      if(itb != basilisk)
      {
        temp_perimeter = budgerigar[itb].get_Perimeter_nodes();
        if(temp_perimeter.size() == 1){budgerigar[itb].set_Perimeter(flowpy);temp_perimeter = budgerigar[itb].get_Perimeter_nodes();}
        for(size_t flude = 0; flude<first_perimeter.size();flude++)
        {
          flowpy.retrieve_current_row_and_col(first_perimeter[flude], this_row, this_col);
          if(this_row>0 && this_row<flowpy.get_NRows()-1 && this_col>0 && this_col<flowpy.get_NCols()-1)
          {
            tested_node = flowpy.retrieve_node_from_row_and_column(this_row,this_col+1);
            if(find(temp_perimeter.begin(), temp_perimeter.end(), tested_node) != temp_perimeter.end()){adjacenty=true;}
            tested_node = flowpy.retrieve_node_from_row_and_column(this_row,this_col-1);
            if(find(temp_perimeter.begin(), temp_perimeter.end(), tested_node) != temp_perimeter.end()){adjacenty=true;}
            tested_node = flowpy.retrieve_node_from_row_and_column(this_row+1,this_col);
            if(find(temp_perimeter.begin(), temp_perimeter.end(), tested_node) != temp_perimeter.end()){adjacenty=true;}
            tested_node = flowpy.retrieve_node_from_row_and_column(this_row-1,this_col);
            if(find(temp_perimeter.begin(), temp_perimeter.end(), tested_node) != temp_perimeter.end()){adjacenty=true;}
          }

          if(!adjacenty){out_perimeter.push_back(temp_perimeter[flude]);}
          adjacenty = false;
        }
      }

    }
  }

  return out_perimeter;
}



/// @brief detect the source nodes in a pixel window around a perimeter, for instance a basin  perimeter
/// @detail It needs a sequence of nodes where it will loop around and gather all the source nodes encountered.
/// @param vector of nodes, Flowinfo object and a JunctionNetwork object and a number of pixel for the window.
/// @return vector of node indices of the new perimeter
/// @author BG
/// @date 10/10/17
vector<int> LSDBasin::get_source_node_from_perimeter(vector<int> perimeter, LSDFlowInfo& flowpy, LSDJunctionNetwork& junky, int pixel_window)
{

  //First, creating a square-shaped mangoose vector (vector that host the pixel window parameter to loop through)
  vector<int> mangoose;
  mangoose.push_back(0);
  for(size_t coati = 1; coati<size_t(pixel_window); coati++)
  {
    mangoose.push_back(coati);
    mangoose.push_back(-coati);
  }
  // mangoose is ready

  //get all the sources
  vector<int> all_sources = junky.get_SourcesVector();

  // creating an empty output vector
  vector<int> selected_sources_nodes;

  //creating the temp variables
  int that_row = 0;
  int that_col = 0;
  int tested_node = 0;

  //loop through the perimeter and neighbooring nodes in the window
  for(size_t coati = 0; coati<perimeter.size();coati++)
  {
    flowpy.retrieve_current_row_and_col(perimeter[coati],that_row,that_col);
    for(size_t hogger =0; hogger<mangoose.size(); hogger++)
    {
      for(size_t hoggest =0; hoggest<mangoose.size(); hoggest++)
      {
        if((that_row + mangoose[hogger]>=0)&&(that_row + mangoose[hogger]< flowpy.get_NRows())&&(that_col+ mangoose[hoggest] >=0)&&(that_col+ mangoose[hoggest]<flowpy.get_NCols()))
        {
          tested_node = flowpy.retrieve_node_from_row_and_column(that_row+mangoose[hogger],that_col+mangoose[hoggest]);
          if(find(all_sources.begin(), all_sources.end(), tested_node) != all_sources.end()){selected_sources_nodes.push_back(tested_node);}
        }


      }
    }
  }

  // now selecting the unique values

  sort( selected_sources_nodes.begin(), selected_sources_nodes.end() );
  selected_sources_nodes.erase( unique( selected_sources_nodes.begin(), selected_sources_nodes.end() ), selected_sources_nodes.end() );

  // it should work
  return selected_sources_nodes;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// THis compile some metrics for each side of a drainage divide across a serie of nodes: min,max,mean,median,std dev ...
// You just have to feed it with a vector of nodes and a template raster. This last can be elevation, slope, a normalized swath and so on
// return a map where the key is a string like "min_in" or ",min_out" as well as "n_node_in"... Ill list it when it will be done.
map<string,float> LSDBasin::get_metrics_both_side_divide(LSDRaster& rasterTemplate, LSDFlowInfo& flowpy, vector<int>& nodes_to_test, map<int,bool>& raster_node_basin)
{

  map<string,float> mapout;
  vector<int> nodes_in, nodes_out;
  vector<float> nodes_in_v, nodes_out_v, nodes_to_test_v;
  int n_nodes_in =0, n_nodes_out = 0,row = 0, col = 0, nnd =0;


  // Now testing which nodes are in  and out of the basin
  for(vector<int>::iterator pudu_deer = nodes_to_test.begin(); pudu_deer != nodes_to_test.end(); pudu_deer++)
  {
    //cout << *pudu_deer << " - metrics calcnode" << endl;
    // get row/col
    if(*pudu_deer != NoDataValue && raster_node_basin[*pudu_deer] != NoDataValue)
    {
      flowpy.retrieve_current_row_and_col(*pudu_deer,row,col);
      if(rasterTemplate.get_data_element(row,col) != NoDataValue )
      {

        if(raster_node_basin[*pudu_deer])
        {
          nodes_in.push_back(*pudu_deer);
          nodes_in_v.push_back(rasterTemplate.get_data_element(row,col));
          n_nodes_in++;
        }
        else
        {
          nodes_out.push_back(*pudu_deer);
          nodes_out_v.push_back(rasterTemplate.get_data_element(row,col));
          n_nodes_out++;
        }
        nodes_to_test_v.push_back(rasterTemplate.get_data_element(row,col));
      }
      else
      {
        nnd++;
      }
    }
    else
    {
      nnd++;
    }
  }
  // cout << "in: " << n_nodes_in << " out: " << n_nodes_out << " nodata: " << nnd << "/" << nodes_to_test.size() << endl;
  mapout["n_nodes_in"] =  (float)n_nodes_in;
  mapout["n_no_data"] = (float)nnd;
  mapout["n_nodes_out"] = (float)n_nodes_out;
  // done

  // compiling the stat using stat_tools
  mapout["mean"] = get_mean_ignore_ndv(nodes_to_test_v, (float)NoDataValue);
  mapout["mean_in"] = get_mean_ignore_ndv(nodes_in_v, (float)NoDataValue);
  mapout["mean_out"] = get_mean_ignore_ndv(nodes_out_v, (float)NoDataValue);
  mapout["median"] = get_median(nodes_to_test_v, (float)NoDataValue);
  mapout["median_in"] = get_median(nodes_in_v, (float)NoDataValue);
  mapout["median_out"] = get_median(nodes_out_v, (float)NoDataValue);
  mapout["max"] = Get_Maximum(nodes_to_test_v, (float)NoDataValue);
  mapout["max_in"] = Get_Maximum(nodes_in_v, (float)NoDataValue);
  mapout["max_out"] = Get_Maximum(nodes_out_v, (float)NoDataValue);
  mapout["min"] = Get_Minimum(nodes_to_test_v, (float)NoDataValue);
  mapout["min_in"] = Get_Minimum(nodes_in_v, (float)NoDataValue);
  mapout["min_out"] = Get_Minimum(nodes_out_v, (float)NoDataValue);
  mapout["StdDev"] = get_standard_deviation(nodes_to_test_v, mapout["mean"], (float)NoDataValue);
  mapout["StdDev_in"] = get_standard_deviation(nodes_in_v, mapout["mean_in"], (float)NoDataValue);
  mapout["StdDev_out"] = get_standard_deviation(nodes_out_v, mapout["mean_out"], (float)NoDataValue);

  return mapout;


}


// This move a square window along the drainage divide and compute the statistics in and out of the basin. Used to compare
// the slope/elevetion/lithology/...
// I'll code a non square version at some point!
// BG - work in progress (on this topic, not on myself, or maybe, what does that even mean?)

void LSDBasin::square_window_stat_drainage_divide(LSDRaster& rasterTemplate, LSDFlowInfo& flowpy, int size_window)
{
  if(DD_preprocessed == false)
  {
    cout << "You need to use the function preprocess_DD_metrics(LSDFlowInfo FlowInfoCorrespondingToThisBasin) before being able to use this function" << endl;
    exit(EXIT_FAILURE);
  }

  // If you call this function, you want to reinitialize the potentially previously calculated map
  stats_around_perimeter_window.clear();


  int row = 0, col = 0,perimeter_index = 0, this_node = 0;
  vector<int> nodes_to_test;
  map<int, bool> raster_node_basin;

  for(int i = 0; i<rasterTemplate.get_NRows();i++)
  //for(map<int,vector<float> >::iterator claus = BasinNodesMapOfXY.begin(); claus != BasinNodesMapOfXY.end(); claus ++)
  {
    for(int j = 0; j<rasterTemplate.get_NCols();j++)
    {
      raster_node_basin[flowpy.retrieve_node_from_row_and_column(i,j)] = false;
    }
   //raster_node_basin[flowpy.get_node_index_of_coordinate_point(claus->second[0],claus->second[1])] = false;
  }

  for(map<int,vector<float> >::iterator oisture = BasinNodesMapOfXY.begin(); oisture != BasinNodesMapOfXY.end();oisture++)
  {
    this_node = flowpy.get_node_index_of_coordinate_point(oisture->second[0],oisture->second[1]);
    raster_node_basin[this_node] = true;
  }

  // cout << "basination done" << endl;


  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Not working after that




  // loop through the perimeter
  for(map<int,vector<float> >::iterator noodle = DD_map.begin(); noodle != DD_map.end(); noodle++)
  {
    // cout << "processing node " << noodle -> first << endl;
    // get the row/col info
    this_node = flowpy.get_node_index_of_coordinate_point(noodle->second[0],noodle->second[1]);
    flowpy.retrieve_current_row_and_col(this_node,row,col);
    // cout << row << "¬¬" << col << endl;
    //loop through the window to gather the wanted node index
    for(int i = row - size_window ; i< row + size_window; i++)
    {
      for(int j = col - size_window ; j< col + size_window; j++)
      {
        // check if the row/col are within the raster
        if(i<rasterTemplate.get_NRows() && j<rasterTemplate.get_NCols() && i >= 0 && j >= 0 )
        {
          nodes_to_test.push_back(flowpy.retrieve_node_from_row_and_column(i,j));
        }
      }
    }
    // ok now I have the number of node to test, I'll just compute the stats
    // cout<< "metrics " << endl;
    if(nodes_to_test.size()>0)
    {
      map<string, float>  temp_map = get_metrics_both_side_divide(rasterTemplate,flowpy, nodes_to_test, raster_node_basin);

      // cout << "done"<< endl;
      temp_map["ID"] = perimeter_index; // I am adding an unique index to then sort the data in python by this index
      temp_map["value_at_DD"] = rasterTemplate.get_data_element(row,col);
      stats_around_perimeter_window[this_node] = temp_map;
      // clearing the vector and start again
    }
    nodes_to_test.clear();
    perimeter_index++;
  }

  // now just processing the distance map
  float dist = 0, last_x = 0 , last_y = 0, this_x = 0, this_y = 0;
  for(vector<int>::iterator titan = Perimeter_nodes.begin(); titan!= Perimeter_nodes.end(); titan ++)
  {
    flowpy.retrieve_current_row_and_col(*titan,row,col);
    flowpy.get_x_and_y_locations(row, col, this_x, this_y);
    if(titan == Perimeter_nodes.begin())
    {
      dist = 0;
    }
    else
    {
      dist = dist + sqrt(pow((this_x-last_x),2)+pow((this_y-last_y),2));
    }
    map_of_dist_perim[*titan] = dist;
    last_x = this_x;
    last_y = this_y;

  }

// Done
  cout << "I have computed statistics around a square window for this basin" << endl;

}


void LSDBasin::write_windowed_stats_around_drainage_divide_csv(string filename, LSDFlowInfo& flowpy)
{
  // these are for extracting element-wise data from the channel profiles.
  cout << "I am now writing your DD stat windowed file:" << endl;
  int this_node, row,col;
  double latitude,longitude, this_x, this_y,last_x,last_y;
  float dist = 0;
  map<string,float> this_map;
  LSDCoordinateConverterLLandUTM Converter;



  // open the data file
  ofstream  data_out;
  data_out.open(filename.c_str());
  data_out << "X,Y,latitude,longitude,distance,value_at_DD,mean,mean_in,mean_out,median,median_in,median_out,StdDev,StdDev_in,StdDev_out,min,min_in,min_out,max,max_in,max_out";
  data_out << endl;


  map<int,map<string,float> >::iterator iter;

  for (iter = stats_around_perimeter_window.begin(); iter != stats_around_perimeter_window.end(); iter++)
  {
    if(iter != stats_around_perimeter_window.begin())
    {
      last_x = this_x;
      last_y = this_y;
    }
    this_node = iter->first;
    this_map = iter->second;
    flowpy.retrieve_current_row_and_col(this_node,row,col);
    flowpy.get_lat_and_long_locations(row, col, latitude, longitude, Converter);
    flowpy.get_x_and_y_locations(row, col, this_x, this_y);

      data_out.precision(9);
      data_out << this_x << ","
                   << this_y << ","
                   << latitude << ","
                   << longitude << ",";
      data_out.precision(5);
      data_out << map_of_dist_perim[iter->first] << ","
               << this_map["value_at_DD"] << ","
               << this_map["mean"] << ","
               << this_map["mean_in"] << ","
               << this_map["mean_out"] << ","
               << this_map["median"] << ","
               << this_map["median_in"] << ","
               << this_map["median_out"] << ","
               << this_map["StdDev"] << ","
               << this_map["StdDev_in"] << ","
               << this_map["StdDev_out"] << ","
               << this_map["min"] << ","
               << this_map["min_in"] << ","
               << this_map["min_out"] << ","
               << this_map["max"] << ","
               << this_map["max_in"] << ","
               << this_map["max_out"]
               << endl;
    }
    data_out.close();
    //"X,Y,latitude,longitude,distance,value_at_DD,mean,mean_in,mean_out,median,median_in,median_out,
    //StdDev,StdDev_in,StdDev_out,min,min_in,min_out,max,max_in,max_out"

}

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This preprocess Drainage Divide metrics to make sure everything is ready to be used for this basin:
// setting the perimeter and calculating various info about it:
// at the moment: x/y coordinates, to deal with later different rasters
// TODO: find a way to sort the perimeter in a vector: each nodes should follow each other
// BG - 26/12/2017 - "hon hon hon" (French Santa, unpublished)
void LSDBasin::preprocess_DD_metrics(LSDFlowInfo flowpy)
{

  // First setting the perimeter
  set_Perimeter(flowpy);
  clean_perimeter(flowpy);

  //organise_perimeter(flowpy);

  // implementing a global map containing vector of info for each nodes atm <x,y> later on distance from origin
  int this_node = 0;
  float this_x = 0, this_y = 0;
  vector<float> info_DD;
  vector<int>::iterator Santa;

  for(Santa = Perimeter_nodes.begin();Santa!=Perimeter_nodes.end();Santa++)
  {
    this_node = *Santa;
    flowpy.get_x_and_y_from_current_node(this_node,this_x,this_y);
    info_DD.push_back(this_x);
    info_DD.push_back(this_y);
    // add some stuffs later

    DD_map[this_node] = info_DD;
    info_DD.clear();
  }

  // Setting a map [nodes_of_basins] = vector[x,y] to make an easier comparison with other raster
  for(Santa = BasinNodes.begin();Santa!=BasinNodes.end();Santa++)
  {
    this_node = *Santa;
    flowpy.get_x_and_y_from_current_node(this_node,this_x,this_y);
    info_DD.push_back(this_x);
    info_DD.push_back(this_y);
    BasinNodesMapOfXY[this_node] = info_DD;
    info_DD.clear();
  }

  // Setting the preprocess checker
  DD_preprocessed = true;
  // done

}

void LSDBasin::organise_perimeter(LSDFlowInfo& flowpy)
{
  // Careful!!! This is a test function, I am definitely trying things here but it might not be ready yet
  vector<int>::iterator YOP = Perimeter_nodes.begin();
  vector<int> row_nodes_to_test, col_nodes_to_test;
  // ### This array will simplify the looping through the perimeter
  static const int arr[] = {-1,0,1};
  vector<int> tester (arr, arr + sizeof(arr) / sizeof(arr[0]) );

  map<int,int> is_done; // check if a node has been processed or not
  int row = 0, col = 0,id = 0, node = 0, n_adj = 0, this_row,this_col;
  float x1 =0 ,x2 = 0, y1 = 0, y2 = 0, last_dist = 0;
  bool ting = true;

  // preprocessing stage to get rid of some points
  print_perimeter_to_csv(flowpy, "/home/boris/Desktop/LSD/capture/sorbas/peritest.csv");


  // Any perimeter nodes should only have 2 neightboors now

  // first, let me indentify the outlet (???), aka the lowest point of the ridge
  node = Perimeter_nodes[0];
  flowpy.retrieve_current_row_and_col(node,row,col);
  is_done[node] = 1;
  map_of_dist_perim[node] = 0;
  Perimeter_nodes_sorted_id[node] = id;
  Perimeter_nodes_sorted.push_back(node);
  // cout << is_done[5260475] << " " << node << " " << 5260475 << endl;
  // exit(EXIT_SUCCESS);


  while(Perimeter_nodes_sorted.size() != 736)
  {

    id++;
    for(vector<int>::iterator YOPL = tester.begin();YOPL!= tester.end() && ting == true ; YOPL++)
    {

      for(vector<int>::iterator LPOY = tester.begin();LPOY != tester.end() && ting == true; LPOY++)
      {
        this_row = row + *YOPL;
        this_col = col + *LPOY;

        if(this_row<flowpy.get_NRows() && this_col<flowpy.get_NCols() && this_row>=0 && this_col >= 0)
        {
          flowpy.get_x_and_y_from_current_node(node,x2,y2);
          // cout << Perimeter_nodes_sorted.size() << endl;
          // cout << node << endl;
          //cout << x2 << " " <<y2 <<endl;

          node = flowpy.retrieve_node_from_row_and_column(this_row,this_col);
          // cout << Perimeter_nodes_map.count(node) << " || " << is_done.count(node) << " || " << x2 << " || " << y2  <<  endl;
          if(Perimeter_nodes_map.count(node) > 0 && is_done.count(node) == 0 && node != NoDataValue && node != -9999 )
          {
            ting = false;
            Perimeter_nodes_sorted.push_back(node);
            is_done[node] = 1;
            Perimeter_nodes_sorted_id[node] = id;
            // distance stuffs
            flowpy.get_x_and_y_from_current_node(Perimeter_nodes_sorted[Perimeter_nodes_sorted.size()-2],x1,y1);
            flowpy.get_x_and_y_from_current_node(node,x2,y2);
            map_of_dist_perim[node] = last_dist + sqrt(pow((x2-x1),2) + pow((y2-y1),2));
            last_dist = map_of_dist_perim[node];
            // cout << map_of_dist_perim[node] << endl;


            flowpy.retrieve_current_row_and_col(node,row,col);


            // debugf
            flowpy.get_x_and_y_from_current_node(node,x2,y2);
            cout << row << " || " << col << endl;
            cout << x2 << " " <<y2 <<endl;
            cout << Perimeter_nodes_sorted.size() << endl;

          }
        }
      }

    }

    ting = true;

  }
  cout << "fupal" << endl;

  // TESTING FUNCTION, DELETE IT AFTERWARDS BORIS!!!!
  Perimeter_nodes = Perimeter_nodes_sorted;
    print_perimeter_to_csv(flowpy, "/home/boris/Desktop/LSD/capture/sorbas/peritest_AFTER.csv");




  // // Now looping from this point while my perimeter_nodes_sorted is not full
  // while(Perimeter_nodes.size()!=Perimeter_nodes_sorted.size())
  // {
  //   // checking all the adjacent nodes (no diagonals on this perimeter)
  //   for(int i = -1; i <=1; i++)
  //   {
  //     for (int j = -1; j<=1;j++)
  //     {
  //       // check the validity of the pointand if were not testing this specific point
  //       if(i<flowpy.get_NRows() && j<flowpy.get_NCols() && i>=0 && j >= 0 && (i !=0 && j!=0))
  //       {
  //         if(is_done[flowpy.get_NodeIndex_from_row_col(row+i,col+j)] != 1)
  //         {
  //           n_adj++;
  //           row_nodes_to_test.push_back(row+i);
  //           col_nodes_to_test.push_back(col+j);
  //         }
  //       }
  //     }
  //   }

  //   // Now we have the number of neighboors

  //   if(n_adj == 1) // easy case
  //   {
  //     node = flowpy.get_NodeIndex_from_row_col(row_nodes_to_test[0],col_nodes_to_test[0]);
  //     is_done[node] = 1;
  //     Perimeter_nodes_sorted.push_back(node);
  //     // distance cauculqtion
  //     flowpy.get_x_and_y_from_current_node(Perimeter_nodes_sorted[Perimeter_nodes_sorted.size()-2],x1,y1);
  //     flowpy.get_x_and_y_from_current_node(node,x2,y2);
  //     map_of_dist_perim[node] = sqrt(pow((x2-x1),2) + pow((y2-y1),2));
  //   }
  //   else if(n_adj>1)
  //   {
  //     vector<int> minindex = Get_Index_Minimum(row_nodes_to_test);
  //     if(minindex.size() == 1)
  //     {
  //       node = flowpy.get_NodeIndex_from_row_col(row_nodes_to_test[minindex[0]],col_nodes_to_test[minindex[0]]);
  //       is_done[node] = 1;
  //       Perimeter_nodes_sorted.push_back(node);
  //       // distance cauculqtion
  //       flowpy.get_x_and_y_from_current_node(Perimeter_nodes_sorted[Perimeter_nodes_sorted.size()-2],x1,y1);
  //       flowpy.get_x_and_y_from_current_node(node,x2,y2);
  //       map_of_dist_perim[node] = sqrt(pow((x2-x1),2) + pow((y2-y1),2));
  //     }

  //   }


  //   // reinitialise everything
  //   flowpy.retrieve_current_row_and_col(node,row,col);
  //   row_nodes_to_test.clear();
  //   col_nodes_to_test.clear();
  //   n_adj = 0;
  // }

}

void LSDBasin::clean_perimeter(LSDFlowInfo& flowpy)
{

  vector<int> light_perimeter; // output perimeter that should be thinned
  // nodes_of_basins; map containing the nodes of the basin (global variable, but not initialized yet)
  int row = 0, col = 0, cptndd = 0, cptndd_tot = 0, this_row = 0, this_col = 0, this_node = 0;
  // tester help to loop through neightboors
  vector<int> tester;
  tester.push_back(-1);
  tester.push_back(1);

  // filling the map of basin nodes to check the nobasin around each node
  for(vector<int>::iterator yo = BasinNodes.begin(); yo!= BasinNodes.end(); yo++)
  {
    nodes_of_basins[*yo] = 1;
  }


  // looping through the the perimeter nodes
  for(vector<int>::iterator uh = Perimeter_nodes.begin();uh != Perimeter_nodes.end(); uh++)
  {
    // the first node is the outlet (I think it is all the time?) so I want it
    if (*uh == BasinNodes[0])
    {
      light_perimeter.push_back(*uh);
    }
    else
    {
      flowpy.retrieve_current_row_and_col(*uh,row,col);
      // looping through the direct row neightboors and incrementing a counter if a direct neighboor is outside the basin
      for(int i = 0; i < 2; i++)
      {
        this_row = row+tester[i]; // tester[0] = -1 and tester 1 = 1
        this_col = col;
        if(this_row>=0 && this_row<get_NRows())
        {
          this_node = flowpy.retrieve_node_from_row_and_column(this_row,this_col);
          if(nodes_of_basins[this_node] !=1) // -> means not in the basin
          {
            cptndd++; // direct neighboors
            cptndd_tot++; // all neighboors (include diagonal)
          }
        }
      }

      // looping through the direct col neightboors and incrementing a counter if a direct neighboor is outside the basin

      for(int i = 0; i < 2; i++)
      {
        this_row = row;
        this_col = col+tester[i];
        if(this_col>=0 && this_col<get_NCols())
        {
          this_node = flowpy.retrieve_node_from_row_and_column(this_row,this_col);
          if(nodes_of_basins[this_node] !=1)
          {
            cptndd++;
            cptndd_tot++;
          }
        }

      }

      for(int i = 0; i<2;i++)
      {
        for(int j =0; j<2 ; j++)
        {
          this_row = row+tester[i];
          this_col = col+tester[j];
          if(this_col>=0 && this_col<get_NCols())
          {
            this_node = flowpy.retrieve_node_from_row_and_column(this_row,this_col);
            if(nodes_of_basins[this_node] !=1)
            {
              cptndd_tot++;
            }
          }
        }
      }
      if(cptndd>0 && cptndd_tot < 6)
      {
        int truc = flowpy.retrieve_node_from_row_and_column(row,col);
        light_perimeter.push_back(truc);
        Perimeter_nodes_map[truc] = 1;
      }

      cptndd = 0;
      cptndd_tot = 0;
    }
  }

  cout << "Cleaned perimeter, you used to have " << Perimeter_nodes.size() << " nodes, now you have " << light_perimeter.size() << " nodes." << endl;

  Perimeter_nodes = light_perimeter;


}


 void LSDBasin::write_elevation_csv(string output_name, LSDFlowInfo& flowpy, LSDRaster& filled)
 {


    // open the file
  ofstream perim_out;
  perim_out.open(output_name.c_str());
  perim_out << "x,y,z" << endl;
  perim_out.precision(5);
  for(vector<int>::iterator Santa = BasinNodes.begin();Santa!=BasinNodes.end();Santa++)
  {
    float this_x,this_y,this_z;
    int this_node;
    this_node = *Santa;
    flowpy.get_x_and_y_from_current_node(this_node,this_x,this_y);
    this_z = filled.get_value_of_point(this_x,this_y);
    perim_out << this_x << ","<< this_y << ","<< this_z << endl;
  }
  perim_out.close();

}
//-----------------------------------------------------------------------------//
void LSDBasin::write_channel_network(string csv_name, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork)
{
  int i, j, JI, stream_order;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  ofstream chan_out;
  chan_out.open(csv_name.c_str());
  chan_out << "JI,stream_order,NI,latitude,longitude" << endl;

  //Loop over every pixel and record it's stream length and basin ID in two vectors
  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    stream_order = JunctionNetwork.get_StreamOrder_of_Node(FlowInfo, BasinNodes[q]);
     if (stream_order != NoDataValue)
     {
       FlowInfo.get_lat_and_long_locations(i, j, latitude, longitude, Converter);
       JI = JunctionNetwork.get_Junction_of_Node(BasinNodes[q], FlowInfo);

       chan_out << JI << ","
                << stream_order << ","
                << BasinNodes[q] << ",";
       chan_out.precision(9);
       chan_out << latitude << ","
                << longitude << ",";
      }
  }
  chan_out.close();
}





// //-------------------------------------------------------------------------//
// // perimeter cleaning Fiona
// // 30/01/18
// //-------------------------------------------------------------------------//
// void LSDBasin::clean_perimeter_fiona(LSDFlowInfo& FlowInfo)
// {
//   // this is another cleaning function. It looks for perimeter nodes which don't border any internal basin nodes and removes them.
//
//   // vector for new perimeter
//   vector<int> thinned_perimeter;
//
//   // loop through each perimeter node and check for internal neighbours
//   for (int i = 0; i < int(Perimeter_nodes.size()); i++)
//   {
//     int this_node = Perimeter_nodes[i];
//
//     // we need to keep the outlet node which should be the first one.
//     if (this_node == BasinNodes[0])
//     {
//       thinned_perimeter.push_back(this_node);
//     }
//     else // not the outlet node
//     {
//
//     }
//
//
//   }
//
//
// }

//----------------------------------------------------------------------------//
// Take points from a LSDSpatialCSVReader and calculate the mean of the points
// within the basin
// FJC 28/09/18
//----------------------------------------------------------------------------//
float LSDBasin::get_basin_mean_from_csv(LSDFlowInfo& FlowInfo, LSDSpatialCSVReader& CSV, string column_name)
{
  // get the latitude, longitude, and data column from the csv
  vector<int> nodes = CSV.get_nodeindices_from_lat_long(FlowInfo);
  vector<float> data = CSV.data_column_to_float(column_name);

  float TotalData = 0;
  int Count = 0;

  // loop through the points and find which ones are in the basin
  for (int i =0; i < int(data.size()); i++)
  {
    if (data[i] != NoDataValue)
    {
      // check if this node is in the basin
      int test = is_node_in_basin(nodes[i]);
      if (test == 1)
      {
        TotalData += data[i];
        Count++;
      }
    }
  }
  float MeanData = TotalData/Count;
  return MeanData;
}

















// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Added lines to separate with LSDCosmoBasin =-=-=-=-=-=-=-=-=-
// I got lost each time I try to find it =-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  +++++++++++++++++
//  +++++++++++++++++
//  +++++++++++++++++
//  +++++++++++++++++
//  +++++++++++++++++
//  COSMOGENIC BASIN
//  +++++++++++++++++
//  +++++++++++++++++
//  +++++++++++++++++
//  +++++++++++++++++
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates a cosmogenic basin. Similar to a normal basin but just has the
// 10Be and 26Al concentrations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::create(int JunctionNumber, LSDFlowInfo& FlowInfo,
                           LSDJunctionNetwork& ChanNet,
                           double N10Be, double delN10Be,
                           double N26Al, double delN26Al)
{

  // The measured 10Be concentration
  measured_N_10Be = N10Be;

  // The measured 26Al concentration
  measured_N_26Al = N26Al;

  // The measured uncertainty in the 10Be concentration
  delN_10Be = delN10Be;

  // The measured uncertainty in the 10Be concentration
  delN_26Al = delN26Al;

  //NO BOUNDS CHECKING ON JunctionNumber

  //setting all of the instance variables for the given junction

  NRows = ChanNet.get_NRows();
  NCols = ChanNet.get_NCols();
  XMinimum = ChanNet.get_XMinimum();
  YMinimum = ChanNet.get_YMinimum();
  DataResolution = ChanNet.get_DataResolution();
  NoDataValue = ChanNet.get_NoDataValue();
  GeoReferencingStrings = ChanNet.get_GeoReferencingStrings();

  Junction = JunctionNumber;

  //cout << "LSDCosmoBasin, line 888, Junction is: " << Junction << endl;

  int basin_outlet = ChanNet.get_Node_of_Junction(Junction);
  BasinNodes = FlowInfo.get_upslope_nodes(basin_outlet);

  //cout << "LSDCosmoBasin, Line 893, basin outlet is: " << basin_outlet << endl;

  NumberOfCells = int(BasinNodes.size());
  Area = NumberOfCells * (DataResolution*DataResolution);

  Beheaded = ChanNet.node_tester(FlowInfo, Junction);

  FlowInfo.retrieve_current_row_and_col(ChanNet.get_Node_of_Junction(Junction), Outlet_i, Outlet_j);

  vector<int> StreamOrderVector = ChanNet.get_StreamOrderVector();

  BasinOrder = StreamOrderVector[Junction];


  int i_max = 0;
  int i_min = 9999999; //a very large number
  int j_max = 0;
  int j_min = 9999999; //a very large number

  int i = 0;
  int j = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);

    if (i > i_max){i_max = i;}
    else if (i < i_min){i_min = i;}
    if (j > j_max){j_max = j;}
    else if (j < j_min){j_min = j;}

  }

  Centroid_i = i_min + ((i_max - i_min)/2);
  Centroid_j = j_min + ((j_max - j_min)/2);   //how do these handle 0.5s ??


  //finished setting all the instance variables


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  This function populates the topographic and production shielding
// It sets the snow shielding to a default of 1
// You need a sheilding raster: we need to implement without one in the near future.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::populate_scaling_vectors(LSDFlowInfo& FlowInfo,
                                               LSDRaster& Elevation_Data,
                                               LSDRaster& T_Shield,
                                               string path_to_atmospheric_data)
{
  int row,col;

  // variables for converting location and elevation
  double this_elevation, this_pressure, this_tshield;

  vector<double> tshield_temp;
  vector<double> prod_temp;
  vector<double> snow_temp;

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();

  // a function for scaling stone production, defaults to 1
  double Fsp = 1.0;

  // the latitude and longitude
  double lat,longitude;

  // decalre converter object
  LSDCoordinateConverterLLandUTM Converter;

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    //exclude NDV from average
    if (Elevation_Data.get_data_element(row,col) != NoDataValue)
    {
      // To get pressure, first get the lat and long
      Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
      //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);

      // now the elevation
      this_elevation = Elevation_Data.get_data_element(row,col);

      // now the pressure
      this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude),
                                        double(this_elevation));

      //cout << "r: " << row << " c: " << col << " lat: " << lat << " long: " << longitude
      //     << " elevation: " << this_elevation << " pressure: " << this_pressure << endl;

      // now get the scaling
      prod_temp.push_back(LSDCRNP.stone2000sp(lat,this_pressure, Fsp));

      // Now get topographic shielding
      this_tshield = double(T_Shield.get_data_element(row,col));
      tshield_temp.push_back(this_tshield);

      // now set the snow sheilding to 1
      snow_temp.push_back(1.0);

    }
    else
    {
      prod_temp.push_back(double(NoDataValue));
      tshield_temp.push_back(double(NoDataValue));
      snow_temp.push_back(double(NoDataValue));
    }
  }

  // set the shielding vectors
  topographic_shielding = tshield_temp;
  production_scaling =  prod_temp;
  snow_shielding = snow_temp;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  This function populates the topographic and production shielding
//  It loads a snow shielding raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::populate_scaling_vectors(LSDFlowInfo& FlowInfo,
                                               LSDRaster& Elevation_Data,
                                               LSDRaster& T_Shield,
                                               LSDRaster& S_Shield,
                                               string path_to_atmospheric_data)
{
  int row,col;

  // variables for converting location and elevation
  double this_elevation, this_pressure, this_tshield, this_sshield;

  vector<double> tshield_temp;
  vector<double> prod_temp;
  vector<double> snow_temp;

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();

  // a function for scaling stone production, defaults to 1
  double Fsp = 1.0;

  // the latitude and longitude
  double lat,longitude;

  // decalre converter object
  LSDCoordinateConverterLLandUTM Converter;

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    //exclude NDV from average
    if (Elevation_Data.get_data_element(row,col) != NoDataValue)
    {
      // To get pressure, first get the lat and long
      Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
      //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);

      // now the elevation
      this_elevation = Elevation_Data.get_data_element(row,col);

      // now the pressure
      this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude),
                                        double(this_elevation));

      //cout << "r: " << row << " c: " << col << " lat: " << lat << " long: " << longitude
      //     << " elevation: " << this_elevation << " pressure: " << this_pressure << endl;


      // now get the scaling
      prod_temp.push_back(LSDCRNP.stone2000sp(lat,this_pressure, Fsp));

      // Now get topographic shielding
      this_tshield = double(T_Shield.get_data_element(row,col));
      tshield_temp.push_back(this_tshield);

      // now get the snow shielding
      this_sshield = double(S_Shield.get_data_element(row,col));
      snow_temp.push_back(this_sshield);

    }
    else
    {
      prod_temp.push_back(double(NoDataValue));
      tshield_temp.push_back(double(NoDataValue));
      snow_temp.push_back(double(NoDataValue));
    }
  }

  // set the shielding vectors
  topographic_shielding = tshield_temp;
  production_scaling =  prod_temp;
  snow_shielding = snow_temp;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function poplates the atmospheric pressure vector.
// It is used for bug-checking and comparison with other cosmo calculators
// Other function do not require this to be called
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::get_atmospheric_pressure(LSDFlowInfo& FlowInfo,
                    LSDRaster& Elevation_Data, string path_to_atmospheric_data)
{
  int row,col;

  // variables for converting location and elevation
  double this_elevation, this_pressure;

  vector<double> pressure_temp;


  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();

  // the latitude and longitude
  double lat,longitude;

  // decalre converter object
  LSDCoordinateConverterLLandUTM Converter;

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {

    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    //exclude NDV from average
    if (Elevation_Data.get_data_element(row,col) != NoDataValue)
    {
      // To get pressure, first get the lat and long
      Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
      //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);

      // now the elevation
      this_elevation = Elevation_Data.get_data_element(row,col);

      // now the pressure
      this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude),
                                        double(this_elevation));

      //cout << "r: " << row << " c: " << col << " lat: " << lat << " long: " << longitude
      //     << " elevation: " << this_elevation << " pressure: " << this_pressure << endl;


      // now get the scaling
      pressure_temp.push_back(this_pressure);

    }
    else
    {
      pressure_temp.push_back(double(NoDataValue));
    }
  }

  // set the pressure vector
  atmospheric_pressure = pressure_temp;

  // print pressure at outlet
  cout << "the pressure at the outlet is: "  << atmospheric_pressure[0] << " mbar" << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function creates the snow and shelf sheilding vectors based on two
// rasters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::populate_snow_and_self_eff_depth_vectors(LSDFlowInfo& FlowInfo,
                                              LSDRaster& snow_eff_depth,
                                              LSDRaster& self_eff_depth)
{
  int row,col;

  // the effective depths at individual nodes
  double this_eff_snow_depth;
  double this_eff_self_depth;

  // temporary vectors that will be copied into the
  vector<double> snow_temp;
  vector<double> self_temp;

  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {
    // get the row and column of the node
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    //exclude NDV from average
    if (snow_eff_depth.get_data_element(row,col) != NoDataValue)
    {
      // get the snow and self shielding
      this_eff_snow_depth = double(snow_eff_depth.get_data_element(row,col));
      this_eff_self_depth = double(self_eff_depth.get_data_element(row,col));

      // add data to the vectors
      snow_temp.push_back(this_eff_snow_depth);
      self_temp.push_back(this_eff_self_depth);
    }
  }

  // update the vectors in the basin object
  self_shield_eff_depth = self_temp;
  snow_shield_eff_depth = snow_temp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function creates the snow and shelf sheilding vectors based on a
// double and a float
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::populate_snow_and_self_eff_depth_vectors(LSDFlowInfo& FlowInfo,
                                              double snow_eff_depth,
                                              LSDRaster& self_eff_depth)
{
  int row,col;

  // the effective depths at individual nodes
  double this_eff_self_depth;

  // temporary vectors that will be copied into the
  vector<double> snow_temp;
  vector<double> self_temp;

  // first put the one element in the snow temp vector
  snow_temp.push_back(snow_eff_depth);

  // now loop through the other vector adding elements
  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {
    // get the row and column of the node
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    //exclude NDV from average
    if (self_eff_depth.get_data_element(row,col) != NoDataValue)
    {
      // get the snow and self shielding
      this_eff_self_depth = double(self_eff_depth.get_data_element(row,col));

      // add data to the vectors
      self_temp.push_back(this_eff_self_depth);
    }
  }

  // update the vectors in the basin object
  self_shield_eff_depth = self_temp;
  snow_shield_eff_depth = snow_temp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function creates the snow and shelf sheilding vectors based on a
// double and a float
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::populate_snow_and_self_eff_depth_vectors(LSDFlowInfo& FlowInfo,
                                              LSDRaster& snow_eff_depth,
                                              double self_eff_depth)
{

  cout << "I am getting effective depths, I have a snow raster but not a self raster" << endl;
  int row,col;

  // the effective depths at individual nodes
  double this_eff_snow_depth;

  // temporary vectors that will be copied into the
  vector<double> snow_temp;
  vector<double> self_temp;

  // first put the one element in the self temp vector
  self_temp.push_back(self_eff_depth);

  // now loop through the other vector adding elements
  for (int q = 0; q < int(BasinNodes.size()); ++q)
  {
    // get the row and column of the node
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    //exclude NDV from average
    if (snow_eff_depth.get_data_element(row,col) != NoDataValue)
    {
      // get the row and column of the node
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      // get the snow and self shielding
      this_eff_snow_depth = double(snow_eff_depth.get_data_element(row,col));

      // add data to the vectors
      snow_temp.push_back(this_eff_snow_depth);
    }
  }

  // update the vectors in the basin object
  self_shield_eff_depth = self_temp;
  snow_shield_eff_depth = snow_temp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function creates the snow and shelf shielding vectors based on a
// double and a float
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::populate_snow_and_self_eff_depth_vectors(double snow_eff_depth,
                                                             double self_eff_depth)
{
  // temporary vectors that will be copied into the
  vector<double> snow_temp;
  vector<double> self_temp;

  // first put the one element in the snow temp vector
  self_temp.push_back(self_eff_depth);
  snow_temp.push_back(snow_eff_depth);

  // update the vectors in the basin object
  self_shield_eff_depth = self_temp;
  snow_shield_eff_depth = snow_temp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This resets the snow and self shielding vectors
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::reset_snow_and_self_eff_depth_vectors()
{
  // temporary vectors that will be copied into the
  vector<double> snow_temp;
  vector<double> self_temp;

  // update the vectors in the basin object
  self_shield_eff_depth = self_temp;
  snow_shield_eff_depth = snow_temp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the erosion rate calculations with formal error analysis
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCosmoBasin::full_CRN_erosion_analysis(double Nuclide_conc, string Nuclide,
                              double Nuclide_conc_err, double prod_uncert_factor,
                              string Muon_scaling)
{
  // the vector for holding the erosion rates and uncertainties
  vector<double> erate_uncert_vec;

  double erate;                   // effective erosion rate g/cm^2/yr
  double erate_external_plus;     // effective erosion rate g/cm^2/yr for AMS uncertainty +
  double erate_external_minus;    // effective erosion rate g/cm^2/yr for AMS uncertainty -
  double dEdExternal;             // change in erosion rate for change in AMS atoms/g
  double External_uncert;         // uncertainty of effective erosion rate g/cm^2/yr for AMS


  double production_uncertainty;  // a lumped production uncertainty value.
                                  // not generally used but needs to be passed
                                  // to the erosion finding routines as a parameter

  // variables for the muon uncertainty
  double average_production_rate; // The average production rate, used in uncertainty
                                  // calculations
  double erate_muon_scheme_schaller;  // erosion rate using schaller scheme
  double erate_muon_scheme_braucher;  // erosion rate using braucher scheme

  double dEdMuonScheme;           // change in erosion rate for change in Muon Scheme
  double Muon_uncert;             // uncertainty of effective erosion rate
                                  // in g/cm^2/yr for different muon schemes

  // variable for the production uncertainty
  double erate_prod_plus;   // erosion rate for positive production uncertainty
  double erate_prod_minus;  // erosion rate for negative production uncertainty

  double dEdProduction;        // change in erosion rate for change in production
  double Prod_uncert;          // uncertainty of effective erosion rate
                               // in g/cm^2/yr for production uncertainty

  double this_prod_difference; // the difference in production for production uncertainty

  // initially we do not modify production rates
  bool is_production_uncertainty_plus_on = false;
  bool is_production_uncertainty_minus_on = false;

  // first get the prediction of the erosion rate
  erate = predict_CRN_erosion(Nuclide_conc, Nuclide, prod_uncert_factor,
                              Muon_scaling, production_uncertainty,
                              average_production_rate,
                              is_production_uncertainty_plus_on,
                              is_production_uncertainty_minus_on);

  double no_prod_uncert = 1.0;    // set the scheme to no production uncertainty
                                  // for the external uncertainty
  // now get the external uncertainty
  erate_external_plus = predict_CRN_erosion(Nuclide_conc+Nuclide_conc_err, Nuclide,
                                             no_prod_uncert, Muon_scaling,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on);
  erate_external_minus = predict_CRN_erosion(Nuclide_conc-Nuclide_conc_err, Nuclide,
                                             no_prod_uncert, Muon_scaling,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on);
  dEdExternal = (erate_external_plus-erate_external_minus)/(2*Nuclide_conc_err);
  External_uncert = fabs(dEdExternal*Nuclide_conc_err);

  //cout << "LSDCosmoBasin, line 1160, erate: " << erate << " and uncertainty: "
  //     << External_uncert << endl;

  // now calculate uncertainty from different muon scaling schemes.
  // The end members are Braucher and Schaller
  string braucher_string = "Braucher";
  string schaller_string = "Schaller";

  // get the difference in the pair
  LSDCRNParameters LSDCRNP;
  int pair_key = 0;       // this is for braucher-schaller
  vector<double> muon_uncert_diff = LSDCRNP.get_uncertainty_scaling_pair(pair_key);

  double this_muon_uncert_dif;
  if(Nuclide == "Be10")
  {
    this_muon_uncert_dif = muon_uncert_diff[0];
  }
  else if (Nuclide == "Al26")
  {
    this_muon_uncert_dif = muon_uncert_diff[1];
  }
  else
  {
    cout << "LINE 1295 LSDBasin you did not supply a valid nuclide, defaulting to 10Be" << endl;
    Nuclide = "Be10";
    this_muon_uncert_dif = muon_uncert_diff[0];
  }

  // now get the muon uncertainty
  erate_muon_scheme_schaller = predict_CRN_erosion(Nuclide_conc, Nuclide,
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on);
  erate_muon_scheme_braucher = predict_CRN_erosion(Nuclide_conc, Nuclide,
                                             no_prod_uncert, braucher_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on);
  dEdMuonScheme = (erate_muon_scheme_schaller-erate_muon_scheme_braucher)/
                  this_muon_uncert_dif;
  Muon_uncert = fabs(dEdMuonScheme*this_muon_uncert_dif);

  //cout << "LSDCosmoBasin, Line 1292, change in scaling production rate: "
  //     << this_muon_uncert_dif << " erate Schal: "
  //     << erate_muon_scheme_schaller << " erate Braucher: "
  //     << erate_muon_scheme_braucher << " and erate uncert: " << Muon_uncert << endl;

  // now get the production uncertainty
  // first set the scaling
  // reset scaling parameters. This is necessary since the F values are
  // reset for local scaling
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

  // now get the uncertainty parameters
  vector<double> prod;
  double prod_plus,prod_minus;
  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[0];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[0];
  }
  else if (Nuclide == "Al26")
  {
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[1];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[1];
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[0];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[0];
  }
  //cout << "Prod plus: " << prod_plus << " prod minus: " << prod_minus << endl;
  this_prod_difference = prod_plus+prod_minus;


  is_production_uncertainty_plus_on = true;
  is_production_uncertainty_minus_on = false;
  erate_prod_plus = predict_CRN_erosion(Nuclide_conc, Nuclide,
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on);

  is_production_uncertainty_plus_on = false;
  is_production_uncertainty_minus_on = true;
  erate_prod_minus = predict_CRN_erosion(Nuclide_conc, Nuclide,
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on);

  dEdProduction = (erate_prod_plus-erate_prod_minus)/
                   this_prod_difference;
  Prod_uncert = fabs(dEdProduction*this_prod_difference);

  //cout << "LSDCosmoBasin, Line 1368, change in production rate for production uncertainty: "
  //     << this_prod_difference << " erate plus: "
  //     << erate_prod_plus << " erate minus: "
  //     << erate_prod_minus << " and erate uncert: " << Prod_uncert << endl;



  // now calculate the total uncertainty
  double total_uncert = sqrt( External_uncert*External_uncert +
                              Muon_uncert*Muon_uncert +
                              Prod_uncert*Prod_uncert);

  erate_uncert_vec.push_back(erate);
  erate_uncert_vec.push_back(External_uncert);
  erate_uncert_vec.push_back(Muon_uncert);
  erate_uncert_vec.push_back(Prod_uncert);
  erate_uncert_vec.push_back(total_uncert);

  return erate_uncert_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function wraps the erosion rate calculations with formal error analysis
// It is the version that includes basin nesting
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCosmoBasin::full_CRN_erosion_analysis_nested(LSDRaster& known_eff_erosion,
                              LSDFlowInfo& FlowInfo, double Nuclide_conc, string Nuclide,
                              double Nuclide_conc_err, double prod_uncert_factor,
                              string Muon_scaling)
{
  // the vector for holding the erosion rates and uncertainties
  vector<double> erate_uncert_vec;

  double erate;                   // effective erosion rate g/cm^2/yr
  double erate_external_plus;     // effective erosion rate g/cm^2/yr for AMS uncertainty +
  double erate_external_minus;    // effective erosion rate g/cm^2/yr for AMS uncertainty -
  double dEdExternal;             // change in erosion rate for change in AMS atoms/g
  double External_uncert;         // uncertainty of effective erosion rate g/cm^2/yr for AMS


  double production_uncertainty;  // a lumped production uncertainty value.
                                  // not generally used but needs to be passed
                                  // to the erosion finding routines as a parameter

  // variables for the muon uncertainty
  double average_production_rate; // The average production rate, used in uncertainty
                                  // calculations
  double erate_muon_scheme_schaller;  // erosion rate using schaller scheme
  double erate_muon_scheme_braucher;  // erosion rate using braucher scheme

  double dEdMuonScheme;           // change in erosion rate for change in Muon Scheme
  double Muon_uncert;             // uncertainty of effective erosion rate
                                  // in g/cm^2/yr for different muon schemes

  // variable for the production uncertainty
  double erate_prod_plus;   // erosion rate for positive production uncertainty
  double erate_prod_minus;  // erosion rate for negative production uncertainty

  double dEdProduction;        // change in erosion rate for change in production
  double Prod_uncert;          // uncertainty of effective erosion rate
                               // in g/cm^2/yr for production uncertainty

  double this_prod_difference; // the difference in production for production uncertainty

  // initially we do not modify production rates
  bool is_production_uncertainty_plus_on = false;
  bool is_production_uncertainty_minus_on = false;

  // first get the prediction of the erosion rate
  erate = predict_CRN_erosion_nested(Nuclide_conc, Nuclide, prod_uncert_factor,
                              Muon_scaling, production_uncertainty,
                              average_production_rate,
                              is_production_uncertainty_plus_on,
                              is_production_uncertainty_minus_on,
                              known_eff_erosion, FlowInfo);

  cout << "Hey Bubba, I got the erosion rate!!!: " << erate << endl << endl << endl;

  double no_prod_uncert = 1.0;    // set the scheme to no production uncertainty
                                  // for the external uncertainty
  // now get the external uncertainty
  erate_external_plus = predict_CRN_erosion_nested(Nuclide_conc+Nuclide_conc_err, Nuclide,
                                             no_prod_uncert, Muon_scaling,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on,
                                             known_eff_erosion, FlowInfo);
  erate_external_minus = predict_CRN_erosion_nested(Nuclide_conc-Nuclide_conc_err, Nuclide,
                                             no_prod_uncert, Muon_scaling,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on,
                                             known_eff_erosion, FlowInfo);
  dEdExternal = (erate_external_plus-erate_external_minus)/(2*Nuclide_conc_err);
  External_uncert = fabs(dEdExternal*Nuclide_conc_err);

  //cout << "LSDCosmoBasin, line 1160, erate: " << erate << " and uncertainty: "
  //     << External_uncert << endl;

  // now calculate uncertainty from different muon scaling schemes.
  // The end members are Braucher and Schaller
  string braucher_string = "Braucher";
  string schaller_string = "Schaller";

  // get the difference in the pair
  LSDCRNParameters LSDCRNP;
  int pair_key = 0;       // this is for braucher-schaller
  vector<double> muon_uncert_diff = LSDCRNP.get_uncertainty_scaling_pair(pair_key);

  double this_muon_uncert_dif;
  if(Nuclide == "Be10")
  {
    this_muon_uncert_dif = muon_uncert_diff[0];
  }
  else if (Nuclide == "Al26")
  {
    this_muon_uncert_dif = muon_uncert_diff[1];
  }
  else
  {
    cout << "LINE 1295 LSDBasin you did not supply a valid nuclide, defaulting to 10Be" << endl;
    Nuclide = "Be10";
    this_muon_uncert_dif = muon_uncert_diff[0];
  }

  // now get the muon uncertainty
  erate_muon_scheme_schaller = predict_CRN_erosion_nested(Nuclide_conc, Nuclide,
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on,
                                             known_eff_erosion, FlowInfo);
  erate_muon_scheme_braucher = predict_CRN_erosion_nested(Nuclide_conc, Nuclide,
                                             no_prod_uncert, braucher_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on,
                                             known_eff_erosion, FlowInfo);
  dEdMuonScheme = (erate_muon_scheme_schaller-erate_muon_scheme_braucher)/
                  this_muon_uncert_dif;
  Muon_uncert = fabs(dEdMuonScheme*this_muon_uncert_dif);

  //cout << "LSDCosmoBasin, Line 1292, change in scaling production rate: "
  //     << this_muon_uncert_dif << " erate Schal: "
  //     << erate_muon_scheme_schaller << " erate Braucher: "
  //     << erate_muon_scheme_braucher << " and erate uncert: " << Muon_uncert << endl;

  // now get the production uncertainty
  // first set the scaling
  // reset scaling parameters. This is necessary since the F values are
  // reset for local scaling
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

  // now get the uncertainty parameters
  vector<double> prod;
  double prod_plus,prod_minus;
  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[0];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[0];
  }
  else if (Nuclide == "Al26")
  {
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[1];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[1];
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_plus();
    prod_plus = prod[0];
    prod = LSDCRNP.set_P0_CRONUS_uncertainty_minus();
    prod_minus = prod[0];
  }
  //cout << "Prod plus: " << prod_plus << " prod minus: " << prod_minus << endl;
  this_prod_difference = prod_plus+prod_minus;


  is_production_uncertainty_plus_on = true;
  is_production_uncertainty_minus_on = false;
  erate_prod_plus = predict_CRN_erosion_nested(Nuclide_conc, Nuclide,
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on,
                                             known_eff_erosion, FlowInfo);

  is_production_uncertainty_plus_on = false;
  is_production_uncertainty_minus_on = true;
  erate_prod_minus = predict_CRN_erosion_nested(Nuclide_conc, Nuclide,
                                             no_prod_uncert, schaller_string,
                                             production_uncertainty,
                                             average_production_rate,
                                             is_production_uncertainty_plus_on,
                                             is_production_uncertainty_minus_on,
                                             known_eff_erosion, FlowInfo);

  dEdProduction = (erate_prod_plus-erate_prod_minus)/
                   this_prod_difference;
  Prod_uncert = fabs(dEdProduction*this_prod_difference);

  //cout << "LSDCosmoBasin, Line 1368, change in production rate for production uncertainty: "
  //     << this_prod_difference << " erate plus: "
  //     << erate_prod_plus << " erate minus: "
  //     << erate_prod_minus << " and erate uncert: " << Prod_uncert << endl;



  // now calculate the total uncertainty
  double total_uncert = sqrt( External_uncert*External_uncert +
                              Muon_uncert*Muon_uncert +
                              Prod_uncert*Prod_uncert);

  erate_uncert_vec.push_back(erate);
  erate_uncert_vec.push_back(External_uncert);
  erate_uncert_vec.push_back(Muon_uncert);
  erate_uncert_vec.push_back(Prod_uncert);
  erate_uncert_vec.push_back(total_uncert);

  cout << endl << endl;
  cout << "erate: " << erate << endl;
  cout << "Ext uncert: " <<  External_uncert << endl;
  cout << "Muon uncert: " <<  Muon_uncert << endl;
  cout << "Prod uncert: " <<  Prod_uncert << endl;
  cout << "total uncert: " <<  total_uncert << endl;



  return erate_uncert_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function is used to get the erosion rate using Newton-Raphson
// method of root finding
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoBasin::predict_CRN_erosion(double Nuclide_conc, string Nuclide,
                                          double prod_uncert_factor,
                                          string Muon_scaling,
                                          double& production_uncertainty,
                                          double& average_production,
                                          bool is_production_uncertainty_plus_on,
                                          bool is_production_uncertainty_minus_on)
{
  // effective erosion rates (in g/cm^2/yr) for running the Newton Raphson
  // iterations
  double erate_guess;
  double eff_erate_guess;
  //double this_eff_erosion_rate;
  //double d_eff_erosion_rate;

  double rho = 2650;  // this is the rock density but in fact it doesn't
                      // really play a role since it is factored into the
                      // apparent erosion to get erosion in mm/yr but the divided
                      // out again.
                      // The value 2650 is used because this is the default in cosmocalc

  // production uncertainty factor is a multiplier that sets the production
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

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

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }


  // now set the scaling parameters
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
    cout << "LSDBasin line 1583, You didn't choos a valid nuclide. Defaulting"
         << " to 10Be." << endl;
    Nuclide = "Be10";
    nuclide_scaling_switches[0] = true;
  }

  // now get the guess from the particle
  // the initial guess just takes scaling from the outlet, and then
  // uses that for the entire basin. This guess will probably be quite
  // far off, but provides a useful starting point

  // the elevation, snow shielding, topographic shielding
  // and production scaling are all independent of the erosion rate
  // and are calculated seperately.
  // IMPORTANT populate scaling vectors must be called in advance!
  if(  production_scaling.size() < 1 )
  {
    cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
         << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
  }

  double total_shielding;
  // if the shielding is based on effective depths
  if (snow_shielding[0] == 0)
  {
    total_shielding =  production_scaling[0]*topographic_shielding[0];
  }
  else
  {
    total_shielding =  production_scaling[0]*topographic_shielding[0]*
                        snow_shielding[0];
  }

  //cout << "LSDBasin line 1128 Prod scaling is: " << production_scaling[0] << endl;

  //cout << "LSDBasin line 1129; total scaling is: " << total_shielding << endl;

  // now recalculate F values to match the total shielding
  LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);
  if (snow_shielding[0] == 0)
  {
    LSDCRNP.set_neutron_scaling(production_scaling[0],topographic_shielding[0],
                             snow_shielding[0]);
  }
  else
  {
    double this_sshield = 1.0;
    LSDCRNP.set_neutron_scaling(production_scaling[0],topographic_shielding[0],
                                this_sshield);
  }

  // at the moment do only the outlet
  bool data_from_outlet_only = false;
  //cout << "LSDBasin line 1739 WARNING YOU ARE ONLY CALCULATING THE OUTLET " << endl;

  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    //cout << "LINE 1134 LSDBasin Nuclide conc is: " << Nuclide_conc << endl;
    eroded_particle.setConc_10Be(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_10Be_neutron_only(rho, LSDCRNP);
    //cout << "Be10, initial erate guess in m/yr with density " << rho << ": " << erate_guess << endl;
  }
  else if (Nuclide == "Al26")
  {
    eroded_particle.setConc_26Al(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_26Al_neutron_only(rho, LSDCRNP);
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    eroded_particle.setConc_10Be(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_10Be_neutron_only(rho, LSDCRNP);
  }

  // convert to  g/cm^2/yr
  eff_erate_guess = 0.1*erate_guess*rho;

  // now using this as the initial guess, use Newton-Raphson to zero in on the
  // correct erosion rate
  double eff_e_new = eff_erate_guess; // the erosion rate upon which we iterate
  double eff_e_change;                // the change in erosion rate between iterations
  double tolerance = 1e-10;           // tolerance for a change in the erosion rate
                                      // between Newton-Raphson iterations
  double eff_e_displace = 1e-6;       // A small displacment in the erosion rate used
                                      // to calculate the derivative
  double N_this_step;                 // the concentration of the nuclide reported this step
  double N_displace;                  // the concentration at the displaced erosion rate
  double N_derivative;                // dN/de derivative for Newton-Raphson
  double f_x;                         // the function being tested by newton raphson
  double f_x_displace;                // the displaced function (for calculating the derivative)

  double this_step_prod_uncert;       // the uncertainty in the production rate
                                      // from this step
  double displace_uncertainty;        // the uncertainty from the displaced calculations
                                      // is not used so a dummy variable is used here

  double this_step_average_production;// the average production rate for this step
  double displace_average_production; // aveage production for the displace step

  do
  {
    // get the new values
    //cout << "Taking a step, eff_e: " << eff_e_new << " data_outlet? " <<  data_from_outlet_only;
    if(self_shield_eff_depth.size() < 1 && snow_shield_eff_depth.size() < 1)
    {
      //cout << "LSDBasin line 1630, You are doing this wihout the effective depth driven shielding" << endl;

      N_this_step = predict_mean_CRN_conc(eff_e_new, Nuclide,prod_uncert_factor,
                                        Muon_scaling,data_from_outlet_only,
                                        this_step_prod_uncert,
                                        this_step_average_production,
                                        is_production_uncertainty_plus_on,
                                        is_production_uncertainty_minus_on);

      // now get the derivative
      N_displace = predict_mean_CRN_conc(eff_e_new+eff_e_displace,Nuclide,
                                       prod_uncert_factor,Muon_scaling,
                                       data_from_outlet_only,displace_uncertainty,
                                       displace_average_production,
                                       is_production_uncertainty_plus_on,
                                       is_production_uncertainty_minus_on);
    }
    else   // if self and snow sheilding are caluclated based on effective depths
    {
      //cout << "LSDBasin line 1649 You are doing this wih the effective depth driven shielding" << endl;

      N_this_step = predict_mean_CRN_conc_with_snow_and_self(eff_e_new, Nuclide,
                                        prod_uncert_factor,
                                        Muon_scaling,data_from_outlet_only,
                                        this_step_prod_uncert,
                                        this_step_average_production,
                                        is_production_uncertainty_plus_on,
                                        is_production_uncertainty_minus_on);
      //cout << " Conc: " << N_this_step << endl;

      // now get the derivative
      N_displace = predict_mean_CRN_conc_with_snow_and_self(eff_e_new+eff_e_displace,Nuclide,
                                       prod_uncert_factor,Muon_scaling,
                                       data_from_outlet_only,displace_uncertainty,
                                       displace_average_production,
                                       is_production_uncertainty_plus_on,
                                       is_production_uncertainty_minus_on);
    }

    f_x =  N_this_step-Nuclide_conc;
    f_x_displace =  N_displace-Nuclide_conc;

    N_derivative = (f_x_displace-f_x)/eff_e_displace;

    if(N_derivative != 0)
    {
      eff_e_new = eff_e_new-f_x/N_derivative;

      // check to see if the difference in erosion rates meet a tolerance
      eff_e_change = f_x/N_derivative;
      //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
    }
    else
    {
      eff_e_change = 0;
    }

  } while(fabs(eff_e_change) > tolerance);

  // replace the production uncertainty
  production_uncertainty = this_step_prod_uncert;

  // replace the average production
  average_production = this_step_average_production;

  return eff_e_new;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function is used to get the erosion rate using Newton-Raphson
// method of root finding
// This allows for nesting because you can fix the erosion rate from
// for certain basin pixels.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoBasin::predict_CRN_erosion_nested(double Nuclide_conc, string Nuclide,
                                          double prod_uncert_factor,
                                          string Muon_scaling,
                                          double& production_uncertainty,
                                          double& average_production,
                                          bool is_production_uncertainty_plus_on,
                                          bool is_production_uncertainty_minus_on,
                                          LSDRaster& eff_erosion_raster,
                                          LSDFlowInfo& FlowInfo)
{
  // effective erosion rates (in g/cm^2/yr) for running the Newton Raphson
  // iterations
  double erate_guess;
  double eff_erate_guess;
  //double this_eff_erosion_rate;
  //double d_eff_erosion_rate;

  double rho = 2650;  // this is the rock density but in fact it doesn't
                      // really play a role since it is factored into the
                      // apparent erosion to get erosion in mm/yr but the divided
                      // out again.
                      // The value 2650 is used because this is the default in cosmocalc

  // production uncertainty factor is a multiplier that sets the production
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

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

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }


  // now set the scaling parameters
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

  // now get the guess from the particle
  // the initial guess just takes scaling from the outlet, and then
  // uses that for the entire basin. This guess will probably be quite
  // far off, but provides a useful starting point

  // the elevation, snow shielding, topographic shielding
  // and production scaling are all independent of the erosion rate
  // and are calculated seperately.
  // IMPORTANT populate scaling vectors must be called in advance!
  if(  production_scaling.size() < 1 )
  {
    cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
         << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
  }

  double total_shielding;
  // if the shielding is based on effective depths
  if (snow_shielding[0] == 0)
  {
    total_shielding =  production_scaling[0]*topographic_shielding[0];
  }
  else
  {
    total_shielding =  production_scaling[0]*topographic_shielding[0]*
                        snow_shielding[0];
  }

  //cout << "LSDBasin line 1128 Prod scaling is: " << production_scaling[0] << endl;

  //cout << "LSDBasin line 1129; total scaling is: " << total_shielding << endl;

  // now recalculate F values to match the total shielding
  LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);
  if (snow_shielding[0] == 0)
  {
    LSDCRNP.set_neutron_scaling(production_scaling[0],topographic_shielding[0],
                             snow_shielding[0]);
  }
  else
  {
    double this_sshield = 1.0;
    LSDCRNP.set_neutron_scaling(production_scaling[0],topographic_shielding[0],
                                this_sshield);
  }

  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {
    //cout << "LINE 1134 LSDBasin Nuclide conc is: " << Nuclide_conc << endl;
    eroded_particle.setConc_10Be(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_10Be_neutron_only(rho, LSDCRNP);
    cout << "Be10, initial erate guess in m/yr with density " << rho << ": " << erate_guess << endl;
  }
  else if (Nuclide == "Al26")
  {
    eroded_particle.setConc_26Al(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_26Al_neutron_only(rho, LSDCRNP);
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    eroded_particle.setConc_10Be(Nuclide_conc);
    erate_guess = eroded_particle.apparent_erosion_10Be_neutron_only(rho, LSDCRNP);
  }

  // convert to  g/cm^2/yr
  eff_erate_guess = 0.1*erate_guess*rho;

  // now using this as the initial guess, use Newton-Raphson to zero in on the
  // correct erosion rate
  double eff_e_new = eff_erate_guess; // the erosion rate upon which we iterate
  double eff_e_change;                // the change in erosion rate between iterations
  double tolerance = 1e-10;           // tolerance for a change in the erosion rate
                                      // between Newton-Raphson iterations
  double eff_e_displace = 1e-6;       // A small displacment in the erosion rate used
                                      // to calculate the derivative
  double N_this_step;                 // the concentration of the nuclide reported this step
  double N_displace;                  // the concentration at the displaced erosion rate
  double N_derivative;                // dN/de derivative for Newton-Raphson
  double f_x;                         // the function being tested by newton raphson
  double f_x_displace;                // the displaced function (for calculating the derivative)

  double this_step_prod_uncert;       // the uncertainty in the production rate
                                      // from this step
  double displace_uncertainty;        // the uncertainty from the displaced calculations
                                      // is not used so a dummy variable is used here

  double this_step_average_production;// the average production rate for this step
  double displace_average_production; // aveage production for the displace step

  // now check if there are unknown erosion rates in basin
  bool there_are_unknowns = are_there_unknown_erosion_rates_in_basin(eff_erosion_raster,FlowInfo);
  if ( there_are_unknowns == false)
  {
    cout << "There are no unknown erosion rates in this basin." << endl;
    cout << " Taking the average of the known erosion rates." << endl;
    eff_e_new = CalculateBasinMean(FlowInfo, eff_erosion_raster);
  }
  else
  {

    // check to see if there are snow and shielding values, if not populate the vecotrs
    if(self_shield_eff_depth.size() < 1 && snow_shield_eff_depth.size() < 1)
    {
      cout << "You don't seem to have populated the snow and self shielding vectors." << endl;
      cout << "Setting these to 0 shielding" << endl;
      double snow_eff_depth = 0;
      double self_eff_depth = 0;
      populate_snow_and_self_eff_depth_vectors(snow_eff_depth, self_eff_depth);
    }

    // now use newton iteration to get the correct cocentration
    do
    {
      N_this_step = predict_mean_CRN_conc_with_snow_and_self_nested(eff_e_new,
                                        eff_erosion_raster,FlowInfo,
                                        Nuclide,
                                        prod_uncert_factor,
                                        Muon_scaling,
                                        this_step_prod_uncert,
                                        this_step_average_production,
                                        is_production_uncertainty_plus_on,
                                        is_production_uncertainty_minus_on);
      //cout << " Conc: " << N_this_step << endl;

      // now get the derivative
      N_displace = predict_mean_CRN_conc_with_snow_and_self_nested(eff_e_new+eff_e_displace,
                                       eff_erosion_raster,FlowInfo,Nuclide,
                                       prod_uncert_factor,Muon_scaling,
                                       displace_uncertainty,
                                       displace_average_production,
                                       is_production_uncertainty_plus_on,
                                       is_production_uncertainty_minus_on);

      f_x =  N_this_step-Nuclide_conc;
      f_x_displace =  N_displace-Nuclide_conc;

      N_derivative = (f_x_displace-f_x)/eff_e_displace;

      if(N_derivative != 0)
      {
        eff_e_new = eff_e_new-f_x/N_derivative;

        // check to see if the difference in erosion rates meet a tolerance
        eff_e_change = f_x/N_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        eff_e_change = 0;
      }
    } while(fabs(eff_e_change) > tolerance);
  }

  // replace the production uncertainty
  production_uncertainty = this_step_prod_uncert;

  // replace the average production
  average_production = this_step_average_production;

  return eff_e_new;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-











//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the concentration of a nuclide as  function of erosion rate
// The erosion rate should be in g/cm^2/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoBasin::predict_mean_CRN_conc(double eff_erosion_rate, string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            bool data_from_outlet_only,
                                            double& production_uncertainty,
                                            double& average_production,
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on)
{
  // production uncertainty factor is a multiplier that sets the production
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

  // these parameters give the average production rate of the entore basin,
  // along with the magnitude of the production uncertainty
  double cumulative_production_rate = 0;
  double average_production_rate;
  double average_production_uncertainty;

  // the average atoms per gram of the nuclide
  double BasinAverage;

  // the total shielding. A product of snow, topographic and production scaling
  double total_shielding;
  double total_shielding_no_uncert;

  // the total atomic concentration of the nuclude in question
  double Total_N = 0;

  int count_samples = 0;

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

  // set a special case if the outlet flag is true
  int end_node;
  if (data_from_outlet_only)
  {
    end_node = 1;
  }
  else
  {
    end_node =  int(BasinNodes.size());
  }

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
    cout << "LSDBasin line 1583, You didn't choos a valid nuclide. Defaulting"
         << " to 10Be." << endl;
    Nuclide = "Be10";
    nuclide_scaling_switches[0] = true;
  }

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }


  // loop through the elevation data
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {
      count_samples++;

      // reset scaling parameters. This is necessary since the F values are
      // reset for local scaling
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

      // the elevation, snow shielding, topographic shielding
      // and production scaling are all independent of the erosion rate
      // and are calculated seperately.
      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // now you need logic to test if you are accounting for self shielding
      if ( self_shielding.size() < 1 )
      {
        total_shielding_no_uncert = production_scaling[q]*topographic_shielding[q]*
                                  snow_shielding[q];
        total_shielding = prod_uncert_factor*total_shielding_no_uncert;
        cumulative_production_rate += total_shielding_no_uncert;
      }
      else
      {
        total_shielding_no_uncert = production_scaling[q]*topographic_shielding[q]*
                                  snow_shielding[q]*self_shielding[q];
        total_shielding = prod_uncert_factor*total_shielding_no_uncert;
        cumulative_production_rate += total_shielding_no_uncert;
      }

      LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

      // get the nuclide concentration from this node
      if (Nuclide == "Be10")
      {

        eroded_particle.update_10Be_SSfull(eff_erosion_rate,LSDCRNP);
        Total_N+=eroded_particle.getConc_10Be();
      }
      else if (Nuclide == "Al26")
      {
        eroded_particle.update_26Al_SSfull(eff_erosion_rate,LSDCRNP);
        Total_N+=eroded_particle.getConc_26Al();
      }
      else
      {
        cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
        cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
        cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
        eroded_particle.update_10Be_SSfull(eff_erosion_rate,LSDCRNP);
        Total_N+=eroded_particle.getConc_10Be();
      }

      //cout << endl << endl << "LINE 2052, total shield: " << total_shielding
      //     << " erosion: " << eff_erosion_rate << " erosion in cm/kyr with rho = 2650: "
      //     << eff_erosion_rate*1e6/2650.0 << " and N: " << eroded_particle.getConc_10Be() << endl;

    }
  }

  BasinAverage = Total_N/double(count_samples);
  average_production_rate = cumulative_production_rate/double(count_samples);
  average_production_uncertainty = average_production_rate*fabs(1-prod_uncert_factor);

  // replace the production uncertanty
  production_uncertainty = average_production_uncertainty;

  // replace the average production rate
  average_production = average_production_rate;

  return BasinAverage;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the concentration of a nuclide as  function of erosion rate
// The erosion rate should be in g/cm^2/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoBasin::predict_mean_CRN_conc_with_snow_and_self(double eff_erosion_rate,
                                            string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            bool data_from_outlet_only,
                                            double& production_uncertainty,
                                            double& average_production,
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on)
{
  // production uncertainty factor is a multiplier that sets the production
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

  // these parameters give the average production rate of the entore basin,
  // along with the magnitude of the production uncertainty
  double cumulative_production_rate = 0;
  double average_production_rate;
  double average_production_uncertainty;

  // the average atoms per gram of the nuclide
  double BasinAverage;

  // the total shielding. A product of snow, topographic and production scaling
  double total_shielding;
  double total_shielding_no_uncert;

  // the total atomic concentration of the nuclude in question
  double Total_N = 0;

  int count_samples = 0;

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

  // set a special case if the outlet flag is true
  int end_node;
  if (data_from_outlet_only)
  {
    end_node = 1;
  }
  else
  {
    end_node =  int(BasinNodes.size());
  }

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

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }

  // parameters for the shielding
  double this_top_eff_depth;
  double this_bottom_eff_depth;

  // loop through the elevation data
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {
      count_samples++;

      // reset scaling parameters. This is necessary since the F values are
      // reset for local scaling
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

      // the elevation, snow shielding, topographic shielding
      // and production scaling are all independent of the erosion rate
      // and are calculated seperately.
      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // now you need logic to test if you are accounting for self shielding
      total_shielding_no_uncert = production_scaling[q]*topographic_shielding[q];
      total_shielding = prod_uncert_factor*total_shielding_no_uncert;
      cumulative_production_rate += total_shielding_no_uncert;

      // scale the F values
      LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

      // check to see if the shielding data exist and if so get the top and bottom
      // effective depths


      // first get snow shielding (as implemented by and effective depth of snow)
      if (snow_shield_eff_depth.size() < 1)
      {
        this_top_eff_depth = 0;
      }
      else if (snow_shield_eff_depth.size() == 1)
      {
        this_top_eff_depth = snow_shield_eff_depth[0];
        //cout << "\n\nSnow shield depth: " <<   this_top_eff_depth << endl;
      }
      else
      {
        this_top_eff_depth = snow_shield_eff_depth[q];
        //cout << "\n\nSnow shield depth: " <<   this_top_eff_depth << endl;
      }

      // now get the self shielding. This is the thickness of the removed
      // layer
      if (self_shield_eff_depth.size() < 1)
      {
        this_bottom_eff_depth = this_top_eff_depth;
      }
      else if (self_shield_eff_depth.size() == 1)
      {
        this_bottom_eff_depth = this_top_eff_depth+self_shield_eff_depth[0];
        //cout << "\n\n Self shield depth: " << self_shield_eff_depth[0]  << endl;
      }
      else
      {
        this_bottom_eff_depth = this_top_eff_depth+self_shield_eff_depth[q];
        //cout << "\n\n Self shield depth: " << self_shield_eff_depth[q]  << endl;
      }

      // get the nuclide concentration from this node
      if (Nuclide == "Be10")
      {
        //cout << "LInE 2271, 10Be" << endl;
        eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                           this_top_eff_depth, this_bottom_eff_depth);
        Total_N+=eroded_particle.getConc_10Be();
      }
      else if (Nuclide == "Al26")
      {
        //cout << "LINE 2278, 26Al" << endl;
        eroded_particle.update_26Al_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                           this_top_eff_depth, this_bottom_eff_depth);
        Total_N+=eroded_particle.getConc_26Al();
      }
      else
      {
        cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
        cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
        cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
        eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                           this_top_eff_depth, this_bottom_eff_depth);
        Total_N+=eroded_particle.getConc_10Be();
      }

      //cout << endl << endl << "LINE 2052, total shield: " << total_shielding
      //     << " erosion: " << eff_erosion_rate << " erosion in cm/kyr with rho = 2650: "
      //     << eff_erosion_rate*1e6/2650.0 << " and N: " << eroded_particle.getConc_10Be() << endl;

    }
  }

  BasinAverage = Total_N/double(count_samples);
  average_production_rate = cumulative_production_rate/double(count_samples);
  average_production_uncertainty = average_production_rate*fabs(1-prod_uncert_factor);

  // replace the production uncertanty
  production_uncertainty = average_production_uncertainty;

  // replace the average production rate
  average_production = average_production_rate;

  return BasinAverage;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the concentration of a nuclide as  function of erosion rate
// The erosion rate should be in g/cm^2/yr
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoBasin::predict_mean_CRN_conc_with_snow_and_self_nested(double eff_erosion_rate,
                                            LSDRaster& known_effective_erosion,
                                            LSDFlowInfo& FlowInfo,
                                            string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            double& production_uncertainty,
                                            double& average_production,
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on)
{
  // production uncertainty factor is a multiplier that sets the production
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

  int end_node;
  end_node =  int(BasinNodes.size());


  // these parameters give the average production rate of the entore basin,
  // along with the magnitude of the production uncertainty
  double cumulative_production_rate = 0;
  double average_production_rate;
  double average_production_uncertainty;

  // the average atoms per gram of the nuclide
  double BasinAverage;

  // the total shielding. A product of snow, topographic and production scaling
  double total_shielding;
  double total_shielding_no_uncert;

  // the total atomic concentration of the nuclude in question
  double Total_N = 0;
  double Total_Mass = 0;

  int count_samples = 0;

  // initiate a particle. We'll just repeatedly call this particle
  // for the sample.
  int startType = 0;
  double Xloc = 0;
  double Yloc = 0;
  double  startdLoc = 0.0;
  double  start_effdloc = 0.0;
  double startzLoc = 0.0;

  int row,col;    // these are for getting the row and column from the know erosion rate raster

  // create a particle at zero depth
  LSDCRNParticle eroded_particle(startType, Xloc, Yloc,
                               startdLoc, start_effdloc, startzLoc);

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;

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

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }

  // parameters for the shielding
  double this_top_eff_depth;
  double this_bottom_eff_depth;

  // loop through the elevation data
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {
      count_samples++;

      // reset scaling parameters. This is necessary since the F values are
      // reset for local scaling
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

      // the elevation, snow shielding, topographic shielding
      // and production scaling are all independent of the erosion rate
      // and are calculated seperately.
      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // now you need logic to test if you are accounting for self shielding
      total_shielding_no_uncert = production_scaling[q]*topographic_shielding[q];
      total_shielding = prod_uncert_factor*total_shielding_no_uncert;
      cumulative_production_rate += total_shielding_no_uncert;

      // scale the F values
      LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

      // check to see if the shielding data exist and if so get the top and bottom
      // effective depths
      // first get snow shielding (as implemented by and effective depth of snow)
      if (snow_shield_eff_depth.size() < 1)
      {
        this_top_eff_depth = 0;
      }
      else if (snow_shield_eff_depth.size() == 1)
      {
        this_top_eff_depth = snow_shield_eff_depth[0];
        //cout << "\n\nSnow shield depth: " <<   this_top_eff_depth << endl;
      }
      else
      {
        this_top_eff_depth = snow_shield_eff_depth[q];
        //cout << "\n\nSnow shield depth: " <<   this_top_eff_depth << endl;
      }

      // now get the self shielding. This is the thickness of the removed
      // layer
      if (self_shield_eff_depth.size() < 1)
      {
        this_bottom_eff_depth = this_top_eff_depth;
      }
      else if (self_shield_eff_depth.size() == 1)
      {
        this_bottom_eff_depth = this_top_eff_depth+self_shield_eff_depth[0];
        //cout << "\n\n Self shield depth: " << self_shield_eff_depth[0]  << endl;
      }
      else
      {
        this_bottom_eff_depth = this_top_eff_depth+self_shield_eff_depth[q];
        //cout << "\n\n Self shield depth: " << self_shield_eff_depth[q]  << endl;
      }

      // get the nuclide concentration from this node
      // get the row and column of the node
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      // get the erosion rate from the raster
      float this_erosion_rate = known_effective_erosion.get_data_element(row,col);

      // check if there is a data value in the raster, if this is the case then
      // calculate the concentration based on the erosion from this raster
      if( this_erosion_rate != NoDataValue)
      {
        //cout << "This pixel at " << row << "," << col << " has an erate of: " << this_erosion_rate << endl;

        if (Nuclide == "Be10")
        {
          //cout << "LInE 2271, 10Be" << endl;
          eroded_particle.update_10Be_SSfull_depth_integrated(this_erosion_rate,LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
           Total_N+=this_erosion_rate*eroded_particle.getConc_10Be();
           Total_Mass+=this_erosion_rate;
        }
        else if (Nuclide == "Al26")
        {
          //cout << "LINE 2278, 26Al" << endl;
          eroded_particle.update_26Al_SSfull_depth_integrated(this_erosion_rate,LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          Total_N+=this_erosion_rate*eroded_particle.getConc_26Al();
          Total_Mass+=this_erosion_rate;
        }
        else
        {
          cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
          cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
          cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
          eroded_particle.update_10Be_SSfull_depth_integrated(this_erosion_rate,LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          Total_N+=this_erosion_rate*eroded_particle.getConc_10Be();
          Total_Mass+=this_erosion_rate;
        }
      }
      else  // this is the erosion rate not previously calculated (i.e., not in the raster)
      {
        if (Nuclide == "Be10")
        {
          //cout << "LInE 2271, 10Be" << endl;
          eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
           Total_N+=eff_erosion_rate*eroded_particle.getConc_10Be();
           Total_Mass+=eff_erosion_rate;
        }
        else if (Nuclide == "Al26")
        {
          //cout << "LINE 2278, 26Al" << endl;
          eroded_particle.update_26Al_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          Total_N+=eff_erosion_rate*eroded_particle.getConc_26Al();
          Total_Mass+=eff_erosion_rate;
        }
        else
        {
          cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
          cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
          cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
          eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                             this_top_eff_depth, this_bottom_eff_depth);
          Total_N+=eff_erosion_rate*eroded_particle.getConc_10Be();
          Total_Mass+=eff_erosion_rate;
        }
      }

    }
  }


  BasinAverage = Total_N/Total_Mass;
  average_production_rate = cumulative_production_rate/double(count_samples);
  average_production_uncertainty = average_production_rate*fabs(1-prod_uncert_factor);

  //cout << "Basin average is: " << BasinAverage << " and effective erosion: " << eff_erosion_rate << endl;


  // replace the production uncertanty
  production_uncertainty = average_production_uncertainty;

  // replace the average production rate
  average_production = average_production_rate;

  return BasinAverage;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the concentration of a nuclide as  function of erosion rate
// It calculates this based on the elevation scaling of the centroid
// The erosion rate should be in g/cm^2/yr
// NOTE: This is extremely inefficient: it is only here as a test to compare
//  how bad it is vs the full production scaling, but if we start using this more
// the averaging needs to be seperated from the concentration calculations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCosmoBasin::predict_mean_CRN_conc_centroid(double eff_erosion_rate, string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            LSDFlowInfo& FlowInfo,
                                            double& production_uncertainty,
                                            double& average_production,
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on)
{
  // production uncertainty factor is a multiplier that sets the production
  // certainty. If it is 1.1, there is 10% production rate uncertainty, or
  // if it is 0.9 there is -10% unvertainty. The reason why it is implemented
  // like this is that this allows gaussian error propigation.
  if (prod_uncert_factor <=0)
  {
    cout << "You have set an unrealistic production uncertainty factor." << endl;
    cout << "Defaulting to 1." << endl;
    prod_uncert_factor = 1;
  }

  // the average atoms per gram of the nuclide
  double AverageTopo;
  double AverageSnow;
  //double AverageProd;
  double AverageSelf;

  // the total shielding. A product of snow, topographic and production scaling
  double total_shielding;

  // the total atomic concentration of the nuclude in question
  double Total_N = 0;

  int count_samples = 0;

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

  // loop through the elevation data, averaging the snow and topo shielding
  double snow_shield_total = 0;
  double topo_shield_total = 0;
  double total_prod_scaling = 0;
  double self_shield_total = 0;
  int centroid_node = 0;   // if the centroid is not in the basin, the 'centroid'
                           // node defaults to the outlet
  int row,col;      // the row and column of the current node
  int end_node = int(BasinNodes.size());
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {
      count_samples++;

      // check to see if this is the centroid
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      if(row == Centroid_i && col == Centroid_j)
      {
        centroid_node = q;
      }

      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }
      snow_shield_total+= snow_shielding[q];
      topo_shield_total+= topographic_shielding[q];
      total_prod_scaling+= production_scaling[q];

      // check for self shielding
      if (self_shielding.size() > 1)
      {
        self_shield_total+= self_shielding[q];
      }
    }
  }

  AverageSnow = snow_shield_total/double(count_samples);
  AverageTopo = topo_shield_total/double(count_samples);
  //AverageProd = total_prod_scaling/double(count_samples);

  if (self_shielding.size() > 1)
  {
    AverageSelf = self_shield_total/double(count_samples);
  }
  else
  {
    AverageSelf = 1.0;
  }

  // at this stage we will try to replicate the basin averaging that goes on in
  // most paper
  // set scaling parameters. This is necessary since the F values are
  // reset for local scaling
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

  // check the production uncertainty bools
  if(is_production_uncertainty_plus_on)
  {
    if(is_production_uncertainty_minus_on)
    {
      cout << "You can't have both plus and minus production uncertainty on" << endl;
      cout << "Setting minus uncertainty to false" << endl;
      is_production_uncertainty_minus_on = false;
    }
  }

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
    cout << "LSDBasin line 1583, You didn't choos a valid nuclide. Defaulting"
         << " to 10Be." << endl;
    Nuclide = "Be10";
    nuclide_scaling_switches[0] = true;
  }


  // now get the shielding. This is based on the average snow sheilding,
  // the average topo shielding, and the production scaling of the centroid
  double total_shielding_no_uncert = AverageSnow*AverageTopo*AverageSelf*
                                     production_scaling[centroid_node];
  total_shielding = prod_uncert_factor*total_shielding_no_uncert;

  // get the uncertanty
  double average_production_uncertainty = total_shielding_no_uncert*(fabs(1-prod_uncert_factor));


  // now recalculate F values to match the total shielding
  LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

  // get the nuclide concentration from this node
  if (Nuclide == "Be10")
  {

    eroded_particle.update_10Be_SSfull(eff_erosion_rate,LSDCRNP);
    Total_N+=eroded_particle.getConc_10Be();
  }
  else if (Nuclide == "Al26")
  {
    eroded_particle.update_26Al_SSfull(eff_erosion_rate,LSDCRNP);
    Total_N+=eroded_particle.getConc_26Al();
  }
  else
  {
    cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
    cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
    cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
    eroded_particle.update_10Be_SSfull(eff_erosion_rate,LSDCRNP);
    Total_N+=eroded_particle.getConc_10Be();
  }

  // replace the production uncertanty number
  production_uncertainty = average_production_uncertainty;

  // replace the average production
  average_production = total_shielding;


  return Total_N;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Check nesting: this loops through an erosion rate raster to see if there
// are unknown erosion rates in the raster. It returns a boolean
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDCosmoBasin::are_there_unknown_erosion_rates_in_basin(LSDRaster& known_erates,LSDFlowInfo& FlowInfo)
{
  bool are_there_unknowns = false;
  int row,col;

  // get the end nod of the basin
  int end_node;
  end_node =  int(BasinNodes.size());

  // loop through the basin nodes. The logic stops as soon as it finds one unknown
  int q = 0;
  while( are_there_unknowns == false && q <= end_node)
  {
    // get the row and column of the node
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

    // get the erosion rate from the raster
    float this_erosion_rate = known_erates.get_data_element(row,col);

    // switch to true if you find an unknown
    if (this_erosion_rate == NoDataValue)
    {
      are_there_unknowns = true;
    }

    q++;

  }
  return are_there_unknowns;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets information of effective elevations for use in
// online calculators
//
// It returns a vector of values:
//  vector<double> parameter_returns;
//  parameter_returns.push_back(AverageProd);
//  parameter_returns.push_back(AverageTopo);
//  parameter_returns.push_back(AverageSelf);
//  parameter_returns.push_back(AverageSnow);
//  parameter_returns.push_back(AverageCombined);
//  parameter_returns.push_back(lat_outlet);
//  parameter_returns.push_back(outlet_pressure);
//  parameter_returns.push_back(outlet_eff_pressure);
//  parameter_returns.push_back(lat_centroid);
//  parameter_returns.push_back(centroid_pressure);
//  parameter_returns.push_back(centroid_eff_pressure);
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCosmoBasin::calculate_effective_pressures_for_calculators(LSDRaster& Elevation,
                        LSDFlowInfo& FlowInfo, string path_to_atmospheric_data)
{

  // the average atoms per gram of the nuclide
  double AverageTopo;
  double AverageProd;
  double AverageCombined;
  double AverageCombinedShielding;
  double AverageSnow;
  double AverageSelf;

  // the number of basin pixels
  int count_samples = 0;

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
  double gamma_spallation = 160;      // in g/cm^2: spallation attentuation depth

  // loop through the elevation data, averaging the snow and topo shielding
  double topo_shield_total = 0;
  double total_prod_scaling = 0;
  double self_shield_total = 0;
  double snow_shield_total = 0;
  double total_combined_scaling = 0;
  double total_combined_shielding = 0;
  double this_snow_shield, this_self_shield;
  int row,col;      // the row and column of the current node
  int end_node = int(BasinNodes.size());
  for (int q = 0; q < end_node; ++q)
  {
    //cout << "node " << q << " of " << end_node << endl;


    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {
      count_samples++;

      // check to see if this is the centroid
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // now get the snow shelding information
      if (snow_shield_eff_depth.size() < 1)
      {
        this_snow_shield = 1;
      }
      else if (snow_shield_eff_depth.size() == 1)
      {
        if (snow_shield_eff_depth[0] != 0)
        {
          this_snow_shield = exp(-snow_shield_eff_depth[0]/gamma_spallation);
        }
        else
        {
          this_snow_shield=1;
        }
      }
      else
      {
        if (snow_shield_eff_depth[q] != 0)
        {
          this_snow_shield = exp(-snow_shield_eff_depth[q]/gamma_spallation);
        }
        else
        {
          this_snow_shield=1;
        }
      }

      // now get the self shielding information
      if (self_shield_eff_depth.size() < 1)
      {
        this_self_shield = 1;
      }
      else if (self_shield_eff_depth.size() == 1)
      {
        if (self_shield_eff_depth[0] != 0)
        {
          this_self_shield = gamma_spallation/self_shield_eff_depth[0]*
                             (1-exp(-self_shield_eff_depth[0]/gamma_spallation));
        }
        else
        {
          this_self_shield = 1;
        }
      }
      else
      {
        if (self_shield_eff_depth[q] != 0)
        {
          this_self_shield = gamma_spallation/self_shield_eff_depth[q]*
                             (1-exp(-self_shield_eff_depth[q]/gamma_spallation));
        }
        else
        {
          this_self_shield=1;
        }
      }

      snow_shield_total += this_snow_shield;
      self_shield_total += this_self_shield;
      topo_shield_total += topographic_shielding[q];
      total_prod_scaling += production_scaling[q];
      total_combined_scaling += topographic_shielding[q]*production_scaling[q]*
                                this_snow_shield*this_self_shield;
      total_combined_shielding += topographic_shielding[q]*
                                this_snow_shield*this_self_shield;


    }
  }

  AverageTopo = topo_shield_total/double(count_samples);
  AverageProd = total_prod_scaling/double(count_samples);
  AverageSelf = self_shield_total/double(count_samples);
  AverageSnow = snow_shield_total/double(count_samples);
  AverageCombined = total_combined_scaling/double(count_samples);
  AverageCombinedShielding = total_combined_shielding/double(count_samples);

  // now find the latitude for both the outlet and the centroid
  // first the outlet
  double lat,longitude;
  double lat_centroid, long_centroid;
  double lat_outlet, long_outlet;
  double this_elevation;
  double centroid_pressure, outlet_pressure;
  double centroid_eff_pressure, outlet_eff_pressure;


  // declare converter object
  LSDCoordinateConverterLLandUTM Converter;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();

  Elevation.get_lat_and_long_locations(Outlet_i, Outlet_j, lat, longitude, Converter);
  lat_outlet = lat;
  long_outlet = longitude;
  Elevation.get_lat_and_long_locations(Centroid_i, Centroid_j, lat, longitude, Converter);
  lat_centroid = lat;
  long_centroid = longitude;

  // get outlet and centroid pressures
  this_elevation = Elevation.get_data_element(Centroid_i, Centroid_j);
  centroid_pressure = LSDCRNP.NCEPatm_2(double(lat_centroid), double(long_centroid),
                                        double(this_elevation));

  this_elevation = Elevation.get_data_element(Outlet_i, Outlet_j);
  outlet_pressure = LSDCRNP.NCEPatm_2(double(lat_outlet), double(long_outlet),
                                        double(this_elevation));

  // now we use newton iteration to calculate the 'effective' pressure for'
  // both the cnetroid and outlet latitutde.
  // First some variables that are used in the newton iteration
  double f_x,f_x_displace;
  double S_displace, S_this_step;
  double P_displace = 0.01;
  double P_change;
  double P_derivative;
  double tolerance = 1e-6;
  double Fsp = 0.978;

  // First for the centroid
  // initial guess is 1000hPa
  double this_P = 1000;
  lat = lat_centroid;
  do
  {
    S_this_step = LSDCRNP.stone2000sp(lat,this_P, Fsp);
    S_displace = LSDCRNP.stone2000sp(lat,this_P+P_displace, Fsp);

    f_x =  S_this_step - AverageProd;
    f_x_displace =  S_displace - AverageProd;

    P_derivative =  (f_x_displace-f_x)/P_displace;

    if(P_derivative != 0)
    {
      //cout << "Pressure before is: " <<this_P << " lat: " << lat;

      this_P = this_P-f_x/P_derivative;

      // check to see if the difference in erosion rates meet a tolerance
      P_change = f_x/P_derivative;
      //cout << " Change is: " << P_change << " target is: " << AverageProd << " and Shielding is: " << S_this_step << endl;

    }
    else
    {
      P_change = 0;
    }
  } while(fabs(P_change) > tolerance);
  centroid_eff_pressure = this_P;

  // Do it again for the outlet
  // initial guess is 1000hPa
  this_P = 1000;
  lat = lat_outlet;
  do
  {
    S_this_step = LSDCRNP.stone2000sp(lat,this_P, Fsp);
    S_displace = LSDCRNP.stone2000sp(lat,this_P+P_displace, Fsp);

    f_x =  S_this_step - AverageProd;
    f_x_displace =  S_displace - AverageProd;

    P_derivative =  (f_x_displace-f_x)/P_displace;

    if(P_derivative != 0)
    {
      this_P = this_P-f_x/P_derivative;

      // check to see if the difference in erosion rates meet a tolerance
      P_change = f_x/P_derivative;
      //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
    }
    else
    {
      P_change = 0;
    }
  } while(fabs(P_change) > tolerance);
  outlet_eff_pressure = this_P;

  vector<double> parameter_returns;
  parameter_returns.push_back(AverageProd);
  parameter_returns.push_back(AverageTopo);
  parameter_returns.push_back(AverageSelf);
  parameter_returns.push_back(AverageSnow);
  parameter_returns.push_back(AverageCombined);
  parameter_returns.push_back(lat_outlet);
  parameter_returns.push_back(outlet_pressure);
  parameter_returns.push_back(outlet_eff_pressure);
  parameter_returns.push_back(lat_centroid);
  parameter_returns.push_back(centroid_pressure);
  parameter_returns.push_back(centroid_eff_pressure);
  parameter_returns.push_back(AverageCombinedShielding);


  return parameter_returns;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets information of effective elevations for use in
// online calculators
//
// This version calculates the locations if there are known erosion rates,
//  used in nesting calculations
//
// It returns a vector of values:
//  vector<double> parameter_returns;
//  parameter_returns.push_back(AverageProd);
//  parameter_returns.push_back(AverageTopo);
//  parameter_returns.push_back(AverageSelf);
//  parameter_returns.push_back(AverageSnow);
//  parameter_returns.push_back(AverageCombined);
//  parameter_returns.push_back(lat_outlet);
//  parameter_returns.push_back(outlet_pressure);
//  parameter_returns.push_back(outlet_eff_pressure);
//  parameter_returns.push_back(lat_centroid);
//  parameter_returns.push_back(centroid_pressure);
//  parameter_returns.push_back(centroid_eff_pressure);
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCosmoBasin::calculate_effective_pressures_for_calculators_nested(LSDRaster& Elevation,
                        LSDFlowInfo& FlowInfo, string path_to_atmospheric_data,
                        LSDRaster& known_eff_erosion_rates)
{

  // make sure the rasters are the same size
  if ( Elevation.does_raster_have_same_dimensions(known_eff_erosion_rates) == false)
  {
    cout << "LSDCosmoBasin::calculate_effective_pressures_for_calculators_nested ERROR!" << endl;
    cout << "The erosion raster and DEM are not the same dimesions." <<endl;
    exit(EXIT_SUCCESS);
  }


  // the average atoms per gram of the nuclide
  double AverageTopo;
  double AverageProd;
  double AverageCombined;
  double AverageCombinedShielding;
  double AverageSnow;
  double AverageSelf;

  // the number of basin pixels
  int count_samples = 0;

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
  double gamma_spallation = 160;      // in g/cm^2: spallation attentuation depth

  // loop through the elevation data, averaging the snow and topo shielding
  double topo_shield_total = 0;
  double total_prod_scaling = 0;
  double self_shield_total = 0;
  double snow_shield_total = 0;
  double total_combined_scaling = 0;
  double total_combined_shielding = 0;
  double this_snow_shield, this_self_shield;
  int row,col;      // the row and column of the current node
  int end_node = int(BasinNodes.size());
  for (int q = 0; q < end_node; ++q)
  {
    // exclude known erosion rate locations from average
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);
    if(known_eff_erosion_rates.get_data_element(row,col) == NoDataValue)
    {
      //exclude NDV from average
      if(topographic_shielding[q] != NoDataValue)
      {
        count_samples++;

        if(  production_scaling.size() < 1 )
        {
          cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
               << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
        }

        // now get the snow shelding information
        if (snow_shield_eff_depth.size() < 1)
        {
          this_snow_shield = 1;
        }
        else if (snow_shield_eff_depth.size() == 1)
        {
          if (snow_shield_eff_depth[0] != 0)
          {
            this_snow_shield = exp(-snow_shield_eff_depth[0]/gamma_spallation);
          }
          else
          {
            this_snow_shield=1;
          }
        }
        else
        {
          if (snow_shield_eff_depth[q] != 0)
          {
            this_snow_shield = exp(-snow_shield_eff_depth[q]/gamma_spallation);
          }
          else
          {
            this_snow_shield=1;
          }
        }

        // now get the self shielding information
        if (self_shield_eff_depth.size() < 1)
        {
          this_self_shield = 1;
        }
        else if (self_shield_eff_depth.size() == 1)
        {
          if (self_shield_eff_depth[0] != 0)
          {
            this_self_shield = gamma_spallation/self_shield_eff_depth[0]*
                               (1-exp(-self_shield_eff_depth[0]/gamma_spallation));
          }
          else
          {
            this_self_shield = 1;
          }
        }
        else
        {
          if (self_shield_eff_depth[q] != 0)
          {
            this_self_shield = gamma_spallation/self_shield_eff_depth[q]*
                               (1-exp(-self_shield_eff_depth[q]/gamma_spallation));
          }
          else
          {
            this_self_shield=1;
          }
        }

        snow_shield_total += this_snow_shield;
        self_shield_total += this_self_shield;
        topo_shield_total += topographic_shielding[q];
        total_prod_scaling += production_scaling[q];
        total_combined_scaling += topographic_shielding[q]*production_scaling[q]*
                                  this_snow_shield*this_self_shield;
        total_combined_shielding += topographic_shielding[q]*
                                  this_snow_shield*this_self_shield;


      }
    }
  }

  AverageTopo = topo_shield_total/double(count_samples);
  AverageProd = total_prod_scaling/double(count_samples);
  AverageSelf = self_shield_total/double(count_samples);
  AverageSnow = snow_shield_total/double(count_samples);
  AverageCombined = total_combined_scaling/double(count_samples);
  AverageCombinedShielding = total_combined_shielding/double(count_samples);

  // now find the latitude for both the outlet and the centroid
  // first the outlet
  double lat,longitude;
  double lat_centroid, long_centroid;
  double lat_outlet, long_outlet;
  double this_elevation;
  double centroid_pressure, outlet_pressure;
  double centroid_eff_pressure, outlet_eff_pressure;


  // declare converter object
  LSDCoordinateConverterLLandUTM Converter;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();

  Elevation.get_lat_and_long_locations(Outlet_i, Outlet_j, lat, longitude, Converter);
  lat_outlet = lat;
  long_outlet = longitude;
  Elevation.get_lat_and_long_locations(Centroid_i, Centroid_j, lat, longitude, Converter);
  lat_centroid = lat;
  long_centroid = longitude;

  // get outlet and centroid pressures
  this_elevation = Elevation.get_data_element(Centroid_i, Centroid_j);
  centroid_pressure = LSDCRNP.NCEPatm_2(double(lat_centroid), double(long_centroid),
                                        double(this_elevation));

  this_elevation = Elevation.get_data_element(Outlet_i, Outlet_j);
  outlet_pressure = LSDCRNP.NCEPatm_2(double(lat_outlet), double(long_outlet),
                                        double(this_elevation));

  // now we use newton iteration to calculate the 'effective' pressure for'
  // both the cnetroid and outlet latitutde.
  // First some variables that are used in the newton iteration
  double f_x,f_x_displace;
  double S_displace, S_this_step;
  double P_displace = 0.01;
  double P_change;
  double P_derivative;
  double tolerance = 1e-6;
  double Fsp = 0.978;

  // First for the centroid
  // initial guess is 1000hPa
  double this_P = 1000;
  lat = lat_centroid;
  do
  {
    S_this_step = LSDCRNP.stone2000sp(lat,this_P, Fsp);
    S_displace = LSDCRNP.stone2000sp(lat,this_P+P_displace, Fsp);

    f_x =  S_this_step - AverageProd;
    f_x_displace =  S_displace - AverageProd;

    P_derivative =  (f_x_displace-f_x)/P_displace;

    if(P_derivative != 0)
    {
      //cout << "Pressure before is: " <<this_P << " lat: " << lat;

      this_P = this_P-f_x/P_derivative;

      // check to see if the difference in erosion rates meet a tolerance
      P_change = f_x/P_derivative;
      //cout << " Change is: " << P_change << " target is: " << AverageProd << " and Shielding is: " << S_this_step << endl;

    }
    else
    {
      P_change = 0;
    }
  } while(fabs(P_change) > tolerance);
  centroid_eff_pressure = this_P;

  // Do it again for the outlet
  // initial guess is 1000hPa
  this_P = 1000;
  lat = lat_outlet;
  do
  {
    S_this_step = LSDCRNP.stone2000sp(lat,this_P, Fsp);
    S_displace = LSDCRNP.stone2000sp(lat,this_P+P_displace, Fsp);

    f_x =  S_this_step - AverageProd;
    f_x_displace =  S_displace - AverageProd;

    P_derivative =  (f_x_displace-f_x)/P_displace;

    if(P_derivative != 0)
    {
      this_P = this_P-f_x/P_derivative;

      // check to see if the difference in erosion rates meet a tolerance
      P_change = f_x/P_derivative;
      //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
    }
    else
    {
      P_change = 0;
    }
  } while(fabs(P_change) > tolerance);
  outlet_eff_pressure = this_P;

  vector<double> parameter_returns;
  parameter_returns.push_back(AverageProd);
  parameter_returns.push_back(AverageTopo);
  parameter_returns.push_back(AverageSelf);
  parameter_returns.push_back(AverageSnow);
  parameter_returns.push_back(AverageCombined);
  parameter_returns.push_back(lat_outlet);
  parameter_returns.push_back(outlet_pressure);
  parameter_returns.push_back(outlet_eff_pressure);
  parameter_returns.push_back(lat_centroid);
  parameter_returns.push_back(centroid_pressure);
  parameter_returns.push_back(centroid_eff_pressure);
  parameter_returns.push_back(AverageCombinedShielding);


  return parameter_returns;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis prints shielding and scaling rasters.
// It uses the full extent of the raster so is storage inefficient.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::print_scaling_and_shielding_rasters(string filename,
                                          LSDFlowInfo& FlowInfo)
{
  // create arrays that will be used to generate rasters
  Array2D<float> Production_Data(NRows,NCols,NoDataValue);
  Array2D<float> Combined_Shielding_Data(NRows,NCols,NoDataValue);
  Array2D<float> Combined_Scaling_Data(NRows,NCols,NoDataValue);

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;
  double gamma_spallation = 160;      // in g/cm^2: spallation attentuation depth

  double this_snow_shield, this_self_shield;
  int row,col;      // the row and column of the current node
  int end_node = int(BasinNodes.size());
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {

      // get the row and column
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // now get the snow shelding information
      if (snow_shield_eff_depth.size() < 1)
      {
        this_snow_shield = 1;
      }
      else if (snow_shield_eff_depth.size() == 1)
      {
        if (self_shield_eff_depth[0] != 0)
        {
          this_snow_shield = exp(-snow_shield_eff_depth[0]/gamma_spallation);
        }
        else
        {
          this_snow_shield=1;
        }
      }
      else
      {
        if (snow_shield_eff_depth[q] != 0)
        {
          this_snow_shield = exp(-snow_shield_eff_depth[q]/gamma_spallation);
        }
        else
        {
          this_snow_shield=1;
        }
      }

      // now get the self shelding information
      if (self_shield_eff_depth.size() < 1)
      {
        this_self_shield = 1;
      }
      else if (self_shield_eff_depth.size() == 1)
      {
        if (self_shield_eff_depth[0] != 0)
        {
          this_self_shield = gamma_spallation/self_shield_eff_depth[0]*
                             (1-exp(-self_shield_eff_depth[0]/gamma_spallation));
        }
        else
        {
          this_self_shield = 1;
        }
      }
      else
      {
        if (self_shield_eff_depth[q] != 0)
        {
          this_self_shield = gamma_spallation/self_shield_eff_depth[q]*
                             (1-exp(-self_shield_eff_depth[q]/gamma_spallation));
        }
        else
        {
          this_self_shield=1;
        }
      }


      Production_Data[row][col] = production_scaling[q];
      Combined_Shielding_Data[row][col] = topographic_shielding[q]*
                                this_snow_shield*this_self_shield;
      Combined_Scaling_Data[row][col] = topographic_shielding[q]*production_scaling[q]*
                                this_snow_shield*this_self_shield;

    }
  }

  // now create and write the rasters
  string bil_ext = "bil";
  cout << "Printing production raster" << endl;
  LSDRaster ProductionRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Production_Data, GeoReferencingStrings);
  string Prod_ext = "_PROD";
  ProductionRaster.write_raster(filename+Prod_ext,bil_ext);

  cout << "Printing CombinedShielding raster" << endl;
  LSDRaster CombinedShieldingRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Combined_Shielding_Data, GeoReferencingStrings);
  string CombShield_ext = "_CSHIELD";
  CombinedShieldingRaster.write_raster(filename+CombShield_ext,bil_ext);

  cout << "Printing CombinedSscaling raster" << endl;
  LSDRaster CombinedScalingRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Combined_Scaling_Data, GeoReferencingStrings);
  string CombScale_ext = "_CSCALE";
  CombinedScalingRaster.write_raster(filename+CombScale_ext,bil_ext);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This returns the combined scaling raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDCosmoBasin::get_combined_scaling_raster(string filename,
                                          LSDFlowInfo& FlowInfo)
{
  // create arrays that will be used to generate rasters
  Array2D<float> Combined_Scaling_Data(NRows,NCols,NoDataValue);

  // now create the CRN parameters object
  LSDCRNParameters LSDCRNP;
  double gamma_spallation = 160;      // in g/cm^2: spallation attentuation depth

  double this_snow_shield, this_self_shield;
  int row,col;      // the row and column of the current node
  int end_node = int(BasinNodes.size());
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {

      // get the row and column
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // now get the snow shelding information
      if (snow_shield_eff_depth.size() < 1)
      {
        this_snow_shield = 1;
      }
      else if (snow_shield_eff_depth.size() == 1)
      {
        if (self_shield_eff_depth[0] != 0)
        {
          this_snow_shield = exp(-snow_shield_eff_depth[0]/gamma_spallation);
        }
        else
        {
          this_snow_shield=1;
        }
      }
      else
      {
        if (self_shield_eff_depth[q] != 0)
        {
          this_snow_shield = exp(-snow_shield_eff_depth[q]/gamma_spallation);
        }
        else
        {
          this_snow_shield=1;
        }
      }

      // now get the self shelding information
      if (self_shield_eff_depth.size() < 1)
      {
        this_self_shield = 1;
      }
      else if (self_shield_eff_depth.size() == 1)
      {
        if (self_shield_eff_depth[0] != 0)
        {
          this_self_shield = gamma_spallation/self_shield_eff_depth[0]*
                             (1-exp(-self_shield_eff_depth[0]/gamma_spallation));
        }
        else
        {
          this_self_shield = 1;
        }
      }
      else
      {
        if (self_shield_eff_depth[q] != 0)
        {
          this_self_shield = gamma_spallation/self_shield_eff_depth[q]*
                             (1-exp(-self_shield_eff_depth[q]/gamma_spallation));
        }
        else
        {
          this_self_shield=1;
        }
      }

      Combined_Scaling_Data[row][col] = topographic_shielding[q]*production_scaling[q]*
                                this_snow_shield*this_self_shield;


    }
  }

  // now create and write the raster
  LSDRaster CombinedScalingRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Combined_Scaling_Data, GeoReferencingStrings);

  return CombinedScalingRaster;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints a raster containing all the concentrations of
// nuclides predicted for the extracted erosion rate
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::print_CRN_conc_raster(string filename,
                                            double eff_erosion_rate, string Nuclide,
                                            string Muon_scaling, LSDFlowInfo& FlowInfo)
{

  // the array to which the data will be printed
  Array2D<float> Conc_Data(NRows,NCols,NoDataValue);

  // the row and col
  int row,col;

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

  double total_shielding;

  // set a special case if the outlet flag is true
  int end_node =  int(BasinNodes.size());

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

  // parameters for the shielding
  double this_top_eff_depth;
  double this_bottom_eff_depth;

  // loop through the elevation data
  for (int q = 0; q < end_node; ++q)
  {

    //exclude NDV from average
    if(topographic_shielding[q] != NoDataValue)
    {

      // get the row and column
      FlowInfo.retrieve_current_row_and_col(BasinNodes[q], row, col);

      // reset scaling parameters. This is necessary since the F values are
      // reset for local scaling
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

      // the elevation, snow shielding, topographic shielding
      // and production scaling are all independent of the erosion rate
      // and are calculated seperately.
      if(  production_scaling.size() < 1 )
      {
        cout << "LSDCosmoBasin, trying to precalculate erosion rate." << endl
             << "Scaling vectors have not been set! You are about to get a seg fault" << endl;
      }

      // get the scaling
      total_shielding = production_scaling[q]*topographic_shielding[q];

      // scale the F values
      LSDCRNP.scale_F_values(total_shielding,nuclide_scaling_switches);

      // first get snow shielding (as implemented by and effective depth of snow)
      if (snow_shield_eff_depth.size() < 1)
      {
        this_top_eff_depth = 0;
      }
      else if (snow_shield_eff_depth.size() == 1)
      {
        this_top_eff_depth = snow_shield_eff_depth[0];
        //cout << "\n\nSnow shield depth: " <<   this_top_eff_depth << endl;
      }
      else
      {
        this_top_eff_depth = snow_shield_eff_depth[q];
        //cout << "\n\nSnow shield depth: " <<   this_top_eff_depth << endl;
      }

      // now get the self shielding. This is the thickness of the removed
      // layer
      if (self_shield_eff_depth.size() < 1)
      {
        this_bottom_eff_depth = this_top_eff_depth;
      }
      else if (self_shield_eff_depth.size() == 1)
      {
        this_bottom_eff_depth = this_top_eff_depth+self_shield_eff_depth[0];
        //cout << "\n\n Self shield depth: " << self_shield_eff_depth[0]  << endl;
      }
      else
      {
        this_bottom_eff_depth = this_top_eff_depth+self_shield_eff_depth[q];
        //cout << "\n\n Self shield depth: " << self_shield_eff_depth[q]  << endl;
      }

      // get the nuclide concentration from this node
      if (Nuclide == "Be10")
      {
        //cout << "LInE 2271, 10Be" << endl;
        eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                           this_top_eff_depth, this_bottom_eff_depth);
        Conc_Data[row][col] = eroded_particle.getConc_10Be();
      }
      else if (Nuclide == "Al26")
      {
        //cout << "LINE 2278, 26Al" << endl;
        eroded_particle.update_26Al_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                           this_top_eff_depth, this_bottom_eff_depth);
        Conc_Data[row][col] = eroded_particle.getConc_26Al();
      }
      else
      {
        cout << "You didn't give a valid nuclide. You chose: " << Nuclide << endl;
        cout << "Choices are 10Be, 26Al.  Note these case sensitive and cannot" << endl;
        cout << "contain spaces or control characters. Defaulting to 10Be." << endl;
        eroded_particle.update_10Be_SSfull_depth_integrated(eff_erosion_rate,LSDCRNP,
                                           this_top_eff_depth, this_bottom_eff_depth);
        Conc_Data[row][col] = eroded_particle.getConc_10Be();
      }
    }
  }

  // now write the raster
  string Conc_ext = "_CONC";
  string bil_ext = "bil";

  cout << "Printing concentration raster" << endl;
  LSDRaster ConcRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Conc_Data, GeoReferencingStrings);
  ConcRaster.write_raster(filename+Conc_ext,bil_ext);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints the information, node by node, in the basin.
// Information is printed as a csv file
// It is used for bug checking
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCosmoBasin::print_particle_csv(string path_to_file, string filename,
                                       LSDFlowInfo& FlowInfo,
                                       LSDRaster& Elevation_Data,
                                       LSDRaster& T_Shield,
                                       string path_to_atmospheric_data)
{
  string full_name = path_to_file+filename+".csv";
  ofstream cosmo_out;
  cosmo_out.open(full_name.c_str());
  cosmo_out.precision(8);

  cosmo_out << "fID,Easting,Northing,Latitude,Longitude,Elevation,Pressure,TopoShield,Production_scaling,Snowshield" << endl;

  // check to see if scaling vecotrs have been made
  int n_nodes = int(BasinNodes.size());

  if (n_nodes != int(production_scaling.size()))
  {
    cout << "LSDCosmoBasin Line 1119, printing node info to csv but am getting shielding first." << endl;
    populate_scaling_vectors(FlowInfo, Elevation_Data, T_Shield, path_to_atmospheric_data);
  }

  // the latitude and longitude
  double lat,longitude;
  float Easting, Northing;
  double this_pressure,this_elevation,this_SShield,this_TShield,this_PShield;
  int row,col;

  // decalre converter object
  LSDCoordinateConverterLLandUTM Converter;

  // get the CRN parameters
  LSDCRNParameters LSDCRNP;

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();

  // now loop through nodes, printing the location and scaling
  for(int n = 0; n < n_nodes; n++)
  {
    FlowInfo.retrieve_current_row_and_col(BasinNodes[n], row, col);

    //exclude NDV from average
    if (Elevation_Data.get_data_element(row,col) != NoDataValue)
    {
      // To get pressure, first get the lat and long
      Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude, Converter);
      Elevation_Data.get_x_and_y_locations(row, col, Easting, Northing);
      //Elevation_Data.get_lat_and_long_locations(row, col, lat, longitude);


      // now the pressure
      this_elevation = Elevation_Data.get_data_element(row,col);
      this_pressure = LSDCRNP.NCEPatm_2(double(lat), double(longitude),
                                        double(this_elevation));

      // get the shielding
      this_TShield = topographic_shielding[n];
      this_SShield = snow_shielding[n];
      this_PShield = production_scaling[n];

      // now print to file
      cosmo_out << n+1 << ","<<Easting<<","<<Northing<<","<<lat<<","<<longitude
                <<","<<this_elevation<<","<<this_pressure<<","<<this_TShield<<","
                <<this_PShield<<","<<this_SShield<< endl;
    }
  }
}


#endif
