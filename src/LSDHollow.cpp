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
#include "LSDHollow.hpp"

using namespace std;
using namespace TNT;

#ifndef LSDHollow_CPP
#define LSDHollow_CPP


void LSDHollow::create(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet){

  //NO BOUNDS CHECKING ON JunctionNumber

  //setting all of the instance variables for the given junction

  NRows = ChanNet.get_NRows();
	NCols = ChanNet.get_NCols();
	XMinimum = ChanNet.get_XMinimum();
	YMinimum = ChanNet.get_YMinimum();
	DataResolution = ChanNet.get_DataResolution();
	NoDataValue = ChanNet.get_NoDataValue();

  Junction = JunctionNumber;
  
  vector <int> JunctionVector = ChanNet.get_JunctionVector();
  vector <int> ReceiverVector = ChanNet.get_ReceiverVector();
  
  LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction, JunctionVector[Junction],
                                                     ReceiverVector[Junction], JunctionVector[ReceiverVector[Junction]], FlowInfo);

  
  int hollow_outlet = StreamLinkVector.get_node_in_channel(0); // get hollow
  HollowNodes = FlowInfo.get_upslope_nodes(hollow_outlet);
                                                                                     
  NumberOfCells = int(HollowNodes.size());
  Area = NumberOfCells * (DataResolution*DataResolution);
  
  Beheaded = ChanNet.node_tester(FlowInfo, Junction);

  FlowInfo.retrieve_current_row_and_col(ChanNet.get_Node_of_Junction(Junction), Outlet_i, Outlet_j);

  int i_max = 0;
  int i_min = 9999999; //a very large number
  int j_max = 0;
  int j_min = 9999999; //a very large number
  
  int i = 0;
  int j = 0;

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
    if (i > i_max){i_max = i;}
    else if (i < i_min){i_min = i;}
    if (j > j_max){j_max = j;}
    else if (j < j_min){j_min = j;}
    
  }
  
  Centroid_i = i_min + ((i_max - i_min)/2);
  Centroid_j = j_min + ((j_max - j_min)/2);   //how do these handle 0.5s ??


  //finished setting all the instance variables
  
  
  // now we set up empty variables to store properties of the hollow
  // these are populated as they are required using the set methods in LSDHollow
    
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
  Perimeter_i = vector<int>(1,NoDataValue);
  Perimeter_j =  vector<int>(1,NoDataValue);
  BasalAge = NoDataValue;
  SoilProduction = NoDataValue;
  CHTMean = NoDataValue;
  Width = NoDataValue;
   
  //finished creating empty variables 

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate mean hollow value.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDHollow::CalculateHollowMean(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float TotalData = 0;
  int CountNDV = 0;
  float HollowAverage;

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
    //exclude NDV from average
    if (Data.get_data_element(i,j) != NoDataValue){
      TotalData += Data.get_data_element(i,j);
    }
    else {
      ++CountNDV;
    }
  }

  HollowAverage = TotalData/(NumberOfCells-CountNDV);

  return HollowAverage;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate max hollow value.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDHollow::CalculateHollowMax(LSDFlowInfo& FlowInfo, LSDRaster Data){

  //could use max_element here? how would that cope with NDVs??

  int i;
  int j;
  float MaxData = 0;
  float CurrentData;

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);
    
    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData > MaxData){
      MaxData = CurrentData;     
    }
  }
    
  return MaxData;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get the raster data passed out as an array of floats in the shape of the hollow.
// SWDG 20/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D<float> LSDHollow::get_Raster_Data_For_Hollow(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  Array2D<float> HollowData(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
    //exclude NDV 
    if (Data.get_data_element(i,j) != NoDataValue){
      HollowData[i][j] = Data.get_data_element(i,j);
    }
  }

  return HollowData;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get the raster data passed out as an array of integers in the shape of the hollow.
// SWDG 20/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D<int> LSDHollow::get_Raster_Data_For_Hollow(LSDFlowInfo& FlowInfo, LSDIndexRaster Data){

  int i;
  int j;
  Array2D<int> HollowData(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
    //exclude NDV 
    if (Data.get_data_element(i,j) != NoDataValue){
      HollowData[i][j] = Data.get_data_element(i,j);
    }
  }

  return HollowData;
  
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the mean hollow aspect. Does not use the normal hollow mean method as angles
// need to be handled differently. 
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDHollow::set_AspectMean(LSDFlowInfo& FlowInfo, LSDRaster Aspect){

  int i;
  int j;
  float avg_r;
  float angle_r;
  float x_component = 0.0;
  float y_component = 0.0;
  int ndv_cell_count = 0;  

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
    if (Aspect.get_data_element(i,j) != NoDataValue){
    
      angle_r = rad(Aspect.get_data_element(i,j));
      x_component += cos(angle_r);
      y_component += sin(angle_r);
  
    }
    else{
      ++ndv_cell_count;
    }
  
  }
    
  x_component = x_component / (HollowNodes.size() - ndv_cell_count);
  y_component = x_component / (HollowNodes.size() - ndv_cell_count);
  avg_r = atan2(y_component, x_component);
  
  AspectMean = deg(avg_r);
   
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the perimeter pixels using a simple edge detection algorithm. This is quite 
// messy and will be improved soon.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDHollow::set_Perimeter(LSDFlowInfo& FlowInfo){

  int i;
  int j;
  vector<int> I;
  vector<int> J;
  int NDVCount = 0;
  Array2D<float> HollowData(NRows, NCols, NoDataValue);

  //create subset arrays for just the hollow data - this should be rolled into its own method.
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
      HollowData[i][j] = HollowNodes[q];
    
  }

  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
      NDVCount = 0;
     
        //count border cells that are NDV
        if (HollowData[i-1][j-1] == NoDataValue){ ++NDVCount; }
        if (HollowData[i][j-1] == NoDataValue){ ++NDVCount; }
        if (HollowData[i+1][j-1] == NoDataValue){ ++NDVCount; }
        if (HollowData[i-1][j] == NoDataValue){ ++NDVCount; }
        if (HollowData[i+1][j] == NoDataValue){ ++NDVCount; }
        if (HollowData[i-1][j+1] == NoDataValue){ ++NDVCount; }
        if (HollowData[i][j+1] == NoDataValue){ ++NDVCount; }
        if (HollowData[i+1][j+1] == NoDataValue){ ++NDVCount; }
        
        if (NDVCount >= 4 && NDVCount < 8){  //increase the first value to get a simpler polygon
          //edge pixel
          I.push_back(i);
          J.push_back(j);
        }
    
  }

  //now have 2 vectors of i and j indexes of every point
  Perimeter_i = I;
  Perimeter_j = J;


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set all of the hollow parameters with one call.
//
// Runs polyfit to get the elevation derivatives, so can be quite memory intensive. Method
// calls all the setters one by one, to populate all the hollow parameters. So a
// hollow can be created and all it's properties set with 2 calls. The BasalAge and SoilProduction
// have default parameters of -9999 as these are rarely used variables.
// SWDG 20/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDHollow::set_All_Parameters(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& CHT,
                                  LSDRaster& Relief, float window_radius,
                                  float SoilProduction, float BasalAge){

  // coefficent matrices for polyfit routine
  Array2D<float> a;
  Array2D<float> b;
  Array2D<float> c;
  Array2D<float> d;
  Array2D<float> e;
  Array2D<float> f;

  Elevation.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster TotalCurv = Elevation.calculate_polyfit_curvature(a,b);
  LSDRaster ProfileCurv = Elevation.calculate_polyfit_profile_curvature(a,b,c,d,e);
  LSDRaster PlanCurv = Elevation.calculate_polyfit_planform_curvature(a,b,c,d,e);
  LSDRaster Aspect = Elevation.calculate_polyfit_aspect(d,e);  
  LSDRaster Slope = Elevation.calculate_polyfit_slope(d,e);
  Array2D<float> FlowDir = Elevation.D_inf_FlowDir(); 
  
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
  set_Perimeter(FlowInfo);
  set_BasalAge(BasalAge);
  set_SoilProduction(SoilProduction);
  set_Width(FlowInfo, FlowDir);
  set_DownslopeLength(FlowInfo, Elevation);
  set_LongProfileLength(FlowInfo);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write integer hollow parameters into the shape of the hollow.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDHollow::write_integer_data_to_LSDIndexRaster(int Param, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<int> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    Output[i][j] = Param;
  }

  LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}                                       

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write real hollow parameters into the shape of the hollow.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDHollow::write_real_data_to_LSDRaster(float Param, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<float> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    Output[i][j] = Param;
  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cookie cut data from an LSDRaster into the shape of the hollow.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDHollow::write_raster_data_to_LSDRaster(LSDRaster Data, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<float> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    Output[i][j] = Data.get_data_element(i,j);
  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cookie cut data from an LSDIndexRaster into the shape of the hollow.
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDHollow::write_raster_data_to_LSDIndexRaster(LSDIndexRaster Data, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<int> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    Output[i][j] = Data.get_data_element(i,j);
  }

  LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Measure the width of a hollow, defined as the width perpendicular to the 
// centroid flow direction. 
// SWDG 19/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDHollow::set_Width(LSDFlowInfo FlowInfo, Array2D<float> FlowDir){
         
  LSDIndexRaster hollow = write_Junction(FlowInfo);
  Array2D<int> HollowArray = hollow.get_RasterData();
    
  float centre_flowdir = FlowDir[Centroid_i][Centroid_j];
  
  float x2;
  float y2;
  vector<int> i_list;
  vector<int> j_list;
  
  int i_new = Centroid_i;
  int j_new = Centroid_j;
  
  i_list.push_back(i_new);
  j_list.push_back(j_new);
  
  float x1 = i_new + 0.5;
  float y1 = j_new - 0.5;
   
  int x_top = 0;
  int y_top = 0;
  
  //get perpendicular flowdirs
  float perp_angle_1 = centre_flowdir - 90;
  float perp_angle_2 = centre_flowdir + 90;
  
  if (perp_angle_1 < 0) {perp_angle_1 = perp_angle_1 + 360;}
   
  if (perp_angle_2 > 360) {perp_angle_2 = perp_angle_2 - 360;} 
  
       
    while (HollowArray[i_new][j_new] != NoDataValue){
      
      x2 = x1 + cos(rad(perp_angle_1)) * DataResolution;
      y2 = y1 - sin(rad(perp_angle_1)) * DataResolution;
        
      i_new = trunc(x2);
      j_new = ceil(y2);  
        
      i_list.push_back(i_new);
      j_list.push_back(j_new);
        
      x1 = x2;
      y1 = y2;

    }
     
    i_new = Centroid_i;
    j_new = Centroid_j;
  
    x1 = i_new + 0.5;
    y1 = j_new - 0.5;
    
    x_top = i_list[i_list.size()-1];
    y_top = j_list[j_list.size()-1];
  
    while (HollowArray[i_new][j_new] != NoDataValue){
    
      x2 = x1 + cos(rad(perp_angle_2)) * DataResolution;
      y2 = y1 - sin(rad(perp_angle_2)) * DataResolution;
    
      i_new = trunc(x2);
      j_new = ceil(y2);  
      
      i_list.push_back(i_new);
      j_list.push_back(j_new);
      
      x1 = x2;
      y1 = y2;

    }
  
  Array2D<int> out(NRows, NCols, NoDataValue);

  for (int q = 0; q < int(i_list.size()); ++q){
    out[i_list[q]][j_list[q]] = 1;
  
  }
  
  int x_bottom = i_list[i_list.size()-1];
  int y_bottom = j_list[j_list.size()-1];
  
  float len = sqrt( ((x_top - x_bottom) * (x_top - x_bottom)) + ((y_top - y_bottom) * (y_top - y_bottom)) );
  
  Width = len;
    
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Measure the downslope length in the hollow, defined as the D8 flow routing distance
// from between the maximum elevation point in the hollow to the minimum elevation, or 
// to the edge of the hollow, whichever comes first.    
//  
// SWDG 20/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDHollow::set_DownslopeLength(LSDFlowInfo FlowInfo, LSDRaster DEM){

  int i;
  int j;
  float MinData = 100000000.0; // a large number
  float MaxData = 0.0; // a small number number
  float CurrentData;
  
  int Min_i = 0;
  int Min_j = 0;
  int Max_i = 0;
  int Max_j = 0;
  
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    CurrentData = DEM.get_data_element(i,j);
    
    if (CurrentData != NoDataValue && CurrentData < MinData){
      MinData = CurrentData;
      Min_i = i;
      Min_j = j;     
    }
    if (CurrentData != NoDataValue && CurrentData > MaxData){
      MaxData = CurrentData;
      Max_i = i;
      Max_j = j;     
    }
  }
  
  //flow from max i,j to min i,j
  
  float root_2 = 1.4142135623;
  float length = 0;
  int node;
  int receiver_node = FlowInfo.retrieve_node_from_row_and_column(Max_i, Max_j);
  int receiver_row = Max_i;
  int receiver_col = Max_j;
  
  LSDIndexRaster Hollow = write_Junction(FlowInfo);
   
  while ((receiver_row != Min_i && receiver_col != Min_j) || Hollow.get_data_element(receiver_row, receiver_col) != NoDataValue){
    
    FlowInfo.retrieve_receiver_information(receiver_node, node, receiver_row, receiver_col);
  
    //update length
    if (FlowInfo.retrieve_flow_length_code_of_node(receiver_node) == 1){ length += DataResolution; }
    else if (FlowInfo.retrieve_flow_length_code_of_node(receiver_node) == 1){ length += (DataResolution * root_2); }
      
    receiver_node = node;
   
  }
  
  DownslopeLength = length;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Measure the long profile length in the hollow, defined as the maximum dimension 
// of the bounding box excluding diagonals.
//  
// SWDG 20/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDHollow::set_LongProfileLength(LSDFlowInfo FlowInfo){

  int i;
  int j; 
  int Min_i = 10000000; //a very large number
  int Min_j = 10000000; //a very large number
  int Max_i = 0;
  int Max_j = 0;
  
  
  LSDIndexRaster Hollow = write_Junction(FlowInfo);
  
  for (int q = 0; q < int(HollowNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(HollowNodes[q], i, j);
    
    if (Hollow.get_data_element(i,j) != NoDataValue){
    
      if (i > Max_i) { Max_i = i; }
      if (i < Min_i) { Min_i = i; }
      if (j > Max_j) { Max_j = j; }
      if (j < Min_j) { Min_j = i; }
       
    }
  }

  if ((Max_i - Min_i) > (Max_j - Min_j)) { LongProfileLength = Max_i - Min_i;}
  if ((Max_i - Min_i) < (Max_j - Min_j)) { LongProfileLength = Max_j - Min_j;}

}

#endif