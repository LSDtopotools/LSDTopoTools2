//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDSwathProfile.cpp
//------------------------------------------------------------------------------
// This code is used to create swath profiles from either raster or sparse data.
// These two options will be coded up into two different classes, in order to
// deal with different types of data from which to construct the profiles.
//
// The generalised swath profile framework for constructing transverse profiles
// is derived from the algorithm described by Hergarten et al. 2013: Generalized
// swath pro?les. This will be extended to also produce generalised longitudinal
// profiles, as well as some other functionality as desired.
//
// The input to the swath profile objects includes the profile itself, which is
// loaded as a PointData structure as described in LSDShapeTools.hpp.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 0.0.1		17/02/2014
// Prerequisite software packages: TNT, PCL and liblas
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

// LSDTopotools objects
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDCloudRaster.hpp"
#include "LSDShapeTools.hpp"
#include "LSDSwathProfile.hpp"
#include "LSDStatsTools.hpp"
#include "LSDSpatialCSVReader.hpp"

// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
//#include <liblas/liblas.hpp>

// TNT
#include "TNT/tnt.h"

using namespace std;
using namespace TNT;

#ifndef LSDSwathProfile_CPP
#define LSDSwathProfile_CPP

//------------------------------------------------------------------------------
// VECTOR GEOMETRY FUNCTIONS
// Functions to solve some geometric problems, such as shortest distance from a
// point to a line etc.

// This function calculates the shortest distance from a point, V_p, to a
// straight line passing through points V_a and V_b, using the vector triple
// product.
float calculate_shortest_distance_to_line(vector<float> V_a, vector<float> V_b, vector<float> V_p)
{
  float d = sqrt(( (V_b[0]-V_a[0])*(V_a[1]-V_p[1]) - (V_b[1]-V_a[1])*(V_a[0]-V_p[0]) )*( (V_b[0]-V_a[0])*(V_a[1]-V_p[1]) - (V_b[1]-V_a[1])*(V_a[0]-V_p[0]) ))
                                                    / sqrt( ((V_b[0]-V_a[0])*(V_b[0]-V_a[0])) + ((V_b[1]-V_a[1])*(V_b[1]-V_a[1])) );
//  cout << "V_a " << V_a[0] << "," << V_a[1] << "\tV_b " << V_b[0] << "," << V_b[1] << "\tV_p " << V_p[0] << "," << V_p[1] << "\td " << d << endl;
  return d;
}
// For a vector line equation V_star = V_a + (V_b-V_a)*t, this function
// calculates t to find the intersection point between the line joining points
// V_a and V_b and the shortest path from a point V_p and this line.
float calculate_t(vector<float> V_a, vector<float> V_b, vector<float> V_p)
{
  float t = -( (V_b[0]-V_a[0])*(V_a[0]-V_p[0]) + (V_b[1]-V_a[1])*(V_a[1]-V_p[1]) ) / ((V_b[0]-V_a[0])*(V_b[0]-V_a[0]) + (V_b[1]-V_a[1])*(V_b[1]-V_a[1]));
  return t;
}
// Calculate disctance between two points
float calculate_distance_between_two_points(vector<float> V_a, vector<float> V_b)
{
  float distance = sqrt((V_a[0]-V_b[0])*(V_a[0]-V_b[0])+(V_a[1]-V_b[1])*(V_a[1]-V_b[1]));
  return distance;
}

// Use cross product to test whether point lies on left or right hand side of
// a line vector in 2D.  Returns true if the point is on the left.
bool test_point_left(vector<float> V_a, vector<float> V_b, vector<float> V_p)
{
  bool test;
  float cross_product = (V_b[0]-V_a[0])*(V_p[1]-V_a[1]) - (V_b[1]-V_a[1])*(V_p[0]-V_a[0]);
  if(cross_product>0) test = true;
  else test = false;
}
//------------------------------------------------------------------------------
// CREATE FUNCTIONS
// the create function. This is default and throws an error
void LSDSwath::create()
{
  cout << "LSDSwathProfile line 102 You need to supply a PointData container of profile coordinates, a raster template and the half width of the profile" << endl;
  exit(EXIT_FAILURE);
}
// This function creates the swath profile template for creating a swath profile
// of a raster dataset.
void LSDSwath::create(PointData& ProfilePoints, LSDRaster& RasterTemplate, float HalfWidth)
{
  NPtsInProfile = ProfilePoints.X.size();   // Number of points in profile
  NoDataValue = RasterTemplate.get_NoDataValue();
  NRows = RasterTemplate.get_NRows();
  NCols = RasterTemplate.get_NCols();
  ProfileHalfWidth = HalfWidth;
  float Resolution = RasterTemplate.get_DataResolution();
  // Loop through profile points, calculating cumulative distance along profile
  // with each iteration
  vector<float> DistanceAlongBaseline_temp(NPtsInProfile,NoDataValue);
  vector<float> BaselineValue_temp(NPtsInProfile,NoDataValue);
  vector<int> BaselineRows_temp(NPtsInProfile,NoDataValue);
  vector<int> BaselineCols_temp(NPtsInProfile,NoDataValue);
  float cumulative_distance = 0;
  DistanceAlongBaseline_temp[0]=cumulative_distance;
  // Retrieve baseline value at each point - sample closest pixel in template raster
  float X_coordinate_shifted_origin = ProfilePoints.X[0] - RasterTemplate.get_XMinimum() - 0.5*Resolution;
  float Y_coordinate_shifted_origin = ProfilePoints.Y[0] - RasterTemplate.get_YMinimum() + 0.5*Resolution;
  // Get row and column of point
  int col_point = int(round(X_coordinate_shifted_origin/Resolution));
  int row_point = (NRows) - int(round(Y_coordinate_shifted_origin/Resolution));
  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows) BaselineValue_temp[0]=NoDataValue;
  else
  {
    BaselineValue_temp[0] = RasterTemplate.get_data_element(row_point,col_point);
    BaselineRows_temp[0] = row_point;
    BaselineCols_temp[0] = col_point;
  }


  for(int i = 1; i<NPtsInProfile; ++i)
  {
    // Calculate cumulative distance
    cumulative_distance += sqrt((ProfilePoints.X[i]-ProfilePoints.X[i-1])*(ProfilePoints.X[i]-ProfilePoints.X[i-1])
                                              +(ProfilePoints.Y[i]-ProfilePoints.Y[i-1])*(ProfilePoints.Y[i]-ProfilePoints.Y[i-1]));
    DistanceAlongBaseline_temp[i]=cumulative_distance;

    // Retrieve baseline value at each point - sample closest pixel in template raster
    X_coordinate_shifted_origin = ProfilePoints.X[i] - RasterTemplate.get_XMinimum() - 0.5*Resolution;
    Y_coordinate_shifted_origin = ProfilePoints.Y[i] - RasterTemplate.get_YMinimum() + 0.5*Resolution;

    // Get row and column of point
    col_point = int(round(X_coordinate_shifted_origin/Resolution));
    row_point = (NRows) - int(round(Y_coordinate_shifted_origin/Resolution));
    if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows) BaselineValue_temp[i]=NoDataValue;
    else
    {
      BaselineValue_temp[i] = RasterTemplate.get_data_element(row_point,col_point);
      BaselineRows_temp[i] = row_point;
      BaselineCols_temp[i] = col_point;
    }
  }

  DistanceAlongBaseline = DistanceAlongBaseline_temp;
  BaselineValue = BaselineValue_temp;

  // Read profile data into a LSDCloud object for querying
  vector<float> zeros(NPtsInProfile,0.0);
  LSDCloudRaster ProfileCloud(ProfilePoints, RasterTemplate);

  // For each point in array, find nearest point along the profile and calculate
  // signed distance to profile.  The convention here is that points lying on
  // left hand side of profile as you traverse from start to finish are
  // considered positive, and those on the right are negative.
  Array2D<float> DistanceToBaseline_temp(NRows,NCols,NoDataValue);
  Array2D<float> ProjectedDistanceAlongBaseline_temp(NRows,NCols,NoDataValue);
  Array2D<float> ProjectedBaselineValue_temp(NRows,NCols,NoDataValue);
  Array2D<float> BaselineValueArray_temp(NRows,NCols,NoDataValue);

  float searchPoint_x,searchPoint_y;
  vector<float> SquaredDistanceToBaseline,temp,temp2;
  vector<int> ProfilePointIndex;
  int K=2;

  // Define bounding box of swath profile
  YMin = ProfileCloud.get_YMin();
  YMax = ProfileCloud.get_YMax();
  XMin = ProfileCloud.get_XMin();
  XMax = ProfileCloud.get_XMax();
  int ColStart = int(floor(XMin/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    cout << flush << "\t\t\t row " << i+1 << "/" << NRows << "\r";
    for(int j = ColStart; j<ColEnd; ++j)
    {
      vector<float> V_p;  // the search point coordinates
      V_p.push_back(float(j)*Resolution);
      V_p.push_back((float(NRows - 1) - float(i))*Resolution);
      // Find nearest two points on profile baseline
      ProfileCloud.NearestNeighbourSearch2D(V_p[0], V_p[1], K, temp, ProfilePointIndex, SquaredDistanceToBaseline);

      if(SquaredDistanceToBaseline.size()<1) cout << "No profile points found - check input files";
      else if(SquaredDistanceToBaseline.size()<2) cout << "Only one point found - check input files";
      else
      {
        // SCENARIO 1
        //----------------------------------------------------------------------
        // Check that two points are adjacenet on baseline - simplest (and
        // usual) scenario)
        vector<float> V_a, V_b; // first, second nearest points on profile
        int V_a_index, V_b_index;
        if(abs(ProfilePointIndex[0]-ProfilePointIndex[1]) <= 1)
        {
//           cout << "SCENARIO 1" << endl;
          // If yes, then find point along cojoining line that minimises this
          // distance.
          // Note that this point is located at the point given by the vector:
          // (V_a + (V_b - V_a)*t).
          // t is calculated using the equation:
          // t = -(V_b-V_a).(V_a-V_p)/(|V_b - V_a|^2)
          V_a_index = ProfilePointIndex[0];
          V_b_index = ProfilePointIndex[1];
          V_a.push_back(ProfileCloud.get_point_x(V_a_index));
          V_a.push_back(ProfileCloud.get_point_y(V_a_index));
          V_b.push_back(ProfileCloud.get_point_x(V_b_index));
          V_b.push_back(ProfileCloud.get_point_y(V_b_index));
        }

        // SCENARIO 2
        //----------------------------------------------------------------------
        // Two nearest points are not adjacent on the profile.  This is quite
        // common on the inside of a bend on a baseline.  At this point, the
        // strategy used is to find all baseline points located at a radius that
        // is equal to that of the nearest baseline point.  Then, the strategy
        // is to search through all located baseline points, and find the
        // nearest point that neighbours one of these baseline points.  The
        // program then searches through to find the closest point on the
        // baseline between those points.
        else
        {
//           cout << "SCENARIO 2" << endl;//" " << sqrt(SquaredDistanceToBaseline[0]) << endl;
          vector<int> ClosestPointIndex;
          float shortest_distance = NoDataValue;
          float distance_to_point;
          //ProfileCloud.RadiusSearch2D(V_p[0], V_p[1], sqrt(SquaredDistanceToBaseline[0]), temp, ClosestPointIndex, temp2);
          ProfileCloud.NearestNeighbourSearch2D(V_p[0], V_p[1], 8, temp, ClosestPointIndex, temp2);
          for(int i_point_iter = 0; i_point_iter < ClosestPointIndex.size(); ++i_point_iter)
          {
            if((ClosestPointIndex[i_point_iter]!=0) && (ClosestPointIndex[i_point_iter]!=NPtsInProfile-1))
            {
              // Check step backwards
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]-1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]-1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]-1;
              }
              // Check step forwards
              V_b_test.clear();
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]+1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]+1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]+1;
              }
            }
            // conditions for edge of profile
            else if(ClosestPointIndex[i_point_iter]==0)
            {
              // Check step forwards only
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]+1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]+1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]+1;
              }
            }
            else if(ClosestPointIndex[i_point_iter]==NPtsInProfile-1)
            {
              // Check step backwards only
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]-1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]-1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]-1;
              }
            }
          }
        }
        // calculate position along baseline vector
        float t = calculate_t(V_a, V_b, V_p);
        float d;
        if(t>0 && t<1)
        // find the distance to the nearest point along the straight line
        // segment that links the two points
        {
//          cout << "test1" << endl;
          d = calculate_shortest_distance_to_line(V_a,V_b,V_p);
          if(d<ProfileHalfWidth)
          {
            ProjectedDistanceAlongBaseline_temp[i][j] = DistanceAlongBaseline[V_a_index]+(DistanceAlongBaseline[V_b_index]-DistanceAlongBaseline[V_a_index])*t;
            // determine which side of the profile the point is on using the cross
            // product
            if(V_a_index < V_b_index)
            {
              if(test_point_left(V_b,V_a,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
            else
            {
              if(test_point_left(V_a,V_b,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
            // Now get baseline value
            BaselineValueArray_temp[i][j] = BaselineValue[V_a_index]+(BaselineValue[V_b_index]-BaselineValue[V_a_index])*t;
          }
          else
          {
            DistanceToBaseline_temp[i][j] = NoDataValue;
            ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
            BaselineValueArray_temp[i][j] = NoDataValue;
          }
        }

        else if(V_a_index==0 || V_a_index==NPtsInProfile-1)
        // Avoid end points
        {
//          cout << "test2" << endl;
          DistanceToBaseline_temp[i][j] = NoDataValue;
          ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
          BaselineValueArray_temp[i][j] = NoDataValue;
        }
        else
        // Select nearest point (e.g. on outer side of bend)
        {
//          cout << "test3" << endl;
          d = calculate_distance_between_two_points(V_a,V_p);
          if(d<ProfileHalfWidth)
          {
            ProjectedDistanceAlongBaseline_temp[i][j] = DistanceAlongBaseline[V_a_index];
            BaselineValueArray_temp[i][j] = BaselineValue[V_a_index];
            // determine which side of the profile the point is on using the cross
            // product.  In this case, take the baseline vector as being formed
            // by straight line that joins the points on either side.
            vector<float> V_c(2,0.0);
            if(V_a_index < V_b_index)
            {
              V_c[0]=ProfileCloud.get_point_x(V_a_index-1);
              V_c[1]=ProfileCloud.get_point_y(V_a_index-1);
              if(test_point_left(V_b,V_c,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
            else
            {
              V_c[0]=ProfileCloud.get_point_x(V_a_index+1);
              V_c[1]=ProfileCloud.get_point_y(V_a_index+1);
              if(test_point_left(V_c,V_b,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
          }
          else
          {
            DistanceToBaseline_temp[i][j] = NoDataValue;
            ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
            BaselineValueArray_temp[i][j] = NoDataValue;
          }
        }
      }
    }
  }
  DistanceToBaselineArray = DistanceToBaseline_temp.copy();
  DistanceAlongBaselineArray = ProjectedDistanceAlongBaseline_temp.copy();
  BaselineValueArray = BaselineValueArray_temp.copy();
  BaselineRows = BaselineRows_temp;
  BaselineCols = BaselineCols_temp;
}

void LSDSwath::create(vector<float>& Y_X_points, LSDRaster& RasterTemplate, float& HalfWidth, float d_space)
{
  cout << "Creation of the swath from WGS coordinates" << endl;
  // initialization of the variables
  vector<vector<float> > points; // contains all the points of the profile
  PointData PD;
  vector<float> unity; // Unity vector
  vector<float> temp_vecta;
  float temp_Y, temp_X, temp_distance_AB;
  float d_point = d_space; // distance between points

  // calculating the unity vector
  temp_distance_AB = sqrt( pow((Y_X_points[3]-Y_X_points[1]),2) + pow((Y_X_points[2]-Y_X_points[0]),2)); // distance between A and B
  temp_X = ((Y_X_points[3]-Y_X_points[1])/temp_distance_AB);
  temp_Y = ((Y_X_points[2]-Y_X_points[0])/temp_distance_AB);
  unity.push_back(temp_Y);
  unity.push_back(temp_X);


  // generation of the points
  for(int i = 0; i < temp_distance_AB; i += d_point)
  {
    temp_X = Y_X_points[1] + unity[1]*i;
    temp_Y = Y_X_points[0] + unity[0]*i;
    temp_vecta.push_back(temp_Y);
    temp_vecta.push_back(temp_X);
    points.push_back(temp_vecta);
    temp_vecta.clear();
    PD.X.push_back(temp_X);
    PD.Y.push_back(temp_Y);
  }
  // pushing the last point
  temp_vecta.push_back(Y_X_points[2]);
  temp_vecta.push_back(Y_X_points[3]);
  points.push_back(temp_vecta);
  temp_vecta.clear();

  // test for development of the function, just to check if the points are written between A and B
  //ofstream output_file;
	//string output_fname = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/Basins/Olt_Basin/DEV_swath_elevations.csv";
	//output_file.open(output_fname.c_str());
	//output_file << "Y,X" << endl;
	//for (int i = 0; i < int(points.size()); ++i)
	//{
	//	output_file << points[i][0] << "," << points[i][1] << endl;
	//}
	//output_file.close();

  // Now run the same constructor as it should be

  PointData ProfilePoints = PD;
  NPtsInProfile = ProfilePoints.X.size();   // Number of points in profile
  NoDataValue = RasterTemplate.get_NoDataValue();
  NRows = RasterTemplate.get_NRows();
  NCols = RasterTemplate.get_NCols();
  ProfileHalfWidth = HalfWidth;
  float Resolution = RasterTemplate.get_DataResolution();
  // Loop through profile points, calculating cumulative distance along profile
  // with each iteration
  vector<float> DistanceAlongBaseline_temp(NPtsInProfile,NoDataValue);
  vector<float> BaselineValue_temp(NPtsInProfile,NoDataValue);
	vector<int> BaselineRows_temp(NPtsInProfile,NoDataValue);
	vector<int> BaselineCols_temp(NPtsInProfile,NoDataValue);
  float cumulative_distance = 0;
  DistanceAlongBaseline_temp[0]=cumulative_distance;
  // Retrieve baseline value at each point - sample closest pixel in template raster
  float X_coordinate_shifted_origin = ProfilePoints.X[0] - RasterTemplate.get_XMinimum();
  float Y_coordinate_shifted_origin = ProfilePoints.Y[0] - RasterTemplate.get_YMinimum();
  // Get row and column of point
  int col_point = int(round(X_coordinate_shifted_origin/Resolution));
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/Resolution));
  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows) BaselineValue_temp[0]=NoDataValue;
  else
  {
    BaselineValue_temp[0] = RasterTemplate.get_data_element(row_point,col_point);
    BaselineRows_temp[0] = row_point;
    BaselineCols_temp[0] = col_point;
  }


  for(int i = 1; i<NPtsInProfile; ++i)
  {
    // Calculate cumulative distance
    cumulative_distance += sqrt((ProfilePoints.X[i]-ProfilePoints.X[i-1])*(ProfilePoints.X[i]-ProfilePoints.X[i-1])
                                              +(ProfilePoints.Y[i]-ProfilePoints.Y[i-1])*(ProfilePoints.Y[i]-ProfilePoints.Y[i-1]));
    DistanceAlongBaseline_temp[i]=cumulative_distance;

    // Retrieve baseline value at each point - sample closest pixel in template raster
    X_coordinate_shifted_origin = ProfilePoints.X[i] - RasterTemplate.get_XMinimum();
    Y_coordinate_shifted_origin = ProfilePoints.Y[i] - RasterTemplate.get_YMinimum();

    // Get row and column of point
    col_point = int(round(X_coordinate_shifted_origin/Resolution));
    row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/Resolution));
    if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows) BaselineValue_temp[i]=NoDataValue;
    else
    {
      BaselineValue_temp[i] = RasterTemplate.get_data_element(row_point,col_point);
      BaselineRows_temp[i] = row_point;
      BaselineCols_temp[i] = col_point;
    }
  }

  DistanceAlongBaseline = DistanceAlongBaseline_temp;
  BaselineValue = BaselineValue_temp;

  // Read profile data into a LSDCloud object for querying
  vector<float> zeros(NPtsInProfile,0.0);
  LSDCloudRaster ProfileCloud(ProfilePoints, RasterTemplate);

  // For each point in array, find nearest point along the profile and calculate
  // signed distance to profile.  The convention here is that points lying on
  // left hand side of profile as you traverse from start to finish are
  // considered positive, and those on the right are negative.
  Array2D<float> DistanceToBaseline_temp(NRows,NCols,NoDataValue);
  Array2D<float> ProjectedDistanceAlongBaseline_temp(NRows,NCols,NoDataValue);
  Array2D<float> ProjectedBaselineValue_temp(NRows,NCols,NoDataValue);
  Array2D<float> BaselineValueArray_temp(NRows,NCols,NoDataValue);

  float searchPoint_x,searchPoint_y;
  vector<float> SquaredDistanceToBaseline,temp,temp2;
  vector<int> ProfilePointIndex;
  int K=2;

  // Define bounding box of swath profile
  YMin = ProfileCloud.get_YMin();
  YMax = ProfileCloud.get_YMax();
  XMin = ProfileCloud.get_XMin();
  XMax = ProfileCloud.get_XMax();
  int ColStart = int(floor(XMin/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    cout << flush << "\t\t\t row " << i+1 << "/" << NRows << "\r";
    for(int j = ColStart; j<ColEnd; ++j)
    {
      vector<float> V_p;  // the search point coordinates
      V_p.push_back(float(j)*Resolution);
      V_p.push_back((float(NRows - 1) - float(i))*Resolution);
      // Find nearest two points on profile baseline
      ProfileCloud.NearestNeighbourSearch2D(V_p[0], V_p[1], K, temp, ProfilePointIndex, SquaredDistanceToBaseline);

      if(SquaredDistanceToBaseline.size()<1) cout << "No profile points found - check input files";
      else if(SquaredDistanceToBaseline.size()<2) cout << "Only one point found - check input files";
      else
      {
        // SCENARIO 1
        //----------------------------------------------------------------------
        // Check that two points are adjacenet on baseline - simplest (and
        // usual) scenario)
        vector<float> V_a, V_b; // first, second nearest points on profile
        int V_a_index, V_b_index;
        if(abs(ProfilePointIndex[0]-ProfilePointIndex[1]) <= 1)
        {
//           cout << "SCENARIO 1" << endl;
          // If yes, then find point along cojoining line that minimises this
          // distance.
          // Note that this point is located at the point given by the vector:
          // (V_a + (V_b - V_a)*t).
          // t is calculated using the equation:
          // t = -(V_b-V_a).(V_a-V_p)/(|V_b - V_a|^2)
          V_a_index = ProfilePointIndex[0];
          V_b_index = ProfilePointIndex[1];
          V_a.push_back(ProfileCloud.get_point_x(V_a_index));
          V_a.push_back(ProfileCloud.get_point_y(V_a_index));
          V_b.push_back(ProfileCloud.get_point_x(V_b_index));
          V_b.push_back(ProfileCloud.get_point_y(V_b_index));
        }

        // SCENARIO 2
        //----------------------------------------------------------------------
        // Two nearest points are not adjacent on the profile.  This is quite
        // common on the inside of a bend on a baseline.  At this point, the
        // strategy used is to find all baseline points located at a radius that
        // is equal to that of the nearest baseline point.  Then, the strategy
        // is to search through all located baseline points, and find the
        // nearest point that neighbours one of these baseline points.  The
        // program then searches through to find the closest point on the
        // baseline between those points.
        else
        {
//           cout << "SCENARIO 2" << endl;//" " << sqrt(SquaredDistanceToBaseline[0]) << endl;
          vector<int> ClosestPointIndex;
          float shortest_distance = NoDataValue;
          float distance_to_point;
          //ProfileCloud.RadiusSearch2D(V_p[0], V_p[1], sqrt(SquaredDistanceToBaseline[0]), temp, ClosestPointIndex, temp2);
          ProfileCloud.NearestNeighbourSearch2D(V_p[0], V_p[1], 8, temp, ClosestPointIndex, temp2);
          for(int i_point_iter = 0; i_point_iter < ClosestPointIndex.size(); ++i_point_iter)
          {
            if((ClosestPointIndex[i_point_iter]!=0) && (ClosestPointIndex[i_point_iter]!=NPtsInProfile-1))
            {
              // Check step backwards
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]-1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]-1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]-1;
              }
              // Check step forwards
              V_b_test.clear();
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]+1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]+1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]+1;
              }
            }
            // conditions for edge of profile
            else if(ClosestPointIndex[i_point_iter]==0)
            {
              // Check step forwards only
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]+1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]+1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]+1;
              }
            }
            else if(ClosestPointIndex[i_point_iter]==NPtsInProfile-1)
            {
              // Check step backwards only
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]-1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]-1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]-1;
              }
            }
          }
        }
        // calculate position along baseline vector
        float t = calculate_t(V_a, V_b, V_p);
        float d;
        if(t>0 && t<1)
        // find the distance to the nearest point along the straight line
        // segment that links the two points
        {
//          cout << "test1" << endl;
          d = calculate_shortest_distance_to_line(V_a,V_b,V_p);
          if(d<ProfileHalfWidth)
          {
            ProjectedDistanceAlongBaseline_temp[i][j] = DistanceAlongBaseline[V_a_index]+(DistanceAlongBaseline[V_b_index]-DistanceAlongBaseline[V_a_index])*t;
            // determine which side of the profile the point is on using the cross
            // product
            if(V_a_index < V_b_index)
            {
              if(test_point_left(V_b,V_a,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
            else
            {
              if(test_point_left(V_a,V_b,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
            // Now get baseline value
            BaselineValueArray_temp[i][j] = BaselineValue[V_a_index]+(BaselineValue[V_b_index]-BaselineValue[V_a_index])*t;
          }
          else
          {
            DistanceToBaseline_temp[i][j] = NoDataValue;
            ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
            BaselineValueArray_temp[i][j] = NoDataValue;
          }
        }

        else if(V_a_index==0 || V_a_index==NPtsInProfile-1)
        // Avoid end points
        {
//          cout << "test2" << endl;
          DistanceToBaseline_temp[i][j] = NoDataValue;
          ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
          BaselineValueArray_temp[i][j] = NoDataValue;
        }
        else
        // Select nearest point (e.g. on outer side of bend)
        {
//          cout << "test3" << endl;
          d = calculate_distance_between_two_points(V_a,V_p);
          if(d<ProfileHalfWidth)
          {
            ProjectedDistanceAlongBaseline_temp[i][j] = DistanceAlongBaseline[V_a_index];
            BaselineValueArray_temp[i][j] = BaselineValue[V_a_index];
            // determine which side of the profile the point is on using the cross
            // product.  In this case, take the baseline vector as being formed
            // by straight line that joins the points on either side.
            vector<float> V_c(2,0.0);
            if(V_a_index < V_b_index)
            {
              V_c[0]=ProfileCloud.get_point_x(V_a_index-1);
              V_c[1]=ProfileCloud.get_point_y(V_a_index-1);
              if(test_point_left(V_b,V_c,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
            else
            {
              V_c[0]=ProfileCloud.get_point_x(V_a_index+1);
              V_c[1]=ProfileCloud.get_point_y(V_a_index+1);
              if(test_point_left(V_c,V_b,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;
            }
          }
          else
          {
            DistanceToBaseline_temp[i][j] = NoDataValue;
            ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
            BaselineValueArray_temp[i][j] = NoDataValue;
          }
        }
      }
    }
  }
  DistanceToBaselineArray = DistanceToBaseline_temp.copy();
  DistanceAlongBaselineArray = ProjectedDistanceAlongBaseline_temp.copy();
  BaselineValueArray = BaselineValueArray_temp.copy();
  BaselineRows = BaselineRows_temp;
  BaselineCols = BaselineCols_temp;

}

//------------------------------------------------------------------------------
// SWATH PROFILE GENERATION
// These routines take a swath profile template, comprising the LSDSwath object,
// and then uses this to construct either transverse (normal to profile) or
// longitudinal (parallel to profile) profiles.

// GET_TRANSVERSE_SWATH_PROFILES
// Function takes a raster and calculates transverse swath profiles, based on
// the swath template in the LSDSwath object.  Note that the input raster at
// present must have the same extent as the original template raster used to
// create the LSDSwath object.
// The function returns a vector container of the desired profiles.  The first
// and second profiles in the container are ALWAYS the mean and standard
// deviation respectively.  The following profiles contain the desired
// percentile profiles indicated in the input vector "desired_percentiles".
void LSDSwath::get_transverse_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles,
       int NormaliseToBaseline)
{
  vector<float> TransverseDistance, RasterValues;
  float Resolution = Raster.get_DataResolution();
  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if((DistanceToBaselineArray[i][j]!=NoDataValue) && (Raster.get_data_element(i,j)!=NoDataValue))
      {
        TransverseDistance.push_back(DistanceToBaselineArray[i][j]);
        // To normalise the profile values, subtract baseline value from raster value
        if (NormaliseToBaseline == 1) RasterValues.push_back(Raster.get_data_element(i,j)-BaselineValueArray[i][j]);
        // Otherwise, just get the raster value
        else RasterValues.push_back(Raster.get_data_element(i,j));
      }
    }
  }
  cout << "test" << endl;
  // Sort values if you need percentiles
  int NumberOfPercentileProfiles = desired_percentiles.size();
  if (NumberOfPercentileProfiles > 0)
  {
    vector<size_t> index_map;
    matlab_float_sort(RasterValues,RasterValues,index_map);
    matlab_float_reorder(TransverseDistance, index_map, TransverseDistance);
  }

  // Bin data
  cout << "BIN DATA" << endl;
  vector< vector<float> > binned_RasterValues;
  bin_data(RasterValues, TransverseDistance, -ProfileHalfWidth, ProfileHalfWidth, BinWidth, mid_points, binned_RasterValues);
  cout << "done" << endl;
  // Produce desired profiles from binned_RasterValues.
  int NBins = mid_points.size();
  vector<float> mean_profile_temp(NBins, NoDataValue);
  vector<float> sd_profile_temp(NBins, NoDataValue);
  vector< vector<float> > output_percentile_profiles_temp;
  vector<int> N_data(NBins,0);
  cout << "\t CALCULATING MEAN AND SD" << endl;
  for(int i = 0; i < NBins; ++i)
  {
    N_data[i] = binned_RasterValues[i].size();
    if (N_data[i] >= 1)
    {
      mean_profile_temp[i]=get_mean(binned_RasterValues[i]);
      sd_profile_temp[i]=get_standard_deviation(binned_RasterValues[i], mean_profile_temp[i]);
    }
    else
    {
      mean_profile_temp[i]=NoDataValue;
      sd_profile_temp[i]=NoDataValue;
    }
  }
  cout << "\t done" << endl;
  cout << "\t CALCULATING PERCENTILES" << endl;
  for(int j = 0; j < NumberOfPercentileProfiles; ++j)
  {
    vector<float> percentile(NBins,NoDataValue);
    for(int i = 0; i < NBins; ++i)
    {
      if (N_data[i] >= 1) percentile[i] = get_percentile(binned_RasterValues[i], desired_percentiles[j]);
      else percentile[i] = NoDataValue;
    }
    output_percentile_profiles_temp.push_back(percentile);
  }
  cout << "\t done" << endl;
  // Load profiles into export vectors
  mean_profile = mean_profile_temp;
  sd_profile = sd_profile_temp;
  output_percentile_profiles = output_percentile_profiles_temp;
}

// GET_LONGITUDINAL_SWATH_PROFILES
// Function takes a raster and calculates longitudinal swath profiles, based on
// the swath template in the LSDSwath object.  Note that the input raster at
// present must have the same extent as the original template raster used to
// create the LSDSwath object.
// The function returns a vector container of the desired profiles.  The first
// and second profiles in the container are ALWAYS the mean and standard
// deviation respectively.  The following profiles contain the desired
// percentile profiles indicated in the input vector "desired_percentiles".
//
// SMM notes: BinWidth is the distance along the swath where you want the points.
void LSDSwath::get_longitudinal_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles,
       int NormaliseToBaseline)
{
  vector<float> LongitudinalDistance, RasterValues;
  float Resolution = Raster.get_DataResolution();
  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if((DistanceAlongBaselineArray[i][j]!=NoDataValue) && (Raster.get_data_element(i,j)!=NoDataValue))
      {
        LongitudinalDistance.push_back(DistanceAlongBaselineArray[i][j]);
        // To normalise the profile values, subtract baseline value from raster value
        if (NormaliseToBaseline == 1) RasterValues.push_back(Raster.get_data_element(i,j)-BaselineValueArray[i][j]);
        // Otherwise, just get the raster value
        else RasterValues.push_back(Raster.get_data_element(i,j));
      }
    }
  }
  // Sort values if you need percentiles
  int NumberOfPercentileProfiles = desired_percentiles.size();
  if (NumberOfPercentileProfiles > 0)
  {
    vector<size_t> index_map;
    matlab_float_sort(RasterValues,RasterValues,index_map);
    matlab_float_reorder(LongitudinalDistance, index_map, LongitudinalDistance);
  }

  // Bin data
  vector< vector<float> > binned_RasterValues;
  float start_point = DistanceAlongBaseline[0];
  float end_point = DistanceAlongBaseline[NPtsInProfile-1];
  bin_data(RasterValues, LongitudinalDistance, start_point, end_point, BinWidth, mid_points, binned_RasterValues);
  // Produce desired profiles from binned_RasterValues.
  int NBins = mid_points.size();
  vector<float> mean_profile_temp(NBins, NoDataValue);
  vector<float> sd_profile_temp(NBins, NoDataValue);
  vector< vector<float> > output_percentile_profiles_temp;
  vector<int> N_data(NBins,0);
  cout << "\t CALCULATING MEAN AND SD" << endl;
  for(int i = 0; i < NBins; ++i)
  {
    N_data[i] = binned_RasterValues[i].size();
    if (N_data[i] >= 1)
    {
      mean_profile_temp[i]=get_mean(binned_RasterValues[i]);
      sd_profile_temp[i]=get_standard_deviation(binned_RasterValues[i], mean_profile_temp[i]);
    }
    else
    {
      mean_profile_temp[i]=NoDataValue;
      sd_profile_temp[i]=NoDataValue;
    }
  }
  cout << "\t done" << endl;
  cout << "\t CALCULATING PERCENTILES" << endl;
  for(int j = 0; j < NumberOfPercentileProfiles; ++j)
  {
    vector<float> percentile(NBins,NoDataValue);
    for(int i = 0; i < NBins; ++i)
    {
      if (N_data[i] >= 1) percentile[i] = get_percentile(binned_RasterValues[i], desired_percentiles[j]);
      else percentile[i] = NoDataValue;
    }
    output_percentile_profiles_temp.push_back(percentile);
  }
  cout << "\t done" << endl;
  // Load profiles into export vectors
  mean_profile = mean_profile_temp;
  sd_profile = sd_profile_temp;
  output_percentile_profiles = output_percentile_profiles_temp;
}

//------------------------------------------------------------------------------
// CREATE RASTER FROM SWATH TOOL
// This function creates a raster of values used by the Swath Profiler
// FJC
// 16/10/15
//
//------------------------------------------------------------------------------
LSDRaster LSDSwath::get_raster_from_swath_profile(LSDRaster& Raster, int NormaliseToBaseline)
{
  vector<float> TransverseDistance, RasterValues;
  Array2D<float> RasterValues_temp(NRows,NCols,NoDataValue);
  float Resolution = Raster.get_DataResolution();
  float XMinimum = Raster.get_XMinimum();
  float YMinimum = Raster.get_YMinimum();
  map<string,string> GeoReferencingStrings = Raster.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if((DistanceToBaselineArray[i][j]!=NoDataValue) && (Raster.get_data_element(i,j)!=NoDataValue) && BaselineValueArray[i][j]!=NoDataValue)
      {
        TransverseDistance.push_back(DistanceToBaselineArray[i][j]);
        // To normalise the profile values, subtract baseline value from raster value
        if (NormaliseToBaseline == 1) RasterValues_temp[i][j] = Raster.get_data_element(i,j)-BaselineValueArray[i][j];
        // Otherwise, just get the raster value
        else RasterValues_temp[i][j] = Raster.get_data_element(i,j);
      }
    }
  }
  LSDRaster SwathRaster(NRows, NCols, XMinimum, YMinimum, Resolution,
                   NoDataValue, RasterValues_temp, GeoReferencingStrings);

  return SwathRaster;
}

//------------------------------------------------------------------------------
// This function uses the swath profile to "fill in" the baseline channel raster value
// with the average value of the pixels along the transverse swath profile
// FJC
// 16/01/17
//
//------------------------------------------------------------------------------
LSDRaster LSDSwath::fill_in_channels_swath(LSDRaster& Raster)
{
	Array2D<float> NewRasterValues(NRows,NCols,NoDataValue);

	float Resolution = Raster.get_DataResolution();
  float XMinimum = Raster.get_XMinimum();
  float YMinimum = Raster.get_YMinimum();
	Array2D<float> RasterValues_temp = Raster.get_RasterData();
	map<string,string> GeoReferencingStrings = Raster.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

	// Go down the baseline and collect values
	vector<float> MeanValues;
	vector<float> ThisPixelValues;

	for (int i = 0; i < int(DistanceAlongBaseline.size()); i++)
	{
		for (int row=RowStart; row<RowEnd; row++)
		{
			for (int col=ColStart; col<ColEnd; col++)
			{
				if (BaselineValueArray[row][col] == BaselineValue[i])
				{
					//this row and col corresponds to this point on the baseline, push back the raster value
					ThisPixelValues.push_back(RasterValues_temp[row][col]);
				}
			}
		}
		// push the vector back to the vector of vectors
		float mean_value = get_mean(ThisPixelValues);
		MeanValues.push_back(mean_value);
		ThisPixelValues.clear();
	}

	// assign the new raster values for the baseline
	for (int i = 0; i < int(DistanceAlongBaseline.size()); i++)
	{
		// get the row and col of the baseline points
		float NewRasterValue = MeanValues[i];
		RasterValues_temp[BaselineRows[i]][BaselineCols[i]] = NewRasterValue;
	}

	// create raster from the array
	NewRasterValues = RasterValues_temp;
	LSDRaster FilledRaster(NRows,NCols,XMinimum,YMinimum,Resolution,NoDataValue,
												 NewRasterValues,GeoReferencingStrings);

	return FilledRaster;

}
//------------------------------------------------------------------------------
// This function takes in a connected components raster and returns the average
// value of another chosen raster and distance along the baseline of each
// connected components patch. The user can choose whether to normalise the
// second raster to the baseline value.
// vector of vectors has the format:
// 0 = patch ID
// 1 = mean distance along the baseline
// 2 = mean raster value for the patch id
// 3 = standard deviation of raster value for patch id
// 4 = standard error of raster value for patch id
// FJC
// 23/01/17
//
//------------------------------------------------------------------------------
vector <vector <float> > LSDSwath::get_connected_components_along_swath(LSDIndexRaster& ConnectedComponents, LSDRaster& RasterTemplate, int NormaliseToBaseline)
{
  vector <vector <float> > MasterVector;
  vector<float> RasterValues;
  vector<float> DistAlongBaseline;
  vector<float> RasterStDevs;
  vector<float> RasterStErrs;

  float Resolution = RasterTemplate.get_DataResolution();
	Array2D<float> RasterValues_temp = RasterTemplate.get_RasterData();
	map<string,string> GeoReferencingStrings = RasterTemplate.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  // get all the unique CC values
  Array2D<int> PatchIDs = ConnectedComponents.get_RasterData();
  vector<int> Unique_Patches = Unique(PatchIDs, NoDataValue);

  for (int i = 0; i < int(Unique_Patches.size()); i++)
  {
    vector<float> raster_values;
    vector<float> DistAlongBaseline_temp;
    int n_observations = 0;
    int target_CC = Unique_Patches[i];
    for (int row=RowStart; row<RowEnd; row++)
    {
      for (int col=ColStart; col<ColEnd; col++)
      {
        // get the connected components value
        int this_CC = ConnectedComponents.get_data_element(row,col);
        if (this_CC == target_CC)
        {
          //push back value and number_observations
          if (NormaliseToBaseline == 1)
          {
            if (RasterValues_temp[row][col] != NoDataValue && BaselineValueArray[row][col] != NoDataValue)
            {
            //normalise to the baseline
              float this_value = RasterValues_temp[row][col] - BaselineValueArray[row][col];
              raster_values.push_back(this_value);
              DistAlongBaseline_temp.push_back(DistanceAlongBaselineArray[row][col]);
            }
          }
          else
          {
            if (RasterValues_temp[row][col] != NoDataValue && BaselineValueArray[row][col] != NoDataValue)
            {
              //just push back
              raster_values.push_back(RasterValues_temp[row][col]);
              DistAlongBaseline_temp.push_back(DistanceAlongBaselineArray[row][col]);
            }
          }
        }
      }
    }
    float mean_value = get_mean(raster_values);
    float stdev = get_standard_deviation(raster_values, mean_value);
    float sterr = get_standard_error(raster_values, stdev);
    float mean_dist = get_mean(DistAlongBaseline_temp);
    RasterValues.push_back(mean_value);
    RasterStDevs.push_back(stdev);
    RasterStErrs.push_back(sterr);
    DistAlongBaseline.push_back(mean_dist);
  }

  vector<float> Unique_Patches_float(Unique_Patches.begin(), Unique_Patches.end());

  // store in the MasterVector
  MasterVector.push_back(Unique_Patches_float);
  MasterVector.push_back(DistAlongBaseline);
  MasterVector.push_back(RasterValues);
  MasterVector.push_back(RasterStDevs);
  MasterVector.push_back(RasterStErrs);

  return MasterVector;
}

//------------------------------------------------------------------------------
// This function takes in a connected components raster
// and returns an array with the distance along the baseline of each
// point in the connected components raster
// The raster_ID argument allows you to name the csv with a meaningful header
// (e.g. elevation)
// FJC
// 28/09/17
//
//------------------------------------------------------------------------------
Array2D<float> LSDSwath::get_BaselineDist_ConnectedComponents(LSDIndexRaster& ConnectedComponents)
{
  // declare baseline dist array
  Array2D<float> BaselineDistance(NRows,NCols,NoDataValue);

  // get some raster info
  float Resolution = ConnectedComponents.get_DataResolution();
	map<string,string> GeoReferencingStrings = ConnectedComponents.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  // get the array of the connected components
  Array2D<int> CC_array = ConnectedComponents.get_RasterData();

  for (int row=RowStart; row<RowEnd; row++)
  {
    for (int col=ColStart; col<ColEnd; col++)
    {
      if (CC_array[row][col] != NoDataValue)
      {
        // get the baseline value
        BaselineDistance[row][col] = DistanceAlongBaselineArray[row][col];
      }
    }
  }

  return BaselineDistance;

}

//------------------------------------------------------------------------------
// This function takes in a connected components raster
// and returns an array with the euclidian distance to the nearest point on
// the baseline for each pixel of the connected components raster
// FJC
// 12/10/17
//
//------------------------------------------------------------------------------
Array2D<float> LSDSwath::get_DistanceToBaseline_ConnectedComponents(LSDIndexRaster& ConnectedComponents)
{
  // declare baseline dist array
  Array2D<float> BaselineDistance(NRows,NCols,NoDataValue);

  // get some raster info
  float Resolution = ConnectedComponents.get_DataResolution();
	map<string,string> GeoReferencingStrings = ConnectedComponents.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  // get the array of the connected components
  Array2D<int> CC_array = ConnectedComponents.get_RasterData();

  for (int row=RowStart; row<RowEnd; row++)
  {
    for (int col=ColStart; col<ColEnd; col++)
    {
      if (CC_array[row][col] != NoDataValue)
      {
        // get the baseline value
        BaselineDistance[row][col] = DistanceToBaselineArray[row][col];
      }
    }
  }

  return BaselineDistance;

}
//------------------------------------------------------------------------------
// This function takes in a raster and returns the mean, min and max values of
// the raster at each point along the swath
// vector of vectors has the format:
// 0 = distance along swath
// 1 = mean value along swath
// 2 = min value along swath
// 3 = max value along swath
// 4 = first quartile (added by BG - > 15/11/2018)
// 5 = third quartile  (added by BG - > 15/11/2018)
// FJC
// 15/02/17
//
//------------------------------------------------------------------------------
vector <vector <float> > LSDSwath::get_RasterValues_along_swath(LSDRaster& RasterTemplate, int NormaliseToBaseline)
{
  // set up vectors
  vector <vector <float> > MasterVector;
  vector<float> DistAlongBaseline;
  vector<float> MeanRasterValues;
  vector<float> MinRasterValues;
  vector<float> MaxRasterValues;
  vector<float> FQRasterValues;
  vector<float> TQRasterValues;

  float Resolution = RasterTemplate.get_DataResolution();
	Array2D<float> RasterValues_temp = RasterTemplate.get_RasterData();
	map<string,string> GeoReferencingStrings = RasterTemplate.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  // loop through each point on the baseline and find all the pixels that are closest to that point
  for (int i = 0; i < int(BaselineValue.size()); i++)
  {
    float this_baselinevalue = BaselineValue[i];
    vector<float> raster_values;
    for (int row=RowStart; row<RowEnd; row++)
    {
      for (int col=ColStart; col<ColEnd; col++)
      {
        //cout << "DistanceAlongBaselineArray: " << DistanceAlongBaselineArray[row][col] << endl;
        // get the points that correspond to this baseline distance
        if (BaselineValueArray[row][col] == this_baselinevalue)
        {
          if (RasterValues_temp[row][col] != NoDataValue)
          {
            if (NormaliseToBaseline == 1)
            {
              // normalise to the baseline values
              float baseline_value = BaselineValue[i];
              float this_value = RasterValues_temp[row][col] - baseline_value;
              raster_values.push_back(this_value);
            }
            else
            {
              // push back corresponding values to a raster
              raster_values.push_back(RasterValues_temp[row][col]);
            }
          }
        }
      }
    }
    if (int(raster_values.size()) > 1)
    {
      float mean_value = get_mean_ignore_ndv(raster_values, NoDataValue);
      float min_value = Get_Minimum(raster_values, NoDataValue);
      float max_value = Get_Maximum(raster_values, NoDataValue);
      float FQ = get_percentile(raster_values, 25);
      float TQ = get_percentile(raster_values, 75);

      DistAlongBaseline.push_back(DistanceAlongBaseline[i]);
      MeanRasterValues.push_back(mean_value);
      MinRasterValues.push_back(min_value);
      MaxRasterValues.push_back(max_value);
      FQRasterValues.push_back(FQ);
      TQRasterValues.push_back(TQ);
      //cout << "Distance: " << DistanceAlongBaseline[i] << " n_raster values: " << raster_values.size() << endl;
    }
    else
    {
      DistAlongBaseline.push_back(DistanceAlongBaseline[i]);
      MeanRasterValues.push_back(NoDataValue);
      MinRasterValues.push_back(NoDataValue);
      MaxRasterValues.push_back(NoDataValue);
      FQRasterValues.push_back(NoDataValue);
      TQRasterValues.push_back(NoDataValue);
    }
  }

  // store in the MasterVector
  MasterVector.push_back(DistAlongBaseline);
  MasterVector.push_back(MeanRasterValues);
  MasterVector.push_back(MinRasterValues);
  MasterVector.push_back(MaxRasterValues);
  MasterVector.push_back(FQRasterValues);
  MasterVector.push_back(TQRasterValues);

  return MasterVector;
}

void LSDSwath::write_RasterValues_along_swath_to_csv(LSDRaster& RasterTemplate, int NormaliseToBaseline, string csv_filename)
{
  vector<vector <float> > MasterVector = get_RasterValues_along_swath(RasterTemplate, NormaliseToBaseline);

  // setup the output csv
  ofstream output_file;
  output_file.open(csv_filename.c_str());
  output_file.precision(8);
  if (!output_file)
  {
     cout << "\n Error opening output csv file. Please check your filename";
     exit(1);
  }
  cout << "Opened the csv" << endl;

  output_file << "DistAlongBaseline,X,Y,latitude,longitude,mean,min,max" << endl;

  double x_loc, y_loc;
  double latitude, longitude;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  int n_points= BaselineValue.size();

  for (int i = 0; i < n_points; i++)
  {
    float this_dist = DistanceAlongBaseline[i];
    if (MasterVector[1][i] != NoDataValue)
    {
      // get the latitude and longitude of the point
      RasterTemplate.get_x_and_y_locations(BaselineRows[i], BaselineCols[i], x_loc, y_loc);
      //cout << "Row: " << row << " Col: " << col << " X: " << x_loc << " Y: " << y_loc << endl;
      RasterTemplate.get_lat_and_long_locations(BaselineRows[i], BaselineCols[i], latitude, longitude, Converter);
      output_file << MasterVector[0][i] << "," << x_loc << "," << y_loc << "," << latitude << "," << longitude << "," << MasterVector[1][i] << "," << MasterVector[2][i] << "," << MasterVector[3][i] << endl;
    }
  }

  output_file.close();

}

//---------------------------------------------------------//
// Take in a swath profile, and get the width of the
// index raster values at every point along the baseline.
// Finds the max and min distance from baseline at each
// point: width = max - min
// Used for getting valley width.
// FJC 21/11/17
//---------------------------------------------------------//
vector<float> LSDSwath::get_widths_along_swath(LSDIndexRaster& RasterForAnalysis)
{
	vector<float> Widths;
	Array2D<int> MaskArray = RasterForAnalysis.get_RasterData();

  float Resolution = RasterForAnalysis.get_DataResolution();
	Array2D<int> RasterValues_temp = RasterForAnalysis.get_RasterData();
	map<string,string> GeoReferencingStrings = RasterForAnalysis.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  // loop through the baseline and find all points in the swath that are closest to this point
  for (int i = 0; i < NPtsInProfile; i++)
  {
    float max_dist = 0;
    float min_dist = 10000000000000000;
    for (int row=RowStart; row<RowEnd; row++)
    {
      for (int col=ColStart; col<ColEnd; col++)
      {
        if (MaskArray[row][col] != NoDataValue)
        {
          float this_dist = DistanceAlongBaselineArray[row][col];
          if (this_dist == DistanceAlongBaseline[i])
          {
            // now check if this is the max or the min for this point.
            if (this_dist > max_dist) { max_dist = this_dist; }
            if (this_dist < min_dist) { min_dist = this_dist; }
          }
        }
      }
    }
    // now get the width at this point (max - min)
    float width = max_dist - min_dist;
    Widths.push_back(width);
  }

  return Widths;
}

//------------------------------------------------------------------------------
// WRITE PROFILES TO FILE
// These routines take a swath profile template, comprising the LSDSwath object,
// and then uses this to construct either transverse (normal to profile) or
// longitudinal (parallel to profile) profiles.
void LSDSwath::write_transverse_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
      string prefix, int NormaliseToBaseline)
{
  string profile_extension = "_trans_profile.txt";
  vector<float> mid_points;
  vector<float> mean_profile;
  vector<float> sd_profile;
  vector< vector<float> > output_percentile_profiles;
  get_transverse_swath_profile(Raster, desired_percentiles, BinWidth, mid_points, mean_profile, sd_profile, output_percentile_profiles, NormaliseToBaseline);
  // Print profiles to file
  string filename = prefix + profile_extension;
  cout << "\t printing profiles to " << filename << endl;
  ofstream ofs;
  ofs.open(filename.c_str());

  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
		exit(EXIT_FAILURE);
  }

  ofs << "Midpoint Mean SD ";
  for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << desired_percentiles[i_perc] << " ";
  ofs << "\n";
  for(int i = 0; i < mid_points.size(); ++i)
  {
    ofs << mid_points[i] << " " << mean_profile[i] << " " << sd_profile[i] << " " ;
    for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << output_percentile_profiles[i_perc][i] << " ";
    ofs << "\n";
  }

  ofs.close();
}
void LSDSwath::write_longitudinal_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
      string prefix, int NormaliseToBaseline)
{
  string profile_extension = "_long_profile.txt";
  vector<float> mid_points;
  vector<float> mean_profile;
  vector<float> sd_profile;
  vector< vector<float> > output_percentile_profiles;
  get_longitudinal_swath_profile(Raster, desired_percentiles, BinWidth, mid_points, mean_profile, sd_profile, output_percentile_profiles, NormaliseToBaseline);
  // Print profiles to file
  string filename = prefix + profile_extension;
  cout << "\t printing profiles to " << filename << endl;
  ofstream ofs;
  ofs.open(filename.c_str());

  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
		exit(EXIT_FAILURE);
  }

  ofs << "Midpoint Mean SD ";
  for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << desired_percentiles[i_perc] << " ";
  ofs << "\n";
  for(int i = 0; i < mid_points.size(); ++i)
  {
    ofs << mid_points[i] << " " << mean_profile[i] << " " << sd_profile[i] << " " ;
    for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << output_percentile_profiles[i_perc][i] << " ";
    ofs << "\n";
  }

  ofs.close();
}

//---------------------------------------------------------------------------//
// Function to print the baseline of the swath profile to a csv. prints the
// distance along swath and the elevation of each point.
// FJC 12/10/17
//---------------------------------------------------------------------------//
void LSDSwath::print_baseline_to_csv(LSDRaster& ElevationRaster, string csv_filename, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet)
{
  Array2D<float> ElevationArray = ElevationRaster.get_RasterData();

  // setup the output csv
  ofstream output_file;
  output_file.open(csv_filename.c_str());
  output_file.precision(8);
  if (!output_file)
  {
     cout << "\n Error opening output csv file. Please check your filename";
     exit(1);
  }
  cout << "Opened the csv" << endl;

  output_file << "DistAlongBaseline,Elevation,X,Y,latitude,longitude,row,col,node,DistFromOutlet" << endl;
  double x_loc, y_loc;
  double latitude, longitude;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  cout << "N baseline rows: " << BaselineRows.size() << " N baseline cols: " << BaselineCols.size() << endl;
  // loop through and get the values for the csv
  for (int i =0; i < NPtsInProfile; i++)
  {
    // get the latitude and longitude of the point
    ElevationRaster.get_x_and_y_locations(BaselineRows[i], BaselineCols[i], x_loc, y_loc);
    //cout << "Row: " << row << " Col: " << col << " X: " << x_loc << " Y: " << y_loc << endl;
    ElevationRaster.get_lat_and_long_locations(BaselineRows[i], BaselineCols[i], latitude, longitude, Converter);
    // get the node
    int this_node = FlowInfo.retrieve_node_from_row_and_column(BaselineRows[i],BaselineCols[i]);
    float DistFromOutlet = DistanceFromOutlet.get_data_element(BaselineRows[i], BaselineCols[i]);
    output_file << DistanceAlongBaseline[i] << "," << BaselineValue[i] << "," << x_loc << "," << y_loc << "," << latitude << "," << longitude << "," << BaselineRows[i] << "," << BaselineCols[i] << "," << this_node << "," << DistFromOutlet << endl;
  }
  output_file.close();
}


// Function to read in a CSV file and snap the points to the baseline.
// Have to do this for the stupid terrace profiles
// FJC 20/03/18
// map< vector<float>> LSDSwath::snap_points_to_baseline(LSDSpatialCSVReader csv_data)
// {
//   // for each row, snap to the nearest baseline.
//
// }



void get_locations_of_points_along_baseline(vector<float> distance_along_baseline)
{
  int N_distances = int(distance_along_baseline.size());
  float this_distance;
  for (int i = 0; i< N_distances; i++)
  {
    this_distance = distance_along_baseline[i];

    // now loop through baseline finding rows and columns

  }
}



LSDRaster LSDSwath::get_raster_DistanceToBaselineArray(LSDRaster& Raster)
{
  vector<float> TransverseDistance, RasterValues;
  Array2D<float> RasterValues_temp(NRows,NCols,NoDataValue);
  float Resolution = Raster.get_DataResolution();
  float XMinimum = Raster.get_XMinimum();
  float YMinimum = Raster.get_YMinimum();
  map<string,string> GeoReferencingStrings = Raster.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if(DistanceToBaselineArray[i][j]!=NoDataValue)
      {
        RasterValues_temp[i][j] = DistanceToBaselineArray[i][j];
      }
    }
  }
  LSDRaster SwathRaster(NRows, NCols, XMinimum, YMinimum, Resolution,
                   NoDataValue, RasterValues_temp, GeoReferencingStrings);

  return SwathRaster;
}


LSDRaster LSDSwath::get_raster_DistanceAlongBaselineArray(LSDRaster& Raster)
{
  vector<float> TransverseDistance, RasterValues;
  Array2D<float> RasterValues_temp(NRows,NCols,NoDataValue);
  float Resolution = Raster.get_DataResolution();
  float XMinimum = Raster.get_XMinimum();
  float YMinimum = Raster.get_YMinimum();
  map<string,string> GeoReferencingStrings = Raster.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if(DistanceAlongBaselineArray[i][j]!=NoDataValue)
      {
        RasterValues_temp[i][j] = DistanceAlongBaselineArray[i][j];
      }
    }
  }
  LSDRaster SwathRaster(NRows, NCols, XMinimum, YMinimum, Resolution,
                   NoDataValue, RasterValues_temp, GeoReferencingStrings);

  return SwathRaster;
}

LSDRaster LSDSwath::get_raster_BaselineValueArray(LSDRaster& Raster)
{
  vector<float> TransverseDistance, RasterValues;
  Array2D<float> RasterValues_temp(NRows,NCols,NoDataValue);
  float Resolution = Raster.get_DataResolution();
  float XMinimum = Raster.get_XMinimum();
  float YMinimum = Raster.get_YMinimum();
  map<string,string> GeoReferencingStrings = Raster.get_GeoReferencingStrings();

  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;

  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if(BaselineValueArray[i][j]!=NoDataValue)
      {
        RasterValues_temp[i][j] = BaselineValueArray[i][j];
      }
    }
  }
  LSDRaster SwathRaster(NRows, NCols, XMinimum, YMinimum, Resolution,
                   NoDataValue, RasterValues_temp, GeoReferencingStrings);

  return SwathRaster;
}

#endif
