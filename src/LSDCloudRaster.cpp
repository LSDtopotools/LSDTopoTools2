//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDCloudRaster.cpp
//------------------------------------------------------------------------------
// This code is based on LSDCloudBase but does not require the liblas software.
// It cannot read in las files but instead uses the cloud functions interfacing
// with LSDRaster objects
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
// Fiona J. Clubb, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 0.0.1		28/03/17
// Prerequisite software packages: TNT, PCL
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

// LSDTopotools objects
#include "LSDRaster.hpp"
#include "LSDCloudRaster.hpp"
#include "LSDShapeTools.hpp"

// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>

// TNT
#include "TNT/tnt.h"

using namespace std;
using namespace TNT;

#ifndef LSDCloudRaster_CPP
#define LSDCloudRaster_CPP

//------------------------------------------------------------------------------
// OPERATOR FUNCTIONS
// Stuff to create a point cloud, and other useful stuff will go here
LSDCloudRaster& LSDCloudRaster::operator=(const LSDCloudRaster& rhs)
{
  if (&rhs != this)
  {
    create(rhs.get_NPts(),rhs.get_XMin(),rhs.get_XMax(),rhs.get_YMin(),rhs.get_YMax(),rhs.get_XOffset(),rhs.get_YOffset(),
           rhs.get_theCloud(),rhs.get_OctreeResolution());
  }
  return *this;
}

// the create function. This is default and throws an error
void LSDCloudRaster::create()
{
	cout << "LSDCloudRaster line 63 Warning you have an empty LSDCloudRaster" << endl;
	exit(EXIT_FAILURE);
}

void LSDCloudRaster::create(int n_pts, float x_min, float x_max, float y_min, float y_max,
              float x_offset, float y_offset, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, const double resolution)
{
  NPts = n_pts;
  pcl::PointCloud<pcl::PointXYZI>::Ptr theCloud;
  if (theCloud->points.size() != NPts)
	{
		cout << "LSDCloudBase line 75 number of points in point cloud is not the same as stated in NPts! Return to Cloud Base!" << endl;
		exit(EXIT_FAILURE);
	}
  // bounding box of point cloud
  XMin = x_min;
  YMin = y_min;
  XMax = x_max;
  YMax = y_max;
  // origin shift
  XOffset = x_offset;
  YOffset = y_offset;

  // Octree spatial structure
  OctreeResolution = resolution;
}

// This creates an LSDCloud from a LSDRaster
void LSDCloudRaster::create(LSDRaster& raster)
{
 	OctreeResolution = raster.get_DataResolution();
	raster_to_cloud(raster);
}
// Read in data from a PointData object (see LSDShapeTools module)
void LSDCloudRaster::create(PointData& point_data, LSDRaster& raster)
{
 	OctreeResolution = raster.get_DataResolution();
	PointData_to_cloud(point_data,raster);
}

// Read in data from an LSDRaster
void LSDCloudRaster::raster_to_cloud(LSDRaster& raster)
{
  XOffset = raster.get_XMinimum();
  YOffset = raster.get_YMinimum();
  NPts = 0;
  vector<float> x_coordinates, y_coordinates, zeta_values;
  for(int i = 0; i < raster.get_NRows(); ++i)
  {
    for(int j = 0; j < raster.get_NCols(); ++j)
    {
      if(raster.get_data_element(i,j) != raster.get_NoDataValue())
      {
        zeta_values.push_back(raster.get_data_element(i,j));
        x_coordinates.push_back(float(j));
        y_coordinates.push_back((float(raster.get_NRows() - 1) - float(i)));
        ++NPts;
      }
    }
  }

  // Generate pointcloud data structure
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZI>);
  theCloud = cloud;

  theCloud->width = zeta_values.size();
  theCloud->height = 1;
  theCloud->points.resize (theCloud->width * theCloud->height);

  // Loop through the .las file, reading the points into the new data structure.
  for (size_t i = 0; i < theCloud->points.size(); ++i)
  {
    theCloud->points[i].x = x_coordinates[i];
    theCloud->points[i].y = y_coordinates[i];
    theCloud->points[i].z = 0;
    theCloud->points[i].intensity = zeta_values[i];
  }
  cloudOctree = (new pcl::octree::OctreePointCloudSearch<pcl::PointXYZI>(OctreeResolution) );
  cloudOctree->setInputCloud(theCloud);
  cloudOctree->addPointsFromInputCloud();
}
// Read in data from a PointData object (see LSDShapeTools module)
void LSDCloudRaster::PointData_to_cloud(PointData& point_data, LSDRaster& raster)
{
  XOffset = raster.get_XMinimum();
  YOffset = raster.get_YMinimum();
  vector<float> x_coordinates, y_coordinates, zeta_values;
  NPts = 0;
  x_coordinates.push_back(point_data.X[0]-XOffset);
  y_coordinates.push_back(point_data.Y[0]-YOffset);
  ++NPts;
  XMin = x_coordinates[0];
  YMin = y_coordinates[0];
  XMax = x_coordinates[0];
  YMax = y_coordinates[0];
  for(int i = 1; i < point_data.X.size(); ++i)
  {
    x_coordinates.push_back(point_data.X[i]-XOffset);
    y_coordinates.push_back(point_data.Y[i]-YOffset);
    ++NPts;
    if(x_coordinates[i] < XMin) XMin = x_coordinates[i];
    if(x_coordinates[i] > XMax) XMax = x_coordinates[i];
    if(y_coordinates[i] < YMin) YMin = y_coordinates[i];
    if(y_coordinates[i] > YMax) YMax = y_coordinates[i];
  }
  // Generate pointcloud data structure
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZI>);
  theCloud = cloud;

  theCloud->width = y_coordinates.size();
  theCloud->height = 1;
  theCloud->points.resize (theCloud->width * theCloud->height);

  // Loop through the .las file, reading the points into the new data structure.
  for (size_t i = 0; i < theCloud->points.size(); ++i)
  {
    theCloud->points[i].x = x_coordinates[i];
    theCloud->points[i].y = y_coordinates[i];
    theCloud->points[i].z = 0;
    theCloud->points[i].intensity = 0;
  }
  cloudOctree = (new pcl::octree::OctreePointCloudSearch<pcl::PointXYZI>(OctreeResolution) );
  cloudOctree->setInputCloud(theCloud);
  cloudOctree->addPointsFromInputCloud();
}

//------------------------------------------------------------------------------
// SEARCHING TOOLS
// a set of tools to conduct search operations on the structured point cloud
// hosted in the LSDCloud object.

// Conducts a focal search (in the xy plane), to find all points within defined
// search radius
void LSDCloudRaster::RadiusSearch2D(float searchPoint_x, float searchPoint_y, float searchRadius,
                              vector<float>& pointValues, vector<int>& pointIndex, vector<float>& pointSquaredDistance)
{
  pcl::PointXYZI searchPoint;
  searchPoint.x = searchPoint_x;
  searchPoint.y = searchPoint_y;
  searchPoint.z = 0;
  searchPoint.intensity = 0;

  if (cloudOctree->radiusSearch (searchPoint, searchRadius, pointIndex, pointSquaredDistance) > 0)
  {
    int NumberOfPoints = pointIndex.size();
    for (int i = 0; i < NumberOfPoints; ++i)
    {
      pointValues.push_back( theCloud->points[ pointIndex[i] ].intensity );
    }
  }
  else
  {
    cout << "No returns in search radius" << endl;
  }
}
// Conducts a focal search (in the xy plane), to find the nearest K points
void LSDCloudRaster::NearestNeighbourSearch2D(float searchPoint_x, float searchPoint_y, int K,
                              vector<float>& pointValues, vector<int>& pointIndex, vector<float>& pointSquaredDistance)
{
  pcl::PointXYZI searchPoint;
  searchPoint.x = searchPoint_x;
  searchPoint.y = searchPoint_y;
  searchPoint.z = 0;
  searchPoint.intensity = 0;

  vector<int> P_i_pointIdxNKNSearch;
  vector<float> P_i_pointNKNSquaredDistance;
  if (cloudOctree->nearestKSearch (searchPoint, K, pointIndex, pointSquaredDistance) > 0)
  {
    int NumberOfPoints = pointIndex.size();
    for (int i = 0; i < NumberOfPoints; ++i)
    {
      pointValues.push_back( theCloud->points[ pointIndex[i] ].intensity );
    }
  }
  else
  {
    cout << "No points to select" << endl;
  }
}


#endif
