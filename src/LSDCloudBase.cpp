//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDCloudBase.cpp
//------------------------------------------------------------------------------
// This code houses the LSDCloud object, and associated functions, designed to
// analyse 3D pointcloud data, such as airborne LiDAR, and interface with the
// raster based LSDTopotools.  Currently reads .las files, but this could be
// expanded in the future to include other input types.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
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
#include "LSDCloudBase.hpp" 
#include "LSDShapeTools.hpp"

// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
#include <liblas/liblas.hpp>

// TNT
#include "TNT/tnt.h"

using namespace std;
using namespace TNT;

#ifndef LSDCloudBase_CPP
#define LSDCloudBase_CPP

//------------------------------------------------------------------------------
// OPERATOR FUNCTIONS
// Stuff to create a point cloud, and other useful stuff will go here
LSDCloud& LSDCloud::operator=(const LSDCloud& rhs)
{
  if (&rhs != this)
  {
    create(rhs.get_NPts(),rhs.get_XMin(),rhs.get_XMax(),rhs.get_YMin(),rhs.get_YMax(),rhs.get_XOffset(),rhs.get_YOffset(),
           rhs.get_theCloud(),rhs.get_OctreeResolution());
  }
  return *this;
}

// the create function. This is default and throws an error
void LSDCloud::create()
{
	cout << "LSDCloud line 64 Warning you have an empty LSDCloud; return to CloudBase!" << endl;
	exit(EXIT_FAILURE);
}

void LSDCloud::create(int n_pts, float x_min, float x_max, float y_min, float y_max,
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

// this creates a point cloud from an input .las file.  There are presently 3
// filter options: 0 -> select all points; 1 -> select vegetation points only;
// 2 -> select ground points only.
void LSDCloud::create(string las_file, LSDRaster& Raster, int filter)
{
	OctreeResolution = Raster.get_DataResolution();
	bool SetOffsets = FALSE;
	YOffset = Raster.get_YMinimum();
  XOffset = Raster.get_XMinimum();
  if(filter == 0)
  {
    read_las_header(las_file, SetOffsets);
    read_las_data(las_file);
  }
  else if(filter == 1)        
  {
    read_las_header(las_file, SetOffsets);
    read_las_data_vegetation(las_file);
  }
  else if(filter == 2)
  {
    read_las_header(las_file, SetOffsets);
    read_las_data_ground(las_file);
  }
  else
  {
  	cout << "LSDCloud line 78 Warning you have an unsupported filter-type; return to CloudBase!" << endl;
  	exit(EXIT_FAILURE);
  }
  //read_raster(filename,filter);
}

void LSDCloud::create(string las_file,int filter,const double octree_resolution)
{
	OctreeResolution = octree_resolution;
	bool SetOffsets = TRUE;
  if(filter == 0)
  {
    read_las_header(las_file,SetOffsets);
    read_las_data(las_file);
  }
  else if(filter == 1)        
  {
    read_las_header(las_file,SetOffsets);
    read_las_data_vegetation(las_file);
  }
  else if(filter == 2)
  {
    read_las_header(las_file,SetOffsets);
    read_las_data_ground(las_file);
  }
  else
  {
  	cout << "LSDCloud line 78 Warning you have an unsupported filter-type; return to CloudBase!" << endl;
  	exit(EXIT_FAILURE);
  }
}

// This creates an LSDCloud from a LSDRaster
void LSDCloud::create(LSDRaster& raster)
{
 	OctreeResolution = raster.get_DataResolution();
	raster_to_cloud(raster);
}
// Read in data from a PointData object (see LSDShapeTools module)
void LSDCloud::create(PointData& point_data, LSDRaster& raster)
{
 	OctreeResolution = raster.get_DataResolution();
	PointData_to_cloud(point_data,raster);
}
//------------------------------------------------------------------------------
// READING .las FILES
// Read in header information only.  If origin shift specified (e.g. using lower
// left corner of a LSDRaster), set SetOffsetsFromTile to FALSE.  If origin
// shift unspecified, use lower left corner of las tile.
void LSDCloud::read_las_header(string in_file, bool SetOffsetsFromTile)
{
  ifstream ifs;
  if (!liblas::Open(ifs, in_file.c_str()))
  {
        std::cout << "Can not open file" << std::endl;;
        exit(1);
  }
  // Set up liblas .las file reader
  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  // Check number of points from header file
  liblas::Header const& header = reader.GetHeader();
//  cout << "\t Points count: " << header.GetPointRecordsCount() << '\n';
  NPts = header.GetPointRecordsCount();
  
  if(SetOffsetsFromTile == TRUE)
  {
    XOffset = header.GetMinX();
    YOffset = header.GetMinY();
  }
  XMin = header.GetMinX() - XOffset;
  YMin = header.GetMinY() - YOffset;
  XMax = header.GetMaxX() - XOffset;
  YMax = header.GetMaxY() - YOffset;
}

// Read in data from .las file.
// Load point cloud, including ground returns and vegetation returns
void LSDCloud::read_las_data(string in_file)
{
  ifstream ifs;
  if (!liblas::Open(ifs, in_file.c_str()))
  {
        std::cout << "Can not open file" << std::endl;;
        exit(1);
  }
  // Set up liblas .las file reader
  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  // Check number of points from header file
  liblas::Header const& header = reader.GetHeader();
  // READ DATA INTO PCL POINT CLOUD
  // I have defined the point structure as x,y,z,intensity.  Note that I have made z=0 to
  // collapse the 3D pointcloud onto a 2D plane, so that I can make 2D range searches on 
  // the x-y plane.  Elevation values are stored in the intensity field
  //----------------------------------------------------------------------------
  // Loop through the .las file, reading the points into the new data structure.
  vector<float> x_coordinates;
  vector<float> y_coordinates;
  vector<float> zeta_values;
    
  while (reader.ReadNextPoint())
  {
    liblas::Point const& p = reader.GetPoint();
    x_coordinates.push_back( p.GetX() - XOffset );
    y_coordinates.push_back( p.GetY() - YOffset );
    zeta_values.push_back( p.GetZ() );
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

  // Generate octree structure
  cloudOctree = (new pcl::octree::OctreePointCloudSearch<pcl::PointXYZI>(OctreeResolution) );
  cloudOctree->setInputCloud(theCloud);
  cloudOctree->addPointsFromInputCloud();
}

void LSDCloud::read_las_data_vegetation(string in_file)
{
  ifstream ifs;
  if (!liblas::Open(ifs, in_file.c_str()))
  {
        std::cout << "Can not open file" << std::endl;;
        exit(1);
  }
  // Set up liblas .las file reader
  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  // Check number of points from header file
  liblas::Header const& header = reader.GetHeader();
  
  // Make list of classes to be selected for by the filter
  vector<liblas::Classification> classes;
  classes.push_back(liblas::Classification(0)); // unclassified/vegetation
  classes.push_back(liblas::Classification(1)); // Unassigned/vegetation
  classes.push_back(liblas::Classification(3)); // low vegetation
  classes.push_back(liblas::Classification(4)); // medium vegetation
  classes.push_back(liblas::Classification(5)); // high vegetation
  
  // Construct filter  
  vector<liblas::FilterPtr> filters;
  liblas::FilterPtr class_filter = liblas::FilterPtr(new liblas::ClassificationFilter(classes));
  
  class_filter->SetType(liblas::FilterI::eInclusion); // eInclusion means to keep the classes that match (eExclusion would throw out those that matched).
  filters.push_back(class_filter);

  // Apply filters to the reader function
  reader.SetFilters(filters);
  
  // Loop through the .las file, reading the points into the new data structure.
  vector<float> x_coordinates;
  vector<float> y_coordinates;
  vector<float> zeta_values;
  
  while (reader.ReadNextPoint())
  {
    liblas::Point const& p = reader.GetPoint();
    x_coordinates.push_back( p.GetX() - XOffset );
    y_coordinates.push_back( p.GetY() - YOffset );
    zeta_values.push_back( p.GetZ() );
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

void LSDCloud::read_las_data_ground(string in_file)
{
  ifstream ifs;
  if (!liblas::Open(ifs, in_file.c_str()))
  {
        std::cout << "Can not open file" << std::endl;;
        exit(1);
  }
  // Set up liblas .las file reader
  liblas::ReaderFactory f;
  liblas::Reader reader = f.CreateWithStream(ifs);
  // Check number of points from header file
  liblas::Header const& header = reader.GetHeader();
  // Make list of classes to be selected for by the filter
  vector<liblas::Classification> classes;
  classes.push_back(liblas::Classification(2)); // Ground
  // Construct filter  
  vector<liblas::FilterPtr> filters;
  liblas::FilterPtr class_filter = liblas::FilterPtr(new liblas::ClassificationFilter(classes));
  
  class_filter->SetType(liblas::FilterI::eInclusion); // eInclusion means to keep the classes that match (eExclusion would throw out those that matched).
  filters.push_back(class_filter);

  // Apply filters to the reader function
  reader.SetFilters(filters);
  
  // Loop through the .las file, reading the points into the new data structure.
  vector<float> x_coordinates;
  vector<float> y_coordinates;
  vector<float> zeta_values;
  
  while (reader.ReadNextPoint())
  {
    liblas::Point const& p = reader.GetPoint();
    x_coordinates.push_back( p.GetX() - XOffset );
    y_coordinates.push_back( p.GetY() - YOffset );
    zeta_values.push_back( p.GetZ() );
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

// Read in data from an LSDRaster
void LSDCloud::raster_to_cloud(LSDRaster& raster)
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
void LSDCloud::PointData_to_cloud(PointData& point_data, LSDRaster& raster)
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
  void LSDCloud::RadiusSearch2D(float searchPoint_x, float searchPoint_y, float searchRadius,
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
  void LSDCloud::NearestNeighbourSearch2D(float searchPoint_x, float searchPoint_y, int K,
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
