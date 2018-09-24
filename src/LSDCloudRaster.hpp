//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDCloudRaster.hpp
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
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDShapeTools.hpp"
// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
using namespace std;
using namespace TNT;

#ifndef LSDCloudRaster_H
#define LSDCloudRaster_H

/// @brief This code houses the LSDCloudRaster object, designed to create clouds from
/// raster objects
/// @author DTM, FJC
/// @date 28/03/17
class LSDCloudRaster
{
  public:
  // Assignment operator.
	LSDCloudRaster& operator=(const LSDCloudRaster& cloud_object);

  LSDCloudRaster()	{ create(); }
  ///@brief create an LSDCloud from a raster.
  ///
  LSDCloudRaster(LSDRaster& raster)	 { create(raster); }

  ///@brief create an LSDCloud from a PointData object (see LSDShapeTools module)
  LSDCloudRaster(PointData& point_data, LSDRaster& raster) { create(point_data,raster); }

  LSDCloudRaster(int n_pts, float x_min, float x_max, float y_min, float y_max, float x_offset, float y_offset, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, const double resolution)
            { create( n_pts, x_min, x_max, y_min, y_max, x_offset, y_offset, cloud, resolution); }


  // get functions
  // these get data elements
  int get_NPts() const {return NPts;}
  float get_XMin() const			{ return XMin; }
  float get_YMin() const			{ return YMin; }
  float get_XMax() const			{ return XMax; }
  float get_YMax() const			{ return YMax; }
  float get_XOffset() const			{ return XOffset; }
  float get_YOffset() const			{ return YOffset; }
  pcl::PointCloud<pcl::PointXYZI>::Ptr get_theCloud() const { return theCloud; }
  float get_OctreeResolution() const   		{ return OctreeResolution; }

  float get_point_x(int index) { return theCloud->points[index].x; }
  float get_point_y(int index) { return theCloud->points[index].y; }

  /// @brief Read in data from an LSDRaster into a LSDCloud object
  ///
  /// @details The point cloud is also structured in an octree for rapid spatial
  /// querying. Uses the PCL library (www.poinclouds.org).
  /// @param LSDRaster -> raster to convert to a LSDCloud
  /// @author DTM
  /// @date 18/02/2014
  void raster_to_cloud(LSDRaster& raster);

  /// @brief Read in data from a PointData (LSDShapeTools) into a LSDCloud
  ///
  /// @details The point cloud is also structured in an octree for rapid spatial
  /// querying. Uses the PCL library (www.poinclouds.org).
  /// @param PointData -> x,y coordinate strucure to load into cloud
  /// @param LSDRaster -> raster used as spatial reference frame
  /// @author DTM
  /// @date 11/04/2014
  void PointData_to_cloud(PointData& point_data, LSDRaster& raster);
  ///---------------------------------------------------------------------------
  /// SEARCH TOOLS

  /// @brief Conducts a focal search (in xy plane), to find all points within
  /// defined search radius
  ///
  /// @details Points are arranged according to proximity to the query point
  /// @param float -> the x coordinate of the query point
  /// @param float -> the y coordinate of the query point
  /// @param float -> the search radius
  /// @param vector<float> -> the values associated with the returned points
  /// @param vector<int> -> the indexes associated with the returned points
  /// @param vector<float> -> the squared distance of the returned points from the query point
  /// @author DTM
  /// @date 18/02/2014
  void RadiusSearch2D(float searchPoint_x, float searchPoint_y, float searchRadius,
                              vector<float>& pointValues, vector<int>& pointIndex, vector<float>& pointSquaredDistance);
  /// @brief Conducts a focal search (in the xy plane), to find the nearest K
  /// points
  ///
  /// @details Points are arranged according to proximity to the query point
  /// @param float -> the x coordinate of the query point
  /// @param float -> the y coordinate of the query point
  /// @param int -> the number of points to find
  /// @param vector<float> -> the values associated with the returned points
  /// @param vector<int> -> the indexes associated with the returned points
  /// @param vector<float> -> the squared distance of the returned points from the query point
  /// @author DTM
  /// @date 18/02/2014
  void NearestNeighbourSearch2D(float searchPoint_x, float searchPoint_y, int K,
                              vector<float>& pointValues, vector<int>& pointIndex, vector<float>& pointSquaredDistance);


	protected:

  pcl::PointCloud<pcl::PointXYZI>::Ptr theCloud;

  // Octree spatial structure
  double OctreeResolution;
  pcl::octree::OctreePointCloudSearch<pcl::PointXYZI> *cloudOctree;

  // GEOREFERENCING
  // bounding box of point cloud
  float XMin;
  float YMin;
  float XMax;
  float YMax;
  // origin shift
  float XOffset;
  float YOffset;

	// metadata
  int NPts;

	private:
  void create();
  void create(LSDRaster& raster);
  void create(int n_pts, float x_min, float x_max, float y_min, float y_max,
              float x_offset, float y_offset, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, const double resolution);
  void create(PointData& point_data, LSDRaster& raster);
};

#endif
