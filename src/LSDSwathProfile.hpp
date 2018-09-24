//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDSwathProfile.hpp
//------------------------------------------------------------------------------
// This code houses the LSDSwath object, used to make swath profiles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 0.0.2		16/02/2018
// Prerequisite software packages: TNT, PCL and liblas
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDCloudRaster.hpp"
#include "LSDShapeTools.hpp"
#include "LSDIndexChannel.hpp"
// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
//#include <liblas/liblas.hpp>
using namespace std;
using namespace TNT;

#ifndef LSDSwathProfile_H
#define LSDSwathProfile_H

/// @brief This code houses the LSDSwath object, used to make swath profiles
/// @author DTM
/// @date 17/02/14
class LSDSwath
{
  public:
  LSDSwath()	{ create(); }
  ///@brief create an LSDSwath using a raster as a template.
  ///
  ///@param PointData ProfilePoints -> coordinates of points forming the
  /// baseline of the profile.
  ///@param LSDRaster RasterTemplate -> a raster dataset that is used as a
  /// template for the swath profile i.e. any rasters that you wish to generate
  /// the profile for should have the same characteristics/spatial extent as the
  /// original template.
  ///@param float ProfileHalfWidth
  ///@author DTM
  ///@date 11/04/2014
  ///
  LSDSwath(PointData& ProfilePoints, LSDRaster& RasterTemplate, float HalfWidth) { create(ProfilePoints, RasterTemplate, HalfWidth); }

  ///@brief create an LSDSwath using a raster as a template.
  ///
  ///@param vector<vector<float> >& Y_X_points -> coordinates of points
  /// forming the extremities of the baseline of the profile.
  ///@param LSDRaster RasterTemplate -> a raster dataset that is used as a
  /// template for the swath profile i.e. any rasters that you wish to generate
  /// the profile for should have the same characteristics/spatial extent as the
  /// original template.
  ///@param float ProfileHalfWidth
  ///@author DTM
  ///@date 28/02/2017
  ///
  LSDSwath(vector<float>& Y_X_points, LSDRaster& RasterTemplate, float& HalfWidth, float d_space) { create(Y_X_points, RasterTemplate, HalfWidth, d_space); }

  /// @brief This moves perpendicular to the swath and gets statisics from bins of distance from the swath
  ///  you can tell it what widths you want and returns a series of swaths at the mean and then percentiles you want
  /// @param Raster The raster from which you extract data
  /// @param desired_percentiles A float vector that contains the percentiles of the data you want to calculate for your output profiles
  /// @param BinWidth The width of each bin
  /// @param mid_points A vector of the mid points, in terms of distance along, of the bins. This is overwritten during computation.
  /// @param mean_profile A vector with the mean profile values. This is overwritten during computation.
  /// @param sd_profile A vector with the standard deviation values. This is overwritten during computation.
  /// @param output_percentile_profiles A vec vector where each vector is the profile at the precentile determined in the desired_percentiles vector.
  ///    This is overwritten during computation.
  /// @param NormaliseToBaseline 0 if you want the raw data, 1 if you want the data normalised to the value at the baseline.
  /// @author DTM
  /// @date 01/01/2015
  void get_transverse_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles,
       int NormaliseToBaseline);

  /// @brief This moves along the swath and gets statisics from bins of distance along the swath
  ///  you can tell it what widths you want and returns a series of swaths at the mean and then percentiles you want
  /// @param Raster The raster from which you extract data
  /// @param desired_percentiles A float vector that contains the percentiles of the data you want to calculate for your output profiles
  /// @param BinWidth The width of each bin
  /// @param mid_points A vector of the mid points, in terms of distance along, of the bins. This is overwritten during computation.
  /// @param mean_profile A vector with the mean profile values. This is overwritten during computation.
  /// @param sd_profile A vector with the standard deviation values. This is overwritten during computation.
  /// @param output_percentile_profiles A vec vector where each vector is the profile at the precentile determined in the desired_percentiles vector.
  ///    This is overwritten during computation.
  /// @param NormaliseToBaseline 0 if you want the raw data, 1 if you want the data normalised to the value at the baseline.
  ///  In general 1 is used to see the value along the baseline relative to that baseline so good for looking at things like
  ///   terraces along a river.
  /// @author DTM
  /// @date 01/01/2015
  void get_longitudinal_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles,
       int NormaliseToBaseline);

  ///@brief create a raster in the shape of the swath profile
  ///@param Raster LSDRaster of interest
  ///@param NormaliseToBaseline if 0 --> raster values; if 1 --> raster values normalised to baseline value
  ///@return raster in shape of swath profile
  ///@author FJC
  ///@date 16/10/15
  LSDRaster get_raster_from_swath_profile(LSDRaster& Raster, int NormaliseToBaseline);

  ///@brief fill in the baseline raster value with the average value of the pixels along the transverse swath profile
  ///@param Raster Raster template
  ///@return LSDRaster with filled in values along the baseline
  ///@author FJC
  ///@date 16/01/17
  LSDRaster fill_in_channels_swath(LSDRaster& Raster);

  ///@brief Get information about connected components along the swath profile
  ///@details This function takes in a connected components raster and returns the average
  /// value of another chosen raster and distance along the baseline of each
  /// connected components patch. The user can choose whether to normalise the
  /// second raster to the baseline value.
  /// vector of vectors has the format:
  /// 0 = patch ID
  /// 1 = mean raster value for the patch id
  /// 2 = mean distance along the baseline
  ///@return vector of vector with patch ids, mean values, and distance along baseline
  ///@author FJC
  ///@date 24/01/17
  vector <vector <float> > get_connected_components_along_swath(LSDIndexRaster& ConnectedComponents, LSDRaster& RasterTemplate, int NormaliseToBaseline);

  ///@brief  This function takes in a raster and returns the mean, min and max values of the raster
  /// at each point along the swath
  /// vector of vectors has the format:
  /// 0 = distance along swath
  /// 1 = mean value along swath
  /// 2 = min value along swath
  /// 3 = max value along swath
  /// if NormaliseToBaseline == 1 then the values will be normalised to the baseline.
  ///@return vector of vectors
  ///@author FJC
  ///@date 15/02/17
  vector <vector <float> > get_RasterValues_along_swath(LSDRaster& RasterTemplate, int NormaliseToBaseline);

  /// @brief  Wrapper to write raster values along swath to csv
  /// @param RasterTemplate raster of values to write
  /// @param NormaliseToBaseline int, if 1 then will be noramlised to the baseline value.
  ///     if 0 then you just get the raw underlying data
  /// @param csv_filename file name of output csv
  ///@author FJC
  ///@date 20/11/17
  void write_RasterValues_along_swath_to_csv(LSDRaster& RasterTemplate, int NormaliseToBaseline, string csv_filename);

  ///@brief  This function takes in a connected components raster and returns an array
  ///@param ConnectedComponents connected components raster
  ///@return array with baseline components
  ///@author FJC
  ///@date 28/09/17
  Array2D<float> get_BaselineDist_ConnectedComponents(LSDIndexRaster& ConnectedComponents);

  ///@brief  This function takes in a connected components raster and returns an array
  /// of the distance to the baseline for each pixl in the raster
  ///@param ConnectedComponents connected components raster
  ///@return array with baseline components
  ///@author FJC
  ///@date 12/10/17
  Array2D<float> get_DistanceToBaseline_ConnectedComponents(LSDIndexRaster& ConnectedComponents);

  ///@brief  This function takes in a raster for analysis and gets the width of pixels in that raster
  ///@param RasterForAnalysis analysis index raster
  ///@return vector with widths along baseline
  ///@author FJC
  ///@date 21/11/17
  vector<float> get_widths_along_swath(LSDIndexRaster& RasterForAnalysis);

  // write profiles to file
  void write_transverse_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix, int NormaliseToBaseline);
  void write_longitudinal_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix, int NormaliseToBaseline);

  ///@brief This prints the baseline to a csv file
  ///@param ElevationRaster the raster
  ///@param csv_filename The name of the file out
  ///@author FJC
  ///@date 12/10/17
  void print_baseline_to_csv(LSDRaster& ElevationRaster, string csv_filename, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet);


  /// @brief This finds locations of points along a baseline
  /// @param A vector of distances along the baseline
  /// @author SMM
  /// @date 16/02/2018
  void get_locations_of_points_along_baseline(vector<float> distance_along_baseline);

  /// @brief This prints the DistanceToBaselineArray to a raster
  /// @param Raster a raster to which the values will be printed to its geometry
  /// @author SMM
  /// @date 16/02/2018
  LSDRaster get_raster_DistanceToBaselineArray(LSDRaster& Raster);

  /// @brief This prints the DistanceAlongBaselineArray to a raster
  /// @param Raster a raster to which the values will be printed to its geometry
  /// @author SMM
  /// @date 16/02/2018
  LSDRaster get_raster_DistanceAlongBaselineArray(LSDRaster& Raster);

  /// @brief This prints the BaselineValueArray to a raster
  /// @param Raster a raster to which the values will be printed to its geometry
  /// @author SMM
  /// @date 16/02/2018
  LSDRaster get_raster_BaselineValueArray(LSDRaster& Raster);



  // get functions
  // these get data elements
  int get_NPtsInProfile() const {return NPtsInProfile;}
  Array2D<float> get_DistanceToBaselineArray() const { return DistanceToBaselineArray; }
  Array2D<float> get_DistanceAlongBaselineArray() const { return DistanceAlongBaselineArray; }
  Array2D<float> get_BaselineValueArray() const { return BaselineValueArray; }
  vector<float> get_DistanceAlongBaseline() const { return DistanceAlongBaseline; }

  float get_XMax() const { return XMax; }
  float get_YMax() const { return YMax; }
  float get_XMin() const { return XMin; }
  float get_YMin() const { return YMin; }
  float get_ProfileHalfWidth() const { return ProfileHalfWidth; }

  vector<int> get_BaselineCols() const { return BaselineCols; }
  vector<int> get_BaselineRows() const { return BaselineRows; }
  vector<float> get_BaselineValue() const { return BaselineValue; }

  protected:

  // Swath template
  vector<float> DistanceAlongBaseline;
  vector<float> BaselineValue;
  vector<int> BaselineRows;  // rows of the baseline points
  vector<int> BaselineCols;  // cols of the baseline points
  Array2D<float> DistanceToBaselineArray;
  Array2D<float> DistanceAlongBaselineArray;
  Array2D<float> BaselineValueArray;

	// metadata
  int NPtsInProfile;
  float ProfileHalfWidth;
  float NoDataValue;
  int NRows;
  int NCols;

  // Bounding Box of profile
  float XMax;
  float XMin;
  float YMax;
  float YMin;

  private:
  void create();
  void create(PointData& ProfilePoints, LSDRaster& RasterTemplate, float ProfileHalfWidth);
  void create(vector<float>& Y_X_points, LSDRaster& RasterTemplate, float& ProfileHalfWidth, float d_space);



};

#endif
