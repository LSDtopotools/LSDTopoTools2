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

#ifndef LSDBasin_H
#define LSDBasin_H

#include <vector>
#include <string>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDParticle.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDSpatialCSVReader.hpp"
using namespace std;
using namespace TNT;

///@brief Object to store information about drainage basins and generate basin average metrics..
class LSDBasin
{

  public:

  /// @brief Default constructor method used to create a basin object.
  ///
  /// @details This is necessary since there is a derived class and with derived
  ///  classes the default constructor of the parent class is automatically called
  /// @author SMM
  /// @date 28/12/14
  LSDBasin()        { create(); }

  /// @brief Constructor method used to create a basin object.
  /// @param Junction outlet junction of the basin to be constructed.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param ChanNet Channel network object.
  /// @author SWDG
  /// @date 11/12/12
  LSDBasin(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet)
                { create(Junction, FlowInfo, ChanNet); }

  LSDBasin(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet, int extended_baslevel_nodes_below_junction)
  {create(JunctionNumber,FlowInfo,ChanNet,extended_baslevel_nodes_below_junction);}


  /// @return Number of rows as an integer.
  int get_NRows() const        { return NRows; }
  /// @return Number of columns as an integer.
  int get_NCols() const        { return NCols; }
  /// @return Minimum X coordinate as an integer.
  float get_XMinimum() const        { return XMinimum; }
  /// @return Minimum Y coordinate as an integer.
  float get_YMinimum() const        { return YMinimum; }
  /// @return Data resolution as an integer.
  float get_DataResolution() const        { return DataResolution; }
  /// @return No Data Value as an integer.
  int get_NoDataValue() const        { return NoDataValue; }

  /// @return Georeferencing information
  map<string,string> get_GeoReferencingStrings() const { return GeoReferencingStrings; }

  /// @return Nodes of the basin.
  vector<int> get_BasinNodes() const { return BasinNodes; }

  /// @return Number of cells in the basin.
  int get_NumberOfCells() const { return NumberOfCells; }

  /// @return Area of basin in spataial units.
  float get_Area() const { return Area; }

  /// @return Junction Number of the basin.
  int get_Junction() const { return Junction; }

  /// @return Stream order of basin.
  int get_BasinOrder() const { return BasinOrder; }

  /// @return Boolean value of whether a basin is beheaded or not.
  bool get_Beheaded() const { return Beheaded; }

  /// @return the node index of the outlet
  int get_Outlet_node() const{ return BasinNodes[0]; }

  /// @return i index of outlet pixel.
  int get_Outlet_i() const { return Outlet_i; }

  /// @return j index of outlet pixel.
  int get_Outlet_j() const { return Outlet_j; }

  ///@return i index of Centroid pixel.
  int get_Centroid_i() const { return Centroid_i; }

  ///@return j index of Centroid pixel.
  int get_Centroid_j() const { return Centroid_j; }

  //Getters of basin parameters

  /// @return Mean slope.
  float get_SlopeMean() const { return SlopeMean; }
  /// @return Mean elevation.
  float get_ElevationMean() const { return ElevationMean; }
  /// @return Mean aspect.
  float get_AspectMean() const { return AspectMean; }
  /// @return Mean relief.
  float get_ReliefMean() const { return ReliefMean; }
  /// @return Mean plan curvature.
  float get_PlanCurvMean() const { return PlanCurvMean; }
  /// @return Mean profile curvature.
  float get_ProfileCurvMean() const { return ProfileCurvMean; }
  /// @return Mean total curvature.
  float get_TotalCurvMean() const { return TotalCurvMean; }
  /// @return Max plan curvature.
  float get_PlanCurvMax() const { return PlanCurvMax; }
  /// @return Max profile curvature.
  float get_ProfileCurvMax() const { return ProfileCurvMax; }
  /// @return Max total curvature.
  float get_TotalCurvMax() const { return TotalCurvMax; }
  /// @return Hillslope length from hilltop flow routing.
  float get_HillslopeLength_HFR() const { return HillslopeLength_HFR; }
  /// @return Hillslope length from boomerang bins.
  float get_HillslopeLength_Binned() const { return HillslopeLength_Binned; }
  /// @return Hillslope length from boomerang splines.
  float get_HillslopeLength_Spline() const { return HillslopeLength_Spline; }
  /// @return Hillslope length from drainage density.
  float get_HillslopeLength_Density() const { return HillslopeLength_Density; }
  /// @return Flow length.
  float get_FlowLength() const { return FlowLength; }
  /// @return Drainage Density.
  float get_DrainageDensity() const { return DrainageDensity; }
  /// @return Basin Perimeter's i index.
  vector<int> get_Perimeter_i() const { return Perimeter_i; }
  /// @return Basin Perimeter's j index.
  vector<int> get_Perimeter_j() const { return Perimeter_j; }
  /// @return Basin Perimiter's node indices
  vector<int> get_Perimeter_nodes() const { return Perimeter_nodes; };
  /// @return Cosmo erosion rate.
  float get_CosmoErosionRate() const { return CosmoErosionRate; }
  /// @return Other eroision rate.
  float get_OtherErosionRate() const { return OtherErosionRate; }
  /// @return Mean hilltop curvature.
  float get_CHTMean() const { return CHTMean; }
  /// @return E* value.
  float get_EStar() const { return EStar; }
  /// @return R* value.
  float get_RStar() const { return RStar; }
  /// @return Number of hilltop pixels in the basin.
  int get_HilltopPx() const { return HilltopPx; }
  /// @return fraction basin that is rock rather than soil mantled
  float get_BedrockFraction() const { return BedrockFraction; };
  float get_Biomass() const { return Biomass; };
  int get_AlternativeIndex() const { return AlternativeIndex; };
  /// @brief Calculate the mean value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the mean of.
  /// @return Mean value.
  /// @author SWDG
  /// @date 11/12/13
  float CalculateBasinMean(LSDFlowInfo& FlowInfo, LSDRaster Data);
  float CalculateBasinMean(LSDFlowInfo& FlowInfo, LSDIndexRaster Data);

  /// @brief Calculate the max value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the max of.
  /// @return Max value.
  /// @author SWDG
  /// @date 11/12/13
  float CalculateBasinMax(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the min value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the minimum of.
  /// @return Min value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinMin(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the median value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the median of.
  /// @return Median value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinMedian(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the percentile value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the percentile of.
  /// @return Percentile value.
  /// @author SWDG
  /// @date 5/2/17
  float CalculateBasinPercentile(LSDFlowInfo& FlowInfo, LSDRaster Data, int Percentile);

  /// @brief Calculate the Standard Deviation of values of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the standard deviation of.
  /// @return Standard deviation value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinStdDev(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the Standard error of values of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the standard error of.
  /// @return Standard error value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinStdError(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the range of values of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the range of.
  /// @return Range value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinRange(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the number of data points of an LSDRaster which fall inside a basin.
  ///
  /// @details Useful for checking that an average value is not just taken from a single data point.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to count.
  /// @return Number of data points.
  /// @author SWDG
  /// @date 17/2/14
  int CalculateNumDataPoints(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Set the mean slope of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Slope Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_SlopeMean(LSDFlowInfo& FlowInfo, LSDRaster Slope){ SlopeMean = CalculateBasinMean(FlowInfo, Slope); }

  /// @brief Set the mean Elevation of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Elevation Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ElevationMean(LSDFlowInfo& FlowInfo, LSDRaster Elevation) { ElevationMean = CalculateBasinMean(FlowInfo, Elevation); }

  /// @brief Set the mean Relief of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Relief Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ReliefMean(LSDFlowInfo& FlowInfo, LSDRaster Relief) { ReliefMean = CalculateBasinMean(FlowInfo, Relief); }

  /// @brief Set the mean PlanCurve of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param PlanCurv Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_PlanCurvMean(LSDFlowInfo& FlowInfo, LSDRaster PlanCurv) { PlanCurvMean = CalculateBasinMean(FlowInfo, PlanCurv); }

  /// @brief Set the mean ProfCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param ProfileCurv Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ProfileCurvMean(LSDFlowInfo& FlowInfo, LSDRaster ProfileCurv) { ProfileCurvMean = CalculateBasinMean(FlowInfo, ProfileCurv); }

  /// @brief Set the mean TotalCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param TotalCurv Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_TotalCurvMean(LSDFlowInfo& FlowInfo, LSDRaster TotalCurv) { TotalCurvMean = CalculateBasinMean(FlowInfo, TotalCurv); }

  /// @brief Set the max PlanCurve of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param PlanCurv Values to find the max of.
  /// @author SWDG
  /// @date 11/12/13
  void set_PlanCurvMax(LSDFlowInfo& FlowInfo, LSDRaster PlanCurv) { PlanCurvMax = CalculateBasinMax(FlowInfo, PlanCurv); }

  /// @brief Set the max ProfCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param ProfileCurv Values to find the max of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ProfileCurvMax(LSDFlowInfo& FlowInfo, LSDRaster ProfileCurv) { ProfileCurvMax = CalculateBasinMax(FlowInfo, ProfileCurv); }

  /// @brief Set the max TotalCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param TotalCurv Values to find the max of.
  /// @author SWDG
  /// @date 11/12/13
  void set_TotalCurvMax(LSDFlowInfo& FlowInfo, LSDRaster TotalCurv) { TotalCurvMax = CalculateBasinMax(FlowInfo, TotalCurv); }

  /// @brief Set the mean hilltop curvature of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param CHT Values to find the mean of.
  /// @author SWDG
  /// @date 12/12/13
  void set_CHTMean(LSDFlowInfo& FlowInfo, LSDRaster CHT) { CHTMean = CalculateBasinMean(FlowInfo, CHT); }

  /// @brief Set the Cosmogenic erosion rate.
  /// @param ErosionRate Erosion rate - No sanity check on this value.
  /// @author SWDG
  /// @date 11/12/13
  void set_CosmoErosionRate(float ErosionRate) { CosmoErosionRate = ErosionRate; }

  /// @brief Set the Other erosion rate.
  /// @param ErosionRate Erosion rate - No sanity check on this value.
  /// @author SWDG
  /// @date 11/12/13
  void set_OtherErosionRate(float ErosionRate) { OtherErosionRate = ErosionRate; }

  /// @brief Calculate E* and R* values for the basin, using hilltop flow routed hillslope lengths.
  /// @param CriticalSlope slope threshold value, typically 0.4.
  /// @author SWDG
  /// @date 12/12/13
  void set_EStar_RStar(float CriticalSlope);

  /// @brief Calculate flow length for the basin using the D8 flow directions.
  /// @param StreamNetwork the channel network.
  /// @param FlowInfo Flowinfo object.
  /// @author SWDG
  /// @date 12/12/13
  void set_FlowLength(LSDIndexRaster& StreamNetwork, LSDFlowInfo& FlowInfo);

  /// @brief Set basin drainage density.
  /// @author SWDG
  /// @date 12/12/13
  void set_DrainageDensity() { DrainageDensity = FlowLength/Area ; }

  /// @brief Set the mean HillslopeLength from hilltop flow routing.
  /// @param FlowInfo Flowinfo object.
  /// @param HillslopeLengths Values to find the mean of.
  /// @author SWDG
  /// @date 12/12/13
  void set_HillslopeLength_HFR(LSDFlowInfo& FlowInfo, LSDRaster HillslopeLengths) { HillslopeLength_HFR = CalculateBasinMean(FlowInfo, HillslopeLengths); }

  /// @brief Set mean HillslopeLengths from boomerang plots from both splines and binned data.
  /// @param Slope LSDRaster of slope.
  /// @param DinfArea D-infinity Flowarea LSDRaster.
  /// @param FlowInfo Flowinfo object.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @param SplineResolution Number of values between each point for the spline curve.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin.
  /// @author SWDG
  /// @date 12/12/13
  void set_HillslopeLengths_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold);

  /// @brief Set the mean HillslopeLength from drainage density.
  /// @author SWDG
  /// @date 12/12/13
  void set_HillslopeLength_Density() { HillslopeLength_Density = (1 / (2 * DrainageDensity)); }

  /// @brief Set the Rock Exposure fraction of the basin
  /// @author DTM
  /// @date 14/07/15
  void set_BedrockFraction(LSDFlowInfo& FlowInfo, LSDRaster RockExposure) { BedrockFraction = CalculateBasinMean(FlowInfo, RockExposure); }

  void set_Biomass(LSDFlowInfo& FlowInfo, LSDRaster BiomassRaster) { Biomass = CalculateBasinMean(FlowInfo, BiomassRaster); }

  void set_AlternativeIndex(LSDFlowInfo& FlowInfo, LSDIndexRaster& AltIndex);


  /// @brief Generate text files containing data to plot boomerangs.
  ///
  /// @details Writes 3 files to the output path, coded with the basin's unique
  /// junction number which can the be read with python and plotted.
  /// @param Slope LSDRaster of slope.
  /// @param DinfArea D-infinity Flowarea LSDRaster.
  /// @param FlowInfo Flowinfo object.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @param SplineResolution Number of values between each point for the spline curve.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin.
  /// @param Path The output path where the data files will be written to, including the final slash.
  /// @author SWDG
  /// @date 12/12/13
  void Plot_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold, string Path);

  /// @brief Set the mean aspect of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Aspect Values to find the mean of.
  /// @author SWDG
  /// @date 17/2/14
  void set_AspectMean(LSDFlowInfo& FlowInfo, LSDRaster Aspect);

  /// @brief Set the perimeter pixels using a simple edge detection algorithm.
  ///
  /// @details This is quite messy and will be improved soon.
  /// @param FlowInfo Flowinfo object.
  /// @author SWDG
  /// @date 12/12/13
  void set_Perimeter(LSDFlowInfo& FlowInfo);

  /// @brief Set the perimeter pixels by passing in
  /// your own vector of perimeter nodes.
  /// @param perimeter_nodes vector of perimeter nodes
  /// @author FJC
  /// @date 26/01/18
  void set_perimeter_from_vector(vector<int> perimeter_nodes) { Perimeter_nodes = perimeter_nodes; }

  /// @brief Prints the perimeter nodes to a csv file
  /// @param FlowInfo the LSDFlowInfo object
  /// @param string perimeter_fname
  /// @author SMM
  /// @date 26/04/2017
  void print_perimeter_to_csv(LSDFlowInfo& FlowInfo, string perimeter_fname);

  /// @brief Prints the perimeter nodes to a csv file plus elevations
  /// @param FlowInfo the LSDFlowInfo object
  /// @param string perimeter_fname
  /// @param perimeter_nodes vector of perimeter nodes that can be passed. Pass an empty vector if you want to use
  /// the default perimeter finder.
  /// @param ElevationRaster elevation raster
  /// @author FJC
  /// @date 10/01/18
  void print_perimeter_hypsometry_to_csv(LSDFlowInfo& FlowInfo, string perimeter_fname, LSDRaster& ElevationRaster);

  /// @brief Orders perimeter nodes from the outlet
  /// @param FlowInfo the LSDFlowInfo object
  /// @author FJC
  /// @date 16/01/18
  vector<int> order_perimeter_nodes(LSDFlowInfo& FlowInfo);

  /// @brief Set the four different hillslope length measurements for the basin.
  /// @param FlowInfo Flowinfo object.
  /// @param HillslopeLengths LSDRaster of hillslope lengths from the hilltop flow routing method.
  /// @param Slope LSDRaster of slope.
  /// @param DinfArea D-infinity Flowarea LSDRaster.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @param SplineResolution Number of values between each point for the spline curve.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin.
  /// @author SWDG
  /// @date 12/12/13
  void set_all_HillslopeLengths(LSDFlowInfo& FlowInfo, LSDRaster& HillslopeLengths, LSDRaster& Slope,
                                LSDRaster& DinfArea, float log_bin_width, int SplineResolution, float bin_threshold);

  /// @brief Set all of the basin parameters with one call.
  ///
  /// @details Runs polyfit to get the elevation derivatives, so can be quite memory intensive. Method
  /// calls all the setters one by one, to populate all the basin parameters. So a
  /// basin can be created and all it's properties set with 2 calls. The erosion rates have default
  /// parameters of -9999 as these are rarely used variables.
  /// @param Elevation LSDRaster of filled elevation values.
  /// @param FlowInfo Flowinfo object.
  /// @param CHT LSDRaster of hilltop curvatures.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param HillslopeLengths LSDRaster of hillslope lengths from the hilltop flow routing method.
  /// @param Relief LSDRaster of the hilltop relief.
  /// @param window_radius Radius in spatial units for the polyft routine.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf. Default value is 0.1.
  /// @param SplineResolution Number of values between each point for the spline curve. Default value is 10000.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin. Default value is 0.
  /// @param CriticalSlope Slope threshold used for E* R* calculations. Default value is 0.4.
  /// @param CosmoErosionRate Erosion rate from cosmo.
  /// @param OtherErosionRate Erosion rate from another source.
  /// @author SWDG
  /// @date 12/12/13
  void set_All_Parameters(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& CHT, LSDIndexRaster& StreamNetwork,
                          LSDRaster& HillslopeLengths, LSDRaster& Relief, float window_radius, float log_bin_width,
                          int SplineResolution, float bin_threshold, float CriticalSlope,
                          float CosmoErosionRate = -9999, float OtherErosionRate = -9999);

  /// @brief Set the count of the number of hilltop pixels in a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Hilltops a raster of hilltop data.
  /// @author SWDG
  /// @date 18/6/15
  void set_HilltopPx(LSDFlowInfo& FlowInfo, LSDRaster Hilltops);

  /// @brief Cookie cut data from an LSDIndexRaster into the shape of the basin.
  /// @param Data LSDIndexRaster data to be written.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDIndexRaster write_raster_data_to_LSDIndexRaster(LSDIndexRaster Data, LSDFlowInfo FlowInfo);

  /// @brief Cookie cut data from an LSDRaster into the shape of the basin.
  /// @param Data LSDRaster data to be written.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_raster_data_to_LSDRaster(LSDRaster& Data, LSDFlowInfo& FlowInfo);

  /// @brief check whether a test node is in the basin or not
  /// @param test_node node to test
  /// @return integer which is 1 if node is in the basin or 0 if it is not
  /// @author FJC
  /// @date 21/02/14
  int is_node_in_basin(int test_node);

  /// @brief This checks the dimensions of the base DEMs along with georeferencing strings
  ///  to see if the other basin is from the same DEM as the current basin
  /// @param other the other LSDBasin
  /// @return true if the other basin if from the same DEM, false otherwise
  /// @author SMM
  /// @date 23/2/2016
  bool are_basins_from_same_base_DEM(LSDBasin& other);

  /// @brief This checks to see if a supplied basin is a subbasin of the other
  /// @param other the other LSDBasin
  /// @return true if the other basin is a subbasin of the first basin
  /// @author SMM
  /// @date 23/2/2016
  bool is_this_a_subbasin(LSDBasin& other);

  /// @brief remove hilltop curvature values that are at the edge of the basin
  /// @param hilltop_curvature raster of CHT
  /// @param FlowInfo Flowinfo object
  /// @return LSDRaster of internal hilltop curvature values
  /// @author FJC
  /// @date 19/03/15
  LSDRaster keep_only_internal_hilltop_curvature(LSDRaster hilltop_curvature, LSDFlowInfo FlowInfo);

  /// @brief Write a real value to an LSDRaster in the shape of the basin.
  /// @param Param real value to be written
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_real_data_to_LSDRaster(float Param, LSDFlowInfo FlowInfo);

  /// @brief Write an integer value to an LSDIndexRaster in the shape of the basin.
  /// @param Param integer value to be written
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDIndexRaster write_integer_data_to_LSDIndexRaster(int Param, LSDFlowInfo FlowInfo);


  /// @brief This function is used to create a single LSDIndexRaster that
  ///  contains all the basins in a DEM. It is done by adding each basin in
  ///  turn (but this only adds a single basin).
  /// @param basin_raster LSDIndexRaster containing the basin codes. It should
  ///  be initialsed using write_integer_data_to_LSDIndexRaster
  /// @param FlowInfo the flowinfo data member
  /// @param drainage_of_other_basins This is a map where the key is an int, which
  ///  has the basin number, and the value is an int, the number of pixels in
  ///  the basin. It is somewhat dangerous because it does not make sure that
  ///  the user is inputting unique basin numbers so some basin numbers could be
  ///  overwritten
  /// @param this_basin_index the index of the basin being added to the raster
  /// @author SMM
  /// @date 18/03/2015
  void add_basin_to_LSDIndexRaster(LSDIndexRaster& basin_raster,
                                           LSDFlowInfo& FlowInfo,
                                           map<int,int>& drainage_of_other_basins,
                                           int this_basin_index);


  /// @brief This function creates a padded, georeferenced raster from a raw
  ///  raster file that is trimmed to the dimension of a basin
  /// @details Includes a padding function so that one can pad either
  ///  to catch streams for CRN analysis or for wider padding to catch
  ///  peaks for topographic shielding (also for CRN analysis)
  /// @param padding_pixels the number of pixels with which to pad the raster
  /// @param FlowInfo an LSDFlowInfo object
  /// @param Raster_Data the raster that gets trimmed. Its georeferencing needs
  ///  to be the same as the basin obect georeferencing.
  /// @author SMM
  /// @date 18/03/2015
LSDRaster TrimPaddedRasterToBasin(int padding_pixels, LSDFlowInfo& FlowInfo,
                                            LSDRaster& Raster_Data);

  /// @brief This function check if two basin are adjacent
  /// @detail return true if the two basin are adjacent with at least one pixel
  ///  TODO add a minimum adjacent pixel parameter
  /// @param LSDBasin another LSDBasin object
  /// @author BG
  /// @date 10/10/2017
  bool is_adjacent(LSDBasin& DifferentBasin, LSDFlowInfo& flowpy);

/// @brief detect the source nodes in a pixel window around a perimeter, for instance a basin  perimeter
/// @detail It needs a sequence of nodes where it will loop around and gather all the source nodes encountered.
/// @param vector of nodes, Flowinfo object and a JunctionNetwork object and a number of pixel for the window.
/// @return vector of node indices of the new perimeter
/// @author BG
/// @date 11/10/17
vector<int> get_source_node_from_perimeter(vector<int> perimeter, LSDFlowInfo& flowpy, LSDJunctionNetwork& junky, int pixel_window);



  /// @brief Write Junction values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of Junction values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDIndexRaster write_Junction(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(Junction, FlowInfo); }

  /// @brief Write NumberOfCells values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of NumberOfCells values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDIndexRaster write_NumberOfCells(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(NumberOfCells, FlowInfo); }

  /// @brief Write BasinOrder values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of BasinOrder values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDIndexRaster write_BasinOrder(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(BasinOrder, FlowInfo); }

  /// @brief Write Area values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of Area values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_Area(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(Area, FlowInfo); }

  /// @brief Write SlopeMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of SlopeMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_SlopeMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(SlopeMean, FlowInfo); }

  /// @brief Write ElevationMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ElevationMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_ElevationMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ElevationMean, FlowInfo); }

  /// @brief Write AspectMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of AspectMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_AspectMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(AspectMean, FlowInfo); }

  /// @brief Write ReliefMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ReliefMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_ReliefMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ReliefMean, FlowInfo); }

  /// @brief Write PlanCurvMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of PlanCurvMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_PlanCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(PlanCurvMean, FlowInfo); }

  /// @brief Write ProfileCurvMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ProfileCurvMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_ProfileCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ProfileCurvMean, FlowInfo); }

  /// @brief Write TotalCurvMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of TotalCurvMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_TotalCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(TotalCurvMean, FlowInfo); }

  /// @brief Write PlanCurvMax values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of PlanCurvMax values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_PlanCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(PlanCurvMax, FlowInfo); }

  /// @brief Write ProfileCurvMax values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ProfileCurvMax values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_ProfileCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ProfileCurvMax, FlowInfo); }

  /// @brief Write TotalCurvMax values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of TotalCurvMax values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_TotalCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(TotalCurvMax, FlowInfo); }

  /// @brief Write HillslopeLength_HFR values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_HFR values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_HillslopeLength_HFR(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_HFR, FlowInfo); }

  /// @brief Write HillslopeLength_Binned values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_Binned values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_HillslopeLength_Binned(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_Binned, FlowInfo); }

  /// @brief Write HillslopeLength_Spline values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_Spline values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_HillslopeLength_Spline(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_Spline, FlowInfo); }

  /// @brief Write HillslopeLength_Density values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_Density values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_HillslopeLength_Density(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_Density, FlowInfo); }

  /// @brief Write FlowLength values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of FlowLength values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_FlowLength(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(FlowLength, FlowInfo); }

  /// @brief Write DrainageDensity values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of DrainageDensity values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_DrainageDensity(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(DrainageDensity, FlowInfo); }

  /// @brief Write CosmoErosionRate values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of CosmoErosionRate values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_CosmoErosionRate(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(CosmoErosionRate, FlowInfo); }

  /// @brief Write OtherErosionRate values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of OtherErosionRate values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_OtherErosionRate(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(OtherErosionRate, FlowInfo); }

  /// @brief Write CHTMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of CHTMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_CHTMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(CHTMean, FlowInfo); }

  /// @brief Write EStar values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of EStar values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_EStar(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(EStar, FlowInfo); }

  /// @brief Write RStar values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of RStar values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13
  LSDRaster write_RStar(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(RStar, FlowInfo); }

  /// @brief Method to merge a vector of LSDRaster basins generated using LSDBasin into
  /// a single LSDRaster for visualisation.
  ///
  /// @details 50% less computationally expesnive than the old method, but still very
  /// inefficient. Does not test for overlaps in the data, will simply overwrite
  /// so that the last value to occupy a cell will be written.
  ///
  /// @param Basins vector of LSDRasters of basins.
  /// @return LSDRaster of all the merged basins.
  /// @author SWDG
  /// @date 07/12/14
  LSDRaster Merge_Basins(vector<LSDRaster> Basins);


  /// @brief Count the number of each unique lithology value contain in the basin from a topographic raster
  /// take a lithologic raster and a topographic raster in argument
  ///
  /// @author BG
  /// @date 17/09/17
  map<int,int> count_unique_values_from_litho_raster(LSDIndexRaster& litho, LSDFlowInfo& topo);


  /// @brief merge and contour the perimeter from a vector of adjacent basins
  /// @detail WARNING There may be 1 pixel-size holes in the perimeter.
  /// @param vector of LSDBasin objects
  /// @author BG
  /// @date 10/10/17
  vector<int> merge_perimeter_nodes_adjacent_basins(vector<LSDBasin> budgerigar, LSDFlowInfo& flowpy);

  /// @brief Compare metrics inside/out of a basin for a given vector of nodes.
  /// @detail Ill detail when it will be done later
  /// @param rasterTemplate and the vector of node to test
  /// @author BG
  /// @date 23/12/17
  map<string,float> get_metrics_both_side_divide(LSDRaster& rasterTemplate, LSDFlowInfo& flowpy, vector<int>& nodes_to_test, map<int,bool>& raster_node_basin);

  /// @brief apply a square window around each perimeter nodes and extract statistics on each sides of the basin.
  /// @detail Ill detail when it will be done later
  /// @param rasterTemplate and the vector of node to test
  /// @author BG
  /// @date 23/12/17
  void square_window_stat_drainage_divide(LSDRaster& rasterTemplate, LSDFlowInfo& flowpy, int size_window);

  /// write the csv file corresponding to the previously calculated windowed stTS
  /// @detail Ill detail when it will be done later
  /// @param
  /// @author BG
  /// @date 23/12/17
  void write_windowed_stats_around_drainage_divide_csv(string full_name, LSDFlowInfo& flowpy);

  /// @brief Preprocess the Drainage Divide tool driver required info
  /// @detail Set the perimeter and set a map containing the corresponding x,y ...
  /// @detail TODO: add distance from origin and other global parameters
  /// @param FlowInfo object corresponding to the original raster where the Basin has been calculated
  /// @author BG
  /// @date 26/12/2017
  void preprocess_DD_metrics(LSDFlowInfo flowpy);

  void organise_perimeter(LSDFlowInfo& flowpy);

  void clean_perimeter(LSDFlowInfo& flowpy);


  /// @brief Write a csv file with X,Y,Z
  /// @param Name of the file and flowinfo object
  /// @author BG
  /// @date true
  void write_elevation_csv(string output_name, LSDFlowInfo& flowpy, LSDRaster& filled);

  /// @brief Write the channel network for this basin to a csv file
  /// @param csv_name
  /// @param FlowInfo
  /// @author FJC
  /// @date 18/08/18
  void write_channel_network(string csv_name, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork);

  /// @brief Get the mean data within a basin from a csv file
  /// @param FlowInfo
  /// @param CSV LSDSpatialCSVReader with the csv data
  /// @param column_name name of the column with the data
  /// @author FJC
  /// @date 29/09/18
  float get_basin_mean_from_csv(LSDFlowInfo& FlowInfo, LSDSpatialCSVReader& CSV, string column_name);

  protected:

  //These instance variables are set at initialisation

  ///Number of rows.
  int NRows;
  ///Number of columns.
  int NCols;
  ///Minimum X coordinate.
  float XMinimum;
  ///Minimum Y coordinate.
  float YMinimum;
  ///Data resolution.
  float DataResolution;
  ///No data value.
  int NoDataValue;
  ///A map of strings for holding georeferencing information
  map<string,string> GeoReferencingStrings;

  /// Junction Number of the basin, serves as a unique ID of a basin.
  int Junction;
  ///Vector of all nodes in basin.
  vector<int> BasinNodes;
  /// Map of nodes in the basin - faster to check if a node is in the basin due to a binary search tree for large maps
  map<int,int> nodes_of_basins;
  /// Number of DEM cells.
  int NumberOfCells;
  /// Area in spatial units of the basin.
  float Area;
  /// Boolean to show if the basin is beheaded or not.
  bool Beheaded;
  /// i index of the outlet pixel
  int Outlet_i;
  /// j index of the outlet pixel
  int Outlet_j;
  /// The Strahler order of the basin
  int BasinOrder;
  ///The i index of the centroid of the basin
  int Centroid_i;
  ///The j index of the centroid of the basin
  int Centroid_j;


  //These instance variables are used to store calculated basin parameters

  /// Mean basin slope.
  float SlopeMean;
  /// Mean basin elevation.
  float ElevationMean;
  /// Mean basin aspect.
  float AspectMean;
  /// Mean basin relief.
  float ReliefMean;
  /// Mean basin planform curvature.
  float PlanCurvMean;
  /// Mean basin profile curvature.
  float ProfileCurvMean;
  /// Mean basin total curvature.
  float TotalCurvMean;
  /// Max basin planform curvature.
  float PlanCurvMax;
  /// Max basin profile curvature.
  float ProfileCurvMax;
  /// Max basin total curvature.
  float TotalCurvMax;
  /// Mean hillslope length from hilltop flow routing.
  float HillslopeLength_HFR;
  /// Mean hillslope length from binned boomerang plots.
  float HillslopeLength_Binned;
  /// Mean hillslope length from spline curves on boomerang plots.
  float HillslopeLength_Spline;
  /// Mean hillslope length from drainage density.
  float HillslopeLength_Density;
  /// Basin flowlength.
  float FlowLength;
  /// Basin drainage density.
  float DrainageDensity;
  /// Basin Perimeter's j index.
  vector<int> Perimeter_i;
  /// Basin Perimeter's j index.
  vector<int> Perimeter_j;
  /// Basin Perimeter's node index
  vector<int> Perimeter_nodes;
  /// Basin Perimeter's node index, sorted by followed order
  vector<int> Perimeter_nodes_sorted;
  /// corresponding map giving an index to the sorted perimeter. Mostly for testing and debugging purposes.
  map<int,int> Perimeter_nodes_sorted_id;
  /// increase the speed of checking whether a node is perimeter or not compare to find in a vector
  map<int,int> Perimeter_nodes_map;
  /// Cosmo erosion rate.
  float CosmoErosionRate;
  /// Other erosion rate.
  float OtherErosionRate;
  /// Mean basin hilltop curvature.
  float CHTMean;
  /// Basin E* value.
  float EStar;
  /// Basin R* value.
  float RStar;
  /// Number of hilltop pixels in a basin.
  float HilltopPx;
  /// Fraction of basin mapped bedrock
  float BedrockFraction;
  /// AGB density
  float Biomass;
  // Alternative index (e.g. lithology)
  int AlternativeIndex;

  // Stuffs for DD purposes
  bool DD_preprocessed;
  // map of stats around the drainage divide
  map<int, map<string, float> > stats_around_perimeter_window;
  // map of distance from the origin of the perimeter, key is the node index
  map<int,float> map_of_dist_perim;
  // map of the perimeter location and metrics information per node
  map<int,vector<float> > DD_map;
  // map of xy location for each basin nodes. It fasten the process even if a bit memory-consuming.
  map<int,vector<float> > BasinNodesMapOfXY;

  private:
  void create();
  void create(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet);
  void create(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet, int extended_baslevel_nodes_below_junction);

};


/// @brief A derived class that is used to compute erosion rates based on
///  concentrations of in-situ cosmogenic nuclides such as 10Be and 26Al
class LSDCosmoBasin: public LSDBasin
{
  public:
    /// @brief Default constructor method used to create a  cosmo basin object.
    /// @param Junction outlet junction of the basin to be constructed.
    /// @param FlowInfo LSDFlowInfo object.
    /// @param ChanNet Channel network object.
    /// @param N10 concentration of 10Be in basin (atoms/g).
    /// @param del_N10 analytical uncertainty of 10Be in basin (atoms/g).
    /// @param N26 concentration of 10Be in basin (atoms/g).
    /// @param del_N26 analytical uncertainty of 10Be in basin (atoms/g).
    /// @author SMM
    /// @date 27/12/14
    LSDCosmoBasin(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet,
             double N10, double del_N10, double N26, double del_N26)
                   { create(Junction, FlowInfo, ChanNet, N10,del_N10,
                               N26, del_N26); }

    /// @brief This function populates the scaling vectors that are used to set
    ///  the production scaling, topographic shielding and snow shielding
    ///  for specific nodes
    ///
    /// @details This is the default scaling that assumes no snow shielding (i.e.,
    ///  snow_shielding == 1). The assumptions are Stone scaling with Fsp = 1
    ///  and that the DEM has a WGS84 ellipsoid
    /// @param FlowInfo the LSDFlowInfo object
    /// @param Elevation_Data the DEM, an LSDRaster object. IMPORTANT!! This needs
    ///  to contain georeferencing information for this to work!!!
    /// @param Topo_Shield an LSDRaster with the topographic shielding
    /// @param path_to_atmospheric_data THis is a path to binary NCEP data.
    /// @author SMM
    /// @date 22/12/2014
    void populate_scaling_vectors(LSDFlowInfo& FlowInfo, LSDRaster& Elevation_Data,
                                  LSDRaster& Topo_Shield, string path_to_atmospheric_data);

    /// @brief This function populates the scaling vectors that are used to set
    ///  the production scaling, topographic shielding and snow shielding
    ///  for specific nodes.
    ///
    /// @details The snow shiedling in inculded in this function
    ///  The assumptions are Stone scaling with Fsp = 1
    ///  and that the DEM has a WGS84 ellipsoid
    /// @param FlowInfo the LSDFlowInfo object
    /// @param Elevation_Data the DEM, an LSDRaster object. IMPORTANT!! This needs
    ///  to contain georeferencing information for this to work!!!
    /// @param Topo_Shield an LSDRaster with the topographic shielding
    /// @param Snow_shield an LSDRaster containing the snow shielding.
    /// @param path_to_atmospheric_data THis is a path to binary NCEP data.
    /// @author SMM
    /// @date 22/12/2014
    void populate_scaling_vectors(LSDFlowInfo& FlowInfo, LSDRaster& Elevation_Data,
                                  LSDRaster& Topo_Shield, LSDRaster& Snow_shield,
                                  string path_to_atmospheric_data);

    /// @brief This function sets the basin vectors for snow and shelf sheilding
    ///  if the effective depth method is to be used.
    ///
    /// @details This calculates snow and self shielding by integrating production
    ///  over the effective depth interval from which sediment comes. Snow
    ///  cover is represented by an effective depth. This method uses
    ///  the full muon profile to predict the CRN concetration
    /// @param FlowInfo the LSDFlowInfo object
    /// @param snow_eff_depth an LSDRaster holding the effective depths (g/cm^2)
    ///  of the snow for each pixel
    /// @param self_eff_depth an LSDRaster holding the effective depths (g/cm^2)
    ///  of self shielding for each pixel. This can implement landsliding,
    ///  since landsliding is effectively a self shielding process
    ///  whereby sediment from depth is exhumed.
    /// @author SMM
    /// @date 23/02/2015
    void populate_snow_and_self_eff_depth_vectors(LSDFlowInfo& FlowInfo,
                                              LSDRaster& snow_eff_depth,
                                              LSDRaster& self_eff_depth);

    /// @brief This function sets the basin vectors for snow and shelf sheilding
    ///  if the effective depth method is to be used.
    ///
    /// @details Overloaded function where snow shielding is simply a double
    /// @param FlowInfo the LSDFlowInfo object
    /// @param snow_eff_depth the shelding depth from snow (over the entire DEM)
    ///  in g/cm^2
    /// @param self_eff_depth an LSDRaster holding the effective depths (g/cm^2)
    ///  of self shielding for each pixel. This can implement landsliding,
    ///  since landsliding is effectively a self shielding process
    ///  whereby sediment from depth is exhumed.
    /// @author SMM
    /// @date 23/02/2015
    void populate_snow_and_self_eff_depth_vectors(LSDFlowInfo& FlowInfo,
                                              double snow_eff_depth,
                                              LSDRaster& self_eff_depth);

    /// @brief This function sets the basin vectors for snow and shelf sheilding
    ///  if the effective depth method is to be used.
    ///
    /// @details Overloaded function where self shielding is simply a double
    /// @param FlowInfo the LSDFlowInfo object
    /// @param snow_eff_depth an LSDRaster holding the effective depths (g/cm^2)
    ///  of the snow for each pixel
    /// @param self_eff_depth the self shelding depth (over the entire DEM)
    ///  in g/cm^2
    ///  whereby sediment from depth is exhumed.
    /// @author SMM
    /// @date 23/02/2015
    void populate_snow_and_self_eff_depth_vectors(LSDFlowInfo& FlowInfo,
                                              LSDRaster& snow_eff_depth,
                                              double self_eff_depth);

    /// @brief This function sets the basin vectors for snow and shelf sheilding
    ///  if the effective depth method is to be used.
    ///
    /// @details Overloaded function where both self and snow shielding
    ///  are simply doubles
    /// @param snow_eff_depth the shelding depth from snow (over the entire DEM)
    ///  in g/cm^2
    /// @param self_eff_depth the self shelding depth (over the entire DEM)
    ///  in g/cm^2
    ///  whereby sediment from depth is exhumed.
    /// @author SMM
    /// @date 23/02/2015
    void populate_snow_and_self_eff_depth_vectors(double snow_eff_depth,
                                                  double self_eff_depth);

    /// @brief this resets the snow and self shielding effective vectors.
    ///
    /// @details resets now and self sheidling vectors to empty vectors
    ///  this allows calucaltions to be made with no snow and self sheilding
    ///  after a self sheilding calculation is made
    /// @author SMM
    /// @date 21/02/2015
    void reset_snow_and_self_eff_depth_vectors();

    /// @brief This is a utility function that populates the atmospheric pressure
    ///  vector. It uses the CRONUS calculator scheme.
    ///
    /// @details The function is mainly used for bug checking
    /// @param FlowInfo the LSDFlowInfo object
    /// @param Elevation_Data the DEM, an LSDRaster object. IMPORTANT!! This needs
    ///  to contain georeferencing information for this to work!!!
    /// @param path_to_atmospheric_data THis is a path to binary NCEP data.
    /// @author SMM
    /// @date 28/01/2015
    void get_atmospheric_pressure(LSDFlowInfo& FlowInfo, LSDRaster& Elevation_Data,
                                  string path_to_atmospheric_data);

    /// @brief This function wraps the erosion rate calculator, and returns
    ///  both the erosion rate as well as the uncertainties
    /// @param Nuclide_conc Concetration of the nuclide
    /// @param Nuclide a string denoting the name of the nuclide (at the moment
    ///  options are 10Be and 26Al)
    /// @param Nuclide_conc_err The instrument error in the nuclide concentration
    /// @param prod_uncert_factor This is a fraction of the total uncertainty
    ///  for the production rates. It is a lumped parameter that can be used
    ///  for just production, or for snow, topo and porduction uncertainty
    /// @param Muon_scaling string that gives the muon scaling scheme
    ///  options are Schaller, Granger and Braucher
    /// @return  a vector of both the erosion rates and the uncertainties of the sample
    /// @author SMM
    /// @date 01/02/2015
    vector<double> full_CRN_erosion_analysis(double Nuclide_conc, string Nuclide,
                            double Nuclide_conc_err, double prod_uncert_factor,
                            string Muon_scaling);

    /// @brief This function wraps the erosion rate calculator, and returns
    ///  both the erosion rate as well as the uncertainties  ^
    /// @param known_eff_erosion a raster containing known effective erosion rates (g/cm2/yr)
    /// @param FlowInfo the LSDFlowInfo object
    /// @param Nuclide_conc Concetration of the nuclide
    /// @param Nuclide a string denoting the name of the nuclide (at the moment
    ///  options are 10Be and 26Al)
    /// @param Nuclide_conc_err The instrument error in the nuclide concentration
    /// @param prod_uncert_factor This is a fraction of the total uncertainty
    ///  for the production rates. It is a lumped parameter that can be used
    ///  for just production, or for snow, topo and porduction uncertainty
    /// @param Muon_scaling string that gives the muon scaling scheme
    ///  options are Schaller, Granger and Braucher
    /// @return  a vector of both the erosion rates and the uncertainties of the sample
    /// @author SMM
    /// @date 11/02/2016
    vector<double> full_CRN_erosion_analysis_nested(LSDRaster& known_eff_erosion,
                              LSDFlowInfo& FlowInfo, double Nuclide_conc, string Nuclide,
                              double Nuclide_conc_err, double prod_uncert_factor,
                              string Muon_scaling);

    /// @brief this uses Newton Raphson iteration to retrieve the erosion rate
    ///  from a basin given a nuclide concentration
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param average_production This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return The effective erosion rate in g/cm^-2/yr
    /// @author SMM
    /// @date 03/01/2015
    double predict_CRN_erosion(double Nuclide_conc, string Nuclide,
                               double prod_uncert_factor,string Muon_scaling,
                               double& production_uncertainty,
                               double& average_production,
                               bool is_production_uncertainty_plus_on,
                               bool is_production_uncertainty_minus_on);

    /// @brief this uses Newton Raphson iteration to retrieve the erosion rate
    ///  from a basin given a nuclide concentration.  This is the nested version.
    ///
    /// @details The nesting can be done by setting the erosion rate of pixels within
    ///  the basin
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param average_production This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return The effective erosion rate in g/cm^-2/yr
    /// @author SMM
    /// @date 03/01/2015
    double predict_CRN_erosion_nested(double Nuclide_conc, string Nuclide,
                                          double prod_uncert_factor,
                                          string Muon_scaling,
                                          double& production_uncertainty,
                                          double& average_production,
                                          bool is_production_uncertainty_plus_on,
                                          bool is_production_uncertainty_minus_on,
                                          LSDRaster& eff_erosion_raster,
                                          LSDFlowInfo& FlowInfo);


    /// @brief this predicts the mean concentration of a nuclide within
    /// a basin
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger
    /// @param data_from_outlet_only boolean that is true of you want
    ///  concentration calculated from the outlet only.
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param production_rate This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 22/12/2014
    double predict_mean_CRN_conc(double eff_erosion_rate, string Nuclide,
                                 double prod_uncert_factor,
                                 string Muon_scaling, bool data_from_outlet_only,
                                 double& production_uncertainty,
                                 double& production_rate,
                                 bool is_production_uncertainty_plus_on,
                                 bool is_production_uncertainty_minus_on);

    /// @brief this predicts the mean concentration of a nuclide within
    ///  a basin. It does a full analyitical solution to account for
    ///  snow and self sheilding
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger
    /// @param data_from_outlet_only boolean that is true of you want
    ///  concentration calculated from the outlet only.
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param production_rate This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 22/12/2014
    double predict_mean_CRN_conc_with_snow_and_self(double eff_erosion_rate, string Nuclide,
                                 double prod_uncert_factor,
                                 string Muon_scaling, bool data_from_outlet_only,
                                 double& production_uncertainty,
                                 double& production_rate,
                                 bool is_production_uncertainty_plus_on,
                                 bool is_production_uncertainty_minus_on);

    /// @brief this predicts the mean concentration of a nuclide within
    ///  a basin. It does a full analyitical solution to account for
    ///  snow and self sheilding. It uses a raster of known erosion rates
    ///  to predict the concentration from the basin: it is primarily used for
    ///  nesting.
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param known_effective_erosion a raster of known effective erosion rates (g/cm^2/yr)
    /// @param FlowInfo the LSDFlowInfo object
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param production_rate This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 06/02/2016
    double predict_mean_CRN_conc_with_snow_and_self_nested(double eff_erosion_rate,
                                            LSDRaster& known_effective_erosion,
                                            LSDFlowInfo& FlowInfo,
                                            string Nuclide,
                                            double prod_uncert_factor, string Muon_scaling,
                                            double& production_uncertainty,
                                            double& average_production,
                                            bool is_production_uncertainty_plus_on,
                                            bool is_production_uncertainty_minus_on);

    /// @brief this predicts the mean concentration of a nuclide within
    ///  a basin, using the production scaling of the centroid
    ///  It replicates the technique used by many authors. This function
    ///  is mainly here to show how far off this method is compared to
    ///  the pixel-by-pixel production scaling
    /// @param eff_erosion rate The erosion rate in g/cm^2/yr
    /// @param Nuclide a string with the nuclide name. At the moment the options are:
    ///   Be10
    ///   Al26
    ///  These are case sensitive
    /// @param prod_uncert_factor production uncertainty factor is a multiplier that sets the production
    ///  certainty. If it is 1.1, there is 10% production rate uncertainty, or
    ///  if it is 0.9 there is -10% unvertainty. The reason why it is implemented
    ///  like this is that this allows gaussian error propigation.
    /// @param Muon_scaling a string that gives the muon scaling scheme.
    ///  options are Schaller, Braucher and Granger.
    /// @param FlowInfo an LSDFlowInfo object
    /// @param production_uncertainty this gives the uncertainty in the production
    ///  rates based on the production_uncert_factor; it is used in gaussian
    ///  error propigation. The parameter is replaced within the function.
    /// @param production_rate This gives the production rate average for the
    ///  basin. It can be used for uncertainty analyis: if the scaling is
    ///  changed the change in this production rate can be used to construct
    ///  the gaussian error propigation terms
    /// @param is_production_uncertainty_plus_on a boolean that is true if the
    ///  production rate uncertainty (+) is switched on
    /// @param is_production_uncertainty_minus_on a boolean that is true if the
    ///  production rate uncertainty (-) is switched on. If the + switch is
    ///  true this parameter defauts to false.
    /// @return the concentration of the nuclide averaged across the DEM.
    /// @author SMM
    /// @date 28/01/2015
    double predict_mean_CRN_conc_centroid(double eff_erosion_rate, string Nuclide,
                                    double prod_uncert_factor, string Muon_scaling,
                                    LSDFlowInfo& FlowInfo, double& production_uncertainty,
                                    double& production_rate,
                                    bool is_production_uncertainty_plus_on,
                                    bool is_production_uncertainty_minus_on);

    /// @breif A function for testing if a known erosion rate raster contains any unknowns within a basin.
    /// @param known_erates a raster of known erosion rates
    /// @param FlowInfo a flow info object
    /// @return a boolean that is true if there are unknown erosion rates and false if not
    /// @author SMM
    /// @date 10/02/2016
    bool are_there_unknown_erosion_rates_in_basin(LSDRaster& known_erates,LSDFlowInfo& FlowInfo);

    /// @brief Prints a csv with information about the nodes in a basin that
    ///  relate to cosmogenic paramters.
    ///
    /// @details the csv file out has the format:
    ///  fID,Easting,Northing,Latitude,Longitude,Elevation,Pressure,...
    ///  TopoShield,Production_scaling,Snowshield
    ///  IMPORTANT: In the normal basin, the outlet is the penultamite node
    ///  in the channel drainaing from the junction. This is so that the basin
    ///  will drain to a node just before the next stram order. However in the
    ///  cosmo basin this is not the case because the cosmo was selected from a
    ///  specific location so we use the closest junction.
    /// @param path_to_file the path to the outfile. Needs a / at the end.
    /// @param filename. This does not have an extension because .csv will be added
    /// @param FlowInfo the LSDFlowInfo object
    /// @param Elevation_Data the DEM, an LSDRaster object. IMPORTANT!! This needs
    ///  to contain georeferencing information for this to work!!!
    /// @param Topo_Shield an LSDRaster with the topographic shielding
    /// @param path_to_atmospheric_data THis is a path to binary NCEP data.
    /// @return the concentration of the nuclide averaged across the DEM
    /// @author SMM
    /// @date 27/12/2014
    void print_particle_csv(string path_to_file, string filename,
                            LSDFlowInfo& FlowInfo, LSDRaster& Elevation_Data,
                            LSDRaster& T_Shield,
                            string path_to_atmospheric_data);


    /// @brief This function gets effective pressure for a basin by averageing
    ///  production rates and then calculating the 'apparent' pressure that can
    ///  produce that production rate. It also returns a number of valriables
    ///  that can be plugged into existing calculators to compare results
    ///  from the LSDCosmo calulcator against other methods
    /// @param Elevation the DEM, an LSDRaster object. IMPORTANT!! This needs
    ///  to contain georeferencing information for this to work!!!
    /// @param FlowInfo the LSDFlowInfo object
    /// @param path_to_atmospheric_data THis is a path to binary NCEP data.
    /// @return a vector of values used in open-source calculators
    ///  vector[0] = AverageProd;
    ///  vector[1] = AverageTopo;
    ///  vector[2] = AverageSelf;
    ///  vector[3] = AverageSnow;
    ///  vector[4] = AverageCombined;
    ///  vector[5] = lat_outlet;
    ///  vector[6] = outlet_pressure;
    ///  vector[7] = outlet_eff_pressure;
    ///  vector[8] = lat_centroid;
    ///  vector[9] = centroid_pressure;
    ///  vector[10] = centroid_eff_pressure;
    /// @author SMM
    /// @date 11/03/2015
    vector<double> calculate_effective_pressures_for_calculators(LSDRaster& Elevation,
                        LSDFlowInfo& FlowInfo, string path_to_atmospheric_data);

    /// @brief This function gets effective pressure for a basin by averageing
    ///  production rates and then calculating the 'apparent' pressure that can
    ///  produce that production rate. It also returns a number of valriables
    ///  that can be plugged into existing calculators to compare results
    ///  from the LSDCosmo calulcator against other methods
    ///  This version allows you to input a know erosion rate raster. Production and other
    ///  scalings are masked anywhere the erosion rates are known.
    /// @param Elevation the DEM, an LSDRaster object. IMPORTANT!! This needs
    ///  to contain georeferencing information for this to work!!!
    /// @param FlowInfo the LSDFlowInfo object
    /// @param path_to_atmospheric_data THis is a path to binary NCEP data.
    /// @param known_eff_erosion An LSDRaster containing the known effective erosion rates
    /// @return a vector of values used in open-source calculators
    ///  vector[0] = AverageProd;
    ///  vector[1] = AverageTopo;
    ///  vector[2] = AverageSelf;
    ///  vector[3] = AverageSnow;
    ///  vector[4] = AverageCombined;
    ///  vector[5] = lat_outlet;
    ///  vector[6] = outlet_pressure;
    ///  vector[7] = outlet_eff_pressure;
    ///  vector[8] = lat_centroid;
    ///  vector[9] = centroid_pressure;
    ///  vector[10] = centroid_eff_pressure;
    /// @author SMM
    /// @date 11/03/2015
    vector<double> calculate_effective_pressures_for_calculators_nested(LSDRaster& Elevation,
                        LSDFlowInfo& FlowInfo, string path_to_atmospheric_data,
                        LSDRaster& known_eff_erosion);

    /// @brief This function prints the production scaling as well as the
    ///  combined scaling and combined shielding on a pixel by pixel basis to
    ///  rasters.
    ///
    /// @details The production rster _PROD has the latitude and elevation adjusted
    ///  scaling using the Stone scheme
    ///   The combined shielding raster _CSHIELD is the product of topographic,
    ///   self and snow shielding, assuming that the snow and self shielding are
    ///   from spallation only
    ///   The combined scaling raster _CSCALE is the combined sheidling values mulitplied by
    ///   the production scaling.
    /// @param filename the name of the file: this is generally the DEM file with
    ///  an extension for the basin name
    /// @param FlowInfo The LSDFlowInfo object
    /// @author SMM
    /// @date 24/03/2015
    void print_scaling_and_shielding_rasters(string filename,LSDFlowInfo& FlowInfo);

    /// @brief returns the combned scaling raster, which has the product of
    ///  the topographic, snow and self shielding, multiplied by production scaling,
    ///  for each pixel.
    ///
    /// @details The snow and self shielding factors are calculate assuming spallation only!
    /// @param filename the name of the file: this is generally the DEM file with
    ///  an extension for the basin name
    /// @param FlowInfo The LSDFlowInfo object
    /// @author SMM
    /// @date 6/02/2016
    LSDRaster get_combined_scaling_raster(string filename,LSDFlowInfo& FlowInfo);


    /// @brief This function is run after the CRN erosion rate have been found
    ///  it prints a raster contiannig the spatially distributed concentration
    ///  of the nuclide predicted for the estimated erosion rate
    /// @param filename the name of the file: this is generally the DEM file with
    ///  an extension for the basin name
    /// @param eff_erosion_rate: the erosion rate calculated for the basin
    /// @param Nuclide a string containing either Al26 or Be10
    /// @param Muon_scaling a string with either Braucher, Schaller or Granger
    /// @param FlowInfo The LSDFlowInfo object
    /// @author SMM
    /// @date 24/03/2015
    void print_CRN_conc_raster(string filename, double eff_erosion_rate, string Nuclide,
                               string Muon_scaling, LSDFlowInfo& FlowInfo);

  protected:
    /// The measured 10Be concentration
    double measured_N_10Be;

    /// The measured 26Al concentration
    double measured_N_26Al;

    /// The measured uncertainty in the 10Be concentration
    double delN_10Be;

    /// The measured uncertainty in the 10Be concentration
    double delN_26Al;

    /// A vector holding the elevations of the data within the basin
    vector<double> snow_shielding;

    /// A vector holding the topographic shielding of nodes within the basin
    vector<double> topographic_shielding;

    /// A vector holding the production scaling of nodes within the basin
    vector<double> production_scaling;

    /// A vector containing self shielding values. This could be due to
    /// a finite thickness of material assumed to be eroding. It can be used
    /// to calculate landsliding adjustment of cosmogenics.
    vector<double> self_shielding;

    /// a vector holding the CRNparticles
    vector<LSDCRNParticle> CRN_particle_vec;

    /// a vector holding the atmospheric pressure
    vector<double> atmospheric_pressure;

    /// this is a vector for holding the effective depth of slf shelding
    /// used in the full muon based self shielding
    /// in g/cm^2
    vector<double> self_shield_eff_depth;

    /// This holds an effective depth of snow for a basin
    /// in g/cm^2
    vector<double> snow_shield_eff_depth;

  private:
    void create(int JunctionNumber, LSDFlowInfo& FlowInfo,
                           LSDJunctionNetwork& ChanNet,
                           double N10Be, double delN10Be,
                           double N26Al, double delN26Al);


};


#endif
