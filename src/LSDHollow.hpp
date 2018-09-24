//beginning of the LSDHollow object

#include <vector>
#include <string>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDHollow_H
#define LSDHollow_H

///@brief Object to store information about hollows and generate hollow average metrics.
class LSDHollow
{

  public:
  
  /// @brief Default constructor method used to create a hollow object.
  /// @param Junction outlet junction of the hollow to be constructed.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param ChanNet Channel network object.
  /// @author SWDG
  /// @date 19/2/14
  LSDHollow(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet)
											{ create(Junction, FlowInfo, ChanNet); }


  /// @return Number of rows as an integer.
	int get_NRows() const				{ return NRows; }
	/// @return Number of columns as an integer.
  int get_NCols() const				{ return NCols; }
  /// @return Minimum X coordinate as an integer.
	float get_XMinimum() const			{ return XMinimum; }
	/// @return Minimum Y coordinate as an integer.
	float get_YMinimum() const			{ return YMinimum; }
	/// @return Data resolution as an integer.
	float get_DataResolution() const	{ return DataResolution; }
	/// @return No Data Value as an integer.
	int get_NoDataValue() const			{ return NoDataValue; }
	
	/// @return Nodes of the hollow.
	vector<int> get_HollowNodes() const { return HollowNodes; }
	
	/// @return Number of cells in the hollow.
	int get_NumberOfCells() const { return NumberOfCells; }
	
  /// @return Area of hollow in spatial units.
  float get_Area() const { return Area; }
  
  /// @return Junction Number of the hollow.
  int get_Junction() const { return Junction; }
    
  /// @return Boolean value of whether a hollow is beheaded or not.										
  bool get_Beheaded() const { return Beheaded; }
  
  /// @return i index of outlet pixel.
  int get_Outlet_i() const { return Outlet_i; }

  /// @return j index of outlet pixel.
  int get_Outlet_j() const { return Outlet_j; }
  
  ///@return i index of Centroid pixel.
  int get_Centroid_i() const { return Centroid_i; }

  ///@return j index of Centroid pixel.
  int get_Centroid_j() const { return Centroid_j; }
  
  //Getters of hollow parameters
  
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
  /// @return Hollow Perimeter's i index.
  vector<int> get_Perimeter_i() const { return Perimeter_i; }
  /// @return Hollow Perimeter's j index.
  vector<int> get_Perimeter_j() const { return Perimeter_j; }
  /// @return Mean hilltop curvature.
  float get_CHTMean() const { return CHTMean; }
  /// @return Soil production rate.
  float get_SoilProduction() const { return SoilProduction; }
  /// @return Basal age of hollow.
  float get_BasalAge() const { return BasalAge; }
  /// @return Hollow Width.
  float get_Width() const { return Width; }
  /// @return Downslope length.
  float get_DownslopeLength() const { return DownslopeLength; }
  /// @return Long Profile Length. 
  float get_LongProfileLength() const {return LongProfileLength; }
  
  /// @brief Calculate the mean value of an LSDRaster which falls inside a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the mean of.
  /// @return Mean value.
  /// @author SWDG
  /// @date 14/2/14
  float CalculateHollowMean(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Calculate the max value of an LSDRaster which falls inside a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the max of.
  /// @return Max value.
  /// @author SWDG
  /// @date 14/2/14
  float CalculateHollowMax(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Return the LSDRaster data from within a hollow as an array.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to extract.
  /// @return Array of values.
  /// @author SWDG
  /// @date 20/2/14  
  Array2D<float> get_Raster_Data_For_Hollow(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Return the LSDIndexRaster data from within a hollow as an array.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to extract.
  /// @return Array of values.
  /// @author SWDG
  /// @date 20/2/14  
  Array2D<int> get_Raster_Data_For_Hollow(LSDFlowInfo& FlowInfo, LSDIndexRaster Data);
  
  /// @brief Set the mean slope of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param Slope Values to find the mean of.
  /// @author SWDG
  /// @date 14/2/14
  void set_SlopeMean(LSDFlowInfo& FlowInfo, LSDRaster Slope){ SlopeMean = CalculateHollowMean(FlowInfo, Slope); }

  /// @brief Set the mean Elevation of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param Elevation Values to find the mean of.
  /// @author SWDG
  /// @date 14/2/14
  void set_ElevationMean(LSDFlowInfo& FlowInfo, LSDRaster Elevation) { ElevationMean = CalculateHollowMean(FlowInfo, Elevation); }

  /// @brief Set the mean Relief of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param Relief Values to find the mean of.
  /// @author SWDG
  /// @date 14/2/14
  void set_ReliefMean(LSDFlowInfo& FlowInfo, LSDRaster Relief) { ReliefMean = CalculateHollowMean(FlowInfo, Relief); }

  /// @brief Set the mean PlanCurve of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param PlanCurv Values to find the mean of.
  /// @author SWDG
  /// @date 14/2/14
  void set_PlanCurvMean(LSDFlowInfo& FlowInfo, LSDRaster PlanCurv) { PlanCurvMean = CalculateHollowMean(FlowInfo, PlanCurv); }

  /// @brief Set the mean ProfCurv of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param ProfileCurv Values to find the mean of.
  /// @author SWDG
  /// @date 14/2/14
  void set_ProfileCurvMean(LSDFlowInfo& FlowInfo, LSDRaster ProfileCurv) { ProfileCurvMean = CalculateHollowMean(FlowInfo, ProfileCurv); }

  /// @brief Set the mean TotalCurv of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param TotalCurv Values to find the mean of.
  /// @author SWDG
  /// @date 14/2/14
  void set_TotalCurvMean(LSDFlowInfo& FlowInfo, LSDRaster TotalCurv) { TotalCurvMean = CalculateHollowMean(FlowInfo, TotalCurv); }

  /// @brief Set the max PlanCurve of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param PlanCurv Values to find the max of.
  /// @author SWDG
  /// @date 14/2/14
  void set_PlanCurvMax(LSDFlowInfo& FlowInfo, LSDRaster PlanCurv) { PlanCurvMax = CalculateHollowMax(FlowInfo, PlanCurv); }

  /// @brief Set the max ProfCurv of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param ProfileCurv Values to find the max of.
  /// @author SWDG
  /// @date 14/2/14
  void set_ProfileCurvMax(LSDFlowInfo& FlowInfo, LSDRaster ProfileCurv) { ProfileCurvMax = CalculateHollowMax(FlowInfo, ProfileCurv); }

  /// @brief Set the max TotalCurv of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param TotalCurv Values to find the max of.
  /// @author SWDG
  /// @date 14/2/14
  void set_TotalCurvMax(LSDFlowInfo& FlowInfo, LSDRaster TotalCurv) { TotalCurvMax = CalculateHollowMax(FlowInfo, TotalCurv); }
 
  /// @brief Set the mean hilltop curvature of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param CHT Values to find the mean of.
  /// @author SWDG
  /// @date 19/2/14
  void set_CHTMean(LSDFlowInfo& FlowInfo, LSDRaster CHT) { CHTMean = CalculateHollowMean(FlowInfo, CHT); }
 
  /// @brief Set the mean aspect of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param Aspect Values to find the mean of.
  /// @author SWDG
  /// @date 19/2/14
  void set_AspectMean(LSDFlowInfo& FlowInfo, LSDRaster Aspect);
  
  /// @brief Set the perimeter pixels using a simple edge detection algorithm. 
  ///
  /// @details This is quite messy and will be improved soon.
  /// @param FlowInfo Flowinfo object.
  /// @author SWDG
  /// @date 19/2/14
  void set_Perimeter(LSDFlowInfo& FlowInfo);
  
  /// @brief Set the Soil production rate.
  /// @param SoilProdRate Soil production rate - No sanity check on this value.
  /// @author SWDG
  /// @date 19/02/14
  void set_SoilProduction(float SoilProdRate) { SoilProduction = SoilProdRate; }
  
  /// @brief Set the basal age.
  /// @param BasalAge Basal age - No sanity check on this value.
  /// @author SWDG
  /// @date 19/02/14
  void set_BasalAge(float basal_age) { BasalAge = basal_age; }
  
  
  
  /// @brief Set all of the hollow parameters with one call.
  ///
  /// @details Runs polyfit to get the elevation derivatives, so can be quite memory intensive. Method
  /// calls all the setters one by one, to populate all the hollow parameters. So a
  /// hollow can be created and all it's properties set with 2 calls. The BasalAge and SoilProduction have default 
  /// parameters of -9999 as these are rarely used variables.
  /// @param Elevation LSDRaster of filled elevation values.
  /// @param FlowInfo Flowinfo object.
  /// @param CHT LSDRaster of hilltop curvatures.
  /// @param Relief LSDRaster of the hilltop relief.
  /// @param window_radius Radius in spatial units for the polyft routine.
  /// @param SoilProduction Soil production rate.
  /// @param BasalAge Basal age of the hollow.
  /// @author SWDG
  /// @date 20/2/14
  void set_All_Parameters(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& CHT,
                                  LSDRaster& Relief, float window_radius,
                                  float SoilProduction = -9999, float BasalAge = -9999);
    
  /// @brief Cookie cut data from an LSDIndexRaster into the shape of the hollow.
  /// @param Data LSDIndexRaster data to be written.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of the data in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDIndexRaster write_raster_data_to_LSDIndexRaster(LSDIndexRaster Data, LSDFlowInfo FlowInfo);
  
  /// @brief Cookie cut data from an LSDRaster into the shape of the hollow.
  /// @param Data LSDRaster data to be written.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of the data in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_raster_data_to_LSDRaster(LSDRaster Data, LSDFlowInfo FlowInfo);
  
  /// @brief Write a real value to an LSDRaster in the shape of the hollow.
  /// @param Param real value to be written
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of the data in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_real_data_to_LSDRaster(float Param, LSDFlowInfo FlowInfo);
 
  /// @brief Write an integer value to an LSDIndexRaster in the shape of the hollow.
  /// @param Param integer value to be written
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of the data in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDIndexRaster write_integer_data_to_LSDIndexRaster(int Param, LSDFlowInfo FlowInfo);
  
  /// @brief Write Junction values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of Junction values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDIndexRaster write_Junction(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(Junction, FlowInfo); }
  
  /// @brief Write NumberOfCells values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of NumberOfCells values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDIndexRaster write_NumberOfCells(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(NumberOfCells, FlowInfo); }
    
  /// @brief Write Area values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of Area values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_Area(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(Area, FlowInfo); }
  
  /// @brief Write SlopeMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of SlopeMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_SlopeMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(SlopeMean, FlowInfo); }
  
  /// @brief Write ElevationMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ElevationMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_ElevationMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ElevationMean, FlowInfo); }
  
  /// @brief Write AspectMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of AspectMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_AspectMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(AspectMean, FlowInfo); }
  
  /// @brief Write ReliefMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ReliefMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_ReliefMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ReliefMean, FlowInfo); }
  
  /// @brief Write PlanCurvMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of PlanCurvMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_PlanCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(PlanCurvMean, FlowInfo); }
  
  /// @brief Write ProfileCurvMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ProfileCurvMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_ProfileCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ProfileCurvMean, FlowInfo); }
  
  /// @brief Write TotalCurvMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of TotalCurvMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_TotalCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(TotalCurvMean, FlowInfo); }
  
  /// @brief Write PlanCurvMax values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of PlanCurvMax values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_PlanCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(PlanCurvMax, FlowInfo); }
  
  /// @brief Write ProfileCurvMax values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ProfileCurvMax values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_ProfileCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ProfileCurvMax, FlowInfo); }
  
  /// @brief Write TotalCurvMax values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of TotalCurvMax values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_TotalCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(TotalCurvMax, FlowInfo); }
    
  /// @brief Write CHTMean values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of CHTMean values in the shape of the hollow.
  /// @author SWDG
  /// @date 14/2/14
  LSDRaster write_CHTMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(CHTMean, FlowInfo); }

  /// @brief Write Width values into the shape of the hollow.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of Width values in the shape of the hollow.
  /// @author SWDG
  /// @date 20/2/14
  LSDRaster write_Width(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(Width, FlowInfo); }
    
  /// @brief Calculate the width of a hollow.
  /// @param FlowInfo Flowinfo object.
  /// @param FlowDir D infinity flow directions.
  /// @author SWDG
  /// @date 19/02/14 
  void set_Width(LSDFlowInfo FlowInfo, Array2D<float> FlowDir);
  
  /// @brief Calculate the downslope length in the hollow, defined as the D8 flow routing distance
  /// from between the maximum elevation point in the hollow to the minimum elevation, or 
  /// to the edge of the hollow, whichever comes first.
  /// @param FlowInfo Flowinfo object.
  /// @param DEM An LSDRaster of elevations. 
  /// @author SWDG
  /// @date 20/02/14  
  void set_DownslopeLength(LSDFlowInfo FlowInfo, LSDRaster DEM);
   
  /// @brief Measure the long profile length in the hollow, defined as the maximum dimension 
  /// of the bounding box excluding diagonals.
  /// @param FlowInfo Flowinfo object. 
  /// @author SWDG
  /// @date 20/02/14  
  void set_LongProfileLength(LSDFlowInfo FlowInfo);
  
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
  /// Junction Number of the hollow, serves as a unique ID of a hollow.
  int Junction;
  ///Vector of all nodes in hollow.
  vector<int> HollowNodes;
  /// Number of DEM cells. 
  int NumberOfCells;
  /// Area in spatial units of the hollow.
  float Area;
  /// Boolean to show if the hollow is beheaded or not.
  bool Beheaded;
  /// i index of the outlet pixel
  int Outlet_i;
  /// j index of the outlet pixel
  int Outlet_j;
  ///The i index of the centroid of the hollow
  int Centroid_i;
  ///The j index of the centroid of the hollow
  int Centroid_j;
  
  
  //These instance variables are used to store calculated hollow parameters 
  
  /// Mean hollow slope.
  float SlopeMean;
  /// Mean hollow elevation.
  float ElevationMean;
  /// Mean hollow aspect.
  float AspectMean;
  /// Mean hollow relief.
  float ReliefMean;
  /// Mean hollow planform curvature.
  float PlanCurvMean;
  /// Mean hollow profile curvature.
  float ProfileCurvMean;
  /// Mean hollow total curvature.
  float TotalCurvMean;
  /// Max hollow planform curvature.
  float PlanCurvMax;
  /// Max hollow profile curvature.
  float ProfileCurvMax;
  /// Max hollow total curvature.
  float TotalCurvMax;
  /// Hollow Perimeter's j index.
  vector<int> Perimeter_i;
  /// Hollow Perimeter's j index.
  vector<int> Perimeter_j;
  /// Mean hollow hilltop curvature.
  float CHTMean;
  /// Soil production rate.
  float SoilProduction;
  /// Basal age of hollow.
  float BasalAge;
  /// Hollow width
  float Width;
  /// Downslope length
  float DownslopeLength;
  /// Long Profile Length
  float LongProfileLength;
  
  private:
	void create(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet);

};

#endif
