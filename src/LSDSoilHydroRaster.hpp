//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDSoilHydroRaster
// Land Surface Dynamics Soil and Hydrologic Raster
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object wraps some functions for manipulating derived soil and hydrology
//  functions that work with rasters
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2015 Simon M. Mudd 2013 5
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


/** @file LSDRaster.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 1.0.0
@brief Main analysis object to interface with other LSD objects.
@details This object contains a diverse range of geomophological
analysis routines which can be used in conjunction with the other objects in
the package.

<b>change log</b>an
MASSIVE MERGE: Starting version 1.0.0 on 15/07/2013

@date 16/07/2013
*/

/**
@mainpage
This  is the documentation for Edinburgh Topographic Analysis Package (ETAP),
incorporating LSDRaster.

These pages will help you get started using this software.

\image html ./logo.png


Tools are included to:
- Generate topographic metrics
- Perform Chi analysis
- And other important science stuff
.

@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDSoilHydroRaster_H
#define LSDSoilHydroRaster_H

#include <string>
#include <vector>
#include <map>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDShapeTools.hpp"
#include "LSDFlowInfo.hpp"
using namespace std;
using namespace TNT;


class LSDSoilHydroRaster: public LSDRaster
{
  public:

    /// @brief The default constructor. Does nothing
    /// @author SMM
    /// @date 11/11/2015
    LSDSoilHydroRaster()     { create(); }

    /// @brief Create a SoilHydroRaster by copying an LSDRaster
    /// @param ThisRaster The LSDRaster to be copied
    /// @author SMM
    /// @date 11/11/2015
    LSDSoilHydroRaster(LSDRaster& ThisRaster)
      { create(ThisRaster); }

    /// @brief Create a SoilHydroRaster by copying an LSDRaster
    /// @param ThisRaster The LSDRaster to be copied
    /// @param value the value which all the data elements take
    /// @author SMM
    /// @date 11/11/2015
    LSDSoilHydroRaster(LSDRaster& ThisRaster, float value)
      { create(ThisRaster, value); }

    /// @brief Creates an LSDSoilHydroRaster object populated with parameter
    /// values stored in OtherRaster.
    ///
    /// @details Fills nodata patches in OtherRaster which are non-nodata in DEM with either a min
    /// or max value taken from OtherRaster, depending on the if min_max is 0 (min) or 1 (max).
    /// @param DEM LSDRaster containing a DEM, used to ensure that every point on the landscape has a parameter value.
    /// @param OtherRaster LSDRaster containing spatial parameter values.
    /// @param min_max Integer flag set to either 0 (min) or 1 (max) to set how the data gaps are filled.
    /// @author SWDG
    /// @date 18/7/16
    LSDSoilHydroRaster(LSDRaster& DEM, LSDRaster& OtherRaster, int min_max)
      { create(DEM, OtherRaster, min_max); }

    /// @brief Create an LSDSoilHydroRaster from memory.
    /// @return LSDRaster
    /// @param nrows An integer of the number of rows.
    /// @param ncols An integer of the number of columns.
    /// @param xmin A float of the minimum X coordinate.
    /// @param ymin A float of the minimum Y coordinate.
    /// @param cellsize A float of the cellsize.
    /// @param ndv An integer of the no data value.
    /// @param data An Array2D of floats in the shape nrows*ncols,
    ///containing the data to be written.
    /// @author SMM
    /// @date 11/11/2015
    LSDSoilHydroRaster(int nrows, int ncols, float xmin, float ymin,
            float cellsize, float ndv, Array2D<float> data)
      { create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }

    /// @brief Create an LSDSoilHydroRaster from memory, includes georeferencing
    /// @return LSDRaster
    /// @param nrows An integer of the number of rows.
    /// @param ncols An integer of the number of columns.
    /// @param xmin A float of the minimum X coordinate.
    /// @param ymin A float of the minimum Y coordinate.
    /// @param cellsize A float of the cellsize.
    /// @param ndv An integer of the no data value.
    /// @param data An Array2D of floats in the shape nrows*ncols,
    /// @param temp_GRS a map of strings containing georeferencing information. Used
    /// mainly with ENVI format files
    ///containing the data to be written.
    /// @author SMM
    /// @date 11/11/2015
    LSDSoilHydroRaster(int nrows, int ncols, float xmin, float ymin,
        float cellsize, float ndv, Array2D<float> data, map<string,string> temp_GRS)
          { create(nrows, ncols, xmin, ymin, cellsize, ndv, data, temp_GRS); }

    /// @brief This sets all NoData pixels to value
    /// @param value The vaile to which all non nodata pixels should be set
    /// @author SMM
    /// @date 11/11/2015
    void SetHomogenousValues(float value);


    /// @brief This calculates snow depth as a function of elevation using
    /// two linear segments after this thesis: http://escholarship.org/uc/item/9zn1c1mk#page-7
    /// @param SlopeAscend The slope on the ascending limb in (g cm^-2 m^-1)
    /// @param SlopeDescend The slope on the descending limb in (g cm^-2 m^-1)
    /// @param PeakElevation The elevation of the peak snowpack (m)
    /// @param PeakSnowpack The peak snowpack in g cm^-2
    /// @param Elevation the Elevation LSDRaster
    /// @author SMM
    /// @date 11/11/2015
    void SetSnowEffDepthBilinear(float SlopeAscend, float SlopeDescend, float PeakElevation,
                                 float PeakSnowpack, LSDRaster& Elevation);


    /// @brief This calculates snow depth as a function of elevation using
    ///   a richards sigmoidal growth model, follwing this paper:
    ///   http://onlinelibrary.wiley.com/doi/10.1002/2015GL063413/epdf
    ///  some reasonable parameters might be:
    ///  MaximumSlope =
    ///  v =
    ///  lambda =
    /// @details The latex code for the equation \text{SWE}=A \left(v \exp \left(\frac{M (v+1)^{\frac{1}{v}+1} (\lambda -\zeta )}{A}+v+1\right)+1\right)^{-1/v}
    /// @param MaximumEffDepth The maximum effective depth (g cm^-2 m^-1)
    /// @param MaximumSlope the maximum slope of the curve. Should probably be less than 0.1
    /// @param v a "shape" parameter that controls how sharp the transitions from no snow
    ///  to snow are. A value between 0.01 and 3 is reasonable.
    /// @param lambda a "location" parameter, which sets where the curve is steep
    ///  this should be set to an elevation where you want snow to really increase
    /// @param Elevation the Elevation LSDRaster
    /// @author SMM
    /// @date 11/11/2015
    void SetSnowEffDepthRichards(float MaximumEffDepth, float MaximumSlope, float v,
                          float lambda, LSDRaster& Elevation);


    /// @brief This is an extremely naivie "landlside" model that is just
    ///  used to generate self shielding rasters to test the CRN basinwide code
    ///  IT SHOULD NOT IN ANY WAY BE USED TO PREDICT SLOPE STABILITY
    /// @param FilledElevation an elevation raster: needs to be filled or the code will crash!
    /// @param initiationPixels and interger giving the number of pixels
    ///  that "initiates" a landslide. Lnadlides then form upslope of pixels meeting
    ///  this threshold
    /// @param MinPixels pixels upslope of the initiation point will not
    ///  become part of the landslide unless they have at least this many contributing
    ///  pixels
    /// @param landslide_thickness The thickness of the landslides. For CRN purposes
    ///   this should be in units g cm^-2
    /// @author SMM
    /// @date 12/11/2015
    void NaiveLandslide(LSDRaster& FilledElevation, int initiationPixels,
                                      int MinPixels, float landslide_thickness);

    /// @brief Calculate h, the soil depth normal to the slope, used in the factor of safety equation.
    ///
    /// @details Call with the soil thickness raster.
    /// @param Slope The arctan of the slope.
    /// @return A raster of soil depth, D.
    /// @author SWDG
    /// @date 13/6/16
    LSDSoilHydroRaster Calculate_h(LSDRaster& Slope);

    /// @brief Calculate w, a hyrdological index, used in the factor of safety equation.
    ///
    /// @details Call with the ratio of recharge to transmissivity.
    /// @param Slope The arctan of the slope.
    /// @param DrainageArea The drainage area raster.
    /// @return A raster of the hydrological index w.
    /// @author SWDG
    /// @date 13/6/16
    LSDSoilHydroRaster Calculate_w(LSDRaster& Slope, LSDRaster& DrainageArea);

    /// @brief Calculate r, the soil to water density ratio, used in the factor of safety equation.
    ///
    /// @details Call with the soil density raster.
    /// @param rhoW A constant value for water density.
    /// @return A raster of the density ratio, r.
    /// @author SWDG
    /// @date 13/6/16
    LSDSoilHydroRaster Calculate_r(float& rhoW);

    /// @brief Calculate C, a cohesion index, used in the factor of safety equation.
    ///
    /// @details Call with the root cohesion raster.
    /// @param Cs The soil cohesion.
    /// @param h The soil thickness, normal to the slope.
    /// @param rhoS The soil density.
    /// @param g The gravitational constant.
    /// @return A raster of cohesion, C.
    /// @author SWDG
    /// @date 13/6/16
    LSDSoilHydroRaster Calculate_C(LSDSoilHydroRaster& Cs, LSDSoilHydroRaster& h, LSDSoilHydroRaster& rhoS, float& g);


    /// @brief Calculate the factor of safety using the sinmap definition.
    ///
    /// @details Call with the cohesion index, C.
    /// @param Slope The arctan of the slope.
    /// @param w The hydrological index.
    /// @param r The density ratio.
    /// @param phi The friction angle of the soil.
    /// @return A raster of factor of safety values.
    /// @author SWDG
    /// @date 13/6/16
    LSDSoilHydroRaster Calculate_sinmap_Fs(LSDRaster& Slope, LSDSoilHydroRaster& w, LSDSoilHydroRaster& r, LSDSoilHydroRaster& phi);

    /// @brief Calculate the the sinmap Stability Index (SI).
    ///
    /// @details This is a wrapper around a lightly modified port of the original sinmap 2.0 implementation.
    /// call with any SoilHydroRaster, it's values are used for identification of NoDataValues.
    /// @param Slope LSDRaster of slope.
    /// @param DrainageArea LSDRaster of drainage area from D-infinity
    /// @param lo_C LSDSoilHydroRaster of lower bound dimensionless cohesion.
    /// @param hi_C LSDSoilHydroRaster of upper bound dimensionless cohesion.
    /// @param lo_phi LSDSoilHydroRaster of lower bound friction angle.
    /// @param hi_phi LSDSoilHydroRaster of upper bound friction angle.
    /// @param lo_RoverT LSDSoilHydroRaster of lower bound ratio of reharge to transmissivity.
    /// @param hi_RoverT LSDSoilHydroRaster of upper bound ratio of reharge to transmissivity.
    /// @param r LSDSoilHydroRaster of water to soil density ratio.
    /// @param lo_FS LSDSoilHydroRaster of minimum factor of safety values.
    /// @param hi_FS LSDSoilHydroRaster of maximum factor of safety values.
    /// @return An LSDSoilHydroRaster coded with stability index values.
    /// @author SWDG
    /// @date 15/6/16
    LSDSoilHydroRaster Calculate_sinmap_SI(LSDRaster Slope, LSDRaster DrainageArea, LSDSoilHydroRaster lo_C, LSDSoilHydroRaster hi_C, LSDSoilHydroRaster lo_phi, LSDSoilHydroRaster hi_phi, LSDSoilHydroRaster lo_RoverT, LSDSoilHydroRaster hi_RoverT, LSDSoilHydroRaster lo_r, LSDSoilHydroRaster hi_r, LSDSoilHydroRaster lo_FS, LSDSoilHydroRaster hi_FS);

    /// @brief Get the maximum value in the SoilHydroRaster, ignoring nodata.
    /// @return The maximum value.
    /// @author SWDG
    /// @date 15/7/16
    float get_maximum_value();

    /// @brief Get the minimum value in the SoilHydroRaster, ignoring nodata.
    /// @return The maximum value.
    /// @author SWDG
    /// @date 15/7/16
    float get_minimum_value();

    /// @brief Calculate the average stability index value for a given set of upslope nodes.
    /// @param Nodes a vector containing all the nodes to be averaged over.
    /// @param FlowInfo an LSDFlowInfo object.
    /// @return The average value.
    /// @author SWDG
    /// @date 25/7/16
    float AverageSI(vector<int>& Nodes, LSDFlowInfo& FlowInfo);

    /// @brief Calculate the average stability index value for a vector containing many sets of upslope nodes.
    /// @param Vector_of_Nodes a vector containing vectors of the nodes to be averaged over.
    /// @param FlowInfo an LSDFlowInfo object.
    /// @return A vector of average values.
    /// @author SWDG
    /// @date 25/7/16
    vector<float> AverageSIs(vector< vector<int> >& Vector_of_Nodes, LSDFlowInfo& FlowInfo);

    /// @brief Write average SI values to spatial units defined by a vector of vectors of nodes.
    /// @param SIs a vector of average values generated by AverageSIs().
    /// @param Vector_of_Nodes a vector containing vectors of the nodes to be written.
    /// @param FlowInfo an LSDFlowInfo object.
    /// @return An LSDSoilHydroRaster with the average value coded to spatial units.
    /// @author SWDG
    /// @date 25/7/16
    LSDSoilHydroRaster WriteAvgSIs(vector<float> SIs, vector< vector<int> >& Vector_of_Nodes, LSDFlowInfo& FlowInfo);

    /// @brief Implementing the terrain shape index as used by Bolstad et al (1998).
    ///
    /// @details From Bolstad et al (1998) Predicting Southern Appalachian overstory
    /// vegetation with digital terrain data which is based on work by McNab (1989)
    /// Terrain Shape Index: Quantifying Effect of Minor Landforms on Tree Height
    ///
    /// These papers use a lot of magic numbers, which turn out to be fudge factors for
    /// working with DEMs with resolutions reported in YARDS. This code will cope with any
    /// stupid spatial unit you choose.
    /// @return An LSDSoilHydroRaster of the TerrainShapeIndex.
    /// @author SWDG
    /// @date 24/1/17
    LSDSoilHydroRaster TerrainShapeIndex();


    /// @brief Method to create a vegetation map based on Bolstad 1998. Calibrated for the Southern Appalachians.
    ///
    /// @details Works best on ~30 meter resolution data. Lidar is too noisy, so resample data before running.
    /// @param DEM LSDRaster of elevlvation values.
    /// @return An LSDIndexRaster coded with the vegetation types: 1 = Cove Hardwood, 2 = Mixed Deciduous, 3 = Xeric Oak-Pine 4 = Northern Hardwood.
    /// @author SWDG
    /// @date 26/1/17
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // Method to create a vegetation map based on Bolstad 1998. Calibrated for the
    // Southern Appalachians.
    //
    // Call with an LSDSoilHydroRaster generated by TerrainShapeIndex()
    // SWDG - 26/1/17
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    LSDIndexRaster BolstadVegetationMap(LSDRaster& DEM);

  protected:

  private:
    void create();
    void create(LSDRaster& OtherRaster);
    void create(LSDRaster& OtherRaster, float value);
    void create(LSDRaster& DEM, LSDRaster& OtherRaster, int min_max);
    void create(int ncols, int nrows, float xmin, float ymin,
                float cellsize, float ndv, Array2D<float> data);
    void create(int ncols, int nrows, float xmin, float ymin,
                float cellsize, float ndv, Array2D<float> data, map<string,string> GRS);

};

#endif
