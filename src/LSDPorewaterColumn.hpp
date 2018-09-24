//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDPorewaterColumn
// Land Surface Dynamics PorewterColumn object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object calculates porewater pressure based on the Iverson 2000 WRR
// model
//
// Developed by:
//  Simon M. Mudd
//  Stuart W.D. Grieve
//
// Copyright (C) 2016 Simon M. Mudd 2013 6
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


/** @file LSDPorewaterColumn.hpp
@author Simon M. Mudd, University of Edinburgh
@author Stuart W. D. Grieve, University of Edinburgh

**/



#ifndef LSDPorewaterColumn_H
#define LSDPorewaterColumn_H

#include <string>
#include <vector>
#include <map>
#include "LSDPorewaterParams.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;


class LSDPorewaterColumn
{
  public:
    /// @brief The default constructor. Does nothing
    /// @author SMM
    /// @date 17/11/2016
    LSDPorewaterColumn()     { create(); }

    /// @brief Create a 
    /// @param ThisRaster The LSDRaster to be copied
    /// @author SMM
    /// @date 17/11/2016
    LSDPorewaterColumn( vector<float> initial_Psi )
      { create(initial_Psi); }

    /// @brief Create a SoilHydroRaster by copying an LSDRaster
    /// @param ThisRaster The LSDRaster to be copied
    /// @author SMM
    /// @date 17/11/2016
    LSDPorewaterColumn( LSDPorewaterParams LSDPP )
      { create(LSDPP); }
      
    /// @brief this prints the pressure to screen
    /// @author SMM
    /// @date 17/11/2016
    void PrintPsiToScreen();
    
    /// @brief This calculates the response function of the porewater to rainfall
    ///  based on Iverson 2000 WRR equation 27e
    /// @param t_star the dimensionless time
    /// @return The response function
    /// @author SMM
    /// @date 17/11/2016
    float CalculateResponseFunction(float t_star);

    /// @brief This calculates transient component of Psi
    /// @param  LSDPP an LSDPorewaterParams object
    /// @param Iz_over_Kz the infilatration rate of this response
    /// @param t_star the dimensionless time
    /// @param T_star the dimensionless storm duration
    /// @return the transient component of Psi for this storm
    /// @author SMM
    /// @date 17/11/2016
    vector<float> CalcualtePsiTransient(LSDPorewaterParams& LSDPP, float Iz_over_Kz, 
                           float t_star, float T_star);

    /// @brief This calculates the transient component of Psi for a given
    ///  rainfall pulse using dimensional tie
    /// @param  LSDPP an LSDPorewaterParams object
    /// @param Iz_over_Kz the infilatration rate of this response
    /// @param t the time in seconds
    /// @param T the storm duration in seconds
    /// @return the transient component of Psi for this storm
    /// @author SMM
    /// @date 17/11/2016
    vector<float> CalculatePsiDimensionalTimeTransient(LSDPorewaterParams& LSDPP, float t, float T, float Iz_over_Kz);

    /// @brief This function takes a time series of rainfall and calculates
    ///  the pore pressures
    /// @param durations storm durations in a vector of floats in seconds
    /// @param intensities storm intensities, in dimensionless Iz_over_Kz
    /// @param LSDPP is an PorewaterParams object
    /// @param t is the dimensional time
    /// @author SMM
    /// @date 17/11/2016
    void CalculatePsiFromTimeSeries(vector<float> durations, vector<float> intensities, 
                                LSDPorewaterParams& LSDPP, float t);


    /// @brief Calcualtes the friction factor component of the factor of safety
    /// @param LSDPP is an PorewaterParams object
    /// @return F_f 
    /// @author SMM
    /// @date 25/11/2016
    float F_f(LSDPorewaterParams& LSDPP);

    /// @brief Calcualtes the friction factor from cohesion component of the factor of safety
    /// @param LSDPP is an PorewaterParams object
    /// @return F_c (a vector of floats)
    /// @author SMM
    /// @date 25/11/2016
    vector<float> F_c(LSDPorewaterParams& LSDPP);

    /// @brief Calcualtes the friction factor from pore pressure component of the factor of safety
    /// @param LSDPP is an PorewaterParams object
    /// @return F_w (a vector of floats)
    /// @author SMM
    /// @date 25/11/2016
    vector<float> F_w(LSDPorewaterParams& LSDPP);

    /// @brief Calcualtes the factor of safety
    /// @param LSDPP is an PorewaterParams object
    /// @return FS (a vector of floats)
    /// @author SMM
    /// @date 26/11/2016
    vector<float> FS(LSDPorewaterParams& LSDPP);
    
    /// @brief This calcualtes the depth of failer of the column. 
    /// @param LSDPP is a PorewaterParams object
    /// @param minimum_depth the minimum depth at which falure can occur
    /// @return depth_of_failure This is -9999 if there is no failure
    /// @author SMM
    /// @date 02/12/2016
    float DepthOfFailure(LSDPorewaterParams& LSDPP, float minimum_depth);

    /// @brief This finds the minimum facotr of safety in a column
    /// @param LSDPP is a PorewaterParams object
    /// @param minimum_depth the minimum depth at which falure can occur
    /// @param depth_of_minFS The depth where the minimum FS occurs
    /// @param minFS The minimum factor of safety value
    /// @author SMM
    /// @date 02/12/2016
    void GetMinFS(LSDPorewaterParams& LSDPP, float minimum_depth, float& depth_of_minFS, float& minFS);
    
    /// @brief This goes through a time series trying to see when the first
    ///  failure occurs and at what depth
    /// @param durations storm durations in a vector of floats in seconds
    /// @param intensities storm intensities, in dimensionless Iz_over_Kz
    /// @param LSDPP is a PorewaterParams object
    /// @param minimum_depth the minimum depth at which falure can occur
    /// @author SMM
    /// @date 02/12/2016
    void ScanTimeseriesForFailure(vector<float> durations, vector<float> intensities,
                                   LSDPorewaterParams& LSDPP, float minimum_depth, 
                                   vector<float> times);
    
    
    /// @return The vector of depths
    vector<float> get_Psi() const { return Psi; }

    /// @return slope in radians
    int get_row() const { return row; }

    /// @return slope in radians
    int get_col() const { return col; }
    
    /// @return slope in radians
    int get_node_index() const { return node_index; }


  protected:
    /// This holds the pressure heads
    vector<float> Psi;
    
    /// This has the information about where the column is
    int row;
    int col;
    int node_index;
  
  private:
  
    void create();
    void create(vector<float> InitialPsi);
    void create(LSDPorewaterParams LSDPP);
  
};

#endif