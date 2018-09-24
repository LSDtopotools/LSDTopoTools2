//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDPorewaterParams
// Land Surface Dynamics PorewterParams object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object interfaces with teh porewater column object
//  In a landsacpe each pixel will have its own pore ressures but 
//  the parameters will be constant (or have similar statistical properties)
//  across a landscape. This structure tries to minimize memory requirements
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



#ifndef LSDPorewaterParams_H
#define LSDPorewaterParams_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;


class LSDPorewaterParams
{
  public:
    /// @brief The default constructor. Does nothing
    /// @author SMM
    /// @date 11/11/2015
    LSDPorewaterParams()     { create(); }

    /// @brief Create a SoilHydroRaster by copying an LSDRaster
    /// @param ThisRaster The LSDRaster to be copied
    /// @author SMM
    /// @date 11/11/2015
    LSDPorewaterParams(string paramfile_path, string paramfile_name)
      { create(paramfile_path,paramfile_name); }
      
      
    /// @return The vector of depths
    vector<float> get_Depths() const { return Depths; }

    /// @return slope in radians
    float get_alpha() const { return alpha; }

    /// @return diffusivity
    float get_D_0() const { return D_0; }
    
    /// @return the saturated hydraulic conductivity
    float get_K_sat() const {return K_sat; }

    /// @return dimensionless diffusivity
    float get_D_hat() const { return D_hat; }

    /// @return depth to water table
    float get_d() const { return d; }

    /// @return friction_angle
    float get_friction_angle() const { return friction_angle; }

    /// @return cohesion
    float get_cohesion() const { return cohesion; }
    
    /// @return weight_of_water
    float get_weight_of_water() const { return weight_of_water; }
    
    /// @return weight_of_soil
    float get_weight_of_soil() const { return weight_of_soil; }
    
    /// @ brief Calculates the beta parameter and sets that data member
    /// @author SMM
    /// @date 17/11/2016
    void calculate_beta();

    /// @brief This calculates dimensionless diffusivity
    /// @author SMM
    /// @date 17/11/2016
    void calculate_D_hat();

    /// @brief converts weeks to seconds
    /// @param weeks a vector of weeks
    /// @return a vector of seconds
    /// @author SMM
    /// @date 25/11/2016
    vector<float> weeks_to_seconds(vector<float> weeks);

    /// @brief converts weeks to seconds
    /// @param weeks weeks
    /// @return seconds
    /// @author SMM
    /// @date 25/11/2016
    float weeks_to_seconds(float weeks);

    /// @brief converts days to seconds
    /// @param weeks a vector of days
    /// @return a vector of seconds
    /// @author SMM
    /// @date 25/11/2016
    vector<float> days_to_seconds(vector<float> days);

    /// @brief converts days to seconds
    /// @param days
    /// @return seconds
    /// @author SMM
    /// @date 25/11/2016
    float days_to_seconds(float days);

    /// @brief converts seconds to weeks
    /// @param seconds a vector of seconds
    /// @return a vector of weeks
    /// @author SMM
    /// @date 25/11/2016
    vector<float> seconds_to_weeks(vector<float> seconds);

    /// @brief converts seconds to weeks
    /// @param seconds
    /// @return weeks
    /// @author SMM
    /// @date 25/11/2016
    float seconds_to_weeks(float seconds);

    /// @brief converts seconds to days
    /// @param weeks a vector of seconds
    /// @return a vector of days
    /// @author SMM
    /// @date 25/11/2016
    vector<float> seconds_to_days(vector<float> seconds);

    /// @brief converts seconds to days
    /// @param seconds
    /// @return days
    /// @author SMM
    /// @date 25/11/2016
    float seconds_to_days(float seconds);



    /// @brief This calculates a steady state pressure profile
    /// @author SMM
    /// @date 17/11/2016
    vector<float> calculate_steady_psi();

    /// @brief This function parses a rainfall file. It only captures the
    ///  rainfall rates, which must have the csv header "rainfall_rate"
    /// @param path the path to the rainfall file
    /// @param filename the name of the file including extension (but this needs to be
    ///  a csv file)
    /// @param intensities a float vector holding rainfall intensities. 
    /// @author SMM
    /// @date 17/11/2016
    void parse_rainfall_file(string path, string filename, vector<float>& intensities);

    /// @brief This function parses a rainfall file derived from MIDAS. 
    /// @param path the path to the rainfall file
    /// @param filename the name of the file including extension (but this needs to be
    ///  a csv file)
    /// @param days the days since 1900. 
    /// @param intensities a float vector holding rainfall intensities. 
    /// @author SMM
    /// @date 22/11/2016
    void parse_MIDAS_rainfall_file(string path, string filename, vector<int>& days, vector<float>& intensities);

    /// @brief This function takes data vectors extracted from a parsed MIDAS
    ///  file (using the parse_MIDAS.py script) and returns an intensity duration
    ///  record
    /// @param days a vector of the days in int format of the rainfall measurements 
    /// @param intensities the intensity measurements: for MIDAS this is mm/day
    /// @param durations this is replaced with the durations in days of the records
    /// @author SMM
    /// @date 23/11/2016
    void parse_MIDAS_duration_intensities(vector<int>& days, vector<float>& intensities, vector<float>& durations);

    /// @brief Thes function prints the parameters to screen
    /// @author SMM
    /// @date 22/11/2016
    void print_parameters_to_screen();


  protected:
    /// This holds the depths in metres
    vector<float> Depths;
    
    /// This holds the slope in radians
    float alpha;
    
    /// The  hydraulic diffusivity in m^2/s
    float D_0;
    
    /// The saturated hydraulic conductivity in m/s
    float K_sat;
    
    /// The dimensionless hydraulic diffusivity
    float D_hat;
    
    /// The depth to saturated water table at steady infiltration in metres
    float d;
    
    /// A parameter that describes the stead state pressure profile
    float beta;
    
    /// The infiltration rate at steady state
    float Iz_over_K_steady;
    
    /// The friction andgle
    float friction_angle;
    
    /// The cohesion
    float cohesion;
    
    /// The weigth of water (density time gravity)
    float weight_of_water;
    
    /// The weight of soil (density times gravity)
    float weight_of_soil;
  
  private:
  
    void create();
    void create(string paramfile_path, string paramfile_name);
  
};

#endif