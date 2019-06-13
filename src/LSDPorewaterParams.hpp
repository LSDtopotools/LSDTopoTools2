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
#include "LSDRaster.hpp"
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

    /// @brief Create a LSDPorewaterParams object using parameter files
    /// @param paramfile_path path to the parameter file
    /// @param paramfile_name The name of the parameter file without the path
    /// @author SMM
    /// @date 11/11/2015
    LSDPorewaterParams(string paramfile_path, string paramfile_name)
      { create(paramfile_path,paramfile_name); }
      
      
    /// @return The vector of depths
    vector<float> get_Depths() const { return Depths; }

    /// @return slope in radians
    float get_alpha() const { return alpha;}

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

    /// @brief A general function to load csv data 
    /// @param filename The name of the csv file INCLUDING path and extension
    /// @return a data map with headers as keys and vectors of strings as values. 
    /// @author SMM
    /// @date 15/01/2019
    map<string, vector<string> > load_csv_data_into_data_map(string filename);

    /// @brief A function for converting string vecs into float vecs. 
    ///  Changes NaN or nan values to -9999 
    /// @param string_vec a vector of strings
    /// @return a vector of floats. 
    /// @author SMM
    /// @date 15/01/2019
    vector<float> convert_string_vec_to_float(vector<string>& string_vec);

    /// @brief A function for converting string vecs into int vecs. 
    ///  Changes NaN or nan values to -9999 
    /// @param string_vec a vector of strings
    /// @return a vector of floats. 
    /// @author SMM
    /// @date 15/01/2019
    vector<int> convert_string_vec_to_int(vector<string>& string_vec);  


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

    /// @brief This function parses a rainfall file. It is more flexible in that it loads
    ///  a csv file and extracts both intensities and durations
    /// @param path the path to the rainfall file
    /// @param filename the name of the file including extension (but this needs to be
    ///  a csv file)
    /// @param durations a float vector holding durations
    /// @param intensities a float vector holding rainfall intensities. 
    /// @author SMM
    /// @date 15/01/2019
    void parse_rainfall_file(string path, string filename, vector<double>& intensities, vector<double>& durations);


    /// @brief Extract duration-intensity from a csv file preprocessed. it ingest data in s and mm/yrs and convert it to the right units (m) and dimentionlessize the time.
    /// @param durations a float vector holding durations
    /// @param intensities a float vector holding rainfall intensities. 
    /// @author BG
    /// @date 01/02/2019
    void get_duration_intensity_from_preprocessed_input(vector<float>& duration_s, vector<float>& intensity_mm);


    ///@brief I will comment after lunch
    ///@author BG
    void combine_time_to_durations_with_intensities( vector<double>& raw_time, vector<double>& raw_intensity, vector<double>& durations, vector<double>& intensities);

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


    /// @brief Return the rainfall csv name
    string get_rainfall_csv_name(){return rainfall_csv_name;};

    string get_path_csv(){return rainfall_csv_path;};

    /// returns the write prefix
    string get_saving_prefix(){return saving_prefix;};

    /// @brief return the slope raster
    LSDRaster get_alpha_raster(){return alpha_raster;}

    /// @brief returns spatially variable alpha. Row col on alpha raster.
    float get_alpha_rowcol(int row, int col) {return alpha_raster.get_data_element(row,col);}

    /// @brief Is full 1D analysis on?
    bool get_full_1D_output() {return full_1D_output;}
    /// @brief Is 2D analysis on?
    bool get_output_2D() {return output_2D;};

    /// @brief Get the number of threads to use for the spatial analysis
    int get_n_threads() {return n_threads;};
    /// @brief get the discrete time (in seconds, same as preprocessed input) of the spatial analysis
    float get_time_of_spatial_analysis() {return time_of_spatial_analysis;};



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
    
    /// A parameter that describes the steady state pressure profile
    /// It is equal to cos^2 alpha - (Iz_over_K_steady) and alpha is the slope angle
    /// see Iverson 2000 page 1902
    float beta;
    
    /// The infiltration rate at steady state. Dimensionless since I_z is in m/s and K is in m/s
    float Iz_over_K_steady;
    
    /// The friction andgle (in radians)
    float friction_angle;
    
    /// The cohesion in Pa
    float cohesion;
    
    /// The weight of water (density times gravity; use SI units (kg/(s*m^2))
    float weight_of_water;
    
    /// The weight of soil (density times gravity; use SI units (kg/(s*m^2))
    float weight_of_soil;

    /// Name of the csv file hosting the rainfall data from param
    string rainfall_csv_name;

    /// path to the csv file
    string rainfall_csv_path;

    /// Saving prefix
    string saving_prefix;

    /// alpha raster
    LSDRaster alpha_raster;

    /// Few bools
    bool full_1D_output;
    bool output_2D;

    /// OMP param
    int n_threads;

    /// Time at which the spatial analysis need to be processed
    float time_of_spatial_analysis;


  
  private:
  
    void create();
    void create(string paramfile_path, string paramfile_name);
  
};

#endif