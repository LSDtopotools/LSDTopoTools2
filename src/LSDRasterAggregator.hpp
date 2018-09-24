//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterAggregator.hpp
// Land Surface Dynamics RasterAggregator
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  that is a very general object for managing aggregated raster calculations
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
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
#include <fstream>
#include <math.h>
#include <iostream>
#include <map>
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
using namespace std;

#ifndef RasterAggregator_HPP
#define RasterAggregator_HPP

/// @brief A general object for holding several different raster layers
class LSDRasterAggregator
{
  public:
  
    /// @brief the default constructor. This doesn't do anything. 
    LSDRasterAggregator()                                { create(); }
    
    /// @brief This constructor requires a filename
    /// @detail The filename and path point to files with the rasters
    /// @param path The path to the files
    /// @param file_prefix the prefix of the files
    /// @author SMM
    /// @date 10/02/2016
    LSDRasterAggregator( string path, string file_prefix)  { create(path,file_prefix); }
    
    /// @brief This function load a csv file containing names of  a DEMs and 
    ///  (possibly) other rasters
    /// @detail The first row of the file contains column headers and is ignored
    ///  Thereafter you have two columns in each row, comma seperated, 
    ///  with the raster type as the first column and the raster filename
    ///  (with full path) as the second. 
    /// @param filename the name of the file
    /// @author SMM
    /// @date 10/02/2016
    void load_raster_filenames(string filename);
    
    /// @brief This loads parameters It creates a map of parameter values
    /// @param filename a string of the full filename
    /// @author SMM
    /// @date 10/02/2016
    void load_parameters(string filename);
   
    /// @brief this function checks the existence and georeferencing of 
    ///  the rasters outlined in the file list
    /// @author SMM
    /// @date 10/02/2016
    void check_rasters();


  protected:
    
    /// the path to the cosmo data. Also used to print results
    string path;
    
    /// the prefix of the parameter files
    string param_name;
    
    /// A map for holding the different raster filenames
    /// The index into the raster filenames is the raster type. 
    /// There are some commone raster types:
    /// DEM is the DEM
    /// ---more to come here---
    map< string,string > raster_filenames;
    
    /// A map holding the parameter values. These are stored as strings
    ///  and converted to the appropriate data type as needed
    map<string,string> parameter_map;
    
    /// The minimum slope for the fill function
    float min_slope;
    
    /// The boundary conditions for flow info calculations
    vector<string> boundary_conditions;
  
  private:
  
    /// @brief the empty create function
    void create();
    
    /// @brief the create function used when calling a file
    /// @param the path to the file. Must have the "/" at the end
    /// @param  file_prefix the prefix (without extension) of the parameter files
    /// @author SMM
    /// @date 10/02/2016
    void create(string path, string file_prefix);
};

/// @brief A derived class that is used to compute erosion rates based on 
///  concentrations of in-situ cosmogenic nuclides such as 10Be and 26Al
class LSDSedimentRouting: public LSDRasterAggregator
{

  public:

    /// @brief This constructor requires a filename
    /// @detail The filename and path point to files with the rasters
    /// @param path The path to the files
    /// @param file_prefix the prefix of the files
    /// @author SMM
    /// @date 10/02/2016
    LSDSedimentRouting( string path, string file_prefix)  { create(path,file_prefix); }
  #
    /// @brief This checks for the parameter values specific to the sediment routing routines
    /// @author SMM
    /// @date 10/02/2016
    void check_parameter_values();
    
    /// @brief Prints the parameter values to screen
    /// @author SMM
    /// @date 10/02/2016
    void print_parameter_values_to_screen();

    /// @brief This calculates the suspended and bedload for a given node in the DEM
    ///  It can be called repeatedly to get the entire DEM
    /// @param node The node at which you want to calculate suspended and bedload
    /// @param FlowInfo the LSDFlowinfo object
    /// @param RasterVec a vector of LSDRasters, the first element is the DEM, 
    ///  the second element is the lithology, the third element is the flow distance
    ///  and the fourth element is the erosion rate
    /// @return a vector containg the bedload and suspended load information
    /// @author SMM
    /// @date 22/03/2016
    vector<float> calculate_suspended_and_bedload(int node, LSDFlowInfo& FlowInfo,
                                                    vector<LSDRaster> RasterVec)


  protected:
  
    /// The number of lithologies
    int N_lithologies;
    
    /// the erosion rate in mm/yr
    float erosion_rate;
    
    /// The erodibiliy coefficients in km^-1
    /// The int is the index into the lithology (coded with an integer)
    map<int,float> erodibility_coefficients;
    
    /// The fertility coefficients: states the fraction in the source material
    ///  that contains zircon
    /// The int is the index into the lithology (coded with an integer)
    map<int,float> fertility_coefficients;
    
    /// The fraction of the source material that is entered as suspended load
    /// The int is the index into the lithology (coded with an integer)
    map<int,float> source_1mm;

    
  private:

    /// @brief the create function used when calling a file
    /// @param the path to the file. Must have the "/" at the end
    /// @param  file_prefix the prefix (without extension) of the parameter files
    /// @author SMM
    /// @date 10/02/2016
    void create(string path, string file_prefix);

};



#endif