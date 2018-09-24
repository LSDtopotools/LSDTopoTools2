//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDAnalysisDriver.hpp
// Land Surface Dynamics Analysis Driver
//
// This object parses parameter files and drives analysis. Its purpose is
// to stop having to write a bunch of .cpp driver functions and instead
// be able to write a parameter files without compiling
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2014 Simon M. Mudd 2014
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
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDStatsTools.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDChannel.hpp"
using namespace std;

#ifndef LSDAnalysisDriver_H
#define LSDAnalysisDriver_H


/// @brief This is a class to manage running LSDTopoTools. It parses a parameter
/// file and then manages running of analyses. 
/// @details The intention of this object is to run analyses via parameter
/// files and not through numerous compiled driver functions. We eventually
/// will want some kind of 'recorder' so that any time this object
/// runs an analysis it gives a full report of what analyses were run so that
/// results are reproducable
class LSDAnalysisDriver
{
	public:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Constructors
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    /// @brief The default constructor. 
    /// @details this asks for a pathname and a filename of the parameter file
    /// It then opens the paramter file and ingests the information
    /// @author SMM
    /// @date 29/07/2014    
    LSDAnalysisDriver()			{ create(); }
	
    /// @brief this constructor just reads the param file given by the path and
    /// filename. You must give the parameter file extension!
    /// @param pname the pathname to the parameter file
    /// @param fname the filename of the parameter file !!INCLUDING EXTENSION!!
    /// @author SMM
    /// @date 29/07/2014	
    LSDAnalysisDriver(string pname, string pfname)			{ create(pname, pfname); }
 
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Main drivers of reading, computation and writing of data 
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-      
    /// @brief This is the main function for parsing the parameter file
    /// @param pathname the path to the paramter file
    /// @param param_fname the name of the parameter file
    /// @author SMM
    /// @date 29/07/2014
    void ingest_data(string pname, string p_fname);
    
    /// @brief This looks through the raster switches, and 
    /// calculates the necessary rasters. 
    /// @author SMM
    /// @date 29/07/2014
    void compute_rasters_from_raster_switches();

    /// @brief This writes rasters based on the analysis switches
    /// @author SMM
    /// @date 29/07/2014
    void write_rasters_from_analysis_switches();

    /// @brief This writes other forms of data based on the analysis switches
    /// @author SMM
    /// @date 29/07/2014
    void write_shapes_from_switches();

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Functions for getting individual datasets
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-      
    /// @brief This simply loads the base raster into the vector of rasters
    /// @author SMM
    /// @date 29/07/2014
    void read_base_raster();       
    
    /// @brief This calculates the fill raster
    /// @author SMM
    /// @date 29/07/2014
    void fill_raster();

    /// @brief This calculates the trimmed and hole filled raster
    /// @author SMM
    /// @date 10/12/2014
    void calculate_trimmed_and_nodata_filled();

    /// @brief This calculates the hillshade raster
    /// @author SMM
    /// @date 29/07/2014
    void calculate_hillshade();

    /// @brief This calculates the drainage area
    /// @author SMM
    /// @date 23/09/2014
    void calculate_drainage_area();
    
    /// @brief This calculates the slope raster
    /// @author SMM
    /// @date 21/09/2014
    void calculate_slope();

    /// @brief This calculates the slope angle raster
    /// @author SMM
    /// @date 13/11/2014
    void calculate_slope_angle();

    /// @brief This calculates the Factor of safety when the hillslopes
    /// are saturated
    /// @author SMM
    /// @date 13/11/2014
    void calculate_FS_sat();

    /// @brief This calculates the various polyfit rasters
    /// @author SMM
    /// @date 21/09/2014
    void calculate_polyfit();

    /// @brief This calculates a curvature threhsold maks
    /// @author SMM
    /// @date 09/09/2016
    void calculate_curvature_mask_threshold();

    /// @brief This calculates a threhsold maks
    /// @author SMM
    /// @date 28/09/2016
    void calculate_mask_threshold();

    /// @brief This calculates the LSDFlowInfo object
    /// @author SMM
    /// @date 29/07/2014
    void calculate_flowinfo();

    /// @brief This calculates the LSDFlowInfo object distance from the outlet
    /// @author SMM
    /// @date 24/09/2014 
    void calculate_flow_distance(); 
    
    /// @brief This calculates the LSDJunctionNetwork object
    /// @return JunctionNetork: the junction network
    /// @author SMM
    /// @date 30/07/2014
    LSDJunctionNetwork calculate_JunctionNetwork();    
        
    /// @brief This gets the nodeindex from the LSDFlowInfo object
    /// IMPORTANT: this is an LSDIndexRaster so it goes into the 
    /// LSDIndexRaster vector!    
    /// @author SMM
    /// @date 29/07/2014
    void calculate_nodeindex();
    
    /// @brief This gets the ContributingPixels from the LSDFlowInfo object
    /// IMPORTANT: this is an LSDIndexRaster so it goes into the 
    /// LSDIndexRaster vector!
    /// @author SMM
    /// @date 30/07/2014
    void calculate_ContributingPixels();

    /// @brief This gets sources for a channel network
    /// @author SMM
    /// @date 30/07/2014
    void calculate_sources();

    /// @brief This calculates Junction index raster
    /// @author SMM
    /// @date 30/07/2014
    void calculate_JunctionIndex(LSDJunctionNetwork& JN);

    /// @brief This calculates stream order raster
    /// @author SMM
    /// @date 30/07/2014
    void calculate_SOArray(LSDJunctionNetwork& JN);
    
    /// @brief This calculates the chi map raster
    /// @author SMM
    /// @date 30/07/2014
    void calculate_chi_map();
    
    /// @brief This retrieves a single thread channel
    /// @author SMM
    /// @date 24/09/2014    
    LSDChannel get_single_thread_channel();


	protected:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Housekeeping functions for making sure the object does not crash
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    

        
    /// @brief This checks to see if boundary condtions have been assigned and 
    /// if not defaults to no flux boundaries
    /// @author SMM
    /// @date 29/07/2014
    void check_boundary_conditions();

    /// @brief This check to see if the filenames, paths and extensions have
    /// been assigned. If not it changes these to defaults
    /// @author SMM
    /// @date 29/07/2014
    void check_file_extensions_and_paths();
    
    /// @brief this adds a slash to the end of the pathname
    /// @author SMM
    /// @date 29/07/2014
    void check_pathname_for_slash();	
 
    /// @brief this adds a slash to the end of the pathname
    /// an overloaded function: this one takes a string and returns the 
    /// string with a slash at the end
    /// @param this_pathname the string containing the pathname
    /// @retun the string with a slash at the end if it hasn't got one already
    /// @author SMM
    /// @date 19/09/2014    
    string check_pathname_for_slash(string this_pathname);   

    /// @brief This checks on the polyfit functions needed
    /// @author SMM
    /// @date 22/09/2014
    void check_polyfit();

    /// @brief this returns the string before the last dot in a string. 
    /// so for example if you gave it paramfile.param it would return paramfile
    /// @param this_string the string you want truncated
    /// @return the truncated string
    /// @author SMM
    /// @date 29/07/2014
    string get_string_before_dot(string this_string);


    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // DATA MEMBERS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
   
    /// This holds the flow info object
    LSDFlowInfo FlowInfo;
    
    /// the just tells the code if the flow info has already been calculated
    bool got_flowinfo;
    
    /// Pointer to junction network object
    //LSDJunctionNetwork JunctionNetwork;
    
    /// tells the code if you've got the junction network
    bool got_JunctionNetwork;
    
    /// tells the code if you've done the polyfit
    bool got_polyfit;
        
    /// the path to the datafiles
    string pathname;
    
    /// the name of the parameter file
    string param_fname;
    		
    /// Extension for reading DEMs. Correspondence with write extensions is checked
    string dem_read_extension;
    
    /// Extension for writing DEMs. Correspondence with read extensions is checked
    string dem_write_extension;
    
    /// Path to files being written. Default is pathname
    string write_path;
    
    /// file prefix of files to be written. Default is the param name prefix 
    string write_fname;

    /// Path to files being read. Default is pathname
    string read_path;
    
    /// file prefix of files to be written. Default is the param name prefix 
    string read_fname;       

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Parameters for various analyses
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    /// the four boundary conditions on the raster for the flow info object 
    vector<string> boundary_conditions;

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Switches for running different analyses
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    /// This map holds all the possible analyses
    map<string,bool> analyses_switches; 

    /// This vector holds various rasters computed during the run
    map<string,LSDRaster> map_of_LSDRasters;
    
    /// as above, but these are index rasters
    map<string,LSDIndexRaster> map_of_LSDIndexRasters;
  
    /// This is a map  container that determines if various rasters are needed for the
    /// analysis. This ensures things like the fill raster are only calculated
    /// once
    map<string,bool> raster_switches;

    /// This is a map that tell where the indices into the raster vecs are    
    //map<string,int> raster_indices;
    
    /// This holds float parameters
    map<string,float> float_parameters;
    
    /// This holds integer parameters
    map<string,int> int_parameters;
    
    /// This holds integer vectors. Can be used to get sources
    map<string, vector<int> > integer_vector_map;
    
    /// This holds names of supporting files, for example files that contain
    /// node of junction indices to be loaded. 
    map<string,string> support_file_names;

    /// This holds names of methods. For example, if the key is drainage_area_method, the string is
    /// the method which is used to calculate drainage area
    map<string,string> method_map;
    

    
  private:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // CREATE FUNCTIONS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
    /// @brief default create function
    /// @details this asks for a pathname and a filename of the parameter file
    /// It then opens the paramter file and ingests the information
    /// @author SMM
    /// @date 29/07/2014
    void create();	

    /// @brief create function
    /// @param pname the pathname to the parameter file
    /// @param fname the filename of the parameter file !!INCLUDING EXTENSION!!
    /// @author SMM
    /// @date 29/07/2014
    void create(string pname, string fname);	




};

#endif
