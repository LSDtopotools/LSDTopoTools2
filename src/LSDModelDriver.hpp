//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDModelDriver
// Land Surface Dynamics Model(s) Driver
//
// This object parses parameter files and drives the models:
//
// LSDRasterModel
// LSDCatchmentModel (coming soon to a repository near you)
//
// Its purpose is to stop having to write a bunch of .cpp driver functions and instead
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
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// This object is written by
/// @author Simon M. Mudd, University of Edinburgh
/// @author David T. Milodowski, University of Edinburgh
/// @author Martin D. Hurst, British Geological Survey
/// @author Fiona Clubb, University of Edinburgh
/// @author Stuart Grieve, University of Edinburgh
/// @author James Jenkinson, University of Edinburgh
/// @author Declan Valters, University of Manchester
///
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
///
/// Version 0.0.1		2015-01-15
///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <string>

#include "LSDAnalysisDriver.hpp"
#include "LSDRasterModel.hpp"
#include "LSDCatchmentModel.hpp"

#ifndef LSDModelDriver_H
#define LSDModelDriver_H



/// @brief This is a class to manage running LSDTopoTools. It parses a parameter
/// file and then manages running of analyses. 
/// @details The intention of this object is to run analyses via parameter
/// files and not through numerous compiled driver functions. We eventually
/// will want some kind of 'recorder' so that any time this object
/// runs an analysis it gives a full report of what analyses were run so that
/// results are reproducable
class LSDModelDriver : public LSDCatchmentModel, LSDRasterModel  //(not sure about this - perhaps keep them totally separate for now?)
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
    /// @author DAV
    /// @date 2015-01-16    
    LSDModelDriver()
    { 
		create(); 
	}
	
    /// @brief this constructor just reads the param file given by the path and
    /// filename. You must give the parameter file extension!
    /// @param pname the pathname to the parameter file
    /// @param fname the filename of the parameter file !!INCLUDING EXTENSION!!
    /// @author DAV
    /// @date 2015-01-16	
    LSDModelDriver(string pname, string pfname)			
    { 
		create(pname, pfname); 
	}
	
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
    
    void select_model_from_model_map();
    
    void run_LSDCatchmentModel_components();
    
    void run_LSDRasterModel_components();
	
	protected:
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
	
	private:
	
	void check_pathname_for_slash();
	
	/// This holds names of the different model cores that can be run from the 
	/// Model Driver
	std::map<string,string> model_map;

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // CREATE FUNCTIONS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
    /// @brief default create function
    /// @details this asks for a pathname and a filename of the parameter file
    /// It then opens the paramter file and ingests the information
    /// @author DAV
    /// @date 2015-01-16
    void create();	

    /// @brief create function
    /// @param pname the pathname to the parameter file
    /// @param fname the filename of the parameter file !!INCLUDING EXTENSION!!
    /// @author DAV
    /// @date 201-01-16
    void create(string pname, string fname);	
	
  
	
};
	
	
#endif	
	
	
	
	
	
	
	
