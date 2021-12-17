//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDParameterParser.hpp
//
// Land Surface Dynamics Parameter Parser Object
//
// An object for keeping track of parameters
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for calculating concntration of environmental tracers, CRNs, TCN, fallout
//  nuclides
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2016 Simon M. Mudd 2016
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
#include <string>
#include <vector>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;


#ifndef LSDParameterParser_HPP
#define LSDParameterParser_HPP

/// A function for dealing with initial ingestion of the parameter files
vector<string> DriverIngestor(int nNumberofArgs,char *argv[]);

/// @brief Object to perform flow routing.
class LSDParameterParser
{
  public:
    /// @brief The create function. This is default and throws an warning.
    /// @author SMM
    /// @date 11/02/16
    LSDParameterParser()        { create(); }

    /// @brief Creates a LSDParamParser object using a path and filename
    /// @param PathName the name of the path
    /// @param FileName the name of the parameter file
    /// @author SMM
    /// @date 3/07/2015
    LSDParameterParser(string PathName, string FileName)
                     { create(PathName,FileName); }

    /// @brief This is the function for reading parameters
    /// @ detail It is identical to the parse_line function in LSDStatsTools but
    ///   has a bigger buffer to cope with long path names
    /// @param infile The file being read
    /// @param pramater a string of the parameter name: returned from the parser
    ///   this is always changed to lower case only!!
    /// @param value A string holding a value. It will need to be further parsed
    ///   depending on the name of the parameter (i.e., to a float, int, string, etc)
    /// @author SMM
    /// @date 02/11/16
    void LSDPP_parse_line(ifstream &infile, string &parameter, string &value);

    /// @brief This takes a parameter file and parses all the parameters into
    ///  a map with string keys and values. Further functions are used to convert these
    ///  into the appropriate data types
    /// @author SMM
    /// @date 03/11/2016
    void LSDPP_parse_file_into_parameter_map(string FullName);

    /// @brief This parses all the required information for File I/O
    /// @detail You need to get the parameter map first
    /// @author SMM
    /// @date 03/11/2016
    void parse_file_IO();

    /// @brief This forces parsing of all the parameters into strings
    /// @detail It is used for copying parameter files
    /// @author SMM
    /// @date 16/03/2017
    void force_parse();

    /// @brief This parses all of the defalt parameter maps
    /// @param default_param_f a map of the default parameters, keys are string, values are floats
    /// @param default_param_i a map of the default parameters, keys are string, values are int
    /// @param default_param_b a map of the default parameters, keys are string, values are bool
    /// @param default_param_s a map of the default parameters, keys are string, values are string
    /// @author SMM
    /// @date 09/11/2016
    void parse_all_parameters(map<string,float> default_map_f,
                      map<string,int> default_map_i, map<string,bool> default_map_b,
                      map<string,string> default_map_s);

    /// @brief This parses all of the defalt parameter maps
    /// Overloaded to include doubles.
    /// @param default_param_f a map of the default parameters, keys are string, values are floats
    /// @param default_param_i a map of the default parameters, keys are string, values are int
    /// @param default_param_b a map of the default parameters, keys are string, values are bool
    /// @param default_param_s a map of the default parameters, keys are string, values are string
    /// @author FJC
    /// @date 20/11/17
    void parse_all_parameters(map<string,float> default_map_f,
                      map<string,int> default_map_i, map<string,bool> default_map_b,
                      map<string,string> default_map_s, map<string,double> default_map_d);

    /// @brief This processes the help information
    /// @param help_map the map of string vectors containing help information. 
    /// @param the file prefix of the help file
    /// @param the version the help file refers to 
    /// @param the appropriate citation for the code
    /// @author SMM
    /// @date 21/07/2021 
    void print_help(map< string, vector<string> > help_map, string file_prefix, string version, string citation);


    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are floats
    /// @author SMM
    /// @date 03/11/2016
    void parse_float_parameters(map<string,float> default_map);

    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are int
    /// @author SMM
    /// @date 03/11/2016
    void parse_int_parameters(map<string,int> default_map);

    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are bool
    /// @author SMM
    /// @date 03/11/2016
    void parse_bool_parameters(map<string,bool> default_map);

    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are strings
    /// @author SMM
    /// @date 03/11/2016
    void parse_string_parameters(map<string,string> default_map);

    /// @brief This function takes a default map and converts it into the parameters
    ///  by comparing the keys to the parsed parameter file
    /// @param default_param a map of the default parameters, keys are string, values are doubles
    /// @author FJC
    /// @date 20/11/2017
    void parse_double_parameters(map<string,double> default_map);

    /// @brief This function takes a parameter from the string map and parses it
    ///  into a vector of strings
    /// @param key the string that is the key into the string map
    /// @return A vector of strings
    /// @author SMM
    /// @date 09/11/2016
    vector<string> parse_string_vector(string key);

    /// @brief This function takes a parameter from the string map and parses it
    ///  into a vector of ints
    /// @param key the string that is the key into the string map
    /// @return A vector of ints
    /// @author SMM
    /// @date 09/11/2016
    vector<int> parse_int_vector(string key);

    /// @brief This function takes a parameter from the string map and parses it
    ///  into a vector of floats
    /// @param key the string that is the key into the string map
    /// @return A vector of floats
    /// @author SMM
    /// @date 09/17/2016
    vector<float> parse_float_vector(string key);

    /// @brief This forces the read and write extensions to bil
    /// @author SMM
    /// @date 02/11/16
    void force_bil_extension();

    /// @return read_extension
    string get_dem_read_extension() const        { return dem_read_extension; }
    /// @return write_extension
    string get_dem_write_extension() const        { return dem_write_extension; }
    /// @return write_path
    string get_write_path() const      { return write_path; }
    /// @return read_path
    string get_read_path() const      { return read_path; }
    /// @return write_fname
    string get_write_fname() const  { return write_fname; }
    /// @return read_fname
    string get_read_fname() const      { return read_fname; }
    /// @return CHeads_file
    string get_CHeads_file() const      { return CHeads_file; }
    /// @return BaseLevelJunctions_file
    string get_BaselevelJunctions_file() const      { return BaselevelJunctions_file; }
    /// @return ChannelSegments_file
    string get_ChannelSegments_file() const      { return ChannelSegments_file; }
    /// @return Floodplain_file
    string get_Floodplain_file() const        { return Floodplain_file; }
    /// @return the boundary conditions
    vector<string> get_boundary_conditions() const     { return boundary_conditions; }
    /// @return the float parameters
    map<string,float> get_float_parameters() const     { return float_parameters; }
    /// @return the float parameters
    map<string,int> get_int_parameters() const     { return int_parameters; }
    /// @return the float parameters
    map<string,bool> get_bool_parameters() const     { return bool_parameters; }
    /// @return the float parameters
    map<string,string> get_string_parameters() const     { return string_parameters; }
    /// @return the double parameters
    map<string,double> get_double_parameters() const     { return double_parameters; }

    /// @brief set the read_fname parameter with a string
    /// @param new_read_fname string containing the new read filename
    ///@author MDH
    ///@date 13/10/2017
    void set_read_fname(string new_read_fname) { read_fname = new_read_fname; }

    /// @brief set the read_fname parameter with a string
    /// @param new_read_fname string containing the new read filename
    ///@author MDH
    ///@date 13/10/2017
    void set_write_fname(string new_write_fname) { write_fname = new_write_fname; }

    /// @brief This checks to see if boundary condtions have been assigned and
    /// if not defaults to no flux boundaries
    /// @author SMM
    /// @date 02/11/2016
    void check_boundary_conditions();

    /// @brief This check to see if the filenames, paths and extensions have
    /// been assigned. If not it changes these to defaults
    /// @author SMM
    /// @date 02/11/2016
    void check_file_extensions_and_paths();


    /// @brief A routine that checks to see if the filename has a path
    ///  If it doesn't, adds the read path
    /// @param this_string The string to check for the path
    /// @author SMM
    /// @date 10/10/2018
    string check_for_path_and_add_read_path_if_required(string this_string);

    /// @brief this returns the string before the last dot in a string.
    /// so for example if you gave it paramfile.param it would return paramfile
    /// @param this_string the string you want truncated
    /// @return the truncated string
    /// @author SMM
    /// @date 29/07/2014
    string get_string_before_dot(string this_string);

    /// @brief This prints your parameters to file so you can check if the
    ///  parameters have been ingested properly
    /// @author SMM
    /// @date 01/11/2016
    void print_parameters();

    /// @brief This prints your parameters to file so you can check if the
    ///  parameters have been ingested properly. You tell it what the file is called
    /// @param fname_prefix the name of the fname (with extension, no path)
    /// @author SMM
    /// @date 16/03/2017
    void print_parameters(string fname_prefix);

    /// @brief This is used to update a parameter file within LSDTopoTools
    ///  and is primarily intended to be used in spawning operations,
    ///  e.g. where you select a bunch of catchemetns and create driver
    ///    functions for each one that can be run on different CPUs.
    /// @param parameter_fname the name of the output parameter file
    /// @param new_read_path the new path.
    /// @param new_read_fname the new read filename (no extension!)
    /// @param new_write_path says what it does on the tin
    /// @param new_write_fname no extension for this file
    /// @param replace_parameters a map of <string,string> that replaces parameters
    /// @author SMM
    /// @date 15/03/2017
    void replace_and_print_parameter_file(string parameter_fname,
                                     string new_read_path, string new_read_fname,
                                     string new_write_path, string new_write_fname,
                                     map<string,string> replace_parameters);

  protected:

    /// The path to the parameter file (used if no read or wirte path is supplied)
    string param_file_path;

    /// The path to the parameter file (used if no read or wirte path is supplied)
    string param_fname;

    /// The exteniosn of the rasters; either bil, asc or flt
    string dem_read_extension;

    /// The exteniosn of the rasters; either bil, asc or flt
    string dem_write_extension;

    /// Path to which files will be written
    string write_path;

    /// Path from which files will be written
    string read_path;

    /// Prefix of files to be written (i.e., no path, no extension)
    string write_fname;

    /// Prefix of files to be read (i.e., no path, no extension)
    string read_fname;

    /// The prefix of the channelheads filename
    string CHeads_file;

    /// The prefix of the baselevel junctions filename
    string BaselevelJunctions_file;

    /// The prefix of the channel segments filename
    string ChannelSegments_file;

    /// The prefix of the floodplain mask filename
    string Floodplain_file;

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    /// the four boundary conditions on the raster for the flow info object
    vector<string> boundary_conditions;

    /// This is the map of all the paramters read from file, in string format
    /// you use the default parameter maps to convert the strings to the
    /// parameter values you want.
    map<string,string> parameter_map;

    /// This is used for bug checking: all defaults used are then passed to this
    /// map, which can be printed and used to check the spelling of the default
    /// parameters
    map<string,string> defaults_used_map;

    /// Also for bug checking: contains the parameters where parameters have been
    /// read from the file
    map<string,string> parameters_read_map;

    /// This holds float parameters
    map<string,float> float_parameters;

    /// This holds integer parameters
    map<string,int> int_parameters;

    /// This holds string parameters
    map<string,string> string_parameters;

    /// This holds bool parameters
    map<string,bool> bool_parameters;

    /// This holds double parameters
    map<string,double> double_parameters;

  private:
    void create();
    void create(string PathName, string FileName);
};

#endif
