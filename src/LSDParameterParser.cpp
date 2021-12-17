//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDParameterParser.cpp
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
#include <fstream>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
#include "LSDParameterParser.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDParameterParser_CPP
#define LSDParameterParser_CPP

// This basic function is not part of the parameter parser object
// but instead is the function called at the beginning of any
// LSDTT driver call
vector<string> DriverIngestor(int nNumberofArgs,char *argv[])
{
  cout << "=========================================================" << endl;
  cout << "|| You have called an LSDTopoTools program.            ||" << endl;
  cout << "|| Prepare to explore topographic data!                ||" << endl;
  cout << "|| You can read the documentation at:                  ||" << endl;
  cout << "=========================================================" << endl;

  string path_name = ".";
  string file_name = "LSDTT_parameters.driver";

  //Test for correct input arguments
  if (nNumberofArgs == 1)
  {
    cout << "You have not given me any arguments. I am going to look" << endl;
    cout << "in this directory for a file with the extension .driver" << endl;
    cout << "I'll use the first one I find (in alphabetical ordering)." << endl;
    cout << "If I don't find one I am going to exit." << endl;
  }
  if (nNumberofArgs == 2)
  {
    cout << "I have one argument. I don't know if this is a directory path" << endl;
    cout << "or a driver filename. I am going to assume it is a directory path" << endl;
    cout << "if it containes the character . or /" << endl;

    string temp_arg = argv[1];
    string s_dot = ".";
    string sl_dot = "/";

    bool this_is_a_path = false;
    bool path_find = false;
    if (temp_arg.find(s_dot) != std::string::npos)
    {
      path_find = true;
    }
    if (temp_arg.find(sl_dot) != std::string::npos)
    {
      path_find = true;
    }

    if (this_is_a_path)
    {
      path_name = temp_arg;
    }
    else
    {
      file_name = temp_arg;
    }
  }
  //Test for correct input arguments
  if (nNumberofArgs==3)
  {
    cout << "I am reading the two arguments you gave me as the path name and the file name." << endl;
    path_name = argv[1];
    file_name = argv[2];
  }
  if (nNumberofArgs>=3)
  {
    cout << "You have provided more than two arguments. " << endl;
    cout << "I only expect 2. I am going to assume you meant" << endl;
    cout << "to give me the first two." << endl;
    path_name = argv[1];
    file_name = argv[2];

  }

  vector<string> path_and_file;
  path_and_file.push_back(path_name);
  path_and_file.push_back(file_name);

  if (path_and_file[0] == "-h" || path_and_file[0] == "--h" || 
      path_and_file[0] == "-H" || path_and_file[0] == "--H" || 
      path_and_file[0] == "--help" || path_and_file[0] == "-help"  ||
      path_and_file[0] == "help" || path_and_file[0] == "Help"  ||
      path_and_file[0] == "HELP" || path_and_file[0] == "-help"  ||
      path_and_file[0] == "-WTF" ||
      path_and_file[1] == "-h" || path_and_file[1] == "--h" || 
      path_and_file[0] == "-H" || path_and_file[0] == "--H" || 
      path_and_file[1] == "--help" || path_and_file[1] == "-help"  ||
      path_and_file[1] == "help" || path_and_file[1] == "Help"  ||
      path_and_file[1] == "HELP" || path_and_file[1] == "-help"  ||
      path_and_file[1] == "-WTF")
  {
    cout << "I am going to print a help file. There will be a .csv and a .html version. " << endl;
    cout << "These files have README in the filename." << endl;
    path_and_file[0] = "./";
    path_and_file[1] = "cry_for_help.txt";

    ofstream ofs;
    ofs.open("./cry_for_help.txt");
    ofs << "# The user has cried for help. " << endl;
    ofs << "# You can find the help file with README in the filename." << endl;
    ofs.close();
  }


  if (path_and_file[0] == "-c" || path_and_file[0] == "--c" || 
      path_and_file[0] == "--C" || path_and_file[0] == "-C"  ||
      path_and_file[0] == "cite" || path_and_file[0] == "Cite"  ||
      path_and_file[0] == "-cite" || path_and_file[0] == "-Cite"  ||
      path_and_file[0] == "--cite" || path_and_file[0] == "--Cite"  ||
      path_and_file[0] == "citation" || path_and_file[0] == "Citation"  ||
      path_and_file[0] == "-citation" || path_and_file[0] == "-Citation"  ||
      path_and_file[0] == "--citation" || path_and_file[0] == "--Citation"  ||
      path_and_file[1] == "-c" || path_and_file[0] == "--c" || 
      path_and_file[1] == "--C" || path_and_file[0] == "-C"  ||
      path_and_file[1] == "cite" || path_and_file[0] == "Cite"  ||
      path_and_file[1] == "-cite" || path_and_file[0] == "-Cite"  ||
      path_and_file[1] == "--cite" || path_and_file[0] == "--Cite"  ||
      path_and_file[1] == "citation" || path_and_file[0] == "Citation"  ||
      path_and_file[1] == "-citation" || path_and_file[0] == "-Citation"  ||
      path_and_file[1] == "--citation" || path_and_file[0] == "--Citation" )
  {
    //cout << "You have chosen to print the citation information. " << endl;
    path_and_file[0] = "./";
    path_and_file[1] = "lsdtt_citation.txt";

    ofstream ofs;
    ofs.open("./lsdtt_citation.txt");
    ofs << "# The user has asked for a citation. " << endl;
    ofs.close();
  }

  if (path_and_file[0] == "-v" || path_and_file[0] == "--v" || 
      path_and_file[0] == "-V" || path_and_file[0] == "--V"  ||
      path_and_file[0] == "version" || path_and_file[0] == "Version"  ||
      path_and_file[0] == "-version" || path_and_file[0] == "-Version"  ||
      path_and_file[0] == "--version" || path_and_file[0] == "--Version"  ||
      path_and_file[1] == "-v" || path_and_file[0] == "--v" || 
      path_and_file[1] == "-V" || path_and_file[0] == "--V"  ||
      path_and_file[1] == "version" || path_and_file[0] == "Version"  ||
      path_and_file[1] == "-version" || path_and_file[0] == "-Version"  ||
      path_and_file[1] == "--version" || path_and_file[0] == "--Version" )
  {
    //cout << "You have chosen to print the version. " << endl;
    path_and_file[0] = "./";
    path_and_file[1] = "lsdtt_version.txt";

    ofstream ofs;
    ofs.open("./lsdtt_version.txt");
    ofs << "# The user has asked for the version. " << endl;
    ofs.close();
  }


  return path_and_file;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create functions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::create()
{
  cout << "I have created an empty parameter parser object. " << endl;
  cout << "Surely you want to give it a filename?" << endl;
}

// This creates using a path and a filename
void LSDParameterParser::create(string PathName, string FileName)
{

  // Make sure the path has an extension
  PathName = FixPath(PathName);
  string FullName = PathName+FileName;

  param_file_path = PathName;
  param_fname = FileName;

  ifstream file_info_in;
  file_info_in.open(FullName.c_str());

  // check if the parameter file exists
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: The parameter file \"" << FullName
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  // now ingest the parameters
  cout << "Parsing the file" << endl;
  LSDPP_parse_file_into_parameter_map(FullName);
  parse_file_IO();

  // make sure the files are okay
  check_boundary_conditions();
  check_file_extensions_and_paths();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Forces the raster extensions to be bil
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParameterParser::force_bil_extension()
{
  cout << "===============================" << endl;
  cout << "WARNING!!! This program requires georeferencing so only ENVI bil format" << endl;
  cout << "Topographic data will be allowed!!" << endl;
  cout << "This is not the same as ESRI bil. In gdal use -of ENVI to ouput to ENVI bil" << endl;
  cout << "===============================" << endl;
  dem_read_extension = "bil";
  dem_write_extension = "bil";
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets a line of the parameter file. Has a long buffer so you can add long path
// names.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParameterParser::LSDPP_parse_line(ifstream &infile, string &parameter, string &value)
{
  char c;
  char buff[1024];
  int pos = 0;
  int word = 0;

  while ( infile.get(c) )
  {
    if (pos >= 1024)
    {
      cout << "Buffer overrun, word too long in parameter line: " << endl;
      string line;
      getline(infile, line);
      cout << "\t" << buff << " ! \n" << line << endl;
      exit(1);
    }
    // preceeding whitespace
    if (c == '#')
    {
      if (word == 0)
      {
        parameter = "NULL";
        value = "NULL";
      }
      if (word == 1)
        value = "NULL";
      word = 2;
    }

    if ((c == ' ' || c == '\t') && pos == 0)
      continue;
    else if ( (c == ':' && word == 0) || ( (c == ' ' || c == '\n' || c == '\t') && word == 1))
    {
      while (buff[pos-1] == ' ' || buff[pos-1] == '\t')
        --pos;    // Trailing whitespace
      buff[pos] = '\0';  // Append Null char
      if (word == 0)
        parameter = buff;  // Assign buffer contents
      else if (word == 1)
        value = buff;
      ++word;
      pos = 0;    // Rewind buffer
    }
    else if ( c == '\n' && word == 0 )
    {
      parameter = "NULL";
      buff[pos] = '\0';
      value = buff;
      ++word;
    }
    else if (word < 2)
    {
      buff[pos] = c;
      ++pos;
    }

    if (c == '\n')
      break;
  }
  if (word == 0)
  {
    parameter = "NULL";
    value = "NULL";
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This reads the parameter file, placing all parameters into a map
// with string values and string key. As you give the parameter parser
// default maps, it will scan these sting and convert them into the correct data type
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParameterParser::LSDPP_parse_file_into_parameter_map(string FullName)
{
  ifstream infile;
  infile.open(FullName.c_str());
  string parameter, value, lower, lower_val;
  string bc;

  cout << "Hello, I am going to parse your LSDTopoTools parameter file for you. " << endl;
  cout << "The parameter filename is: " << FullName << endl;

  // this will hold all the parameter values.
  map<string,string> temp_parameters;

  // now ingest parameters
  while (infile.good())
  {
    LSDPP_parse_line(infile, parameter, value);
    lower = parameter;
    //if (parameter == "NULL")
    //  continue;
    //for (unsigned int i=0; i<parameter.length(); ++i)
    //{
    //  lower[i] = tolower(parameter[i]);
    //}

    cout << "parameter is: " << lower << " and value is: " << value << endl;

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);

    temp_parameters[lower] = value;
  }

  parameter_map = temp_parameters;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This uses the parameter map to get file input and output
void LSDParameterParser::parse_file_IO()
{
  cout << endl << endl << endl << "----------------------" << endl;
  cout << "Parsing the file I/O" << endl;

  bool found_cheads = false;
  if(parameter_map.find("dem read extension") != parameter_map.end())
  {
    dem_read_extension = parameter_map["dem read extension"];
    // get rid of any control characters from the end (if param file was made in DOS)
    dem_read_extension = RemoveControlCharactersFromEndOfString(dem_read_extension);
  }
  if(parameter_map.find("dem write extension") != parameter_map.end())
  {
    dem_write_extension = parameter_map["dem write extension"];
    // get rid of any control characters from the end (if param file was made in DOS)
    dem_write_extension = RemoveControlCharactersFromEndOfString(dem_write_extension);
  }
  if(parameter_map.find("write path") != parameter_map.end())
  {
    write_path = parameter_map["write path"];
    // get rid of any control characters from the end (if param file was made in DOS)
    write_path = RemoveControlCharactersFromEndOfString(write_path);
  }
  if(parameter_map.find("write fname") != parameter_map.end())
  {
    write_fname = parameter_map["write fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    write_fname = RemoveControlCharactersFromEndOfString(write_fname);
    //cout << "Got the write name, it is: "  << write_fname << endl;
  }
  if(parameter_map.find("read path") != parameter_map.end())
  {
    read_path = parameter_map["read path"];
    // get rid of any control characters from the end (if param file was made in DOS)
    read_path = RemoveControlCharactersFromEndOfString(read_path);
    //cout << "Got the write name, it is: "  << write_fname << endl;
  }
  else
  {
    cout << "I did not find a read path so I am assuming the file is in this current directory." << endl;
    read_path = "./";
  }

  if(parameter_map.find("read fname") != parameter_map.end())
  {
    read_fname = parameter_map["read fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    read_fname = RemoveControlCharactersFromEndOfString(read_fname);
    //cout << "Got the write name, it is: "  << write_fname << endl;
  }

  if(parameter_map.find("CHeads_file") != parameter_map.end())
  {
    cout << endl << endl << endl;
    cout << "===========================" << endl;
    cout << "I have a channel heads file using the parameter string CHeads_file." << endl;
    cout << "If you have another parameter channel heads fname it will overwrite this parameter." << endl;
    cout << "This option is used for backwards compatibility but the channel heads fname is preferred." << endl;
    cout << "===========================" << endl;
    cout << endl << endl << endl;

    CHeads_file = parameter_map["CHeads_file"];
    // get rid of any control characters from the end (if param file was made in DOS)
    CHeads_file = RemoveControlCharactersFromEndOfString(CHeads_file);
    cout << "Got the channel heads name, it is: " << CHeads_file << endl;
    if(found_cheads)
    {
      cout << "Warning, channel head file is being overwritten--you have it twice in parameter file." << endl;
    }
    found_cheads = true;
  }



  if(parameter_map.find("channel heads fname") != parameter_map.end())
  {
    cout << "I found a channel heads fname which is" << parameter_map["channel heads fname"] << endl;
    CHeads_file = parameter_map["channel heads fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    CHeads_file = RemoveControlCharactersFromEndOfString(CHeads_file);
    //cout << "Got the channel heads name, it is: " << CHeads_file << endl;
    if(found_cheads)
    {
      cout << "Warning, channel head file is being overwritten--you have it twice in parameter file." << endl;
    }
    found_cheads = true;
  }
  else
  {
    if(found_cheads == false)
    {
      CHeads_file = "NULL";
    }
    else
    {
      cout << "I didn't find a channel heads fname but I have a previous channel heads from the CHeads_file" << endl;
      cout << "I will use that. It is: " << CHeads_file << endl;
    }
  }



  if(parameter_map.find("baselevel junctions fname") != parameter_map.end())
  {
    BaselevelJunctions_file = parameter_map["baselevel junctions fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    BaselevelJunctions_file = RemoveControlCharactersFromEndOfString(BaselevelJunctions_file);
  }
  else
  {
    BaselevelJunctions_file = "NULL";
  }
  if(parameter_map.find("channel segments fname") != parameter_map.end())
  {
    ChannelSegments_file = parameter_map["channel segments fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    ChannelSegments_file = RemoveControlCharactersFromEndOfString(ChannelSegments_file);
  }
  else
  {
    ChannelSegments_file = "NULL";
  }
  if(parameter_map.find("floodplain fname") != parameter_map.end())
  {
    Floodplain_file = parameter_map["floodplain fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    Floodplain_file = RemoveControlCharactersFromEndOfString(Floodplain_file);
  }
  else
  {
    Floodplain_file = "NULL";
  }
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This forces parsing of everything in the file
// Used for copying parameter files
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::force_parse()
{
  cout << "Forcing parsing of the parameter file" << endl;
  for( map<string, string >::iterator it = parameter_map.begin(); it != parameter_map.end(); ++it)
  {
    string key = it->first;
    cout << "Key is: " <<it->first << "\n";
    if(key != "CHeads file" && key != "read fname" && key != "write fname" &&
        key != "read path" && key != "write path" && key != "read extension" &&
        key != "write extension" && key != "NULL")
    {
      parameters_read_map[it->first] = it->second;
    }
  }

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses all the default parameter maps.
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::parse_all_parameters(map<string,float> default_map_f,
                      map<string,int> default_map_i, map<string,bool> default_map_b,
                      map<string,string> default_map_s)
{
  parse_float_parameters(default_map_f);
  parse_int_parameters(default_map_i);
  parse_bool_parameters(default_map_b);
  parse_string_parameters(default_map_s);

}




//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses all the default parameter maps.
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::parse_all_parameters(map<string,float> default_map_f,
                      map<string,int> default_map_i, map<string,bool> default_map_b,
                      map<string,string> default_map_s,
                      map<string,double> default_map_d)
{
  parse_float_parameters(default_map_f);
  parse_int_parameters(default_map_i);
  parse_bool_parameters(default_map_b);
  parse_string_parameters(default_map_s);
  parse_double_parameters(default_map_d);

}


void LSDParameterParser::print_help(map< string, vector<string> > help_map, string help_prefix, string version, string citation)
{
  
  // Open the help file
  ofstream help_ofs;
  string help_fname = write_path+help_prefix+".csv";

  help_ofs.open(help_fname.c_str());
  vector<string> this_string_vec;

  help_ofs << "name,type,default,description,guidance" << endl;
  help_ofs <<"version,"<< version << ",citation,"<< citation << ", please also look at statements printed by code for citations to specific routines." << endl;

  // get the parameters
  vector<string> these_keys = extract_keys(parameters_read_map);
  //cout << "The number of read keys are: " << these_keys.size() << endl;

  vector<string> more_keys = extract_keys(defaults_used_map);
  //cout << "The number of defaults are: " << more_keys.size() << endl;

  these_keys.insert( these_keys.end(), more_keys.begin(), more_keys.end() );
  //cout << "The updated size is: " << these_keys.size() << endl;
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    //cout << "The key is: " <<  these_keys[i] << endl;
    if ( help_map.find( these_keys[i] ) == help_map.end() ) 
    {
      help_ofs << these_keys[i] << ",not available,not available,not available,not available" << endl;
    } 
    else 
    {
      this_string_vec = help_map[ these_keys[i] ];
      help_ofs << these_keys[i] << "," << this_string_vec[0] << ","
               << this_string_vec[1] << "," << this_string_vec[2] << ","
               << this_string_vec[3] << endl;
    }

  }

  help_ofs.close();


  // Now we do the html file
  ofstream help_ofs_html;
  string help_fname_html = write_path+help_prefix+".html";
  help_ofs_html.open(help_fname_html.c_str());

  cout << "Printing the help files to: " << endl;
  cout << help_fname << endl;
  cout << help_fname_html << endl;

  help_ofs_html << "<!DOCTYPE html>\n<html>\n<head>\n\t<title>Help information for " << help_prefix << "</title>";
  help_ofs_html << "<style>table, th, td {\n\tborder: 1px solid black;\n\tborder-collapse: collapse;\n}\n</style>\n</head>\n<body>" << endl << endl;

  help_ofs_html << "<h1>Help information for " << help_prefix << "</h1>" << endl << endl;
  help_ofs_html << "<p>You are using version " << version << " of this software.</p>" << endl; 
  help_ofs_html << "<p>If the use of this software results in a publication, please cite <a href=\"" << citation << "\">" << citation << "</a>. In addition please also look at statements printed by the code for citations to specific routines and algorithms.</p>" << endl; 

  help_ofs_html << "<h2>Making the parameter file</h2>" << endl;
  help_ofs_html << "<p>You call lsdtt command line tools with a parameter file. This file can be made in a text editor. ";
  help_ofs_html << "The programs read ENVI bil format rasters projected into UTM (the WGS84 version). You must convert and project your DEM before using these tools. "; 
  help_ofs_html << "You should have a line of your parameter file that reads:.</p>" << endl; 
  help_ofs_html << "<p>read fname: FILE_PREFIX</p>" << endl; 
  help_ofs_html << "<p>where FILE_PREFIX is the name of the raster without the .bil extension</p>" << endl;
  help_ofs_html << "<p>The rest of the parameter file has the format:</p>\n<p>PARAMETER_NAME: VALUE</p><p>That is, there needs to be a colon and a space after the parameter name." << endl;
  help_ofs_html << "The default parameter value will be used if you don't give the parameter file the name and value of a parameter</p>" << endl << endl;

  // And now for the table. 
  // First the header
  help_ofs_html << "<h2>Parameter table</h2>" << endl;
  help_ofs_html << "<table>\n\t<colgroup>\n\t\t<col span=\"1\" style=\"width: 15%;\">"
                << "\n\t\t<col span=\"1\" style=\"width: 10%;\">"
                << "\n\t\t<col span=\"1\" style=\"width: 15%;\">"
                << "\n\t\t<col span=\"1\" style=\"width: 30%;\">"
                << "\n\t\t<col span=\"1\" style=\"width: 30%;\">"
                << "\n\t</colgroup>" << endl;
  
  help_ofs_html << "\n\t<tbody>\n\t<tr>\n\t\t<th>param name</th>\n\t\t<th>type</th>\n\t\t<th>default</th>\n\t\t<th>description</th>\n\t\t<th>guidance</th>\n\t</tr>" << endl;

  // And now for the rest of it. 
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    //cout << "The key is: " <<  these_keys[i] << endl;
    if ( help_map.find( these_keys[i] ) == help_map.end() ) 
    {
      help_ofs_html << "\t<tr>\n\t\t<td>" << these_keys[i] << "</td>\n\t\t<td>not available</td>\n\t\t<td>not available</td>\n\t\t<td>not available</td>\n\t\t<td>not available</td>\n\t</tr>" << endl;
    } 
    else 
    {
      this_string_vec = help_map[ these_keys[i] ];

      help_ofs_html << "\t<tr>\n\t\t<td>" << these_keys[i] << "</td>\n\t\t<td>" << this_string_vec[0] << "</td>\n\t\t<td>" 
                    << this_string_vec[1] << "</td>\n\t\t<td>" << this_string_vec[2] << "</td>\n\t\t<td>"
                    << this_string_vec[3] <<"</td>\n\t</tr>" << endl;
    }

  }

  // and now the end
  help_ofs_html << "\n\t<tbody>\n\t</table>\n\t</body>\n</html>" << endl;

  help_ofs_html.close();

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// These two functions takes a map of defualt parameters and returns the parameters for the
// current implementation
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::parse_float_parameters(map<string,float> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    //cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      //cout << "Found key " << these_keys[i];

      // convert the value to float
      float_parameters[these_keys[i]] = atof(parameter_map[these_keys[i]].c_str());
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
      //cout << " it is: " << parameter_map[these_keys[i]] << " check: " << float_parameters[these_keys[i]] << endl;

    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      float_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = dtoa(default_map[these_keys[i]]);

    }

  }
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Parse double parameters (for coords) FJC 20/11/17
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::parse_double_parameters(map<string,double> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    //cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      //cout << "Found key " << these_keys[i];

      // convert the value to float
      double_parameters[these_keys[i]] = atof(parameter_map[these_keys[i]].c_str());
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
      //cout << " it is: " << parameter_map[these_keys[i]] << " check: " << double_parameters[these_keys[i]] << endl;

    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      double_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = dtoa(default_map[these_keys[i]]);

    }

  }
}

void LSDParameterParser::parse_int_parameters(map<string,int> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    //cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      //cout << "Found the key" << endl;
      // convert the value to float
      int_parameters[these_keys[i]] = atoi(parameter_map[these_keys[i]].c_str());
      //cout << " it is: " << parameter_map[these_keys[i]] << " check: " << int_parameters[these_keys[i]] << endl;
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      int_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = itoa(default_map[these_keys[i]]);
    }
  }
}


void LSDParameterParser::parse_bool_parameters(map<string,bool> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    //cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      // convert the value to bool
      string value = parameter_map[these_keys[i]];
      bool temp_bool = (value == "true" || value== "True" || value == "TRUE" || value== "T" || value== "t") ? true : false;
      bool_parameters[these_keys[i]] = temp_bool;
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      bool_parameters[these_keys[i]] = default_map[these_keys[i]];
      if (default_map[these_keys[i]] == true)
      {
        defaults_used_map[these_keys[i]] = "true";
      }
      else
      {
        defaults_used_map[these_keys[i]] = "false" ;
      }

    }
  }
}

void LSDParameterParser::parse_string_parameters(map<string,string> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    //cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      string_parameters[these_keys[i]] = parameter_map[these_keys[i]];
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      string_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = default_map[these_keys[i]];
    }
  }
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses a vector of strings
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<string> LSDParameterParser::parse_string_vector(string key)
{
  string this_string = string_parameters[key];

  // reset the string vec
  vector<string> this_string_vec;

  // create a stringstream
  stringstream ss(this_string);

  // import the data, using a comma to separate
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_string_vec.push_back( substr );
  }

  return this_string_vec;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses a vector of ints
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDParameterParser::parse_int_vector(string key)
{
  string this_string = string_parameters[key];

  // reset the string vec
  vector<int> this_int_vec;

  // create a stringstream
  stringstream ss(this_string);

  // import the data, using a comma to separate
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_int_vec.push_back( atoi(substr.c_str()) );
  }

  return this_int_vec;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses a vector of ints
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDParameterParser::parse_float_vector(string key)
{
  cout << "I am going to parse a float vector for you!" << endl;

  string this_string = string_parameters[key];

  // reset the string vec
  vector<float> this_float_vec;

  // create a stringstream
  stringstream ss(this_string);

  // import the data, using a comma to separate
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_float_vec.push_back( atof(substr.c_str()) );
  }

  return this_float_vec;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This checks boundary conditions  for flow info
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::check_boundary_conditions()
{
  if( int(boundary_conditions.size()) != 4)
  {
    //cout << "Boundary conditions not assigned! Defaulting to no flux."  << endl;
    vector<string> temp_bc(4);
    for (int i = 0; i< 4; i++)
    {
      temp_bc[i] = "n";
    }
    boundary_conditions = temp_bc;
  }

  //for (int i =0; i< 4; i++)
  //{
  //  cout << "Boundary["<<i<<"]: "<<boundary_conditions[i]<< endl;
  //}
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function checks the file extensions for reading and writing DEMs
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::check_file_extensions_and_paths()
{

  // first check the extensions
  if (dem_read_extension != "asc"  && dem_read_extension != "flt" && dem_read_extension != "bil" &&
      dem_write_extension != "asc"  && dem_write_extension != "flt" && dem_write_extension != "bil")
  {
    //cout << "LSDParameterParser Raster file extension not assigned! Defaulting to bil format." << endl;
    //cout << "You entered: " << dem_read_extension << "!" <<endl;
    dem_read_extension = "bil";
    dem_write_extension = "bil";
  }
  else
  {
    if (dem_read_extension != "asc"  && dem_read_extension != "flt" && dem_read_extension != "bil")
    {
      //cout << "DEM read extension not assigned, defaulting to write extension." << endl;
      dem_read_extension = dem_write_extension;
    }
    else
    {
      //cout << "DEM write extension not assigned, defaulting to read extension." << endl;
      dem_write_extension = dem_read_extension;
    }
  }

  // now check the paths
  //cout << "Write path length is: " << write_path.length() << endl;
  if (write_path.length() == 0)
  {
    write_path = param_file_path;
    if (read_path.length() != 0)
    {
      write_path = read_path;
    }
  }

  //cout << "CHECKING NAMES, Write fname is: " << write_fname << endl;
  //cout << "The write fname length is " << write_fname.length() << endl;
  if (write_fname.length() == 0)
  {
    if (read_fname.length() != 0)
    {
      write_fname = read_fname;
    }
    write_fname = get_string_before_dot(param_fname);
    //cout << "Write fname not assigned, defaulting to name of parameter file." << endl;
    //cout << "The write fname is: " << write_fname << endl;
  }

  // now check the path
  //cout << "Read path length is: " << read_path.length() << endl;
  if (read_path.length() == 0)
  {
    read_path = write_path;
  }
  if (read_fname.length() == 0)
  {
    read_fname = get_string_before_dot(param_fname);
    //cout << "Read fname not assigned, defaulting to name of parameter file." << endl;
    //cout << "The read fname is: " << read_fname << endl;
  }

  // make sure the read and write paths have the slash at the end
  write_path = FixPath(write_path);
  read_path = FixPath(read_path);

  cout << "The full read fname is:\n " << read_path+read_fname << endl;
  cout << "The full write fname is:\n " << write_path+write_fname << endl;
  //cout << "The read and write extensions are: " << dem_read_extension
  //     << " " << dem_write_extension << endl;

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function checks filenames to see if they include a path. 
// If not it add the read path
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string LSDParameterParser::check_for_path_and_add_read_path_if_required(string this_string)
{
  string sl_dot = "/";
  string new_string;
  
  cout << "The string to check is: " << this_string << endl;
  if (this_string == "NULL" || this_string == "null" || this_string == "Null")
  {
    new_string = "NULL";
  }
  else
  {
    if (this_string.find(sl_dot) != std::string::npos)
    {
      cout << "This filename includes a path. I am not going to modify it." << endl;
      new_string = this_string;
    } 
    else
    {
      cout << "This finlename doesn't have a path. I am adding the read path." << endl;
      new_string = read_path+this_string;
      cout << "The new filename is: " << new_string << endl;
    }
  } 

  return new_string;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function strips the text after the final dot in a string
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string LSDParameterParser::get_string_before_dot(string this_string)
{
  string cut_string;
  unsigned found = this_string.find_last_of(".");
  cut_string = this_string.substr(0,found);
  return cut_string;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This prints parameters read to file, so you can make sure your parameters
// have ingested properly
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::print_parameters()
{
  string fname = write_path+write_fname+"_ingestedParam.param";
  ofstream params_out;
  params_out.open(fname.c_str());

  params_out << "# Here are the paramters ingested and set by the parameter parser" << endl;
  params_out << "# The file names and paths are: " << endl;
  params_out << "read path: " << read_path << endl;
  params_out << "read fname: " << read_fname << endl;
  params_out << "write path: " << write_path << endl;
  params_out << "write fname: " << write_fname << endl;
  params_out << "CHeads file: " << CHeads_file << endl;

  params_out << "read extension: " << dem_read_extension << endl;
  params_out << "write extension: " << dem_write_extension << endl << endl;

  params_out << "# ===================================="  << endl;
  params_out << "# Now for parameters read from file." << endl;
  params_out << "# If an expected parameter is not here check your spelling." << endl;

  vector<string> empty_vec;
  vector<string> these_keys = extract_keys(parameters_read_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << parameters_read_map[these_keys[i]]  << endl;
  }


  params_out << endl << "# ===================================="  << endl;
  params_out << "# Now for the default parameters." << endl;

  these_keys = empty_vec;
  these_keys = extract_keys(defaults_used_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << defaults_used_map[these_keys[i]]  << endl;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This prints parameters read to file, so you can make sure your parameters
// have ingested properly
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::print_parameters(string fname_prefix)
{
  string fname = write_path+fname_prefix;
  ofstream params_out;
  params_out.open(fname.c_str());

  params_out << "# Here are the paramters ingested and set by the parameter parser" << endl;
  params_out << "# The file names and paths are: " << endl;
  params_out << "read path: " << read_path << endl;
  params_out << "read fname: " << read_fname << endl;
  params_out << "write path: " << write_path << endl;
  params_out << "write fname: " << write_fname << endl;
  params_out << "CHeads file: " << CHeads_file << endl;

  params_out << "read extension: " << dem_read_extension << endl;
  params_out << "write extension: " << dem_write_extension << endl << endl;

  params_out << "# ===================================="  << endl;
  params_out << "# Now for parameters read from file." << endl;
  params_out << "# If an expected parameter is not here check your spelling." << endl;

  vector<string> empty_vec;
  vector<string> these_keys = extract_keys(parameters_read_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << parameters_read_map[these_keys[i]]  << endl;
  }


  params_out << endl << "# ===================================="  << endl;
  params_out << "# Now for the default parameters." << endl;

  these_keys = empty_vec;
  these_keys = extract_keys(defaults_used_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << defaults_used_map[these_keys[i]]  << endl;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This prints parameters read to file, so you can make sure your parameters
// have ingested properly
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParameterParser::replace_and_print_parameter_file(string parameter_fname,
                                     string new_read_path, string new_read_fname,
                                     string new_write_path, string new_write_fname,
                                     map<string,string> replace_parameters)
{
  string fname = write_path+parameter_fname;
  ofstream params_out;
  params_out.open(fname.c_str());

  params_out << "# This is an adjusted parameter file" << endl;
  params_out << "# The file names and paths are: " << endl;
  params_out << "read path: " << new_read_path << endl;
  params_out << "read fname: " << new_read_fname << endl;
  params_out << "write path: " << new_write_path << endl;
  params_out << "write fname: " << new_write_fname << endl;
  params_out << "CHeads file: " << CHeads_file << endl;

  params_out << "read extension: " << dem_read_extension << endl;
  params_out << "write extension: " << dem_write_extension << endl << endl;

  params_out << "# ===================================="  << endl;
  params_out << "# Now for parameters read from file." << endl;
  params_out << "# If an expected parameter is not here check your spelling." << endl;

  vector<string> empty_vec;
  vector<string> these_keys = extract_keys(parameters_read_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    // If the key is contained in the replace parameters, use the replace parameter
    if(replace_parameters.find(these_keys[i]) != replace_parameters.end())
    {
      cout << "I found a replace parameter!" << endl;
      params_out << these_keys[i] << ": " << replace_parameters[these_keys[i]]  << endl;
    }
    else
    {
      params_out << these_keys[i] << ": " << parameters_read_map[these_keys[i]]  << endl;
    }
  }

  params_out << endl << "# ===================================="  << endl;
  params_out << "# Now for the default parameters." << endl;

  these_keys = empty_vec;
  these_keys = extract_keys(defaults_used_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    // If the key is contained in the replace parameters, use the replace parameter
    if(replace_parameters.find(these_keys[i]) != replace_parameters.end())
    {
      cout << "I found a replace parameter!" << endl;
      params_out << these_keys[i] << ": " << replace_parameters[these_keys[i]]  << endl;
    }
    else
    {
      params_out << these_keys[i] << ": " << defaults_used_map[these_keys[i]]  << endl;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#endif
