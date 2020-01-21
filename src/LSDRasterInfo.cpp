//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterInfo
// Land Surface Dynamics Raster Information
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for examining the georeferencing and existence of rasters
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

#ifndef LSDRasterInfo_CPP
#define LSDRasterInfo_CPP

#include <fstream>
#include <string>
#include <map>
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDRasterInfo.hpp"
using namespace std;


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Basic create function, does nothing
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterInfo::create()
{}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Create function that reads the header of the file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterInfo::create(string filename, string extension)
{
  read_header(filename, extension);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Create function that reads from a raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterInfo::create(LSDRaster& Raster)
{
      ///Number of rows.
  NRows = Raster.get_NRows();
  NCols = Raster.get_NCols();
  XMinimum = Raster.get_XMinimum();
  YMinimum = Raster.get_YMinimum();
  DataResolution = Raster.get_DataResolution();
  NoDataValue = Raster.get_NoDataValue();
  GeoReferencingStrings = Raster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Create function that reads from an index raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterInfo::create(LSDIndexRaster& IRaster)
{
      ///Number of rows.
  NRows = IRaster.get_NRows();
  NCols = IRaster.get_NCols();
  XMinimum = IRaster.get_XMinimum();
  YMinimum = IRaster.get_YMinimum();
  DataResolution = IRaster.get_DataResolution();
  NoDataValue = IRaster.get_NoDataValue();
  GeoReferencingStrings = IRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// The equality operator. Used for testing if the rasters have the 
//  same dimensions and georeferencing
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterInfo::operator==(LSDRasterInfo& other) 
{
  map<string,string> DRS = other.get_GeoReferencingStrings();
  string mi_str_key = "ENVI_map_info";
  string cs_str_key = "ENVI_coordinate_system";
  
  bool is_georef_and_dimensions_same;

  if (NRows == other.get_NRows() &&
      NCols == other.get_NCols() &&
      XMinimum == other.get_XMinimum() &&
      YMinimum == other.get_YMinimum() &&
      DataResolution == other.get_DataResolution() &&
      NoDataValue == other.get_NoDataValue() &&
      GeoReferencingStrings[mi_str_key] == DRS[mi_str_key] &&
      GeoReferencingStrings[cs_str_key] == DRS[cs_str_key])
  {
    is_georef_and_dimensions_same = true;
  }
  else
  {
    is_georef_and_dimensions_same = false;
  }

  return is_georef_and_dimensions_same;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// The inequality operator. Used for testing if the rasters have the 
//  same dimensions and georeferencing
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterInfo::operator!=(LSDRasterInfo& other) 
{
  map<string,string> DRS = other.get_GeoReferencingStrings();
  string mi_str_key = "ENVI_map_info";
  string cs_str_key = "ENVI_coordinate_system";
  
  bool is_georef_and_dimensions_same;

  if (NRows == other.get_NRows() &&
      NCols == other.get_NCols() &&
      XMinimum == other.get_XMinimum() &&
      YMinimum == other.get_YMinimum() &&
      DataResolution == other.get_DataResolution() &&
      NoDataValue == other.get_NoDataValue() &&
      GeoReferencingStrings[mi_str_key] == DRS[mi_str_key] &&
      GeoReferencingStrings[cs_str_key] == DRS[cs_str_key])
  {
    is_georef_and_dimensions_same = false;
  }
  else
  {
    is_georef_and_dimensions_same = true;
  }

  return is_georef_and_dimensions_same;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function reads a DEM header
// One has to provide both the filename and the extension
// the '.' between the filename and extension is not included
// for example, if the full filename is test.asc
// then
// filename = "test"
// and
// ext = "asc"
// The full filename could also be "test.01.asc"
// so filename would be "test.01"
// and ext would again be "asc"
//
// SMM 2015, copied from the LSDRaster version
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterInfo::read_header(string filename, string extension)
{
  string string_filename;
  string dot = ".";
  string_filename = filename+dot+extension;
  //cout << "\n\nLoading LSDRasterInfo, the filename is " << string_filename << endl;

  if (extension == "asc")
  {
    // open the data file
    ifstream data_in(string_filename.c_str());

    if( data_in.fail() )
    {
      cout << "\nFATAL ERROR: the header file \"" << string_filename
         << "\" doesn't exist" << std::endl;
      exit(EXIT_FAILURE);
    }

    //Read in raster data
    string str;    // a temporary string for discarding text

    // read the georeferencing data and metadata
    data_in >> str >> NCols;
                cout << "NCols: " << NCols << " str: " << endl;
    data_in >> str >> NRows;
    //cout << "NRows: " << NRows << " str: " << endl;
    data_in >> str >> XMinimum >> str >> YMinimum
          >> str >> DataResolution
          >> str >> NoDataValue;

    //cout << "Loading asc file; NCols: " << NCols << " NRows: " << NRows << endl
    //     << "X minimum: " << XMinimum << " YMinimum: " << YMinimum << endl
    //     << "Data Resolution: " << DataResolution << " and No Data Value: "
    //     << NoDataValue << endl;
  }
  else if (extension == "flt")
  {
    // float data (a binary format created by ArcMap) has a header file
    // this file must be opened first
    string header_filename;
    string header_extension = "hdr";
    header_filename = filename+dot+header_extension;

    ifstream ifs(header_filename.c_str());
    if( ifs.fail() )
    {
      cout << "\nFATAL ERROR: the header file \"" << header_filename
         << "\" doesn't exist" << std::endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      string str;
      ifs >> str >> NCols;
      //cout << "NCols: " << NCols << " str: " << endl;
      ifs >> str >> NRows;
      //cout << "NRows: " << NRows << " str: " << endl;
      ifs >> str >> XMinimum >> str >> YMinimum
          >> str >> DataResolution
          >> str >> NoDataValue;
    }
    ifs.close();

    //cout << "Loading flt file; NCols: " << NCols << " NRows: " << NRows << endl
    //     << "X minimum: " << XMinimum << " YMinimum: " << YMinimum << endl
    //     << "Data Resolution: " << DataResolution << " and No Data Value: "
    //     << NoDataValue << endl;

  }
  else if (extension == "bil")
  {
    // float data (a binary format created by ArcMap) has a header file
    // this file must be opened first
    string header_filename;
    string header_extension = "hdr";
    header_filename = filename+dot+header_extension;
    int NoDataExists = 0;
    int DataType = 4;     // default is float data

    ifstream ifs(header_filename.c_str());
    if( ifs.fail() )
    {
      cout << "\nFATAL ERROR: the header file \"" << header_filename
           << "\" doesn't exist" << std::endl;
      exit(EXIT_FAILURE);
    }
    else
    {  
      string str;
      ifs >> str;
      if (str != "ENVI")
      {
        cout << "\nFATAL ERROR: this is not an ENVI header file!, first line is: " 
             << str << endl;
        exit(EXIT_FAILURE);       
      }
      else
      {
        // the the rest of the lines
        int NChars = 5000; // need a big buffer beacause of the projection string
        char* thisline = new char[NChars]; //char thisline[NChars]; MSVC/Clang uniformisation
        vector<string> lines;
        while( ifs.getline(thisline, NChars) )
        {
          lines.push_back(thisline);
        }
        
        // now loop through and get the number of rows
        int counter = 0;
        int NLines = int(lines.size());
        int this_NRows = 0;
        size_t found; 
        string str_find = "lines";
        while (counter < NLines)
        {
          found = lines[counter].find(str_find); 
          if (found!=string::npos)
          {
            // get the data using a stringstream
            istringstream iss(lines[counter]);
            iss >> str >> str >> str;
            this_NRows = atoi(str.c_str());
            //cout << "NRows = " << this_NRows << endl;
            NRows = this_NRows;
            
            // advance to the end so you move on to the new loop
            counter = lines.size();            
          }
          else
          {
            counter++;
          }
        }

        // get the number of columns
        counter = 0;
        int this_NCols = 0;
        str_find = "samples";
        while (counter < NLines)
        {
          found = lines[counter].find(str_find); 
          if (found!=string::npos)
          {
            // get the data using a stringstream
            istringstream iss(lines[counter]);
            iss >> str >> str >> str;
            this_NCols = atoi(str.c_str());
            //cout << "NCols = " << this_NCols << endl;
            NCols = this_NCols;
            
            // advance to the end so you move on to the new loop            
            counter = lines.size();    
          }
          else
          {
            counter++;
          }
        }  
                     
        // get data ignore value
        counter = 0;
        float this_NoDataValue = 0;
        str_find = "data ignore value";
        while (counter < NLines)
        {
          found = lines[counter].find(str_find); 
          if (found!=string::npos)
          {
            // get the data using a stringstream
            istringstream iss(lines[counter]);
            iss >> str >> str >> str >> str >> str;
            this_NoDataValue = atoi(str.c_str());
            //cout << "NCols = " << this_NCols << endl;
            NoDataValue = this_NoDataValue;
            
            NoDataExists = 1;   // set this to true
            
            // advance to the end so you move on to the new loop            
            counter = lines.size();    
          }
          else
          {
            counter++;
          }
        }   

        // get data type
        counter = 0;

        str_find = "data type";
        while (counter < NLines)
        {
          found = lines[counter].find(str_find); 
          if (found!=string::npos)
          {
            // get the data using a stringstream
            istringstream iss(lines[counter]);
            iss >> str >> str >> str >> str >> str;
            DataType = atoi(str.c_str());
            cout << "Data Type = " << DataType << endl;
                     
            // advance to the end so you move on to the new loop            
            counter = lines.size();    
          }
          else
          {
            counter++;
          }
        }  
        
        // get the map info
        counter = 0;
        string this_map_info = "empty";
        str_find = "map info";
        while (counter < NLines)
        {
          found = lines[counter].find(str_find); 
          if (found!=string::npos)
          {
            //cout << "Found map info on line " << counter << '\n';  
      
            // now split the line 
            size_t start_pos;
            size_t end_pos;
            string open_curly_bracket = "{";
            string closed_curly_bracket = "}";
            start_pos = lines[counter].find(open_curly_bracket);
            end_pos = lines[counter].find(closed_curly_bracket);
            //cout << "startpos: " << start_pos << " and end pos: " << end_pos << endl;
            string info_str = lines[counter].substr(start_pos+1, end_pos-start_pos-1);
            //cout << "\nThe map info string is:\n" << info_str << endl;
            string mi_key = "ENVI_map_info";
            GeoReferencingStrings[mi_key] = info_str;

            // now parse the string
            vector<string> mapinfo_strings;
            istringstream iss(info_str);
            while( iss.good() )
            {
              string substr;
              getline( iss, substr, ',' );
              mapinfo_strings.push_back( substr );
            }
            XMinimum = atof(mapinfo_strings[3].c_str());	          	          
            float YMax = atof(mapinfo_strings[4].c_str());
            	       	   	          
            DataResolution = atof(mapinfo_strings[5].c_str());
            
            // get Y minium
            YMinimum = YMax - NRows*DataResolution;	          
            
            //using a string comparison as float(X) != float(X) in many cases due to floating point math
            // http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm  - SWDG	          
            if (mapinfo_strings[5] != mapinfo_strings[6]) 
            {
              cout << "Warning! Loading ENVI DEM, but X and Y data spacing are different!" << endl;
            }

            //cout << "Xmin: " << XMinimum << " YMin: " << YMinimum << " spacing: " 
            //     << DataResolution << endl;

            counter = lines.size();
          }
          else
          {
            counter++;
          }
        }     

        // get the projection string
        counter = 0;
        string this_coordinate_system_string = "empty";
        str_find = "coordinate system string";
        while (counter < NLines)
        {
          found = lines[counter].find(str_find); 
          if (found!=string::npos)
          {
            //cout << "Found coordinate system string on line " << counter << '\n';  

            // now split the line 
            size_t start_pos;
            size_t end_pos;
            string open_curly_bracket = "{";
            string closed_curly_bracket = "}";
            start_pos = lines[counter].find(open_curly_bracket);
            end_pos = lines[counter].find(closed_curly_bracket);
            //cout << "startpos: " << start_pos << " and end pos: " << end_pos << endl;
            string csys_str = lines[counter].substr(start_pos+1, end_pos-start_pos-1);
            //cout << "\nThe coordinate system string is:\n" << csys_str << endl;
            string cs_key = "ENVI_coordinate_system";
            GeoReferencingStrings[cs_key] = csys_str;
            counter = lines.size();
          }
          else
          {
            counter++;
          }
        }          
      }         
    }
    ifs.close(); 
     
    // this is the array into which data is fed
    if (NoDataExists == 0)
    {
      NoDataValue = -9999;
      
      if(extension == "bil")
      {
        string header_filename;
        string header_extension = "hdr";
        header_filename = filename+dot+header_extension;
        ofstream file_out;
        file_out.open(header_filename.c_str(), std::ios_base::app);
        file_out << "data ignore value = -9999" << endl;
        file_out.close(); 
      }
      
    }
    
    

    //cout << "Loading ENVI bil file; NCols: " << NCols << " NRows: " << NRows << endl
    //   << "X minimum: " << XMinimum << " YMinimum: " << YMinimum << endl
    //     << "Data Resolution: " << DataResolution << " and No Data Value: "
    //     << NoDataValue << endl;

  }
  else
  {
    cout << "You did not enter and approprate extension!" << endl
          << "You entered: " << extension << " options are .flt, .asc and .bil" << endl;
    exit(EXIT_FAILURE);
  }


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDRasterInfo::get_UTM_information(int& UTM_zone, bool& is_North)
{

  // set up strings and iterators
  map<string,string>::iterator iter;

  //check to see if there is already a map info string
  string mi_key = "ENVI_map_info";
  iter = GeoReferencingStrings.find(mi_key);
  if (iter != GeoReferencingStrings.end() )
  {
    string info_str = GeoReferencingStrings[mi_key] ;

    // now parse the string
    vector<string> mapinfo_strings;
    istringstream iss(info_str);
    while( iss.good() )
    {
      string substr;
      getline( iss, substr, ',' );
      mapinfo_strings.push_back( substr );
    }
    UTM_zone = atoi(mapinfo_strings[7].c_str());
    //cout << "Line 1041, UTM zone: " << UTM_zone << endl;
    //cout << "LINE 1042 LSDRaster, N or S: " << mapinfo_strings[7] << endl;

    // find if the zone is in the north
    string n_str = "n";
    string N_str = "N";
    is_North = false;
    size_t found = mapinfo_strings[8].find(N_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    //cout << "is_North is: " << is_North << endl;

  }
  else
  {
    UTM_zone = NoDataValue;
    is_North = false;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Checks to see is a point is in the raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDRasterInfo::check_if_point_is_in_raster(float X_coordinate,float Y_coordinate)
{
  bool is_in_raster = true;

  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1)
  {
    is_in_raster = false;
  }

  return is_in_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


void LSDRasterInfo::print_raster_information()
{
  cout << "NRows: " << NRows << endl;
  cout << "NCols: "  << NCols << endl;
  cout << "XMinimum: "  << XMinimum << endl;
  cout << "YMinimum: "  << YMinimum << endl;
  cout << "DataResolution: " << DataResolution << endl;
  cout << "NoDataValue: " << NoDataValue << endl;
  cout << "ENVI_map_info: " << GeoReferencingStrings["ENVI_map_info"] << endl;
  cout << "ENVI_coordinate_system: " << GeoReferencingStrings["ENVI_coordinate_system"] << endl;
}

#endif
