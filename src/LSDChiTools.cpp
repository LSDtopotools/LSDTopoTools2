//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiTools
// Land Surface Dynamics ChiTools object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for performing various analyses in chi space
//
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//  Boris Gailleton
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDChiTools.cpp
// LSDChiTools object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDChiTools_CPP
#define LSDChiTools_CPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <queue>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDChiTools.hpp"
#include "LSDBasin.hpp"
#include "LSDChiNetwork.hpp"
#include "LSDMostLikelyPartitionsFinder.hpp"

// Testing some timings here
#include <chrono>

#ifdef _OPENMP 
  #include <thread>
#endif

using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDRaster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDRaster& ThisRaster)
{
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDRaster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDIndexRaster& ThisRaster)
{
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDFlowInfo
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDFlowInfo& ThisFI)
{
  NRows = ThisFI.get_NRows();
  NCols = ThisFI.get_NCols();
  XMinimum = ThisFI.get_XMinimum();
  YMinimum = ThisFI.get_YMinimum();
  DataResolution = ThisFI.get_DataResolution();
  NoDataValue = ThisFI.get_NoDataValue();
  GeoReferencingStrings = ThisFI.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDFlowInfo
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDJunctionNetwork& ThisJN)
{
  NRows = ThisJN.get_NRows();
  NCols = ThisJN.get_NCols();
  XMinimum = ThisJN.get_XMinimum();
  YMinimum = ThisJN.get_YMinimum();
  DataResolution = ThisJN.get_DataResolution();
  NoDataValue = ThisJN.get_NoDataValue();
  GeoReferencingStrings = ThisJN.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// the create function. This is default and throws an error
// SMM 2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChiTools::create()
{
  // cout << "LSDChiTools line 64 Warning you have an empty LSDChiTools!" << endl;
  // exit(EXIT_FAILURE);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This resets all the data maps
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::reset_data_maps()
{
  map<int,float> empty_map;
  vector<int> empty_vec;

  M_chi_data_map = empty_map;
  b_chi_data_map = empty_map;
  elev_data_map = empty_map;
  chi_data_map = empty_map;
  flow_distance_data_map = empty_map;
  drainage_area_data_map = empty_map;
  node_sequence = empty_vec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc)
{

  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;

  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
// Same as above but with floats
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc)
{

  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;

  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to convert a node position with a row and column to a lat
// and long coordinate
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_lat_and_long_locations(int row, int col, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter)
{
  // get the x and y locations of the node
  double x_loc,y_loc;
  get_x_and_y_locations(row, col, x_loc, y_loc);

  // get the UTM zone of the node
  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);
  //cout << endl << endl << "Line 1034, UTM zone is: " << UTM_zone << endl;


  if(UTM_zone == NoDataValue)
  {
    lat = NoDataValue;
    longitude = NoDataValue;
  }
  else
  {
    // set the default ellipsoid to WGS84
    int eId = 22;

    double xld = double(x_loc);
    double yld = double(y_loc);

    // use the converter to convert to lat and long
    double Lat,Long;
    Converter.UTMtoLL(eId, yld, xld, UTM_zone, is_North, Lat, Long);


    lat = Lat;
    longitude = Long;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to convert a x/y position that is not necessarly a node into lat/long
// and long coordinate
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_lat_and_long_locations_from_coordinate(float X, float Y, double& lat,
                   double& longitude, LSDCoordinateConverterLLandUTM Converter)
{

  // get the UTM zone of the node
  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);
  //cout << endl << endl << "Line 1034, UTM zone is: " << UTM_zone << endl;


  if(UTM_zone == NoDataValue)
  {
    lat = NoDataValue;
    longitude = NoDataValue;
  }
  else
  {
    // set the default ellipsoid to WGS84
    int eId = 22;

    double xld = double(X);
    double yld = double(Y);

    // use the converter to convert to lat and long
    double Lat,Long;
    Converter.UTMtoLL(eId, yld, xld, UTM_zone, is_North, Lat, Long);


    lat = Lat;
    longitude = Long;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_UTM_information(int& UTM_zone, bool& is_North)
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
    if (found!=string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=string::npos)
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
// This updates the chi values using a new chi raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::update_chi_data_map(LSDFlowInfo& FlowInfo, LSDRaster& Chi_coord)
{
  if (chi_data_map.size() == 0)
  {
    cout << "Trying to update chi but you have not run the automator yet to" << endl;
    cout << "organise the sources and channels. LSDChiTools::update_chi_data_map" << endl;
  }
  else
  {
    int n_nodes = int(node_sequence.size());
    int this_node,row,col;
    float updated_chi;
    for(int node = 0; node<n_nodes; node++)
    {
      this_node = node_sequence[node];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      updated_chi = Chi_coord.get_data_element(row,col);
      chi_data_map[this_node] = updated_chi;
    }
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This updates the chi values by calculating them directly from the FlowInfo object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::update_chi_data_map(LSDFlowInfo& FlowInfo, float A_0, float movern)
{
  if (chi_data_map.size() == 0)
  {
    cout << "Trying to update chi but you have not run the automator yet to" << endl;
    cout << "organise the sources and channels. LSDChiTools::update_chi_data_map" << endl;
  }
  else
  {
    cout << "WARNING: update_chi_data_map" << endl;
    cout<< "Chi is being calculated from outlet nodes, so interior basins may have different starting chi values." << endl;

    float thresh_area_for_chi = 0;  // this gets chi in all nodes. Not much slower and avoids errors
    LSDRaster this_chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

    int n_nodes = int(node_sequence.size());
    int this_node,row,col;
    float updated_chi;
    for(int node = 0; node<n_nodes; node++)
    {
      this_node = node_sequence[node];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      updated_chi = this_chi_coordinate.get_data_element(row,col);
      chi_data_map[this_node] = updated_chi;
    }
  }
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This updates the chi values by calculating them directly from the FlowInfo object
// The outlet_node_from_basin_key_map map is generated by
// get_outlet_node_from_basin_key_map()
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::update_chi_data_map_for_single_basin(LSDFlowInfo& FlowInfo, float A_0, float movern,
                                     int minimum_contributing_pixels, int basin_key,
                                     map<int,int> outlet_node_from_basin_key_map)
{
  if (chi_data_map.size() == 0)
  {
    cout << "Trying to update chi but you have not run the automator yet to" << endl;
    cout << "organise the sources and channels. LSDChiTools::update_chi_data_map_for_single_basin" << endl;
  }
  else
  {
    // get the source node of the basin in question
    int outlet_node_index;
    if ( outlet_node_from_basin_key_map.find(basin_key) == outlet_node_from_basin_key_map.end() )
    {
      cout << "Fatal error: you have selected a basin key (" << basin_key << " that is not amongst the basins." << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      outlet_node_index = outlet_node_from_basin_key_map[basin_key];
    }

    // now get the chi coordinate updated for this basin only
    //cout << "Outlet node index is: " << outlet_node_index << endl;
    map<int,float> new_chi_map = FlowInfo.get_upslope_chi_from_single_starting_node(outlet_node_index , movern, A_0, minimum_contributing_pixels);

    map<int,float>::iterator iter = new_chi_map.begin();
    while(iter != new_chi_map.end())
    {
      //cout << "node is: " << iter->first << endl;
      chi_data_map[iter->first] = iter->second;
      iter++;
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This updates the chi values by calculating them directly from the FlowInfo object
// The outlet_node_from_basin_key_map map is generated by
// get_outlet_node_from_basin_key_map()
// same as above but includes discharge raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::update_chi_data_map_for_single_basin(LSDFlowInfo& FlowInfo, float A_0, float movern,
                                     int minimum_contributing_pixels, int basin_key,
                                     map<int,int> outlet_node_from_basin_key_map, LSDRaster& Discharge)
{
  if (chi_data_map.size() == 0)
  {
    cout << "Trying to update chi but you have not run the automator yet to" << endl;
    cout << "organise the sources and channels. LSDChiTools::update_chi_data_map_for_single_basin" << endl;
  }
  else
  {
    // get the source node of the basin in question
    int outlet_node_index;
    if ( outlet_node_from_basin_key_map.find(basin_key) == outlet_node_from_basin_key_map.end() )
    {
      cout << "Fatal error: you have selected a basin key (" << basin_key << " that is not amongst the basins." << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      outlet_node_index = outlet_node_from_basin_key_map[basin_key];
    }

    // now get the chi coordinate updated for this basin only
    map<int,float> new_chi_map = FlowInfo.get_upslope_chi_from_single_starting_node(outlet_node_index , movern, A_0, minimum_contributing_pixels, Discharge);

    map<int,float>::iterator iter = new_chi_map.begin();
    while(iter != new_chi_map.end())
    {
      chi_data_map[iter->first] = iter->second;
      iter++;
    }
  }
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname,
                                 float A_0, float m_over_n, float area_threshold)
{

  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());

  chi_map_csv_out.precision(9);

  float chi_coord;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  chi_map_csv_out << "latitude,longitude,chi" << endl;

  LSDRaster Chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0, area_threshold);

  float NDV = Chi.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      chi_coord =  Chi.get_data_element(row,col);

      if (chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out << latitude << "," << longitude  << "," << chi_coord << endl;
      }
    }
  }

  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2. You feed it the chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname,
                                 LSDRaster& chi_coord)
{

  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());



  float this_chi_coord;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  chi_map_csv_out << "latitude,longitude,chi" << endl;

  float NDV = chi_coord.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      this_chi_coord = chi_coord.get_data_element(row,col);

      if (this_chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out.precision(9);
        chi_map_csv_out << latitude << "," << longitude  << ",";
        chi_map_csv_out.precision(5);
        chi_map_csv_out << this_chi_coord << endl;
      }
    }
  }

  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2. You feed it the chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname,
                                 LSDRaster& chi_coord, LSDIndexRaster& basin_raster)
{

  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());



  float this_chi_coord;
  int this_basin_number;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  chi_map_csv_out << "latitude,longitude,chi,basin_junction" << endl;

  float NDV = chi_coord.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      this_chi_coord = chi_coord.get_data_element(row,col);
      this_basin_number = basin_raster.get_data_element(row,col);


      if (this_chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out.precision(9);
        chi_map_csv_out << latitude << "," << longitude  << ",";
        chi_map_csv_out.precision(5);
        chi_map_csv_out << this_chi_coord << ",";
        chi_map_csv_out << this_basin_number << endl;
      }
    }
  }

  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a csv with basin_ID, all the litho count and the total
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::simple_litho_basin_to_csv(LSDFlowInfo& FlowInfo, string csv_slbc_fname,
                                  map<int,map<int,int> > map_slbc)
{
  cout << "I am writing a simple basin/litho csv file: " << csv_slbc_fname << endl;
  // opening the stream to write the csv file
  ofstream csv_out;
  csv_out.open(csv_slbc_fname.c_str());
  // writing the header
  csv_out << "basin_id,";
  for(map<int,int>::iterator it = map_slbc[0].begin(); it !=map_slbc[0].end();++it)
  {
    csv_out << it->first << ",";
  }
  csv_out << "total"<< endl;
  // preparing the writing
  int total_temp = 0;
  map<int,int> tempmapcsv; // temporary map to avoid mapception confusions
  // writing, first loop through basins
  for(map<int,map<int,int> >::iterator it1 = map_slbc.begin(); it1 !=map_slbc.end();++it1)
  {
    csv_out << it1->first<<","; // writing the basin ID
    tempmapcsv = it1->second; // Getting the map of values
    // Second loop through the map of litho
    for(map<int,int>::iterator it = tempmapcsv.begin(); it !=tempmapcsv.end();++it)
    {
      csv_out << it->second << ","; // writing the value
      total_temp += it->second; // temporary temp
    }
    // writing and reinitializing the total
    csv_out << total_temp << endl;
    total_temp = 0;
  }
  // Done, closing the file stream
  csv_out.close();
  cout << "I am done writing the simple litho file." << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a csv with basin_ID, all the litho count and the total
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::extended_litho_basin_to_csv(LSDFlowInfo& FlowInfo, string csv_slbc_fname,
                                  map<int,map<int,int> > map_slbc)
{
  cout << "I am writing a simple basin/litho csv file: " << csv_slbc_fname << endl;
  // opening the stream to write the csv file
  ofstream csv_out;
  csv_out.open(csv_slbc_fname.c_str());
  // writing the header
  csv_out << "method," << "basin_id,";
  for(map<int,int>::iterator it = map_slbc[0].begin(); it !=map_slbc[0].end();++it)
  {
    csv_out << it->first << ",";
  }
  csv_out << "total"<< endl;

  // Some stats on it
  map<int,map<int,float> > percentmap;
  map<int,float> tempmappercent;
  map<int,int> tempmapcsv1;
  int total_temp;
  float tvalue;

  for(map<int,map<int,int> >::iterator it1 = map_slbc.begin(); it1 !=map_slbc.end();++it1)
  {
    total_temp = 0;
    // Second loop through the map of litho
    tempmapcsv1 = it1->second;
    for(map<int,int>::iterator it = tempmapcsv1.begin(); it !=tempmapcsv1.end();++it)
    {
      total_temp += it->second;
    }
    for(map<int,int>::iterator it = tempmapcsv1.begin(); it !=tempmapcsv1.end();++it)
    {
      tvalue = it->second;
      tempmappercent[it->first] = tvalue*100/float(total_temp);
    }
    percentmap[it1->first] = tempmappercent;
    tempmappercent.clear();
    total_temp  = 0;
  }


  // preparing the writing
  total_temp = 0;
  map<int,int> tempmapcsv; // temporary map to avoid mapception confusions
  // writing, first loop through basins
  for(map<int,map<int,int> >::iterator it1 = map_slbc.begin(); it1 !=map_slbc.end();++it1)
  {
    // writing the raw count
    csv_out << "count,";
    csv_out << it1->first<<","; // writing the basin ID
    tempmapcsv = it1->second; // Getting the map of values
    // Second loop through the map of litho
    for(map<int,int>::iterator it = tempmapcsv.begin(); it !=tempmapcsv.end();++it)
    {
      csv_out << it->second << ","; // writing the value
      total_temp += it->second; // temporary temp
    }
    // writing and reinitializing the total
    csv_out << total_temp << endl;
    total_temp = 0;

    // writing the percentage
    csv_out << "percentage,";
    csv_out << it1->first<<","; // writing the basin ID
    tempmappercent = percentmap[it1->first]; // Getting the map of values
    // Second loop through the map of litho
    float percentemptot =0;
    for(map<int,float>::iterator it = tempmappercent.begin(); it !=tempmappercent.end();++it)
    {
      csv_out << it->second << ","; // writing the value
      percentemptot += it->second; // temporary temp
    }
    // writing and reinitializing the total
    csv_out << percentemptot << endl; // should be 100 if everything is ok
    percentemptot = 0;
  }
  // Done, closing the file stream
  csv_out.close();
  cout << "I am done writing the simple litho file." << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function creates the data structures for keeping track of
//  the channel network but only maps the chi coordinate.
// Mainly used for calculating m/n ratio.
// DOES NOT segment the chi-elevation profiles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator_chi_only(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate)
{
  // These elements access the chi data
  vector< vector<float> > chi_coordinates;
  vector< vector<int> > chi_node_indices;

  // these are for the individual channels
  vector<float> these_chi_coordinates;
  vector<int> these_chi_node_indices;

  // these are maps that will store the data
  map<int,float> chi_coord_map;
  map<int,float> elev_map;
  map<int,float> area_map;
  map<int,float> flow_distance_map;
  vector<int> node_sequence_vec;

  // these are vectors that will store information about the individual nodes
  // that allow us to map the nodes to specific channels during data visualisation

  // These two maps have each node in the channel (the index)
  // linked to a key (either the baselevel key or source key)
  map<int,int> these_source_keys;
  map<int,int> these_baselevel_keys;

  // These two maps link keys, which are incrmented by one, to the
  // junction or node of the baselevel or source
  map<int,int> this_key_to_source_map;
  map<int,int> this_key_to_baselevel_map;

  // these are for working with the FlowInfo object
  int this_node,row,col;
  int this_base_level, this_source_node;

  vector<int> empty_vec;
  ordered_baselevel_nodes = empty_vec;
  ordered_source_nodes = empty_vec;


  // get the number of channels
  int source_node_tracker = -1;
  int baselevel_tracker = -1;
  int ranked_source_node_tracker = -1;
  int n_channels = int(source_nodes.size());
  for(int chan = 0; chan<n_channels; chan++)
  {
    //cout << "Sampling channel " << chan+1 << " of " << n_channels << endl;

    // get the base level
    this_base_level = baselevel_node_of_each_basin[chan];
    //cout << "Got the base level" << endl;

    // If a key to this base level does not exist, add one.
    if ( this_key_to_baselevel_map.find(this_base_level) == this_key_to_baselevel_map.end() )
    {
      baselevel_tracker++;

      // this resets the ranked source node tracker
      ranked_source_node_tracker = -1;

      //cout << "Found a new baselevel. The node is: " << this_base_level << " and key is: " << baselevel_tracker << endl;
      this_key_to_baselevel_map[this_base_level] = baselevel_tracker;
      ordered_baselevel_nodes.push_back(this_base_level);

      // Get the node of the source to the mainstem. This works because the
      // mainstem is always the first channel listed in a basin.
      source_node_of_mainstem_map[baselevel_tracker] = source_nodes[chan];
    }

    // now add the source tracker
    source_node_tracker++;
    ranked_source_node_tracker++;

    // get the source node
    this_source_node = source_nodes[chan];

    // add the node to the trackers so that we can trace individual basin nodes
    // for m over n calculations
    ordered_source_nodes.push_back(this_source_node);
    source_nodes_ranked_by_basin.push_back(ranked_source_node_tracker);

    // now add the source node to the data map
    this_key_to_source_map[this_source_node] = source_node_tracker;

    //cout << "The source key is: " << source_node_tracker << " and basin key is: " << baselevel_tracker << endl;

    // get this particular channel (it is a chi network with only one channel)
    LSDChiNetwork ThisChiChannel(FlowInfo, source_nodes[chan], outlet_nodes[chan],
                                Elevation, FlowDistance, DrainageArea,chi_coordinate);

    // okay the ChiNetwork has all the data about the m vales at this stage.
    // Get these vales and print them to a raster
    chi_coordinates = ThisChiChannel.get_chis();
    chi_node_indices = ThisChiChannel.get_node_indices();

    // now get the number of channels. This should be 1!
    int n_channels = int(chi_coordinates.size());
    if (n_channels != 1)
    {
      cout << "Whoa there, I am trying to make a chi map but something seems to have gone wrong with the channel extraction."  << endl;
      cout << "I should only have one channel per look but I have " << n_channels << " channels." << endl;
    }

    // now get chi coordantes
    these_chi_coordinates = chi_coordinates[0];
    these_chi_node_indices = chi_node_indices[0];

    //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;


    int n_nodes_in_channel = int(these_chi_coordinates.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {

      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;

      // only take the nodes that have not been found
      if (chi_coord_map.find(this_node) == chi_coord_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);

        //cout << "This is a new node; " << this_node << endl;
        chi_coord_map[this_node] = these_chi_coordinates[node];
        elev_map[this_node] = Elevation.get_data_element(row,col);
        area_map[this_node] = DrainageArea.get_data_element(row,col);
        flow_distance_map[this_node] = FlowDistance.get_data_element(row,col);
        node_sequence_vec.push_back(this_node);

        these_source_keys[this_node] = source_node_tracker;
        these_baselevel_keys[this_node] = baselevel_tracker;
      }
      else
      {
        //cout << "I already have node: " << this_node << endl;
      }
    }
  }

  //cout << "I am all finished segmenting the channels!" << endl;

  // set the object data members
  elev_data_map = elev_map;
  chi_data_map = chi_coord_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_sequence_vec;

  source_keys_map = these_source_keys;
  baselevel_keys_map = these_baselevel_keys;
  key_to_source_map = this_key_to_source_map;
  key_to_baselevel_map = this_key_to_baselevel_map;
  //cout << "BUG TRACKER" << endl; exit(EXIT_FAILURE);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating segments from all sources in a DEM
// The sources and their outlets are supplied by the source and outlet nodes
// vectors. These are generated from the LSDJunctionNetwork function
// get_overlapping_channels
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int target_nodes,
                                    int n_iterations, int skip,
                                    int minimum_segment_length, float sigma)
{
  // cout << "Starting here!!!" << endl;
  // cout << "Sources:" << source_nodes.size() << endl;
  // cout << "outlet_nodes:" << outlet_nodes.size() << endl;
  // cout << "baselevel nodes:" << baselevel_node_of_each_basin.size() << endl;

  // for(size_t i=0; i< source_nodes.size(); i++)
  // {
  //   std::cout << "source_nodes:" << source_nodes[i];
  //   std::cout << "outlet_nodes:" << outlet_nodes[i];
  //   std::cout << "baselevel_node_of_each_basin:" << baselevel_node_of_each_basin[i] << std::endl;
  // }

  // IMPORTANT THESE PARAMETERS ARE NOT USED BECAUSE CHI IS CALCULATED SEPARATELY
  // However we need to give something to pass to the Monte carlo functions
  // even through they are not used (they are inherited)
  float A_0 = 1;
  float m_over_n = 0.5;

  // These elements access the chi data
  vector< vector<float> > chi_m_means;
  vector< vector<float> > chi_b_means;
  vector< vector<float> > chi_coordinates;
  vector< vector<int> > chi_node_indices;

  // these are for the individual channels
  vector<float> these_chi_m_means;
  vector<float> these_chi_b_means;
  vector<float> these_chi_coordinates;
  vector<int> these_chi_node_indices;

  // these are maps that will store the data
  map<int,float> m_means_map;
  map<int,float> b_means_map;
  map<int,float> chi_coord_map;
  map<int,float> elev_map;
  map<int,float> area_map;
  map<int,float> flow_distance_map;
  vector<int> node_sequence_vec;

  // these are vectors that will store information about the individual nodes
  // that allow us to map the nodes to specific channels during data visualisation

  // These two maps have each node in the channel (the index)
  // linked to a key (either the baselevel key or source key)
  map<int,int> these_source_keys;
  map<int,int> these_baselevel_keys;

  // These two maps link keys, which are incrmented by one, to the
  // junction or node of the baselevel or source
  map<int,int> this_key_to_source_map;
  map<int,int> this_key_to_baselevel_map;

  // these are for working with the FlowInfo object
  int this_node,row,col;
  int this_base_level, this_source_node;

  // get the number of channels
  int source_node_tracker = -1;
  int baselevel_tracker = -1;
  int ranked_source_node_tracker = -1;
  int n_channels = int(source_nodes.size());
  for(int chan = 0; chan<n_channels; chan++)
  {
    // cout << "Sampling channel " << chan+1 << " of " << n_channels << endl;

    // get the base level
    this_base_level = baselevel_node_of_each_basin[chan];
    // cout << "Got the base level: " << this_base_level << endl;

    // If a key to this base level does not exist, add one.
    if ( this_key_to_baselevel_map.find(this_base_level) == this_key_to_baselevel_map.end() )
    {
      baselevel_tracker++;
      // this resets the ranked source node tracker
      ranked_source_node_tracker = -1;
      // cout << "Found a new baselevel. The node is: " << this_base_level << " and key is: " << baselevel_tracker << endl;
      this_key_to_baselevel_map[this_base_level] = baselevel_tracker;
      ordered_baselevel_nodes.push_back(this_base_level);

      // Get the node of the source to the mainstem. This works because the
      // mainstem is always the first channel listed in a basin.
      source_node_of_mainstem_map[baselevel_tracker] = source_nodes[chan];
    }

    // now add the source tracker
    source_node_tracker++;
    ranked_source_node_tracker++;

    // get the source node
    this_source_node = source_nodes[chan];

    // add the node to the trackers so that we can trace individual basin nodes
    // for m over n calculations
    ordered_source_nodes.push_back(this_source_node);
    source_nodes_ranked_by_basin.push_back(ranked_source_node_tracker);

    // now add the source node to the data map
    this_key_to_source_map[this_source_node] = source_node_tracker;

    // cout << "The source key is: " << source_node_tracker << " and basin key is: " << baselevel_tracker << endl;

    // get this particular channel (it is a chi network with only one channel)
    LSDChiNetwork ThisChiChannel(FlowInfo, source_nodes[chan], outlet_nodes[chan],
                                Elevation, FlowDistance, DrainageArea,chi_coordinate);

    // split the channel
    // cout << "Splitting channels" << endl;
    ThisChiChannel.split_all_channels(A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma);

    // monte carlo sample all channels
    // cout << "Entering the monte carlo sampling" << endl;
    ThisChiChannel.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations, skip, minimum_segment_length, sigma);

    // okay the ChiNetwork has all the data about the m vales at this stage.
    // Get these vales and print them to a raster
    chi_m_means = ThisChiChannel.get_m_means();
    chi_b_means = ThisChiChannel.get_b_means();
    chi_coordinates = ThisChiChannel.get_chis();
    chi_node_indices = ThisChiChannel.get_node_indices();

    // now get the number of channels. This should be 1!
    int n_channels = int(chi_m_means.size());
    if (n_channels != 1)
    {
      cout << "Whoa there, I am trying to make a chi map but something seems to have gone wrong with the channel extraction."  << endl;
      cout << "I should only have one channel per look but I have " << n_channels << " channels." << endl;
    }

    // now get the m_means out
    these_chi_m_means = chi_m_means[0];
    these_chi_b_means = chi_b_means[0];
    these_chi_coordinates = chi_coordinates[0];
    these_chi_node_indices = chi_node_indices[0];

    // cout << "I have " << these_chi_m_means.size() << " nodes." << endl;


    int n_nodes_in_channel = int(these_chi_m_means.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {

      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;

      // only take the nodes that have not been found
      if (m_means_map.find(this_node) == m_means_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);

        //cout << "This is a new node; " << this_node << endl;
        m_means_map[this_node] = these_chi_m_means[node];
        b_means_map[this_node] = these_chi_b_means[node];
        chi_coord_map[this_node] = these_chi_coordinates[node];
        elev_map[this_node] = Elevation.get_data_element(row,col);
        area_map[this_node] = DrainageArea.get_data_element(row,col);
        flow_distance_map[this_node] = FlowDistance.get_data_element(row,col);
        node_sequence_vec.push_back(this_node);

        these_source_keys[this_node] = source_node_tracker;
        these_baselevel_keys[this_node] = baselevel_tracker;

      }
      else
      {
        //cout << "I already have node: " << this_node << endl;
      }
    }
  }

  //cout << "I am all finished segmenting the channels!" << endl;

  // set the object data members
  M_chi_data_map =m_means_map;
  b_chi_data_map = b_means_map;
  elev_data_map = elev_map;
  chi_data_map = chi_coord_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_sequence_vec;

  source_keys_map = these_source_keys;
  baselevel_keys_map = these_baselevel_keys;
  key_to_source_map = this_key_to_source_map;
  key_to_baselevel_map = this_key_to_baselevel_map;

  // get the fitted elevations
  calculate_segmented_elevation(FlowInfo);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// Theoretically this function will only be compiled/used/callable if openMP is activated during compilation
// Let me (Boris) know if there is any conflict with that, but that's unlikely.
// B.G
// THis is a multithreading attempt on getting m_chi, the function is nearly identical to the serial one
// I changed variable initialization to get each thread what they need privately, and (so far) I preconstructed the vectors affected by several thread at a time with their final size.
// If it works, It would mean that our routines can easily be parallel
#ifdef _OPENMP
struct sort_my_chinetworks
{
  // Vector of river node (pointer of)
  LSDChiNetwork* CN;
  // Elevation at eh base of the river
  int chansize;
};

// These are the operator used to sort the river per grater/lower elevation in the priority queue
bool operator>( const sort_my_chinetworks& lhs, const sort_my_chinetworks& rhs )
{
  return lhs.chansize > rhs.chansize;
}
bool operator<( const sort_my_chinetworks& lhs, const sort_my_chinetworks& rhs )
{
  return lhs.chansize < rhs.chansize;
};
//#######################################################################################################
//#######################################################################################################


void LSDChiTools::chi_map_automator(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int target_nodes,
                                    int n_iterations, int skip,
                                    int minimum_segment_length, float sigma, int nthreads)
{

  // IMPORTANT THESE PARAMETERS ARE NOT USED BECAUSE CHI IS CALCULATED SEPARATELY
  // However we need to give something to pass to the Monte carlo functions
  // even through they are not used (they are inherited)
  float A_0 = 1;
  float m_over_n = 0.5;

  // These elements access the chi data
  vector< vector<float> > chi_m_means;
  vector< vector<float> > chi_b_means;
  vector< vector<float> > chi_coordinates;
  vector< vector<int> > chi_node_indices;

  // these are for the individual channels
  vector<float> these_chi_m_means;
  vector<float> these_chi_b_means;
  vector<float> these_chi_coordinates;
  vector<int> these_chi_node_indices;

  // these are maps that will store the data
  map<int,float> m_means_map;
  map<int,float> b_means_map;
  map<int,float> chi_coord_map;
  map<int,float> elev_map;
  map<int,float> area_map;
  map<int,float> flow_distance_map;
  vector<int> node_sequence_vec;

  vector<LSDChiNetwork*> full_chi_network;

  // these are vectors that will store information about the individual nodes
  // that allow us to map the nodes to specific channels during data visualisation

  // These two maps have each node in the channel (the index)
  // linked to a key (either the baselevel key or source key)
  map<int,int> these_source_keys;
  map<int,int> these_baselevel_keys;

  // These two maps link keys, which are incrmented by one, to the
  // junction or node of the baselevel or source
  map<int,int> this_key_to_source_map;
  map<int,int> this_key_to_baselevel_map;

  // these are for working with the FlowInfo object
  int this_node,row,col;
  int this_base_level, this_source_node;

  // get the number of channels
  int source_node_tracker = -1;
  int baselevel_tracker = -1;
  int ranked_source_node_tracker = -1;
  int n_channels = int(source_nodes.size());
  for(int chan = 0; chan<n_channels; chan++)
  {
    //cout << "Sampling channel " << chan+1 << " of " << n_channels << endl;

    // get the base level
    this_base_level = baselevel_node_of_each_basin[chan];
    //cout << "Got the base level" << endl;

    // If a key to this base level does not exist, add one.
    if ( this_key_to_baselevel_map.find(this_base_level) == this_key_to_baselevel_map.end() )
    {
      baselevel_tracker++;
      // this resets the ranked source node tracker
      ranked_source_node_tracker = -1;
      //cout << "Found a new baselevel. The node is: " << this_base_level << " and key is: " << baselevel_tracker << endl;
      this_key_to_baselevel_map[this_base_level] = baselevel_tracker;
      ordered_baselevel_nodes.push_back(this_base_level);

      // Get the node of the source to the mainstem. This works because the
      // mainstem is always the first channel listed in a basin.
      source_node_of_mainstem_map[baselevel_tracker] = source_nodes[chan];
    }

    // now add the source tracker
    source_node_tracker++;
    ranked_source_node_tracker++;

    // get the source node
    this_source_node = source_nodes[chan];

    // add the node to the trackers so that we can trace individual basin nodes
    // for m over n calculations
    ordered_source_nodes.push_back(this_source_node);
    source_nodes_ranked_by_basin.push_back(ranked_source_node_tracker);

    // now add the source node to the data map
    this_key_to_source_map[this_source_node] = source_node_tracker;

    //cout << "The source key is: " << source_node_tracker << " and basin key is: " << baselevel_tracker << endl;

    // get this particular channel (it is a chi network with only one channel)
    // I am storing in the heap and not the stack for testing purposes
    full_chi_network.push_back(new LSDChiNetwork(FlowInfo, source_nodes[chan], outlet_nodes[chan],
                                Elevation, FlowDistance, DrainageArea,chi_coordinate));

    // split the channel
    //cout << "Splitting channels" << endl;
    // ThisChiChannel.split_all_channels(A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma);
    // full_chi_network.push_back(ThisChiChannel);
    // monte carlo sample all channels
    //cout << "Entering the monte carlo sampling" << endl;
    // ThisChiChannel.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations, skip, minimum_segment_length, sigma);

    // okay the ChiNetwork has all the data about the m vales at this stage.
    // Get these vales and print them to a raster
    // chi_m_means = ThisChiChannel.get_m_means();
    // chi_b_means = ThisChiChannel.get_b_means();
    // chi_coordinates = ThisChiChannel.get_chis();
    chi_node_indices = full_chi_network.back()->get_node_indices();

    // now get the number of channels. This should be 1!
    int n_channels = int(chi_node_indices.size());
    if (n_channels != 1)
    {
      cout << "Whoa there, I am trying to make a chi map but something seems to have gone wrong with the channel extraction."  << endl;
      cout << "I should only have one channel per look but I have " << n_channels << " channels." << endl;
    }

    // now get the m_means out
    // these_chi_m_means = chi_m_means[0];
    // these_chi_b_means = chi_b_means[0];
    // these_chi_coordinates = chi_coordinates[0];
    these_chi_node_indices = chi_node_indices[0];

    //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;


    int n_nodes_in_channel = int(these_chi_node_indices.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {

      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;

      // only take the nodes that have not been found
      if (m_means_map.find(this_node) == m_means_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);

        //cout << "This is a new node; " << this_node << endl;
        // m_means_map[this_node] = these_chi_m_means[node];
        // b_means_map[this_node] = these_chi_b_means[node];
        // chi_coord_map[this_node] = these_chi_coordinates[node];
        elev_map[this_node] = Elevation.get_data_element(row,col);
        area_map[this_node] = DrainageArea.get_data_element(row,col);
        flow_distance_map[this_node] = FlowDistance.get_data_element(row,col);
        node_sequence_vec.push_back(this_node);

        these_source_keys[this_node] = source_node_tracker;
        these_baselevel_keys[this_node] = baselevel_tracker;

      }
      else
      {
        //cout << "I already have node: " << this_node << endl;
      }
    }
  }

  // Last step: sorting my rivers per size: I want my largest one to be dynamically allocated to the different thread first!
  priority_queue< sort_my_chinetworks, vector<sort_my_chinetworks>, greater<sort_my_chinetworks> > myp;
  for (int i=0; i<int(full_chi_network.size()); i++)
  {
    // GEtting each chi_network
    LSDChiNetwork* this_CN = full_chi_network[i];
    vector<int> these = this_CN->get_node_indices()[0];
    int this_chansize = these.size();
    sort_my_chinetworks this_sort;
    this_sort.CN = this_CN;
    this_sort.chansize = this_chansize;
    myp.push(this_sort);
  }

  // I have all the channel oredered, I need to reconstruct the vector
  // full_chi_network.clear();
  vector<LSDChiNetwork*> full_chi_network2;
  while(myp.size()>0)
  {
    sort_my_chinetworks this_sort = myp.top(); // Getting the first element
    LSDChiNetwork* this_CN = this_sort.CN;
    myp.pop(); // Getting rid of the top element
    if(this_sort.chansize>30)
      full_chi_network2.push_back(this_CN);
    // cout << "Size: " << this_sort.chansize << endl;
  }
  full_chi_network = full_chi_network2;
  // I need to reverse it lol, otherwise my largest river are processed last 
  reverse(full_chi_network.begin(),full_chi_network.end());
  // DOne :)

  cout<< "I will now proceed to the Monte Carlo iterative scheme: I have to process " << full_chi_network.size() << "rivers!" << endl;
  cout << "I will test " << n_iterations << " combination of nodes for EACH of these rivers. This will take time!" << endl;
  cout << "I will let you know when I'll be done so you can contact my masters if I crash." << endl;

  //------------ KEEP THAT PART -------------------------------------
  //------------ OpenMP test, crashes on windows for some reasons ---
  //------------ Testing manual threasing below in case -------------
  // #pragma omp parallel num_threads(nthreads)
  // {
    
  //   #pragma omp for schedule(dynamic,2)
  //   for (int i=0; i<int(full_chi_network.size()); i++)
  //   {
  //     // ThisChiChannel.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations, skip, minimum_segment_length, sigma);
  //     float tA_0 = A_0; float  tm_over_n = m_over_n; int tn_iterations = n_iterations; int tskip = skip;int ttarget_nodes = target_nodes;  int tminimum_segment_length = minimum_segment_length; float tsigma = sigma;
  //     LSDChiNetwork* DatChan = full_chi_network[i];
  //     cout<< "--[" << i << "]--" << endl;
  //     DatChan->split_all_channels(tA_0, tm_over_n, tn_iterations, tskip, ttarget_nodes, tminimum_segment_length, tsigma);
  //     DatChan->monte_carlo_sample_river_network_for_best_fit_after_breaks(tA_0, tm_over_n, tn_iterations, tskip, tminimum_segment_length, tsigma);
  //   }
  //   #pragma omp barrier
  // }

  //--------- End of onpenMP
  //------------------------


  //------ c++11 thread test
  unsigned int n_possible_threads = thread::hardware_concurrency();
  if(int(n_possible_threads)>nthreads)
    n_possible_threads = nthreads; // Reducing your number of threads
  else if(n_possible_threads == 0)
  {
    cout << "Erm, I did not manage to understand how many concurrency thread you can have. I am defaulting to 2." << endl;
    n_possible_threads = 2; // Some program cannot get the number of threads from the computer, Therefore I am defaulting to 2
  }
  // Now creating the different vector of pointers for all the threads
  vector<vector<LSDChiNetwork*> > multithreading_tasks(n_possible_threads);
  for(size_t tec=0;tec<n_possible_threads;tec++)
  {
    vector<LSDChiNetwork*> empty_vec;
    multithreading_tasks[tec] = empty_vec;
  }

  // giving one river at each thead containers
  int comptathread = 0;
  for(size_t tec=0; tec<full_chi_network.size();tec++)
  {
    multithreading_tasks[comptathread].push_back(full_chi_network[tec]);
    // Preparing the next round
    if(comptathread<n_possible_threads - 1)
      comptathread++;
    else
      comptathread=0;
  }

  // OK let's multithread, note that the vector of threads is minus 1 to still use the main one
  vector<std::thread> running_threads;
  for(size_t tec=1; tec<n_possible_threads;tec++) // looping though my tasks to multithread and leaving the first one for my main
  {
    running_threads.push_back(std::thread (&LSDChiTools::internal_function_multi_ksn, this, std::ref(multithreading_tasks[tec]),A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma));
  }
  // ALso running for my main thread
  this->internal_function_multi_ksn(multithreading_tasks[0] ,A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma);
  // Once done with main thread, waiting for others
  for(size_t tec=0; tec< running_threads.size(); tec++)
    running_threads[tec].join();

  // cout.clear(); // reenabling cout
  cout << "Monte Carlo samplig over!" << endl;

  // Getting the data out of this shit
  for(size_t i=0; i<full_chi_network.size(); i++)
  { 

    chi_m_means = full_chi_network[i]->get_m_means();
    chi_b_means = full_chi_network[i]->get_b_means();
    chi_coordinates = full_chi_network[i]->get_chis(); 
    chi_node_indices = full_chi_network[i]->get_node_indices();

    if(chi_m_means.size()==0)
      continue;

    these_chi_m_means = chi_m_means[0];
    these_chi_b_means = chi_b_means[0];
    these_chi_coordinates = chi_coordinates[0];
    these_chi_node_indices = chi_node_indices[0];
    //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;

    int this_node,row,col;
    int n_nodes_in_channel = int(these_chi_m_means.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {


      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;

      // only take the nodes that have not been found
      if (m_means_map.find(this_node) == m_means_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);

        //cout << "This is a new node; " << this_node << endl;
        m_means_map[this_node] = these_chi_m_means[node];
        b_means_map[this_node] = these_chi_b_means[node];
        chi_coord_map[this_node] = these_chi_coordinates[node];

      }
      else
      {
        //cout << "I already have node: " << this_node << endl;
      }
    }
  }


  // cout << "Done with stuff here" << endl;

  // set the object data members
  M_chi_data_map = m_means_map;
  b_chi_data_map = b_means_map;
  elev_data_map = elev_map;
  chi_data_map = chi_coord_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_sequence_vec;

  source_keys_map = these_source_keys;
  baselevel_keys_map = these_baselevel_keys;
  key_to_source_map = this_key_to_source_map;
  key_to_baselevel_map = this_key_to_baselevel_map;

  // The following loop removes my chi networks from memory to avoid bad leaks
  for(size_t tec=0; tec<full_chi_network.size(); tec++)
    delete full_chi_network[tec];

  // delete full_chi_network;
  // delete full_chi_network2;
  // for(size_t tec=0; tec<full_chi_network2.size(); tec++)
  //   delete full_chi_network2[tec];

  // get the fitted elevations
  calculate_segmented_elevation(FlowInfo);

}

void LSDChiTools::internal_function_multi_ksn(vector<LSDChiNetwork*>& me_vec_of_chi_network, float tA_0 , float  tm_over_n , int tn_iterations, int tskip , int ttarget_nodes , int tminimum_segment_length , float tsigma )
{
  for (size_t tec=0; tec<me_vec_of_chi_network.size(); tec++)
  {
    LSDChiNetwork* DatChan = me_vec_of_chi_network[tec];
    DatChan->split_all_channels(tA_0, tm_over_n, tn_iterations, tskip, ttarget_nodes, tminimum_segment_length, tsigma);
    // cout << "||" <<  std::this_thread::get_id() << "||" << tec << "/" <<  me_vec_of_chi_network.size()  << "||" << endl;
    DatChan->monte_carlo_sample_river_network_for_best_fit_after_breaks(tA_0, tm_over_n, tn_iterations, tskip, tminimum_segment_length, tsigma);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// Attempt 1: failed. I am keeping it for one commit just in case
// void LSDChiTools::chi_map_automator(LSDFlowInfo& FlowInfo,
//                                     vector<int> source_nodes,
//                                     vector<int> outlet_nodes,
//                                     vector<int> baselevel_node_of_each_basin,
//                                     LSDRaster& Elevation, LSDRaster& FlowDistance,
//                                     LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
//                                     int target_nodes,
//                                     int n_iterations, int skip,
//                                     int minimum_segment_length, float sigma, int nthreads)
// {

//   if(nthreads>omp_get_max_threads())
//   {
//     cout << "FATAL_ERROR::You are trying to run m_chi routines in parallel with more processor than I own. I'll go on strike if you do that." << endl;
//     cout << "Alright let's be nice, I'll use as much power as I can then: " << omp_get_max_threads() << " threads (part of cpu)." << endl;
//     nthreads = omp_get_max_threads();
//     // exit(EXIT_FAILURE);
//   }

//   // IMPORTANT THESE PARAMETERS ARE NOT USED BECAUSE CHI IS CALCULATED SEPARATELY
//   // However we need to give something to pass to the Monte carlo functions
//   // even through they are not used (they are inherited)
//   float A_0 = 1;
//   float m_over_n = 0.5;

//   // NUmber of channels to process
//   int n_channels = int(source_nodes.size());




//   // these are maps that will store the data
//   map<int,float> m_means_map;
//   map<int,float> b_means_map;
//   map<int,float> chi_coord_map;
//   map<int,float> elev_map;
//   map<int,float> area_map;
//   map<int,float> flow_distance_map;
//   vector<vector<int> > these_node_sequence_vec(n_channels);

//   // these are vectors that will store information about the individual nodes
//   // that allow us to map the nodes to specific channels during data visualisation

//   // These two maps have each node in the channel (the index)
//   // linked to a key (either the baselevel key or source key)
//   map<int,int> these_source_keys;
//   map<int,int> these_baselevel_keys;

//   // These two maps link keys, which are incrmented by one, to the
//   // junction or node of the baselevel or source
//   map<int,int> this_key_to_source_map;
//   map<int,int> this_key_to_baselevel_map;

//   vector<vector<int> > temp_ordered_baselevel_nodes(n_channels);
//   vector<vector<int> > temp_ordered_source_nodes(n_channels);
//   vector<vector<int> > temp_source_nodes_ranked_by_basin(n_channels);



//   // get the number of channels
//   int source_node_tracker = -1;
//   int baselevel_tracker = -1;
//   int ranked_source_node_tracker = -1;
//   #pragma omp parallel num_threads(nthreads)
//   {
//     // It means that I want my threads to process one river each at a time: because node are disorganise there is no point doing it by chunks
//     #pragma omp for schedule (dynamic,1) 
//     for(int chan = 0; chan<n_channels; chan++)
//     {
//       int thisT = omp_get_thread_num();
//       cout<<"THREAD NUMBER " << thisT << endl;
//       // these are for the individual channels
//       // these are for working with the FlowInfo object
//       int this_node,row,col;
//       int this_base_level, this_source_node;
//         // These elements access the chi data -> one vector per channel
//       vector< vector<float> > chi_m_means;
//       vector< vector<float> > chi_b_means;
//       vector< vector<float> > chi_coordinates;
//       vector< vector<int> > chi_node_indices;
//       vector<float> these_chi_m_means;
//       vector<float> these_chi_b_means;
//       vector<float> these_chi_coordinates;
//       vector<int> these_chi_node_indices;
//       vector<int> this_temp_ordered_baselevel_nodes;
//       vector<int> this_temp_ordered_source_nodes;
//       vector<int> this_temp_source_nodes_ranked_by_basin;

//       // get the base level
//       this_base_level = baselevel_node_of_each_basin[chan];

//       // If a key to this base level does not exist, add one.
//       if ( this_key_to_baselevel_map.find(this_base_level) == this_key_to_baselevel_map.end() )
//       {
//         // cout << "BL tracking DEBUG 1 || T" << thisT << endl;
//         baselevel_tracker++;
//         // this resets the ranked source node tracker
//         ranked_source_node_tracker = -1;
//         //cout << "Found a new baselevel. The node is: " << this_base_level << " and key is: " << baselevel_tracker << endl;
//         this_key_to_baselevel_map[this_base_level] = baselevel_tracker;
//         this_temp_ordered_baselevel_nodes.push_back(this_base_level);

//         // Get the node of the source to the mainstem. This works because the
//         // mainstem is always the first channel listed in a basin.
//         // source_node_of_mainstem_map[baselevel_tracker] = source_nodes[chan];
//         // cout << "BL tracking DEBUG 2  || T" << thisT << endl;

//       }

//       // now add the source tracker
//       source_node_tracker++;
//       ranked_source_node_tracker++;

//       // get the source node
//       this_source_node = source_nodes[chan];

//       // add the node to the trackers so that we can trace individual basin nodes
//       // for m over n calculations
//       this_temp_ordered_source_nodes.push_back(this_source_node);
//       this_temp_source_nodes_ranked_by_basin.push_back(ranked_source_node_tracker);

//       // now add the source node to the data map
//       this_key_to_source_map[this_source_node] = source_node_tracker;

//       //cout << "The source key is: " << source_node_tracker << " and basin key is: " << baselevel_tracker << endl;

//       // get this particular channel (it is a chi network with only one channel)
//       cout.setstate(ios_base::failbit); // disabling cout messages
//       LSDChiNetwork ThisChiChannel(FlowInfo, source_nodes[chan], outlet_nodes[chan],
//                                   Elevation, FlowDistance, DrainageArea,chi_coordinate);

//       // split the channel
//       //cout << "Splitting channels" << endl;
//       ThisChiChannel.split_all_channels(A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma);

//       // monte carlo sample all channels
//       //cout << "Entering the monte carlo sampling" << endl;
//       ThisChiChannel.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations, skip, minimum_segment_length, sigma);
//       cout.clear(); // reenabling cout
//       // okay the ChiNetwork has all the data about the m vales at this stage.
//       // Get these vales and print them to a raster
//       chi_m_means = ThisChiChannel.get_m_means();
//       chi_b_means = ThisChiChannel.get_b_means();
//       chi_coordinates = ThisChiChannel.get_chis();
//       chi_node_indices = ThisChiChannel.get_node_indices();

//       // now get the number of channels. This should be 1!
//       int n_channels = int(chi_m_means.size());
//       if (n_channels != 1)
//       {
//         cout << "Whoa there, I am trying to make a chi map but something seems to have gone wrong with the channel extraction."  << endl;
//         cout << "I should only have one channel per look but I have " << n_channels << " channels." << endl;
//       }

//       // now get the m_means out
//       these_chi_m_means = chi_m_means[0];
//       these_chi_b_means = chi_b_means[0];
//       these_chi_coordinates = chi_coordinates[0];
//       these_chi_node_indices = chi_node_indices[0];

//       //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;

//       temp_ordered_baselevel_nodes[chan] = this_temp_ordered_baselevel_nodes;
//       temp_ordered_source_nodes[chan] = this_temp_ordered_source_nodes;
//       temp_source_nodes_ranked_by_basin[chan] = this_temp_source_nodes_ranked_by_basin;


//       int n_nodes_in_channel = int(these_chi_m_means.size());
//       vector<int> temp_node_river(n_nodes_in_channel);
//       for (int node = 0; node< n_nodes_in_channel; node++)
//       {

//         this_node =  these_chi_node_indices[node];
//         //cout << "This node is " << this_node << endl;

//         // only take the nodes that have not been found
//         if (m_means_map.find(this_node) == m_means_map.end() )
//         {
//           FlowInfo.retrieve_current_row_and_col(this_node,row,col);

//           //cout << "This is a new node; " << this_node << endl;
//           m_means_map[this_node] = these_chi_m_means[node];
//           b_means_map[this_node] = these_chi_b_means[node];
//           chi_coord_map[this_node] = these_chi_coordinates[node];
//           elev_map[this_node] = Elevation.get_data_element(row,col);
//           area_map[this_node] = DrainageArea.get_data_element(row,col);
//           flow_distance_map[this_node] = FlowDistance.get_data_element(row,col);
//           temp_node_river.push_back(this_node);

//           these_source_keys[this_node] = source_node_tracker;
//           these_baselevel_keys[this_node] = baselevel_tracker;

//         }
//         else
//         {
//           //cout << "I already have node: " << this_node << endl;
//         }
//       }
//       these_node_sequence_vec[chan] = temp_node_river;
//     }
//   }

//   // flattening the node sequence
//   vector<int> node_sequence_vec;
//   for(size_t i=0; i< these_node_sequence_vec.size(); i++)
//   {
//     vector<int> this_vecnode = these_node_sequence_vec[i];
//     for (size_t j=0; j<this_vecnode.size();j++)
//       node_sequence_vec.push_back(this_vecnode[j]);
//   }

//   // Flattening other stuff
//   for(size_t i=0; i<temp_ordered_baselevel_nodes.size();i++)
//   {
//     vector<int> this_stuff = temp_ordered_baselevel_nodes[i], this_stuff2 =  temp_source_nodes_ranked_by_basin[i], this_stuff3 = temp_ordered_source_nodes[i];

//     for(size_t j = 0; j<this_stuff.size(); j++)
//     {
//       ordered_baselevel_nodes.push_back(this_stuff[j]);
//     }
//     for(size_t j = 0; j<this_stuff2.size(); j++)
//     {
//       source_nodes_ranked_by_basin.push_back(this_stuff2[j]);
//     }
//     for(size_t j = 0; j<this_stuff3.size(); j++)
//     {
//       ordered_source_nodes.push_back(this_stuff3[j]);
//     }
//   }

//   cout << "I am all finished segmenting the channels!" << endl;

//   // set the object data members
//   M_chi_data_map =m_means_map;
//   b_chi_data_map = b_means_map;
//   elev_data_map = elev_map;
//   chi_data_map = chi_coord_map;
//   flow_distance_data_map = flow_distance_map;
//   drainage_area_data_map = area_map;
//   node_sequence = node_sequence_vec;

//   source_keys_map = these_source_keys;
//   baselevel_keys_map = these_baselevel_keys;
//   key_to_source_map = this_key_to_source_map;
//   key_to_baselevel_map = this_key_to_baselevel_map;

//   // get the fitted elevations
//   calculate_segmented_elevation(FlowInfo);

// }
// //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#endif 

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Hacky version of the chi_automator when applying a filtering to the topography
// Reason is that the entire knicpoint algorithm is built on map preprocessed by this function, but I need them before doing the segmentation
// I'll make that more elegant when time will allow it
// BG
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_pre_automator(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int target_nodes,
                                    int n_iterations, int skip,
                                    int minimum_segment_length, float sigma)
{
  vector< vector<int> > chi_node_indices;
  vector<int> these_chi_node_indices;
  vector<int> node_sequence_vec;
  map<int,int> these_source_keys;
  map<int,int> these_baselevel_keys;
  map<int,int> this_key_to_source_map;
  map<int,int> this_key_to_baselevel_map;

  // these are for working with the FlowInfo object
  int this_node,row,col;
  int this_base_level, this_source_node;

  // get the number of channels
  int source_node_tracker = -1;
  int baselevel_tracker = -1;
  int ranked_source_node_tracker = -1;
  int n_channels = int(source_nodes.size());
  for(int chan = 0; chan<n_channels; chan++)
  {
    //cout << "Sampling channel " << chan+1 << " of " << n_channels << endl;

    // get the base level
    this_base_level = baselevel_node_of_each_basin[chan];
    //cout << "Got the base level" << endl;

    // If a key to this base level does not exist, add one.
    if ( this_key_to_baselevel_map.find(this_base_level) == this_key_to_baselevel_map.end() )
    {
      baselevel_tracker++;
      // this resets the ranked source node tracker
      ranked_source_node_tracker = -1;
      //cout << "Found a new baselevel. The node is: " << this_base_level << " and key is: " << baselevel_tracker << endl;
      this_key_to_baselevel_map[this_base_level] = baselevel_tracker;
      ordered_baselevel_nodes.push_back(this_base_level);

      // Get the node of the source to the mainstem. This works because the
      // mainstem is always the first channel listed in a basin.
      source_node_of_mainstem_map[baselevel_tracker] = source_nodes[chan];
    }

    // now add the source tracker
    source_node_tracker++;
    ranked_source_node_tracker++;

    // get the source node
    this_source_node = source_nodes[chan];

    // add the node to the trackers so that we can trace individual basin nodes
    // for m over n calculations
    ordered_source_nodes.push_back(this_source_node);
    source_nodes_ranked_by_basin.push_back(ranked_source_node_tracker);

    // now add the source node to the data map
    this_key_to_source_map[this_source_node] = source_node_tracker;

    //cout << "The source key is: " << source_node_tracker << " and basin key is: " << baselevel_tracker << endl;

    // get this particular channel (it is a chi network with only one channel)
    LSDChiNetwork ThisChiChannel(FlowInfo, source_nodes[chan], outlet_nodes[chan],
                                Elevation, FlowDistance, DrainageArea,chi_coordinate);
    chi_node_indices = ThisChiChannel.get_node_indices();
    these_chi_node_indices = chi_node_indices[0];

    //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;


    int n_nodes_in_channel = int(these_chi_node_indices.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {

      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;

      FlowInfo.retrieve_current_row_and_col(this_node,row,col);

      node_sequence_vec.push_back(this_node);

      these_source_keys[this_node] = source_node_tracker;
      these_baselevel_keys[this_node] = baselevel_tracker;

    }
  }

 
  node_sequence = node_sequence_vec;
  source_keys_map = these_source_keys;
  baselevel_keys_map = these_baselevel_keys;
  key_to_source_map = this_key_to_source_map;
  key_to_baselevel_map = this_key_to_baselevel_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is a much more rudimentary version that mimics the
// channel steepness caluclations.
// chi needs to be calculated outside of the function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator_rudimentary(LSDFlowInfo& FlowInfo,
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    vector<int> baselevel_node_of_each_basin,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance,
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate,
                                    int regression_nodes)
{

  // the data is stored in maps, for easier testing if a node has been
  // visited.
  // You might consider having these as data elements in the object so you don't
  // have to pass them
  map<int,float> gradient_data_map;
  map<int,float> intercept_data_map;
  map<int,float> R2_data_map;
  map<int,float> chi_coordinate_data_map;
  map<int,float> elevation_data_map;
  map<int,float> flow_distance_map;
  map<int,float> area_map;
  vector<int> node_order;

  // check if the number of nodes are odd .If not add 1
  if (regression_nodes % 2 == 0)
  {
    cout << "Hello user. You need an odd number of regression nodes." << endl;
    regression_nodes = regression_nodes+1;
    cout << " Changing your regression nodes to " << regression_nodes << endl;
  }

  // now get the midpoint
  int mp_nodes = (regression_nodes-1)/2;

  //cout << "The number of mp nodes is: " << mp_nodes << endl;

  // these keep track of the beginning and ending nodes of a given channel
  int channel_start_node;
  int channel_end_node;
  float channel_end_elevation;

  // vectors for holding the chi elevation data
  vector<float> chi_vec;
  vector<float> elev_vec;
  vector<float> empty_vec;

  // these are extracted from the channel segments using linear regression
  float intercept,gradient,R_squared;

  //float this_chi;
  //float this_elev;
  int this_mp_node;
  int this_end_node;
  int this_start_node;

  // these are for getting information out of the FlowInfo object
  int row,col, this_node;
  int r_node, r_row,r_col;          // reciever row and column.

  // The way this works is that it starts at the top of a channel. It then works
  // its way down and find the node that is the midpoint and the node that is the
  // end point. The midpoint node is where the data will be recorded.
  // It then puts the data from the start node to the end node into a vector
  // and performs a linear regression of this vector. The regression data from these
  // vectors are recorded at the nodes.
  // We then want to cover all the nodes with data so what happens if some nodes
  // do not become midpoints?
  // We could start at the top and get the first midpoint.
  // From there we can work our way down checking if the top of the regression segment
  // is more than one node down from the end point...

  // get the number of channels
  int n_channels = int(source_nodes.size());
  // now loop through the channels
  for(int chan = 0; chan<n_channels; chan++)
  {
    channel_start_node = source_nodes[chan];
    channel_end_node = outlet_nodes[chan];

    // Get the elevation of the end node as a secondary check of the ending of the channel
    // segment
    FlowInfo.retrieve_current_row_and_col(channel_end_node,row,col);
    channel_end_elevation = Elevation.get_data_element(row,col);

    // reset the flag for ending the channel
    bool is_end_of_channel = false;

    // set the segment start node to the channel start node
    this_start_node = channel_start_node;

    // now retrieve the midpoint node
    this_node = channel_start_node;
    for(int n = 0; n<mp_nodes; n++)
    {
      FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
      this_node = r_node;
    }
    this_mp_node = this_node;
    this_node = r_node;

    // now go down one step
    FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
    this_node = r_node;

    // now get the end node
    for(int n = 0; n<mp_nodes; n++)
    {
      FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
      this_node = r_node;
    }
    this_end_node = this_node;

    //================================================
    // This loop is for bug checking
    //this_node = this_start_node;
    //do
    //{
    //  // get the elevation and chi vectors by following the flow
    //  cout << "This node is: " << this_node << endl;
    //  FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    //  FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
    //  this_node = r_node;
    //}
    //while(this_node != this_end_node);
    //
    //cout << "And the midpoint node was: " << this_mp_node << endl;
    //================================================

    // we search down the channel, collecting linear regressions at the
    // midpoint of the intervals
    while (is_end_of_channel == false)
    {
      // get a vector of chi and elevation from the start node to the end node
      chi_vec = empty_vec;
      elev_vec = empty_vec;

      // copy the data elements into the vecotrs. This is a little stupid
      // because one might just use a deque to pop the first element
      // and push the last, but the linear regression takes vectors,
      // not deques so you would have to copy the deque element-wise anyway
      // If you wanted, you could speed this up by implementing a linear regression
      // of deques, but that will need to wait for another day.
      this_node = this_start_node;
      do
      {
        // get the elevation and chi vectors by following the flow
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        chi_vec.push_back(chi_coordinate.get_data_element(row,col));
        elev_vec.push_back(Elevation.get_data_element(row,col));

        FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
        this_node = r_node;
      } while(this_node != this_end_node);

      // do a linear regression on the segment
      least_squares_linear_regression(chi_vec,elev_vec, intercept, gradient, R_squared);

      // now add the intercept and gradient data to the correct node
      // only take data that has not been calculated before
      // The channels are in order of descending length so data from
      // longer channels take precidence.
      if (gradient_data_map.find(this_mp_node) == gradient_data_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_mp_node,row,col);
        gradient_data_map[this_mp_node] = gradient;
        intercept_data_map[this_mp_node] = intercept;
        R2_data_map[this_mp_node] = R_squared;
        chi_coordinate_data_map[this_mp_node] = chi_coordinate.get_data_element(row,col);
        elevation_data_map[this_mp_node] = Elevation.get_data_element(row,col);
        flow_distance_map[this_mp_node] = FlowDistance.get_data_element(row,col);
        area_map[this_mp_node] = DrainageArea.get_data_element(row,col);
        node_order.push_back(this_mp_node);
      }
      else
      {
        is_end_of_channel = true;
      }

      // now move all the nodes down one
      FlowInfo.retrieve_receiver_information(this_start_node,r_node,r_row,r_col);
      this_start_node = r_node;

      FlowInfo.retrieve_receiver_information(this_mp_node,r_node,r_row,r_col);
      this_mp_node = r_node;

      FlowInfo.retrieve_receiver_information(this_end_node,r_node,r_row,r_col);
      this_end_node = r_node;

      // check if we are at the end of the channel
      if (this_end_node == channel_end_node)
      {
        is_end_of_channel = true;
      }
      // also check if the end node is lower elevation than the end node,
      // just to try and stop the channel passing the end node
      FlowInfo.retrieve_current_row_and_col(this_end_node,row,col);
      if (channel_end_elevation > Elevation.get_data_element(row,col))
      {
        is_end_of_channel = true;
      }
    }          // This finishes the regression segment loop
  }            // This finishes the channel and resets channel start and end nodes

  // set the data objects
  M_chi_data_map = gradient_data_map;
  b_chi_data_map = intercept_data_map;
  elev_data_map = elevation_data_map;
  chi_data_map = chi_coordinate_data_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_order;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the basins and an additional file that has basin centroids
// and labelling information for plotting
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDChiTools::get_basin_raster(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Junctions)
{
  int N_Juncs = Junctions.size();
  LSDCoordinateConverterLLandUTM Converter;


  // Get some data members for holding basins and the raster
  vector<LSDBasin> AllTheBasins;
  map<int,int> drainage_of_other_basins;
  LSDIndexRaster BasinMasterRaster;

  //cout << "I am trying to print basins, found " << N_BaseLevelJuncs << " base levels." << endl;
  // Loop through the junctions
  for(int BN = 0; BN<N_Juncs; BN++)
  {
    //cout << "Getting basin " << BN << " and the junction is: "  << BaseLevelJunctions[BN] << endl;
    LSDBasin thisBasin(Junctions[BN],FlowInfo, JunctionNetwork);
    //cout << "...got it!" << endl;
    AllTheBasins.push_back(thisBasin);

    // This is required if the basins are nested--test the code which numbers
    // to be overwritten by a smaller basin
    drainage_of_other_basins[Junctions[BN]] = thisBasin.get_NumberOfCells();

  }

  // now loop through everything again getting the raster
  if (N_Juncs > 0)     // this gets the first raster
  {
    BasinMasterRaster = AllTheBasins[0].write_integer_data_to_LSDIndexRaster(Junctions[0], FlowInfo);
  }

  // now add on the subsequent basins
  for(int BN = 1; BN<N_Juncs; BN++)
  {
    AllTheBasins[BN].add_basin_to_LSDIndexRaster(BasinMasterRaster, FlowInfo,
                              drainage_of_other_basins, Junctions[BN]);
  }

  return BasinMasterRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// extract lithology data per basins from a litho map, a FlowInfo and a junction networks
// ongoing work
// BG - 15/09/2017
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,map<int,int> > LSDChiTools::get_basin_lithocount(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork, LSDIndexRaster& litho,
                               vector<int> Junctions)
{
  int N_Juncs = Junctions.size();
  LSDCoordinateConverterLLandUTM Converter;
  map<int,map<int,int> > lithocount;

  // Get some data members for holding basins and the raster
  vector<LSDBasin> AllTheBasins;
  map<int,int> drainage_of_other_basins;
  LSDIndexRaster BasinMasterRaster;

  //cout << "I am trying to print basins, found " << N_BaseLevelJuncs << " base levels." << endl;
  // Loop through the junctions
  for(int BN = 0; BN<N_Juncs; BN++)
  {
    //cout << "Getting basin " << BN << " and the junction is: "  << BaseLevelJunctions[BN] << endl;
    LSDBasin thisBasin(Junctions[BN],FlowInfo, JunctionNetwork);
    //cout << "...got it!" << endl;
    AllTheBasins.push_back(thisBasin);

    // This is required if the basins are nested--test the code which numbers
    // to be overwritten by a smaller basin
    drainage_of_other_basins[Junctions[BN]] = thisBasin.get_NumberOfCells();
    lithocount[BN] = thisBasin.count_unique_values_from_litho_raster(litho, FlowInfo);
  }

  return lithocount;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to tag channels with a segment number
// It decides on segments if the M_Chi value has changed so should only be used
// with chi networks that have used a skip of 0 and a monte carlo itertions of 1
// This data is used by other routines to look at the spatial distribution of
// hillslope-channel coupling.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::segment_counter(LSDFlowInfo& FlowInfo, float maximum_segment_length)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node;
  int segment_counter = 0;
  map<int,int> this_segment_counter_map;
  float last_M_chi, this_M_chi, last_flow_length, this_flow_length, segment_length;

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    this_node = node_sequence[0];
    last_M_chi =  M_chi_data_map[this_node];
    last_flow_length = flow_distance_data_map[this_node];

    for (int n = 0; n< n_nodes; n++)
    {

      // Get the M_chi and flow_length from the current node
      this_node = node_sequence[n];
      this_M_chi = M_chi_data_map[this_node];
      this_flow_length = flow_distance_data_map[this_node];

      // update the current segment length
      segment_length = fabs(this_flow_length-last_flow_length);

      // If the M_chi has changed, increment the segment counter
      if (this_M_chi != last_M_chi)
      {
        segment_counter++;
        last_M_chi = this_M_chi;
        last_flow_length = this_flow_length;
      }
      // Else if the segment length has exceeded the maximum length of segments increment the segment counter
      else if (segment_length >= maximum_segment_length)
      {
        segment_counter++;
        last_M_chi = this_M_chi;
        last_flow_length = this_flow_length;
      }

      // Print the segment counter to the data map
      this_segment_counter_map[this_node]  = segment_counter;
    }
  }
  segment_counter_map = this_segment_counter_map;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to tag channels with a segment number
// It decides on segments if the M_Chi value has changed so should only be used
// with chi networks that have used a skip of 0 and a monte carlo itertions of 1
// This data is used by other routines to look at the spatial distribution of
// hillslope-channel coupling. Generates an LSDIndexRaster of the channel network
// indexed by segment numbers for feeding to the hilltop flow routing analysis.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDChiTools::segment_mapping(LSDFlowInfo& FlowInfo, float maximum_segment_length)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node, row, col;
  int segment_counter = 0;
  map<int,int> this_segment_counter_map;
  float last_M_chi, this_M_chi, last_flow_length, this_flow_length, segment_length;

  //declare empty array for raster generation
  Array2D<int> SegmentedStreamNetworkArray(NRows,NCols,-9999);

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    //get the node
    this_node = node_sequence[0];
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);

    last_M_chi =  M_chi_data_map[this_node];
    last_flow_length = flow_distance_data_map[this_node];

    for (int n = 0; n< n_nodes; n++)
    {

      // Get the M_chi from the current node
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      this_M_chi = M_chi_data_map[this_node];
      this_flow_length = flow_distance_data_map[this_node];

      // update the current segment length
      segment_length = fabs(this_flow_length-last_flow_length);

      // If the M_chi has changed, increment the segment counter
      if (this_M_chi != last_M_chi)
      {
        segment_counter++;
        last_M_chi = this_M_chi;
        last_flow_length = this_flow_length;
      }
      // Else if the segment length has exceeded the maximum length of segments increment the segment counter
      else if (segment_length >= maximum_segment_length)
      {
        segment_counter++;
        last_M_chi = this_M_chi;
        last_flow_length = this_flow_length;
      }

      // Print the segment counter to the data map and raster
      this_segment_counter_map[this_node]  = segment_counter;
      SegmentedStreamNetworkArray[row][col] = segment_counter;
    }
  }
  segment_counter_map = this_segment_counter_map;

  return LSDIndexRaster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,SegmentedStreamNetworkArray,GeoReferencingStrings);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to tag channels with a segment number
// It decides on segments if the M_Chi value has changed so should only be used
// with chi networks that have used a skip of 0 and a monte carlo itertions of 1
// This data is used by other routines to look at the spatial distribution of
// hillslope-channel coupling.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::segment_counter_knickpoint(LSDFlowInfo& FlowInfo, float threshold_knickpoint, float threshold_knickpoint_length)
{
  cout << "Deprecated function, do not use it anymore. I am keeping it for a time just to check if someone still need it" << endl;
  exit(EXIT_FAILURE);
  /*// these are for extracting element-wise data from the channel profiles.
  //int abs_threshhold_knickpoint = abs (threshold_knickpoint);
  int this_node = 0;
  int segment_counter_knickpoint = 0; // count the number of knickpoints
  int segment_counter = 0; // count the number of segments
  map<int,float> this_segment_counter_knickpoint_map;
  map<int,int> this_segment_counter_map;
  map<int,int> this_segment_knickpoint_sign_map;
  map<int,int> this_segment_length_map;
  float last_M_chi, this_M_chi;
  float delta_m = 0; // difference between last and new m_chi
  int knickpoint_sign = 0; // sign of the knickpoint: + =1 and - = -1
  float temp_delta_m = 0; // debugging stuff
  float this_segment_length = 0;
  int last_node = 0;
  int n_nodes_segment = 0;
  float x1_temp =0;
  float y1_temp =0;
  float x2_temp =0;
  float y2_temp =0;
  bool same_channel;
  float max_knickpoint_value =0;
  int new_knickpoint_counter =0;


  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    this_node = node_sequence[0];
    last_M_chi =  M_chi_data_map[this_node];

    for (int n = 0; n< n_nodes; n++)
    {

      // set the nodes number and keep information about the previous one
      if(n>0)
      {
        last_node = this_node;
      }
      this_node = node_sequence[n];
      // Get the M_chi from the current node
      this_M_chi = M_chi_data_map[this_node];
      // increment the segment node counter
      n_nodes_segment++;

      // If the M_chi has changed, do stuffs
      if (this_M_chi != last_M_chi)
      {
        segment_counter++; // increment the  segment counter
        delta_m= last_M_chi/this_M_chi; // Ratio between last and new chi steepness
        if(delta_m<1){knickpoint_sign = -1;} else {knickpoint_sign = 1;} // Assign the knickpoint sign value
        //delta_m = abs(delta_m); // required to detect all the knickpoints
        //if(delta_m > temp_delta_m) {temp_delta_m = delta_m;} // debugging stuff
        // now checking if the difference of m_chi between two segment is not due to a channel change
          // first retrieving xy coordinates
        FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
        FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);

        // Then check if the distance betweenthe two is more than 2 nodes (distance between two points via pytagore or thing like this)
        if (sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp)) > (2*FlowInfo.get_DataResolution()))
        {
          same_channel = false;
        }
        // done



        // Check if the threshold is (over)reached
        //if(delta_m > abs_threshhold_knickpoint && same_channel) if we are are using the threshold value

        // SMM NOTE: WHAT IS GOING ON HERE?
        // Is this supposed to be if(same_channel) ?

        if(true) // useless thing
        {
          segment_counter_knickpoint++; // number of knickpoints
          this_segment_counter_knickpoint_map[this_node] = delta_m; // adding the knickpoint value
          if(delta_m>max_knickpoint_value){max_knickpoint_value=delta_m; cout << max_knickpoint_value <<"||"<<this_segment_counter_knickpoint_map[this_node] << endl; } // assign the new knickpoint max value
        }
        //this_segment_length = n_nodes_segment * FlowInfo.get_DataResolution(); // getting the length of the segment using the resolution * number of nodes
        same_channel = true; // Set back the same channel parameter to true
        // now assign the segment lenght to all the point of the segment
        for(int i = n_nodes_segment ; i > 0 ; i--)
        {
          this_segment_length_map[node_sequence[n - i]] =  this_segment_length;
        }

        last_M_chi = this_M_chi;
        this_segment_knickpoint_sign_map[this_node] = knickpoint_sign; // assign the segment polarity
        n_nodes_segment = 0; // reinitialyse the # of node for the next segment
        this_segment_length = 0; // same
      }
      else
      {
        // incrementing the segment length
        if (n>0)
        {
          // Now calculating the distance between the two last nodes

          // Retrieving the x/y coordinates of the last two nodes
          FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
          FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
          // calculate and increment the distance from the last node
          this_segment_length += sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
          //cout << "distance : " << this_segment_length;

        }
      }


      // Print the segment counter to the data map
      this_segment_counter_map[this_node]  = segment_counter;
    }





    // Just Ignore this process for now!
    if(threshold_knickpoint_length>0)
    {


      // calculating equivalent between nodes_of_knickpoints and nodes to be able to properly navigate inside the map
      int nodes_of_knickpoints [segment_counter_knickpoint];
      int knickpoint_id = 0;
      int last_knickpoint_id = 0;
      for(int i = 0; i< n_nodes; i++)
      {
        this_node = node_sequence[i];
        if(this_segment_counter_knickpoint_map.count(this_node)) // supposed to test if this map has a value assigned for this key
        {
          nodes_of_knickpoints[knickpoint_id] = this_node;
          last_knickpoint_id = knickpoint_id; // stock the last knickpoint_id
          knickpoint_id++;
        }
      }
      // setting up the calculation for the length threshold
      bool still_processing = true;
      bool still_processing_total = true;
      float old_max_knickpoint_value = max_knickpoint_value;
      float distance_to_process_down = threshold_knickpoint_length/2; // represents the threshold down the knickpoint
      float distance_to_process_up = threshold_knickpoint_length/2; // represents the threshold up the knickpoint
      int number_of_erase = 0;

      while (still_processing_total)
      {
        still_processing = true;
        vector <int> knickpoint_to_delete;
        knickpoint_id = 0;

        map<int, float>::iterator it = this_segment_counter_knickpoint_map.begin();
        while(still_processing && it != this_segment_counter_knickpoint_map.end())
        {
          //cout<<it->first<<" :: "<<it->second<<endl;
          if(it->second == max_knickpoint_value)
          {
            distance_to_process_down = threshold_knickpoint_length/2;
            distance_to_process_up = threshold_knickpoint_length/2;
            // let's know test the distance before and beyond this nodes
            for(int g = 1; distance_to_process_down>0 || distance_to_process_up >0 ;g++)
            {

              //Calculate the case down
              if(this_segment_counter_knickpoint_map.count(nodes_of_knickpoints[knickpoint_id - g]) && distance_to_process_down>0)
              {
                //cout << "down exists" << endl;
                this_node = nodes_of_knickpoints[knickpoint_id - g];
                last_node = nodes_of_knickpoints[knickpoint_id - g+1];
                FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
                float temp_distance = sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                if(temp_distance < distance_to_process_down)
                {
                  distance_to_process_down -= temp_distance;
                  knickpoint_to_delete.push_back(this_node);
                } else{distance_to_process_down = 0;}
              }
              else{distance_to_process_down = 0;}
              // calculate the case up
              if(this_segment_counter_knickpoint_map.count(nodes_of_knickpoints[knickpoint_id + g]) && distance_to_process_up> 0)
              {
                //cout << "up exists" << endl;
                this_node = nodes_of_knickpoints[knickpoint_id + g];
                last_node = nodes_of_knickpoints[knickpoint_id + g-1];
                FlowInfo.get_x_and_y_from_current_node(last_node, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node, x2_temp, y2_temp);
                float temp_distance = sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                if(temp_distance < distance_to_process_down)
                {
                  distance_to_process_down-=temp_distance;
                  knickpoint_to_delete.push_back(this_node);
                } else{distance_to_process_down = 0;}
              }
              else {distance_to_process_up = 0;}
            }
            // Now I have to erase everything I planned to
            for (vector<int>::iterator it2 = knickpoint_to_delete.begin() ; it2 != knickpoint_to_delete.end(); it2++)
            {
              this_segment_counter_knickpoint_map.erase(*it2);
              //cout << *it2 <<endl;
              number_of_erase++;
            }

            it = this_segment_counter_knickpoint_map.begin();
            knickpoint_to_delete.clear();

            // calculate the new maximum
            old_max_knickpoint_value = max_knickpoint_value;
            max_knickpoint_value = 0;
            for(int j = 0; j< segment_counter_knickpoint; j++)
            {
              //cout << "Am I reaching this point?? j: "<< j << endl;
              this_node = nodes_of_knickpoints[j];
              if(this_segment_counter_knickpoint_map.count(this_node))
              {
                if(this_segment_counter_knickpoint_map[this_node] > max_knickpoint_value && this_segment_counter_knickpoint_map[this_node] < old_max_knickpoint_value && j!=knickpoint_id )
                  {
                    max_knickpoint_value = this_segment_counter_knickpoint_map[this_node];

                  }
              }
            }
            // reset the internal loop
            still_processing = false;


          }
          else
          {

            it++;
            knickpoint_id++;
            cout <<"application of the length threshold: " << it->first << "||" << number_of_erase <<"||"<<old_max_knickpoint_value<<"||" <<max_knickpoint_value <<endl;
            if(it->first == nodes_of_knickpoints[last_knickpoint_id] || max_knickpoint_value == 0)
            {
              still_processing_total = false;
            }
          }
        }
      }
    }
    // alright, trying a new stuff for the length threshold this time











    


    if(false)
    {
      // now sorting and calculating the knickpoint values by length (hopefully)
      int nodes_of_knickpoints [segment_counter_knickpoint];
      map<int,vector <float> > temp_knickpoint_map;
      //map <int,float,float>
      int knickpoint_id = 0;
      for(int i = 0; i< n_nodes; i++)
      {
        this_node = node_sequence[i];
        if(this_segment_counter_knickpoint_map.count(this_node)) // supposed to test if this map has a value assigned for this key
        {
          nodes_of_knickpoints[knickpoint_id] = this_node;
          FlowInfo.get_x_and_y_from_current_node(this_node, x1_temp, y1_temp);
          float this_node_float = this_node;
          vector<float> temp_vecta (3);
          temp_vecta.push_back(this_node);
          temp_vecta.push_back(x1_temp);
          temp_vecta.push_back(y1_temp);
          temp_knickpoint_map[knickpoint_id] = temp_vecta;
          knickpoint_id++;
          cout << temp_vecta[1] << endl;

          //cout << this_node << " || " << nodes_of_knickpoints[knickpoint_id] << endl;
        }
      }
    }




    // this part begin the calculation for the length threshold to erase too dense areas and set new variables
    float distance_to_process_down = threshold_knickpoint_length/2; // represents the threshold
    float distance_to_process_up = threshold_knickpoint_length/2; // represents the threshold
    float distance_to_substract = 0; // temp variable that remove the distance already processed from the threshold
    bool still_processing = true; // tell the for loop when to reloop fromn the beginning
    bool still_processing_total = false; // tell the main loop when to stop
    int old_max = 0; // used to stock the previous maximum data
    new_knickpoint_counter = segment_counter_knickpoint; // will return the number of knickpoints after deletion of the oldest
    //number_of_nodes_to_investigate_length = threshold_knickpoint_length;
    // beginning the calculation
    cout << "beginning the length calculation stuffs" << endl;




    // old method, will erase if I found another way
    while(still_processing_total) // This will be shutted down when the last node will have been processed
    {
      still_processing = true; // Required to trigger next loop
      for(int i = 0; i< segment_counter_knickpoint && still_processing; i++)
      {
        this_node = nodes_of_knickpoints[i]; // go through the knickpoints
        if(this_segment_counter_knickpoint_map.count(this_node)) // check if this node still exists
        {
          if(  this_segment_counter_knickpoint_map[this_node]==max_knickpoint_value) // Check if it is the current maximum value
          {
            for(int g = 1; distance_to_process_down > 0 || distance_to_process_up > 0 ; g++) // If so, test the adjacent nodes in order to delete the required ones
            {
              // now getting the coordinates of the wanted nodes
              cout<<"aue coute"<<endl;
              if(this_segment_length_map.count(this_node-g) && this_segment_length_map.count(this_node-g+1))
              {

                FlowInfo.get_x_and_y_from_current_node(this_node-g, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node-g+1, x2_temp, y2_temp);


                // Check if it exists and if it is on the same river
                if(this_segment_counter_knickpoint_map.count(this_node-g) && distance_to_process_down >0)
                {
                  if(sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp))<= distance_to_process_down)
                  {
                    this_segment_counter_knickpoint_map.erase(this_node-g);
                    this_segment_knickpoint_sign_map.erase(this_node-g);
                    new_knickpoint_counter--;
                    cout << "something to test blablabla" << endl;
                    distance_to_process_down -= sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                  }
                  else
                  {
                    distance_to_process_down = 0;
                  }
                }
              }
              else
              {
                distance_to_process_down = 0;
              }
              if(this_segment_length_map.count(this_node+g) && this_segment_length_map.count(this_node+g-1))
              {
                FlowInfo.get_x_and_y_from_current_node(this_node+g, x1_temp, y1_temp);
                FlowInfo.get_x_and_y_from_current_node(this_node+g-1, x2_temp, y2_temp);
                if(this_segment_counter_knickpoint_map.count(this_node+g) && distance_to_process_up>0)
                {
                  if(sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp))<= distance_to_process_up)
                  {
                    cout << "before erase :" << this_segment_knickpoint_sign_map[this_node+g] << endl;
                    this_segment_counter_knickpoint_map.erase(this_node+g);
                    this_segment_knickpoint_sign_map.erase(this_node+g);
                    new_knickpoint_counter--;
                    cout << "after erase :" << this_segment_knickpoint_sign_map[this_node+g] << endl;
                    distance_to_process_up -= sqrt((x2_temp - x1_temp)*(x2_temp - x1_temp)+(y2_temp - y1_temp)*(y2_temp - y1_temp));
                  }
                  else
                  {
                    distance_to_process_up=0;
                  }
                }
              }
              else
              {
                distance_to_process_up=0;
              }
            }

            old_max = max_knickpoint_value;
            max_knickpoint_value = 0;
            for(int j = 0; j< segment_counter_knickpoint; j++)
            {
              //cout << "Am I reaching this point?? j: "<< j << endl;
              this_node = nodes_of_knickpoints[j];
              if(this_segment_counter_knickpoint_map.count(this_node))
              {
                if(this_segment_counter_knickpoint_map[this_node] > max_knickpoint_value && this_segment_counter_knickpoint_map[this_node] <= old_max)
                  {
                    max_knickpoint_value =this_segment_counter_knickpoint_map[this_node];
                    //if(this_segment_counter_knickpoint_map[this_node] == max_knickpoint_value) {cout<< "biiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiite"<< endl;}
                  }
              }
            }
            still_processing = false;
            distance_to_process_up = threshold_knickpoint_length/2;
            distance_to_process_down = threshold_knickpoint_length/2;
          }
        }
        if(i == segment_counter_knickpoint-1)
        {
          still_processing = false;
          still_processing_total = false;
        }
        else{still_processing = true;}
      }
    }
  }


  cout << "segment_counter_knickpoint is   " << new_knickpoint_counter << "/" << segment_counter << " delta max is " << temp_delta_m << endl;
  // print everything in the public/protected variables
  segment_counter_knickpoint_map = this_segment_counter_knickpoint_map;
  segment_knickpoint_sign_map = this_segment_knickpoint_sign_map;
  segment_length_map = this_segment_length_map;
  */
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_previous_mchi_for_all_sources(LSDFlowInfo& Flowinfo)
{
  // setting the variables for extracting the knickpoints
  
   // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    int working_source = source_keys_map[node_sequence[0]]; // This is the working source key, the one you already have extracted the information
    int current_source = source_keys_map[node_sequence[0]]; // this is the currently tested source_key that will become the working key is different than previous key
    int starting_node_of_source_key = get_ending_node_of_source(Flowinfo ,working_source); // This store the starting node of the river with this source_key
    
    // Now getting the receiving node of the river with this source
    int receiving_node_of_source_key, temp_row, temp_col;
    Flowinfo.retrieve_receiver_information(starting_node_of_source_key, receiving_node_of_source_key, temp_row,temp_col);

    //Finally getting the 
    float m_chi_receiving_river = M_chi_data_map[receiving_node_of_source_key];
    
    // Creating temp data_map to save everything
    map<int,int> this_map_source_key_receiver; // 
    map<int,float> this_map_source_key_receiver_mchi;
    this_map_source_key_receiver[working_source] = source_keys_map[receiving_node_of_source_key];
    this_map_source_key_receiver_mchi[working_source] = m_chi_receiving_river;

    // done initializing, let's do it for all the rivers

    for (int n = 0; n< n_nodes; n++)
    {
      current_source = source_keys_map[node_sequence[n]];
      //cout << source_keys_map[node_sequence[starting_node_of_source_key]] << "||" << node_sequence[source_keys_map[receiving_node_of_source_key]] << endl;
      if(current_source != working_source && current_source != -9999)
      {
        // cout << "changing sources" << endl;
        working_source = current_source;
        starting_node_of_source_key = get_ending_node_of_source(Flowinfo, working_source);
        Flowinfo.retrieve_receiver_information(starting_node_of_source_key, receiving_node_of_source_key, temp_row,temp_col);
        m_chi_receiving_river = M_chi_data_map[receiving_node_of_source_key];
        this_map_source_key_receiver[working_source] = source_keys_map[receiving_node_of_source_key];
        this_map_source_key_receiver_mchi[working_source] = m_chi_receiving_river;
      }  

    }

    // generalizing the maps
    map_source_key_receiver = this_map_source_key_receiver;
    map_source_key_receiver_mchi = this_map_source_key_receiver_mchi;


  }

}

LSDRaster LSDChiTools::prefilter_river_topography(LSDRaster& filled_topography, LSDFlowInfo& FlowInfo, double lambda_TVD)
{
  // cout << "Work in progress, But first Lunch" << endl;
    // preparing the needed iterators
  vector<int> ordered_source_keys;
  Array2D<float> next_topo = filled_topography.get_RasterData();
  set_map_of_source_and_node(FlowInfo,5);


  // getiing the SK in the right order to make sure I am processing it bottom to top
  int n_nodes = (node_sequence.size());
  int last_SK = -9999;
  int this_SK = 0;
  for(int n=0; n<n_nodes ;n++)
  {
    this_SK = source_keys_map[node_sequence[n]];
    if(this_SK != last_SK)
    {
      ordered_source_keys.push_back(this_SK);
      last_SK = this_SK;
    }
  }

  // Initializing some variables
  vector<int> vecnode;

  // Looping through all the sources key
  float this_base_elevation = 0, retrend_bae_elevation = 0;
  for(vector<int>::iterator tit = ordered_source_keys.begin(); tit != ordered_source_keys.end(); tit++)
  {

    this_SK = *tit;
    vecnode = map_node_source_key[this_SK]; // getting the vector of river node to check
    // exit(EXIT_FAILURE);
    
    vector<double> detrend(vecnode.size()), denoised(vecnode.size()), retrend(vecnode.size());
    int fnode = vecnode.back();
    int recnode,rrow, rcol,row,col;

    // First step is to determine our boundary condition -> the base elevation that will be retrended
    FlowInfo.retrieve_current_row_and_col(fnode,row,col);
    FlowInfo.retrieve_receiver_information(fnode,recnode,rrow, rcol);
    if(fnode == recnode || filled_topography.get_data_element(row,col) == NoDataValue)
    {
      // If baselevel
      this_base_elevation = filled_topography.get_data_element(row,col);

    }
    else
    {
      // Else getting the receiver
      this_base_elevation = filled_topography.get_data_element(rrow,rcol);
      retrend_bae_elevation = next_topo[row][col];

    }

    // Going through the vector of nodes
    size_t rit_a = vecnode.size() - 1;
    int last_node; 
    float last_elev;
    for(vector<int>::reverse_iterator ti2 = vecnode.rbegin(); ti2 != vecnode.rend(); ++ti2)
    {
      int this_node = *ti2;

      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      // FlowInfo.retrieve_receiver_information(this_node,recnode,rrow, rcol);

      if(ti2 == vecnode.rbegin())
      {
        detrend[rit_a] = filled_topography.get_data_element(row,col) - this_base_elevation;
      }
      else
      {

        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        detrend[rit_a] = filled_topography.get_data_element(row,col) - last_elev;
        // cout << "D:" << detrend[rit_a] << endl;
        // cout << "LE:" << last_elev << endl;
      }
      last_node = this_node;
      last_elev = filled_topography.get_data_element(row,col);
      rit_a--;
    }
    // exit(EXIT_FAILURE);


    // Denoising
    if(vecnode.size()>1)
    {
      denoised = TV1D_denoise_v2(detrend, lambda_TVD);
    }
    else
    {
      denoised = detrend;
    }

    // Now retrending

    for(int rit_b = int(vecnode.size()) - 1 ; rit_b >= 0; --rit_b)
    {

      if(rit_b == vecnode.size() - 1)
      {
        // cout << "Premier node: " << this_base_elevation << endl;
        int tnode = vecnode[rit_b], trow=0,tcol=0;
        FlowInfo.retrieve_receiver_information(tnode,recnode,trow, tcol);
        retrend_bae_elevation = next_topo[trow][tcol];
        retrend[rit_b] = retrend_bae_elevation + denoised[rit_b];
      }
      else
      {
        retrend[rit_b] = retrend[rit_b+1] + denoised[rit_b];
      }
      int this_node = vecnode[rit_b];


      FlowInfo.retrieve_current_row_and_col(this_node,row,col);

      if(row>=0 && col>=0 && row<NRows && col<NCols)
      {

        next_topo[row][col] = retrend[rit_b];
      }
      else{cout<<"DEBUGWARNING::Unreal Node here ?!" << endl;}
    }
    cout << "done" << endl;

  }

  
  LSDRaster ThisRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, next_topo, GeoReferencingStrings);

  return ThisRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// Rather than recreating, changing and messing wiht knickpoints functions during developent,
// Eveything will now be controlled from this function calling the adapted and up-to-date function
// Leaving me time to develop what will be the final cleanest one and making easier the subdivision
// in loads of little functions rather than one big script-like one. OBJECT ORIENTED POWER ˁ˚ᴥ˚ˀ
// BG 
void LSDChiTools::ksn_knickpoint_automator(LSDFlowInfo& FlowInfo, string OUT_DIR, string OUT_ID, float MZS_th, float lambda_TVD,int stepped_combining_window,int window_stepped, float n_std_dev, int kp_node_search)
{

  cout << "Getting ready for the knickpoint detection algorithm ...";
  // The first preprocessing step is to preselect the river we want to process
  // Potentially data selection function to be added here, exempli gratia lenght threshold for tributaries
  // this first function fill a map[source key] = vector<node for this rive including the receiver node>
  set_map_of_source_and_node(FlowInfo,5);

  // We need to set the d(Segmented_elevation) to get the stepped knickpoints
  derive_the_segmented_elevation();
  cout << " OK" << endl ;

  // vector<int> dsaf = FlowInfo.get_DonorStackVector();

  // We preprocessed the metrics we need, let's get into the Data!
  // cout << "Experimentaion on TVD on elevation here" << endl;
  // TVD_on_segelev(FlowInfo);
  // cout << " JUMANGI" << endl ;


  cout << " Denoising the ksn or mchi and the differential segmenting elevation (Total Variation Denoising adapted from Condat, 2013) ..." << endl;
  // Applying the Total_variation_denoising on m_chi.
  // This is really efficient Algorithm, I am denoising the b_chi as well, for testing purposes
  TVD_on_my_ksn(lambda_TVD);
  cout << " OK" << endl ;


  cout << " Extracting general metrics for rivers ...";
  // This will increment maps with source keys as key and various metrics such as river length, Chi lenght...
  compute_basic_matrics_per_source_keys(FlowInfo);
  cout << " OK" << endl ;



  // main function that increment the map_of_knickpoints by detecting the changes in ksn within rivers
  // /!\ Contain a cout statement
  cout << "Detecting raw ksn knickpoints and stepped knickpoints for each source ..." << endl;
  ksn_knickpoint_detection_new(FlowInfo);
  cout << " OK" << endl;


  // As proud as I was to have implement that, I unfortunately don't need it anymore...
  // // Now dealing with outlier detection
  // // first calculating the KDE
  // cout << "Kernel Density Estimation per river (Deprecated, but I keep it for some debugging purposes) ...";
  // ksn_kp_KDE();
  // cout << " OK" << endl ;
  // I just keep it for latter purposes, RIP AGU method.

  // Processing the knickpoints to combine the composite knickpoints
  cout << "Combining ksn knickpoints ..." << endl;
  ksn_knickpoints_combining(FlowInfo,kp_node_search);
  cout << " OK" << endl ;

  cout << "Getting the stepped_knickpoints ..." << endl;
  stepped_knickpoints_detection_v2(FlowInfo,window_stepped,n_std_dev);
  // stepped_knickpoints_combining(FlowInfo, stepped_combining_window);
  cout << " OK" << endl ;


  // Ok let's detect oultliers here
  cout << "Generating some stats ...";
  ksn_knickpoint_outlier_automator(FlowInfo, MZS_th);
  cout << " OK" << endl ;


  //printing the raw ksn knickpoint file
  string this_name = OUT_DIR + OUT_ID + "_ksnkp_raw.csv";
  cout << "Printing data into csv files ..." << endl;
  cout << "Raw..." << endl;
  print_raw_ksn_knickpoint(FlowInfo, this_name);
  this_name = OUT_DIR + OUT_ID + "_ksnkp_SK.csv";
  cout << "SK..." << endl;
  print_bandwidth_ksn_knickpoint(this_name);
  this_name = OUT_DIR + OUT_ID + "_ksnkp_mchi.csv";
  cout << "mchi/ksn..." << endl;
  print_mchisegmented_knickpoint_version(FlowInfo, this_name);
  this_name = OUT_DIR + OUT_ID + "_ksnkp.csv";
  cout << "knickpoints..." << endl;
  print_final_ksn_knickpoint(FlowInfo, this_name);
  cout << " OK" << endl ;

  // Old to keep
  // this_name = OUT_DIR + OUT_ID + "_ksnkp_SEGELEV.csv";
  // print_mchisegmented_knickpoint_version_test_on_segmented_elevation(FlowInfo, this_name);


}


// Rather than recreating, changing and messing wiht knickpoints functions during developent,
// Eveything will now be controlled from this function calling the adapted and up-to-date function
// Leaving me time to develop what will be the final cleanest one and making easier the subdivision
// in loads of little functions rather than one big script-like one. OBJECT ORIENTED POWER ˁ˚ᴥ˚ˀ
// BG 
void LSDChiTools::ksn_knickpoint_automator_no_file(LSDFlowInfo& FlowInfo, float MZS_th, float lambda_TVD,int stepped_combining_window,int window_stepped, float n_std_dev, int kp_node_search)
{

  cout << "Getting ready for the knickpoint detection algorithm ...";
  // The first preprocessing step is to preselect the river we want to process
  // Potentially data selection function to be added here, exempli gratia lenght threshold for tributaries
  // this first function fill a map[source key] = vector<node for this rive including the receiver node>
  set_map_of_source_and_node(FlowInfo,5);

  // We need to set the d(Segmented_elevation) to get the stepped knickpoints
  derive_the_segmented_elevation();
  cout << " OK" << endl ;

  // vector<int> dsaf = FlowInfo.get_DonorStackVector();

  // We preprocessed the metrics we need, let's get into the Data!
  // cout << "Experimentaion on TVD on elevation here" << endl;
  // TVD_on_segelev(FlowInfo);
  // cout << " JUMANGI" << endl ;


  cout << " Denoising the ksn or mchi and the differential segmenting elevation (Total Variation Denoising adapted from Condat, 2013) ..." << endl;
  // Applying the Total_variation_denoising on m_chi.
  // This is really efficient Algorithm, I am denoising the b_chi as well, for testing purposes
  TVD_on_my_ksn(lambda_TVD);
  cout << " OK" << endl ;


  cout << " Extracting general metrics for rivers ...";
  // This will increment maps with source keys as key and various metrics such as river length, Chi lenght...
  compute_basic_matrics_per_source_keys(FlowInfo);
  cout << " OK" << endl ;



  // main function that increment the map_of_knickpoints by detecting the changes in ksn within rivers
  // /!\ Contain a cout statement
  cout << "Detecting raw ksn knickpoints and stepped knickpoints for each source ..." << endl;
  ksn_knickpoint_detection_new(FlowInfo);
  cout << " OK" << endl;


  // As proud as I was to have implement that, I unfortunately don't need it anymore...
  // // Now dealing with outlier detection
  // // first calculating the KDE
  // cout << "Kernel Density Estimation per river (Deprecated, but I keep it for some debugging purposes) ...";
  // ksn_kp_KDE();
  // cout << " OK" << endl ;
  // I just keep it for latter purposes, RIP AGU method.

  // Processing the knickpoints to combine the composite knickpoints
  cout << "Combining ksn knickpoints ..." << endl;
  ksn_knickpoints_combining(FlowInfo,kp_node_search);
  cout << " OK" << endl ;

  cout << "Getting the stepped_knickpoints ..." << endl;
  stepped_knickpoints_detection_v2(FlowInfo,window_stepped,n_std_dev);
  // stepped_knickpoints_combining(FlowInfo, stepped_combining_window);
  cout << " OK" << endl ;


  // Ok let's detect oultliers here
  cout << "Generating some stats ...";
  ksn_knickpoint_outlier_automator(FlowInfo, MZS_th);
  cout << " OK" << endl ;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// New function for the knickpoint detection
// save the difference and ratio between m_chi values of each segments
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::ksn_knickpoint_detection_new(LSDFlowInfo& FlowInfo)
{

  // preparing the needed iterators
  map<int,vector<int> >::iterator SK;
  // Initializing some variables
  int this_SK;
  vector<int> vecnode;

  // Looping through all the sources key
  for(SK = map_node_source_key.begin(); SK != map_node_source_key.end(); SK++)
  {
    this_SK = SK->first;
    // cout << "this SK:" << this_SK << endl;
    vecnode = SK->second; // getting the vector of river node to check
    // cout << "vecnode gottend:" << vecnode.size() << endl;

    if(vecnode.size()>5)
      ksn_knickpoint_raw_river(this_SK,vecnode);
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Detect the knickpoint in one river and increment the global map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::ksn_knickpoint_raw_river(int SK, vector<int> vecnode)
{
  // Setting the iterator(s)
  vector<int>::iterator node = vecnode.begin(); // first node of the river -> the source

  // Setting function variables
  int last_node = *node; // last node is the first node
  node++; // switching to the second node
  int this_node = *node; // this node is the second one 
  // vector that will contain the nodes having a knickpoint
  vector<int> vecdif;
  // Bunch of floats
  //float dkdc = 0;     // Don't seem to need these (SMM)
  float dchi = 0;     // Don't seem to need these (SMM)
  float dksn = 0;
  float this_ksn = TVD_m_chi_map[this_node];
  float last_ksn = TVD_m_chi_map[last_node]; // Setting last and this ksn
  // chi stuff
  float this_chi = chi_data_map[this_node];
  float last_chi = chi_data_map[last_node];
  // elevation to get like

    //REPLY TO REVIEW VARIABLES
  // vector<int> node_to_imp;
  // vector<double> temp_cx; 

  // Looping through the nodes from the second one
  for( ; node != vecnode.end(); node++) // the first ";" is normal: it states that I have no initial conditions 
  {
    // initializing the variables for this run
    this_node = *node;
    this_ksn = TVD_m_chi_map[this_node];
    this_chi = chi_data_map[this_node];
    // dsegelev = segelev_diff[this_node];



    // if ksn has change, Implementing a raw knickpoint, quantifying it with delta ksn
    if((this_ksn != last_ksn) && this_ksn != -9999 && last_ksn != -9999 )
    {

      // deta ksn from bottom to top
      dksn = last_ksn - this_ksn;

      dchi = this_chi - last_chi;
      // saving the value in the global raw_ksn_kp map
      raw_ksn_kp_map[this_node] =  dksn;
      // RPLY TO REVIEW CX STUFF
      Cx_TVD_map[this_node] = dksn / dchi;
      // saving the node to get the vector of ksn knickpoints per source, for later grouping purpose for example
      vecdif.push_back(this_node);

      
    }
    // setting the next last variables
    last_node = this_node;
    last_ksn = this_ksn;
    last_chi = this_chi;

  }
  // vector<double> temp_cx_2;

  // //Reply to review
  // if(node_to_imp.size()>1)
  // {
  //   temp_cx_2 = TV1D_denoise_v2(temp_cx, 1);
  //   for(size_t urugay=0;urugay<node_to_imp.size();urugay++)
  //   {
  //     this_node = node_to_imp[urugay];
  //     Cx_2TVD_map[this_node] = 
  //   }
  // }


  // implementing the global map
  map_node_source_key_kp[SK] = vecdif;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Detection of stepped kncikpoints - NEW
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::stepped_knickpoints_detection_v2(LSDFlowInfo& Flowinfo, int window, float n_std_dev)
{
    // preparing the needed iterators
  map<int,vector<int> >::iterator SK;
  // Initializing some variables
  // int this_SK;             // This doens't seem to be used (SMM)
  vector<int> vecnode;
  int HW = int(window/2);

  // Looping through all the sources key
  for(SK = map_node_source_key.begin(); SK != map_node_source_key.end(); SK++)
  {
    
    // the source key
    //this_SK = SK->first;       // This doens't seem to be used (SMM)
    
    // the vector containing the river nodes
    vecnode = SK->second;
    // first check if the size is over the window
    if(int(vecnode.size()) > window)
    { 
      // first step is to get the windowed stats, the moving window will gather informations across the vector of nodes 
      map<string,vector<float> > these_stats = get_windowed_stats_for_knickpoints(vecnode,HW);

      for(size_t it = 0; it< vecnode.size(); it++)
      {
        int this_node = vecnode[it];
        float this_std = these_stats["std_dev"][it], this_mean = these_stats["mean"][it];
        mean_for_kp[this_node] = this_mean;
        std_for_kp[this_node] = this_std;

        if(segelev_diff[this_node] >= (n_std_dev * this_std))
        {
          kp_segdrop[this_node] = segelev_diff[this_node];
        }
      }
    }
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Apply a moving window on the river nodes to get the mean and std-dev on each nodes
// HW is the HalfWidth
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<string,vector<float> > LSDChiTools::get_windowed_stats_for_knickpoints(vector<int> vecnode,int HW)
{
  // first get the general vector of value
  vector<float> vecval, window_vecval_mean, window_vecval_std;
  int window = HW*2;

  // looping through the node of the river and getting the corresponding delta segmented elevation
  for(vector<int>::iterator it = vecnode.begin(); it != vecnode.end(); it ++)
  {
    int this_node = *it;
    if(chi_data_map[this_node] > 0){
      vecval.push_back(segelev_diff[this_node]);
    }
  }
  // Vecval is now implemented with the values

  // Now moving window along the delta segelev profile
  for(size_t it = 0;it<vecnode.size(); it++)
  {
    // we are in the first nodes boundary: we cannot loop before so the same value will be applied everywhere in this case
    if(int(it)<HW)
    {
      vector<float> this_vecval;
      // getting the value for the window
      for(int o = 0; o < window ; o++ )
      {
        this_vecval.push_back(vecval[o]);
      }
      // getting the stats
      float this_mean = get_mean(this_vecval);
      // implemeting
      window_vecval_mean.push_back(this_mean);
      float this_std = get_standard_deviation(this_vecval,this_mean);
      window_vecval_std.push_back(this_std);

      this_vecval.clear();    
    }
    // we are in the middle of the river nodes
    else if (it < (vecnode.size()-HW))
    {
      vector<float> this_vecval;
      for(size_t o = it; o < (it+window) ; o++ )
      {
        this_vecval.push_back(vecval[o]);
      }
      float this_mean =get_mean(this_vecval);
      window_vecval_mean.push_back(this_mean);
      float this_std = get_standard_deviation(this_vecval,this_mean);
      window_vecval_std.push_back(this_std);

      this_vecval.clear();    
    }
    //we are in the last part
    else if (it >= (vecnode.size()-HW))
    {
      vector<float> this_vecval;
      for(size_t o = (vecnode.size()-window); o < (vecnode.size()) ; o++ )
      {
        this_vecval.push_back(vecval[o]);
      }
      float this_mean = get_mean(this_vecval);
      window_vecval_mean.push_back(this_mean);
      float this_std = get_standard_deviation(this_vecval,this_mean);
      window_vecval_std.push_back(this_std);

      this_vecval.clear();    
    }

  }

  // generating the outputs
  map<string,vector<float> > map_of_stats;
  map_of_stats["mean"] = window_vecval_mean;
  map_of_stats["std_dev"] = window_vecval_std;

  return map_of_stats;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Detection of stepped kncikpoints - old version
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::raw_stepped_knickpoint_detection(int SK, vector<int> vecnode)
{
  // Setting the iterator(s)
  vector<int>::iterator node = vecnode.begin(); // first node of the river -> the source

  // Setting function variables
  int last_node = *node; // last node is the first node
  node++; // switching to the second node
  int this_node = *node; // this node is the second one 
  // vector that will contain the nodes having a knickpoint
  vector<int> vecdif;
  // Bunch of floats
  float this_TVD_b_chi = TVD_b_chi_map[this_node], last_TVD_b_chi = TVD_b_chi_map[last_node]; // Setting last and this segmented elevation change

  // Looping through the nodes from the second one
  for( ; node != vecnode.end(); node++) // the first ";" is normal: it states that I have no initial conditions 
  {
    // initializing the variables for this run
    this_node = *node;
    this_TVD_b_chi = TVD_b_chi_map[this_node];
   
    // if b_chi has change, Implementing a raw knickpoint and calculating the d|ksn|/dchi
    if(this_TVD_b_chi != last_TVD_b_chi || TVD_m_chi_map[last_node] != TVD_m_chi_map[this_node])
    {
      raw_delta_segelev_from_TVDb_chi[this_node] = segelev_diff[this_node]; 
      // testing something here
      // float this_segdiff = (TVD_m_chi_map[last_node] * chi_data_map[last_node] + TVD_b_chi_map[last_node]) - (TVD_m_chi_map[this_node] * chi_data_map[this_node] + TVD_b_chi_map[this_node]);
      // raw_delta_segelev_from_TVDb_chi[this_node] = this_segdiff;
      vecdif.push_back(this_node);
    }
    // setting the next last variables
    last_node = this_node;
    last_TVD_b_chi = this_TVD_b_chi;
  }

  // implementing the global map
  map_node_source_key_kp_stepped[SK] = vecdif;
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  Group the adjacent local knickpoints                  =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

void LSDChiTools::ksn_knickpoints_combining(LSDFlowInfo& Flowinfo, int kp_node_search)
{
  // this function will combine the composite knickpoints

  // First looping through the source keys
  map<int,vector<int> >::iterator henri;
  // saaving an ID for each knickpoints
  id_kp = 0;

  for(henri = map_node_source_key_kp.begin(); henri != map_node_source_key_kp.end(); henri ++)
  {
    // Now looping through the node for each rivers
    int this_SK = henri -> first;
    // Getting the vector of node per river AND vector of ksn knickpoints nodes per river
    vector<int> vecnode_kp = henri->second, vecnode_river = map_node_source_key[this_SK];
    if(vecnode_kp.size()>0)
    {
      // getting the groups of vector depending on your combination window
      vector<vector<int> > grouped_kp = group_local_kp(vecnode_kp,vecnode_river,Flowinfo, kp_node_search);

      // We have the group of vector now lets run through it to get the requested values
      for(vector<vector<int> >::iterator vlad = grouped_kp.begin(); vlad != grouped_kp.end(); vlad ++)
      {
        // Investigating the first node
        vector<int> this_vecnode = *vlad;

        // Easy case: non composite knickpoint, let's just record the same info thatn the raw detection
        if(this_vecnode.size() == 1)
        {
          int this_node = this_vecnode[0];// node of the knickpoint
          ksn_kp_map[this_node] = raw_ksn_kp_map[this_node]; // delta ksn of the knickpoint
          sharpness_ksn_length[this_node] = 0; // sharpness = 0 as the knickpoint is a point
          ksn_extent[this_node] = make_pair(this_node, this_node); // the "extent" nodes are the same

          // Getting the coordinate of the centroid of the kp -> just regular x
          float this_x = 0,this_y = 0;
          Flowinfo.get_x_and_y_from_current_node(this_node, this_x, this_y);
          ksn_centroid[this_node] = make_pair(this_x,this_y); // centroid does not move as well

          // finally getting the ID of my kp and its location (It is not relevant in this case, but we need to implement the maps)
          ksn_kp_ID[this_node] = id_kp;
          nearest_node_centroid_kp[this_node] = this_node;
          flow_distance_kp_centroid_map[this_node] = flow_distance_data_map[this_node];
          id_kp ++;
        }

        else if(this_vecnode.size() > 1)
        {
          // harder case: several knickpoints, let's go step by step
          // the identifying node of the kp is the first one
          int this_node = this_vecnode[0];
          
          // summing the segelev and dksn for this group of node
          ksn_kp_map[this_node] = get_dksn_from_composite_kp(this_vecnode); // gobal ksn value of the knickpoint, simple sum

          // Sharpness is the width of the combined knickpoints
          sharpness_ksn_length[this_node] = get_kp_sharpness_length(this_vecnode, Flowinfo);
          ksn_extent[this_node] = make_pair(this_vecnode[0], this_vecnode.back()); // the extent nodes are the extreme of the grouping vector


          // The neirest node from the centre of the knickpoint (in regards to chi distance)
          int nearnode = get_ksn_centroid_coordinates(Flowinfo, this_vecnode, vecnode_river,ksn_kp_map[this_node]); // get the weighted x and y of the centroid.

          flow_distance_kp_centroid_map[this_node] = flow_distance_data_map[nearnode];
          // getting the Flow dist at the weighted distance
          float this_x = 0,this_y = 0;
          Flowinfo.get_x_and_y_from_current_node(nearnode, this_x, this_y);
          // implementing the global maps
          ksn_centroid[this_node] = make_pair(this_x,this_y);
          nearest_node_centroid_kp[this_node] = nearnode;
          ksn_kp_ID[this_node] = id_kp;
          id_kp++;
        }
      }
    }

  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  get the delta ksn value for a composite knickpoint    =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

float LSDChiTools::get_dksn_from_composite_kp(vector<int> vecnode)
{

  float out_value = 0;
  for(vector<int>::iterator bob = vecnode.begin(); bob!= vecnode.end(); bob++)
  {
    int this_node = *bob;
    out_value += raw_ksn_kp_map[this_node];
  }

  return out_value;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  get the delta segelev value for a composite knickpoint - OLD    =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float LSDChiTools::get_dseg_drop_from_composite_kp(vector<int> vecnode)
{
  // I'll get the segelev augmentation NORMALIZED by the chi distance
  // I have to find a way to normalize
  float out_value = 0; //, chi_min = 0, chi_max =0;
  for(vector<int>::iterator bob = vecnode.begin(); bob!= vecnode.end(); bob++)
  {
    int this_node = *bob;
    out_value += raw_delta_segelev_from_TVDb_chi[this_node];
  }

  return out_value;
}

float LSDChiTools::get_kp_sharpness_length(vector<int> vecnode, LSDFlowInfo& Flowinfo)
{
  float total_distance = 0;
  int last_node = vecnode[0];

  for(vector<int>::iterator gog = vecnode.begin(); gog!= vecnode.end(); gog++)
  {
    int this_node = *gog;
    total_distance += abs(flow_distance_data_map[last_node] - flow_distance_data_map[this_node]);
  }

  return total_distance;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  get the centroid node of grouped knickpoints          =
//                  Current version                       =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDChiTools::get_ksn_centroid_coordinates(LSDFlowInfo& Flowinfo, vector<int> vecnode,vector<int> vecnode_river, float total_dksn)
{

  pair<float,float> out_pair; // x,y coordinates

  vector<int> this_vecnode_river = get_vecnode_river_from_extent(vecnode[0], vecnode.back(), vecnode_river);

  // first, let's get the chi, the dksn and the weighted distance of each nodes
  vector<float> chi_vec , dksn_vec;
  
  // looping through each nodes to group, and getting the info in the same order that the original vector of nodes
  for(size_t iter = 0; iter < this_vecnode_river.size(); iter++)
  {
    chi_vec.push_back(chi_data_map[this_vecnode_river[iter]]);
    if (raw_ksn_kp_map.count(this_vecnode_river[iter]) != 0)
    {
      dksn_vec.push_back(abs(raw_ksn_kp_map[this_vecnode_river[iter]])); // I am using the raw map, that host the values before combining
    }
    else
    {
      dksn_vec.push_back(0);
    }
  }


  // Applying a barycenter-like centering of the data to get the chi barycentre
  float M = total_dksn, chi_center = 0, sum_of_mi_ri = 0;
  // Summing the weighted chi distance by dksn

  for(size_t iter = 0; iter < this_vecnode_river.size(); iter++)
  {
    float mi = dksn_vec[iter], ri = chi_vec[iter];
    sum_of_mi_ri += mi*ri;
    // cout << sum_of_mi_ri << " || " << ri << " || " << mi <<  " || " << M << endl;
  }

  // applying the weight
  chi_center = abs(sum_of_mi_ri/M);
  // cout << chi_center << endl;

  // *cymbal and trumbet noise* we have our chi center

  // then looping from bottom to top of nodes to get the boundary nodes of this last
  bool found_it = false;
  int ninf = 0, nsup = 0, this_node = 0, last_node = 0 ;
  for(vector<int>::iterator hibou = this_vecnode_river.begin(); found_it == false; hibou ++)
  {
    this_node = *hibou;
    // cout << this_node << " || " << this_vecnode_river.back() << endl ;
    if((chi_data_map[this_node] > chi_center )|| ( this_node == this_vecnode_river.back()))
    {
      found_it = true;
      ninf = last_node;
      nsup = this_node;
    }
    last_node = this_node;
  }

  // now getting the distance between nodes
  int centre_node = 0;

  float up_chi_diff = abs(chi_data_map[nsup] - chi_center);
  float down_chi_diff = abs(chi_data_map[ninf] - chi_center);

  if(up_chi_diff<=down_chi_diff)
  {
    centre_node = nsup;
  }
  else
  {
    centre_node = ninf;
  }

  return centre_node;

}


vector<int> LSDChiTools::get_vecnode_river_from_extent(int first_node, int last_node, vector<int> vecnode_river)
{

  // first_node is bottom node and last node is the upper node inb term of elevation
        // cout << endl << elev_data_map[vecnode_river[0]] << " || " << elev_data_map[vecnode_river.back()] << endl;
        // cout <<  elev_data_map[first_node] << " || " << elev_data_map[last_node] << endl;


        // exit(EXIT_FAILURE);

  vector<int> out_node;
  bool pushit = false;

  // looping through the river nodes
  for(vector<int>::iterator uter = vecnode_river.begin(); uter != vecnode_river.end(); uter++)
  {
    int this_node = *uter;
    // Checking if I need to begin the push-in
    if(this_node == first_node){pushit  = true;}
    // push when relevant
    if(pushit){out_node.push_back(this_node);}
    // Checking if I need to stop pushing
    if(this_node == last_node){pushit = false;}
  }
  // cout << endl << out_node.size() << endl;


  return out_node;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  This method combine the main stepped variation        =
//                  Old version                           =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

void LSDChiTools::stepped_knickpoints_combining(LSDFlowInfo& Flowinfo, int kp_node_search)
{

// this function will combine the composite knickpoints

  // First looping through the source keys
  map<int,vector<int> >::iterator henri;

  for(henri = map_node_source_key_kp_stepped.begin(); henri != map_node_source_key_kp_stepped.end(); henri ++)
  {
    // Now looping through the node for each rivers
    int this_SK = henri -> first;
    vector<int> vecnode_kp = henri->second, vecnode_river = map_node_source_key[this_SK];
    // cout << this_SK << endl;
    if(vecnode_kp.size()>0)
    {
      // getting the groups of vector, each vector contains the knickpoints nodes to group
      vector<vector<int> > grouped_kp = group_local_kp(vecnode_kp,vecnode_river,Flowinfo, kp_node_search);
      // We have the group of vector now lets run through it to get the requested values
      for(vector<vector<int> >::iterator vlad = grouped_kp.begin(); vlad != grouped_kp.end(); vlad ++)
      {
        // getting the vector of node we need: the nodes containing 
        vector<int> this_vecnode = *vlad;
        // Easy case: non composite knickpoint, let's just record the same info thatn the raw detection
        if(this_vecnode.size() == 1)
        {
          int this_node = this_vecnode[0];// node of the knickpoint
          kp_segdrop[this_node] = raw_delta_segelev_from_TVDb_chi[this_node];
          sharpness_stepped_length[this_node] = 0; // sharpness = 0 as the knickpoint is a point
          stepped_extent[this_node] = make_pair(this_node, this_node); // the extent nodes are the same
          // Getting the coordinate of the centroid of the kp -> just regular x
          float this_x = 0,this_y = 0;
          Flowinfo.get_x_and_y_from_current_node(this_node, this_x, this_y);
          stepped_centroid[this_node] = make_pair(this_x,this_y);
          
          nearest_node_centroid_kp_stepped[this_node] = this_node;
          flow_distance_stepped_kp_centroid_map[this_node] = flow_distance_data_map[this_node];
        }
        else if(this_vecnode.size() > 1)
        {
          // Second case: several knickpoints to group, let's go step by step
          // the identifying node of the kp is the first one
          int this_node = this_vecnode[0];
          
          // summing the segelev for this group of node
          kp_segdrop[this_node] = get_dseg_drop_from_composite_kp(this_vecnode); // gobal segmentation elevation value of the knickpoint

          
          // Sharpness is the width of the combined knickpoints
          sharpness_stepped_length[this_node] = get_kp_sharpness_length(this_vecnode, Flowinfo);
          stepped_extent[this_node] = make_pair(this_vecnode[0], this_vecnode.back()); // the extent nodes are the extreme of the DD


          // The neirest node from the centre of the knickpoint (in regards to chi distance)
          
          // TO SORT, NOT ADAPTED TO THAT  -> TODO NEXT WEEK BORIS
          int nearnode = get_ksn_centroid_coordinates(Flowinfo, this_vecnode, vecnode_river,kp_segdrop[this_node]); // get the weighted x and y of the centroid.

          flow_distance_stepped_kp_centroid_map[this_node] = flow_distance_data_map[nearnode];
          // getting the 
          float this_x = 0,this_y = 0;
          Flowinfo.get_x_and_y_from_current_node(nearnode, this_x, this_y);

          stepped_centroid[this_node] = make_pair(this_x,this_y);
          nearest_node_centroid_kp_stepped[this_node] = nearnode;



        }
      }
    }
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  Group the adjacent local knickpoints from two vectors =
//                  Old version                           =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

vector<vector<int> > LSDChiTools::old_group_local_kp(vector<int> vecnode_kp, vector<int> vecnode_river,LSDFlowInfo& Flowinfo)
{
  // creating a map of neightboors for each nodes of the river 
  size_t it;
  map<int,pair<int,int> > neightboors;
  for(it = 0; it< vecnode_river.size(); it ++)
  {
    int this_node = vecnode_river[it];
    pair<int,int> this_pair;
    if(it == 0)
    {
      this_pair = make_pair(-9999, vecnode_river[it+1]);
    }
    else if(it == (vecnode_river.size() - 1))
    {
      this_pair = make_pair(vecnode_river[it-1],-9999);
    }
    else
    {
      this_pair = make_pair(vecnode_river[it-1],vecnode_river[it+1]);
    }
    neightboors[this_node] = this_pair;
  }
  // Done

  // Now looping through each nodes containing knickpoints to check if they have a direct neighboor
  vector<vector<int> > out_vector;
  vector<int> this_vecnode;
  if(vecnode_kp.size()>0)
  {
    // cout << "DEBUG test 1" << endl;
    // dealing with the first node
    this_vecnode.push_back(vecnode_kp[0]);
    if(vecnode_kp[0] != neightboors[0].second )
    {
      out_vector.push_back(this_vecnode);
      this_vecnode.clear();
    }
    // cout << "DEBUG test 2" << endl;

    // other nodes

    for(it = 1; it < vecnode_kp.size(); it++)
    {
      // which one is our working node
      int this_node = vecnode_kp[it], next_node = vecnode_kp[it+1], last_node = vecnode_kp[it-1] ;
      this_vecnode.push_back(this_node);
      if(it < (vecnode_kp.size()-1))
      {
        if( next_node != neightboors[this_node].second && last_node != neightboors[this_node].first)
        {
          out_vector.push_back(this_vecnode);
          this_vecnode.clear();
        }   
      }
      else 
      {
        out_vector.push_back(this_vecnode);
        this_vecnode.clear();
      }
    }

    // cout << "DEBUG test 3" << endl;

  }

  return out_vector;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  Group the adjacent local knickpoints from two vectors =
//                  New version                           =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

vector<vector<int> > LSDChiTools::group_local_kp(vector<int> vecnode_kp, vector<int> vecnode_river, LSDFlowInfo& Flowinfo, int kp_node_search)
{

  // pixel window to check on the knickpoints
  int HW = kp_node_search;
  // cout << "DEBUG_1" << endl;
  // getting the index of each knickpoint in the node vector
  size_t iced_t =0;
  vector<int> corresponding_index;
  for(; iced_t < vecnode_kp.size(); iced_t++)
  {
    int idx =0;
    bool found_it = false;
    for (size_t yh =0; found_it == false; yh++)
    {
      if(vecnode_kp[iced_t] == vecnode_river[yh])
      {
        idx = yh; // corresponding index
        found_it = true;
      }
    }
    corresponding_index.push_back(idx);
  }


  // Now I have the corresponding index
  // cout << "DEBUG_2 || " << corresponding_index.back() << endl;

  // Now creating a vector of number of node between this kp node and the following
  vector<int> n_node_to_next;

  for(size_t it = 0; it < corresponding_index.size()-1; it++)
  {
    n_node_to_next.push_back(corresponding_index[it+1] - corresponding_index[it]);
  }
  n_node_to_next.push_back(0);

  // cout << (int)vecnode_river.size() << " || " << (int)n_node_to_next.size() << " || " << (int)corresponding_index.size() << " || " << (int)vecnode_kp.size() << endl;

  // I got the number of node in between a knickpoint and the next
  // cout << "DEBUG_3" << endl;

  vector<vector<int> > out_vector;
  vector<int> this_vec;
  int counting_these_fucking_knickpoints = 0;

  for(size_t it = 0; it < vecnode_kp.size(); it++)
  {

    // saving the node
    bool save_the_raster = true;
    // cout << "DEBUG_3.1" << endl;
    int this_idx = corresponding_index[it];
    // cout << "DEBUG_3.2" << endl;
    this_vec.push_back(vecnode_kp[it]);
    counting_these_fucking_knickpoints++;
    // cout << "THIS KP NODE = " << vecnode_kp[it] << endl;

    if(n_node_to_next[it] <= HW)
    {
      // cout << "DEBUG_3.3" << endl;
      // cout << this_idx << " || " << n_node_to_next[it] << " || " << vecnode_river.size()-1 << endl;

      if((this_idx+n_node_to_next[it])> int(vecnode_river.size()-1))
      {
        cout << "GP_KP::FATAL ERROR, This is a weird error, try to rerun the analysis, if it persists contact B.G." << endl;
        exit(EXIT_FAILURE);
      }

      if((this_idx < int(vecnode_river.size()-1)) && (this_idx+n_node_to_next[it]<(int(vecnode_river.size()-1))) && it != vecnode_kp.size()-1)
      {
        // cout << n_node_to_next[it] << endl;
        // cout << "DEBUG_4" << endl;
        float this_kp = raw_ksn_kp_map[vecnode_river[this_idx]], next_kp = raw_ksn_kp_map[vecnode_river[this_idx+n_node_to_next[it]]] ;
        // cout << "DEBUG_5" << endl;
        // Check if they are both the same polarity
        // cout << this_idx << " || " << n_node_to_next[it] << " || " << vecnode_river.size()-1 << endl;
        if( (this_kp > 0 && next_kp > 0) || (this_kp <0 && next_kp<0) )

        {
          // cout << "DEBUG_5.5" << endl;
          save_the_raster = false;
        }
      }
    }
    // if not of these, I am saving this vector of node and clearing it
    if(save_the_raster == true)
    {
      // this should not happen, it helps me to debug
      if(this_vec.size() == 0)
      {
        cout <<"FATAL ERROR: void vector while combining (ERROR #45)" << endl;
        exit(EXIT_FAILURE);
      }
      out_vector.push_back(this_vec);
      this_vec.clear();
    }
  }

  // DEBUGGING to keep
  // int sumdfsdfa = 0;
  // for(vector<vector<int> >::iterator vlad = out_vector.begin(); vlad != out_vector.end(); vlad ++)
  // {
  //   vector<int> gyuyg = *vlad;
  //   sumdfsdfa += gyuyg.size();
  // }
  // cout<< "out: " << sumdfsdfa << " || in: " << vecnode_kp.size() << " || in_2: " << counting_these_fucking_knickpoints  << endl;  
  return out_vector;

}

/// Function to test TVD on elevation
/// This test is after submission to answer to WS comments
/// There is great potential in suck a method, let see how sensitive to lambda elevation is
/// The trickiest part would be to adpat lambda, I am afraid that the TVD might just try to flatten the thing
/// I'll try to comment this function on the go
/// Boris
void  LSDChiTools::TVD_on_segelev(LSDFlowInfo& Flowinfo)
{
  // Set the variables
  map<int,vector<int> >::iterator chirac;
  vector<int> this_vec;
  int this_river = 0;

  // Looping through the source key and getting the associated vector of river nodes
  for(chirac = map_node_source_key.begin(); chirac != map_node_source_key.end() ; chirac ++)
  {
    this_vec = chirac->second; // the vector of node
    vector<float> gros_test; // Debugging purposes
    // this next function Apply the TVD on the vector. It directly save the results in a global map and return general informations for debugging
    cout << "Denoising river #" << this_river << " with size " << this_vec.size() << endl;

    // Checking if the river isn't too small
    if(this_vec.size()>20){

      // I need a vector per data to TVD and for the corresponding resuslts
      vector<double> this_segelev_to_TVD, the_TVDed_vector_segelev, this_elev_to_TVD, the_TVDed_vector_elev, this_delev_to_TVD, the_TVDed_vector_delev, retrended_elevation, this_delev, this_chi, m_chi_z_from_retrended,curv_chi_from_retrended;
      // Pferd means Horse in German
      size_t pferd;
      // iterator through the node vector
      vector<int>::iterator vlad = this_vec.begin();

      // Vgetting the data to denoise
      for(; vlad != this_vec.end() ; vlad++)
      {
        int this_node = *vlad;
        this_segelev_to_TVD.push_back(double(segmented_elevation_map[this_node]));
        this_elev_to_TVD.push_back(elev_data_map[this_node]);
        this_delev_to_TVD.push_back(elev_data_map[this_node]); // storing the elevation before actually getting it 
        this_chi.push_back(chi_data_map[this_node]);

      }

      // detrening the elevation here: Just a first order derivation
      for(pferd = 0; pferd<this_vec.size()-1; pferd++)
      {
        this_delev_to_TVD[pferd] = this_elev_to_TVD[pferd] - this_elev_to_TVD[pferd + 1] ;
      }

      // This is just a Q&D way to initialize the vectors to the right format. Not proud of that but it's working.
      this_delev = this_delev_to_TVD;
      m_chi_z_from_retrended = this_delev_to_TVD;
      curv_chi_from_retrended = this_delev_to_TVD;

      // Dealing with the base node of each river: its derivation compare to its local receiver
      int receiver;
      Flowinfo.retrieve_receiver_information(this_vec[this_vec.size()-1],receiver);
      // And applying it
      this_delev_to_TVD[this_vec.size()-1] = elev_data_map[this_vec[this_vec.size()-1]] - elev_data_map[receiver];

      // Getting the base elevation of this particular river
      float base_elevation = elev_data_map[receiver]; // elevation at level -1

      // range of lambda to test
      static const double arr[] = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5,6,7,8,9,10,20};//,10,20,25,50,100,500};
      vector<double> lambdac(arr, arr + sizeof(arr) / sizeof(arr[0]));

      // OK now I am running the denoising on each vector and ...
      double this_lambda = 0;
      for(vector<double>::iterator jacques=lambdac.begin();jacques!=lambdac.end();jacques++)
      {

        this_lambda = *jacques; // This simply is the lambda currently tested
        the_TVDed_vector_segelev = TV1D_denoise_v2(this_segelev_to_TVD, this_lambda);
        the_TVDed_vector_elev = TV1D_denoise_v2(this_elev_to_TVD, this_lambda);
        the_TVDed_vector_delev = TV1D_denoise_v2(this_delev_to_TVD, this_lambda);

        // ... retrending the elevation
        retrended_elevation.clear(); // safety first but probably useless
        retrended_elevation = the_TVDed_vector_elev;
        double last_retrended_elev = 0;
        for (size_t jaguar=1;jaguar<=this_vec.size();jaguar++)
        {
          // inverting the iterator
          size_t djag = this_vec.size() - jaguar;
          // if this is the base I am detrending in regardfs to its receiver
          if(djag == this_vec.size()-1)
          {
            retrended_elevation[djag] = base_elevation + the_TVDed_vector_delev[djag];
          }
          else
          {
            // otherwise I detrend from my last detrending
            retrended_elevation[djag] = last_retrended_elev + the_TVDed_vector_delev[djag];
          }
          // Done, saving last step
          last_retrended_elev = retrended_elevation[djag];
        }



        // Now I just have to save the results into global maps
        int last_node = -9999; // forcing segfault in case it happens 
        for(pferd=0; pferd<this_vec.size(); pferd++)
        {
          int this_node = this_vec[pferd];

          pair<double,int> tkey(this_lambda,this_node);

          pair<double,int> last_key;

          TVD_segelev[tkey] = the_TVDed_vector_segelev[pferd];
          TVD_elev[tkey] = the_TVDed_vector_elev[pferd];
          TVD_delev[tkey] = the_TVDed_vector_delev[pferd];
          delev[tkey] = this_delev[pferd];
          TVD_retrend[tkey] = retrended_elevation[pferd];

          last_key = tkey;
        }
        // cout << "GSDJKFSDFSDF" << endl;
        // Calculating here the first derivative from chi-retrended profile (-> m_chi)
        for (size_t jaguar=1;jaguar<=this_vec.size();jaguar++)
        {
          // inverting the iterator
          size_t djag = this_vec.size() - jaguar;
          // if this is the base I am deriving in regardfs to its receiver
          if(djag == this_vec.size()-1)
          {
            pair<double,int> tpppp(this_lambda,receiver);
            // cout << djag << endl;
            // cout << TVD_retrend[tpppp] << endl;
            // cout << this_chi[djag] << endl;
            // cout << chi_data_map[receiver] << endl;
            double this_elev_to_use = 0;
            if(TVD_retrend.count(tpppp) == 0)
            {
              this_elev_to_use = elev_data_map[receiver];
            }
            else
            {
              this_elev_to_use = TVD_retrend[tpppp];
            }

            m_chi_z_from_retrended[djag] = abs((retrended_elevation[djag] - this_elev_to_use));///(this_chi[djag]-chi_data_map[receiver]));
          }
          else
          {
            // otherwise I detrend from my last detrending
            m_chi_z_from_retrended[djag] = abs((retrended_elevation[djag] - retrended_elevation[djag+1]));///(this_chi[djag]-this_chi[djag+1]));
          }
          // Done, saving last step
        }

        //Saving the map
        for(pferd=0; pferd<this_vec.size(); pferd++)
        {
          int this_node = this_vec[pferd];

          pair<double,int> tkey(this_lambda,this_node);

          pair<double,int> last_key;

          globmap_mchi_z_from_retrend[tkey] = m_chi_z_from_retrended[pferd];

          last_key = tkey;
        }

        // Calculating Now the curvature from mchi-from-retrended profile (-> C_chi)
        for (size_t jaguar=1;jaguar<=this_vec.size();jaguar++)
        {
          // inverting the iterator
          size_t djag = this_vec.size() - jaguar;
          // if this is the base I am deriving in regardfs to its receiver
          if(djag == this_vec.size()-1)
          {
            pair<double,int> tpppp(this_lambda,receiver);
            double this_mchi_to_use = 0;

            this_mchi_to_use = globmap_mchi_z_from_retrend[tpppp];

            curv_chi_from_retrended[djag] = (m_chi_z_from_retrended[djag] - this_mchi_to_use);///(this_chi[djag]-chi_data_map[receiver]);
          }
          else
          {
            // otherwise I detrend from my last detrending
            // cout << djag << endl;
            curv_chi_from_retrended[djag] = (m_chi_z_from_retrended[djag] - m_chi_z_from_retrended[djag+1]);///(this_chi[djag]-this_chi[djag+1]);
          }
          // Done, saving last step
        }

        //Saving the map
        for(pferd=0; pferd<this_vec.size(); pferd++)
        {
          int this_node = this_vec[pferd];

          pair<double,int> tkey(this_lambda,this_node);

          pair<double,int> last_key;

          globmap_cchi_from_retrend[tkey] = curv_chi_from_retrended[pferd];

          last_key = tkey;
        }



      }

    }

  // Next river
  this_river++ ;

  }

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  Apply a Total Variation Denoising filter on the data  =
//    Coded in LSDStatTools adapted from Condat L.2013    =
//            DOI: 10.1109/LSP.2013.2278339               =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

void  LSDChiTools::TVD_on_my_ksn( float lambda)
{
  // Set the variables
  map<int,vector<int> >::iterator chirac;
  vector<int> this_vec;
  int this_river = 0;
  // int max_node = 500;
  // Looping through the source key and getting the associated vector of river nodes
  for(chirac = map_node_source_key.begin(); chirac != map_node_source_key.end() ; chirac ++)
  {
    // cout << "Denoising river #" << this_river << endl;
    // int this_SK = chirac->first; // we don't need it, but maybe at some points
    this_vec = chirac->second; // the vector of node
    vector<float> gros_test; // Debugging purposes
    // this next function Apply the TVD on the vector. It directly save the results in a global map and return general informations for debugging
    cout << "Denoising river #" << this_river << "with size " << this_vec.size() << "\r";

    if(this_vec.size()>20){

      gros_test = TVD_this_vec(this_vec, lambda);
    }
  this_river++ ;

  }

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  Apply a Total Variation Denoising filter on a specific vector =
//    Coded in LSDStatTools adapted from Condat L.2013            =
//            DOI: 10.1109/LSP.2013.2278339                       =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float>  LSDChiTools::TVD_this_vec(vector<int> this_vec, float lambda)
{

  // Creating the containers I will need for the denoising
  vector<double> this_val_mchi, this_val_bchi, this_val_segelev; //TVD_segelev_diff
  // this iterator will iterate through the vector of node you want to denoise
  vector<int>::iterator chirac = this_vec.begin();

  // Loop through the data to gather the vector we want to TVD
  for( ; chirac != this_vec.end() ; chirac++)
  {
    int this_node = *chirac;
    if(chi_data_map[this_node] > 0){ //  Checking if there are no data


      this_val_mchi.push_back((double)M_chi_data_map[this_node]);
      // this_val_bchi.push_back((double)b_chi_data_map[this_node]);
      // this_val_segelev.push_back((double)segelev_diff[this_node]);
      // NOTE: We recast everythin to double, floating points somehow generate bugs, in rare cases.
    }


  }
  // Calling the actual denoising coded in Stat tools
  double clambda = lambda;
  vector<float> outtemp;
  if(this_val_mchi.size()>0)
  {
    vector<double> this_val_mchi_TVDed = TV1D_denoise_v2(this_val_mchi, clambda);
  
    // Our data is Denoised, yaaaay. Let's implement the global map to save it.
    for(size_t plo = 0; plo < this_vec.size() ; plo++ )
    {
      int this_node = this_vec[plo];
      TVD_m_chi_map[this_node] = (float)this_val_mchi_TVDed[plo];
    }
  
    // Formatting a debugging vector that I sometime use. Ignore that.
    vector<float> outtemp(this_val_mchi_TVDed.begin(),this_val_mchi_TVDed.end());

  }

  return outtemp;
  // Done, EZPZ
}

/// Test function to adapat the denoising size but not working at the moment.
void  LSDChiTools::TVD_this_vec_v2(vector<int> this_vec, float lambda, int max_node, string type)
{

  // Creating the containers I will need for the denoising
  vector<double> this_val_mchi, this_val_bchi, this_val_segelev; //TVD_segelev_diff
  // this iterator will iterate through the vector of node you want to denoise
  vector<int>::iterator chirac = this_vec.begin();


  for( ; chirac != this_vec.end() ; chirac++)
  {
    int this_node = *chirac;
    if(chi_data_map[this_node] > 0){
      this_val_mchi.push_back((double)M_chi_data_map[this_node]);
      this_val_bchi.push_back((double)b_chi_data_map[this_node]);
      this_val_segelev.push_back((double)segelev_diff[this_node]);
    }

  }

  // Calling the actual denoising coded in Stat tools
  double clambda = lambda;
  vector<double> this_val_mchi_TVDed = TV1D_denoise_v2(this_val_mchi, clambda);
  // vector<double> this_val_bchi_TVDed = TV1D_denoise_v2(this_val_bchi, dlamda);
  vector<double> this_val_segelev_TVDed = TV1D_denoise_v2(this_val_segelev, 5);


  // vector<double> this_val_TVDed_Corrected = correct_TVD_vec(this_val_TVDed);
  vector<int> vecnode_to_save;
  int beginning_index;
  if(type == "begin")
  {
    int slider = int(0.6*max_node);
    // vector<int>::const_iterator cons = this_vec.begin() + slider;
    vector<int> temp_vecta(this_vec.begin(), this_vec.begin()+ slider );
    vecnode_to_save = temp_vecta;
    beginning_index = 0;
  }
  else if(type == "middle")
  {
    int beg = int(0.2*max_node);
    vector<int> temp_vecta(this_vec.begin()+beg,this_vec.end()-beg);
    vecnode_to_save = temp_vecta;
    beginning_index = int(0.2*max_node);

  }
  else
  {
    int slider = int(0.4 * max_node);
    vector<int> temp_vecta(this_vec.begin()+ slider,this_vec.end());
    vecnode_to_save = temp_vecta;
    beginning_index = int(0.4 * max_node);

  }

  for(vector<int>::iterator plo = vecnode_to_save.begin(); plo != vecnode_to_save.end() ; plo++ )
  {
    // to switch to intermediate after test
    int this_node = *plo;
    TVD_m_chi_map[this_node] = (float)this_val_mchi_TVDed[beginning_index];
    // TVD_b_chi_map[this_node] = (float)this_val_bchi_TVDed[beginning_index];
    // TVD_segelev_diff[this_node] = (float)this_val_segelev_TVDed[beginning_index];

    // TVD_m_chi_map_non_corrected[this_node] = (float)this_val_TVDed[plo];
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Correct some abberations in the TVD      =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> LSDChiTools::correct_TVD_vec(vector<double> this_val)
{
  // looping through the nodes, if the value of the nex node exactly equals the value of the previous node
  // then this is an artifact of the TVD process for some reasons

  size_t tuile = 1;
  vector<double> this_val_out;
  // first elements unchanged


  this_val_out.push_back(this_val[0]);
  size_t si = this_val.size();


  for(;tuile < this_val.size() - 1 ; tuile++)
  {
    // if(this_val[tuile-1] +2 < this_val[tuile]) {cout << this_val[tuile-1] << " || " << this_val[tuile] << " || " << this_val[tuile+1] << endl;}
    if(this_val[tuile-1] == this_val[tuile+1])
    {
      
      this_val_out.push_back(this_val[tuile-1]);
    }
    else
    {
      // cout << this_val[tuile-1] << " || " << this_val[tuile+1] << endl ;
      this_val_out.push_back(this_val[tuile]);
    }
  }

  // last elements unchanged as well

  this_val_out.push_back(this_val[si-1]);


  return this_val_out;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// lump the m_chi to detect outliers        =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

void LSDChiTools::lump_my_ksn(int n_nodlump)
{

  // Set the variables
  map<int,vector<int> >::iterator chirac;
  vector<int> this_vec;
  for(chirac = map_node_source_key.begin(); chirac != map_node_source_key.end() ; chirac ++)
  {
    this_vec = chirac->second;
    size_t testi = 2*n_nodlump;
    if(this_vec.size() > testi)
    {
      lump_this_vec(this_vec,n_nodlump);
    }
    else
    {
      cout << "ignoring lumping on source " << chirac->first << ": not enough nodes." << endl;
    }
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// lump m_chi for these specific nodes      =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

void LSDChiTools::lump_this_vec(vector<int> this_vec, int n_nodlump)
{
  // Set the variables
  float this_mean = 0;
  vector<int>::iterator tnode;
  vector<float> this_val;
  int this_node;

  for(tnode = this_vec.begin(); tnode != this_vec.end() ; tnode ++)
  {
    this_node = *tnode;
    this_val.push_back(M_chi_data_map[this_node]);
  }

  
  
  for(size_t op = 0 ; op < this_vec.size() ; op ++)
  {
    this_node = this_vec[op];
    if(int(op) < int(n_nodlump))
    {
      vector<float>::const_iterator beg = this_val.begin(), en = this_val.begin()+ n_nodlump + op;
      vector<float> tvec(beg ,en );
      this_mean = get_mean_ignore_ndv(tvec, NoDataValue);
    }
    else if(op < this_vec.size() - n_nodlump)
    {
      vector<float>::const_iterator beg = this_val.begin()+ op, en = this_val.begin()+ op + n_nodlump;
      vector<float> tvec(beg ,en );

      this_mean = get_mean_ignore_ndv(tvec, NoDataValue);

    }
    else
    {
      vector<float>::const_iterator beg = this_val.begin()+ op, en  = this_val.end();
      vector<float> tvec(beg ,en );

      this_mean = get_mean_ignore_ndv(tvec , NoDataValue);

    }

    lumped_m_chi_map[this_node] = this_mean;
  }


}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// derive the segmented elevation to get it ready for the stepped knickpints quantification =
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDChiTools::derive_the_segmented_elevation()
{

  // Set the variables
  map<int,vector<int> >::iterator chirac;
  vector<int> this_vec;
  vector<int>::iterator nonode;
  int last_node = 0, this_node = 0;
  
  // looping through each river
  for(chirac = map_node_source_key.begin(); chirac != map_node_source_key.end() ; chirac ++)
  {
    this_vec = chirac->second;

    // looping through each node in the river if the river has at least two nodes
    if(this_vec.size() > 1)
    {
      for(nonode = this_vec.begin(); nonode != this_vec.end(); nonode++)
      {
        this_node = *nonode;
        if(nonode != this_vec.begin())
        {
          // i derive if this is not the last node
          this_node = *nonode;
          segelev_diff[this_node] = (segmented_elevation_map[last_node] - segmented_elevation_map[this_node]); // ---> to add if we derive it to chi (chi_data_map[last_node] - chi_data_map[this_node]); 
          segelev_diff_second[this_node] = segelev_diff_second[last_node] - segelev_diff_second[this_node];
        }
        else
        {
          // the first derivative is 0
          segelev_diff[this_node] = 0;
          segelev_diff_second[this_node] = 0;
          

        }
        // incrementing the rest
        last_node = this_node;
      }
    }

  }

}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate the KDE over the knickpoint map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::ksn_kp_KDE()
{
  // first setting the main iterator
  map<int,vector<int> >::iterator jacques;

  // function variables
  int this_SK = 0;
  vector<int> vecnode;

  for(jacques = map_node_source_key_kp.begin(); jacques != map_node_source_key_kp.end(); jacques++)
  {
    this_SK = jacques->first;
    vecnode = jacques->second;
    if(vecnode.size()>0)
     {KDE_vec_node_mchi(vecnode,this_SK);}

  }

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate the KDE using the mchi value corresponding to a vector of node index
// BG
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDChiTools::KDE_vec_node_mchi(vector<int> vecnode, int SK)
{
  // setting the iterator
  vector<int>::iterator valachie = vecnode.begin();

  // first getting the corresponding vector of values
  vector<float> veksn, veKDE;
  int this_node = 0;
  for(;valachie != vecnode.end(); valachie++)
  {
    this_node = *valachie;
    veksn.push_back(raw_dksndchi_kp_map[this_node]);
  }

  // now getting the KDE corresponding vector
  pair<float,vector<float> > pagul = auto_KDE(veksn);

  // incrementing the bandwidth map
  KDE_bandwidth_per_source_key[SK] = pagul.first;

  // dealing with retrieving the KDE per nodes
  veKDE = pagul.second;
  float this_KDE = 0;
  for(size_t uip = 0; uip < veKDE.size(); uip++)
  {
    this_node = vecnode[uip];
    this_KDE = veKDE[uip];
    raw_KDE_kp_map[this_node] = this_KDE;
  }

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Automate the outlier detection - ATM I am testing a bunch ou method
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::ksn_knickpoint_outlier_automator(LSDFlowInfo& FlowInfo, float MZS_th)
{

  // Ok looping through river
  map<int,vector<int> >::iterator salazar;
  vector<int> vecnode,vecoutlier_MZS_dkdc;
  vector<float> vecval;
  for(salazar = map_node_source_key_kp.begin(); salazar!= map_node_source_key_kp.end(); salazar++)
  {

    vecnode = salazar->second;
    if(vecnode.size()>0)
    {
      vecval = get_value_from_map_and_node(vecnode,raw_dksndchi_kp_map);
      vecoutlier_MZS_dkdc = is_outlier_MZS(vecval, NoDataValue, MZS_th);

      for(size_t hi = 0; hi < vecnode.size(); hi++)
      {
        map_outlier_MZS_dksndchi[vecnode[hi]] = vecoutlier_MZS_dkdc[hi];
      }
    }

  }

  // now dealing with the values after combining the knickpoints
  vecval.clear();
  vecnode.clear();

  // getting all the final knickpoints
  for(map<int,float>::iterator gorilla = ksn_kp_map.begin(); gorilla != ksn_kp_map.end(); gorilla ++)
  {
    vecnode.push_back(gorilla->first);
    vecval.push_back(gorilla->second);
  }

  vector<int> vecoutlier_MZS_combined = is_outlier_MZS(vecval, NoDataValue, MZS_th);

  for(size_t hi = 0; hi < vecnode.size(); hi++)
  {
    map_outlier_MZS_combined[vecnode[hi]] = vecoutlier_MZS_combined[hi];
  }




}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// write a file with the source_key, basinkey and the associated bandwidth
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDChiTools::print_bandwidth_ksn_knickpoint(string filename)
{
  
    // open the data file
  ofstream  file_out;
  file_out.open(filename.c_str());
  file_out << "source_key,basin_key,length,chi" << endl;

  int this_source_key, this_basin_key, this_node = 0;
  float this_bandwidth = 0;
  map<int,vector<int> >::iterator OL;

  for(OL = map_node_source_key.begin(); OL !=  map_node_source_key.end() ; OL++)
  {
    vector<int> albert = OL->second;
    if(albert.size()>0)
    {
      this_source_key = OL->first;
      this_node = OL->second[0];
      this_basin_key = baselevel_keys_map[this_node];
      this_bandwidth = KDE_bandwidth_per_source_key[this_source_key];
      file_out << this_source_key << ","
               << this_basin_key << ","
               << map_flow_length_source_key[this_source_key] << ","
               << map_chi_length_source_key[this_source_key] << endl;
    }
  }
  file_out.close();

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// write a file with the raw ksn knickpoints
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDChiTools::print_raw_ksn_knickpoint(LSDFlowInfo& FlowInfo, string filename)
{
  // these are for extracting element-wise data from the channel profiles.

  int this_node,row,col;
  float this_kp;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "longitude,latitude,elevation,flow_distance,chi,drainage_area,delta_ksn,dksn/dchi,KDE,basin_key,out_MZS,Cx,source_key";

  chi_data_out << endl;

  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    map<int,float>::iterator iter;

    for (iter = raw_ksn_kp_map.begin(); iter != raw_ksn_kp_map.end(); iter++)
    {
        this_node = iter->first;
        this_kp = iter->second;
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);



        chi_data_out.precision(9);
        chi_data_out << latitude << ","
                     << longitude << ",";
        chi_data_out.precision(5);
        chi_data_out << elev_data_map[this_node] << ","
                     << flow_distance_data_map[this_node] << ","
                     << chi_data_map[this_node] << ","
                     << drainage_area_data_map[this_node] << ","
                     << this_kp << ","
                     << raw_dksndchi_kp_map[this_node] << ","
                     << raw_KDE_kp_map[this_node] << ","
                     << baselevel_keys_map[this_node]<< ","
                     << map_outlier_MZS_dksndchi[this_node] << ","
                     << Cx_TVD_map[this_node] << ","
                     << source_keys_map[this_node];

        chi_data_out << endl;
    }
  }

  chi_data_out.close();
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// write a file with the final knickpoint output
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDChiTools::print_final_ksn_knickpoint(LSDFlowInfo& FlowInfo, string filename)
{
  // these are for extracting element-wise data from the channel profiles.

  int this_node,row,col;
  float this_kp;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // Map to check if a node has already been processed or not when adding the stepped component later
  map<int,bool> is_done;
  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "node,X,Y,latitude,longitude,elevation,flow_distance,chi,drainage_area,delta_ksn,delta_segelev,sharpness,sign,out,Cx,basin_key,source_key";

  chi_data_out << endl;

  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    map<int,float>::iterator iter;

    for (iter = ksn_kp_map.begin(); iter != ksn_kp_map.end(); iter++)
    {
        this_node = iter->first;
        this_kp = iter->second;
        int nearnode = nearest_node_centroid_kp[this_node];
        // int nearnode_stepped = nearest_node_centroid_kp_stepped[this_node];    // not used (SMM)
        float this_segelev = 0;

        if(kp_segdrop.count(nearnode) == 1)
        {
          is_done[nearnode] = true;
          this_segelev = kp_segdrop[nearnode];
        }
        // get the centroid location
        // float this_x = ksn_centroid[this_node].first, this_y = ksn_centroid[this_node].second; // BUGGY AT THE MOMENT
        float this_x = 0, this_y =0;
        FlowInfo.retrieve_current_row_and_col(nearnode,row,col);
        get_x_and_y_locations(row, col, this_x, this_y);

        get_lat_and_long_locations_from_coordinate(this_x, this_y, latitude, longitude, Converter);

        // Just adding sign column for plotting purposes
        int this_sign = 0;
        if(this_kp>0){this_sign = 1;}else{this_sign=(-1);}

        chi_data_out << this_node << ",";
        chi_data_out << this_x << ",";
        chi_data_out << this_y << ",";
        chi_data_out.precision(9);
        chi_data_out << latitude << ","
                     << longitude << ",";
        chi_data_out.precision(5);
        // NOTE: The IDENTIFYING NODE on the knickpoint map is the first of a knickpoint group - On the global map it provide the nearest node to get the coordinates of the centroid
        chi_data_out << elev_data_map[nearnode] << "," // NOTE -> nearnode is the centroid node, however the knickpoint info are stored in this node. I am still working on the centroidisation of the node
                     << flow_distance_data_map[nearnode] << ","
                     << chi_data_map[nearnode] << "," // NOTE -> nearnode is the centroid node, however the knickpoint info are stored in this node. I am still working on the centroidisation of the node
                     << drainage_area_data_map[nearnode] << "," // NOTE -> nearnode is the centroid node, however the knickpoint info are stored in this node. I am still working on the centroidisation of the node
                     << this_kp << ","
                     << this_segelev << ","
                     << sharpness_ksn_length[this_node] << ","
                     << this_sign << ","
                     << map_outlier_MZS_combined[this_node] << ","
                     << Cx_TVD_map[this_node] << ","
                     << baselevel_keys_map[this_node]<< ","
                     << source_keys_map[this_node];

        chi_data_out << endl;
    }

    // I have printed the ksn knickpoints to the file, and the stepped knickpoints when located at the same location. I have recorded the one I've already written.
    // Now I'll complete the writing with the knickpoints that are only stepped

    for (iter = kp_segdrop.begin(); iter != kp_segdrop.end(); iter++)
    {
        this_node = iter->first;
        float this_segelev = iter->second;
        // int nearnode = nearest_node_centroid_kp_stepped[this_node];
        float this_kp = 0; // All the ksn knickpoints have been written, these one only have a segelev component

        // if(chi_data_map[this_node] == 0)
        // {
        //   cout << "This node is screwed" <<endl;
        // }
        // if(chi_data_map[nearnode] == 0)
        // {
        //   cout << "nearnode is screwed" <<endl ;
        // }

        if(is_done.count(this_node) != 1 && chi_data_map[this_node] != 0)
        {
          // cout << "Tbg 45b" << endl;
          // get the centroid location
          // float this_x = ksn_centroid[this_node].first, this_y = ksn_centroid[this_node].second; // BUGGY AT THE MOMENT
          float this_x = 0, this_y =0;
          FlowInfo.retrieve_current_row_and_col(this_node,row,col);
          get_x_and_y_locations(row, col, this_x, this_y);

          get_lat_and_long_locations_from_coordinate(this_x, this_y, latitude, longitude, Converter);

          // Just adding sign column for plotting purposes
          int this_sign = 1;
          chi_data_out << this_node << ",";
          chi_data_out << this_x << ",";
          chi_data_out << this_y << ",";
          chi_data_out.precision(9);
          chi_data_out << latitude << ","
                       << longitude << ",";
          chi_data_out.precision(5);
          // NOTE: The IDENTIFYING NODE on the knickpoint map is the first of a knickpoint group - On the global map it provide the nearest node to get the coordinates of the centroid
          chi_data_out << elev_data_map[this_node] << "," // NOTE -> this_node is the centroid node, however the knickpoint info are stored in this node. I am still working on the centroidisation of the node
                       << flow_distance_data_map[this_node] << ","
                       << chi_data_map[this_node] << "," // NOTE -> this_node is the centroid node, however the knickpoint info are stored in this node. I am still working on the centroidisation of the node
                       << drainage_area_data_map[this_node] << "," // NOTE -> this_node is the centroid node, however the knickpoint info are stored in this node. I am still working on the centroidisation of the node
                       << this_kp << ","
                       << this_segelev << ","
                       << sharpness_ksn_length[this_node] << ","
                       << this_sign << ","
                       << map_outlier_MZS_combined[this_node] << ","
                       << baselevel_keys_map[this_node]<< ","
                       << source_keys_map[this_node];

          chi_data_out << endl;
        }
    }

  }
  // done, let me close the file correctly
  chi_data_out.close();
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file for the knickpoint algorithm
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_mchisegmented_knickpoint_version_test_on_segmented_elevation(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row, col;
  double latitude,longitude;
  double this_x,this_y;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // test to see if there is segment numbering
  bool have_segments = false;
  if( segment_counter_map.size() == node_sequence.size())
  {
    have_segments = true;
  }

  // test to see if the fitted elevations have been calculated
  bool have_segmented_elevation = false;
  have_segmented_elevation = true;
  


  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "node,X,Y,row,col,latitude,longitude,chi,elevation,flow_distance,drainage_area,lambda,TVD_segelev,TVD_elev,TVD_delev,TVD_retrend,delev,mchi_from_retrend,c_chi_from_retrend,source_key,basin_key";


  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (map<pair<double,int>, double>::iterator gala = TVD_segelev.begin(); gala!=TVD_segelev.end();gala++)
    {
      pair<double,int> this_pair;
      this_pair = gala->first;
      this_node = this_pair.second;
      // pregetting the data
      double this_lambda, this_TVDELEVTEST, this_TVD_elev, this_TVD_delev, this_TVD_retrend, this_delev, this_mchi_from_retrend, this_C_chi;
      this_lambda = this_pair.first;
      this_TVDELEVTEST = gala->second;
      this_TVD_elev = TVD_elev[this_pair];
      this_TVD_delev = TVD_delev[this_pair];
      this_TVD_retrend = TVD_retrend[this_pair];
      this_delev = delev[this_pair];
      this_mchi_from_retrend = globmap_mchi_z_from_retrend[this_pair];
      this_C_chi = globmap_cchi_from_retrend[this_pair];

      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);
      get_x_and_y_locations(row, col, this_x, this_y);
      int ksnkp = 0, segelevkp = 0; // 0 = no knickpoint, -1 negative, 1 positive
      // checking if there is a raw knickpoint here
      if(raw_ksn_kp_map.count(this_node)!=0) 
      {
        if(raw_ksn_kp_map[this_node]>0)
        {
          ksnkp = 1;
        }
        else if(raw_ksn_kp_map[this_node]<0)
        {
          ksnkp = -1;
        }

        if(segelev_diff[this_node]>0)
        {
          segelevkp = 1;
        }
        if(segelev_diff[this_node]>0)
        {
          segelevkp = -1;
        }

      }

      chi_data_out << this_node << ","
                   << this_x << ","
                   << this_y << ","
                   << row << ","
                   << col << ",";

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << this_lambda << ","
                   << this_TVDELEVTEST << ","
                   << this_TVD_elev << ','
                   << this_TVD_delev << ','
                   << this_TVD_retrend << ','
                   << this_delev << ','
                   << this_mchi_from_retrend << ','
                   << this_C_chi << ','
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// Creating a new series of knickpoint extraction

void LSDChiTools::generate_knickpoint_overview(LSDFlowInfo& FlowInfo, LSDRaster& Elevation, LSDIndexRaster& FlowAcc , float movern, float tA0, string filename)
{

  // Experimental function to respond to WS
  cout << "Building knickpoint Overview metrics" << endl;
  // Firstly I am getting the vector of separated flows
  vector< vector<int> > vector_flow_route = FlowInfo.get_vectors_of_flow(Elevation);
  // ;// = FlowInfo.get_vectors_of_flow(Elevation);
  // cout << "DEBUG::Here " << endl;

  // Going through it with a verry relevant iterator name (No incidence on code comprehension so:
  Array2D<float> tdenoised_elevation(NRows,NCols,NoDataValue);
  Array2D<float> tkp_overview_chi_coordinates(NRows,NCols,NoDataValue);
  Array2D<float> tCx_TVD_map(NRows,NCols,NoDataValue);
  Array2D<float> tMx_TVD_map(NRows,NCols,NoDataValue);

  for(size_t h=0;h<NRows;h++)
  {
    for(size_t g=0;g<NCols;g++)
    {
      if(Elevation.get_data_element(h,g) != NoDataValue)
      {
        tdenoised_elevation[h][g] = 0;
        tkp_overview_chi_coordinates[h][g] = 0;
        tCx_TVD_map[h][g] = 0;
        tMx_TVD_map[h][g] = 0;
      }
    }

  }

  cout << "Processing the metrics: I have to process " << vector_flow_route.size() << " flow vectors." << endl;
  
  int cpt = 1;
  double lambda_TVD = 10;
  for(vector< vector<int> >::iterator Hund_und_Katze =  vector_flow_route.begin(); Hund_und_Katze != vector_flow_route.end(); Hund_und_Katze++)
  {
    // if(cpt % 10 == 0)
    // {
    //   cout << cpt << "/" << vector_flow_route.size() << " processed" << "\r";
    // }
    // cpt++;


    // this vector of connected node
    vector<int> nodevec = *Hund_und_Katze;
    vector<double> elevec;
    vector<int> Vrow,Vcol;
    vector<int> Vx,Vy;
    // getting the elevation

    for(size_t i=0; i<nodevec.size(); i++)
    {
      int row=0, col=0, this_node = nodevec[i];
      double x=0,y=0;
      // Getting the info
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      FlowInfo.get_x_and_y_locations(row,col,x,y);
      // Implementing the vectors
      elevec.push_back(Elevation.get_data_element(row, col));
      Vrow.push_back(row);
      Vcol.push_back(col);      
      Vx.push_back(x);
      Vy.push_back(y);
    }
    // Now I have the info I need

    // TVD on my elevation:

    // first need to detrend
    vector<double> this_detrend_elevec;
    for(size_t i=0; i<nodevec.size(); i++)
    {
      int this_node = nodevec[i];
      int row = Vrow[i];
      int col = Vcol[i];
      double this_detrend = 0;
      double last_elev = 0;
      if(i==0)
      {
        int receiver_id, receiver_row, receiver_col;
        FlowInfo.retrieve_receiver_information(this_node,receiver_id,receiver_row,receiver_col);
        last_elev = Elevation.get_data_element(receiver_row, receiver_col);
      }
      else
      {
        last_elev = elevec[i-1];
      }
      this_detrend_elevec.push_back(elevec[i] - last_elev);
      // cout << this_detrend_elevec.back() << endl;
      if(this_detrend_elevec.back()<0)
      {
        cout<< "Knickpoint overview failure: got contradictory elevations" << endl;
        exit(EXIT_FAILURE);
      }
    }  
    
    vector<double> this_denoised_delevec;
    if(this_detrend_elevec.size()>2)
    {
      this_denoised_delevec = TV1D_denoise_v2(this_detrend_elevec, lambda_TVD);
    }
    else
    {
      // No denoising if too small -> median of values
      double med = get_median(this_detrend_elevec);
      for(size_t yu=0; yu< this_detrend_elevec.size(); yu++)
      {
        this_denoised_delevec.push_back(med);
      }
    }

    // Now I need to retrend:
    vector<double> this_retrend_elevec = this_denoised_delevec;
    for(size_t i=0; i<nodevec.size(); i++)
    {
      int this_node = nodevec[i];
      int row = Vrow[i];
      int col = Vcol[i];
      double this_detrend = 0;
      double last_elev = 0;
      if(i==0)
      {
        int receiver_id, receiver_row, receiver_col;

        FlowInfo.retrieve_receiver_information(this_node,receiver_id,receiver_row,receiver_col);
        if(receiver_id == this_node || tdenoised_elevation[receiver_row][receiver_col] == NoDataValue)
        {
          // Base level node
          tdenoised_elevation[row][col] = 0;
          this_denoised_delevec[i] = 0;
          last_elev = 0;
        }
        else
        {

          last_elev = tdenoised_elevation[receiver_row][receiver_col];
        }
      }
      else
      {

        last_elev = this_retrend_elevec[i-1];
      }
      this_retrend_elevec[i] = this_detrend_elevec[i] + last_elev;
      tdenoised_elevation[row][col] = this_detrend_elevec[i] + last_elev;

    }  

    // Lets reiterate through that shit and implement the metrics we need
    // Thanks to the clever node ordering, I can calculate chi and stuffs confidently
    // if all of that takes time I'll break it into several vector of nodes from their stream order or other things
    // I also need two more additionnal vector!
    vector<double> delta_chi, mchi_vec, mchi_vec_tvd;
    for(size_t i=0; i<nodevec.size(); i++)
    {
      int this_node = nodevec[i];
      int row = Vrow[i];
      int col = Vcol[i];
      
      // receiver info
      vector<double> steepness;
      int receiver_id, receiver_row, receiver_col;
      FlowInfo.retrieve_receiver_information(this_node,receiver_id,receiver_row,receiver_col);
      if(receiver_id == this_node || tdenoised_elevation[receiver_row][receiver_col] == NoDataValue)
      {
        // I am a baselevel
        tkp_overview_chi_coordinates[row][col] = 0;
        delta_chi.push_back(0);
        mchi_vec.push_back(0);
      }
      else
      {
        double rx,ry;
        // getting the receiver/donors XY
        FlowInfo.get_x_and_y_locations(receiver_row,receiver_col,rx,ry);
        double dx = sqrt(pow(Vx[i]-rx, 2) + pow(Vy[i]-ry, 2));
        // Calculating Chi
        double tDA = FlowAcc.get_data_element(row,col);
        double last_DA = FlowAcc.get_data_element(receiver_row,receiver_col);
        double last_chi = tkp_overview_chi_coordinates[receiver_row][receiver_col];
        double tchi = last_chi + pow((((tA0/last_DA)+(tA0/tDA))/2)*dx,movern);
        tkp_overview_chi_coordinates[row][col] = tchi; // got it
        double dchi = tchi - last_chi;
        // double tdelev =
        mchi_vec.push_back((tdenoised_elevation[row][col] - tdenoised_elevation[receiver_row][receiver_col])/dchi); 
        delta_chi.push_back(dchi);
        
      }

      // TVD on my Mchi
      if(mchi_vec.size()>1) 
      {
        mchi_vec_tvd = TV1D_denoise_v2(mchi_vec, lambda_TVD);
      }
      else
      { 
        mchi_vec_tvd=mchi_vec;
      }

      // getting the rest
      for(size_t i=0; i<nodevec.size(); i++)
      {
        int this_node = nodevec[i];
        int row = Vrow[i];
        int col = Vcol[i];
        
        // receiver info
        vector<double> steepness;
        int receiver_id, receiver_row, receiver_col;
        FlowInfo.retrieve_receiver_information(this_node,receiver_id,receiver_row,receiver_col);
        if(receiver_id == this_node || tdenoised_elevation[receiver_row][receiver_col] == NoDataValue)
        {
          // I am a baselevel
          tMx_TVD_map[row][col] = 0;
          tCx_TVD_map[row][col] = 0;
        }
        else
        {
          // Need to get everything now
          tMx_TVD_map[row][col] = mchi_vec_tvd[i];
          tCx_TVD_map[row][col] = (tMx_TVD_map[row][col] - tMx_TVD_map[receiver_row][receiver_col])/delta_chi[i];  
        }
      }
    }



      
  }

  cout << endl;
  cout << "Done with calculation, now writing to file the different rasters ";
  // Done with all the workable nodes

  // Let's write the rasters now
  LSDRaster Rdenoised_elevation(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, tdenoised_elevation,GeoReferencingStrings);
  Rdenoised_elevation.write_raster(filename+"_denoised_dem","bil");
  cout << " Denoised elevation OK ";

  LSDRaster Rkp_overview_chi_coordinates(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, tkp_overview_chi_coordinates,GeoReferencingStrings);
  Rkp_overview_chi_coordinates.write_raster(filename+"_chi_glob","bil");
  cout << " Chi global raster OK ";

  LSDRaster RCx_TVD_map(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, tCx_TVD_map,GeoReferencingStrings);
  RCx_TVD_map.write_raster(filename+"_Cchi_TVD","bil");
  cout << " Curvature from TVD elevation OK ";

  LSDRaster RMx_TVD_map(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, tMx_TVD_map,GeoReferencingStrings);
  RMx_TVD_map.write_raster(filename+"_Mchi_TVD","bil");
  cout << " Chi-z gradient from TVD elevation OK " << endl;

        
     
}






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file for the knickpoint algorithm
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_mchisegmented_knickpoint_version(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row, col;
  double latitude,longitude;
  double this_x,this_y;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // test to see if there is segment numbering
  bool have_segments = false;
  if( segment_counter_map.size() == node_sequence.size())
  {
    have_segments = true;
  }

  // test to see if the fitted elevations have been calculated
  bool have_segmented_elevation = false;
  have_segmented_elevation = true;
  


  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "node,X,Y,row,col,latitude,longitude,chi,elevation,flow_distance,drainage_area,m_chi,lumped_ksn,TVD_ksn,TVD_segelev_diff,b_chi,ksnkp,segelevkp,source_key,basin_key";
  if(have_segmented_elevation)
  {
    chi_data_out << ",segmented_elevation,mean_segdiff,std_segdiff,segdiff";
  }
  if (have_segments)
  {
    chi_data_out << ",segment_number";
    cout << "I added the segment number in the csv file"<< endl;
  }
  chi_data_out << endl;




  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);
      get_x_and_y_locations(row, col, this_x, this_y);
      int ksnkp = 0, segelevkp = 0; // 0 = no knickpoint, -1 negative, 1 positive
      // checking if there is a raw knickpoint here
      if(raw_ksn_kp_map.count(this_node)!=0) 
      {
        if(raw_ksn_kp_map[this_node]>0)
        {
          ksnkp = 1;
        }
        else if(raw_ksn_kp_map[this_node]<0)
        {
          ksnkp = -1;
        }

        if(segelev_diff[this_node]>0)
        {
          segelevkp = 1;
        }
        if(segelev_diff[this_node]>0)
        {
          segelevkp = -1;
        }

      }

      chi_data_out << this_node << ","
                   << this_x << ","
                   << this_y << ","
                   << row << ","
                   << col << ",";

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << M_chi_data_map[this_node] << ","
                   << lumped_m_chi_map[this_node] << ","
                   << TVD_m_chi_map[this_node] << ","
                   << TVD_segelev_diff[this_node] << ","
                   << b_chi_data_map[this_node] << ","
                   << ksnkp << ","
                   << segelevkp << ","                  
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      if(have_segmented_elevation)
      {
        chi_data_out << "," << segmented_elevation_map[this_node];
        chi_data_out << "," << mean_for_kp[this_node];
        chi_data_out << "," << std_for_kp[this_node];
        chi_data_out << "," << segelev_diff[this_node];


      }
      if (have_segments)
      {
        chi_data_out << "," << segment_counter_map[this_node];
      }
      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


void LSDChiTools::set_map_of_source_and_node(LSDFlowInfo& FlowInfo, int n_node_downstream)
{
  // find the number of nodes
  // setting the initilal condition
  int n_nodes = (node_sequence.size()), last_SK = source_keys_map[node_sequence[0]], this_SK = source_keys_map[node_sequence[0]], this_node = node_sequence[0], temp_receiver_node = 0, last_node = 0;
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    // setting temporary vector that I may need
    vector<int> temp_node_SK;
    temp_node_SK.push_back(this_node);

    // Now looping through all the node, gathering the one in the same river
    for (int n = 0; n< n_nodes; n++)
    {
      // getting each nodes information
      this_node = node_sequence[n];
      this_SK = source_keys_map[this_node];
      if(this_SK == last_SK)
      {
        // If the source key is the same than the previous one ---> incrementing the vector of node for each river
        temp_node_SK.push_back(this_node);
      }
      else
      {
        // Investigating the knickpoints between the tributaries and the main stem. I am saving that in a different file.
        // if different source key: first getting the receiving node
        vector<int> vector_of_nodes_dowstream_for_this_source_key; 
        for(int i = 0; i<n_node_downstream; i++)
        {
          FlowInfo.retrieve_receiver_information(last_node,temp_receiver_node);
          // pushing it back
          if(temp_receiver_node != -9999 || temp_receiver_node != NoDataValue || temp_receiver_node != 0)
          {
            vector_of_nodes_dowstream_for_this_source_key.push_back(temp_receiver_node);
          }
        }
        // saving this source key
        map_node_source_key[last_SK] = temp_node_SK;
        map_source_key_vecnode_of_receiver[last_SK] = vector_of_nodes_dowstream_for_this_source_key;
        // clearing the vector for the next source key and saving the current node in the new river
        temp_node_SK.clear();
        temp_node_SK.push_back(this_node);
      }
      // saving the last node info for next loop
      last_SK = this_SK;
      last_node = this_node;
    }
  }
// Done
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Just get the length, number of nodes and other basic metrics about the river
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::compute_basic_matrics_per_source_keys(LSDFlowInfo& FlowInfo)
{

  // setting the variables
  map<int,vector<int> >::iterator valachie;
  int this_SK =0;
  float dist = 0 , chi_dist = 0;
  vector<int> vecval;

  // Looping through each source key and getting the length of river in meters and chi space.
  // note: the first and last node of the vector are the extremes of each rivers plus the receiver node.
  for(valachie = map_node_source_key.begin(); valachie != map_node_source_key.end(); valachie++)
  {
    this_SK = valachie->first;
    vecval = valachie->second;
    chi_dist = abs(chi_data_map[vecval[0]] - chi_data_map[vecval.back()]);
    dist = abs(flow_distance_data_map[vecval[0]] - flow_distance_data_map[vecval.back()]);
    map_flow_length_source_key[this_SK] = dist;
    map_chi_length_source_key[this_SK] = chi_dist;
  }
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// medium old function for the knickpoint detection
// save the difference and ratio between m_chi values of each segments
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::ksn_knickpoint_detection(LSDFlowInfo& FlowInfo)
{

  // setting the variables for extracting the knickpoints
  int this_node = 0; // Node to investigate
  map<int,float> this_kickpoint_diff_map; // map of the delta k_sn, key is node number from Flowinfo
  map<int,float> this_knickpoint_rad; // map of the radian angle, key is node number from Flowinfo
  map<int,float> this_kickpoint_ratio_map; // map of the ratio k_sn, key is node number from Flowinfo
  map<int,int> this_knickpoint_sign_map; // map of the sign k_sn, key is node number from Flowinfo
  float last_M_chi, this_M_chi, last_M_atan, this_M_atan; // Floating storage of each last/new knickpoint value
  float delta_mchi = 0; // difference between last and new m_chi
  float ratio_mchi = 0; // ratio between last and new m_chi
  float delta_atan = 0; // difference for the slope in radian
  int knickpoint_sign = 0; // sign of the knickpoint: + =1 and - = -1
  int last_node = 0; // store the last investigated node to investigate
  int number_of_0 = 0; // debug storage of the number of knickpoint ratio recasted to avoid x/0
  int n_knp = 0; // hum ... probably the number of knickpoint, I have no recognition of this
  float max_elev = 0; // used to recast the M_chi value for angle determination 
  float max_chi = 0; // used to recast the M_chi value for angle determination
  float natural_coeff = 0;  // used to recast the M_chi value for angle determination
  map<int,vector<int> > this_node_kp_per_source_key; // this map store the node of each river knickpoint, the key is the source_key
  map<int,float> this_cumul_ksn_map; // This map store the cumulative ksn for each rivers
  map<int,float> this_cumul_rksn_map; // This map store the cumulative ratio ksn for each rivers
  map<int,float> this_cumul_rad_map; // This map store the cumulative ksn for each rivers
  map<pair<int,int>, float> this_knickzone_raw_cumul_ksn_map; //this map store the raw cumulative value of each knickpoints
  map<pair<int,int>, float> this_knickzone_raw_cumul_rksn_map; //this map store the raw cumulative value of each knickpoints
  map<pair<int,int>, float> this_knickzone_raw_cumul_rad_map; //this map store the raw cumulative value of each knickpoints
  map<pair<int,int>, vector<int> > this_knickzone_list_of_nodes; // This map stores the nodes of each knickzones


  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    // Preprocessing step, I am getting the maximum of chi values, in order to recast the calculation of the angle.
    // short explanation: the angle is calculated using arctan to get a rad of the segment. ksn is obtained by calculating Mchi with A0 = 1
    // This usually give to Chi values between 0 and 20~30ish depending on the landscape. We keep this value for the ksn knickpoints calculation
    // However, to get a "natural" angle, we recast the Mchi to correspond to a Chi value comparable to the elevation, thus using a A0 to get a maximum chi similar to the maximum elevation.

    // First we want to get the maximum elevation and the maximum chi and the nodes per rivers
    set_map_of_source_and_node(FlowInfo, 20);
    for (int n = 0; n< n_nodes; n++)
    {
      if(elev_data_map[node_sequence[n]] > max_elev)
      {
        max_elev = elev_data_map[node_sequence[n]];
      }
      if(chi_data_map[node_sequence[n]] > max_chi)
      {
        max_chi = chi_data_map[node_sequence[n]];
      }
    }

    // natural coeff is then maxz/maxChi (= A0_new^(-m/n))
    natural_coeff = max_elev/max_chi;
    natural_coeff = 1/(natural_coeff);


    // At this point, you have the required condition to launch the knickpoint analysis
    cout << "I am now extracting the knickpoint dataset" << endl;

    // Initializing the first node
    this_node = node_sequence[0];
    last_M_chi =  M_chi_data_map[this_node];

    // Engaging the loop through the rivers nodes
    for (int n = 0; n< n_nodes; n++)
    {
      // set the nodes number and keep information about the previous one
      last_node = this_node;
      this_node = node_sequence[n];


      // Get the M_chi from the current node
      this_M_chi = M_chi_data_map[this_node];
      last_M_chi = M_chi_data_map[last_node];

      // recasting if negative M_Chi. These negative values are artifact for really flat segments, recasting it to 0 is then inconsequential
      if(this_M_chi < 0 && n>0){this_M_chi = 0;} // getting rid of the negative values because we don't want it, I don't want the n = 0 to avoid detecting fake knickpoint if the first value is actually negative

      // If the M_chi has changed I increment the knickpoints, I also check if the two point are on the same channel to avoid stange unrelated knickpoints
      if (this_M_chi != last_M_chi && source_keys_map[this_node] == source_keys_map[last_node])
      {
        // If this condition is satisfied, we change segment and the knickpoint will be saved
        // -> first thing to do is to save the node into the map of knickpoint per river
        this_node_kp_per_source_key[source_keys_map[this_node]].push_back(this_node);
       

        // Calculation of the arctan of the M_Chi to get angle of M_chi segment, We want the absolute value, atan can have some sign issue. Uses of atan2 solve this quadrant issue, however require the x/y value that would be painful to get here
        // We use the natural coeff to adjust M_chi
        last_M_atan = abs(atan(natural_coeff * last_M_chi));// if you want degrees *180/M_PI;
        this_M_atan = abs(atan(natural_coeff * this_M_chi));// if you want degrees *180/M_PI;

        // dealing with division/0 and calculation of the ratio
        if(this_M_chi == 0)
        {
          ratio_mchi = -9999; // correspond to +infinite
          number_of_0++;
        }
        else
        {
          ratio_mchi = last_M_chi/this_M_chi; // Ratio between last and new chi steepness
        }

        // calculation of the delta between the two segments referred as diff
        delta_mchi = (last_M_chi-this_M_chi); // diff between last and new chi steepness - note that there is no minus because we loop up to bottom

        // Determination of the sign
        if(delta_mchi<=0){knickpoint_sign = -1;} else {knickpoint_sign = 1;} // Assign the knickpoint sign value

        // Delta between the radian angle 
        delta_atan = (last_M_atan - this_M_atan); // note that the minus is because we loop up to bottom, in fact no, no minus, my bad

        // Allocate the values to local maps
        this_kickpoint_diff_map[this_node] = delta_mchi;
        this_kickpoint_ratio_map[this_node] = ratio_mchi;
        this_knickpoint_rad[this_node] = delta_atan;
        this_knickpoint_sign_map[this_node] = knickpoint_sign;
        n_knp ++;
        // end of the loop
      }
    }

    // Now calculating the cumulation of angle and ksn, we could incorporate that in the last loop, but it is clearer that way I think. The loss of efficienty isn't too bad

    // Iteration over each river, 
    map<int,vector<int> >::iterator marten; //I quite like the map<...>::iterator
    int begining_node = 0, ending_node = 0, last_node_tapir = 0; // temporaries integers
    vector<int> node_to_implement_for_the_knickzone_cumulation;
    vector<int> node_per_knickzones; 
    pair<int,int> temp_knickzone (0,0);
    for(marten = (this_node_kp_per_source_key.begin());marten != this_node_kp_per_source_key.end();marten++)
    {

      for(vector<int>::reverse_iterator tapir = marten->second.rbegin(); tapir!= marten->second.rend(); ++tapir) // The iterator can be reversed to loop backward !! I love iterators. 
      {
        // Now looping through the knickpoint of each river, *tapir is the pointer that refers to node number. /!\ Note that I am loooing using a reverse_iterator to go from the bottom to the top of each river
        
        // now calculating the cumulation for each methods
        // # First case, the node is the first of the river, we want to save the value as the first to be cumulated
        if(tapir == marten->second.rbegin())
        {
        // Getting the requested ksn values
        ksn_cumul_knickpoint_map[*tapir] = this_kickpoint_diff_map[*tapir]; // cumulating the ksn value
        // Getting the requested ratio ksn values
        rksn_cumul_knickpoint_map[*tapir] = this_kickpoint_ratio_map[*tapir]; // cumulating the rksn value
        // Getting the radian values
        rad_cumul_knickpoint_map[*tapir] = this_knickpoint_rad[*tapir]; // cumulating the rad value 

        node_to_implement_for_the_knickzone_cumulation.push_back(*tapir); // saving the node for completion
        begining_node = *tapir; //this will be a beginning node whatever happense

        }
        // if this is the last element of the river, we save everything and reinitialize for the following
        else if(tapir == marten->second.rend()-1)
        {
          if(this_knickpoint_sign_map[last_node_tapir] == this_knickpoint_sign_map[*tapir])
          {
            ksn_cumul_knickpoint_map[*tapir] += this_kickpoint_diff_map[*tapir];
            rksn_cumul_knickpoint_map[*tapir] += this_kickpoint_ratio_map[*tapir];
            rad_cumul_knickpoint_map[*tapir] += this_knickpoint_rad[*tapir];
            ending_node = *tapir;
            temp_knickzone = make_pair (begining_node,ending_node);

            node_to_implement_for_the_knickzone_cumulation.push_back(*tapir);
            this_knickzone_list_of_nodes[temp_knickzone] = node_to_implement_for_the_knickzone_cumulation; // saving the list of nodes per knickzones

          }
          else
          {
            ksn_cumul_knickpoint_map[*tapir] += this_kickpoint_diff_map[*tapir];
            rksn_cumul_knickpoint_map[*tapir] += this_kickpoint_ratio_map[*tapir];
            rad_cumul_knickpoint_map[*tapir] += this_knickpoint_rad[*tapir];
            temp_knickzone = make_pair (begining_node,ending_node);
            this_knickzone_list_of_nodes[temp_knickzone] = node_to_implement_for_the_knickzone_cumulation;

            node_to_implement_for_the_knickzone_cumulation.clear();
            begining_node = *tapir;
            ending_node = *tapir;
            temp_knickzone = make_pair (begining_node,ending_node);
            this_knickzone_list_of_nodes[temp_knickzone] = node_to_implement_for_the_knickzone_cumulation;

          }

          // reinitializing vector
          node_to_implement_for_the_knickzone_cumulation.clear();
        }
        // # Other Case, the signs are the same, so we just cumulate the values and save the node to change
        else if(this_knickpoint_sign_map[last_node_tapir] == this_knickpoint_sign_map[*tapir])
        {
          ksn_cumul_knickpoint_map[*tapir] += this_kickpoint_diff_map[*tapir]; // cumulating the ksn valu
          rksn_cumul_knickpoint_map[*tapir] += this_kickpoint_ratio_map[*tapir]; // cumulating the rksn valu
          rad_cumul_knickpoint_map[*tapir] += this_knickpoint_rad[*tapir]; // cumulating the rad value
          node_to_implement_for_the_knickzone_cumulation.push_back(*tapir); // saving the node for completion
        }
        // # Finally, the sign are differents so we save the cumulative value for each nodes of the raw knickzone then reinitialize everything for the following
        else
        {

          // This also end a knickzone
          ending_node = last_node_tapir;
          temp_knickzone = make_pair (begining_node,ending_node);
          this_knickzone_list_of_nodes[temp_knickzone] = node_to_implement_for_the_knickzone_cumulation; // saving the list of nodes per knickzones

          begining_node = *tapir;
          // reinitializing vector
          node_to_implement_for_the_knickzone_cumulation.clear();
          node_to_implement_for_the_knickzone_cumulation.push_back(*tapir);
        } 

        // Here lies the remianing of my beloved derivation code that diseapeared after realizing how useless it was

        // storing the last chi value for the derivative
        last_node_tapir = *tapir;
      }
    }
  }

  // print everything in the public/protected maps -> saving the calculated data in the system
  ksn_ratio_knickpoint_map = this_kickpoint_ratio_map;
  ksn_diff_knickpoint_map = this_kickpoint_diff_map;
  ksn_sign_knickpoint_map = this_knickpoint_sign_map;
  ksn_rad_knickpoint_map = this_knickpoint_rad;
  ksn_cumul_knickpoint_map = this_cumul_ksn_map;
  rksn_cumul_knickpoint_map = this_cumul_rksn_map;
  rad_cumul_knickpoint_map = this_cumul_rad_map;

  knickzone_weighting_completion(this_knickzone_list_of_nodes);

  //  cout << "I finished to detect the knickpoints, you have " << n_knp << " knickpoints, thus " << number_of_0 << " ratios are switched to -9999 due to 0 divisions." << endl;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Loop through all the knickzones to weight all the different combinations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::knickzone_weighting_completion(map<pair<int,int>, vector<int> > mapofnode)
{
  // Good resolution, I'll try to better comment my codes now then
  // This function attend to weight and save all the combinations of knickzones previously calculated
  // The idea is to apply an outlier to each possibles knickzones to select the most probable ones
  // I'll begin by declaring all the variables I'll need for this purpose
  int begining_node = 0, ending_node = 0, this_knickzone_ID = 0;
  vector<int> working_nodes; //  temporary vector to store the working nodes
  float chi_size = 0, weighted_sum_ksn = 0, weighted_sum_rksn = 0, weighted_sum_rad = 0, weighter_coeff = 1, weighter = 0, ksn_sum = 0, rksn_sum = 0, rad_sum = 0 ; 
  pair<int,int> temp_pair;


  // Message to check
  cout<< "I am now testing and writing all the possible knickzones combinations, weighted by their chi lenght. The execution time depends on the number of detected knickzones and can take a while." << endl;

  // then I'll loop throught the map:
  // first through the list of knickzones
  map<pair<int,int>, vector<int> >::iterator Dwarf_epauletted_fruit_bat; // iterator is a tool to loop through map, efficiently and cleanly, with a random animal name because it is monday morning and life need random animal names on a monday morning
  for(Dwarf_epauletted_fruit_bat = mapofnode.begin(); Dwarf_epauletted_fruit_bat!=mapofnode.end();Dwarf_epauletted_fruit_bat++)
  {

    working_nodes = Dwarf_epauletted_fruit_bat->second;
    //Ok, now I am looping through the nodes of each knickzones
    for(size_t ity = 0; ity!= working_nodes.size();ity ++)
    {
      //I am letting this typical debug statement to check if I am in a ascending mode or not. It can be consufing depending How I looped before 
      // cout << elev_data_map[Dwarf_epauletted_fruit_bat->second[ity]] << "||" << ity << endl;
      // I am in an ascending node mode

      // Just A quick note why I am note using iterator for this vector iteration, I find iterator really useful but slightly less clear when we want to use previous or next element in a vector or array
      // Note that I don't know yet If I will use that but anyway let's code

      // The difficulty here is too loop through all the knickzone combinations:
      // We need to loop through all the i to n combinations with 0 <= i <= n 

      // Each knickzone will be store in the global knickzone_WP_ksn maps

      begining_node = working_nodes[ity];
      ending_node = working_nodes[working_nodes.size()-1];

      // Checking if the knickzone does have 
      if(begining_node != ending_node)
      {
        chi_size = chi_data_map[ending_node]-chi_data_map[begining_node]; 
        weighter_coeff = 1;
        // second iteration through the knickzone
        for(size_t frutbat = ity; frutbat != working_nodes.size(); frutbat ++)
        {
          ending_node = working_nodes[frutbat];
          weighter = (exp((-(chi_data_map[ending_node]-chi_data_map[begining_node])) / (chi_size * weighter_coeff) ) );
          weighted_sum_ksn += (ksn_diff_knickpoint_map[ending_node] * weighter);
          weighted_sum_rksn += (ksn_ratio_knickpoint_map[ending_node] * weighter);
          weighted_sum_rad += (ksn_rad_knickpoint_map[ending_node] * weighter);
          ksn_sum += ksn_diff_knickpoint_map[ending_node];
          rksn_sum += ksn_ratio_knickpoint_map[ending_node];
          rad_sum += ksn_rad_knickpoint_map[ending_node];

          temp_pair = make_pair (begining_node,ending_node);
          knickzone_WP_ksn[temp_pair] = weighted_sum_ksn;
          knickzone_WP_rksn[temp_pair] = weighted_sum_rksn;
          knickzone_WP_rad[temp_pair] = weighted_sum_rad;
          knickzone_raw_cumul_ksn[temp_pair] = ksn_sum;
          knickzone_raw_cumul_rksn[temp_pair] = rksn_sum;
          knickzone_raw_cumul_rad[temp_pair] = rad_sum;
          knickzone_ID[temp_pair] = this_knickzone_ID;

        }

      }
      else
      {
        temp_pair = make_pair (begining_node,ending_node);
        knickzone_WP_ksn[temp_pair] = ksn_diff_knickpoint_map[ending_node];
        knickzone_WP_rksn[temp_pair] = ksn_ratio_knickpoint_map[ending_node];
        knickzone_WP_rad[temp_pair] = ksn_rad_knickpoint_map[ending_node];
        knickzone_raw_cumul_ksn[temp_pair] = ksn_diff_knickpoint_map[ending_node];
        knickzone_raw_cumul_rksn[temp_pair] = ksn_ratio_knickpoint_map[ending_node];
        knickzone_raw_cumul_rad[temp_pair] = ksn_rad_knickpoint_map[ending_node];
        knickzone_ID[temp_pair] = this_knickzone_ID;
      }

      weighted_sum_rad = 0;
      weighted_sum_ksn = 0;
      weighted_sum_rksn = 0;
      ksn_sum = 0;
      rksn_sum = 0;
      rad_sum = 0;
      

    }
  this_knickzone_ID ++;
  }
  cout << "I am done testing all your knickzones combinations" << endl;
} 

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file - knickpoint version
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_knickzone_to_csv(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  cout << "I am now writing your ksn knickpoint file:" << endl;
  int A_node, B_node, row,col;
  double Alatitude,Alongitude,Blatitude,Blongitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "Alatitude,Alongitude,Blatitude,Blongitude,Aelevation,Belevation,Aflow_distance,Bflow_distance,Achi,Bchi,Adrainage_area,Bdrainage_area,ksn,rksn,sign,rad,Wgksn,Wgrksn,Wgrad,source_key,basin_key,knickzone_key,lenght";

  chi_data_out << endl;

  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    map<pair<int,int>,float>::iterator iter;

    for (iter = knickzone_raw_cumul_ksn.begin(); iter != knickzone_raw_cumul_ksn.end(); iter++)
    {
        A_node = iter->first.first;
        B_node = iter->first.second;
        FlowInfo.retrieve_current_row_and_col(A_node,row,col);
        get_lat_and_long_locations(row, col, Alatitude, Alongitude, Converter);
        FlowInfo.retrieve_current_row_and_col(B_node,row,col);
        get_lat_and_long_locations(row, col, Blatitude, Blongitude, Converter);
        // cout << "printing node " << this_node << " with diff " << ksn_diff_knickpoint_map[this_node] << endl; 
        
        chi_data_out.precision(9);
        chi_data_out << Alatitude << ","
                     << Alongitude << ","
                     << Blatitude << ","
                     << Blongitude << ",";
        chi_data_out.precision(5);
        chi_data_out << elev_data_map[A_node] << ","
                     << elev_data_map[B_node] << ","
                     << flow_distance_data_map[A_node] << ","
                     << flow_distance_data_map[B_node] << ","
                     << chi_data_map[A_node] << ","
                     << chi_data_map[B_node] << ","
                     << drainage_area_data_map[A_node] << ","
                     << drainage_area_data_map[B_node] << ","
                     << knickzone_raw_cumul_ksn[iter->first] << ","
                     << knickzone_raw_cumul_rksn[iter->first] << ","
                     << ksn_sign_knickpoint_map[A_node] << ","
                     << knickzone_raw_cumul_rad[iter->first] << ","
                     << knickzone_WP_ksn[iter->first] << ","
                     << knickzone_WP_rksn[iter->first] << ","
                     << knickzone_WP_rad[iter->first] << ","
                     << source_keys_map[A_node] << ","
                     << baselevel_keys_map[A_node]<< ","
                     << knickzone_ID[iter->first] << ","
                     << (chi_data_map[B_node]-chi_data_map[A_node]);
        chi_data_out << endl;
    }
  }

  chi_data_out.close();
  cout << "I am done, your file is:" << endl;
  cout << filename << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file - knickpoint version
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_knickpoint_to_csv(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  cout << "I am now writing your ksn knickzone file:" << endl;
  int this_node, row,col;
  double latitude,longitude, this_x, this_y;;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "Y,X,latitude,longitude,elevation,flow_distance,chi,drainage_area,ksn,rksn,sign,rad,cumul_ksn,cumul_rksn,cumul_rad,source_key,basin_key";

  chi_data_out << endl;




  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    map<int,float>::iterator iter;

    for (iter = ksn_diff_knickpoint_map.begin(); iter != ksn_diff_knickpoint_map.end(); iter++)
    {
        this_node = iter->first;
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        get_x_and_y_locations(row, col, this_x, this_y);
        
        chi_data_out.precision(9);
        chi_data_out << this_y << ","
                     << this_x << ","
                     << latitude << ","
                     << longitude << ",";
        chi_data_out.precision(5);
        chi_data_out << elev_data_map[this_node] << ","
                     << flow_distance_data_map[this_node] << ","
                     << chi_data_map[this_node] << ","
                     << drainage_area_data_map[this_node] << ","
                     << ksn_diff_knickpoint_map[this_node] << ","
                     << ksn_ratio_knickpoint_map[this_node] << ","
                     << ksn_sign_knickpoint_map[this_node] << ","
                     << ksn_rad_knickpoint_map[this_node] << ","
                     << ksn_cumul_knickpoint_map[this_node] << ","
                     << rksn_cumul_knickpoint_map[this_node] << ","
                     << rad_cumul_knickpoint_map[this_node] << ","
                     << source_keys_map[this_node] << ","
                     << baselevel_keys_map[this_node];
        chi_data_out << endl;
    }
  }

  chi_data_out.close();
  cout << "I am done, your file is:" << endl;
  cout << filename << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file - knickpoint version
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_segmented_gradient_to_csv(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  cout << "I am now writing your segmented gradient file:" << endl;
  int this_node, row,col;
  double latitude,longitude, this_x, this_y;;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "Y,X,latitude,longitude,elevation,flow_distance,drainage_area,segmented_gradient,source_key,basin_key";

  chi_data_out << endl;




  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    map<int,float>::iterator iter;

    for (iter = M_chi_data_map.begin(); iter != M_chi_data_map.end(); iter++)
    {
        this_node = iter->first;
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        get_x_and_y_locations(row, col, this_x, this_y);
        
        chi_data_out.precision(9);
        chi_data_out << this_y << ","
                     << this_x << ","
                     << latitude << ","
                     << longitude << ",";
        chi_data_out.precision(5);
        chi_data_out << elev_data_map[this_node] << ","
                     << flow_distance_data_map[this_node] << ","
                     << drainage_area_data_map[this_node] << ","
                     << M_chi_data_map[this_node] << ","
                     << source_keys_map[this_node] << ","
                     << baselevel_keys_map[this_node];
        chi_data_out << endl;
    }
  }

  chi_data_out.close();
  cout << "I am done, your file is:" << endl;
  cout << filename << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the fitted elevations: It uses m_chi and b_chi
// data to get the fitted elevation of the channel points.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_segmented_elevation(LSDFlowInfo& FlowInfo)
{
  // these are for extracting element-wise data from the channel profiles.
  int this_node;
  map<int,float> this_segmented_elevation_map;
  float this_M_chi, this_b_chi, this_chi, this_segemented_elevation;

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot calculate segments since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {

      // Get the M_chi and b_chi from the current node
      this_node = node_sequence[n];
      this_M_chi = M_chi_data_map[this_node];
      this_b_chi = b_chi_data_map[this_node];
      this_chi = chi_data_map[this_node];

      // calculate elevations simply based on the fact that we are fitting segments
      // with the equation z = M_chi*chi+b_chi
      this_segemented_elevation = this_M_chi*this_chi+this_b_chi;

      // Print the segment counter to the data map
      this_segmented_elevation_map[this_node]  = this_segemented_elevation;
    }
  }
  segmented_elevation_map = this_segmented_elevation_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This renumbers sources by base level. Each base level node has
// a number of sources and these get an incremental value which is used
// for the combination vector. These serve as maps between the source keys
// and the MLE values
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::baselevel_and_source_splitter(vector<int>& n_sources_for_baselevel,
                                                vector<int>& index_into_sources_vec)
{

  // You need to loop through all the sources, for each baselevel you get all the number rankings
  int n_sources = int(ordered_source_nodes.size());
  int baselevel_node;

  // this vector contains the index into the ordered source node vector
  // of the starting point for a given baselelvel node
  vector<int> starting_index_of_source_for_baselevel_node;
  vector<int> n_sources_each_baselevel;

  // loop through allsources, tracking where the baselevel nodes change.
  int this_baselevel_node = -1;
  int n_sources_this_baselevel;
  for (int i = 0; i< n_sources; i++)
  {
    // get the baselevel node of each of the sources
    baselevel_node = baselevel_keys_map[ ordered_source_nodes[i] ];

    if(baselevel_node != this_baselevel_node)
    {
      if (this_baselevel_node != -1)
      {
        n_sources_each_baselevel.push_back(n_sources_this_baselevel);
      }


      starting_index_of_source_for_baselevel_node.push_back(i);
      this_baselevel_node = baselevel_node;

      n_sources_this_baselevel = 0;
    }

    n_sources_this_baselevel++;
  }

  // and now get the last number of baselelvel nodes
  n_sources_each_baselevel.push_back(n_sources_this_baselevel);

  // now print out the results
  //cout << endl << endl << "============" << endl;
  bool print_for_debug = false;
  if(print_for_debug)
  {
    for(int i = 0; i< n_sources; i++)
   {
      cout << "Source number is: " << ordered_source_nodes[i] << " and baselelvel: " << baselevel_keys_map[ ordered_source_nodes[i] ] << endl;
    }

    int n_bl = int(starting_index_of_source_for_baselevel_node.size());
    cout << endl << endl << "============" << endl;
    cout << "n_bl: " << n_bl << endl;
    for(int i = 0; i< n_bl; i++)
    {
      cout << "Baselevel node is: " << ordered_baselevel_nodes[i] << " n sources: "
          << n_sources_each_baselevel[i] << " start_index: "
          << starting_index_of_source_for_baselevel_node[i] << endl;
    }

    cout << endl << endl << "============" << endl;
    cout << "Let me get the numbereing for you by basin" << endl;
    for(int i = 0; i< n_sources; i++)
    {
      cout << "Renumbered source key is: " << source_nodes_ranked_by_basin[i] << endl;
    }

  }

  // replace the two vectors
  n_sources_for_baselevel = n_sources_each_baselevel,
  index_into_sources_vec = starting_index_of_source_for_baselevel_node;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is used to debug the basin and source indexing
void LSDChiTools::print_basin_and_source_indexing_to_screen()
{
  cout << endl << endl << "=============" << endl << "printing for debugging" << endl;
  int n_sources = int(ordered_source_nodes.size());
  for(int i = 0; i< n_sources; i++)
  {
    cout << "Source number is: " << ordered_source_nodes[i] << " and baselelvel: " << baselevel_keys_map[ ordered_source_nodes[i] ] << endl;
  }

  vector<int> n_sources_for_baselevel;
  vector<int> index_into_sources_vec;
  baselevel_and_source_splitter(n_sources_for_baselevel,index_into_sources_vec);

  int n_bl = int(index_into_sources_vec.size());
  cout << endl << endl << "============" << endl;
  cout << "n_bl: " << n_bl << endl;
  for(int i = 0; i< n_bl; i++)
  {
    cout << "Baselevel node is: " << ordered_baselevel_nodes[i] << " n sources: "
        << n_sources_for_baselevel[i] << " start_index: "
        << index_into_sources_vec[i] << endl;
  }

  cout << endl << endl << "============" << endl;
  cout << "Let me get the numbering for you by basin" << endl;
  for(int i = 0; i< n_sources; i++)
  {
    cout << "Renumbered source key is: " << source_nodes_ranked_by_basin[i] << endl;
  }

  // now get the source keys
  cout << endl << endl << "========" << endl << "sources and keys" << endl;
  map<int,int>::iterator iter;
  iter = key_to_source_map.begin();
  while(iter != key_to_source_map.end())
  {
    cout << "source is: " << iter->first << " and key is: " << iter->second << endl;
    iter++;
  }

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::test_segment_collinearity(LSDFlowInfo& FlowInfo, int reference_channel, int test_channel,
                                             float sigma)
{
  // The way this works is that one of the segments (delineated by its source number)
  // is taken as a reference, and then all other segments are compared to how
  // closely they match this segment. If the chi value of the segment being tested
  // is greater than the maximum chi of the reference segment or less than the
  // minimum chi of the reference segment the data point is ignored.

  float MLE = 1;
  // first get the source node of the reference channel
  if ( reference_channel >= int(key_to_source_map.size()) || test_channel >= int(key_to_source_map.size()) )
  {
    cout << "LSDChiTools::test_segment_collinearity One of the sources is not in the channel network. Source is: " << reference_channel << endl;
  }
  else
  {
    vector<float> elev_data_chan0;
    vector<float> chi_data_chan0;
    get_chi_elevation_data_of_channel(FlowInfo, reference_channel, chi_data_chan0, elev_data_chan0);

    vector<float> elev_data_chan1;
    vector<float> chi_data_chan1;
    get_chi_elevation_data_of_channel(FlowInfo, test_channel, chi_data_chan1, elev_data_chan1);

    vector<float> residuals = project_data_onto_reference_channel(chi_data_chan0, elev_data_chan0,
                                 chi_data_chan1,elev_data_chan1);
    MLE = calculate_MLE_from_residuals(residuals, sigma);

  }
  return MLE;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::test_segment_collinearity_using_points(LSDFlowInfo& FlowInfo, int reference_channel, int test_channel,
                                             float sigma, vector<float> chi_distances_to_test)
{
  // The way this works is that one of the segments (delineated by its source number)
  // is taken as a reference, and then all other segments are compared to how
  // closely they match this segment. If the chi value of the segment being tested
  // is greater than the maximum chi of the reference segment or less than the
  // minimum chi of the reference segment the data point is ignored.

  float MLE = 1;
  // first get the source node of the reference channel
  if ( reference_channel >= int(key_to_source_map.size()) || test_channel >= int(key_to_source_map.size()) )
  {
    cout << "LSDChiTools::test_segment_collinearity One of the sources is not in the channel network. Source is: " << reference_channel << endl;
  }
  else
  {
    vector<float> elev_data_chan0;
    vector<float> chi_data_chan0;
    get_chi_elevation_data_of_channel(FlowInfo, reference_channel, chi_data_chan0, elev_data_chan0);

    vector<float> elev_data_chan1;
    vector<float> chi_data_chan1;
    get_chi_elevation_data_of_channel(FlowInfo, test_channel, chi_data_chan1, elev_data_chan1);

    vector<float> residuals = project_points_onto_reference_channel(chi_data_chan0, elev_data_chan0,
                                 chi_data_chan1,elev_data_chan1, chi_distances_to_test);

    // look through the residuals and remove any that are nodata
    int n_resid = int(residuals.size());
    vector<float> valid_residuals;
    for(int i = 0; i<n_resid; i++)
    {
      if (residuals[i] != -9999)
      {
        valid_residuals.push_back(residuals[i]);
      }
    }
    if (valid_residuals.size() != 0)
    {
      // only update MLE if there are residuals.
      MLE = calculate_MLE_from_residuals(residuals, sigma);
    }

  }
  return MLE;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDChiTools::retrieve_all_residuals_by_basin(LSDFlowInfo& FlowInfo, bool only_use_mainstem_as_reference,
                                                 int baselevel_key)
{
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }


  // run the splitter
  // this gets the starting index of the sources for each basin.
  // It means that the channel numbers are linked to the channels in the basin
  vector<int> start_node_for_baselelvel;
  vector<int> n_sources_in_basin;
  baselevel_and_source_splitter(n_sources_in_basin, start_node_for_baselelvel);
  int channel_offset = start_node_for_baselelvel[baselevel_key];
  int n_channels = n_sources_in_basin[baselevel_key];



  // now get all the possible two pair combinations of these channels
  bool zero_indexed = true;   // this is just because the channels are numbered from zero
  int k = 2;                  // We want combinations of 2 channels

  // the vec vec holds a vector of each possible combination of channels
  // each vector has two elements in it: the first and second channel in the comibination
  vector< vector<int> > combo_vecvev = combinations(n_channels, k, zero_indexed);
  vector<float> elev_data_chan0;
  vector<float> chi_data_chan0;
  vector<float> elev_data_chan1;
  vector<float> chi_data_chan1;
  vector<float> residuals;
  vector<float> cumulative_residuals;

  // Drop out if there is only a single channel in the basin
  if (n_channels == 1)
  {
    cout << "This basin only has one channel." << endl;
    residuals.push_back(0.0);
    return residuals;
  }


  vector<float> MLEs;
  vector<int> MLE_index;   // the index into the combo_vecvec that is used to
                           // tell which combinations have MLE values

  int last_ref_channel = -1;

  int n_combinations = int(combo_vecvev.size());
  vector<int> this_combo;
  int chan0,chan1;

  if (only_use_mainstem_as_reference)
  {
    n_combinations = n_channels-1;
  }

  // we loop through the different combinations in the vecvec
  for (int combo = 0; combo < n_combinations; combo++)
  {
    this_combo = combo_vecvev[combo];

    // you need to map these combinations onto the particular channels of this basin
    // These channels refere to the source keys
    chan0 = this_combo[0]+channel_offset;
    chan1 = this_combo[1]+channel_offset;
    //cout << "chan0 is: " << chan0 << "  and chan1 is: " << chan1 << " and combo 0 is: " << this_combo[0] <<endl;


    // only get the reference channel if the channel has changed.
    // This collects the chi-elevation data of the reference channel
    if (last_ref_channel != chan0)
    {
      get_chi_elevation_data_of_channel(FlowInfo, chan0, chi_data_chan0, elev_data_chan0);
    }

    // This gets the chi-elevation data of the test channel. Again, the chan1
    // parameter referes to the source key.
    get_chi_elevation_data_of_channel(FlowInfo, chan1, chi_data_chan1, elev_data_chan1);

    // Now return the residuals between the reference channel and test channel.
    // Each node in the test channel gets a residual, it is projected to a
    // linear fit between nodes on the reference channel
    residuals = project_data_onto_reference_channel(chi_data_chan0, elev_data_chan0,
                                 chi_data_chan1,elev_data_chan1);

    // append the residual vector
    cumulative_residuals.insert(cumulative_residuals.end(), residuals.begin(), residuals.end());

  }

  return cumulative_residuals;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::test_all_segment_collinearity_by_basin(LSDFlowInfo& FlowInfo, bool only_use_mainstem_as_reference,
                                                 int baselevel_key,
                                                 vector<int>& reference_source, vector<int>& test_source,
                                                 vector<float>& MLE_values, vector<float>& RMSE_values,
                                                 float sigma)
{
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }


  // run the splitter
  // this gets the starting index of the sources for each basin.
  // It means that the channel numbers are linked to the channels in the basin
  vector<int> start_node_for_baselelvel;
  vector<int> n_sources_in_basin;
  baselevel_and_source_splitter(n_sources_in_basin, start_node_for_baselelvel);
  int channel_offset = start_node_for_baselelvel[baselevel_key];
  int n_channels = n_sources_in_basin[baselevel_key];

  // Drop out if there is only a single channel in the basin
  if (n_channels == 1)
  {
    cout << "This basin only has one channel." << endl;
    return 1.0;
  }


  //cout << "Let me check the basin indexing for you." << endl;
  //int n_bl =  int(start_node_for_baselelvel.size());
  //for(int i = 0; i< n_bl; i++)
  //{
  //  cout << "Baselevel node is: " << ordered_baselevel_nodes[i] << " n sources: "
  //      << n_sources_in_basin[i] << " start_index: "
  //      << start_node_for_baselelvel[i] << endl;
  //}
  //cout << "Basin key is: " << baselevel_key << " This basin has " << n_channels << " sources and start node is "<< channel_offset <<endl;

  // placeholder vectors: will replace the passed vectors
  vector<int> this_reference_source;
  vector<int> this_test_source;
  vector<float> these_MLE_values;
  vector<float> these_RMSE_values;

  // now get all the possible two pair combinations of these channels
  bool zero_indexed = true;   // this is just because the channels are numbered from zero
  int k = 2;                  // We want combinations of 2 channels

  // the vec vec holds a vector of each possible combination of channels
  // each vector has two elements in it: the first and second channel in the comibination
  vector< vector<int> > combo_vecvev = combinations(n_channels, k, zero_indexed);
  vector<float> elev_data_chan0;
  vector<float> chi_data_chan0;
  vector<float> elev_data_chan1;
  vector<float> chi_data_chan1;
  vector<float> residuals;

  int n_residuals;



  vector<float> MLEs;
  vector<int> MLE_index;   // the index into the combo_vecvec that is used to
                           // tell which combinations have MLE values

  int last_ref_channel = -1;

  int n_combinations = int(combo_vecvev.size());
  vector<int> this_combo;
  int chan0,chan1;

  if (only_use_mainstem_as_reference)
  {
    n_combinations = n_channels-1;
  }

  // we loop through the different combinations in the vecvec
  for (int combo = 0; combo < n_combinations; combo++)
  {
    this_combo = combo_vecvev[combo];

    // you need to map these combinations onto the particular channels of this basin
    // These channels refere to the source keys
    chan0 = this_combo[0]+channel_offset;
    chan1 = this_combo[1]+channel_offset;
    //cout << "chan0 is: " << chan0 << "  and chan1 is: " << chan1 << " and combo 0 is: " << this_combo[0] <<endl;


    // only get the reference channel if the channel has changed.
    // This collects the chi-elevation data of the reference channel
    if (last_ref_channel != chan0)
    {
      get_chi_elevation_data_of_channel(FlowInfo, chan0, chi_data_chan0, elev_data_chan0);
    }

    // This gets the chi-elevation data of the test channel. Again, the chan1
    // parameter referes to the source key.
    get_chi_elevation_data_of_channel(FlowInfo, chan1, chi_data_chan1, elev_data_chan1);

    // Now return the residuals between the reference channel and test channel.
    // Each node in the test channel gets a residual, it is projected to a
    // linear fit between nodes on the reference channel
    residuals = project_data_onto_reference_channel(chi_data_chan0, elev_data_chan0,
                                 chi_data_chan1,elev_data_chan1);
    n_residuals = int(residuals.size());
    //cout << "Basin: " << baselevel_key << " The number of residuals are: " << n_residuals << endl;

    // Now get the MLE and RMSE for this channel pair. It only runs if
    // there are residuals. Otherwise it means that the channels are non-overlapping
    if (n_residuals > 0)
    {
      float MLE1 = calculate_MLE_from_residuals(residuals, sigma);
      float RMSE = calculate_RMSE_from_residuals(residuals);
      last_ref_channel = chan0;
      //cout << "MLE: " << MLE1 << " and RMSE: " << RMSE << endl;

      // If we are only using the mainstem channel, we only use the first channel
      // as a reference channel. The first channel is denoted by this_combo[0] == 0
      //cout << "The use only mainstem is: " << only_use_mainstem_as_reference << endl;
      if (only_use_mainstem_as_reference)
      {
        //cout << "Checking the combination, combo 0 is: " << this_combo[0] << endl;
        if (this_combo[0] > 0)
        {
          // skip to the last node
          //cout << "I am skipping to the LAST NODE in the combinations" << endl;
          combo = n_combinations;
        }
        else
        {
          these_MLE_values.push_back(MLE1);
          these_RMSE_values.push_back(RMSE);
          this_reference_source.push_back(chan0);
          this_test_source.push_back(chan1);
        }
      }
      else
      {
        these_MLE_values.push_back(MLE1);
        these_RMSE_values.push_back(RMSE);
        this_reference_source.push_back(chan0);
        this_test_source.push_back(chan1);
      }
    }
    else
    {
        these_MLE_values.push_back(1.0);
        these_RMSE_values.push_back(0.0);
        this_reference_source.push_back(chan0);
        this_test_source.push_back(chan1);
    }
  }

  MLE_values = these_MLE_values;
  RMSE_values = these_RMSE_values;
  reference_source = this_reference_source;
  test_source = this_test_source;

  float tot_MLE = 1;
  for (int res = 0; res < int(these_MLE_values.size()); res++)
  {
    tot_MLE = tot_MLE*these_MLE_values[res];
  }


  //cout << "Let me tell you all about the MLE values " << endl;
  // for debugging
  bool print_results = false;
  if(print_results)
  {
    for (int res = 0; res < int(MLE_values.size()); res++)
    {
      cout << "reference_source: " << reference_source[res] << " "
           << "test_source: " << test_source[res] << " "
           << "MLE_values: " << MLE_values[res] << " "
           << "RMSE_values: " << RMSE_values[res] << endl;
    }
  }

  //cout << "N_residuals: " << MLE_values.size() << endl << endl;

  //cout << "Total MLE is: " << tot_MLE << endl;
  return tot_MLE;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::test_all_segment_collinearity_by_basin_using_points(LSDFlowInfo& FlowInfo, bool only_use_mainstem_as_reference,
                                                 int baselevel_key,
                                                 vector<int>& reference_source, vector<int>& test_source,
                                                 vector<float>& MLE_values, vector<float>& RMSE_values,
                                                 float sigma, vector<float> chi_fractions_for_testing)
{
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << " and sigma is " << sigma << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }


  // run the splitter
  // this gets the starting index of the sources for each basin.
  // It means that the channel numbers are linked to the channels in the basin
  vector<int> start_node_for_baselelvel;
  vector<int> n_sources_in_basin;
  baselevel_and_source_splitter(n_sources_in_basin, start_node_for_baselelvel);
  int channel_offset = start_node_for_baselelvel[baselevel_key];
  int n_channels = n_sources_in_basin[baselevel_key];

  // Drop out if there is only a single channel in the basin
  if (n_channels == 1)
  {
    cout << "This basin only has one channel." << endl;
    return 1.0;
  }

  // we need to know the length of the mainstem.
  int MS_channel = channel_offset;
  vector<float> MS_chi;
  vector<float> MS_elev;
  get_chi_elevation_data_of_channel(FlowInfo, MS_channel, MS_chi, MS_elev);

  int n_MS_nodes = int(MS_chi.size());
  float MS_length = MS_chi[0]-MS_chi[n_MS_nodes-1];
  //cout << "LSDChiTools::test_all_segment_collinearity_by_basin_using_points The mainstem has a length of: " << MS_length << endl;
  vector<float> chi_test_distances;
  int n_frac = int(chi_fractions_for_testing.size());
  float this_distance;
  for(int f = 0; f<n_frac; f++)
  {
    this_distance = chi_fractions_for_testing[f]*MS_length;
    chi_test_distances.push_back(this_distance);
    //cout << "The frac is: " << chi_fractions_for_testing[f] << " and distance: " << chi_test_distances[f] << endl;
  }




  //cout << "Let me check the basin indexing for you." << endl;
  //int n_bl =  int(start_node_for_baselelvel.size());
  //for(int i = 0; i< n_bl; i++)
  //{
  //  cout << "Baselevel node is: " << ordered_baselevel_nodes[i] << " n sources: "
  //      << n_sources_in_basin[i] << " start_index: "
  //      << start_node_for_baselelvel[i] << endl;
  //}
  //cout << "Basin key is: " << baselevel_key << " This basin has " << n_channels << " sources and start node is "<< channel_offset <<endl;

  // placeholder vectors: will replace the passed vectors
  vector<int> this_reference_source;
  vector<int> this_test_source;
  vector<float> these_MLE_values;
  vector<float> these_RMSE_values;

  // now get all the possible two pair combinations of these channels
  bool zero_indexed = true;   // this is just because the channels are numbered from zero
  int k = 2;                  // We want combinations of 2 channels

  // the vec vec holds a vector of each possible combination of channels
  // each vector has two elements in it: the first and second channel in the comibination
  vector< vector<int> > combo_vecvev = combinations(n_channels, k, zero_indexed);
  vector<float> elev_data_chan0;
  vector<float> chi_data_chan0;
  vector<float> elev_data_chan1;
  vector<float> chi_data_chan1;
  vector<float> residuals;

  int n_residuals;

  vector<float> MLEs;
  vector<int> MLE_index;   // the index into the combo_vecvec that is used to
                           // tell which combinations have MLE values

  int last_ref_channel = -1;

  int n_combinations = int(combo_vecvev.size());
  vector<int> this_combo;
  int chan0,chan1;

  if (only_use_mainstem_as_reference)
  {
    n_combinations = n_channels-1;
  }

  // we loop through the different combinations in the vecvec
  //cout << "LSDChiTools::test_all_segment_collinearity_by_basin_using_points This basin has " << n_combinations << " combinations." << endl;
  for (int combo = 0; combo < n_combinations; combo++)
  {
    this_combo = combo_vecvev[combo];

    // you need to map these combinations onto the particular channels of this basin
    // These channels refere to the source keys
    chan0 = this_combo[0]+channel_offset;
    chan1 = this_combo[1]+channel_offset;
    //cout << "chan0 is: " << chan0 << "  and chan1 is: " << chan1 << " and combo 0 is: " << this_combo[0] <<endl;


    // only get the reference channel if the channel has changed.
    // This collects the chi-elevation data of the reference channel
    if (last_ref_channel != chan0)
    {
      get_chi_elevation_data_of_channel(FlowInfo, chan0, chi_data_chan0, elev_data_chan0);
    }

    // This gets the chi-elevation data of the test channel. Again, the chan1
    // parameter referes to the source key.
    get_chi_elevation_data_of_channel(FlowInfo, chan1, chi_data_chan1, elev_data_chan1);

    // Now return the residuals between the reference channel and test channel.
    // Each node in the test channel gets a residual, it is projected to a
    // linear fit between nodes on the reference channel
    residuals = project_points_onto_reference_channel(chi_data_chan0, elev_data_chan0,
                                 chi_data_chan1,elev_data_chan1, chi_test_distances);
    n_residuals = int(residuals.size());
    //cout << "LSDChiTools::test_all_segment_collinearity_by_basin_using_points Basin: " << baselevel_key << " The number of residuals are: " << n_residuals << endl;
    //for (int i = 0; i< n_residuals; i++)
    //{
    //  cout << "residual[" << i<<"]: " << residuals[i] << endl;
    //}


    // Now get the MLE and RMSE for this channel pair. It only runs if
    // there are residuals. Otherwise it means that the channels are non-overlapping
    if (n_residuals > 0)
    {
      float MLE1 = calculate_MLE_from_residuals(residuals, sigma);
      float RMSE = calculate_RMSE_from_residuals(residuals);
      last_ref_channel = chan0;
      //cout << "sigma is: " << sigma << " MLE: " << MLE1 << " and RMSE: " << RMSE << endl;

      // If we are only using the mainstem channel, we only use the first channel
      // as a reference channel. The first channel is denoted by this_combo[0] == 0
      //cout << "The use only mainstem is: " << only_use_mainstem_as_reference << endl;
      if (only_use_mainstem_as_reference)
      {
        //cout << "Checking the combination, combo 0 is: " << this_combo[0] << endl;
        if (this_combo[0] > 0)
        {
          // skip to the last node
          //cout << "I am skipping to the LAST NODE in the combinations" << endl;
          combo = n_combinations;
        }
        else
        {
          these_MLE_values.push_back(MLE1);
          these_RMSE_values.push_back(RMSE);
          this_reference_source.push_back(chan0);
          this_test_source.push_back(chan1);
        }
      }
      else
      {
        these_MLE_values.push_back(MLE1);
        these_RMSE_values.push_back(RMSE);
        this_reference_source.push_back(chan0);
        this_test_source.push_back(chan1);
      }
    }
    else
    {
        these_MLE_values.push_back(1.0);
        these_RMSE_values.push_back(0.0);
        this_reference_source.push_back(chan0);
        this_test_source.push_back(chan1);
    }
  }

  MLE_values = these_MLE_values;
  RMSE_values = these_RMSE_values;
  reference_source = this_reference_source;
  test_source = this_test_source;

  float tot_MLE = 1;
  for (int res = 0; res < int(these_MLE_values.size()); res++)
  {
    tot_MLE = tot_MLE*these_MLE_values[res];
  }


  //cout << "Let me tell you all about the MLE values " << endl;
  // for debugging
  bool print_results = false;
  if(print_results)
  {
    for (int res = 0; res < int(MLE_values.size()); res++)
    {
      cout << "reference_source: " << reference_source[res] << " "
           << "test_source: " << test_source[res] << " "
           << "MLE_values: " << MLE_values[res] << " "
           << "RMSE_values: " << RMSE_values[res] << endl;
    }
  }

  //cout << "N_residuals: " << MLE_values.size() << endl << endl;

  //cout << "Total MLE is: " << tot_MLE << endl;
  return tot_MLE;

}









//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::test_collinearity_by_basin_disorder(LSDFlowInfo& FlowInfo,
                                                 int baselevel_key)
{
  float disorder_stat = -9999;
  
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  
  // Drop out if there is only a single channel in the basin
  //if (n_channels == 1)
  //{
  //  cout << "This basin only has one channel." << endl;
  //  return 1.0;
  //}


  // This is a brute force way to get the complete chi data map
  
  vector<float> this_basin_chi;
  vector<float> this_basin_elevation;
  int n_nodes = int(node_sequence.size());
  int this_node;
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    if (baselevel_keys_map[this_node] == baselevel_key)
    {

      this_basin_chi.push_back(chi_data_map[this_node]);
      this_basin_elevation.push_back(elev_data_map[this_node]);
    }
  }
  
  // now sort these vectors
    // initiate the sorted vectors
  vector<float> chi_sorted;
  vector<float> elev_sorted;
  vector<size_t> index_map;

  // sort the vectors
  matlab_float_sort(this_basin_elevation, elev_sorted, index_map);
  matlab_float_reorder(this_basin_chi, index_map, chi_sorted);
  
  // now calculate disorder
  float this_delta_chi = 0;
  float sum_delta_chi = 0;
  
  int n_nodes_this_basin = int(chi_sorted.size());
  float chi_max = chi_sorted[n_nodes_this_basin-1];
  float chi_min = chi_sorted[0];
  // Sorry Simon, the following cout slows down the code when many rivers are tested: terminal output overflow
  // cout << "My first guess of minimum chi is: " << chi_min << endl;

  for(int i = 0; i<n_nodes_this_basin-1; i++)
  {
    this_delta_chi = fabs(chi_sorted[i+1]-chi_sorted[i]);
    sum_delta_chi+=this_delta_chi;
    if(chi_sorted[i] > chi_max)
    {
      chi_max = chi_sorted[i];
    }
    if(chi_sorted[i] < chi_min)
    {
      chi_min = chi_sorted[i];
    }
  }
  if(chi_sorted[n_nodes_this_basin-1] > chi_max)
  {
    chi_max = chi_sorted[n_nodes_this_basin-1];
  }
  float chi_range = chi_max-chi_min;
  
  disorder_stat = (sum_delta_chi - chi_range)/chi_range;

  //cout << "Let me tell you all about the MLE values " << endl;
  // for debugging
  bool print_results = false;
  if(print_results)
  {
    cout << "I need to code this up." << endl;
  }

  return disorder_stat;

}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDChiTools::test_collinearity_by_basin_disorder_with_uncert(LSDFlowInfo& FlowInfo,
                                                 int baselevel_key)
{
  vector<float> disorder_stat_vec;
  
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  // this map will hold references to the sources. There are two of them
  // because the first has the sources as the keys and the second has the
  // combination index as the keys. 
  map<int,int> sources_are_keys;
  map<int,int> comboindex_are_keys;
  
  vector<float> this_basin_chi;
  vector<float> this_basin_elevation;
  vector<int> this_basin_source;
  int n_nodes = int(node_sequence.size());
  int this_node;
  int comboindex = 0;     // this is used to store an index into the combinations
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    if (baselevel_keys_map[this_node] == baselevel_key)
    {

      this_basin_chi.push_back(chi_data_map[this_node]);
      this_basin_elevation.push_back(elev_data_map[this_node]);
      this_basin_source.push_back(source_keys_map[this_node]);
      
      // if the key doesn't exist, add a source key counter
      if ( sources_are_keys.find( source_keys_map[this_node] ) == sources_are_keys.end() )
      {
        sources_are_keys[ source_keys_map[this_node] ] = comboindex;
        comboindex++;
      }
    }
  }
  
  // now do the second map by inverting the first map
  for(map<int,int>::iterator iter =sources_are_keys.begin(); iter != sources_are_keys.end(); ++iter)
  {
    int k =  iter->first;
    int v = iter->second;
    comboindex_are_keys[v] = k;
  }
  
  int n_sources = int(comboindex_are_keys.size());
  //cout << "The number of channels are: " << n_sources << endl;
  //cout << "The source node is: " << comboindex_are_keys[0] << endl;
  
  
  // get the combinations
  //cout << "Let me get some combinations for you." << endl;
  int n_elements = n_sources-1;
  int n_in_each_combo = 3;
  
  if (n_sources < n_in_each_combo)
  {
    cout << "Not enough channels in this basin! I am returning a nodata vector." << endl;
    disorder_stat_vec.push_back(-9999);
  }
  else
  {
    // this gets the combinations
    bool zero_indexed = false;
    vector< vector<int> > combo_vecvec = combinations(n_elements, n_in_each_combo, zero_indexed);
    
    bool print_combinations = false;
    if(print_combinations)
    {
      for (int i = 0; i< int(combo_vecvec.size()); i++)
      {
        for (int j = 0; j< int(combo_vecvec[i].size()); j++)
        {
          cout << combo_vecvec[i][j] << " ";
        }
        cout << endl;
      }
    }
    int n_combinations = int(combo_vecvec.size());
  
    // now you enter the combinations loop
    //int n_data_points = int(this_basin_source.size());
    //cout << "Checking the combinations" << endl;
    for (int combo = 0; combo<n_combinations; combo++)
    {
      vector<int> these_combos = combo_vecvec[combo];
      vector<int> these_combo_sources;
      
      // add the trunk channel source
      these_combo_sources.push_back(comboindex_are_keys[0]); 
      //cout << 0 <<":" << comboindex_are_keys[0] << " ";
      
      for (int source_key = 0 ; source_key < int(these_combos.size()); source_key++)
      {
        these_combo_sources.push_back( comboindex_are_keys[ these_combos[source_key] ] );
        //cout << these_combos[source_key] <<":" << comboindex_are_keys[ these_combos[source_key] ] << " ";
      }
      //cout << endl;
      
      
      // now we sample and sort these vectors
      vector<float> chi_combo;
      vector<float> elev_combo;
      vector<float> source_combo;
      
      for(int node = 0; node< int(this_basin_chi.size()); node++)
      {
        if(find(these_combo_sources.begin(), these_combo_sources.end(), this_basin_source[node]) != these_combo_sources.end()) 
        {
          chi_combo.push_back(this_basin_chi[node]);
          elev_combo.push_back(this_basin_elevation[node]);
        }
      }
        
      // now you have the elevation and chi for this combination, sort them
      // initiate the sorted vectors
      vector<float> chi_sorted;
      vector<float> elev_sorted;
      vector<size_t> index_map;    
  
      // sort the vectors
      matlab_float_sort(elev_combo, elev_sorted, index_map);
      matlab_float_reorder(chi_combo, index_map, chi_sorted);
        
      // now calculate disorder
      float this_delta_chi = 0;
      float sum_delta_chi = 0;
      
      int n_nodes_this_basin = int(chi_sorted.size());
      float chi_max = chi_sorted[n_nodes_this_basin-1];
      float chi_min = chi_sorted[0];
      // Sorry Simon, the following cout slows down the code when many rivers are tested: terminal output overflow
      // cout << "My first guess of minimum chi is: " << chi_min << endl;
      for(int i = 0; i<n_nodes_this_basin-1; i++)
      {
        this_delta_chi = fabs(chi_sorted[i+1]-chi_sorted[i]);
        sum_delta_chi+=this_delta_chi;
        if(chi_sorted[i] > chi_max)
        {
          chi_max = chi_sorted[i];
        }
        if(chi_sorted[i] < chi_min)
        {
          chi_min = chi_sorted[i];
        }
      }
      if(chi_sorted[n_nodes_this_basin-1] > chi_max)
      {
        chi_max = chi_sorted[n_nodes_this_basin-1];
      }
      float chi_range = chi_max-chi_min;
        
      float disorder_stat = (sum_delta_chi - chi_range)/chi_range;
      disorder_stat_vec.push_back(disorder_stat);
    }
  }


  // for debugging
  bool print_results = false;
  if(print_results)
  {
    cout << "Disorder stats for all the combinations are: " << endl;
    for(int i = 0; i< int(disorder_stat_vec.size()); i++)
    {
      cout << disorder_stat_vec[i] << " ";
    }
    cout << endl;
    cout << "I'm returning the disorder stat vector." << endl;
  }
  

  return disorder_stat_vec;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function reproduce the bug from the disorder metric linked to chi_min
// Used to assess the extent of the bug and make sure it is sorted
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDChiTools::TEST_FUNCTION_OLD_DISORDER_DO_NOT_USE(LSDFlowInfo& FlowInfo,
                                                 int baselevel_key)
{
  vector<float> disorder_stat_vec;
  
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" << endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  // this map will hold references to the sources. There are two of them
  // because the first has the sources as the keys and the second has the
  // combination index as the keys. 
  map<int,int> sources_are_keys;
  map<int,int> comboindex_are_keys;
  
  vector<float> this_basin_chi;
  vector<float> this_basin_elevation;
  vector<int> this_basin_source;
  int n_nodes = int(node_sequence.size());
  int this_node;
  int comboindex = 0;     // this is used to store an index into the combinations
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    if (baselevel_keys_map[this_node] == baselevel_key)
    {

      this_basin_chi.push_back(chi_data_map[this_node]);
      this_basin_elevation.push_back(elev_data_map[this_node]);
      this_basin_source.push_back(source_keys_map[this_node]);
      
      // if the key doesn't exist, add a source key counter
      if ( sources_are_keys.find( source_keys_map[this_node] ) == sources_are_keys.end() )
      {
        sources_are_keys[ source_keys_map[this_node] ] = comboindex;
        comboindex++;
      }
    }
  }
  
  // now do the second map by inverting the first map
  for(map<int,int>::iterator iter =sources_are_keys.begin(); iter != sources_are_keys.end(); ++iter)
  {
    int k =  iter->first;
    int v = iter->second;
    comboindex_are_keys[v] = k;
  }
  
  int n_sources = int(comboindex_are_keys.size());
  //cout << "The number of channels are: " << n_sources << endl;
  //cout << "The source node is: " << comboindex_are_keys[0] << endl;
  
  
  // get the combinations
  //cout << "Let me get some combinations for you." << endl;
  int n_elements = n_sources-1;
  int n_in_each_combo = 3;
  
  if (n_sources < n_in_each_combo)
  {
    cout << "Not enough channels in this basin! I am returning a nodata vector." << endl;
    disorder_stat_vec.push_back(-9999);
  }
  else
  {
    // this gets the combinations
    bool zero_indexed = false;
    vector< vector<int> > combo_vecvec = combinations(n_elements, n_in_each_combo, zero_indexed);
    
    bool print_combinations = false;
    if(print_combinations)
    {
      for (int i = 0; i< int(combo_vecvec.size()); i++)
      {
        for (int j = 0; j< int(combo_vecvec[i].size()); j++)
        {
          cout << combo_vecvec[i][j] << " ";
        }
        cout << endl;
      }
    }
    int n_combinations = int(combo_vecvec.size());
  
    // now you enter the combinations loop
    //int n_data_points = int(this_basin_source.size());
    //cout << "Checking the combinations" << endl;
    for (int combo = 0; combo<n_combinations; combo++)
    {
      vector<int> these_combos = combo_vecvec[combo];
      vector<int> these_combo_sources;
      
      // add the trunk channel source
      these_combo_sources.push_back(comboindex_are_keys[0]); 
      //cout << 0 <<":" << comboindex_are_keys[0] << " ";
      
      for (int source_key = 0 ; source_key < int(these_combos.size()); source_key++)
      {
        these_combo_sources.push_back( comboindex_are_keys[ these_combos[source_key] ] );
        //cout << these_combos[source_key] <<":" << comboindex_are_keys[ these_combos[source_key] ] << " ";
      }
      //cout << endl;
      
      
      // now we sample and sort these vectors
      vector<float> chi_combo;
      vector<float> elev_combo;
      vector<float> source_combo;
      
      for(int node = 0; node< int(this_basin_chi.size()); node++)
      {
        if(find(these_combo_sources.begin(), these_combo_sources.end(), this_basin_source[node]) != these_combo_sources.end()) 
        {
          chi_combo.push_back(this_basin_chi[node]);
          elev_combo.push_back(this_basin_elevation[node]);
        }
      }
        
      // now you have the elevation and chi for this combination, sort them
      // initiate the sorted vectors
      vector<float> chi_sorted;
      vector<float> elev_sorted;
      vector<size_t> index_map;    
  
      // sort the vectors
      matlab_float_sort(elev_combo, elev_sorted, index_map);
      matlab_float_reorder(chi_combo, index_map, chi_sorted);
        
      // now calculate disorder
      float this_delta_chi = 0;
      float sum_delta_chi = 0;
      
      int n_nodes_this_basin = int(chi_sorted.size());
      float corrected_chi_max = chi_sorted[n_nodes_this_basin-1];
      float corrected_chi_min = chi_sorted[0];
      float chi_min = 10000;
      float chi_max = 0;
      if(chi_min != corrected_chi_min)  
        cout << "BUGDETECTED:: CHI MIN=" << corrected_chi_min <<  " BUT CHI_MIN_USED="  << chi_min << endl;

      for(int i = 0; i<n_nodes_this_basin-1; i++)
      {
        this_delta_chi = fabs(chi_sorted[i+1]-chi_sorted[i]);
        sum_delta_chi+=this_delta_chi;
        if(chi_sorted[i] > chi_max)
        {
          chi_max = chi_sorted[i];
        }
        if(chi_sorted[i] < chi_min)
        {
          chi_min = chi_sorted[i];
        }
      }
      if(chi_sorted[n_nodes_this_basin-1] > chi_max)
      {
        chi_max = chi_sorted[n_nodes_this_basin-1];
      }
      float chi_range = chi_max-chi_min;
        
      float disorder_stat = (sum_delta_chi - chi_range)/chi_range;
      disorder_stat_vec.push_back(disorder_stat);
    }
  }


  // for debugging
  bool print_results = false;
  if(print_results)
  {
    cout << "Disorder stats for all the combinations are: " << endl;
    for(int i = 0; i< int(disorder_stat_vec.size()); i++)
    {
      cout << disorder_stat_vec[i] << " ";
    }
    cout << endl;
    cout << "I'm returning the disorder stat vector." << endl;
  }
  

  return disorder_stat_vec;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment, It saves the source keys involved in each combination to get an idea of spatial distributions - B.G 2019
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map< vector<int>, vector<float> > LSDChiTools::test_collinearity_by_basin_disorder_with_uncert_retain_key(LSDFlowInfo& FlowInfo,
                                                 int baselevel_key)
{
  vector<float> disorder_stat_vec;
  vector<vector<int> > disorder_stat_vec_associated_SK;
  map< vector<int>, vector<float> > output;


  
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  // this map will hold references to the sources. There are two of them
  // because the first has the sources as the keys and the second has the
  // combination index as the keys. 
  map<int,int> sources_are_keys;
  map<int,int> comboindex_are_keys;
  
  vector<float> this_basin_chi;
  vector<float> this_basin_elevation;
  vector<int> this_basin_source;
  int n_nodes = int(node_sequence.size());
  int this_node;
  int comboindex = 0;     // this is used to store an index into the combinations
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    if (baselevel_keys_map[this_node] == baselevel_key)
    {

      this_basin_chi.push_back(chi_data_map[this_node]);
      this_basin_elevation.push_back(elev_data_map[this_node]);
      this_basin_source.push_back(source_keys_map[this_node]);
      
      // if the key doesn't exist, add a source key counter
      if ( sources_are_keys.find( source_keys_map[this_node] ) == sources_are_keys.end() )
      {
        sources_are_keys[ source_keys_map[this_node] ] = comboindex;
        comboindex++;
      }
    }
  }
  
  // now do the second map by inverting the first map
  for(map<int,int>::iterator iter =sources_are_keys.begin(); iter != sources_are_keys.end(); ++iter)
  {
    int k =  iter->first;
    int v = iter->second;
    comboindex_are_keys[v] = k;
  }
  
  int n_sources = int(comboindex_are_keys.size());
  //cout << "The number of channels are: " << n_sources << endl;
  //cout << "The source node is: " << comboindex_are_keys[0] << endl;
  
  
  // get the combinations
  //cout << "Let me get some combinations for you." << endl;
  int n_elements = n_sources-1;
  int n_in_each_combo = 3;
  
  if (n_sources < n_in_each_combo)
  {
    cout << "Not enough channels in this basin! I am returning a nodata vector." << endl;
    disorder_stat_vec.push_back(-9999);
    disorder_stat_vec_associated_SK.push_back({-9999});
  }
  else
  {
    // this gets the combinations
    bool zero_indexed = false;
    vector< vector<int> > combo_vecvec = combinations(n_elements, n_in_each_combo, zero_indexed);
    
    bool print_combinations = false;
    if(print_combinations)
    {
      for (int i = 0; i< int(combo_vecvec.size()); i++)
      {
        for (int j = 0; j< int(combo_vecvec[i].size()); j++)
        {
          cout << combo_vecvec[i][j] << " ";
        }
        cout << endl;
      }
    }
    int n_combinations = int(combo_vecvec.size());
    
    // now you enter the combinations loop
    //int n_data_points = int(this_basin_source.size());
    //cout << "Checking the combinations" << endl;
    for (int combo = 0; combo<n_combinations; combo++)
    {
      vector<int> these_combos = combo_vecvec[combo];
      vector<int> these_combo_sources;
      
      // add the trunk channel source
      these_combo_sources.push_back(comboindex_are_keys[0]); 
      
      for (int source_key = 0 ; source_key < int(these_combos.size()); source_key++)
      {
        these_combo_sources.push_back( comboindex_are_keys[ these_combos[source_key] ] );
      }
      disorder_stat_vec_associated_SK.push_back(these_combo_sources);
            
      
      // now we sample and sort these vectors
      vector<float> chi_combo;
      vector<float> elev_combo;
      vector<float> source_combo;
      
      for(int node = 0; node< int(this_basin_chi.size()); node++)
      {
        if(find(these_combo_sources.begin(), these_combo_sources.end(), this_basin_source[node]) != these_combo_sources.end()) 
        {
          chi_combo.push_back(this_basin_chi[node]);
          elev_combo.push_back(this_basin_elevation[node]);
        }
      }


        
      // now you have the elevation and chi for this combination, sort them
      // initiate the sorted vectors
      vector<float> chi_sorted;
      vector<float> elev_sorted;
      vector<size_t> index_map; 
 
      // sort the vectors
      matlab_float_sort(elev_combo, elev_sorted, index_map);
      matlab_float_reorder(chi_combo, index_map, chi_sorted);


      // now calculate disorder
      float this_delta_chi = 0;
      float sum_delta_chi = 0;
      
      int n_nodes_this_basin = int(chi_sorted.size());
      float chi_max = chi_sorted[n_nodes_this_basin-1];
      float chi_min = chi_sorted[0];
      // Sorry Simon, the following cout slows down the code when many rivers are tested: terminal output overflow
      // cout << "My first guess of minimum chi is: " << chi_min << endl;
      for(int i = 0; i<n_nodes_this_basin-1; i++)
      {
        this_delta_chi = fabs(chi_sorted[i+1]-chi_sorted[i]);
        sum_delta_chi += this_delta_chi;
        if(chi_sorted[i] > chi_max)
        {
          chi_max = chi_sorted[i];
        }
        if(chi_sorted[i] < chi_min)
        {
          chi_min = chi_sorted[i];
        }
      }
      if(chi_sorted[n_nodes_this_basin-1] > chi_max)
      {
        chi_max = chi_sorted[n_nodes_this_basin-1];
      }
      float chi_range = chi_max-chi_min;
        
      float disorder_stat = (sum_delta_chi - chi_range)/chi_range;
      disorder_stat_vec.push_back(disorder_stat);
      // std::cout << "distat is " << chi_sorted.size()  << std::endl;
    }
  }



  // for debugging
  bool print_results = false;
  if(print_results)
  {
    cout << "Disorder stats for all the combinations are: " << endl;
    for(int i = 0; i< int(disorder_stat_vec.size()); i++)
    {
      cout << disorder_stat_vec[i] << " ";
    }
    cout << endl;
    cout << "I'm returning the disorder stat vector." << endl;
  }
  
  for(size_t i=0; i<disorder_stat_vec_associated_SK.size();i++)
  {
    if(output.find( disorder_stat_vec_associated_SK[i] ) == output.end())
    {
      output[ disorder_stat_vec_associated_SK[i] ] = {disorder_stat_vec[i]};
    }
    else
    {
      output[ disorder_stat_vec_associated_SK[i] ].push_back(disorder_stat_vec[i]);
    }
  }

  return output;

}

void LSDChiTools::precombine_sources_for_disorder_with_uncert_opti(LSDFlowInfo& FlowInfo,int baselevel_key, 
  map<int,int>& sources_are_keys, map<int,int>& comboindex_are_keys, vector<vector<int> >& combo_vecvec, vector<int>& nodes_in_basin
  , vector<int>& this_basin_source, int n_in_each_combo)
{

  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  int n_nodes = int(node_sequence.size());
  int this_node;
  int comboindex = 0;     // this is used to store an index into the combinations
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];
    this_basin_source.push_back(source_keys_map[this_node]); 

    if (baselevel_keys_map[this_node] == baselevel_key)
    {
      nodes_in_basin.push_back(this_node);
      // if the key doesn't exist, add a source key counter
      if ( sources_are_keys.find( source_keys_map[this_node] ) == sources_are_keys.end() )
      {
        sources_are_keys[ source_keys_map[this_node] ] = comboindex;     
        comboindex++;
      }
    }
  }
  
  // now do the second map by inverting the first map
  for(map<int,int>::iterator iter =sources_are_keys.begin(); iter != sources_are_keys.end(); ++iter)
  {
    int k =  iter->first;
    int v = iter->second;
    comboindex_are_keys[v] = k;
  }
  
  int n_sources = int(comboindex_are_keys.size());
  //cout << "The number of channels are: " << n_sources << endl;
  //cout << "The source node is: " << comboindex_are_keys[0] << endl;
  
  
  // get the combinations
  //cout << "Let me get some combinations for you." << endl;
  int n_elements = n_sources-1;
  vector< vector<int> > temp_combo_vecvec;
  
  if (n_sources < n_in_each_combo)
  {
    cout << "Not enough channels in this basin! I am aborting, the rest of your code might break." << endl;
    return;
  }
  else
  {
    // this gets the combinations
    bool zero_indexed = false;
    temp_combo_vecvec = combinations(n_elements, n_in_each_combo, zero_indexed);
  }


  combo_vecvec = temp_combo_vecvec;
 // Done
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment, It saves the source keys involved in each combination to get an idea of spatial distributions - B.G 2019
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDChiTools::opti_collinearity_by_basin_disorder_with_uncert_retain_key(LSDFlowInfo& FlowInfo,
          int baselevel_key, map<int,int>& sources_are_keys, map<int,int>& comboindex_are_keys, vector<float>& this_basin_chi,
          vector<float>& this_basin_elevation, vector<int>& this_basin_source, vector<int>& this_basin_river_node, vector< vector<int> >& combo_vecvec, vector<int>& n_pix_comb)
{

  // get the number of combinations
  int n_combinations = int(combo_vecvec.size());

  
  // Initilising the output and global related vectors
  vector<float> disorder_stat_vec(combo_vecvec.size()); 
  
  
  // Checking that the requested basin does exists
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }
  
  // Getting number of node to check some boundaries
  int n_nodes = int(node_sequence.size());
  // Getting the number of source keys involved  
  int n_sources = int(comboindex_are_keys.size());

  // get the combinations
  //cout << "Let me get some combinations for you." << endl;
  int n_elements = n_sources-1;
  // int n_in_each_combo = 3;
  
  // Preprocessing river size here, not forcing it not avoid reprocessing
  this->generate_river_size(false);

  // I am adapting the function from another one, During this process I keep so artifacts of the original structure that I may forget to remove. Ignore and forgive. YOLO.
  if(true)
  {
    size_t id_disorders = 0;
    

    // Going through the combination of all sources
    for (int combo = 0; combo<n_combinations; combo++)
    {


      // THis hosts the current combination of sources
      vector<int> these_combos = combo_vecvec[combo];
      // temp vector to add the next one
      vector<int> these_combo_sources;
      
      // add the trunk channel source
      these_combo_sources.push_back(comboindex_are_keys[0]); 
      
      // std::cout << "bk:" <<baselevel_key <<" ||";

      // THis gets the combination right into        
      for (int source_key = 0 ; source_key < int(these_combos.size()); source_key++)
      {
        these_combo_sources.push_back( comboindex_are_keys[ these_combos[source_key] ] );
      }


      

      // begin = std::chrono::steady_clock::now();

      // predefining the size of the vectors
      size_t this_size = 0;
      for (int source_key = 0 ; source_key < int(these_combo_sources.size()); source_key++)
      {
        this_size += size_t(source_key_to_size[these_combo_sources[source_key]]);
      }


      // now we sample and sort these vectors
      vector<float> chi_combo(this_size);
      vector<float> elev_combo(this_size);

      // Getting the current number of pixels for normalisation
      n_pix_comb[id_disorders] = int(this_size);

      
      // THE BIG CHANGE FROM PREVIOUS ALGORITHMS IS THAT IT TAKES ALREADY SORTED CHI_COMBOS AND ELEV_COMBOS VECTORS
      // THEREFORE THERE IS NO NEED TO RESORT THEM AND IT SAVES A LOT OF TIME
      size_t ti = 0;
      for(int node = 0; node< int(this_basin_chi.size()); node++)
      {
        if(find(these_combo_sources.begin(), these_combo_sources.end(), this_basin_source[node]) != these_combo_sources.end()) 
        {
          chi_combo[ti] = this_basin_chi[node];
          elev_combo[ti] = this_basin_elevation[node];
          ti++;
        }
      }


      // now you have the elevation and chi for this combination, sort them
      // initiate the sorted vectors as references to the other ones for comparison purposes with previous algorithms
      vector<float>& chi_sorted = chi_combo;
      vector<float>& elev_sorted = elev_combo;

      // now calculate disorder
      float this_delta_chi = 0;
      float sum_delta_chi = 0;
      
      int n_nodes_this_basin = int(chi_sorted.size());
      float chi_max = chi_sorted[n_nodes_this_basin-1];
      float chi_min = chi_sorted[0];
      // Sorry Simon, the following cout slows down the code when many rivers are tested: terminal output overflow
      // cout << "My first guess of minimum chi is: " << chi_min << endl;
      for(int i = 0; i<n_nodes_this_basin-1; i++)
      {
        this_delta_chi = fabs(chi_sorted[i+1]-chi_sorted[i]);
        sum_delta_chi+=this_delta_chi;
        if(chi_sorted[i] > chi_max)
        {
          chi_max = chi_sorted[i];
        }
        if(chi_sorted[i] < chi_min)
        {
          chi_min = chi_sorted[i];
        }
      }
      if(chi_sorted[n_nodes_this_basin-1] > chi_max)
      {
        chi_max = chi_sorted[n_nodes_this_basin-1];
      }
      float chi_range = chi_max-chi_min;
        
      float disorder_stat = (sum_delta_chi - chi_range)/chi_range;
      disorder_stat_vec[id_disorders] = disorder_stat;
      id_disorders++;
      // std::cout << "distat is " << chi_sorted.size()  << std::endl;

    }


  }
  
  return disorder_stat_vec;

}


void LSDChiTools::generate_river_size( bool force)
{
  // Checking if this has already been generated or if you want to force regeneration
  if(source_key_to_size.size()>0 && force == false)
    return;
  
  map<int,int> this_riv_key;

  for(size_t i=0; i<node_sequence.size(); i++)
  {
    int this_node = node_sequence[i];
    // cout << this_node << endl;
    int SK = source_keys_map[this_node];
    if(this_riv_key.count(SK) == 0)
      this_riv_key[SK] = 1;
    else
      this_riv_key[SK] = this_riv_key[SK] + 1;
  }


  source_key_to_size = this_riv_key;

  return;
  // Done

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment, It saves the source keys involved in each combination to get an idea of spatial distributions - B.G 2019
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map< vector<int>, vector<float> > LSDChiTools::TEST_FUNCTION_OLD_DISORDER_DO_NOT_USE_KEY(LSDFlowInfo& FlowInfo,
                                                 int baselevel_key)
{
  vector<float> disorder_stat_vec;
  vector<vector<int> > disorder_stat_vec_associated_SK;
  map< vector<int>, vector<float> > output;
  
  //cout << "Testing the segment collinearity for basin key " << baselevel_key << endl;
  // get some information about the number of basins
  int n_basins = int(ordered_baselevel_nodes.size());
  if (baselevel_key >= n_basins)
  {
    cout << "Fatal error LSDChiTools::test_all_segment_collinearity_by_basin_disorder" <<endl;
    cout << "You have selected a basin that doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  // this map will hold references to the sources. There are two of them
  // because the first has the sources as the keys and the second has the
  // combination index as the keys. 
  map<int,int> sources_are_keys;
  map<int,int> comboindex_are_keys;
  
  vector<float> this_basin_chi;
  vector<float> this_basin_elevation;
  vector<int> this_basin_source;
  int n_nodes = int(node_sequence.size());
  int this_node;
  int comboindex = 0;     // this is used to store an index into the combinations
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    if (baselevel_keys_map[this_node] == baselevel_key)
    {

      this_basin_chi.push_back(chi_data_map[this_node]);
      this_basin_elevation.push_back(elev_data_map[this_node]);
      this_basin_source.push_back(source_keys_map[this_node]);
      
      // if the key doesn't exist, add a source key counter
      if ( sources_are_keys.find( source_keys_map[this_node] ) == sources_are_keys.end() )
      {
        sources_are_keys[ source_keys_map[this_node] ] = comboindex;
        comboindex++;
      }
    }
  }
  
  // now do the second map by inverting the first map
  for(map<int,int>::iterator iter =sources_are_keys.begin(); iter != sources_are_keys.end(); ++iter)
  {
    int k =  iter->first;
    int v = iter->second;
    comboindex_are_keys[v] = k;
  }
  
  int n_sources = int(comboindex_are_keys.size());
  //cout << "The number of channels are: " << n_sources << endl;
  //cout << "The source node is: " << comboindex_are_keys[0] << endl;
  
  
  // get the combinations
  //cout << "Let me get some combinations for you." << endl;
  int n_elements = n_sources-1;
  int n_in_each_combo = 3;
  
  if (n_sources < n_in_each_combo)
  {
    cout << "Not enough channels in this basin! I am returning a nodata vector." << endl;
    disorder_stat_vec.push_back(-9999);
    disorder_stat_vec_associated_SK.push_back({-9999});
  }
  else
  {
    // this gets the combinations
    bool zero_indexed = false;
    vector< vector<int> > combo_vecvec = combinations(n_elements, n_in_each_combo, zero_indexed);
    
    bool print_combinations = false;
    if(print_combinations)
    {
      for (int i = 0; i< int(combo_vecvec.size()); i++)
      {
        for (int j = 0; j< int(combo_vecvec[i].size()); j++)
        {
          cout << combo_vecvec[i][j] << " ";
        }
        cout << endl;
      }
    }
    int n_combinations = int(combo_vecvec.size());
  
    // now you enter the combinations loop
    //int n_data_points = int(this_basin_source.size());
    //cout << "Checking the combinations" << endl;
    for (int combo = 0; combo<n_combinations; combo++)
    {
      vector<int> these_combos = combo_vecvec[combo];
      vector<int> these_combo_sources;
      
      // add the trunk channel source
      these_combo_sources.push_back(comboindex_are_keys[0]); 

      //cout << 0 <<":" << comboindex_are_keys[0] << " ";
      
      for (int source_key = 0 ; source_key < int(these_combos.size()); source_key++)
      {
        these_combo_sources.push_back( comboindex_are_keys[ these_combos[source_key] ] );
        //cout << these_combos[source_key] <<":" << comboindex_are_keys[ these_combos[source_key] ] << " ";
      }
      //cout << endl;
      disorder_stat_vec_associated_SK.push_back(these_combo_sources);
            
      
      // now we sample and sort these vectors
      vector<float> chi_combo;
      vector<float> elev_combo;
      vector<float> source_combo;
      
      for(int node = 0; node< int(this_basin_chi.size()); node++)
      {
        if(find(these_combo_sources.begin(), these_combo_sources.end(), this_basin_source[node]) != these_combo_sources.end()) 
        {
          chi_combo.push_back(this_basin_chi[node]);
          elev_combo.push_back(this_basin_elevation[node]);
        }
      }
        
      // now you have the elevation and chi for this combination, sort them
      // initiate the sorted vectors
      vector<float> chi_sorted;
      vector<float> elev_sorted;
      vector<size_t> index_map;    
  
      // sort the vectors
      matlab_float_sort(elev_combo, elev_sorted, index_map);
      matlab_float_reorder(chi_combo, index_map, chi_sorted);
        
      // now calculate disorder
      float this_delta_chi = 0;
      float sum_delta_chi = 0;
      
      int n_nodes_this_basin = int(chi_sorted.size());
      float corrected_chi_max = chi_sorted[n_nodes_this_basin-1];
      float corrected_chi_min = chi_sorted[0];
      float chi_min = 10000;
      float chi_max = 0;
      // if(chi_min != corrected_chi_min)
      // {
      //   cout << "BUGDETECTED:: CHI MIN=" << corrected_chi_min <<  " BUT CHI_MIN_USED="  << chi_min << endl;
      // }

      for(int i = 0; i<n_nodes_this_basin-1; i++)
      {
        this_delta_chi = fabs(chi_sorted[i+1]-chi_sorted[i]);
        sum_delta_chi+=this_delta_chi;
        if(chi_sorted[i] > chi_max)
        {
          chi_max = chi_sorted[i];
        }
        if(chi_sorted[i] < chi_min)
        {
          chi_min = chi_sorted[i];
        }
      }
      if(chi_sorted[n_nodes_this_basin-1] > chi_max)
      {
        chi_max = chi_sorted[n_nodes_this_basin-1];
      }
      float chi_range = chi_max-chi_min;
        
      float disorder_stat = (sum_delta_chi - chi_range)/chi_range;
      disorder_stat_vec.push_back(disorder_stat);
    }
  }


  // for debugging
  bool print_results = false;
  if(print_results)
  {
    cout << "Disorder stats for all the combinations are: " << endl;
    for(int i = 0; i< int(disorder_stat_vec.size()); i++)
    {
      cout << disorder_stat_vec[i] << " ";
    }
    cout << endl;
    cout << "I'm returning the disorder stat vector." << endl;
  }
  
  for(size_t i=0; i<disorder_stat_vec_associated_SK.size();i++)
  {
    if(output.find( disorder_stat_vec_associated_SK[i] ) == output.end())
    {
      output[ disorder_stat_vec_associated_SK[i] ] = {disorder_stat_vec[i]};
    }
    else
    {
      output[ disorder_stat_vec_associated_SK[i] ].push_back(disorder_stat_vec[i]);
    }
  }

  return output;

}






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_median_residuals(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix)
{
  vector< vector<float> > median_vecvec;
  vector< vector<float> > Q1_vecvec;
  vector< vector<float> > Q3_vecvec;
  vector<int> outlet_jns;

  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_median_residuals" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;
  //float thresh_area_for_chi = 0;      // This just gets chi from all pixels.

  string filename_bstats = file_prefix+"_movern_residuals_median.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  string filename_bstats_Q1 = file_prefix+"_movern_residuals_Q1.csv";
  ofstream stats_by_basin_out_Q1;
  stats_by_basin_out_Q1.open(filename_bstats_Q1.c_str());

  string filename_bstats_Q3 = file_prefix+"_movern_residuals_Q3.csv";
  ofstream stats_by_basin_out_Q3;
  stats_by_basin_out_Q3.open(filename_bstats_Q3.c_str());

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // calculate chi
    update_chi_data_map(FlowInfo, A_0, movern[i]);

    // The stats for the residuals
    vector<float> median_values, Q1_values, Q3_values;

    vector<int> all_basin_keys;

    // now run the collinearity test
    vector<float> these_residuals;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      these_residuals = retrieve_all_residuals_by_basin(FlowInfo, only_use_mainstem_as_reference,
                                                        basin_key);
      // BUg HERE WHEN FEEDING VOID VECTOR - TO FIX - Boris                                                  
      vector<float> these_stats = calculate_descriptive_stats(these_residuals);
      Q1_values.push_back(these_stats[1]);
      median_values.push_back(these_stats[2]);
      Q3_values.push_back(these_stats[3]);

      cout << "basin: " << basin_key << " and median residual is: " << these_stats[2] << endl;
    }

    // add the data to the vecvecs

    median_vecvec.push_back(median_values);
    Q1_vecvec.push_back(Q1_values);
    Q3_vecvec.push_back(Q3_values);
  }


  // Print the median data
  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",m_over_n="<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<median_vecvec[i][basin_key];
    }
    stats_by_basin_out << endl;
  }
  stats_by_basin_out.close();


  // Print the Q1 data
  stats_by_basin_out_Q1 << "basin_key,outlet_jn";
  stats_by_basin_out_Q1.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out_Q1 << ",m_over_n="<<movern[i];
  }
  stats_by_basin_out_Q1 << endl;
  stats_by_basin_out_Q1.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out_Q1 << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out_Q1 << "," <<Q1_vecvec[i][basin_key];
    }
    stats_by_basin_out_Q1 << endl;
  }
  stats_by_basin_out_Q1.close();


  // Print the Q3 data
  stats_by_basin_out_Q3 << "basin_key,outlet_jn";
  stats_by_basin_out_Q3.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out_Q3 << ",m_over_n="<<movern[i];
  }
  stats_by_basin_out_Q3 << endl;
  stats_by_basin_out_Q3.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out_Q3 << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out_Q3 << "," <<Q3_vecvec[i][basin_key];
    }
    stats_by_basin_out_Q3 << endl;
  }
  stats_by_basin_out_Q3.close();

}











//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma)
{
  // these vectors store all the values which are then used for printing
  vector< vector<float> > RMSE_vecvec;
  vector< vector<float> > MLE_vecvec;
  vector< vector<float> > total_MLE_vecvec;
  vector<int> reference_keys;
  vector<int> test_keys;
  vector<int> outlet_jns;

  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;
  //float thresh_area_for_chi = 0;      // This just gets chi from all pixels.

  string filename_bstats = file_prefix+"_basinstats.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_fullstats.csv";
    ofstream movern_stats_out;
    movern_stats_out.open(filename_fullstats.c_str());

    // calculate chi
    update_chi_data_map(FlowInfo, A_0, movern[i]);

    // these are the vectors that will hold the information about the
    // comparison between channels.
    // the _all vectors are one of all the basins
    // reference source is the source key of the reference channel
    vector<int> reference_source, all_reference_source;
    // test source is the source key of the test channel
    vector<int> test_source, all_test_source;
    // MLE the maximum liklihood estimator
    vector<float> MLE_values, all_MLE_values;
    // RMSE is the root mean square error
    vector<float> RMSE_values, all_RMSE_values;

    vector<float> tot_MLE_vec;
    // basin keys
    vector<int> all_basin_keys;

    // now run the collinearity test
    float tot_MLE;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      tot_MLE = test_all_segment_collinearity_by_basin(FlowInfo, only_use_mainstem_as_reference,
                                  basin_key,
                                  reference_source, test_source, MLE_values, RMSE_values, sigma);
      // concatenate the vectors to the "all" vectors
      all_reference_source.insert(all_reference_source.end(), reference_source.begin(), reference_source.end() );
      all_test_source.insert(all_test_source.end(), test_source.begin(), test_source.end() );
      all_MLE_values.insert(all_MLE_values.end(), MLE_values.begin(), MLE_values.end() );
      all_RMSE_values.insert(all_RMSE_values.end(), RMSE_values.begin(), RMSE_values.end() );
      all_basin_keys.insert(all_basin_keys.end(), reference_source.size(), basin_key);

      tot_MLE_vec.push_back(tot_MLE);
      cout << "basin: " << basin_key << " and tot_MLE: " << tot_MLE << endl;
    }

    // add the data to the vecvecs
    MLE_vecvec.push_back(all_MLE_values);
    RMSE_vecvec.push_back(all_RMSE_values);
    total_MLE_vecvec.push_back(tot_MLE_vec);
    reference_keys = all_reference_source;
    test_keys = all_test_source;

    // now print the data to the file
    movern_stats_out << "basin_key,reference_source_key,test_source_key,MLE,RMSE" << endl;
    int n_rmse_vals = int(all_RMSE_values.size());
    for(int i = 0; i<n_rmse_vals; i++)
    {
      movern_stats_out << all_basin_keys[i] << ","
                       << all_reference_source[i] << ","
                       << all_test_source[i] << ","
                       << all_MLE_values[i] << ","
                       << all_RMSE_values[i] << endl;
    }
    movern_stats_out.close();

  }

  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",m_over_n = "<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<total_MLE_vecvec[i][basin_key];
    }
    stats_by_basin_out << endl;
  }

  stats_by_basin_out.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge(LSDFlowInfo& FlowInfo,
                        LSDJunctionNetwork& JN, float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix,
                        LSDRaster& Discharge, float sigma)
{
  // these vectors store all the values which are then used for printing
  vector< vector<float> > RMSE_vecvec;
  vector< vector<float> > MLE_vecvec;
  vector< vector<float> > total_MLE_vecvec;
  vector<int> reference_keys;
  vector<int> test_keys;
  vector<int> outlet_jns;


  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;
  //float thresh_area_for_chi = 0;      // This just gets chi from all pixels.

  string filename_bstats = file_prefix+"_basinstats.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_fullstats.csv";
    ofstream movern_stats_out;
    movern_stats_out.open(filename_fullstats.c_str());

    // calculate chi
    float area_threshold = 0;

    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    // these are the vectors that will hold the information about the
    // comparison between channels.
    // the _all vectors are one of all the basins
    // reference source is the source key of the reference channel
    vector<int> reference_source, all_reference_source;
    // test source is the source key of the test channel
    vector<int> test_source, all_test_source;
    // MLE the maximum liklihood estimator
    vector<float> MLE_values, all_MLE_values;
    // RMSE is the root mean square error
    vector<float> RMSE_values, all_RMSE_values;

    vector<float> tot_MLE_vec;
    // basin keys
    vector<int> all_basin_keys;
    // now run the collinearity test
    float tot_MLE;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      tot_MLE = test_all_segment_collinearity_by_basin(FlowInfo, only_use_mainstem_as_reference,
                                  basin_key,
                                  reference_source, test_source, MLE_values, RMSE_values, sigma);
      // concatenate the vectors to the "all" vectors
      all_reference_source.insert(all_reference_source.end(), reference_source.begin(), reference_source.end() );
      all_test_source.insert(all_test_source.end(), test_source.begin(), test_source.end() );
      all_MLE_values.insert(all_MLE_values.end(), MLE_values.begin(), MLE_values.end() );
      all_RMSE_values.insert(all_RMSE_values.end(), RMSE_values.begin(), RMSE_values.end() );
      all_basin_keys.insert(all_basin_keys.end(), reference_source.size(), basin_key);
        
      tot_MLE_vec.push_back(tot_MLE);
      cout << "basin: " << basin_key << " and tot_MLE: " << tot_MLE << endl;
    }

    // add the data to the vecvecs
    MLE_vecvec.push_back(all_MLE_values);
    RMSE_vecvec.push_back(all_RMSE_values);
    total_MLE_vecvec.push_back(tot_MLE_vec);
    reference_keys = all_reference_source;
    test_keys = all_test_source;

    // now print the data to the file
    movern_stats_out << "basin_key,reference_source_key,test_source_key,MLE,RMSE" << endl;
    int n_rmse_vals = int(all_RMSE_values.size());
    for(int i = 0; i<n_rmse_vals; i++)
    {
      movern_stats_out << all_basin_keys[i] << ","
                       << all_reference_source[i] << ","
                       << all_test_source[i] << ","
                       << all_MLE_values[i] << ","
                       << all_RMSE_values[i] << endl;
    }
    movern_stats_out.close();

  }

  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",m_over_n = "<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<total_MLE_vecvec[i][basin_key];
    }
    stats_by_basin_out << endl;
  }

  stats_by_basin_out.close();

}

















//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_points(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma,
                        vector<float> chi_fractions_vector)
{
  // these vectors store all the values which are then used for printing
  vector< vector<float> > RMSE_vecvec;
  vector< vector<float> > MLE_vecvec;
  vector< vector<float> > total_MLE_vecvec;
  vector<int> reference_keys;
  vector<int> test_keys;
  vector<int> outlet_jns;

  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;
  //float thresh_area_for_chi = 0;      // This just gets chi from all pixels.

  string filename_bstats = file_prefix+"_basinstats.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_fullstats.csv";
    ofstream movern_stats_out;
    movern_stats_out.open(filename_fullstats.c_str());

    // calculate chi
    update_chi_data_map(FlowInfo, A_0, movern[i]);

    // these are the vectors that will hold the information about the
    // comparison between channels.
    // the _all vectors are one of all the basins
    // reference source is the source key of the reference channel
    vector<int> reference_source, all_reference_source;
    // test source is the source key of the test channel
    vector<int> test_source, all_test_source;
    // MLE the maximum liklihood estimator
    vector<float> MLE_values, all_MLE_values;
    // RMSE is the root mean square error
    vector<float> RMSE_values, all_RMSE_values;

    vector<float> tot_MLE_vec;
    // basin keys
    vector<int> all_basin_keys;

    // now run the collinearity test
    float tot_MLE;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      tot_MLE = test_all_segment_collinearity_by_basin_using_points(FlowInfo, only_use_mainstem_as_reference,
                                  basin_key,
                                  reference_source, test_source, MLE_values, RMSE_values, sigma,
                                  chi_fractions_vector);
      // concatenate the vectors to the "all" vectors
      all_reference_source.insert(all_reference_source.end(), reference_source.begin(), reference_source.end() );
      all_test_source.insert(all_test_source.end(), test_source.begin(), test_source.end() );
      all_MLE_values.insert(all_MLE_values.end(), MLE_values.begin(), MLE_values.end() );
      all_RMSE_values.insert(all_RMSE_values.end(), RMSE_values.begin(), RMSE_values.end() );
      all_basin_keys.insert(all_basin_keys.end(), reference_source.size(), basin_key);

      tot_MLE_vec.push_back(tot_MLE);
      cout << "basin: " << basin_key << " and tot_MLE: " << tot_MLE << endl;
    }

    // add the data to the vecvecs
    MLE_vecvec.push_back(all_MLE_values);
    RMSE_vecvec.push_back(all_RMSE_values);
    total_MLE_vecvec.push_back(tot_MLE_vec);
    reference_keys = all_reference_source;
    test_keys = all_test_source;

    // now print the data to the file
    movern_stats_out << "basin_key,reference_source_key,test_source_key,MLE,RMSE" << endl;
    int n_rmse_vals = int(all_RMSE_values.size());
    for(int i = 0; i<n_rmse_vals; i++)
    {
      movern_stats_out << all_basin_keys[i] << ","
                       << all_reference_source[i] << ","
                       << all_test_source[i] << ","
                       << all_MLE_values[i] << ","
                       << all_RMSE_values[i] << endl;
    }
    movern_stats_out.close();

  }

  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",m_over_n = "<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<total_MLE_vecvec[i][basin_key];
    }
    stats_by_basin_out << endl;
  }

  stats_by_basin_out.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_points(LSDFlowInfo& FlowInfo,
                        LSDJunctionNetwork& JN, float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix,
                        LSDRaster& Discharge, float sigma,
                        vector<float> chi_fractions_vector)
{
  // these vectors store all the values which are then used for printing
  vector< vector<float> > RMSE_vecvec;
  vector< vector<float> > MLE_vecvec;
  vector< vector<float> > total_MLE_vecvec;
  vector<int> reference_keys;
  vector<int> test_keys;
  vector<int> outlet_jns;


  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;
  //float thresh_area_for_chi = 0;      // This just gets chi from all pixels.

  string filename_bstats = file_prefix+"_basinstats.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_fullstats.csv";
    ofstream movern_stats_out;
    movern_stats_out.open(filename_fullstats.c_str());

    // calculate chi
    float area_threshold = 0;

    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    // these are the vectors that will hold the information about the
    // comparison between channels.
    // the _all vectors are one of all the basins
    // reference source is the source key of the reference channel
    vector<int> reference_source, all_reference_source;
    // test source is the source key of the test channel
    vector<int> test_source, all_test_source;
    // MLE the maximum liklihood estimator
    vector<float> MLE_values, all_MLE_values;
    // RMSE is the root mean square error
    vector<float> RMSE_values, all_RMSE_values;

    vector<float> tot_MLE_vec;

    // now run the collinearity test
    float tot_MLE;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      tot_MLE = test_all_segment_collinearity_by_basin_using_points(FlowInfo, only_use_mainstem_as_reference,
                                  basin_key,
                                  reference_source, test_source, MLE_values, RMSE_values, sigma,
                                  chi_fractions_vector);
      // concatenate the vectors to the "all" vectors
      all_reference_source.insert(all_reference_source.end(), reference_source.begin(), reference_source.end() );
      all_test_source.insert(all_test_source.end(), test_source.begin(), test_source.end() );
      all_MLE_values.insert(all_MLE_values.end(), MLE_values.begin(), MLE_values.end() );
      all_RMSE_values.insert(all_RMSE_values.end(), RMSE_values.begin(), RMSE_values.end() );

      tot_MLE_vec.push_back(tot_MLE);
      cout << "basin: " << basin_key << " and tot_MLE: " << tot_MLE << endl;
    }

    // add the data to the vecvecs
    MLE_vecvec.push_back(all_MLE_values);
    RMSE_vecvec.push_back(all_RMSE_values);
    total_MLE_vecvec.push_back(tot_MLE_vec);
    reference_keys = all_reference_source;
    test_keys = all_test_source;

    // now print the data to the file
    movern_stats_out << "reference_source_key,test_source_key,MLE,RMSE" << endl;
    int n_rmse_vals = int(all_RMSE_values.size());
    for(int i = 0; i<n_rmse_vals; i++)
    {
      movern_stats_out << all_reference_source[i] << ","
                       << all_test_source[i] << ","
                       << all_MLE_values[i] << ","
                       << all_RMSE_values[i] << endl;
    }
    movern_stats_out.close();

  }

  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",m_over_n = "<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<total_MLE_vecvec[i][basin_key];
    }
    stats_by_basin_out << endl;
  }

  stats_by_basin_out.close();

}






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_points_MC(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma,
                        int n_fracs,
                        int MC_iterations,
                        float max_frac)
{
  // these vectors store all the values which are then used for printing
  vector< vector<float> > min_MLE, first_quartile_MLE, median_MLE, third_quartile_MLE, max_MLE;

  // open the file that contains the basin stats
  string filename_bstats = file_prefix+"_points_MC_basinstats.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());


  // get the see for the random number generator
  long seed = time(NULL);

  vector< vector<float> > RMSE_vecvec;
  vector< vector<float> > MLE_vecvec;
  vector< vector<float> > total_MLE_vecvec;
  vector<int> reference_keys;
  vector<int> test_keys;
  vector<int> outlet_jns;

  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // calculate chi
    update_chi_data_map(FlowInfo, A_0, movern[i]);

    // get some vecvecs for storing information about the sources, MLEs, etc
    // for each iteration
    vector<int> master_reference_sources, master_test_sources, master_basin_keys;
    vector< vector<float> > total_MLE_iteration;


    for (int iteration = 0; iteration < MC_iterations; iteration++)
    {

      // these are the vectors that will hold the information about the
      // comparison between channels.
      // the _all vectors are one of all the basins
      // reference source is the source key of the reference channel
      vector<int> reference_source, all_reference_source;
      // test source is the source key of the test channel
      vector<int> test_source, all_test_source;
      // MLE the maximum liklihood estimator
      vector<float> MLE_values, all_MLE_values;
      // RMSE is the root mean square error
      vector<float> RMSE_values, all_RMSE_values;

      vector<float> tot_MLE_vec;
      // basin keys
      vector<int> all_basin_keys;

      // now run the collinearity test
      float tot_MLE;

      // here we need the logic for the chi_fractions
      float partition_size = max_frac/float(n_fracs);
      vector<float> chi_fractions_vector;
      for(int partition = 0; partition < n_fracs; partition++)
      {
        float this_part = ran3(&seed)*partition_size;
        //cout << "This part is: " << this_part << endl;

        chi_fractions_vector.push_back( float(partition)*partition_size + this_part );
        //cout << "Partion["<< partition << "]: " << chi_fractions_vector[partition] << endl;
      }


      for(int basin_key = 0; basin_key<n_basins; basin_key++)
      {
        //cout << "The sigma is: " << sigma << endl;
        tot_MLE = test_all_segment_collinearity_by_basin_using_points(FlowInfo, only_use_mainstem_as_reference,
                                    basin_key,
                                    reference_source, test_source, MLE_values, RMSE_values, sigma,
                                    chi_fractions_vector);
        // concatenate the vectors to the "all" vectors
        all_reference_source.insert(all_reference_source.end(), reference_source.begin(), reference_source.end() );
        all_test_source.insert(all_test_source.end(), test_source.begin(), test_source.end() );
        all_MLE_values.insert(all_MLE_values.end(), MLE_values.begin(), MLE_values.end() );
        all_RMSE_values.insert(all_RMSE_values.end(), RMSE_values.begin(), RMSE_values.end() );
        all_basin_keys.insert(all_basin_keys.end(), reference_source.size(), basin_key);

        tot_MLE_vec.push_back(tot_MLE);
        //cout << "basin: " << basin_key << " and tot_MLE: " << tot_MLE << endl;
      }

      // so we end up with vectors containing the MLE values for each source,
      // the RMSE, the reference and test sources, and the basin keys.

      // add the data to the vecvecs
      MLE_vecvec.push_back(all_MLE_values);
      RMSE_vecvec.push_back(all_RMSE_values);
      total_MLE_vecvec.push_back(tot_MLE_vec);
      reference_keys = all_reference_source;
      test_keys = all_test_source;

      // if this is the first iteration, push the reference, test, and basin keys vector
      // into the master vector
      if(iteration == 0)
      {
        master_reference_sources = all_reference_source;
        master_test_sources = all_test_source;
        master_basin_keys = master_basin_keys;
      }
      else
      {
        // check to see if the sources vectors are the same size
        if (master_reference_sources != all_reference_source)
        {
          cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_points_MC" << endl;
          cout << "We have a problem with the source vector, buddy!" << endl;
        }
      }

      // in this space we will need logic for storing all the MLE data and RMSE data
      // for each iteration and each basin. I'll save that work for later.


      // now get the total MLE vec for each iteration
      total_MLE_iteration.push_back(tot_MLE_vec);
    }

    // now for this m/n, we need to calculate the stats of the MLEs
    // this is a little slow since we need to reorganise the data in vecotrs
    // by repeatedly looping over the vecvec
    vector<float> this_movern_min_MLE;
    vector<float> this_movern_first_quartile_MLE;
    vector<float> this_movern_median_MLE;
    vector<float> this_movern_third_quartile_MLE;
    vector<float> this_movern_max_MLE;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      vector<float> this_basins_MLE;
      for(int iteration = 0; iteration < MC_iterations; iteration ++)
      {
        this_basins_MLE.push_back( total_MLE_iteration[iteration][basin_key] );
      }

      vector<float> descriptive_stats_MLE =  calculate_descriptive_stats(this_basins_MLE);

      this_movern_min_MLE.push_back(descriptive_stats_MLE[0]);
      this_movern_first_quartile_MLE.push_back(descriptive_stats_MLE[1]);
      this_movern_median_MLE.push_back(descriptive_stats_MLE[2]);
      this_movern_third_quartile_MLE.push_back(descriptive_stats_MLE[3]);
      this_movern_max_MLE.push_back(descriptive_stats_MLE[4]);
    }

    // push the data for this m/n into the vecvecs
    min_MLE.push_back(this_movern_min_MLE);
    first_quartile_MLE.push_back(this_movern_first_quartile_MLE);
    median_MLE.push_back(this_movern_median_MLE);
    third_quartile_MLE.push_back(this_movern_third_quartile_MLE);
    max_MLE.push_back(this_movern_max_MLE);

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_pointsMC.csv";
    ofstream movern_stats_out;
    movern_stats_out.open(filename_fullstats.c_str());

    // now print the data to the file
    movern_stats_out << "basin_number,minimum_MLE,first_quartile_MLE,median_MLE,third_quartile_MLE,maximum_MLE" << endl;
    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      movern_stats_out << basin_key << ","
                       << this_movern_min_MLE[basin_key] << ","
                       << this_movern_first_quartile_MLE[basin_key] << ","
                       << this_movern_median_MLE[basin_key] << ","
                       << this_movern_third_quartile_MLE[basin_key] << ","
                       << this_movern_max_MLE[basin_key] << endl;
    }
    movern_stats_out.close();
  }


  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",median_MLE_m_over_n="<<movern[i];
  }
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",FQ_MLE_m_over_n="<<movern[i];
  }
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",TQ_MLE_m_over_n="<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<median_MLE[i][basin_key];
    }
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<first_quartile_MLE[i][basin_key];
    }
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<third_quartile_MLE[i][basin_key];
    }
    stats_by_basin_out << endl;
  }

  stats_by_basin_out.close();

}













//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function test the collinearity of all segments compared to a reference
// segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_points_MC(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        bool only_use_mainstem_as_reference,
                        string file_prefix, float sigma,
                        int n_fracs,
                        int MC_iterations,
                        float max_frac, LSDRaster& Discharge)
{
  // these vectors store all the values which are then used for printing
  vector< vector<float> > min_MLE, first_quartile_MLE, median_MLE, third_quartile_MLE, max_MLE;

  // open the file that contains the basin stats
  string filename_bstats = file_prefix+"_points_MC_basinstats_Q.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  // get the see for the random number generator
  long seed = time(NULL);

  vector< vector<float> > RMSE_vecvec;
  vector< vector<float> > MLE_vecvec;
  vector< vector<float> > total_MLE_vecvec;
  vector<int> reference_keys;
  vector<int> test_keys;
  vector<int> outlet_jns;

  //cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern" << endl;
  //cout << "I am defaulting to A_0 = 1." << endl;
  vector<float> movern;
  float A_0 = 1;

  int n_basins = int(ordered_baselevel_nodes.size());

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }

  cout << endl << endl << "==========================" << endl;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    cout << "i: " << i << " and m over n: " << movern[i] << " ";

    // calculate chi
    float area_threshold = 0;
    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    // get some vecvecs for storing information about the sources, MLEs, etc
    // for each iteration
    vector<int> master_reference_sources, master_test_sources, master_basin_keys;
    vector< vector<float> > total_MLE_iteration;

    for (int iteration = 0; iteration < MC_iterations; iteration++)
    {

      // these are the vectors that will hold the information about the
      // comparison between channels.
      // the _all vectors are one of all the basins
      // reference source is the source key of the reference channel
      vector<int> reference_source, all_reference_source;
      // test source is the source key of the test channel
      vector<int> test_source, all_test_source;
      // MLE the maximum liklihood estimator
      vector<float> MLE_values, all_MLE_values;
      // RMSE is the root mean square error
      vector<float> RMSE_values, all_RMSE_values;

      vector<float> tot_MLE_vec;
      // basin keys
      vector<int> all_basin_keys;

      // now run the collinearity test
      float tot_MLE;

      // here we need the logic for the chi_fractions
      float partition_size = max_frac/float(n_fracs);
      vector<float> chi_fractions_vector;
      for(int partition = 0; partition < n_fracs; partition++)
      {
        float this_part = ran3(&seed)*partition_size;
        //cout << "This part is: " << this_part << endl;

        chi_fractions_vector.push_back( float(partition)*partition_size + this_part );
        //cout << "Partion["<< partition << "]: " << chi_fractions_vector[partition] << endl;
      }


      for(int basin_key = 0; basin_key<n_basins; basin_key++)
      {
        //cout << "The sigma is: " << sigma << endl;
        tot_MLE = test_all_segment_collinearity_by_basin_using_points(FlowInfo, only_use_mainstem_as_reference,
                                    basin_key,
                                    reference_source, test_source, MLE_values, RMSE_values, sigma,
                                    chi_fractions_vector);
        // concatenate the vectors to the "all" vectors
        all_reference_source.insert(all_reference_source.end(), reference_source.begin(), reference_source.end() );
        all_test_source.insert(all_test_source.end(), test_source.begin(), test_source.end() );
        all_MLE_values.insert(all_MLE_values.end(), MLE_values.begin(), MLE_values.end() );
        all_RMSE_values.insert(all_RMSE_values.end(), RMSE_values.begin(), RMSE_values.end() );
        all_basin_keys.insert(all_basin_keys.end(), reference_source.size(), basin_key);

        tot_MLE_vec.push_back(tot_MLE);
        //cout << "basin: " << basin_key << " and tot_MLE: " << tot_MLE << endl;
      }

      // so we end up with vectors containing the MLE values for each source,
      // the RMSE, the reference and test sources, and the basin keys.

      // add the data to the vecvecs
      MLE_vecvec.push_back(all_MLE_values);
      RMSE_vecvec.push_back(all_RMSE_values);
      total_MLE_vecvec.push_back(tot_MLE_vec);
      reference_keys = all_reference_source;
      test_keys = all_test_source;

      // if this is the first iteration, push the reference, test, and basin keys vector
      // into the master vector
      if(iteration == 0)
      {
        master_reference_sources = all_reference_source;
        master_test_sources = all_test_source;
        master_basin_keys = master_basin_keys;
      }
      else
      {
        // check to see if the sources vectors are the same size
        if (master_reference_sources != all_reference_source)
        {
          cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_points_MC" << endl;
          cout << "We have a problem with the source vector, buddy!" << endl;
        }
      }

      // in this space we will need logic for storing all the MLE data and RMSE data
      // for each iteration and each basin. I'll save that work for later.


      // now get the total MLE vec for each iteration
      total_MLE_iteration.push_back(tot_MLE_vec);
    }

    // now for this m/n, we need to calculate the stats of the MLEs
    // this is a little slow since we need to reorganise the data in vecotrs
    // by repeatedly looping over the vecvec
    vector<float> this_movern_min_MLE;
    vector<float> this_movern_first_quartile_MLE;
    vector<float> this_movern_median_MLE;
    vector<float> this_movern_third_quartile_MLE;
    vector<float> this_movern_max_MLE;

    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      vector<float> this_basins_MLE;
      for(int iteration = 0; iteration < MC_iterations; iteration ++)
      {
        this_basins_MLE.push_back( total_MLE_iteration[iteration][basin_key] );
      }

      vector<float> descriptive_stats_MLE =  calculate_descriptive_stats(this_basins_MLE);

      this_movern_min_MLE.push_back(descriptive_stats_MLE[0]);
      this_movern_first_quartile_MLE.push_back(descriptive_stats_MLE[1]);
      this_movern_median_MLE.push_back(descriptive_stats_MLE[2]);
      this_movern_third_quartile_MLE.push_back(descriptive_stats_MLE[3]);
      this_movern_max_MLE.push_back(descriptive_stats_MLE[4]);
    }

    // push the data for this m/n into the vecvecs
    min_MLE.push_back(this_movern_min_MLE);
    first_quartile_MLE.push_back(this_movern_first_quartile_MLE);
    median_MLE.push_back(this_movern_median_MLE);
    third_quartile_MLE.push_back(this_movern_third_quartile_MLE);
    max_MLE.push_back(this_movern_max_MLE);

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_pointsMC_Q.csv";
    ofstream movern_stats_out;
    movern_stats_out.open(filename_fullstats.c_str());

    // now print the data to the file
    movern_stats_out << "basin_number,minimum_MLE,first_quartile_MLE,median_MLE,third_quartile_MLE,maximum_MLE" << endl;
    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      movern_stats_out << basin_key << ","
                       << this_movern_min_MLE[basin_key] << ","
                       << this_movern_first_quartile_MLE[basin_key] << ","
                       << this_movern_median_MLE[basin_key] << ","
                       << this_movern_third_quartile_MLE[basin_key] << ","
                       << this_movern_max_MLE[basin_key] << endl;
    }
    movern_stats_out.close();
  }


  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",median_MLE_m_over_n="<<movern[i];
  }
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",FQ_MLE_m_over_n="<<movern[i];
  }
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",TQ_MLE_m_over_n="<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<median_MLE[i][basin_key];
    }
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<first_quartile_MLE[i][basin_key];
    }
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<third_quartile_MLE[i][basin_key];
    }
    stats_by_basin_out << endl;
  }

  stats_by_basin_out.close();

}









//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This functions test the goodness of fit for the m/n ratio using the 
// disorder method propsoed by Hergarten et al 2016
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_disorder(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        string file_prefix, bool use_uncert)
{
  cout << "I am now entering the disorder loop." << endl;

  int n_basins = int(ordered_baselevel_nodes.size());
  
  // some parameters to be stored in vectors
  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_disorder" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<int> outlet_jns;
  vector<float> movern;
  float A_0 = 1;  
  
  vector< vector<float> > disorder_vecvec;
  

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }
  
  
  // We do the full disorder statistic first
  cout << "I am calculating the disorder statistic!" << endl;
  
  vector<float> emptyvec;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    //cout << "i: " << i << " and m over n: " << movern[i] << " ";

    vector<float> these_disorders;

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_fullstats_disorder.csv";
    //ofstream movern_stats_out;
    //movern_stats_out.open(filename_fullstats.c_str());

    // calculate chi
    float area_threshold = 0;

    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                 area_threshold);
    update_chi_data_map(FlowInfo, this_chi);

    // these are the vectors that will hold the information about the disorder by basin
    vector<float> tot_MLE_vec;

    // now run the collinearity test
    float disorder_stat;

    //for(int basin_key = 0; basin_key<1; basin_key++)
    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      //disorder_stat = test_all_segment_collinearity_by_basin_using_disorder(FlowInfo, basin_key);
      disorder_stat = test_collinearity_by_basin_disorder(FlowInfo, basin_key);
      cout << "basin: " << basin_key << " and m/n is: " << movern[i] << " and disorder stat is: " << disorder_stat << endl;
      these_disorders.push_back(disorder_stat);
    }
    disorder_vecvec.push_back(these_disorders);
    
  }
  
  // now we need to loop through the vecvec getting the minimum disorder for each basin
  map<int,float> best_fit_movern_disorder_map;
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    vector<float> these_disorders;
    float min_disorder = 1000000000000; // a big number since the disorder needs to be smaller than this.  
    best_fit_movern_disorder_map[basin_key] = -9999; 
    for(int i = 0; i< n_movern; i++)
    {
      // check to see if this disorder is a minimum for this basin
      if (disorder_vecvec[i][basin_key] < min_disorder)
      {
        min_disorder = disorder_vecvec[i][basin_key];
        best_fit_movern_disorder_map[basin_key] = movern[i]; 
      }
    }
  }
  
  if(file_prefix != "EXPLICITELY_DO_NOT_SAVE_THE_OUTPUT")
  {
    // open the file that contains the basin stats
    string filename_bstats = file_prefix+"_disorder_basinstats.csv";
    ofstream stats_by_basin_out;
    stats_by_basin_out.open(filename_bstats.c_str());

    stats_by_basin_out << "basin_key,outlet_jn";
    stats_by_basin_out.precision(4);
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << ",m_over_n = "<<movern[i];
    }
    stats_by_basin_out << endl;
    stats_by_basin_out.precision(9);
    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
      for(int i = 0; i< n_movern; i++)
      {
        stats_by_basin_out << "," <<disorder_vecvec[i][basin_key];
      }
      stats_by_basin_out << endl;
    }
    stats_by_basin_out.close();
  }
  
  
  // Now if the use uncertainty flag is true, calculate the disorder statistics. 
  if (use_uncert)
  {
    // first initiate the vectors for holding the combination vecs.
    // the key is the basin number. Each element will hold the m over n value
    // with the lowest disorder
    map<int, vector<float> > best_fit_movern_for_basins;
    map<int, vector<float> > lowest_disorder_for_basins;

    for(int i = 0; i< n_movern; i++)
    {
      // get the m over n value
      float this_movern = float(i)*delta_movern+start_movern;
      movern.push_back( this_movern );
      cout << "i: " << i << " and m over n: " << movern[i] << endl;

      // calculate chi
      float area_threshold = 0;
  
      LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                   area_threshold);
      update_chi_data_map(FlowInfo, this_chi);
      
      // now loop through basins
      //for(int basin_key = 0; basin_key<1; basin_key++)
      for(int basin_key = 0; basin_key<n_basins; basin_key++)
      {
        cout << "Testing uncert the basin key is " << basin_key << endl;
        vector<float> disorder_stats = test_collinearity_by_basin_disorder_with_uncert(FlowInfo, basin_key);
        
        // if this is the first m over n value, then initiate the vectors for this basin key
        if (i == 0)
        {
          lowest_disorder_for_basins[basin_key] = disorder_stats;
          
          int n_combos_this_basin = int(disorder_stats.size());
          vector<float> best_fit_movern;
          for(int bf = 0; bf < n_combos_this_basin; bf++)
          {
            best_fit_movern.push_back( this_movern );
          }
          best_fit_movern_for_basins[basin_key] = best_fit_movern;
        }
        else
        {
          // loop through all the combos and get the best fit movern
          vector<float> existing_lowest_disorder = lowest_disorder_for_basins[basin_key];
          vector<float> existing_best_fit_movern = best_fit_movern_for_basins[basin_key];
          int n_combos_this_basin = int(disorder_stats.size());
          
          for(int bf = 0; bf < n_combos_this_basin; bf++)
          {
            if (existing_lowest_disorder[bf] > disorder_stats[bf] )
            {
              existing_lowest_disorder[bf] = disorder_stats[bf];
              existing_best_fit_movern[bf] = this_movern;
            }
          }
          lowest_disorder_for_basins[basin_key] = existing_lowest_disorder;
          best_fit_movern_for_basins[basin_key] = existing_best_fit_movern;
        }
      }  // end basin loop
    }    // end m/n loop
    
    if(file_prefix != "EXPLICITELY_DO_NOT_SAVE_THE_OUTPUT") // This is a condition in the rare case you don't want to save the file, just to calculate.
    {
      // open the outfile
      string filename_fullstats = file_prefix+"_fullstats_disorder_uncert.csv";
      ofstream stats_by_basin_out;
      stats_by_basin_out.open(filename_fullstats.c_str());
    
      stats_by_basin_out << "basin_key,N_combinations,minimum,first_quartile,median,third_quartile,maximum,mean,standard_deviation,standard_error,MAD, best_fit_for_all_tribs" << endl;
      stats_by_basin_out.precision(8);
      for(int basin_key = 0; basin_key<n_basins; basin_key++)
      {
        vector<float> these_movern = best_fit_movern_for_basins[basin_key];
        int n_combinations =  int(these_movern.size());
        
        vector<float> these_stats = calculate_descriptive_stats(these_movern);
        stats_by_basin_out << basin_key << ",";
        stats_by_basin_out << n_combinations << ",";
        stats_by_basin_out << these_stats[0] <<",";
        stats_by_basin_out << these_stats[1] <<",";
        stats_by_basin_out << these_stats[2] <<",";
        stats_by_basin_out << these_stats[3] <<",";
        stats_by_basin_out << these_stats[4] <<",";
        stats_by_basin_out << these_stats[5] <<",";
        stats_by_basin_out << these_stats[6] <<",";
        stats_by_basin_out << these_stats[7] <<",";
        stats_by_basin_out << these_stats[8] <<",";
        stats_by_basin_out << best_fit_movern_disorder_map[basin_key] << endl;
      }
      stats_by_basin_out.close();
    }
  }
}









//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This functions test the goodness of fit for the concavity using the 
// disorder method propsoed by Hergarten et al 2016
// Uses a discharge raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_with_discharge_using_disorder(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JN,
                        float start_movern, float delta_movern, int n_movern,
                        string file_prefix, bool use_uncert, LSDRaster& Discharge)
{
  cout << "I am now entering the disorder loop." << endl;

  int n_basins = int(ordered_baselevel_nodes.size());
  
  // some parameters to be stored in vectors
  cout << "LSDChiTools::calculate_goodness_of_fit_collinearity_fxn_movern_using_disorder" << endl;
  cout << "I am defaulting to A_0 = 1." << endl;
  vector<int> outlet_jns;
  vector<float> movern;
  float A_0 = 1;  
  
  vector< vector<float> > disorder_vecvec;
  

  // get the outlet junction of each basin key
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    int outlet_node = ordered_baselevel_nodes[basin_key];
    int outlet_jn = JN.get_Junction_of_Node(outlet_node, FlowInfo);
    outlet_jns.push_back(outlet_jn);
  }
  
  
  // We do the full disorder statistic first
  cout << "I am calculating the disorder statistic!" << endl;
  
  vector<float> emptyvec;
  for(int i = 0; i< n_movern; i++)
  {
    // get the m over n value
    movern.push_back( float(i)*delta_movern+start_movern );
    //cout << "i: " << i << " and m over n: " << movern[i] << " ";

    vector<float> these_disorders;

    // open the outfile
    string filename_fullstats = file_prefix+"_"+dtoa(movern[i])+"_fullstatsq_disorder.csv";
    //ofstream movern_stats_out;
    //movern_stats_out.open(filename_fullstats.c_str());

    // calculate chi
    float area_threshold = 0;
    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    // these are the vectors that will hold the information about the disorder by basin
    vector<float> tot_MLE_vec;

    // now run the collinearity test
    float disorder_stat;

    //for(int basin_key = 0; basin_key<1; basin_key++)
    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      //disorder_stat = test_all_segment_collinearity_by_basin_using_disorder(FlowInfo, basin_key);
      disorder_stat = test_collinearity_by_basin_disorder(FlowInfo, basin_key);
      cout << "basin: " << basin_key << " and m/n is: " << movern[i] << " and disorder stat is: " << disorder_stat << endl;
      these_disorders.push_back(disorder_stat);
    }
    disorder_vecvec.push_back(these_disorders);
    
  }
  
  // now we need to loop through the vecvec getting the minimum disorder for each basin
  map<int,float> best_fit_movern_disorder_map;
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    vector<float> these_disorders;
    float min_disorder = 1000000000000; // a big number since the disorder needs to be smaller than this.  
    best_fit_movern_disorder_map[basin_key] = -9999; 
    for(int i = 0; i< n_movern; i++)
    {
      // check to see if this disorder is a minimum for this basin
      if (disorder_vecvec[i][basin_key] < min_disorder)
      {
        min_disorder = disorder_vecvec[i][basin_key];
        best_fit_movern_disorder_map[basin_key] = movern[i]; 
      }
    }
  }
  

  // open the file that contains the basin stats
  string filename_bstats = file_prefix+"_disorder_basinstats.csv";
  ofstream stats_by_basin_out;
  stats_by_basin_out.open(filename_bstats.c_str());

  stats_by_basin_out << "basin_key,outlet_jn";
  stats_by_basin_out.precision(4);
  for(int i = 0; i< n_movern; i++)
  {
    stats_by_basin_out << ",m_over_n = "<<movern[i];
  }
  stats_by_basin_out << endl;
  stats_by_basin_out.precision(9);
  for(int basin_key = 0; basin_key<n_basins; basin_key++)
  {
    stats_by_basin_out << basin_key << "," << outlet_jns[basin_key];
    for(int i = 0; i< n_movern; i++)
    {
      stats_by_basin_out << "," <<disorder_vecvec[i][basin_key];
    }
    stats_by_basin_out << endl;
  }
  stats_by_basin_out.close();
  
  
  // Now if the use uncertainty flag is true, calculate the disorder statistics. 
  if (use_uncert)
  {
    // first initiate the vectors for holding the combination vecs.
    // the key is the basin number. Each element will hold the m over n value
    // with the lowest disorder
    map<int, vector<float> > best_fit_movern_for_basins;
    map<int, vector<float> > lowest_disorder_for_basins;

    for(int i = 0; i< n_movern; i++)
    {
      // get the m over n value
      float this_movern = float(i)*delta_movern+start_movern;
      movern.push_back( this_movern );
      cout << "i: " << i << " and m over n: " << movern[i] << endl;

      // calculate chi
      float area_threshold = 0;
  
      LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern[i], A_0,
                                   area_threshold);
      update_chi_data_map(FlowInfo, this_chi);
      
      // now loop through basins
      //for(int basin_key = 0; basin_key<1; basin_key++)
      for(int basin_key = 0; basin_key<n_basins; basin_key++)
      {
        cout << "Testing uncert the basin key is " << basin_key << endl;
        vector<float> disorder_stats = test_collinearity_by_basin_disorder_with_uncert(FlowInfo, basin_key);
        
        // if this is the first m over n value, then initiate the vectors for this basin key
        if (i == 0)
        {
          lowest_disorder_for_basins[basin_key] = disorder_stats;
          
          int n_combos_this_basin = int(disorder_stats.size());
          vector<float> best_fit_movern;
          for(int bf = 0; bf < n_combos_this_basin; bf++)
          {
            best_fit_movern.push_back( this_movern );
          }
          best_fit_movern_for_basins[basin_key] = best_fit_movern;
        }
        else
        {
          // loop through all the combos and get the best fit movern
          vector<float> existing_lowest_disorder = lowest_disorder_for_basins[basin_key];
          vector<float> existing_best_fit_movern = best_fit_movern_for_basins[basin_key];
          int n_combos_this_basin = int(disorder_stats.size());
          
          for(int bf = 0; bf < n_combos_this_basin; bf++)
          {
            if (existing_lowest_disorder[bf] > disorder_stats[bf] )
            {
              existing_lowest_disorder[bf] = disorder_stats[bf];
              existing_best_fit_movern[bf] = this_movern;
            }
          }
          lowest_disorder_for_basins[basin_key] = existing_lowest_disorder;
          best_fit_movern_for_basins[basin_key] = existing_best_fit_movern;
        }
      }  // end basin loop
    }    // end m/n loop
    

    // open the outfile
    string filename_fullstats = file_prefix+"_fullstats_disorder_uncert.csv";
    ofstream stats_by_basin_out;
    stats_by_basin_out.open(filename_fullstats.c_str());
  
    stats_by_basin_out << "basin_key,N_combinations,minimum,first_quartile,median,third_quartile,maximum,mean,standard_deviation,standard_error,MAD, best_fit_for_all_tribs" << endl;
    stats_by_basin_out.precision(8);
    for(int basin_key = 0; basin_key<n_basins; basin_key++)
    {
      vector<float> these_movern = best_fit_movern_for_basins[basin_key];
      int n_combinations =  int(these_movern.size());
      
      vector<float> these_stats = calculate_descriptive_stats(these_movern);
      stats_by_basin_out << basin_key << ",";
      stats_by_basin_out << n_combinations << ",";
      stats_by_basin_out << these_stats[0] <<",";
      stats_by_basin_out << these_stats[1] <<",";
      stats_by_basin_out << these_stats[2] <<",";
      stats_by_basin_out << these_stats[3] <<",";
      stats_by_basin_out << these_stats[4] <<",";
      stats_by_basin_out << these_stats[5] <<",";
      stats_by_basin_out << these_stats[6] <<",";
      stats_by_basin_out << these_stats[7] <<",";
      stats_by_basin_out << these_stats[8] <<",";
      stats_by_basin_out << best_fit_movern_disorder_map[basin_key] << endl;
    }
    stats_by_basin_out.close();
  }
}










//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This drives the m/n MCMC analysis
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::MCMC_driver(LSDFlowInfo& FlowInfo, int minimum_contributing_pixels, float sigma,
                                 float movern_minimum, float movern_maximum,
                                 int N_chain_links,
                                 string OUT_DIR, string OUT_ID,
                                 bool use_points)
{
  // so we loop through the basins
  float this_dmovern_stddev = 0.1;
  int n_basins = int(key_to_baselevel_map.size());
  for (int basin_key = 0; basin_key < n_basins; basin_key++)
  {
    cout << "Running MCMC on basin: " << basin_key << endl;

    float this_sigma = MCMC_for_movern_tune_sigma(FlowInfo, minimum_contributing_pixels,
                               this_dmovern_stddev,
                               movern_minimum, movern_maximum,
                               basin_key, use_points);

    //float this_dmovern_sigma = MCMC_for_movern_tune_dmovern(FlowInfo, minimum_contributing_pixels, sigma,
    //                             movern_minimum, movern_maximum, basin_key);
    //
    // now run the chain
    string chain_file = OUT_DIR+OUT_ID+"_Basin"+itoa(basin_key)+"_chain.csv";
    bool printChain = true;
    float accept = MCMC_for_movern(chain_file, printChain, FlowInfo,
                                 minimum_contributing_pixels,
                                 N_chain_links, this_sigma, this_dmovern_stddev,
                                 movern_minimum,movern_maximum,basin_key, use_points);
    cout << "The final acceptance rate was: " << accept << endl;
  }



}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This will tune the dmovern_stddev until the acceptance rate is ~0.25
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::MCMC_for_movern_tune_dmovern(LSDFlowInfo& FlowInfo, int minimum_contributing_pixels, float sigma,
                                 float movern_minimum, float movern_maximum,
                                 int basin_key, bool use_points)
{

  float min_acceptance_rate = 0.2;
  float max_acceptance_rate = 0.33;
  float this_acceptance_rate = 1.0;

  string ChainFname = "NULL";
  bool printChain = false;
  int NIterations = 1000;


  float range = movern_maximum-movern_minimum;
  float max_dmovern = range/3;
  int n_steps = 0;
  float this_dmovern_stddev = 0.2;
  while( n_steps < 10)
  {

    // NOTE: It would probably make more sense to use Newton-Raphson to get the
    // correct dmovern but we will use this stupid method for now.

    cout << "Looking at the acceptance rate. The current dmovern_stddev is: " << this_dmovern_stddev << endl;
    this_acceptance_rate = MCMC_for_movern(ChainFname, printChain, FlowInfo, minimum_contributing_pixels, NIterations, sigma, this_dmovern_stddev,
                     movern_minimum, movern_maximum, basin_key, use_points);

    cout << "The acceptance rate is: " << this_acceptance_rate << " and the m/n stddev is: " << this_dmovern_stddev << endl;

    if (this_acceptance_rate > max_acceptance_rate)
    {
      this_dmovern_stddev = this_dmovern_stddev*1.45;
    }
    else if (this_acceptance_rate < min_acceptance_rate)
    {
      this_dmovern_stddev = this_dmovern_stddev*0.77;
    }
    else
    {
      n_steps = 10;
    }
    // Limit the maximum stddec of dmovern to one third of the range in m/n
    if (this_dmovern_stddev > max_dmovern)
    {
      this_dmovern_stddev = max_dmovern;
      n_steps = 10;
    }

    // we don't want this to go on forever.
    n_steps++;
  }

  return this_dmovern_stddev;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This will tune the dmovern_stddev until the acceptance rate is ~0.25
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::MCMC_for_movern_tune_sigma(LSDFlowInfo& FlowInfo, int minimum_contributing_pixels,
                                              float dmovern_stddev,
                                              float movern_minimum, float movern_maximum,
                                              int basin_key, bool use_points)
{
  // This tries to get the correct acceptance rate by tuning sigma.
  // You give it a stddev
  // Remember: increasing sigma makes all the MLE values go up so increases
  // the acceptance rate.



  float min_acceptance_rate = 0.2;
  float max_acceptance_rate = 0.33;
  float this_acceptance_rate = 1.0;

  string ChainFname = "NULL";
  bool printChain = false;
  int NIterations = 2500;
  float this_sigma = 2000;

  if (use_points)
  {
    this_sigma = 100;
  }

  int n_steps = 0;

  while( n_steps < 20)
  {

    // NOTE: It would probably make more sense to use Newton-Raphson to get the
    // correct dmovern but we will use this stupid method for now.

    cout << "Looking at the acceptance rate. I've changed sigma to: " << this_sigma << endl;
    this_acceptance_rate = MCMC_for_movern(ChainFname, printChain, FlowInfo, minimum_contributing_pixels, NIterations, this_sigma, dmovern_stddev,
                     movern_minimum, movern_maximum, basin_key, use_points);

    cout << "The acceptance rate is: " << this_acceptance_rate << " and the m/n stddev is: " << this_sigma << endl;

    if (this_acceptance_rate > max_acceptance_rate)
    {
      this_sigma = this_sigma*0.77;    // these factors are arbitray: they just need to increase or decrease the acceptance rate.
    }
    else if (this_acceptance_rate < min_acceptance_rate)
    {
      this_sigma = this_sigma*1.45;
    }
    else
    {
      n_steps = 20;
    }

    // we don't want this to go on forever.
    n_steps++;
  }

  return this_sigma;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This runs an MCMC chain on the goodness of fit for a basin
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDChiTools::MCMC_for_movern(string ChainFname, bool printChain, LSDFlowInfo& FlowInfo,
                                 int minimum_contributing_pixels, int NIterations,
                                 float sigma, float dmovern_stddev,
                                 float movern_minimum, float movern_maximum,
                                 int basin_key, bool use_points)
{

  map<int,int> outlet_node_from_basin_key_map = get_outlet_node_from_basin_key_map();

  vector<float> chi_upslope_fracs;
  float start_frac = 0.4;
  float dfrac = 0.025;
  for(int i = 0; i< 11; i++)
  {
    chi_upslope_fracs.push_back(start_frac - float(i)*dfrac);
  }


  //Declarations
  float LastLikelihood;               //Last accepted likelihood
  float NewLikelihood;                //New likelihood
  float LikelihoodRatio;              //Ratio between last and new likelihoods
  float AcceptanceProbability;        //New iteration is accepted if likelihood ratio exceeds

  int NAccepted = 0;      //count accepted parameters
  int NRejected = 0;      //count rejected parameters

  // these are the vectors that will hold the information about the
  // comparison between channels.
  // the _all vectors are one of all the basins
  // reference source is the source key of the reference channel
  vector<int> reference_source, all_reference_source;
  // test source is the source key of the test channel
  vector<int> test_source, all_test_source;
  // MLE the maximum liklihood estimator
  vector<float> MLE_values, all_MLE_values;
  // RMSE is the root mean square error
  vector<float> RMSE_values, all_RMSE_values;

  // you need a seed for the random number generator
  long seed = time(NULL);

  // these are numbers for the change in   chi
  float gauss_mean = 0;
  float std_dev = dmovern_stddev;
  float gauss_minimum = -3.0*std_dev;   // the gaussian function implemented in statstools truncates at 3 standard deviations
  bool allowNegative = true;

  ofstream ChainFileOut(ChainFname.c_str());
  if(printChain)
  {
    ChainFileOut  << "i,movern_New,movern_Old,NewLikelihood,LastLikelihood,NAccepted,NRejected" << endl;
  }


  // Caluclate initial chi
  float A_0 = 1;
  float dmovern, movern_old, movern_new, reflect;
  movern_old = 0.5;
  movern_new = 0.5;
  update_chi_data_map(FlowInfo, A_0, movern_new);

  // get the initial MLE
  bool only_use_mainstem_as_reference = true;
  if (use_points)
  {
    LastLikelihood = test_all_segment_collinearity_by_basin_using_points(FlowInfo, only_use_mainstem_as_reference,
                                  basin_key,reference_source, test_source, MLE_values, RMSE_values, sigma,
                                  chi_upslope_fracs);
  }
  else
  {
    LastLikelihood = test_all_segment_collinearity_by_basin(FlowInfo, only_use_mainstem_as_reference,
                                  basin_key,reference_source, test_source, MLE_values, RMSE_values, sigma);
  }

  int print_interval = 100;

  // Now do the metropolis hastings algorithm
  for (int j = 0; j<NIterations; j++)
  {
    if (j%print_interval == 0)
    {
      cout << "Iteration " << j << " of " << NIterations << endl;
    }

    // Vary the movern value
    dmovern = getGaussianRandom(gauss_minimum, gauss_mean, allowNegative);
    movern_new = movern_old + dmovern;
    //cout << "dmovern is: " << dmovern << " and New m over n is: " <<  movern_new << endl;
    // reflect the data if necessary
    if ( movern_new < movern_minimum)
    {
      reflect = movern_minimum - movern_new;
      movern_new = reflect+movern_minimum;
    }
    if ( movern_new > movern_maximum)
    {
      reflect = movern_new- movern_maximum;
      movern_new = movern_maximum - reflect;
    }

    //cout << "After reflect New m over n is: " <<  movern_new << endl;

    // run the model with the new parameters
    update_chi_data_map_for_single_basin(FlowInfo, A_0, movern_new,
                                     minimum_contributing_pixels, basin_key,
                                     outlet_node_from_basin_key_map);
    //cout << "Got chi " << endl;
    if (use_points)
    {
      NewLikelihood = test_all_segment_collinearity_by_basin_using_points(FlowInfo, only_use_mainstem_as_reference,
                                    basin_key,reference_source, test_source, MLE_values, RMSE_values, sigma,
                                    chi_upslope_fracs);
    }
    else
    {
      NewLikelihood = test_all_segment_collinearity_by_basin(FlowInfo, only_use_mainstem_as_reference,
                                    basin_key,reference_source, test_source, MLE_values, RMSE_values, sigma);
    }
    // get the likelihood ratio
    LikelihoodRatio = NewLikelihood/LastLikelihood;

    //cout << "Ratio: " << LikelihoodRatio << "New MLE: " << NewLikelihood << " and old MLE: " << LastLikelihood << endl;

    // get the acceptance probability (this is set up so that occasional
    // guesses that are worse than the lst one get accepted so that
    // the chain can visit all of parameter space)
    AcceptanceProbability = ran3(&seed);

    // if accepted
    if (LikelihoodRatio > AcceptanceProbability)
    {
      NAccepted++;
    }
    else
    {
      NRejected++;
    }

    // Now print to the chain file
    if(printChain)
    {
      ChainFileOut  << j << "," << movern_new << "," << movern_old << ","
                    << NewLikelihood << "," << LastLikelihood << ","
                    << NAccepted << "," << NRejected << endl;

    // reset m over n old
    if (LikelihoodRatio > AcceptanceProbability)
    {
      LastLikelihood = NewLikelihood;
      movern_old = movern_new;
    }
    }
  }


  ChainFileOut.close();


  float accept_percent = float(NAccepted)/ float(NAccepted+NRejected);
  return accept_percent;



}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a series of simple profiles (chi-elevation) as a function of
// movern
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_profiles_as_fxn_movern(LSDFlowInfo& FlowInfo, string filename, float start_movern, float delta_movern, int n_movern)
{
  float A_0 = 1;
  float this_movern;

  vector<float> movern_values;
  vector< vector<float> > chi_vecvec;
  vector<float> empty_vec;
  vector<float> this_chi_vec;
  int this_node;
  int n_nodes = int(node_sequence.size());

  // loop through m over n values
  for(int i = 0; i< n_movern; i++)
  {

    this_movern =  float(i)*delta_movern+start_movern;
    update_chi_data_map(FlowInfo, A_0, this_movern);

    cout << "m/n is: " << this_movern << endl;

    movern_values.push_back(this_movern);
    this_chi_vec = empty_vec;

    // now get the chi values for each node and push them into the chi_vecvec
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      this_chi_vec.push_back(chi_data_map[this_node]);
    }
    chi_vecvec.push_back(this_chi_vec);
  }
  cout << "Okay, I've got all the chi values in the vecvec." << endl;
  // okay, we are done getting all the chi values, now add these into the file

  ofstream chi_csv_out;
  cout << "Running the printing for movern. Filename is: " << filename << endl;
  chi_csv_out.open(filename.c_str());
  chi_csv_out << "source_key,basin_key,elevation";
  for (int i = 0; i< n_movern; i++)
  {
    chi_csv_out << ",m_over_n = " << movern_values[i];
  }
  chi_csv_out << endl;

  // now loop through all the nodes
  chi_csv_out.precision(5);
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    chi_csv_out << source_keys_map[this_node] << ","
                 << baselevel_keys_map[this_node] << ","
                 << elev_data_map[this_node];

    for (int i = 0; i< n_movern; i++)
    {
      chi_csv_out << "," << chi_vecvec[i][n];
    }
    chi_csv_out << endl;
  }
  chi_csv_out.close();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a series of simple profiles (chi-elevation) as a function of
// movern
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_profiles_as_fxn_movern_with_burned_raster(LSDFlowInfo& FlowInfo,
                                         string filename, float start_movern,
                                         float delta_movern, int n_movern,
                                         LSDRaster& BurnRaster, string burned_column_name)
{
  float A_0 = 1;
  float this_movern;

  vector<float> movern_values;
  vector< vector<float> > chi_vecvec;
  vector<float> empty_vec;
  vector<float> this_chi_vec;
  int this_node;
  int n_nodes = int(node_sequence.size());

  // loop through m over n values
  for(int i = 0; i< n_movern; i++)
  {

    this_movern =  float(i)*delta_movern+start_movern;
    update_chi_data_map(FlowInfo, A_0, this_movern);

    cout << "m/n is: " << this_movern << endl;

    movern_values.push_back(this_movern);
    this_chi_vec = empty_vec;

    // now get the chi values for each node and push them into the chi_vecvec
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      this_chi_vec.push_back(chi_data_map[this_node]);
    }
    chi_vecvec.push_back(this_chi_vec);
  }
  cout << "Okay, I've got all the chi values in the vecvec." << endl;
  // okay, we are done getting all the chi values, now add these into the file

  ofstream chi_csv_out;
  cout << "Running the printing for movern. Filename is: " << filename << endl;
  chi_csv_out.open(filename.c_str());
  chi_csv_out << "source_key,basin_key,elevation";
  for (int i = 0; i< n_movern; i++)
  {
    chi_csv_out << ",m_over_n = " << movern_values[i];
  }
  chi_csv_out << "," << burned_column_name << endl;

  // now loop through all the nodes
  chi_csv_out.precision(5);
  int curr_row,curr_col;
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);

    chi_csv_out << source_keys_map[this_node] << ","
                 << baselevel_keys_map[this_node] << ","
                 << elev_data_map[this_node];

    for (int i = 0; i< n_movern; i++)
    {
      chi_csv_out << "," << chi_vecvec[i][n];
    }
    chi_csv_out << "," << BurnRaster.get_data_element(curr_row,curr_col);
    chi_csv_out << endl;
  }
  chi_csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a series of simple profiles (chi-elevation) as a function of
// movern
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_profiles_as_fxn_movern_with_secondary_raster(LSDFlowInfo& FlowInfo,
                                         string filename, float start_movern,
                                         float delta_movern, int n_movern,
                                         LSDRaster& BurnRaster, LSDRaster& SecondaryBurnRaster, string burned_column_name, string secondary_burned_column_name)
{
  float A_0 = 1;
  float this_movern;

  vector<float> movern_values;
  vector< vector<float> > chi_vecvec;
  vector<float> empty_vec;
  vector<float> this_chi_vec;
  int this_node;
  int n_nodes = int(node_sequence.size());

  // loop through m over n values
  for(int i = 0; i< n_movern; i++)
  {

    this_movern =  float(i)*delta_movern+start_movern;
    update_chi_data_map(FlowInfo, A_0, this_movern);

    cout << "m/n is: " << this_movern << endl;

    movern_values.push_back(this_movern);
    this_chi_vec = empty_vec;

    // now get the chi values for each node and push them into the chi_vecvec
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      this_chi_vec.push_back(chi_data_map[this_node]);
    }
    chi_vecvec.push_back(this_chi_vec);
  }
  cout << "Okay, I've got all the chi values in the vecvec." << endl;
  // okay, we are done getting all the chi values, now add these into the file

  ofstream chi_csv_out;
  cout << "Running the printing for movern. Filename is: " << filename << endl;
  chi_csv_out.open(filename.c_str());
  chi_csv_out << "source_key,basin_key,elevation";
  for (int i = 0; i< n_movern; i++)
  {
    chi_csv_out << ",m_over_n = " << movern_values[i];
  }
  chi_csv_out << "," << burned_column_name;
  chi_csv_out << "," << secondary_burned_column_name << endl;

  // now loop through all the nodes
  chi_csv_out.precision(5);
  int curr_row,curr_col;
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);

    chi_csv_out << source_keys_map[this_node] << ","
                 << baselevel_keys_map[this_node] << ","
                 << elev_data_map[this_node];

    for (int i = 0; i< n_movern; i++)
    {
      chi_csv_out << "," << chi_vecvec[i][n];
    }
    chi_csv_out << "," << BurnRaster.get_data_element(curr_row,curr_col);
    chi_csv_out << "," << SecondaryBurnRaster.get_data_element(curr_row,curr_col);
    chi_csv_out << endl;
  }
  chi_csv_out.close();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a series of simple profiles (chi-elevation) as a function of
// movern
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_profiles_as_fxn_movern_with_discharge(LSDFlowInfo& FlowInfo, string filename,
           float start_movern, float delta_movern, int n_movern, LSDRaster& Discharge)
{
  float A_0 = 1.0;
  float this_movern;

  vector<float> movern_values;
  vector< vector<float> > chi_vecvec;
  vector<float> empty_vec;
  vector<float> this_chi_vec;
  int this_node;
  int n_nodes = int(node_sequence.size());

  // loop through m over n values
  for(int i = 0; i< n_movern; i++)
  {

    this_movern =  float(i)*delta_movern+start_movern;

    // calculate chi
    float area_threshold = 0;

    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(this_movern,A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    cout << "m/n is: " << this_movern << endl;

    movern_values.push_back(this_movern);
    this_chi_vec = empty_vec;

    // now get the chi values for each node and push them into the chi_vecvec
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      this_chi_vec.push_back(chi_data_map[this_node]);
    }
    chi_vecvec.push_back(this_chi_vec);
  }
  cout << "Okay, I've got all the chi values in the vecvec." << endl;
  // okay, we are done getting all the chi values, now add these into the file

  ofstream chi_csv_out;
  cout << "Running the printing for movern. Filename is: " << filename << endl;
  chi_csv_out.open(filename.c_str());
  chi_csv_out << "source_key,basin_key,elevation";
  for (int i = 0; i< n_movern; i++)
  {
    chi_csv_out << ",m_over_n = " << movern_values[i];
  }
  chi_csv_out << endl;

  // now loop through all the nodes
  chi_csv_out.precision(5);
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    chi_csv_out << source_keys_map[this_node] << ","
                 << baselevel_keys_map[this_node] << ","
                 << elev_data_map[this_node];

    for (int i = 0; i< n_movern; i++)
    {
      chi_csv_out << "," << chi_vecvec[i][n];
    }
    chi_csv_out << endl;
  }
  chi_csv_out.close();
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a series of simple profiles (chi-elevation) as a function of
// movern
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_profiles_as_fxn_movern_with_discharge_and_burned_raster(LSDFlowInfo& FlowInfo, string filename,
           float start_movern, float delta_movern, int n_movern, LSDRaster& Discharge,
           LSDRaster& BurnRaster, string burned_column_name)
{
  float A_0 = 1.0;
  float this_movern;

  vector<float> movern_values;
  vector< vector<float> > chi_vecvec;
  vector<float> empty_vec;
  vector<float> this_chi_vec;
  int this_node;
  int n_nodes = int(node_sequence.size());
  int curr_row,curr_col;

  // loop through m over n values
  for(int i = 0; i< n_movern; i++)
  {

    this_movern =  float(i)*delta_movern+start_movern;

    // calculate chi
    float area_threshold = 0;

    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(this_movern,A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    cout << "m/n is: " << this_movern << endl;

    movern_values.push_back(this_movern);
    this_chi_vec = empty_vec;

    // now get the chi values for each node and push them into the chi_vecvec
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      this_chi_vec.push_back(chi_data_map[this_node]);
    }
    chi_vecvec.push_back(this_chi_vec);
  }
  cout << "Okay, I've got all the chi values in the vecvec." << endl;
  // okay, we are done getting all the chi values, now add these into the file

  ofstream chi_csv_out;
  cout << "Running the printing for movern. Filename is: " << filename << endl;
  chi_csv_out.open(filename.c_str());
  chi_csv_out << "source_key,basin_key,elevation";
  for (int i = 0; i< n_movern; i++)
  {
    chi_csv_out << ",m_over_n = " << movern_values[i];
  }
  chi_csv_out << "," << burned_column_name << endl;

  // now loop through all the nodes
  chi_csv_out.precision(5);
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    // get the row and column to use with the raster for burning
    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);

    chi_csv_out << source_keys_map[this_node] << ","
                 << baselevel_keys_map[this_node] << ","
                 << elev_data_map[this_node];

    for (int i = 0; i< n_movern; i++)
    {
      chi_csv_out << "," << chi_vecvec[i][n];
    }
    chi_csv_out << "," << BurnRaster.get_data_element(curr_row,curr_col);
    chi_csv_out << endl;
  }
  chi_csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a series of simple profiles (chi-elevation) as a function of
// movern
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void LSDChiTools::print_profiles_as_fxn_movern_with_discharge_and_secondary_raster(LSDFlowInfo& FlowInfo, string filename,
           float start_movern, float delta_movern, int n_movern, LSDRaster& Discharge,
           LSDRaster& BurnRaster,LSDRaster& SecondaryBurnRaster, string burned_column_name,string secondary_burned_column_name)
{
  float A_0 = 1.0;
  float this_movern;

  vector<float> movern_values;
  vector< vector<float> > chi_vecvec;
  vector<float> empty_vec;
  vector<float> this_chi_vec;
  int this_node;
  int n_nodes = int(node_sequence.size());
  int curr_row,curr_col;

  // loop through m over n values
  for(int i = 0; i< n_movern; i++)
  {

    this_movern =  float(i)*delta_movern+start_movern;

    // calculate chi
    float area_threshold = 0;

    LSDRaster this_chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(this_movern,A_0,
                                 area_threshold, Discharge);
    update_chi_data_map(FlowInfo, this_chi);

    cout << "m/n is: " << this_movern << endl;

    movern_values.push_back(this_movern);
    this_chi_vec = empty_vec;

    // now get the chi values for each node and push them into the chi_vecvec
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      this_chi_vec.push_back(chi_data_map[this_node]);
    }
    chi_vecvec.push_back(this_chi_vec);
  }
  cout << "Okay, I've got all the chi values in the vecvec." << endl;
  // okay, we are done getting all the chi values, now add these into the file

  ofstream chi_csv_out;
  cout << "Running the printing for movern. Filename is: " << filename << endl;
  chi_csv_out.open(filename.c_str());
  chi_csv_out << "source_key,basin_key,elevation";
  for (int i = 0; i< n_movern; i++)
  {
    chi_csv_out << ",m_over_n = " << movern_values[i];
  }
  chi_csv_out << "," << burned_column_name;
  chi_csv_out << "," << secondary_burned_column_name << endl;
  // now loop through all the nodes
  chi_csv_out.precision(5);
  for (int n = 0; n< n_nodes; n++)
  {
    this_node = node_sequence[n];

    // get the row and column to use with the raster for burning
    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);

    chi_csv_out << source_keys_map[this_node] << ","
                 << baselevel_keys_map[this_node] << ","
                 << elev_data_map[this_node];

    for (int i = 0; i< n_movern; i++)
    {
      chi_csv_out << "," << chi_vecvec[i][n];
    }
    chi_csv_out << "," << BurnRaster.get_data_element(curr_row,curr_col);
    chi_csv_out << "," << SecondaryBurnRaster.get_data_element(curr_row,curr_col);
    chi_csv_out << endl;
  }
  chi_csv_out.close();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the source node of a given key
// Returns a map were the key is the basin key and the value is the
// node index of the outlet node
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int,int> LSDChiTools::get_outlet_node_from_basin_key_map()
{
  map<int,int> outlet_node_from_basin_key;
  map<int,int>::iterator iter = key_to_baselevel_map.begin();

  while(iter != key_to_baselevel_map.end())
  {
    int outlet_node  = iter->first;
    int basin_key = iter->second;

    outlet_node_from_basin_key[basin_key] = outlet_node;
    iter++;
  }
  return outlet_node_from_basin_key;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the source node of a given key
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChiTools::get_source_from_source_key(int source_key)
{
  // first get the source node of the reference channel
  int source_node = -9999;
  map<int,int>::iterator iter = key_to_source_map.begin();
  while(iter != key_to_source_map.end())
  {
    if (iter->second == source_key)
    {
      source_node = iter->first;
      //cout << "I found the source key AWESOME! The source key is: " << source_node << endl;
    }
    iter++;
  }

  if ( source_node == -9999 )
  {
    cout << "LSDChiTools::get_starting_node_of_source " << endl;
    cout << "FATAL ERROR: This source is not in the channel network. Source is: " << source_key << endl;
    exit(EXIT_FAILURE);
  }

  return source_node;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Given a source key, find the index of the starting node in the node sequence
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChiTools::get_starting_node_of_source(int source_key)
{
  int source_node = get_source_from_source_key(source_key);
  int this_starting_node = -9999;

  // Find the node sequence index of this node;
  int this_ns_node = -1;

  bool not_starting_node = true;
  while (not_starting_node)
  {
    this_ns_node++;
    if(node_sequence[this_ns_node] == source_node)
    {
      not_starting_node = false;
      this_starting_node = this_ns_node;
    }
  }
  //cout << "The starting node in the sequence is: " << this_starting_node << endl;

  return this_starting_node;
}


int LSDChiTools::get_ending_node_of_source(LSDFlowInfo& FlowInfo,int source_key)
{
  int starting_node_index = get_starting_node_of_source(source_key);
  int node_indenter = starting_node_index;
  int temp1, temp2, temp_node_index = starting_node_index, temp_node = node_sequence[starting_node_index];
  int n = node_sequence.size();

  

  while(source_keys_map[temp_node] == source_key && temp_node_index < n)
  {

    node_indenter ++;
    temp_node_index = node_indenter;
    FlowInfo.retrieve_receiver_information(node_sequence[temp_node_index], temp_node, temp1, temp2);

    

  }
  int ending_node = temp_node;


  return ending_node;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets the number of channels
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDChiTools::get_number_of_channels()
{
  int n_channels = int(key_to_source_map.size());
  return n_channels;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Extract the elevation and chi data from a channel
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_chi_elevation_data_of_channel(LSDFlowInfo& FlowInfo, int source_key,
                                vector<float>& chi_data, vector<float>& elevation_data)
{
  int n_channels = int(key_to_source_map.size());
  if (source_key>=n_channels)
  {
    cout << "LSDChiTools::get_chi_elevation_data_of_channel FATAL ERROR" << endl;
    cout << "This source key is not in channel network" << endl;
    exit(EXIT_FAILURE);
  }

  // get source node
  int starting_source = get_source_from_source_key(source_key);
  //cout << "LSDChiTools::get_chi_elevation_data_of_channel, Starting source is: " << starting_source << endl;

  vector<float> this_chi;
  vector<float> this_elevation;

  // add the source to the chi elevation vectors
  this_chi.push_back(chi_data_map[starting_source]);
  this_elevation.push_back(elev_data_map[starting_source]);

  //cout << "Starting chi is: " << chi_data_map[starting_source] << endl;

  // now work downstream until you get to a different source or
  // a baselevel node
  bool is_end = false;
  int current_node = starting_source;
  int receiver_node,receiver_row,receiver_col;
  int this_source_key;
  while(is_end == false)
  {
    FlowInfo.retrieve_receiver_information(current_node,receiver_node, receiver_row,receiver_col);

    if(current_node == receiver_node)
    {
      // this is a baselelvel node, simply switch on the is_end boolean
      is_end = true;
    }
    else
    {
      this_source_key = source_keys_map[receiver_node];
      if (this_source_key != source_key)
      {
        //cout << "I made it to the end of this channel" << endl;
      }
      else
      {
        this_chi.push_back(chi_data_map[receiver_node]);
        this_elevation.push_back(elev_data_map[receiver_node]);
      }
    }
    // increment the node downstream
    current_node = receiver_node;

  }

  elevation_data = this_elevation;
  chi_data = this_chi;
  //cout << "Starting chi is: " << chi_data[0] << " and ending chi is: " << chi_data[n_nodes-1] << endl;

  // For debugging
  //int n_nodes = int(elevation_data.size());
  //for(int i= 0; i<n_nodes; i++)
  //{
  //  cout << chi_data[i] << "," << elevation_data[i] << endl;
  //}
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Project data onto a reference chi-elevation profile
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDChiTools::project_data_onto_reference_channel(vector<float>& reference_chi,
                                 vector<float>& reference_elevation, vector<float>& trib_chi,
                                 vector<float>& trib_elevation)
{
  // How this works is that you take the tributary elevations and then
  // determine the elevation on the reference at the same chi. This is done by
  // interpolating the elevation as a linear fit between the two adjacent chi
  // points on the reference channel.
  vector<float> joint_chi;
  vector<float> trib_joint_elev;
  vector<float> ref_joint_elev;
  vector<float> residuals;


  float this_chi;
  float max_ref_chi;
  float min_ref_chi;
  float max_trib_chi;
  float min_trib_chi;

  int n_ref_nodes = int(reference_chi.size());
  int n_trib_nodes = int(trib_chi.size());
  if (n_ref_nodes <= 1 || n_trib_nodes <= 1)
  {
    cout << "LSDChiTools::project_data_onto_reference_channel WARNING" << endl;
    cout << "The reference channel has 1 or zero nodes." << endl;
    return residuals;
  }
  else
  {
    max_ref_chi = reference_chi[0];
    min_ref_chi = reference_chi[n_ref_nodes-1];

    max_trib_chi = trib_chi[0];
    min_trib_chi = trib_chi[n_trib_nodes-1];
  }

  // test to see if there is overlap
  if(min_trib_chi > max_ref_chi || max_trib_chi < min_ref_chi)
  {
    cout << "LSDChiTools::project_data_onto_reference_channel These channels do not overlap." << endl;
    return residuals;
  }

  // The reference chis monotonically decrease so we will keep track of what
  // indices the bounding chi points are.
  int start_ref_index = 0;
  int end_ref_index = 1;

  // begin by ramping up to the first node within the reference vector
  float this_trib_chi = trib_chi[0];
  int this_node=0;
  if(this_trib_chi > max_ref_chi)
  {
    // the node in the trib is outside of the reference frame. Increment
    // the index until it is.
    while(this_node < n_trib_nodes && trib_chi[this_node] <= max_ref_chi)
    {
      this_node++;
    }
  }

  this_chi = trib_chi[this_node];
  float ref_chi_upstream = reference_chi[start_ref_index];
  float ref_chi_downstream =  reference_chi[end_ref_index];

  //cout << "The number of trib nodes is: " << n_trib_nodes << endl;
  // now ramp up the the start ref index and end ramp index
  bool found_joint_chi;
  for (int i = this_node; i<n_trib_nodes; i++)
  {
    found_joint_chi = false;
    //cout << "trib node: " << i << endl;
    // Get the chi for this node
    this_chi = trib_chi[i];

    // now test if it is between the upstream and downstream chi coordinates in the reference vector
    if (this_chi < ref_chi_upstream && this_chi > ref_chi_downstream)
    {
      // It is between these reference chi values!
      //cout << "FOUND CHI This chi is: " << this_chi << " and bounds are: " << ref_chi_upstream << "," << ref_chi_downstream << endl;
      found_joint_chi = true;
    }
    else
    {
      // we didn't find the chi, we need to move through the reference vector to find
      // the chi value
      bool found_ref_nodes = false;
      while (end_ref_index < n_ref_nodes && found_ref_nodes == false)
      {
        start_ref_index++;
        end_ref_index++;
        ref_chi_upstream = reference_chi[start_ref_index];
        ref_chi_downstream =  reference_chi[end_ref_index];
        if (this_chi < ref_chi_upstream && this_chi > ref_chi_downstream)
        {
          found_ref_nodes = true;
          found_joint_chi = true;
          //cout << "FOUND CHI This chi is: " << this_chi << " and bounds are: " << ref_chi_upstream << "," << ref_chi_downstream << endl;
        }
      }
      // There is different logic if we reached the end of the reference vector
      if (end_ref_index == n_ref_nodes-1)
      {
        //cout << "I am at the end of the reference vector" << endl;
        i = n_trib_nodes-1;
      }
    }

    if(found_joint_chi)
    {
      // we need to calculate the eleavtion on the reference vector
      float dist_ref = ref_chi_upstream-ref_chi_downstream;
      float chi_frac = (this_chi-ref_chi_downstream)/dist_ref;
      float joint_elev = chi_frac*(reference_elevation[start_ref_index]-reference_elevation[end_ref_index])
                             +reference_elevation[end_ref_index];
      joint_chi.push_back(this_chi);
      trib_joint_elev.push_back(trib_elevation[i]);
      ref_joint_elev.push_back(joint_elev);
      residuals.push_back(trib_elevation[i]-joint_elev);
    }

  }

  //for(int i = 0; i<int(joint_chi.size()); i++)
  //{
  //  cout << "residual[" << i << "]: "<< residuals[i] << endl;
  //}

  // this section is for debugging
  bool print_for_debugging = false;
  if (print_for_debugging)
  {
    ofstream chans_out;
    chans_out.open("Test_project.csv");
    chans_out << "chi,elev,channel_code" << endl;
    for(int i = 0; i<n_ref_nodes; i++)
    {
      chans_out << reference_chi[i] << "," << reference_elevation[i] <<",0" <<endl;
    }
    for(int i = 0; i<n_trib_nodes; i++)
    {
      chans_out << trib_chi[i] << "," << trib_elevation[i] <<",1" <<endl;
    }
    for(int i = 0; i< int(joint_chi.size()); i++)
    {
      chans_out << joint_chi[i] << "," << trib_joint_elev[i] <<",2" <<endl;
    }
    for(int i = 0; i<int(joint_chi.size()); i++)
    {
      chans_out << joint_chi[i] << "," << ref_joint_elev[i] <<",3" <<endl;
    }
    chans_out.close();
  }


  return residuals;

}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Project data onto a reference chi-elevation profile
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDChiTools::project_points_onto_reference_channel(vector<float>& reference_chi,
                                 vector<float>& reference_elevation, vector<float>& trib_chi,
                                 vector<float>& trib_elevation,
                                 vector<float> chi_distances_to_test)
{
  //cout << endl << endl << "=========="  << endl << "projecting points" << endl;

  // How this works is that you take the tributary elevations and then
  // determine the elevation on the reference at the same chi. This is done by
  // interpolating the elevation as a linear fit between the two adjacent chi
  // points on the reference channel.
  vector<float> joint_chi;
  vector<float> trib_joint_elev;
  vector<float> ref_joint_elev;
  vector<float> residuals;


  //float this_chi;
  float max_ref_chi;
  float min_ref_chi;
  float max_trib_chi;
  float min_trib_chi;

  int n_ref_nodes = int(reference_chi.size());
  int n_trib_nodes = int(trib_chi.size());
  if (n_ref_nodes <= 1 || n_trib_nodes <= 1)
  {
    //cout << "LSDChiTools::project_points_onto_reference_channel WARNING" << endl;
    //cout << "The reference channel has 1 or zero nodes." << endl;
    return residuals;
  }
  else
  {
    max_ref_chi = reference_chi[0];
    min_ref_chi = reference_chi[n_ref_nodes-1];

    max_trib_chi = trib_chi[0];
    min_trib_chi = trib_chi[n_trib_nodes-1];
  }

  //cout << "Max ref:  " << max_ref_chi << " min ref:  " << min_ref_chi << endl;
  //cout << "Max test: " << max_trib_chi << " min test:  " << min_trib_chi << endl;


  // test to see if there is overlap
  if(min_trib_chi > max_ref_chi || max_trib_chi < min_ref_chi)
  {
    //cout << "LSDChiTools::project_points_onto_reference_channel These channels do not overlap." << endl;
    return residuals;
  }

  // get the length of the trib
  float length_of_trib = max_trib_chi-min_trib_chi;

  // get the distances upstream of the base of the tributary.
  vector<float> chi_upstream;
  for(int i = 0; i<n_trib_nodes; i++)
  {
    chi_upstream.push_back(trib_chi[i]-min_trib_chi);
  }

  int n_points = int(chi_distances_to_test.size());
  // now we search down the chi_distances vector to get the distances along
  // the main stem corresponding to the correct place on the reference channel
  // The reference chis monotonically decrease so we will keep track of what
  // indices the bounding chi points are.

  for(int i = 0; i<n_points; i++ )
  {
    if (chi_distances_to_test[i] > length_of_trib)
    {
      residuals.push_back(-9999);
    }
    else
    {
      // need to get to the node in the trib that is less than the chi to test
      bool got_there_yet = false;
      int this_trib_node = 0;
      while ( (got_there_yet == false)  && this_trib_node!= n_trib_nodes)
      {
        //cout << "chi distance: " << chi_upstream[this_trib_node] << endl;
        if(chi_distances_to_test[i] > chi_upstream[this_trib_node])
        {
          //cout << "The distance is greater than the distance of this node!" << endl;
          got_there_yet = true;
        }
        else
        {
          this_trib_node++;
        }
      }
      //cout << "The distance needed is: " << chi_distances_to_test[i] << endl;
      //cout << "Found the distance: " <<  chi_upstream[this_trib_node] << endl;
      //cout << "The chi location for testing on the tributary is: " << trib_chi[this_trib_node] << endl;

      // we have either got the node or reached the end of the tributary
      if (this_trib_node == n_trib_nodes)
      {
        // we got to the end without finding the correct node. Return nodata
        residuals.push_back(-9999);
      }
      else
      {
        // now you need to find the bounding nodes on the mainstem
        bool found_joint_chi;
        int start_ref_index = 0;
        int end_ref_index = 1;
        float ref_chi_upstream = reference_chi[start_ref_index];
        float ref_chi_downstream =  reference_chi[end_ref_index];
        float this_chi = trib_chi[this_trib_node];

        // now test if it is between the upstream and downstream chi coordinates in the reference vector
        if (this_chi  < ref_chi_upstream && this_chi > ref_chi_downstream)
        {
          // It is between these reference chi values!
          //cout << "FOUND CHI This chi is: " << this_chi << " and bounds are: " << ref_chi_upstream << "," << ref_chi_downstream << endl;
          found_joint_chi = true;
        }
        else
        {
          // we didn't find the chi, we need to move through the reference vector to find
          // the chi value
          bool found_ref_nodes = false;
          while (end_ref_index < n_ref_nodes && found_ref_nodes==false)
          {
            start_ref_index++;
            end_ref_index++;
            ref_chi_upstream = reference_chi[start_ref_index];
            ref_chi_downstream =  reference_chi[end_ref_index];
           if (this_chi  < ref_chi_upstream && this_chi > ref_chi_downstream)
           {

             found_ref_nodes = true;
             found_joint_chi = true;
             //cout << "FOUND CHI This chi is: " <<this_chi  << " and bounds are: " << ref_chi_upstream << "," << ref_chi_downstream << endl;
           }
          }
          // There is different logic if we reached the end of the reference vector
          if (end_ref_index == n_ref_nodes-1)
          {
            //cout << "I am at the end of the reference vector" << endl;
            i = n_trib_nodes-1;
          }
        }

        if(found_joint_chi)
        {
          // we need to calculate the eleavtion on the reference vector
          float dist_ref = ref_chi_upstream-ref_chi_downstream;
          float chi_frac = (this_chi-ref_chi_downstream)/dist_ref;
          float joint_elev = chi_frac*(reference_elevation[start_ref_index]-reference_elevation[end_ref_index])
                             +reference_elevation[end_ref_index];
          joint_chi.push_back(this_chi);
          trib_joint_elev.push_back(trib_elevation[this_trib_node]);
          ref_joint_elev.push_back(joint_elev);
          residuals.push_back(trib_elevation[this_trib_node]-joint_elev);

          //cout << "trib chi: " << this_chi << " elev trib: " << trib_elevation[this_trib_node] << " test elev: " << joint_elev << endl;

        }
        else
        {
          residuals.push_back(-9999);
        }
      }
    }
  }

  // get rid of nodata residuals
  int n_residuals = int(residuals.size());
  vector<float> thinned_residuals;
  for(int i = 0; i<n_residuals; i++)
  {
    if (residuals[i] != -9999)
    {
      thinned_residuals.push_back(residuals[i]);
    }
  }
  residuals = thinned_residuals;



  // this section is for debugging
  bool print_for_debugging = false;
  if (print_for_debugging)
  {
    ofstream chans_out;
    chans_out.open("Test_project.csv");
    chans_out << "chi,elev,channel_code" << endl;
    for(int i = 0; i<n_ref_nodes; i++)
    {
      chans_out << reference_chi[i] << "," << reference_elevation[i] <<",0" <<endl;
    }
    for(int i = 0; i<n_trib_nodes; i++)
    {
      chans_out << trib_chi[i] << "," << trib_elevation[i] <<",1" <<endl;
    }
    for(int i = 0; i< int(joint_chi.size()); i++)
    {
      chans_out << joint_chi[i] << "," << trib_joint_elev[i] <<",2" <<endl;
    }
    for(int i = 0; i<int(joint_chi.size()); i++)
    {
      chans_out << joint_chi[i] << "," << ref_joint_elev[i] <<",3" <<endl;
    }
    chans_out.close();
  }


  return residuals;

}







//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This does slope-area analysis
// The strategy is to loop through each source
// progress down each but do not cross tributary junctions
// The data is not at all the points since we are averaging across multiple
// pixels.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_slope_area_data(LSDFlowInfo& FlowInfo, float vertical_interval,
                                  vector<int>& midpoint_nodes, vector<float>& slopes)
{
  // set verbose to true if you want to print the data as you go along
  bool verbose = false;

  // these are the vectors that will hold all the data
  // Note we don't need to keep track of all the area, elevation, location, etc.
  // because these are all keyed to to the node index of the midpoint.
  // The only thing we need to retain is the slope.
  vector<float> SA_slope;
  vector<int> SA_midpoint_node;

  // used to calculate the S-A data
  float half_interval = vertical_interval*0.5;
  //float midpoint_area;
  float upstream_elevation;
  float upstream_flow_distance;
  float downstream_elevation = 0;  // need to set these to avoid compiler warning
  float downstream_flow_distance = 0;
  float midpoint_node;
  float target_end_interval_elevation;
  float target_midpoint_interval_elevation;

  // we loop through every source
  int n_sources = int(ordered_source_nodes.size());
  if(verbose)
  {
    cout << "Super times, you are going to do a S-A analysis." << endl;
    cout << "Number of sources is: " << n_sources << endl;
  }


  int top_interval_node;
  int search_node;
  int this_source_node;
  int this_source_key;
  int row,col;
  for(int s = 0; s<n_sources; s++)
  {
    top_interval_node = ordered_source_nodes[s];
    this_source_node = top_interval_node;
    this_source_key = source_keys_map[this_source_node];
    // now trace downstream until you first get to the midpoint,
    // and then to the final node.
    bool is_this_final_node = false;
    bool is_end_interval;
    bool is_midpoint_interval;
    int last_node;

    midpoint_node = 0;    // this will be updated later

    if (verbose)
    {
      cout << "=====================" << endl;
      cout << "I am starting on a channel with source node " << this_source_node << " and key: " << this_source_key << endl;
    }

    // now, start at the top node of each channel and
    // work down. The algorithm looks downstream until it hits
    // the midpoint, and then continues until it hits
    // the end point
    // if it encounters the end of the channel before it hits
    // the end point the loop exits
    while (is_this_final_node == false)
    {
      // get the upstream elevation and flow distance
      upstream_elevation = elev_data_map[top_interval_node];
      upstream_flow_distance = flow_distance_data_map[top_interval_node];

      target_end_interval_elevation = upstream_elevation-vertical_interval;
      target_midpoint_interval_elevation = upstream_elevation-half_interval;

      if(verbose)
      {
        cout << "Looking for a midpoint for top interval node: " << top_interval_node << endl;
        cout << "Top z: " << upstream_elevation << " and top fd: " << upstream_flow_distance << endl;
      }

      // this gets the receiver (placed into the seach node)
      FlowInfo.retrieve_receiver_information(top_interval_node,
                     search_node, row, col);

      // check to see if this is the last element or in a tributary
      if (search_node == top_interval_node || this_source_key != source_keys_map[search_node])
      {
        is_this_final_node = true;
      }
      else
      {
        // reset midpoint and end flags
        is_end_interval = false;
        is_midpoint_interval = false;

        // now work downstream
        while ( is_this_final_node == false &&  is_end_interval == false)
        {
          //cout << "search_node: " << search_node << " elev: " << elevations[chan][search_node]
          //     << " and target mp, end: " << target_mp_interval_elevations << " " << target_end_interval_elevations << endl;

          // see if search node is the midpoint node
          if ( elev_data_map[search_node] <= target_midpoint_interval_elevation &&  is_midpoint_interval == false)
          {
            //midpoint_area = drainage_area_data_map[search_node];
            midpoint_node = search_node;

            // set midpoint flag so it doens't collect downstream nodes
            is_midpoint_interval = true;

            if(verbose)
            {
              cout << endl << endl << "=========" << endl << "I found the midpoint!" << endl;
            }
          }

          // see if the search node is the end node
          if (elev_data_map[search_node] <= target_end_interval_elevation)
          {
            downstream_elevation = elev_data_map[search_node];
            downstream_flow_distance = flow_distance_data_map[search_node];

            // make sure the code knows this is the end, the only end, my friend.
            is_end_interval = true;;
          }

          // now move downstream
          last_node = search_node;
          FlowInfo.retrieve_receiver_information(last_node,
                     search_node, row, col);

          // test is this is the end
          if (search_node == last_node || this_source_key != source_keys_map[search_node])
          {
            is_this_final_node = true;
          }

        }

        // if is_end_interval is true, that means it found the end interval.
        // record the information. If it didn't reach the end flag that means
        // the previous node hit the final node in the channel before finding the
        // end interval.
        if (is_end_interval)
        {
          float slope = (upstream_elevation-downstream_elevation)/
                       (upstream_flow_distance-downstream_flow_distance);


          // It the slope is zero or less than zero we ignore the data
          if(slope > 0)
          {
            // record the data
            SA_midpoint_node.push_back(midpoint_node);
            SA_slope.push_back(slope);
          }
        }             // if statement for recording data on S-A to data containers
      }               // end check if tributary only has one node
      // now we need to update the top interval node
      int old_top_interval = top_interval_node;
      FlowInfo.retrieve_receiver_information(old_top_interval,
                     top_interval_node, row, col);
      if(verbose)
      {
        cout << "Resetting the top interval node, old:  " << old_top_interval << " and new: " << top_interval_node << endl;
      }

    }                 // check if this is the final node of the source trib
  }                   // end sources loop (at this point we go to the next source)

  midpoint_nodes = SA_midpoint_node;
  slopes = SA_slope;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function bootstraps the S-A data to get confidence intervals on the
// best fit m/n ratio
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::bootstrap_slope_area_data(LSDFlowInfo& FlowInfo,
                                          vector<int>& SA_midpoint_node,
                                          vector<float>& SA_slope,
                                          int N_iterations,
                                          float bootstrap_keep_data_prob,
                                          string filename)
{

  // open the file
  ofstream bootstrap_out;
  bootstrap_out.open(filename.c_str());
  bootstrap_out << "basin_key,MS_minimum,MS_first_quartile,MS_median,MS_third_quartile,MS_maximum,MS_MAD,All_minimum,All_first_quartile,All_median,All_third_quartile,All_maximum,All_MAD" << endl;


  // we will store the data in maps where the key is the source node
  // This is because we will bin the data by source.
  vector<float> empty_vec;
  map< int, vector<float> > log_slope_map;
  map< int, vector<float> > log_area_map;
  map< int, int > basin_key_of_this_source_map;

  map< int, vector<float> > log_slope_map_by_basin;
  map< int, vector<float> > log_area_map_by_basin;

  int this_node;
  int this_source_key;
  int this_basin_key;

  float l10_DA, l10_S;

  int n_nodes = int(SA_midpoint_node.size());
  if (n_nodes <= 0)
  {
    cout << "Trying to print SA data but there doesn't seem to be any." << endl;
    cout << "Have you run the automator and gathered the sources yet?" << endl;
  }
  else
  {
    // get the data vectors out
    for (int n = 0; n< n_nodes; n++)
    {
      // get the source node
      this_node = SA_midpoint_node[n];
      this_source_key = source_keys_map[this_node];
      this_basin_key = baselevel_keys_map[this_node];
      //cout << "This source key is: " << this_source_key << endl;

      // see if we have a vector for that source node
      if( log_area_map.find(this_source_key) == log_area_map.end())
      {
        log_area_map[this_source_key] = empty_vec;
        log_slope_map[this_source_key] = empty_vec;
      }
      // see if we have a vector for that source node
      if( log_area_map_by_basin.find(this_basin_key) == log_area_map_by_basin.end())
      {
        log_area_map_by_basin[this_basin_key] = empty_vec;
        log_slope_map_by_basin[this_basin_key] = empty_vec;
      }
      // check if we have the basin of this source
      if (basin_key_of_this_source_map.find(this_source_key) == basin_key_of_this_source_map.end() )
      {
        basin_key_of_this_source_map[this_source_key] = baselevel_keys_map[this_node];
      }

      // add to this source's log S, log A data. We will later use these to bin
      if(drainage_area_data_map[this_node] > 0 && SA_slope[n] > 0)
      {
        l10_DA =  log10(drainage_area_data_map[this_node]);
        l10_S = log10(SA_slope[n]);

        log_area_map[this_source_key].push_back( l10_DA );
        log_slope_map[this_source_key].push_back( l10_S );

        log_area_map_by_basin[this_basin_key].push_back( l10_DA );
        log_slope_map_by_basin[this_basin_key].push_back( l10_S );
      }
      //log_area_map[this_source_key].push_back( drainage_area_data_map[this_node] );
      //log_slope_map[this_source_key].push_back( SA_slope[n] );
    }
  }
  // Okay, so at this point we have data maps that have the slope and areas in maps
  // with the key as the source, and a "basin_key_of_this_source" where
  // the key is this source number and the value is the basin.
  // We can loop through the basins or all sources to get at the
  // basin-based data


  map<int, vector<float> > log_area_mainstem;
  map<int, vector<float> > log_slope_mainstem;

  //map<int, vector<float> > mainstem_summary_stats_by_basin;
  //map<int, vector<float> > all_summary_stats_by_basin;

  // This loops through the basin keys, getting out the mainstem sources
  map<int,int>::iterator iter = source_node_of_mainstem_map.begin();
  while (iter != source_node_of_mainstem_map.end())
  {
    // get the basin
    int basin_key = iter->first;
    // get the mainstem source
    int mainstem_source_node = iter->second;
    int this_source_key =  source_keys_map[mainstem_source_node];

    cout << "Bootstrapping, basin_key is: " << basin_key << endl;

    // bootstrap the main stem
    vector<float> this_log_area_mainstem = log_area_map[this_source_key];
    vector<float> this_log_slope_mainstem = log_slope_map[this_source_key];
    vector<float> MS_summary = bootstrap_linear_regression(this_log_area_mainstem, this_log_slope_mainstem, N_iterations,bootstrap_keep_data_prob);


    vector<float> this_log_area_all = log_area_map_by_basin[basin_key];
    vector<float> this_log_slope_all = log_slope_map_by_basin[basin_key];
    vector<float> All_summary = bootstrap_linear_regression(this_log_area_all, this_log_slope_all, N_iterations,bootstrap_keep_data_prob);


    // print results to file
    bootstrap_out << basin_key << "," << MS_summary[0] << "," << MS_summary[1] << ","
                  << MS_summary[2] << "," << MS_summary[3] << "," << MS_summary[4]
                  << "," << MS_summary[8] << ","
                  << All_summary[0] << "," << All_summary[1] <<","
                  << All_summary[2] << "," << All_summary[3] << "," << All_summary[4]
                  << "," << All_summary[8] << endl;

    iter++;
  }

  // Now go through all the data, collating data for each basin




}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This bins the data
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::bin_slope_area_data(LSDFlowInfo& FlowInfo,
                                          vector<int>& SA_midpoint_node,
                                          vector<float>& SA_slope,
                                          float log_bin_width,
                                          string filename)
{

  vector<float> empty_vec;

  // we will store the data in maps where the key is the source node
  // This is because we will bin the data by source.
  //map< int, vector<float> > slope;
  //map< int, vector<float> > area;
  map< int, vector<float> > log_slope_map;
  map< int, vector<float> > log_area_map;
  map< int, int > basin_key_of_this_source_map;

  int this_node;
  int this_source_key;

  int n_nodes = int(SA_midpoint_node.size());
  if (n_nodes <= 0)
  {
    cout << "Trying to print SA data but there doesn't seem to be any." << endl;
    cout << "Have you run the automator and gathered the sources yet?" << endl;
  }
  else
  {
    // get the data vectors out
    for (int n = 0; n< n_nodes; n++)
    {
      // get the source node
      this_node = SA_midpoint_node[n];
      this_source_key = source_keys_map[this_node];
      //cout << "This source key is: " << this_source_key << endl;

      // see if we have a vector for that source node
      if( log_area_map.find(this_source_key) == log_area_map.end())
      {
        log_area_map[this_source_key] = empty_vec;
        log_slope_map[this_source_key] = empty_vec;
      }
      // check if we have the basin of this source
      if (basin_key_of_this_source_map.find(this_source_key) == basin_key_of_this_source_map.end() )
      {
        basin_key_of_this_source_map[this_source_key] = baselevel_keys_map[this_node];
      }

      // add to this source's log S, log A data. We will later use these to bin
      if(drainage_area_data_map[this_node] > 0 && SA_slope[n] > 0)
      {
        log_area_map[this_source_key].push_back( log10(drainage_area_data_map[this_node]) );
        log_slope_map[this_source_key].push_back( log10(SA_slope[n]) );
      }
    }
  }

  // now we bin the data
  vector<float> midpoints_output;
  vector<float> MeanX_output;
  vector<float> MedianX_output;
  vector<float> StandardDeviationX_output;
  vector<float> StandardErrorX_output;
  vector<float> MADX_output;

  vector<float> MeanY_output;
  vector<float> MinimumY_output;
  vector<float> FirstQuartileY_output;
  vector<float> MedianY_output;
  vector<float> ThirdQuartileY_output;
  vector<float> MaximumY_output;
  vector<float> StandardDeviationY_output;
  vector<float> StandardErrorY_output;
  vector<float> MADY_output;
  vector<int> number_observations_output;
  float NoDataValue = -9999;


  // these are the vectors holding all the compiled information
  vector<int> binned_basin_keys;
  vector<int> binned_source_keys;
  vector<float> binned_logA_means;
  vector<float> binned_logA_midpoints;
  vector<float> binned_logA_medians;
  vector<float> binned_logA_stdErr;
  vector<float> binned_logA_MAD;

  vector<float> binned_logS_means;
  vector<float> binned_logS_medians;
  vector<float> binned_logS_FirstQuartile;
  vector<float> binned_logS_ThirdQuartile;
  vector<float> binned_logS_StandardDeviation;
  vector<int> binnned_NObvs;

  // loop through all the source nodes
  map<int, vector<float> >::iterator it;
  for(it = log_area_map.begin(); it != log_area_map.end(); ++it)
  {
    this_source_key =  it->first;

    //cout << "The source key is: " << this_source_key << endl;

    // extract the log S-log A data for this source
    vector<float> log_area = log_area_map[this_source_key];
    vector<float> log_slope = log_slope_map[this_source_key];

    // this gets the binned data for this particular tributary
    bin_data(log_area, log_slope, log_bin_width, midpoints_output, MeanX_output,
             MedianX_output, StandardDeviationX_output, StandardErrorX_output,
             MADX_output, MeanY_output, MinimumY_output,FirstQuartileY_output,
             MedianY_output, ThirdQuartileY_output, MaximumY_output,
             StandardDeviationY_output, StandardErrorY_output, MADY_output,
             number_observations_output, NoDataValue);

    // now we need to add this information to the master vectors
    int n_bins = int(midpoints_output.size());
    int n_Obvs;
    for(int i = 0; i<n_bins; i++)
    {
      n_Obvs = number_observations_output[i];

      // only record data if there are enough observations
      if(n_Obvs > 0)
      {
        binned_basin_keys.push_back(basin_key_of_this_source_map[this_source_key]);
        binned_source_keys.push_back(this_source_key);
        binned_logA_means.push_back(MeanX_output[i]);
        binned_logA_midpoints.push_back(midpoints_output[i]);
        binned_logA_medians.push_back(MedianX_output[i]);
        binned_logA_stdErr.push_back(StandardErrorX_output[i]);
        binned_logA_MAD.push_back(MADX_output[i]);

        binned_logS_means.push_back(MeanY_output[i]);
        binned_logS_medians.push_back(MedianY_output[i]);
        binned_logS_FirstQuartile.push_back(FirstQuartileY_output[i]);
        binned_logS_ThirdQuartile.push_back(ThirdQuartileY_output[i]);
        binned_logS_StandardDeviation.push_back(StandardDeviationY_output[i]);
        binnned_NObvs.push_back(n_Obvs);
      }
    }
  }


  // now print to file
  int n_data_points = int(binnned_NObvs.size());
  ofstream  binned_out;
  binned_out.open(filename.c_str());
  binned_out << "basin_key,source_key,midpoints_log_A,mean_log_A,median_log_A,logA_stdErr,logA_MAD,mean_log_S,median_log_S,logS_FirstQuartile,logS_ThirdQuartile,logS_stdDev,n_observations" << endl;
  for(int i = 0; i<n_data_points; i++)
  {
    binned_out << binned_basin_keys[i] << ","
               << binned_source_keys[i] << ","
               << binned_logA_midpoints[i] << ","
               << binned_logA_means[i] << ","
               << binned_logA_medians[i] << ","
               << binned_logA_stdErr[i] << ","
               << binned_logA_MAD[i] << ","
               << binned_logS_means[i] << ","
               << binned_logS_medians[i] << ","
               << binned_logS_FirstQuartile[i] << ","
               << binned_logS_ThirdQuartile[i] << ","
               << binned_logS_StandardDeviation[i] << ","
               << binnned_NObvs[i] << endl;
  }
  binned_out.close();
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This segments the binned slope data
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::segment_binned_slope_area_data(LSDFlowInfo& FlowInfo,
                                          vector<int>& SA_midpoint_node,
                                          vector<float>& SA_slope,
                                          float log_bin_width,
                                          int minimum_segment_length,
                                          string filename)
{


  vector<float> empty_vec;

  // we will store the data in maps where the key is the source node
  // This is because we will bin the data by source.
  //map< int, vector<float> > slope;
  //map< int, vector<float> > area;
  map< int, vector<float> > log_slope_map;
  map< int, vector<float> > log_area_map;
  map< int, int > basin_key_of_this_source_map;

  int this_node;
  int this_source_key;

  int n_nodes = int(SA_midpoint_node.size());
  if (n_nodes <= 0)
  {
    cout << "Trying to print SA data but there doesn't seem to be any." << endl;
    cout << "Have you run the automator and gathered the sources yet?" << endl;
  }
  else
  {
    // get the data vectors out
    for (int n = 0; n< n_nodes; n++)
    {
      // get the source node
      this_node = SA_midpoint_node[n];
      this_source_key = source_keys_map[this_node];
      //cout << "This source key is: " << this_source_key << endl;

      // see if we have a vector for that source node
      if( log_area_map.find(this_source_key) == log_area_map.end())
      {
        log_area_map[this_source_key] = empty_vec;
        log_slope_map[this_source_key] = empty_vec;
      }
      // check if we have the basin of this source
      if (basin_key_of_this_source_map.find(this_source_key) == basin_key_of_this_source_map.end() )
      {
        basin_key_of_this_source_map[this_source_key] = baselevel_keys_map[this_node];
      }

      // add to this source's log S, log A data. We will later use these to bin
      log_area_map[this_source_key].push_back( log10(drainage_area_data_map[this_node]) );
      log_slope_map[this_source_key].push_back( log10(SA_slope[n]) );
    }
  }

  // now we bin the data
  vector<float> midpoints_output;
  vector<float> MeanX_output;
  vector<float> MedianX_output;
  vector<float> StandardDeviationX_output;
  vector<float> StandardErrorX_output;
  vector<float> MADX_output;

  vector<float> MeanY_output;
  vector<float> MinimumY_output;
  vector<float> FirstQuartileY_output;
  vector<float> MedianY_output;
  vector<float> ThirdQuartileY_output;
  vector<float> MaximumY_output;
  vector<float> StandardDeviationY_output;
  vector<float> StandardErrorY_output;
  vector<float> MADY_output;
  vector<int> number_observations_output;
  float NoDataValue = -9999;


  // these are the vectors holding all the compiled information
  vector<int> binned_basin_keys;
  vector<int> binned_source_keys;
  vector<float> binned_logA_means;
  vector<float> binned_logA_midpoints;
  vector<float> binned_logA_medians;
  vector<float> binned_logA_stdErr;
  vector<float> binned_logA_MAD;

  vector<float> binned_logS_means;
  vector<float> binned_logS_medians;
  vector<float> binned_logS_FirstQuartile;
  vector<float> binned_logS_ThirdQuartile;
  vector<float> binned_logS_StandardDeviation;
  vector<int> binnned_NObvs;

  // this holds the source numbers for each basin
  vector<int> basins_with_data;
  int last_basin = 0; // I am initializing last basin to 0 to avoid warnings. 13/11/2017 (if something break after that.) Boris
  map<int, vector<int> > basin_and_sources_map;

  // loop through all the source nodes
  map<int, vector<float> >::iterator it;
  for(it = log_area_map.begin(); it != log_area_map.end(); ++it)
  {
    this_source_key =  it->first;

    //cout << "The source key is: " << this_source_key << endl;
    // Initiate the last basin
    //last_basin = basin_key_of_this_source_map[this_source_key];

    // extract the log S-log A data for this source
    vector<float> log_area = log_area_map[this_source_key];
    vector<float> log_slope = log_slope_map[this_source_key];

    // this gets the binned data for this particular tributary
    bin_data(log_area, log_slope, log_bin_width, midpoints_output, MeanX_output,
             MedianX_output, StandardDeviationX_output, StandardErrorX_output,
             MADX_output, MeanY_output, MinimumY_output,FirstQuartileY_output,
             MedianY_output, ThirdQuartileY_output, MaximumY_output,
             StandardDeviationY_output, StandardErrorY_output, MADY_output,
             number_observations_output, NoDataValue);

    // now we need to add this information to the master vectors
    int n_bins = int(midpoints_output.size());
    int n_Obvs;
    for(int i = 0; i<n_bins; i++)
    {
      n_Obvs = number_observations_output[i];

      // only record data if there are enough observations
      if(n_Obvs > 0)
      {
        binned_basin_keys.push_back(basin_key_of_this_source_map[this_source_key]);
        binned_source_keys.push_back(this_source_key);
        binned_logA_means.push_back(MeanX_output[i]);
        binned_logA_midpoints.push_back(midpoints_output[i]);
        binned_logA_medians.push_back(MedianX_output[i]);
        binned_logA_stdErr.push_back(StandardErrorX_output[i]);
        binned_logA_MAD.push_back(MADX_output[i]);

        binned_logS_means.push_back(MeanY_output[i]);
        binned_logS_medians.push_back(MedianY_output[i]);
        binned_logS_FirstQuartile.push_back(FirstQuartileY_output[i]);
        binned_logS_ThirdQuartile.push_back(ThirdQuartileY_output[i]);
        binned_logS_StandardDeviation.push_back(StandardDeviationY_output[i]);
        binnned_NObvs.push_back(n_Obvs);

        // we need to collect a list of basins. This is rather inefficient but
        // in grand scheme of things this is far from the rate limiting step
        if (basins_with_data.size() == 0)
        {
          // This is for the first basin.
          basins_with_data.push_back(basin_key_of_this_source_map[this_source_key]);
          last_basin = basin_key_of_this_source_map[this_source_key];
        }
        else
        {
          // After the first basin, we just check to see if the basin has changed
          if (basin_key_of_this_source_map[this_source_key] !=  last_basin)
          {
            basins_with_data.push_back(basin_key_of_this_source_map[this_source_key]);
            last_basin = basin_key_of_this_source_map[this_source_key];
          }
        }
      }
    }
  }



  // open the file
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << "basin_key,median_log_A,median_log_S,logS_FirstQuartile,logS_ThirdQuartile,segment_number,segmented_log_S,segment_slope,segment_intercept,segement_R2,segment_Durbin_Watson" << endl;



  // we neeed to get the data for individual basins
  // please forgive me but this will be an extremely rudimentary and slow algoritms since
  // I want to finish soon to get beer.
  int n_tot_SA_nodes =  int(binned_basin_keys.size());
  int n_basins_in_set = int(basins_with_data.size());
  int this_basin;
  for(int i = 0; i< n_basins_in_set; i++)
  {
    this_basin = basins_with_data[i];

    // now get the data. For now we'll use the mean area and median slope since this
    // is a prototype.
    vector<float> area_data;
    vector<float> slope_data;
    vector<float> std_dev_data;
    vector<float> FirstQuartile_data;
    vector<float> ThirdQuartile_data;

    int mainstem_source = -9999;

    // loop through all the nodes collecting only nodes that are in the correct
    // basin and with the main stem (later versions will include tribs)
    for (int n = 0; n< n_tot_SA_nodes; n++)
    {
      if (binned_basin_keys[n] == this_basin)
      {

        // This looks for the first trib in basin, which is the mainstem
        if (mainstem_source == -9999)
        {
          mainstem_source = binned_source_keys[n];
        }

        if (binned_source_keys[n] == mainstem_source)
        {
          area_data.push_back(binned_logA_medians[n]);
          slope_data.push_back(binned_logS_medians[n]);
          std_dev_data.push_back(binned_logS_StandardDeviation[n]);
          FirstQuartile_data.push_back(binned_logS_FirstQuartile[n]);
          ThirdQuartile_data.push_back(binned_logS_ThirdQuartile[n]);
        }
      }
    }

    // Okay, now we partition the data
    // For now
    LSDMostLikelyPartitionsFinder Partitioner(minimum_segment_length,area_data, slope_data);

    // Partition the data
    //float sigma = 0.0005;  // this is a placeholder. Later we can use slope uncertainties. NOW USING MEASURED ERROR
    // We use the standard error of the S values as the sigma in partitioner.
    //cout << "This basin is: " << this_basin << endl;
    Partitioner.best_fit_driver_AIC_for_linear_segments(std_dev_data);

    // Now we extract all the data from the partitions
    vector<float> sigma_values;
    sigma_values.push_back(1);
    int node = 0;
    vector<float> b_values;
    vector<float> m_values;
    vector<float> r2_values;
    vector<float> DW_values;
    vector<float> fitted_y;
    vector<int> seg_lengths;
    float this_MLE;
    int this_n_segments;
    int this_n_nodes;
    float this_AIC;
    float this_AICc;

    // These are some functions that I am using to figure out what the most likeley partitioner is doing
    //Partitioner.print_to_screen_most_likeley_segment_lengths();

    Partitioner.get_data_from_best_fit_lines(node, sigma_values,
                      b_values, m_values,r2_values, DW_values, fitted_y,seg_lengths,
                      this_MLE,  this_n_segments,  this_n_nodes,
                      this_AIC,  this_AICc);



    // We need to make some new vectors for storing the relevant segment number
    vector<int> seg_number;
    vector<float> seg_m;
    vector<float> seg_b;
    vector<float> seg_r2;
    vector<float> seg_DW;
    for (int seg = 0; seg< this_n_segments; seg++)
    {
      for(int seg_node =0; seg_node< seg_lengths[seg]; seg_node++)
      {
        seg_number.push_back(seg);
        seg_m.push_back( m_values[seg]);
        seg_b.push_back( b_values[seg]);
        seg_r2.push_back( r2_values[seg]);
        seg_DW.push_back( DW_values[seg]);
      }
    }


    // Now we print the data
    for (int sn = 0; sn < int(area_data.size()); sn++)
    {
      outfile << this_basin << "," << area_data[sn] << "," << slope_data[sn] <<  ","
              << FirstQuartile_data[sn] << "," << ThirdQuartile_data[sn] << ","
              << seg_number[sn] << "," << fitted_y[sn] << ","
              << seg_m[sn] << "," << seg_b[sn] << "," << seg_r2[sn] << ","
              << seg_DW[sn] << endl;
    }
  }

  outfile.close();



}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print S-A data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_slope_area_data_to_csv(LSDFlowInfo& FlowInfo,
                                          vector<int>& SA_midpoint_node,
                                          vector<float>& SA_slope,
                                          string filename)
{
  // open the data file
  int row,col;
  LSDCoordinateConverterLLandUTM Converter;
  int n_nodes = int(SA_midpoint_node.size());
  double latitude,longitude;
  ofstream  SA_out;
  SA_out.open(filename.c_str());
  cout << "Opening the data file: " << filename << endl;
  SA_out << "latitude,longitude,chi,elevation,flow_distance,drainage_area,slope,source_key,basin_key" << endl;
  int this_node;
  if (n_nodes <= 0)
  {
    cout << "Trying to print SA data but there doesn't seem to be any." << endl;
    cout << "Have you run the automator and gathered the sources yet?" << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = SA_midpoint_node[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      SA_out.precision(9);
      SA_out << latitude << ","
             << longitude << ",";
      SA_out.precision(5);
      SA_out << chi_data_map[this_node] << ","
             << elev_data_map[this_node] << ","
             << flow_distance_data_map[this_node] << ","
             << drainage_area_data_map[this_node] << ","
             << SA_slope[n] << ","
             << source_keys_map[this_node] << ","
             << baselevel_keys_map[this_node];
      SA_out << endl;
    }
  }

  SA_out.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print chi maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_chi_data_map_to_csv(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,chi,elevation,flow_distance,drainage_area,source_key,basin_key" << endl;
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print chi maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_chi_data_map_to_csv_for_single_basin(LSDFlowInfo& FlowInfo, string filename, int basin_key)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,chi,elevation,flow_distance,drainage_area,source_key,basin_key" << endl;
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];

      if (baselevel_keys_map[this_node] == basin_key)
      {

        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);

        chi_data_out.precision(9);
        chi_data_out << latitude << ","
                   << longitude << ",";
        chi_data_out.precision(5);
        chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

        chi_data_out << endl;
      }
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_full(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row, col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // test to see if there is segment numbering
  bool have_segments = false;
  if( segment_counter_map.size() == node_sequence.size())
  {
    have_segments = true;
  }

  // test to see if the fitted elevations have been calculated
  bool have_segmented_elevation = false;
  if( segmented_elevation_map.size() == node_sequence.size())
  {
    have_segmented_elevation = true;
  }


  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "node,row,col,latitude,longitude,chi,elevation,flow_distance,drainage_area,m_chi,b_chi,source_key,basin_key";
  if(have_segmented_elevation)
  {
    chi_data_out << ",segmented_elevation";
  }
  if (have_segments)
  {
    chi_data_out << ",segment_number";
    cout << "I added the segment number in the csv file"<< endl;
  }
  chi_data_out << endl;




  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out << this_node << ","
                   << row << ","
                   << col << ",";
      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << ","
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      if(have_segmented_elevation)
      {
        chi_data_out << "," << segmented_elevation_map[this_node];
      }
      if (have_segments)
      {
        chi_data_out << "," << segment_counter_map[this_node];
      }
      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_intersources_mchi_map( string filename)
{

  // open the data file
  ofstream  file_out_;
  file_out_.open(filename.c_str());
  file_out_ << "source_key,receiving_source_key,m_chi,chi";
  file_out_ << endl;

 
  for (map<int,int>::iterator alpaca = map_source_key_receiver.begin(); alpaca != map_source_key_receiver.end() ; alpaca++)
  {
    file_out_.precision(5);
    file_out_ << alpaca->first << ","
                 << alpaca->second << ","
                 << map_source_key_receiver_mchi[alpaca->first]
                 << M_chi_data_map[alpaca->first]
                 << chi_data_map[alpaca->first];
    file_out_ << endl;
  }
  

  file_out_.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Development function to Print data maps to file including knickpoints
// BG
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_full_knickpoints(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = (node_sequence.size());

  // test to see if there is segment numbering
  bool have_segments = false;
  if( segment_counter_map.size() == node_sequence.size())
  {
    have_segments = true;
  }

  // test to see if the fitted elevations have been calculated
  bool have_segmented_elevation = false;
  if( segmented_elevation_map.size() == node_sequence.size())
  {
    have_segmented_elevation = true;
  }


  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,chi,elevation,flow_distance,drainage_area,m_chi,b_chi,source_key,basin_key";
  if(have_segmented_elevation)
  {
    chi_data_out << ",segmented_elevation";
  }
  if (have_segments)
  {
    chi_data_out << ",segment_number";
  }
  chi_data_out << ",knickpoints,knickpoint_sign,segment_length"; // add the knickpoint col
  chi_data_out << endl;

  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << ","
                   << source_keys_map[this_node] << ","
                   << baselevel_keys_map[this_node];

      if(have_segmented_elevation)
      {
        chi_data_out << "," << segmented_elevation_map[this_node];
      }
      if (have_segments)
      {
        chi_data_out << "," << segment_counter_map[this_node];
      }

      chi_data_out << "," << segment_counter_knickpoint_map[this_node];
      chi_data_out << "," << segment_knickpoint_sign_map[this_node];
      chi_data_out << "," << segment_length_map[this_node];
      chi_data_out << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_source_keys(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col, key;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  map<int, int>::iterator it;

  // open the data file
  ofstream  source_keys_out;
  source_keys_out.open(filename.c_str());
  source_keys_out << "latitude,longitude,source_node,source_key" << endl;

  // loop through the source key map
  for ( it = key_to_source_map.begin(); it != key_to_source_map.end(); it++ )
  {
    key = it->second;
    this_node = it->first;
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    source_keys_out.precision(9);
    source_keys_out << latitude << ","
                   << longitude << "," << this_node << ",";
    source_keys_out.precision(5);
    source_keys_out << key << endl;
  }

  source_keys_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_baselevel_keys(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, this_junc,row,col, key;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  map<int, int>::iterator it;

  // open the data file
  ofstream  baselevel_keys_out;
  baselevel_keys_out.open(filename.c_str());
  baselevel_keys_out << "latitude,longitude,baselevel_node,baselevel_junction,baselevel_key" << endl;

  // loop through the source
  for ( it = key_to_baselevel_map.begin(); it != key_to_baselevel_map.end(); it++ )
  {
    key = it->second;
    this_node = it->first;
    this_junc = JunctionNetwork.get_Junction_of_Node(this_node, FlowInfo);
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    baselevel_keys_out.precision(9);
    baselevel_keys_out << latitude << ","
                   << longitude << "," << this_node << ",";
    baselevel_keys_out.precision(5);
    baselevel_keys_out << this_junc <<"," << key << endl;
  }

  baselevel_keys_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints the basins and an additional file that has basin centroids
// and labelling information for plotting
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_basins(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& JunctionNetwork,
                               vector<int> Junctions, string base_filename)
{
  int N_Juncs = Junctions.size();
  LSDCoordinateConverterLLandUTM Converter;


  // Get some data members for holding basins and the raster
  vector<LSDBasin> AllTheBasins;
  map<int,int> drainage_of_other_basins;
  LSDIndexRaster BasinMasterRaster;

  string basin_raster_name = base_filename+"_AllBasins";
  string basin_info_name = base_filename+"_AllBasinsInfo.csv";

  ofstream basin_info_out;
  basin_info_out.open(basin_info_name.c_str());
  basin_info_out << "latitude,longitude,outlet_latitude,outlet_longitude,outlet_junction,basin_key" << endl;

  // Make sure the full lat-long information is printed
  basin_info_out.precision(9);

  // These store row and column information for converting the outlet and centroid to
  // latitude and longitude
  int centroid_i, centroid_j, outlet_i, outlet_j;
  double out_lat,out_long, cen_lat, cen_long;

  //cout << "I am trying to print basins, found " << N_BaseLevelJuncs << " base levels." << endl;
  // Loop through the junctions
  for(int BN = 0; BN<N_Juncs; BN++)
  {
    //cout << "Getting basin " << BN << " and the junction is: "  << Junctions[BN] << endl;
    LSDBasin thisBasin(Junctions[BN],FlowInfo, JunctionNetwork);
    //cout << "...got it!" << endl;
    AllTheBasins.push_back(thisBasin);

    // This is required if the basins are nested--test the code which numbers
    // to be overwritten by a smaller basin
    drainage_of_other_basins[Junctions[BN]] = thisBasin.get_NumberOfCells();


    // We need to see if we can find the basin key
    int basin_key = -9999;

    // need to node index of this junction
    int node_of_junction =  JunctionNetwork.get_Node_of_Junction( Junctions[BN] );
    if ( key_to_baselevel_map.find( node_of_junction) != key_to_baselevel_map.end() )
    {
      basin_key = key_to_baselevel_map[node_of_junction];
    }

    // get the centroid and outlet locations
    centroid_i = thisBasin.get_Centroid_i();
    centroid_j = thisBasin.get_Centroid_j();

    outlet_i = thisBasin.get_Outlet_i();
    outlet_j = thisBasin.get_Outlet_j();

    // Find the latitude and longitude of the outlet and centroid
    get_lat_and_long_locations(centroid_i, centroid_j, cen_lat, cen_long, Converter);
    get_lat_and_long_locations(outlet_i, outlet_j, out_lat, out_long, Converter);

    basin_info_out << cen_lat << "," << cen_long << "," << out_lat << "," << out_long << "," << Junctions[BN] << "," << basin_key << endl;
  }
  basin_info_out.close();

  // now loop through everything again getting the raster
  if (N_Juncs > 0)     // this gets the first raster
  {
    BasinMasterRaster = AllTheBasins[0].write_integer_data_to_LSDIndexRaster(Junctions[0], FlowInfo);
  }

  // now add on the subsequent basins
  for(int BN = 1; BN<N_Juncs; BN++)
  {
    AllTheBasins[BN].add_basin_to_LSDIndexRaster(BasinMasterRaster, FlowInfo,
                              drainage_of_other_basins, Junctions[BN]);
  }


  // We need to use bil format since there is georeferencing
  string raster_ext = "bil";
  // print the basin raster
  BasinMasterRaster.write_raster(basin_raster_name, raster_ext);
  cout << "Finished with exporting basins!" << endl;

}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_basic(LSDFlowInfo& FlowInfo, string filename)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;

  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,m_chi,b_chi" << endl;

  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter);

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(6);
      chi_data_out << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << "," << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

///@brief Return a map of useful maps
///@author B.G.
///@date 21/11/2018
map<string, map<int,float> > LSDChiTools::get_data_maps()
{
  map<string, map<int,float> > output;
  output["flow_distance"]=flow_distance_data_map;
  output["elevation"] = elev_data_map;
  output["chi"] = chi_data_map;
  output["DA"] = drainage_area_data_map;
  output["m_chi"] = M_chi_data_map;

  map<int,float> tempmapsk;
  for(map<int,int>::iterator it= source_keys_map.begin(); it!= source_keys_map.end(); it++)
  {
    tempmapsk[it->first]=float(it->second);
  }
  
  output["SK"] = tempmapsk;

  return output;
}


///@brief return vectors of integer data calculated by Mudd et al., 2014 JGR.
///@author B.G.
///@date 21/11/2018
map<string, vector<int> > LSDChiTools::get_integer_vecdata_for_m_chi(LSDFlowInfo &Flowinfo)
{
  size_t size = node_sequence.size();
  vector<int> arrnodeID(size), arrrow(size), arrcol(size), arrsource_key(size), arrbasin_key(size);
  for (size_t n=0; n<size; n++)
  {
    int this_node = node_sequence[n];
    arrnodeID[n] = this_node;
    int trow,tcol;
    Flowinfo.retrieve_current_row_and_col(this_node,trow,tcol);
    arrrow[n] = trow; arrcol[n] = tcol;
    arrsource_key[n] =  source_keys_map[this_node];
    arrbasin_key[n] =  baselevel_keys_map[this_node];
  }

  map<string, vector<int> > output;
  output["nodeID"] = arrnodeID;
  output["row"] = arrrow;
  output["col"] = arrcol;
  output["source_key"] = arrsource_key;
  output["basin_key"] = arrbasin_key;

  return output;

}

///@brief return vectors of integer data calculated by Mudd et al., 2014 JGR.
///@author B.G.
///@date 21/11/2018
map<string, vector<int> > LSDChiTools::get_integer_vecdata_for_knickpoint_analysis(LSDFlowInfo &Flowinfo)
{

  vector<int> arrnodeID, arrrow, arrcol, arrsource_key, arrbasin_key, arrsign;
  int cpt = -1, this_node = 0, row,col;
  map<int,bool> is_done;
  map<int,float>::iterator iter;
  for (iter = ksn_kp_map.begin(); iter != ksn_kp_map.end(); iter++)
  {
    this_node = iter->first;
    float this_kp = iter->second;
    int nearnode = nearest_node_centroid_kp[this_node];
    // int nearnode_stepped = nearest_node_centroid_kp_stepped[this_node];    // not used (SMM)
    float this_segelev = 0;

    if(kp_segdrop.count(nearnode) == 1)
    {
      is_done[nearnode] = true;
      this_segelev = kp_segdrop[nearnode];
    }
    // get the centroid location
      // float this_x = ksn_centroid[this_node].first, this_y = ksn_centroid[this_node].second; // BUGGY AT THE MOMENT
    float this_x = 0, this_y =0;
    Flowinfo.retrieve_current_row_and_col(nearnode,row,col);
    get_x_and_y_locations(row, col, this_x, this_y);

     // Just adding sign column for plotting purposes
    int this_sign = 0;
    if(this_kp>0){this_sign = 1;}else{this_sign=(-1);}

    arrnodeID.push_back(nearnode);
    arrrow.push_back(row); 
    arrcol.push_back(col);
    arrsource_key.push_back(source_keys_map[nearnode]);
    arrbasin_key.push_back(baselevel_keys_map[nearnode]);
    arrsign.push_back(this_sign);    
  }

  // Now implementing the lonely step knickpoints
  for (iter = kp_segdrop.begin(); iter != kp_segdrop.end(); iter++)
  {
    cpt++;
    this_node = iter->first;
    float this_segelev = iter->second;
    // int nearnode = nearest_node_centroid_kp_stepped[this_node];
    float this_kp = 0; // All the ksn knickpoints have been written, these one only have a segelev component

    // if(chi_data_map[this_node] == 0)
    // {
    //   cout << "This node is screwed" <<endl;
    // }
    // if(chi_data_map[nearnode] == 0)
    // {
    //   cout << "nearnode is screwed" <<endl ;
    // }

    if(is_done.count(this_node) != 1 && chi_data_map[this_node] != 0)
    {
      // cout << "Tbg 45b" << endl;
      // get the centroid location
      // float this_x = ksn_centroid[this_node].first, this_y = ksn_centroid[this_node].second; // BUGGY AT THE MOMENT
      float this_x = 0, this_y =0;
      Flowinfo.retrieve_current_row_and_col(this_node,row,col);
      get_x_and_y_locations(row, col, this_x, this_y);
      arrnodeID.push_back(this_node);
      arrrow.push_back(row);        
      arrcol.push_back(col);
      arrsource_key.push_back(source_keys_map[this_node]);
      arrbasin_key.push_back(baselevel_keys_map[this_node]);
      arrsign.push_back(1);
    }
  }



  map<string, vector<int> > output;
  output["nodeID"] = arrnodeID;
  output["row"] = arrrow;
  output["col"] = arrcol;
  output["source_key"] = arrsource_key;
  output["basin_key"] = arrbasin_key;
  output["sign"] = arrsign;

  return output;

}

///@brief return vectors of integer data calculated by Mudd et al., 2014 JGR.
///@author B.G.
///@date 21/11/2018
map<string, vector<float> > LSDChiTools::get_float_vecdata_for_knickpoint_analysis(LSDFlowInfo &Flowinfo)
{

  vector<int> arrnodeID;
  vector<float> this_delta, this_x_coord, this_y_coord, this_elevation, this_drainage_area, this_flow_distance, this_chi, this_segelev;

  int cpt = -1, this_node = 0, row,col;
  map<int,bool> is_done;
  map<int,float>::iterator iter;
  for (iter = ksn_kp_map.begin(); iter != ksn_kp_map.end(); iter++)
  {
    this_node = iter->first;
    float this_kp = iter->second;
    int nearnode = nearest_node_centroid_kp[this_node];
    // int nearnode_stepped = nearest_node_centroid_kp_stepped[this_node];    // not used (SMM)
    float this_segelev_val = 0;

    if(kp_segdrop.count(nearnode) == 1)
    {
      is_done[nearnode] = true;
      this_segelev_val = kp_segdrop[nearnode];
    }
    // get the centroid location
      // float this_x = ksn_centroid[this_node].first, this_y = ksn_centroid[this_node].second; // BUGGY AT THE MOMENT
    float this_x = 0, this_y =0;
    Flowinfo.retrieve_current_row_and_col(nearnode,row,col);
    get_x_and_y_locations(row, col, this_x, this_y);

     // Just adding sign column for plotting purposes
    int this_sign = 0;
    if(this_kp>0){this_sign = 1;}else{this_sign=(-1);}

    this_delta.push_back(this_kp);
    this_x_coord.push_back(this_x);
    this_y_coord.push_back(this_y);
    this_elevation.push_back(elev_data_map[nearnode]);
    this_drainage_area.push_back(drainage_area_data_map[nearnode]);
    this_flow_distance.push_back(flow_distance_data_map[nearnode]);
    this_chi.push_back(chi_data_map[nearnode]);
    this_segelev.push_back(this_segelev_val);
  }

  // Now implementing the lonely step knickpoints
  for (iter = kp_segdrop.begin(); iter != kp_segdrop.end(); iter++)
  {
    cpt++;
    this_node = iter->first;
    float this_segelev_val = iter->second;
    // int nearnode = nearest_node_centroid_kp_stepped[this_node];
    float this_kp = 0; // All the ksn knickpoints have been written, these one only have a segelev component

    // if(chi_data_map[this_node] == 0)
    // {
    //   cout << "This node is screwed" <<endl;
    // }
    // if(chi_data_map[nearnode] == 0)
    // {
    //   cout << "nearnode is screwed" <<endl ;
    // }

    if(is_done.count(this_node) != 1 && chi_data_map[this_node] != 0)
    {
      // cout << "Tbg 45b" << endl;
      // get the centroid location
      // float this_x = ksn_centroid[this_node].first, this_y = ksn_centroid[this_node].second; // BUGGY AT THE MOMENT
      float this_x = 0, this_y =0;
      Flowinfo.retrieve_current_row_and_col(this_node,row,col);
      get_x_and_y_locations(row, col, this_x, this_y);
      this_delta.push_back(this_kp);
      this_x_coord.push_back(this_x);
      this_y_coord.push_back(this_y);
      this_elevation.push_back(elev_data_map[this_node]);
      this_drainage_area.push_back(drainage_area_data_map[this_node]);
      this_flow_distance.push_back(flow_distance_data_map[this_node]);
      this_chi.push_back(chi_data_map[this_node]);
      this_segelev.push_back(this_segelev_val);
    }
  }



  map<string, vector<float> > output;
  output["delta_ksn"] = this_delta;
  output["x"] = this_x_coord;
  output["y"] = this_y_coord;
  output["elevation"] = this_elevation;
  output["flow_distance"] = this_flow_distance;
  output["drainage_area"] = this_drainage_area;
  output["chi"] = this_chi;
  output["delta_zseg"] = this_segelev;


  return output;

}

//     chi_data_out << elev_data_map[this_node] << ","
//                  << flow_distance_data_map[this_node] << ","
//                  << chi_data_map[this_node] << ","
//                  << drainage_area_data_map[this_node] << ","
//                  << ksn_diff_knickpoint_map[this_node] << ","
//                  << ksn_ratio_knickpoint_map[this_node] << ","
//                  << ksn_sign_knickpoint_map[this_node] << ","
//                  << ksn_rad_knickpoint_map[this_node] << ","
//                  << ksn_cumul_knickpoint_map[this_node] << ","
//                  << rksn_cumul_knickpoint_map[this_node] << ","
//                  << rad_cumul_knickpoint_map[this_node] << ","
//                  << source_keys_map[this_node] << ","
//                  << baselevel_keys_map[this_node];
//     chi_data_out << endl;

// ///@brief return vectors of integer data calculated by Mudd et al., 2014 JGR.
// ///@author B.G.
// ///@date 21/11/2018
// map<string, vector<float> > LSDChiTools::get_float_vecdata_for_knickpoint_analysis(LSDFlowInfo &Flowinfo)
// {
//   // Overall size
//   size_t size = ksn_diff_knickpoint_map.size() + ;
//   // Getting the data ready
//   vector<int> arrnodeID(size);
//   vector<float> this_delta(size), this_ratio_ksn(size), this_x_coord(size), this_y_coord(size), this_elevation(size), this_drainage_area(size), this_flow_distance(size), this_chi(size), this_segelev(size);

//   size_T cpt = -1;
//   for (iter = ksn_diff_knickpoint_map.begin(); iter != ksn_diff_knickpoint_map.end(); iter++)
//   {
//     cpt++;
//     // getting node info
//     int this_node = node_sequence[n];
//     arrnodeID[n] = this_node;
//     int trow,tcol;
//     float this_x, this_y;
//     Flowinfo.retrieve_current_row_and_col(this_node,trow,tcol);
//     Flowinfo.get_x_and_y_from_current_node(this_node, this_x, this_y);

//     // feeding the vectors
//     this_delta[n] = ksn_diff_knickpoint_map[this_node];
//     this_ratio_ksn[n] = ksn_ratio_knickpoint_map[this_node];
//     this_segelev[n]
//     this_x_coord[n] = this_x;
//     this_y_coord[n] = this_y;
//     this_elevation[n] = elev_data_map[this_node];
//     this_drainage_area[n] = drainage_area_data_map[this_node];
//     this_flow_distance[n] = flow_distance_data_map[this_node];
//     this_chi[n] = chi_data_map[this_node];


//   }

//   map<string, vector<float> > output;
//   output["m_chi"] = this_m_chi;
//   output["b_chi"] = this_b_chi;
//   output["x"] = this_x_coord;
//   output["y"] = this_y_coord;
//   output["elevation"] = this_elevation;
//   output["drainage_area"] = this_drainage_area;
//   output["flow_distance"] = this_flow_distance;
//   output["chi"] = this_chi;

//   return output;

// }


///@brief return vectors of integer data calculated by Mudd et al., 2014 JGR.
///@author B.G.
///@date 21/11/2018
map<string, vector<float> > LSDChiTools::get_float_vecdata_for_m_chi(LSDFlowInfo &Flowinfo)
{
  // Overall size
  size_t size = node_sequence.size();
  // Getting the data ready
  vector<int> arrnodeID(size);
  vector<float> this_m_chi(size), this_b_chi(size), this_x_coord(size), this_y_coord(size), this_elevation(size), this_drainage_area(size), this_flow_distance(size), this_chi(size), this_segelev(size);

  // Looping through each nodes
  for (size_t n=0; n<size; n++)
  {
    // getting node info
    int this_node = node_sequence[n];
    arrnodeID[n] = this_node;
    int trow,tcol;
    float this_x, this_y;
    Flowinfo.retrieve_current_row_and_col(this_node,trow,tcol);
    Flowinfo.get_x_and_y_from_current_node(this_node, this_x, this_y);

    // feeding the vectors
    this_m_chi[n] = M_chi_data_map[this_node];
    this_b_chi[n] = b_chi_data_map[this_node];
    this_x_coord[n] = this_x;
    this_y_coord[n] = this_y;
    this_elevation[n] = elev_data_map[this_node];
    this_drainage_area[n] = drainage_area_data_map[this_node];
    this_flow_distance[n] = flow_distance_data_map[this_node];
    this_chi[n] = chi_data_map[this_node];
    this_segelev[n] = segmented_elevation_map[this_node];

  }

  map<string, vector<float> > output;
  output["m_chi"] = this_m_chi;
  output["b_chi"] = this_b_chi;
  output["x"] = this_x_coord;
  output["y"] = this_y_coord;
  output["elevation"] = this_elevation;
  output["drainage_area"] = this_drainage_area;
  output["flow_distance"] = this_flow_distance;
  output["chi"] = this_chi;
  output["segmented_elevation"] = this_segelev;

  return output;

}

///@brief Return the node sequence
///@author B.G.
///@date 21/11/2018
vector<int> LSDChiTools::get_vectors_of_node()
{
  vector<int> output;
  output = node_sequence;
  return output;
}

void LSDChiTools::update_chi_vector_for_opti_disorder_with_uncert(vector<int>& sorted_nodes, vector<float>& that_chi)
{
  for(size_t i=0;i<sorted_nodes.size();i++)
  {
    int that_node = sorted_nodes[i];
    float tthat_chi;
    if(chi_data_map.count(that_node) == 0)
      tthat_chi = 0;
    else
      tthat_chi = chi_data_map[that_node];
    that_chi[i] = tthat_chi;
  }
}

void LSDChiTools::update_elevation_vector_for_opti_disorder_with_uncert(vector<int>& sorted_nodes, vector<float>& that_elev)
{
  for(size_t i=0;i<sorted_nodes.size();i++)
  {
    int that_node = sorted_nodes[i];
    float tthat_elev;
    if(elev_data_map.count(that_node) == 0)
      tthat_elev = 0;
    else
      tthat_elev = elev_data_map[that_node];

    that_elev[i] = tthat_elev;
  }

}

#endif
