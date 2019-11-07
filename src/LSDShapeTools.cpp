//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDShapeTools
// Land Surface Dynamics Shapefile tools
//
// A collection of routines for maipulating the binary ESRI shapefile format
// for use within the Edinburgh Land Surface Dynamics group topographic toolbox
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

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------


#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <functional> // For string splitter
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>
#include <utility>
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
using namespace std;

#ifndef ShapeTools_CPP
#define ShapeTools_CPP





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to test the Byte order of the system.
// Returns a boolean value where true is little endian.
//
// SWDG 11/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
bool SystemEndiannessTest(){

  int TestInt = 1;  //this is stored as 4 bytes in memory
  int ReconstructedTestInt;
  int ReconstructedTestIntSwapped;

  char * TestBytes = (char *) &TestInt;  //convert each byte of the int into a char

  memcpy(&ReconstructedTestInt, TestBytes, 4);

  BYTE temp = TestBytes[0];
  TestBytes[0] = TestBytes[3];
  TestBytes[3] = temp;

  temp = TestBytes[1];
  TestBytes[1] = TestBytes[2];
  TestBytes[2] = temp;

  memcpy(&ReconstructedTestIntSwapped, TestBytes, 4);

  if (ReconstructedTestInt == 1){
    //cout << "Little Endian" << endl;
    return true;
  }
  else if (ReconstructedTestIntSwapped == 1){
    //cout << "Big Endian" << endl;
    return false;
  }
  else{
    cout << "Unable to determine endianness of system." << endl;
    exit(EXIT_FAILURE);
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to swap the byte order of a word in memory. Used if the system's byte order
// does not match the data in the shapefile.
//
// SWDG 11/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void ByteSwap(int length, void * ByteData){

    BYTE temp;

    for (int i = 0; i< length/2; ++i){

      temp = ((BYTE *) ByteData)[i];
      ((BYTE *) ByteData)[i] = ((BYTE *) ByteData)[length-i-1];
      ((BYTE *) ByteData)[length-i-1] = temp;

    }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to get the size of the binary file being loaded.
//
// Taken from http://www.dreamincode.net/forums/topic/170054-understanding-and-reading-binary-files-in-c/
//
// SWDG 10/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
long getFileSize(FILE *file){
  long lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, 2);
  lEndPos = ftell(file);
  fseek(file, lCurPos, 0);
  return lEndPos;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to load an ESRI ShapeFile.
//
// Only works for X,Y point shapefiles at present and it's behaviour is totally undefined
// if you pass in any other type of file.
//
// In future this will be rebuilt into a full class that can support shapefiles of
// different types.
//
// Built in part from:
// http://www.dreamincode.net/forums/topic/170054-understanding-and-reading-binary-files-in-c/
//
// SWDG 13/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
PointData LoadShapefile(string Filename){

  BYTE *ByteData;      // Pointer to our buffered data
  FILE *file = NULL;    // File pointer

  // Open the file in binary mode using the "rb" format string
  // This also checks if the file exists and/or can be opened for reading correctly
  if ((file = fopen(Filename.c_str(), "rb")) == NULL){
    cout << "Could not open specified file" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "File opened successfully" << endl;
  }

  // Get the size of the file in bytes
  long fileSize = getFileSize(file);

  // Allocate space in the buffer for the whole file
  ByteData = new BYTE[fileSize];

  // Read the file in to the buffer
  fread(ByteData, fileSize, 1, file);

  //Declare variables used in the method
  int FileLength;
  int ShapeType;
  double Xmin;
  double Ymin;
  double Xmax;
  double Ymax;
  double Zmin;
  double Zmax;
  double Mmin;
  double Mmax;
  int RecordLength;
  int NoOfRecords;
  PointData Points;
  double TempX;
  double TempY;

  bool LittleEndian = SystemEndiannessTest(); // Get byteorder of the system

  // If system byte order is Little Endian
  if (LittleEndian == true){

    // Get the length of the file
    ByteSwap(4, ByteData+24);
    memcpy(&FileLength, ByteData+24, 4);

    // Get type of shape in file (not currently used) see
    memcpy(&ShapeType, ByteData+32, 4);

    // Get Georeferencing data (not currently used)
    memcpy(&Xmin, ByteData+36, 8);
    memcpy(&Ymin, ByteData+44, 8);
    memcpy(&Xmax, ByteData+52, 8);
    memcpy(&Ymax, ByteData+60, 8);
    memcpy(&Zmin, ByteData+68, 8);
    memcpy(&Zmax, ByteData+76, 8);
    memcpy(&Mmin, ByteData+84, 8);
    memcpy(&Mmax, ByteData+92, 8);

    // Get the length of the first record (All records are the same length for points)
    ByteSwap(4, ByteData+104);
    memcpy(&RecordLength, ByteData+104, 4);

    //Calculate the number of records in the file
    NoOfRecords = (FileLength-50)/(RecordLength+4);  // FileLength - 50(the length in words of the header)
                                                     // divided by RecordLength+4 (4 is the length in words of the record header)

    if (NoOfRecords == 0){
      cout << "Empty Shapefile. No Data to read!\n" << endl;
      exit(EXIT_FAILURE);
    }

    //Read all of the records into the Points structure
    for (int q = 0; q < NoOfRecords; ++q){

      memcpy(&TempX, ByteData+112+(q * ((RecordLength+4)*2)), 8); // RecordLength+4*2 == 28 (The total length of a record)
      memcpy(&TempY, ByteData+120+(q * ((RecordLength+4)*2)), 8);

      Points.X.push_back(TempX);
      Points.Y.push_back(TempY);

    }

  }
  // If system byte order is Big Endian
  else{

    // Get the length of the file
    memcpy(&FileLength, ByteData+24, 4);

    // Get type of shape in file (not currently used) see
    ByteSwap(4, ByteData+32);
    memcpy(&ShapeType, ByteData+32, 4);

    // Get Georeferencing data (not currently used)
    ByteSwap(8, ByteData+36);
    memcpy(&Xmin, ByteData+36, 8);
    ByteSwap(8, ByteData+44);
    memcpy(&Ymin, ByteData+44, 8);
    ByteSwap(8, ByteData+52);
    memcpy(&Xmax, ByteData+52, 8);
    ByteSwap(8, ByteData+60);
    memcpy(&Ymax, ByteData+60, 8);
    ByteSwap(8, ByteData+68);
    memcpy(&Zmin, ByteData+68, 8);
    ByteSwap(8, ByteData+76);
    memcpy(&Zmax, ByteData+76, 8);
    ByteSwap(8, ByteData+84);
    memcpy(&Mmin, ByteData+84, 8);
    ByteSwap(8, ByteData+92);
    memcpy(&Mmax, ByteData+92, 8);

    // Get the length of the first record (All records are the same length for points)
    memcpy(&RecordLength, ByteData+104, 4);

    //Calculate the number of records in the file
    NoOfRecords = (FileLength-50)/(RecordLength+4);  // FileLength - 50(the length in words of the header)
                                                     // divided by RecordLength+4 (4 is the length in words of the record header)

    if (NoOfRecords == 0){
      cout << "Empty Shapefile. No Data to read!\n" << endl;
      exit(EXIT_FAILURE);
    }

    //Read all of the records into the Points structure
    for (int q = 0; q < NoOfRecords; ++q){

      ByteSwap(8, ByteData+112+(q * ((RecordLength+4)*2)));
      memcpy(&TempX, ByteData+112+(q * ((RecordLength+4)*2)), 8); // RecordLength+4*2 == 28 (The total length of a record)
      ByteSwap(8, ByteData+120+(q * ((RecordLength+4)*2)));
      memcpy(&TempY, ByteData+120+(q * ((RecordLength+4)*2)), 8);

      Points.X.push_back(TempX);
      Points.Y.push_back(TempY);

    }

  }

  // Close the file and return the point data
  fclose(file);
  return Points;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to load an ESRI polyline Shapefile.
//
// Only works for polyline shapefiles at present and it's behaviour is totally undefined
// if you pass in any other type of file.
//
// In future this will be rebuilt into a full class that can support shapefiles of
// different types.
//
// Returns a vector of points. So that each item in the vector represents a single polyline.
//
// Built in part from:
// http://www.dreamincode.net/forums/topic/170054-understanding-and-reading-binary-files-in-c/
//
// SWDG 17/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<PointData> LoadPolyline(string Filename){

  BYTE *ByteData;      // Pointer to our buffered data
  FILE *file = NULL;    // File pointer

  // Open the file in binary mode using the "rb" format string
  // This also checks if the file exists and/or can be opened for reading correctly
  if ((file = fopen(Filename.c_str(), "rb")) == NULL){
    cout << "Could not open specified file" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "File opened successfully" << endl;
  }

  // Get the size of the file in bytes
  long fileSize = getFileSize(file);

  // Allocate space in the buffer for the whole file
  ByteData = new BYTE[fileSize];

  // Read the file in to the buffer
  fread(ByteData, fileSize, 1, file);

  //Declare variables used in the method
  int FileLength;
  int ShapeType;
  double Xmin;
  double Ymin;
  double Xmax;
  double Ymax;
  double Zmin;
  double Zmax;
  double Mmin;
  double Mmax;
  int RecordLength;
  PointData Points;
  double TempX;
  double TempY;
  vector<PointData> Polylines;
  int Skip;
  int shapetype;
  int numparts;
  int numpoints;

  bool LittleEndian = SystemEndiannessTest(); // Get byteorder of the system

  // If system byte order is Little Endian
  if (LittleEndian == true){

    // Get the length of the file
    ByteSwap(4, ByteData+24);
    memcpy(&FileLength, ByteData+24, 4);

    if (FileLength == 50){
      cout << "Empty Shapefile. No Data to read!\n" << endl;
      exit(EXIT_FAILURE);
    }

    // Get type of shape in file (not currently used) see
    memcpy(&ShapeType, ByteData+32, 4);

    // Get Georeferencing data (not currently used)
    memcpy(&Xmin, ByteData+36, 8);
    memcpy(&Ymin, ByteData+44, 8);
    memcpy(&Xmax, ByteData+52, 8);
    memcpy(&Ymax, ByteData+60, 8);
    memcpy(&Zmin, ByteData+68, 8);
    memcpy(&Zmax, ByteData+76, 8);
    memcpy(&Mmin, ByteData+84, 8);
    memcpy(&Mmax, ByteData+92, 8);

    // Get the length of the first record
    ByteSwap(4, ByteData+104);
    memcpy(&RecordLength, ByteData+104, 4);
    Skip = RecordLength*2 + 8 + 100;

    memcpy(&shapetype, ByteData+108, 4);
    memcpy(&numparts, ByteData+144, 4);
    memcpy(&numpoints, ByteData+148, 4);  // number of XY pairs in the first polyline

    for (int q = 0; q < numpoints; ++q){

      memcpy(&TempX, ByteData+156+(q*16), 8);
      memcpy(&TempY, ByteData+164+(q*16), 8);

      Points.X.push_back(TempX);
      Points.Y.push_back(TempY);

    }

    Polylines.push_back(Points);

    while (Skip < (FileLength*2)){

      ByteSwap(4, ByteData+Skip+4);
      memcpy(&RecordLength, ByteData+Skip+4, 4);
      memcpy(&shapetype, ByteData+Skip+8, 4);
      memcpy(&numparts, ByteData+Skip+44, 4);
      memcpy(&numpoints, ByteData+Skip+48, 4);  // number of XY pairs in the polyline

      PointData Points;
      for (int w = 0; w < numpoints; ++w){

        memcpy(&TempX, ByteData+Skip+56+(w*16), 8);
        memcpy(&TempY, ByteData+Skip+64+(w*16), 8);

        Points.X.push_back(TempX);
        Points.Y.push_back(TempY);

      }

      Polylines.push_back(Points);
      Skip += (RecordLength*2) +8;

    }

  }
  // If system byte order is Big Endian
  else{

    // Get the length of the file
    memcpy(&FileLength, ByteData+24, 4);

    if (FileLength == 50){
      cout << "Empty Shapefile. No Data to read!\n" << endl;
      exit(EXIT_FAILURE);
    }

    // Get type of shape in file (not currently used) see
    ByteSwap(4, ByteData+32);
    memcpy(&ShapeType, ByteData+32, 4);

    // Get Georeferencing data (not currently used)

    ByteSwap(8, ByteData+36);
    memcpy(&Xmin, ByteData+36, 8);
    ByteSwap(8, ByteData+44);
    memcpy(&Ymin, ByteData+44, 8);
    ByteSwap(8, ByteData+52);
    memcpy(&Xmax, ByteData+52, 8);
    ByteSwap(8, ByteData+60);
    memcpy(&Ymax, ByteData+60, 8);
    ByteSwap(8, ByteData+68);
    memcpy(&Zmin, ByteData+68, 8);
    ByteSwap(8, ByteData+76);
    memcpy(&Zmax, ByteData+76, 8);
    ByteSwap(8, ByteData+84);
    memcpy(&Mmin, ByteData+84, 8);
    ByteSwap(8, ByteData+92);
    memcpy(&Mmax, ByteData+92, 8);

    // Get the length of the first record

    memcpy(&RecordLength, ByteData+104, 4);
    Skip = RecordLength*2 + 8 + 100;

    ByteSwap(8, ByteData+108);
    memcpy(&shapetype, ByteData+108, 4);
    ByteSwap(8, ByteData+144);
    memcpy(&numparts, ByteData+144, 4);
    ByteSwap(8, ByteData+148);
    memcpy(&numpoints, ByteData+148, 4);  // number of XY pairs in the first polyline

    for (int q = 0; q < numpoints; ++q){

      ByteSwap(8, ByteData+156+(q*16));
      memcpy(&TempX, ByteData+156+(q*16), 8);
      ByteSwap(8, ByteData+164+(q*16));
      memcpy(&TempY, ByteData+164+(q*16), 8);

      cout << TempX << " " << TempY << endl;
      Points.X.push_back(TempX);
      Points.Y.push_back(TempY);

    }

    Polylines.push_back(Points);

    cout << "---------" << endl;

    while (Skip < (FileLength*2)){

      ByteSwap(8, ByteData+Skip+4);
      memcpy(&RecordLength, ByteData+Skip+4, 4);
      ByteSwap(8, ByteData+Skip+8);
      memcpy(&shapetype, ByteData+Skip+8, 4);
      ByteSwap(8, ByteData+Skip+44);
      memcpy(&numparts, ByteData+Skip+44, 4);
      ByteSwap(8, ByteData+Skip+48);
      memcpy(&numpoints, ByteData+Skip+48, 4);  // number of XY pairs in the polyline

      PointData Points;
      for (int w = 0; w < numpoints; ++w){

        ByteSwap(8, ByteData+Skip+56+(w*16));
        memcpy(&TempX, ByteData+Skip+56+(w*16), 8);
        ByteSwap(8, ByteData+Skip+64+(w*16));
        memcpy(&TempY, ByteData+Skip+64+(w*16), 8);

        cout << TempX << " " << TempY << endl;
        Points.X.push_back(TempX);
        Points.Y.push_back(TempY);

      }

      Polylines.push_back(Points);
      Skip += (RecordLength*2) +8;

      cout << "------" << endl;

    }

  }

  // Close the file and return the point data
  fclose(file);
  return Polylines;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to convert an IndexChannelTree to a PointData object.
//
// Returns a vector of points.
// DTM 11/07/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
PointData LoadChannelTree(string Filename, int multistem_option, int trib_number)
{
  if((multistem_option != 0) && (multistem_option != 1) && (multistem_option !=2))
  {
    cout << "multistem_option must be 0 (mainstem only), 1 (all tributaries) or 2 (specify channel number).  Setting mainstem only default" << endl;
    multistem_option = 0;
  }

  ifstream channel_data_in;
  channel_data_in.open(Filename.c_str());

  if( channel_data_in.fail() )
  {
    cout << "\nFATAL ERROR: the channel network file \"" << Filename << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  PointData Points;

  int channel_number;
  int receiver_cnumber;
  int recevier_cnode;

  int node;
  int row;
  int col;

  float flow_dist;
  float elev;
  float drain_area;

  //int last_cn = 0;    // this is 1 if this is the first node in a channel
  int last_receiver_node = -1;
  //int last_receiver_channel;

  float XMinimum,YMinimum,DataResolution,NoDataValue;
  int NRows,NCols;

  channel_data_in >> NRows >> NCols >> XMinimum >> YMinimum >> DataResolution >> NoDataValue;
  float x,y;
  while( channel_data_in >> channel_number >> receiver_cnumber >> recevier_cnode
                         >> node >> row >> col >> flow_dist >> elev >> drain_area)
  {
    // get the receiver_channel and receiver node for the first channel (these will be recursive)
    if (last_receiver_node == -1)
    {
      last_receiver_node = recevier_cnode;
      //last_receiver_channel = receiver_cnumber;
    }
    // now load everything into the PointData object :-)

    if(multistem_option == 0)
    {
      if(channel_number == 0)
      {
        x = XMinimum + float(col)*DataResolution + 0.5*DataResolution;
        y = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
        Points.X.push_back(x);
        Points.Y.push_back(y);
      }
    }
    else if(multistem_option == 1)
    {
      x = XMinimum + float(col)*DataResolution + 0.5*DataResolution;
      y = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
      Points.X.push_back(x);
      Points.Y.push_back(y);
    }
    else if(multistem_option == 2)
    {
      if(channel_number == trib_number)
      {
        x = XMinimum + float(col)*DataResolution + 0.5*DataResolution;
        y = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
        Points.X.push_back(x);
        Points.Y.push_back(y);
      }
    }
    else
    {
      if(channel_number == 0)
      {
        x = XMinimum + float(col)*DataResolution + 0.5*DataResolution;
        y = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
        Points.X.push_back(x);
        Points.Y.push_back(y);
      }
    }
  }
  return Points;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get point data from vectors of X and Y coordinates
// FJC 17/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
PointData get_point_data_from_coordinates(vector<float>& X_coordinates, vector<float>& Y_coordinates)
{
  // convert to doubles
  vector<double> x(X_coordinates.begin(), X_coordinates.end());
  vector<double> y(Y_coordinates.begin(), Y_coordinates.end());


  PointData Points;
  Points.X = x;
  Points.Y = y;

  return Points;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get point data from vectors of X and Y coordinates
// FJC 17/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
PointData get_point_data_from_coordinates(vector<double>& X_coordinates, vector<double>& Y_coordinates)
{
  PointData Points;
  Points.X = X_coordinates;
  Points.Y = Y_coordinates;

  return Points;
}






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// These are functions for the OGC well know text object
// It is used to parse OGC WKT format strings in order to quantify a coordinate system
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This just populates the UTM coordinate system
void LSDOGCWKTCRSReader::create()
{
  populate_UTM_data();
}

// This reads an OGC WKT file and also populates the UTM projections
void LSDOGCWKTCRSReader::create(string filename)
{
  read_OGCWKT_file(filename);

  cout << "I've got the file, now for UTM" << endl;
  populate_UTM_data();
}

// Does the same thing as above but you can enter the filename and the path separately
void LSDOGCWKTCRSReader::create(string path, string filename)
{
    string fname = FixPath(path)+ filename;
    create(fname);
}

// This reads a file
void LSDOGCWKTCRSReader::read_OGCWKT_file(string filename)
{

  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load OGC WKT CRS data but the file " << endl
         << filename << endl
         << "doesn't exist!" << endl;
    exit(EXIT_FAILURE);
  }

  // initiate the string to hold the lines file
  string line_from_file;
  vector<string> this_string_vec;
  vector<string> lines;
  string temp_string;

  // get all the lines from the file
  while( ifs.good() )
  {
    getline(ifs, line_from_file);
    lines.push_back(line_from_file);
  }
  ifs.close();

  string PROJCS = "PROJCS";
  int n_lines = int(lines.size());
  for(int i = 0; i< n_lines; i++)
  {
    cout << lines[i] << endl << "========" << endl;
    // make sure the line is actually a projection
    if(lines[i].find(PROJCS) != string::npos)
    {
      parse_OGCWKT_line(lines[i]);
    }
  }

  //cout << endl << endl << endl << "========================" << endl;
  //cout << "Let me parse the first line." << endl;
  //parse_OGCWKT_line(lines[0]);

  //cout << endl << endl << endl << "========================" << endl;
  //cout << "Let me parse the third line." << endl;
  //parse_OGCWKT_line(lines[2]);

  //cout << endl << endl << endl << "========================" << endl;
  //cout << "Let me parse the fourth line." << endl;
  //parse_OGCWKT_line(lines[3]);

}


// This populates all the WGS 84 UTM coordinate systems
// This is rather wasteful in terms of memory and probably better to do it on the fly, 
// but I'm worried about mixing file reading and creating projections on the fly
void LSDOGCWKTCRSReader::populate_UTM_data()
{
  // There are 60 UTM zones, and each has a north and south, so in total there are 120 projections. 
  // The northern zones have EPSG codes 32601-32660 whereas southern have 32701-32760 

  int EPSG;
  string EPSG_type = "Transverse_Mercator";
  LSDReferenceEllipsoid RefEllipse("WGS84");

  map<string,double> UTM;
  UTM["k_0"] = 0.9996;

    
  // First do the northern zones   
  for (int i = 1; i<=60; i++)
  {
    EPSG = 32600+i;
    map<string,int> this_UTM_int;
    this_UTM_int["zone"] = i;
    this_UTM_int["is_north"] = 1; 

    LSDProjectionInfo ThisProjection(EPSG,EPSG_type,this_UTM_int,UTM,RefEllipse);

    projections_map[EPSG] = ThisProjection;
  }

  // Now the sourthern zones
  // First do the northern zones   
  for (int i = 1; i<=60; i++)
  {
    EPSG = 32700+i;
    map<string,int> this_UTM_int;
    this_UTM_int["zone"] = i;
    this_UTM_int["is_north"] = 1; 

    LSDProjectionInfo ThisProjection(EPSG,EPSG_type,this_UTM_int,UTM,RefEllipse);

    projections_map[EPSG] = ThisProjection;
  }  

}


// This parses a single line of OGC WKT text. The entire projection definition
// needs to be on the same line. 
// For formate see here: https://www.opengeospatial.org/standards/wkt-crs
void LSDOGCWKTCRSReader::parse_OGCWKT_line(string line)
{

  // First we get rid of the spaces. They are bad
  line = RemoveSpaces(line);

  // Now parse the projection
  string formatted_projection = parse_projection(line);

  // Parse the ellipsoid
  string formatted_ellipsoid = parse_ellipsoid(line);

  // Now parse the EPSG
  int EPSG = parse_EPSG(line);
  
  // These are parameters that will be fed into an LSDProjectionInfo Object
  string EPSG_type;
  LSDReferenceEllipsoid RefEllipse(formatted_ellipsoid);

  map<string,double> DoubleParams;
  map<string,int> IntParams;

  if (formatted_projection == "Lambert_Conformal_Conic_2SP" || 
      formatted_projection == "Albers_Conic_Equal_Area" ||
      formatted_projection == "Transverse_Mercator")
  {
    EPSG_type = formatted_projection;
    // Now we parse the parameters
    if (formatted_projection == "Lambert_Conformal_Conic_2SP" || 
        formatted_projection == "Albers_Conic_Equal_Area")
    {
      cout << "This is a conic projection. Let me get the parameters." << endl;
      DoubleParams = parse_conic_parameters(line);
      

    }
    else if (formatted_projection == "Transverse_Mercator")
    {
      cout << "This is a mercator projection." << endl;
    }

    //cout << "Let me make the projection" << endl;
    LSDProjectionInfo ThisProjection(EPSG,EPSG_type,IntParams,DoubleParams,RefEllipse);

    //cout << "Now I will add it to the map." << endl;
    projections_map[EPSG] = ThisProjection;
  }
  else
  {
    cout << "The projection type " << formatted_projection << " is not currently supported." << endl;
    cout << "Options are (this is case sensitive): "  << endl;
    cout << "Lambert_Conformal_Conic_2SP" << endl;
    cout << "Albers_Conic_Equal_Area" << endl;
    cout << "Transverse_Mercator" << endl;
  }

  // At the end of this logic you will now have an entry into the projection map with the EPSG code
  // as the key


}



// This function extract the ellipsoid from an OGC WKT line. 
// It looks for specific patterns for common ellipsoids. 
string LSDOGCWKTCRSReader::parse_ellipsoid(string line)
{
  
  string newline = line = RemoveSpaces(line);
  
  // Now get the kind of projection
  string keyword = "ELLIPSOID";
  if (line.find(keyword) == string::npos)
  {
    // ELLIPSOID WAS NOT FOUND. Change to SPHEROID
    keyword = "SPHEROID";

    if (line.find(keyword) == string::npos)
    {
      keyword = "NULL";
    }
  }
  string formatted_projection; 
  // Now parse the actual line
  if (keyword != "NULL")
  {
    string proj_string = parse_keyword(newline, keyword);

    unsigned first = proj_string.find_first_of("[")+2;
    unsigned last = proj_string.find_first_of(",")-1;

    string projection = proj_string.substr(first, last-first);
    //cout << "The ellipoid is: " << projection << endl;
    

    // Now parse. We constrain common examples
    if (projection == "WGS84" || projection == "WGS_84" || projection == "wgs_84" || 
        projection == "wgs84" || projection == "WGS_1984" || projection == "WGS1984" || 
        projection == "wgs1984" || projection == "wgs_1984")
    {
      formatted_projection = "WGS84";
    }
    else if (projection == "WGS72" || projection == "WGS_72" || projection == "wgs_72" || 
        projection == "wgs72" || projection == "WGS_1972" || projection == "WGS1972" || 
        projection == "wgs1972" || projection == "wgs_1972")
    {
      formatted_projection = "WGS72";
    }
    else if (projection == "GRS80" || projection == "GRS_80" || projection == "grs_80" || 
        projection == "grs80" || projection == "GRS_1980" || projection == "GRS1980" || 
        projection == "grs1980" || projection == "grs_1980")
    {
      formatted_projection = "GRS80";
    }
    else if (projection == "Clarke_1866" || projection == "Clarke1866" || projection == "CLARKE1866" || 
             projection == "clarke1866" || projection == "CLARKE_1866" || projection == "clarke_1866" || 
             projection == "clarke 1866" || projection == "CLARKE 1866" || projection == "Clarke 1866" || 
             projection == "Clarke_66" || projection == "Clarke66" || projection == "CLARKE66" || 
             projection == "clarke66" || projection == "CLARKE_66" || projection == "clarke_66" )
    {
      formatted_projection = "Clarke1866";
    }
    else
    {
      cout << "Warning, the ellipsoid is not one of WGS84, WGS72, Clarke1866, or GRS80." << endl;
      formatted_projection = projection; 
    }
  }
  else
  {
    cout << "Warning, I did not find an ellipsoid or spheroid code so I am assuming WGS84" << endl;
    formatted_projection = "WGS84";  
  }
  return formatted_projection;
}


// This parses the projection from an OGC WKT line
string LSDOGCWKTCRSReader::parse_projection(string line)
{
  
  string newline = line = RemoveSpaces(line);
  
  // Now get the kind of projection
  string keyword = "PROJECTION";
  string proj_string = parse_keyword(newline, keyword);
  

  unsigned first = proj_string.find_first_of("[")+2;
  unsigned last = proj_string.find_first_of("]")-1;

  string projection = proj_string.substr(first, last-first);

  
  string formatted_projection; 

  // Now parse. We constrain common examples
  if (projection == "Lambert_Conformal_Conic_2SP" || projection == "Lambert_Conformal_Conic" ||
     projection == "lambert_conformal_conic_2sp" || projection == "lambert_conformal_conic" ||
     projection == "LAMBERT_CONFORMAL_CONIC_2SP" || projection == "LAMBERT_CONFORMAL_CONIC" ||
     projection == "LCC" || projection == "lcc")
  {
    formatted_projection = "Lambert_Conformal_Conic_2SP";
  }
  else if (projection == "Transverse_Mercator" || projection == "transverse_mercator" ||
           projection == "TRANSVERSE_MERCATOR" || projection == "UTM" ||
           projection == "utm")
  {
    formatted_projection = "Transverse_Mercator";
  }
  else if (projection == "Albers_Conic_Equal_Area" || projection == "albers_conic_equal_area" ||
           projection == "ALBERS_CONIC_EQUAL_AREA")
  {
    formatted_projection = "Albers_Conic_Equal_Area";
  }
  else
  {
    formatted_projection = projection; 
  }
  //cout << "The projection is: " << formatted_projection << endl;
  return formatted_projection;
}

// This parses the EPSG code from and OGC WKT line
int LSDOGCWKTCRSReader::parse_EPSG(string line)
{
  // This parsing relies on the fact that the EPSG for the projection is the last EPSG code to appear
  unsigned last_EPSG = line.rfind("EPSG");
  //cout << "The position is: " << last_EPSG << endl;
  //cout << "The string length is: " << line.length() << endl;

  // The +8 is to cut out the actual "ESPG",". So you should have the EPSG code without the leading "
  // And then it finishes with a "
  string EPSG_part = line.substr(last_EPSG+7); 
  //cout << "EPSG part is: " << EPSG_part << endl;

  // Now get the actual EPSG string
  last_EPSG = EPSG_part.find_first_of("]")-1;
  string EPSG_str =  EPSG_part.substr(0,last_EPSG);
  //cout << "EPSG_string is: " << EPSG_str << endl;

  int EPSG_code = atoi(EPSG_str.c_str());
  return EPSG_code;

}

// This parses the parameters for a conic projection
// It works for lambert conformal conic and albers equal area conic
map<string,double> LSDOGCWKTCRSReader::parse_conic_parameters(string line)
{
  // we will be reading in these variables
  map<string,double> param_map;
  unsigned first, last;
  string keyword, params_str, params_str2;
  bool found;

  // Get the first parallel (known as phi_1 in Snyder 1987)
  found = false;
  keyword = "Standard_Parallel_1";
  if (line.find(keyword) == string::npos)
  {
    keyword = "standard_parallel_1";
    if (line.find(keyword) == string::npos)
    {
      cout << "Standard parallel 1 not found! Defaulting to 0." << endl;
      param_map["phi_1"] = 0;
    }
    else 
    {
      found = true;
    }
  } 
  else
  {    
    found = true;
  }  
    
  if (found)
  {  
    first = line.rfind(keyword)+1;
    params_str =  line.substr(first);
    first = params_str.find_first_of(",")+1;
    last = params_str.find_first_of("]");
    string params_str2 = params_str.substr(first,last-first);

    param_map["phi_1"] = atof(params_str2.c_str());
  }

  // Now get the second standard parallel (known as phi_1 in Snyder 1987)
  keyword = "Standard_Parallel_2";
  if (line.find(keyword) == string::npos)
  {
    keyword = "standard_parallel_2";
    if (line.find(keyword) == string::npos)
    {
      cout << "Standard parallel 2 not found! Defaulting to 0." << endl;
      param_map["phi_2"] = 0;
    }
    else 
    {
      found = true;
    }
  } 
  else
  {    
    found = true;
  }
  if (found)  
  {
    first = line.rfind(keyword)+1;
    params_str =  line.substr(first);
    first = params_str.find_first_of(",")+1;
    last = params_str.find_first_of("]");
    params_str2 = params_str.substr(first,last-first);

    param_map["phi_2"] = atof(params_str2.c_str());
  }


  // Now get the Latitude_Of_Origin (also known as phi_0 in Snyder 1987)
  keyword = "Latitude_Of_Origin";
  found = false;
  if (line.find(keyword) == string::npos) 
  {
    keyword = "latitude_of_center";
    if (line.find(keyword) == string::npos) 
    {
      //cout << "Warning! Origin latitude not found!" << endl;
      param_map["phi_0"] = 0;
    }
    else
    {
      found = true;
    }
  }
  else
  {
    found = true;
  }

  if (found)
  {
    first = line.rfind(keyword)+1;
    params_str =  line.substr(first);
    first = params_str.find_first_of(",")+1;
    last = params_str.find_first_of("]");
    params_str2 = params_str.substr(first,last-first);

    //cout << "The parameter string is: " << params_str2 << endl;
    param_map["phi_0"] = atof(params_str2.c_str());
  }

  // Now get the Longitude_Of_Origin (also known as central meridan or lambda_0 in Snyder 1987)
  keyword = "Longitude_Of_Origin";
  found = false;
  if (line.find(keyword) == string::npos) 
  {
    keyword = "longitude_of_center";
    if (line.find(keyword) == string::npos) 
    {
      keyword = "Central_Meridian";
      if (line.find(keyword) == string::npos) 
      {
        cout << "Warning! Origin latitude not found!" << endl;
        param_map["lambda_0"] = 0;
      }
      else
      {
        found = true;
      }
    }
    else
    {
      found = true;
    }
  }
  else
  {
    found = true;
  }

  if (found)
  {
    first = line.rfind(keyword)+1;
    params_str =  line.substr(first);
    first = params_str.find_first_of(",")+1;
    last = params_str.find_first_of("]");
    params_str2 = params_str.substr(first,last-first);

    //cout << "The parameter string is: " << params_str2 << endl;
    param_map["lambda_0"] = atof(params_str2.c_str());
  }


  cout << "Central latitude / phi_0 is: " << param_map["phi_0"] << endl;
  cout << "Central longitude / lambda_0 is: " << param_map["lambda_0"] << endl;
  cout << "Standard parallel 1 / phi_1 is: " << param_map["phi_1"] << endl;
  cout << "Standard parallel 2 / phi_2 is: " << param_map["phi_2"] << endl;
  return param_map;
}

// This is for parsing a specific, user-defined keyword
string LSDOGCWKTCRSReader::parse_keyword(string line, string keyword)
{
  // this is meant to parse brackets in a WKT file
  // we can use this recursively
  int first = line.find(keyword);
  string test_extract;

  if (first == -1)
  {
    //cout << "Keyword not found" << endl;
  }
  else
  {
    //cout << "I found the keyword. Now I will try to close the brackets" << endl;
    string newline = line.substr(first);
    //cout << "The new line is: " << endl;
    //cout << newline << endl;

    //const char* char_str = newline.c_str();

    unsigned start = newline.find_first_of('[');
    unsigned index = start;
    unsigned end;

    char open_bracket = '[';
    char closed_bracket = ']';



    int counter = 1;
    while (counter > 0)
    {
      index++;

      // Check bounds
      if( index >= newline.length() )
      {
        //cout << "Coming to the end of the line. I did not find a closing square bracket!" << endl;
        counter = 0;
      }
      else
      {
        const char c = newline[index];
        if(c == open_bracket)
        {
          counter++;
          //cout << "Found" << endl;
        }
        else if( c == closed_bracket)
        {
          counter--;
        }
      }

      if (counter == 0)
      {
        //cout << "Found a closing bracket! The character is: " << newline[index] << endl;
        end = index;
      }

    }
    //cout << "Start is: " << start << " and end is: " << end << endl;
    test_extract = newline.substr(start, end-start+1);


  }
  //cout << "test extract is: " << test_extract << endl;
  return test_extract;

}








//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function, sets up some vectors for holding
// ellipses and datums
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCoordinateConverterLLandUTM::create()
{
  // this sets up the ellipsoid and datum vectors
  vector<LSDEllipsoid> Ellipse_data_temp;

  // declare the names of the ellipsoids
  char t00[] = "Airy1830";
  char t01[] = "AiryModified";
  char t02[] = "AustralianNational";
  char t03[] = "Bessel1841Namibia";
  char t04[] = "Bessel1841";
  char t05[] = "Clarke1866";
  char t06[] = "Clarke1880";
  char t07[] = "EverestIndia1830";
  char t08[] = "EverestSabahSarawak";
  char t09[] = "EverestIndia1956";
  char t10[] = "EverestMalaysia1969";
  char t11[] = "EverestMalay_Sing";
  char t12[] = "EverestPakistan";
  char t13[] = "Fischer1960Modified";
  char t14[] = "Helmert1906";
  char t15[] = "Hough1960";
  char t16[] = "Indonesian1974";
  char t17[] = "International1924";
  char t18[] = "Krassovsky1940";
  char t19[] = "GRS80";
  char t20[] = "SouthAmerican1969";
  char t21[] = "WGS72";
  char t22[] = "WGS84";

  // add the ellipsoids to the vector
  LSDEllipsoid E00( 0, t00,    6377563.396,  299.3249646);
  LSDEllipsoid E01( 1, t01,    6377340.189,  299.3249646);
  LSDEllipsoid E02( 2, t02,  6378160,  298.25);
  LSDEllipsoid E03( 3, t03,  6377483.865,  299.1528128);
  LSDEllipsoid E04( 4, t04,    6377397.155,  299.1528128);
  LSDEllipsoid E05( 5, t05,    6378206.4,  294.9786982);
  LSDEllipsoid E06( 6, t06,    6378249.145,  293.465);
  LSDEllipsoid E07( 7, t07,  6377276.345,  300.8017);
  LSDEllipsoid E08( 8, t08,  6377298.556,  300.8017);
  LSDEllipsoid E09( 9, t09,  6377301.243,  300.8017);
  LSDEllipsoid E10(10, t10,  6377295.664,  300.8017);  //Dana has no datum that uses this LSDEllipsoid E00!
  LSDEllipsoid E11(11, t11,  6377304.063,  300.8017);
  LSDEllipsoid E12(12, t12,  6377309.613,  300.8017);
  LSDEllipsoid E13(13, t13,  6378155,  298.3);
  LSDEllipsoid E14(14, t14,    6378200,  298.3);
  LSDEllipsoid E15(15, t15,    6378270,  297);
  LSDEllipsoid E16(16, t16,    6378160,  298.247);
  LSDEllipsoid E17(17, t17,  6378388,  297);
  LSDEllipsoid E18(18, t18,    6378245,  298.3);
  LSDEllipsoid E19(19, t19,      6378137,  298.257222101);
  LSDEllipsoid E20(20, t20,  6378160,  298.25);
  LSDEllipsoid E21(21, t21,      6378135,  298.26);
  LSDEllipsoid E22(22, t22,      6378137,  298.257223563);

  Ellipse_data_temp.push_back(E00);
  Ellipse_data_temp.push_back(E01);
  Ellipse_data_temp.push_back(E02);
  Ellipse_data_temp.push_back(E03);
  Ellipse_data_temp.push_back(E04);
  Ellipse_data_temp.push_back(E05);
  Ellipse_data_temp.push_back(E06);
  Ellipse_data_temp.push_back(E07);
  Ellipse_data_temp.push_back(E08);
  Ellipse_data_temp.push_back(E09);
  Ellipse_data_temp.push_back(E10);
  Ellipse_data_temp.push_back(E11);
  Ellipse_data_temp.push_back(E12);
  Ellipse_data_temp.push_back(E13);
  Ellipse_data_temp.push_back(E14);
  Ellipse_data_temp.push_back(E15);
  Ellipse_data_temp.push_back(E16);
  Ellipse_data_temp.push_back(E17);
  Ellipse_data_temp.push_back(E18);
  Ellipse_data_temp.push_back(E19);
  Ellipse_data_temp.push_back(E20);
  Ellipse_data_temp.push_back(E21);
  Ellipse_data_temp.push_back(E22);

  //names for ellipsoidId's
  int eClarke1866 = 5;
  int eGRS80 = 19;
  int eWGS72 = 21;
  int eWGS84 = 22;


  // now for the datum
  vector<LSDDatum> Datum_data_temp;

  // initiate the datum names
  char T00[] = "NAD27_AK";
  char T01[] = "NAD27_AK_AleutiansE";
  char T02[] = "NAD27_AK_AleutiansW";
  char T03[] = "NAD27_Bahamas";
  char T04[] = "NAD27_Bahamas_SanSalv";
  char T05[] = "NAD27_AB_BC";
  char T06[] = "NAD27_MB_ON";
  char T07[] = "NAD27_NB_NL_NS_QC";
  char T08[] = "NAD27_NT_SK";
  char T09[] = "NAD27_YT";
  char T10[] = "NAD27_CanalZone";
  char T11[] = "NAD27_Cuba";
  char T12[] = "NAD27_Greenland";
  char T13[] = "NAD27_Carribean";
  char T14[] = "NAD27_CtrlAmerica";
  char T15[] = "NAD27_Canada";
  char T16[] = "NAD27_ConUS";
  char T17[] = "NAD27_ConUS_East";
  char T18[] = "NAD27_ConUS_West";
  char T19[] = "NAD27_Mexico";
  char T20[] = "NAD83_AK";
  char T21[] = "NAD83_AK_Aleutians";
  char T22[] = "NAD83_Canada";
  char T23[] = "NAD83_ConUS";
  char T24[] = "NAD83_Hawaii";
  char T25[] = "NAD83_Mexico_CtrlAmerica";
  char T26[] = "WGS72";
  char T27[] = "WGS84";

  LSDDatum D00(0, T00,      eClarke1866,  -5,  135,  172); //NAD27 for Alaska Excluding Aleutians
  LSDDatum D01( 1, T01,  eClarke1866,  -2,  152,  149); //NAD27 for Aleutians East of 180W
  LSDDatum D02( 2, T02,  eClarke1866,  2,  204,  105); //NAD27 for Aleutians West of 180W
  LSDDatum D03( 3, T03,    eClarke1866,  -4,  154,  178); //NAD27 for Bahamas Except SanSalvadorIsland
  LSDDatum D04( 4, T04,  eClarke1866,  1,  140,  165); //NAD27 for Bahamas SanSalvadorIsland
  LSDDatum D05( 5, T05,    eClarke1866,  -7,  162,  188); //NAD27 for Canada Alberta BritishColumbia
  LSDDatum D06( 6, T06,    eClarke1866,  -9,  157,  184); //NAD27 for Canada Manitoba Ontario
  LSDDatum D07( 7, T07,    eClarke1866,  -22,  160,  190); //NAD27 for Canada NewBrunswick Newfoundland NovaScotia Quebec
  LSDDatum D08( 8, T08,    eClarke1866,  4,  159,  188); //NAD27 for Canada NorthwestTerritories Saskatchewan
  LSDDatum D09( 9, T09,      eClarke1866,  -7,  139,  181); //NAD27 for Canada Yukon
  LSDDatum D10(10, T10,    eClarke1866,  0,  125,  201); //NAD27 for CanalZone (ER: is that Panama??)
  LSDDatum D11(11, T11,      eClarke1866,  -9,  152,  178); //NAD27 for Cuba
  LSDDatum D12(12, T12,    eClarke1866,  11,  114,  195); //NAD27 for Greenland (HayesPeninsula)
  LSDDatum D13(13, T13,    eClarke1866,  -3,  142,  183); //NAD27 for Antigua Barbados Barbuda Caicos Cuba DominicanRep GrandCayman Jamaica Turks
  LSDDatum D14(14, T14,    eClarke1866,  0,  125,  194); //NAD27 for Belize CostaRica ElSalvador Guatemala Honduras Nicaragua
  LSDDatum D15(15, T15,    eClarke1866,  -10,  158,  187); //NAD27 for Canada
  LSDDatum D16(16, T16,    eClarke1866,  -8,  160,  176); //NAD27 for CONUS
  LSDDatum D17(17, T17,    eClarke1866,  -9,  161,  179); //NAD27 for CONUS East of Mississippi Including Louisiana Missouri Minnesota
  LSDDatum D18(18, T18,    eClarke1866,  -8,  159,  175); //NAD27 for CONUS West of Mississippi Excluding Louisiana Missouri Minnesota
  LSDDatum D19(19, T19,    eClarke1866,  -12,  130,  190); //NAD27 for Mexico
  LSDDatum D20(20, T20,      eGRS80,    0,  0,  0); //NAD83 for Alaska Excluding Aleutians
  LSDDatum D21(21, T21,    eGRS80,    -2,  0,  4); //NAD83 for Aleutians
  LSDDatum D22(22, T22,    eGRS80,    0,  0,  0); //NAD83 for Canada
  LSDDatum D23(23, T23,    eGRS80,    0,  0,  0); //NAD83 for CONUS
  LSDDatum D24(24, T24,    eGRS80,    1,  1,  -1); //NAD83 for Hawaii
  LSDDatum D25(25, T25,  eGRS80,    0,  0,  0); //NAD83 for Mexico CentralAmerica
  LSDDatum D26(26, T26,      eWGS72,    0,  0,  0); //WGS72 for world
  LSDDatum D27(27, T27,      eWGS84,    0,  0,  0); //WGS84 for world

  Datum_data_temp.push_back(D00);
  Datum_data_temp.push_back(D01);
  Datum_data_temp.push_back(D02);
  Datum_data_temp.push_back(D03);
  Datum_data_temp.push_back(D04);
  Datum_data_temp.push_back(D05);
  Datum_data_temp.push_back(D06);
  Datum_data_temp.push_back(D07);
  Datum_data_temp.push_back(D08);
  Datum_data_temp.push_back(D09);
  Datum_data_temp.push_back(D10);
  Datum_data_temp.push_back(D11);
  Datum_data_temp.push_back(D12);
  Datum_data_temp.push_back(D13);
  Datum_data_temp.push_back(D14);
  Datum_data_temp.push_back(D15);
  Datum_data_temp.push_back(D16);
  Datum_data_temp.push_back(D17);
  Datum_data_temp.push_back(D18);
  Datum_data_temp.push_back(D19);
  Datum_data_temp.push_back(D20);
  Datum_data_temp.push_back(D21);
  Datum_data_temp.push_back(D22);
  Datum_data_temp.push_back(D23);
  Datum_data_temp.push_back(D24);
  Datum_data_temp.push_back(D25);
  Datum_data_temp.push_back(D26);
  Datum_data_temp.push_back(D27);

  Ellipsoids = Ellipse_data_temp;
  Datums = Datum_data_temp;

  RADIANS_PER_DEGREE = M_PI/180.0;
  DEGREES_PER_RADIAN = 180.0/M_PI;

  /** Useful constants **/
  TWOPI = 2.0 * M_PI;
  HALFPI = M_PI / 2.0;

  // Grid granularity for rounding UTM coordinates to generate MapXY.
  grid_size = 100000.0;    // 100 km grid

  // WGS84 Parameters
  WGS84_A=6378137.0;    // major axis
  WGS84_B=6356752.31424518;  // minor axis
  WGS84_F=0.0033528107;    // ellipsoid flattening
  WGS84_E=0.0818191908;    // first eccentricity
  WGS84_EP=0.0820944379;    // second eccentricity

  // UTM Parameters
  UTM_K0=0.9996;      // scale factor
  UTM_FE=500000.0;    // false easting
  UTM_FN_N=0.0;           // false northing, northern hemisphere
  UTM_FN_S=10000000.0;    // false northing, southern hemisphere
  UTM_E2=(WGS84_E*WGS84_E);  // e^2
  UTM_E4=(UTM_E2*UTM_E2);    // e^4
  UTM_E6=(UTM_E4*UTM_E2);    // e^6
  UTM_EP2=(UTM_E2/(1-UTM_E2));  // e'^2

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// converts LatLong to UTM coords
// 3/22/95: by ChuckGantz chuck.gantz@globalstar.com, from USGS Bulletin 1532.
// Lat and Long are in degrees;
// North latitudes and East Longitudes are positive.
//
// Minor modifications for our objects by SMM, 07/12/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCoordinateConverterLLandUTM::LLtoUTM(int eId, double Lat, double Long,
             double& Northing, double& Easting, int& Zone)
{

  double UTMEasting, UTMNorthing;

  double a = WGS84_A;
  double eccSquared = UTM_E2;
  double k0 = UTM_K0;

  double LongOrigin;
  double eccPrimeSquared;
  double N, T, C, A, M;

        //Make sure the longitude is between -180.00 .. 179.9
  double LongTemp = (Long+180)-int((Long+180)/360)*360-180;

  double LatRad = Lat*RADIANS_PER_DEGREE;
  double LongRad = LongTemp*RADIANS_PER_DEGREE;
  double LongOriginRad;
  int    ZoneNumber;

  ZoneNumber = int((LongTemp + 180)/6) + 1;

  if( Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 )
  {
    ZoneNumber = 32;
  }
        // Special zones for Svalbard
  if( Lat >= 72.0 && Lat < 84.0 )
  {
    if(      LongTemp >= 0.0  && LongTemp <  9.0 ) ZoneNumber = 31;
    else if( LongTemp >= 9.0  && LongTemp < 21.0 ) ZoneNumber = 33;
    else if( LongTemp >= 21.0 && LongTemp < 33.0 ) ZoneNumber = 35;
    else if( LongTemp >= 33.0 && LongTemp < 42.0 ) ZoneNumber = 37;
   }
        // +3 puts origin in middle of zone
  LongOrigin = (ZoneNumber - 1)*6 - 180 + 3;
  LongOriginRad = LongOrigin * RADIANS_PER_DEGREE;

  //compute the UTM Zone from the latitude and longitude
  cout << "Zone is  " << ZoneNumber << endl;
  Zone = ZoneNumber;

  eccPrimeSquared = (eccSquared)/(1-eccSquared);

  N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad));
  T = tan(LatRad)*tan(LatRad);
  C = eccPrimeSquared*cos(LatRad)*cos(LatRad);
  A = cos(LatRad)*(LongRad-LongOriginRad);

  M = a*((1 - eccSquared/4 - 3*eccSquared*eccSquared/64
                - 5*eccSquared*eccSquared*eccSquared/256) * LatRad
               - (3*eccSquared/8 + 3*eccSquared*eccSquared/32
                  + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
               + (15*eccSquared*eccSquared/256
                  + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad)
               - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad));

  UTMEasting = (double)
          (k0*N*(A+(1-T+C)*A*A*A/6
                 + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
           + 500000.0);

  UTMNorthing = (double)
          (k0*(M+N*tan(LatRad)
               *(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                 + (61-58*T+T*T+600*C-330*eccPrimeSquared)*A*A*A*A*A*A/720)));

  if(Lat < 0)
  {
    //10000000 meter offset for southern hemisphere
    UTMNorthing += 10000000.0;
  }

  Northing = UTMNorthing;
  Easting = UTMEasting;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// converts LatLong to UTM coords
// 3/22/95: by ChuckGantz chuck.gantz@globalstar.com, from USGS Bulletin 1532.
// Lat and Long are in degrees;
// North latitudes and East Longitudes are positive.
//
// Minor modifications for our objects by SMM, 07/12/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCoordinateConverterLLandUTM::LLtoUTM_ForceZone(int eId, double Lat, double Long,
             double& Northing, double& Easting, int Zone)
{

  double UTMEasting, UTMNorthing;

  double a = WGS84_A;
  double eccSquared = UTM_E2;
  double k0 = UTM_K0;

  double LongOrigin;
  double eccPrimeSquared;
  double N, T, C, A, M;

        //Make sure the longitude is between -180.00 .. 179.9
  double LongTemp = (Long+180)-int((Long+180)/360)*360-180;

  double LatRad = Lat*RADIANS_PER_DEGREE;
  double LongRad = LongTemp*RADIANS_PER_DEGREE;
  double LongOriginRad;
  int ZoneNumber;

  ZoneNumber = int((LongTemp + 180)/6) + 1;

  if( Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 )
  {
    ZoneNumber = 32;
  }
        // Special zones for Svalbard
  if( Lat >= 72.0 && Lat < 84.0 )
  {
    if(      LongTemp >= 0.0  && LongTemp <  9.0 ) ZoneNumber = 31;
    else if( LongTemp >= 9.0  && LongTemp < 21.0 ) ZoneNumber = 33;
    else if( LongTemp >= 21.0 && LongTemp < 33.0 ) ZoneNumber = 35;
    else if( LongTemp >= 33.0 && LongTemp < 42.0 ) ZoneNumber = 37;
  }

  if (ZoneNumber != Zone)
  {
    cout << "WARNING: the point is located in zone " << ZoneNumber << endl
         << "but you are forcing the points into Zone " << Zone << endl;
   }

  // +3 puts origin in middle of zone
  LongOrigin = (Zone - 1)*6 - 180 + 3;
  LongOriginRad = LongOrigin * RADIANS_PER_DEGREE;

  eccPrimeSquared = (eccSquared)/(1-eccSquared);

  N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad));
  T = tan(LatRad)*tan(LatRad);
  C = eccPrimeSquared*cos(LatRad)*cos(LatRad);
  A = cos(LatRad)*(LongRad-LongOriginRad);

  M = a*((1 - eccSquared/4 - 3*eccSquared*eccSquared/64
                - 5*eccSquared*eccSquared*eccSquared/256) * LatRad
               - (3*eccSquared/8 + 3*eccSquared*eccSquared/32
                  + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
               + (15*eccSquared*eccSquared/256
                  + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad)
               - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad));

  UTMEasting = (double)
          (k0*N*(A+(1-T+C)*A*A*A/6
                 + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
           + 500000.0);

  UTMNorthing = (double)
          (k0*(M+N*tan(LatRad)
               *(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                 + (61-58*T+T*T+600*C-330*eccPrimeSquared)*A*A*A*A*A*A/720)));

  if(Lat < 0)
  {
    //10000000 meter offset for southern hemisphere
    UTMNorthing += 10000000.0;
  }

  Northing = UTMNorthing;
  Easting = UTMEasting;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// converts UTM coords to LatLong;  3/22/95: by ChuckGantz chuck.gantz@globalstar.com, from USGS Bulletin 1532.
// Lat and Long are in degrees;  North latitudes and East Longitudes are positive.
//
// Minor modifications for our objects by SMM, 07/12/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCoordinateConverterLLandUTM::UTMtoLL(int eId, double UTMNorthing, double UTMEasting,
                                             int UTMZone, bool isNorth,
                                             double& Lat, double& Long)
{

  double k0 = UTM_K0;
  double a = WGS84_A;
  double eccSquared = UTM_E2;
  double eccPrimeSquared;
  double e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared));
  double N1, T1, C1, R1, D, M;
  double LongOrigin;
  double mu, phi1Rad;
  double x, y;
  int ZoneNumber;

  x = UTMEasting - 500000.0; //remove 500,000 meter offset for longitude
  y = UTMNorthing;

  ZoneNumber = UTMZone;
  if(isNorth == false)
  {
    //cout << "Line 1010, you are in the Southern hemisphere!"<< endl;
    //remove 10,000,000 meter offset used for southern hemisphere
    y -= 10000000.0;
  }

        //+3 puts origin in middle of zone
  LongOrigin = (ZoneNumber - 1)*6 - 180 + 3;
  eccPrimeSquared = (eccSquared)/(1-eccSquared);

  M = y / k0;
  mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64
                   -5*eccSquared*eccSquared*eccSquared/256));

  phi1Rad = mu + ((3*e1/2-27*e1*e1*e1/32)*sin(2*mu)
                        + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
                        + (151*e1*e1*e1/96)*sin(6*mu));

  N1 = a/sqrt(1-eccSquared*sin(phi1Rad)*sin(phi1Rad));
  T1 = tan(phi1Rad)*tan(phi1Rad);
  C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad);
  R1 = a*(1-eccSquared)/pow(1-eccSquared*sin(phi1Rad)*sin(phi1Rad), 1.5);
  D = x/(N1*k0);

  Lat = phi1Rad - ((N1*tan(phi1Rad)/R1)
                         *(D*D/2
                           -(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24
                           +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared
                             -3*C1*C1)*D*D*D*D*D*D/720));

  Lat = Lat * DEGREES_PER_RADIAN;

  Long = ((D-(1+2*T1+C1)*D*D*D/6
                 +(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)
                 *D*D*D*D*D/120)
                / cos(phi1Rad));
  Long = LongOrigin + Long * DEGREES_PER_RADIAN;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//==============================================================================
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Converts from british national grid to lat long
// From hannah fry's website:
// http://www.hannahfry.co.uk/blog//2011/10/10/converting-british-national-grid-to-latitude-and-longitude
//
// see also
// https://github.com/chrisveness/geodesy/blob/master/osgridref.js
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCoordinateConverterLLandUTM::BNGtoLL(double Northing, double Easting, double& Lat, double& Long)
{
  // The Airy 180 semi-major and semi-minor axes used for OSGB36 (m)
  double a = 6377563.396;
  double b =  6356256.909;

  double F0 = 0.9996012717;    // scale factor on the central meridian
  double lat0 = 49*M_PI/180;   // Latitude of true origin (radians)
  double lon0 = -2*M_PI/180;      // Longtitude of true origin and central meridian (radians)

  // Northing & easting of true origin (m)
  double N0 = -100000;
  double E0 = 400000;

  double e2 = 1 - (b*b)/(a*a);  // eccentricity squared
  double n = (a-b)/(a+b);

  // Initialise the iterative variables
  double lat = lat0;
  double M = 0;
  double M1,M2,M3,M4;

  // Loop until you get the right latitude
  while(  (Northing - N0 - M) > 0.00001)
  {
    lat = (Northing-N0-M)/(a*F0) + lat;
    M1 = (1 + n + (5/4)*n*n + (5/4)*n*n*n) * (lat-lat0);
    M2 = (3*n + 3*n*n + (21/8)*n*n*n) * sin(lat-lat0) * cos(lat+lat0);
    M3 = ((15/8)*n*n + (15/8)*n*n*n) * sin(2*(lat-lat0)) * cos(2*(lat+lat0));
    M4 = (35/24)*n*n*n * sin(3*(lat-lat0)) * cos(3*(lat+lat0));

    // calculate the meridional arc
    M = b * F0 * (M1 - M2 + M3 - M4);
  }

  // transverse radius of curvature
  double nu = a*F0/sqrt(1-e2*sin(lat)*sin(lat));

  // meridional radius of curvature
  double rho = a*F0*(1-e2)*pow((1-e2*sin(lat)*sin(lat)),-1.5);
  double eta2 = nu/rho-1;
  double nu3 = nu*nu*nu;
  double nu5 = nu*nu*nu*nu*nu;
  double nu7 = nu*nu*nu*nu*nu*nu*nu;    // BATMAN!!
  double tanlat2 = tan(lat)*tan(lat);
  double tanlat4 = tanlat2*tanlat2;

  double secLat = 1./cos(lat);
  double VII = tan(lat)/(2*rho*nu);
  double VIII = tan(lat)/(24*rho*nu3)*(5+3*tanlat2+eta2-9*tanlat2*eta2);
  double IX = tan(lat)/(720*rho*nu5)*(61+90*tanlat2+45*tanlat4);
  double X = secLat/nu;
  double XI = secLat/(6*nu3)*(nu/rho+2*tanlat2);
  double XII = secLat/(120*nu5)*(5+28*tanlat2+24*tanlat4);
  double XIIA = secLat/(5040*nu7)*(61+662*tanlat2+1320*tanlat4+720*tanlat2*tanlat4);
  double dE = Easting-E0;

  lat = lat - VII*dE*dE + VIII*dE*dE*dE*dE - IX*dE*dE*dE*dE*dE*dE;
  double lon = lon0 + X*dE - XI*dE*dE*dE + XII*dE*dE*dE*dE*dE - XIIA*dE*dE*dE*dE*dE*dE*dE;

  // convert to degrees
  lat = lat*180/M_PI;
  lon = lon*180/M_PI;

  Lat = lat;
  Long = lon;


}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// converts LatLongHt in datum dIn, to LatLongHt in datum dTo;
// 2002dec: by Eugene Reimer, from PeterDana equations.
// Lat and Long params are in degrees;
// North latitudes and East longitudes are positive;  Height is in meters;
// ==This approach to Datum-conversion is a waste of time;
// to get acceptable accuracy a large table is needed -- see NADCON, NTv2...
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCoordinateConverterLLandUTM::DatumConvert(int dIn, double LatIn,
                  double LongIn, double HtIn,
                  int dTo,  double& LatTo, double& LongTo, double& HtTo)
{
  double a,ee,N,X,Y,Z,EE,p,b,t;

  //--transform to XYZ, using the "In" ellipsoid
  //LongIn += 180;
  a = Ellipsoids[Datums[dIn].eId].EquatorialRadius;
  ee= Ellipsoids[Datums[dIn].eId].eccSquared;
  N = a / sqrt(1 - ee*sin(rad(LatIn))*sin(rad(LatIn)));
  X = (N + HtIn) * cos(rad(LatIn)) * cos(rad(LongIn));
  Y = (N + HtIn) * cos(rad(LatIn)) * sin(rad(LongIn));
  Z = (N*(1-ee) + HtIn) * sin(rad(LatIn));

  //--apply delta-terms dX dY dZ
  //cout<<"\tX:" <<X <<" Y:" <<Y <<" Z:" <<Z;    //==DEBUG
  X+= Datums[dIn].dX - Datums[dTo].dX;
  Y+= Datums[dIn].dY - Datums[dTo].dY;
  Z+= Datums[dIn].dZ - Datums[dTo].dZ;
  //cout<<"\tX:" <<X <<" Y:" <<Y <<" Z:" <<Z;    //==DEBUG

  //--transform back to LatLongHeight, using the "To" ellipsoid
  a = Ellipsoids[Datums[dTo].eId].EquatorialRadius;
  ee= Ellipsoids[Datums[dTo].eId].eccSquared;
  EE= ee/(1-ee);
  p = sqrt(X*X + Y*Y);
  b = a*sqrt(1-ee);
  t = atan(Z*a/(p*b));
  LatTo = atan((Z+EE*b*sin(t)*sin(t)*sin(t)) / (p-ee*a*cos(t)*cos(t)*cos(t)));
  LongTo= atan(Y/X);
  HtTo  = p/cos(LatTo) - a/sqrt(1-ee*sin(LatTo)*sin(LatTo));
  LatTo = deg(LatTo);
  LongTo = deg(LongTo);
  LongTo -= 180;
}




////==============================================================================
////=========================================================
////=========================================================
////=========================================================
////=========================================================
//// YYY    YYY   YYY    YYY     OOOO
////  YYY  YYY     YYY  YYY     OOOOOO
////   YYYYYY       YYYYYY    OOO    OOO
////    YYYY         YYYY     OOO    OOO
////    YYYY         YYYY     OOO    OOO
////    YYYY         YYYY     OOO    OOO
////    YYYY         YYYY       OOOOOO
////    YYYY         YYYY        OOOO
//// BELOW HERE IS ALL THE NEW GEOREFERNCING STUFF
//// Developed by SMM in 2019
//==============================================================================
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// New Reference Ellipsoid object
//==================================================================
void LSDReferenceEllipsoid::create()
{
  //cout << "Setting your ellipsoid to WGS84"  << endl;
  string Name = "WGS84";
  set_ellipsoid_parameters(Name);
}


void LSDReferenceEllipsoid::create(string Name)
{
  cout << "Setting your ellipsoid"  << endl;
  set_ellipsoid_parameters(Name);
}


// This sets the ellipsoid parameters from a list of common ellipsoids
// TODO (SMM 19/01/2019) Add a martian coordinate system
void LSDReferenceEllipsoid::set_ellipsoid_parameters(string Name)
{
  // This stored ellipsoid parameters
  map<string, pair<double,double> > ellipsoid_map;

  // declare the names of the ellipsoids
  // The nmubers can be found here: https://en.wikipedia.org/wiki/Earth_ellipsoid
  ellipsoid_map["Airy1830"] = make_pair(6377563.396,  299.3249646);
  ellipsoid_map["AiryModified"] = make_pair(6377340.189,  299.3249646);
  ellipsoid_map["AustralianNational"] = make_pair(6378160,  298.25);
  ellipsoid_map["Bessel1841Namibia"] = make_pair(6377483.865,  299.1528128);
  ellipsoid_map["Bessel1841"] = make_pair(6377397.155,  299.1528128);
  ellipsoid_map["Clarke1866"] = make_pair(6378206.4,  294.9786982);
  ellipsoid_map["Clarke1880"] = make_pair(6378249.145,  293.465);
  ellipsoid_map["EverestIndia1830"] = make_pair(6377276.345,  300.8017);
  ellipsoid_map["EverestSabahSarawak"] = make_pair(6377298.556,  300.8017);
  ellipsoid_map["EverestIndia1956"] = make_pair(6377298.556,  300.8017);
  ellipsoid_map["EverestMalaysia1969"] = make_pair(6377301.243,  300.8017);
  ellipsoid_map["EverestMalay_Sing"] = make_pair(6377295.664,  300.8017);
  ellipsoid_map["EverestPakistan"] = make_pair(6377304.063,  300.8017);
  ellipsoid_map["Fischer1960Modified"] = make_pair(6378155,  298.3);
  ellipsoid_map["Helmert1906"] = make_pair(6378200,  298.3);
  ellipsoid_map["Hough1960"] = make_pair(6378270,  297);
  ellipsoid_map["Indonesian1974"] = make_pair(6378160,  298.247);
  ellipsoid_map["International1924"] = make_pair(6378388,  297);
  ellipsoid_map["Krassovsky1940"] = make_pair(6378245,  298.3);
  ellipsoid_map["GRS80"] = make_pair(6378137,  298.257222101);
  ellipsoid_map["SouthAmerican1969"] = make_pair(6378160,  298.25);
  ellipsoid_map["WGS72"] = make_pair(6378135,  298.26);
  ellipsoid_map["WGS84"] = make_pair(6378137,  298.257223563);
  ellipsoid_map["GRS_1967_Truncated"] = make_pair(6378160,  298.25);
  ellipsoid_map["Mars_2000_IAU_IAG"] = make_pair(3396190.0,  169.89444722361179);

  if ( ellipsoid_map.find(Name) == ellipsoid_map.end() )
  {
    cout << "I did not find the ellipsoid " << Name << ". I am defaulting to WGS84." << endl;
    Name = "WGS84";
  }
  else
  {
    cout.precision(12);
    //cout << "I did find the ellipsoid. The parameters are: " << ellipsoid_map[Name].first << " " << ellipsoid_map[Name].second << endl;
  }

  EquatorialRadius = ellipsoid_map[Name].first;
  
  // now calculate all the subsidiary data elements. 
  f_inverse = ellipsoid_map[Name].second;
  f = 1/ellipsoid_map[Name].second;
  
  // See Snyder 1987 for details: https://pubs.er.usgs.gov/publication/pp1395
  // Page 13 of that document
  eccSquared = 2*f - f*f;
  ecc = sqrt(eccSquared);

}

// Placeholder function to make sure the ellipsoid is valid. 
bool LSDReferenceEllipsoid::check_ellipsoid()
{
  bool It_fucking_works_since_we_checked_during_creation = true;
  return It_fucking_works_since_we_checked_during_creation;

}



//==================================================================
//==================================================================
// The projection info object
//==================================================================

// Default constructor. Assumes a default projections
void LSDProjectionInfo::create()
{
  //cout << "You a creating a projection info object with no information." << endl;
  //cout << "I will assume a default projection" << endl;
  //cout << "This is a lambert conformal conical with parameters from Snyder 1987 page 296" << endl;
  //int EPSG_code = 0;
  //set_projection(EPSG_code);

}

// Constructor where you give the EPSG code
void LSDProjectionInfo::create(int EPSG_code)
{
  //cout << "You a creating a projection with EPSG code " << EPSG_code << endl;
  set_projection(EPSG_code);

}

/// This constructor just assigns all the data members from values supplied by the user. 
void LSDProjectionInfo::create(int tEPSG_code, string tEPSG_type, 
                      map<string,int> tIntProjParameters, map<string,double> tProjParameters, 
                      LSDReferenceEllipsoid tRefEllipsoid)
{
  EPSG = tEPSG_code;
  EPSG_type = tEPSG_type;
  IntProjParameters = tIntProjParameters;
  ProjParameters = tProjParameters;
  ReferenceEllipsoid = tRefEllipsoid;


}

/// A placeholder function that will eventually check for validity of the projection
bool LSDProjectionInfo::check_projection()
{
  // Currently no checking is implemented. To be done once we implement the WKT file reading functionality
  bool this_needs_to_be_updated = true;
  return this_needs_to_be_updated;
}

// This sets the projection. For the EPSG code to work its detauls must be 
// entered into the file below. 
// In the future we might consider ingesting xml EPSG details 
// in order to streamline this. 
void LSDProjectionInfo::set_projection(int EPSG_code)
{

  // We need a way to do this from WKT. 
  // See http://spatialreference.org/ref/esri/north-america-lambert-conformal-conic/


  // This is a series of projection parameters
  map<int, map<string, double> > EPSG_list;

  // This is a series of projection parameters as integers
  map<int, map<string, int> > EPSG_int_list;

  // This can (for the time being) have the different types of projection for 
  // which we intend to write code to perform forward and inverse transformations. 
  // The types will be:
  // LCC (Lamber Conformal Conical)
  // AEA (Albers Equal Area)
  // UTM (Universal Transverse Mercator)
  // BNG (British National Grid)
  map<int, string> EPSG_type_map;

  

  // Conical projections include the origin latitude and longitude (phi_0, lambda_0)
  // and the standard parallels (phi_1 and phi_2)
  // These parameters apply to both lambert conformal conical and albers equal area projections
  // See Snyder 1987 page 101 and page 106
  // THESE ARE REPORTED IN DEGREES!!!!!! 
  map<string,double> Conical;
  Conical["phi_0"] = 23;
  Conical["lambda_0"] = -96;        // West is negative
  Conical["phi_1"] = 33;
  Conical["phi_2"] = 45; 
  
  // A placeholder!! Not a real EPSG.
  // This comes from the snyder worked example on page 296
  EPSG_list[0] = Conical;   
  EPSG_type_map[0] = "Lambert_Conformal_Conic_2SP";     // this is not random, it comes from standard WKT of projections. 
                                                    // see, e.g.,  http://spatialreference.org/ref/esri/north-america-lambert-conformal-conic/ogcwkt/       

  // Now for the example on page 101
  Conical["phi_0"] = 23;
  Conical["lambda_0"] = -96;        // West is negative
  Conical["phi_1"] = 29.5;
  Conical["phi_2"] = 45.5;  

  EPSG_list[1] = Conical;   
  EPSG_type_map[1] = "Albers_Conic_Equal_Area";     // this is not random, it comes from standard WKT of projections. 
                                                    // see, e.g., http://spatialreference.org/ref/esri/102001/ogcwkt/)     

  
  map<string,double> UTM;
  UTM["k_0"] = 0.9996;

  map<string,int> UTM_int;
  UTM_int["zone"] =34;
  UTM_int["is_north"] = 1;

  EPSG_list[32634] = UTM;
  EPSG_int_list[32634] = UTM_int;
  EPSG_type_map[32634] = "Transverse_Mercator";

  // now zone 18 (for the example)
  UTM_int["zone"] =18;
  UTM_int["is_north"] = 1;
  EPSG_list[32618] = UTM;
  EPSG_int_list[32618] = UTM_int;
  EPSG_type_map[32618] = "Transverse_Mercator";

  
  if (EPSG_list.find(EPSG_code) == EPSG_list.end())
  {
    cout << "I did not find this EPSG code. " << endl;
    cout << "I am going to set to a default LCC projection" << endl;
    EPSG = 0;
    ProjParameters = EPSG_list[0]; 
    EPSG_type = "Lambert_Conformal_Conic_2SP";
  } 
  else
  {
    cout << "I found this code!! Hooray!" << endl;
    EPSG = EPSG_code;
    ProjParameters = EPSG_list[EPSG_code];
    IntProjParameters = EPSG_int_list[EPSG_code];
    EPSG_type = EPSG_type_map[EPSG_code];
  }

}

// This simply prints the projection parameters to screen.
void LSDProjectionInfo::print_projection_parameters()
{
  cout << "Your projection has the EPSG code: " << EPSG << endl;
  cout << "The type is: " << EPSG_type << endl;
  cout << "The parameters are: " << endl;

  vector<string> keys =  extract_keys(ProjParameters);
  for (int i = 0; i< int(keys.size()); i++)
  {
    cout << keys[i] << ": " << ProjParameters[ keys[i] ] << endl;
  }
}




//==================================================================
//==================================================================
// The new coordinate converter object
//==================================================================
void LSDCoordinateConverter::create()
{
  cout << "You are creating an empty coordinate converter. " << endl;
  cout << "I am going to assume you want WGS84 Ellipse and WGS84 datum." << endl;

  LSDReferenceEllipsoid AnEllipsoid("WGS84");
  RefEllipsoid = AnEllipsoid;

  int default_EPSG = 0;
  LSDProjectionInfo Proj(default_EPSG);
  ProjParams = Proj;

  RADIANS_PER_DEGREE = M_PI/180.0;
  DEGREES_PER_RADIAN = 180.0/M_PI;

}

// In this create function you need an EPSG code. 
// If it is not included then it will go to the default EPSG code.
void LSDCoordinateConverter::create(int EPSG_code)
{
  cout << "You did not give me an ellipsoid. I will use WGS84." << endl;  
  LSDReferenceEllipsoid AnEllipsoid("WGS84");
  RefEllipsoid = AnEllipsoid;

  // Now get the projection
  LSDProjectionInfo Proj(EPSG_code);
  ProjParams = Proj;

  RADIANS_PER_DEGREE = M_PI/180.0;
  DEGREES_PER_RADIAN = 180.0/M_PI;


}

// In this create function you need an EPSG code and an ellipsoid. 
// If these are not included then it will go to the default EPSG code
// and the default ellipsoid, which is WGS84.
void LSDCoordinateConverter::create(int EPSG_code, string Ellipsoid_name)
{
  cout << "I am going to create an ellipsoid using " << Ellipsoid_name << endl;
  LSDReferenceEllipsoid AnEllipsoid(Ellipsoid_name);
  RefEllipsoid = AnEllipsoid;

  // Now get the projection
  LSDProjectionInfo Proj(EPSG_code);
  ProjParams = Proj;

  RADIANS_PER_DEGREE = M_PI/180.0;
  DEGREES_PER_RADIAN = 180.0/M_PI;

}

// This converts latitude-longitude data to easting and northing data
pair<double,double> LSDCoordinateConverter::Convert_LL_to_EN(pair<double,double> Lat_Long)
{
  // check the ellipsoid and the coordinate system
  bool does_it_work_ellipse = RefEllipsoid.check_ellipsoid();
  bool does_it_work_EPSG = ProjParams.check_projection();
  
  // the easting and northing pair
  pair<double,double> E_N;

  if (does_it_work_ellipse && does_it_work_EPSG)
  {
    // figure out which type of projection it is and proceed accordingly
    if (ProjParams.get_EPSG_type() == "Albers_Conic_Equal_Area")
    {
      cout << "You want Albers Equal Area Conic, super! I've coded that up." << endl;
      E_N = AEAC_forward(Lat_Long);
    }
    else if (ProjParams.get_EPSG_type() == "Transverse_Mercator")
    {
      cout << "You want transverse mercator, super! I've coded that up." << endl;
      E_N = UTM_forward(Lat_Long);
    }
    else if (ProjParams.get_EPSG_type() == "Lambert_Conformal_Conic_2SP")
    {
      cout << "You want Lambert Conformal Conic, super! I've coded that up." << endl;
      E_N = LCC_forward(Lat_Long);
    }
    else
    {
      cout << "I didn't find a valid type. Can you double check?" << endl;
      E_N = make_pair(-9999,-9999);
    }
  }
  else
  {
    cout << "Your coordinate system or ellipsoid doesn't work! " << endl;
    cout << "This is caused by refering to either a reference ellipsoid or EPSG code that is not in our database. " << endl;
    E_N = make_pair(-9999,-9999);
  }

  return E_N;


}


// This converts latitude-longitude data to easting and northing data
pair<double,double> LSDCoordinateConverter::Convert_EN_to_LL(pair<double,double> Easting_Northing)
{
  // check the ellipsoid and the coordinate system
  bool does_it_work_ellipse = RefEllipsoid.check_ellipsoid();
  bool does_it_work_EPSG = ProjParams.check_projection();
  
  // the easting and northing pair
  pair<double,double> Lat_Long;

  if (does_it_work_ellipse && does_it_work_EPSG)
  {
    // figure out which type of projection it is and proceed accordingly
    if (ProjParams.get_EPSG_type() == "Albers_Conic_Equal_Area")
    {
      cout << "You want albers equal area conic, super! I've coded that up." << endl;
      Lat_Long = AEAC_inverse(Easting_Northing);
    }
    else if (ProjParams.get_EPSG_type() == "Transverse_Mercator")
    {
      //cout << "Not coded transverse mercator yet" << endl;
      cout << "You want transverse mercator, super! I've coded that up." << endl;
      Lat_Long = UTM_inverse(Easting_Northing);
    }
    else if (ProjParams.get_EPSG_type() == "Lambert_Conformal_Conic_2SP")
    {
      cout << "You want Lambert Conformal Conic, super! I've coded that up." << endl;
      Lat_Long = LCC_inverse(Easting_Northing);
    }
    else
    {
      cout << "I didn't find a valid type. Can you double check?" << endl;
      Lat_Long = make_pair(-9999,-9999);
    }
  }
  else
  {
    cout << "Your coordinate system or ellipsoid doesn't work! " << endl;
    cout << "This is caused by refering to either a reference ellipsoid or EPSG code that is not in our database. " << endl;
    Lat_Long = make_pair(-9999,-9999);
  }

  return Lat_Long;


}

//=================================================================================================
//=================================================================================================
// LAMBERT
// This does a lambert conical conformal forward conversion. 
pair<double,double> LSDCoordinateConverter::LCC_forward(pair<double,double> Lat_Long)
{
  // get all the necessary parameters
  double ecc = RefEllipsoid.get_ecc();
  double eccSquared = RefEllipsoid.get_eccSquared();
  double a = RefEllipsoid.get_EquatorialRadius();

  // the project parameters are in a map. 
  // The checking of the map should happen before this stage when you call the EPSG code. 
  map<string, double> PP = ProjParams.get_ProjParameters();

  // get coordinates and point location, converting to radians
  double phi_0_rad = PP["phi_0"]*RADIANS_PER_DEGREE;
  double lambda_0_rad = PP["lambda_0"]*RADIANS_PER_DEGREE;
  double phi_1_rad = PP["phi_1"]*RADIANS_PER_DEGREE;
  double phi_2_rad = PP["phi_2"]*RADIANS_PER_DEGREE;

  //Make sure the longitude is between -180.00 .. 179.9
  double Long = Lat_Long.second;
  double LongTemp = (Long+180)-int((Long+180)/360)*360-180;

  double phi_rad = Lat_Long.first*RADIANS_PER_DEGREE;
  double lambda_rad = LongTemp*RADIANS_PER_DEGREE;

  // everything below this line comes from equations on page 107 and 108 of Snyder 1987.
  // The worked example is on page 296. 
  // This is a very inefficient way of doing this because there are constants being calculated
  // here that could be computed only one time rather than for every point. 
  // However that would require a more complicated object structure and we do not expect conversion
  // to these coordinates to be rate limiting since it only occurs either at the end of analysis during
  // file printing or at the beginning of the analysis. 
  // SMM 21/01/2019
  double m_1 = cos( phi_1_rad ) / (sqrt(1 - eccSquared*( sin( phi_1_rad ) )*( sin( phi_1_rad ) )));
  double m_2 = cos( phi_2_rad ) / (sqrt(1 - eccSquared*( sin( phi_2_rad ) )*( sin( phi_2_rad ) )));

  double t_1_top = tan(M_PI*0.25 - phi_1_rad*0.5);
  double t_1_bottom = pow(((1 - ecc*sin(phi_1_rad))/(1 + ecc*sin(phi_1_rad))),ecc*0.5);
  double t_1 = t_1_top/t_1_bottom;

  double t_2_top = tan(M_PI*0.25 - phi_2_rad*0.5);
  double t_2_bottom = pow(((1 - ecc*sin(phi_2_rad))/(1 + ecc*sin(phi_2_rad))),ecc*0.5);
  double t_2 = t_2_top/t_2_bottom;

  double t_0_top = tan(M_PI*0.25 - phi_0_rad*0.5);
  double t_0_bottom = pow(((1 - ecc*sin(phi_0_rad))/(1 + ecc*sin(phi_0_rad))),ecc*0.5);
  double t_0 = t_0_top/t_0_bottom;    

  double t_top = tan(M_PI*0.25 - phi_rad*0.5);
  double t_bottom = pow(((1 - ecc*sin(phi_rad))/(1 + ecc*sin(phi_rad))),ecc*0.5);
  double t = t_top/t_bottom;    

  double n = (log(m_1)- log(m_2))/(log(t_1)-log(t_2));
  double F = m_1/ ( n*pow(t_1,n) );
  
  double rho_0 = a*F*(pow(t_0,n));
  double rho = a*F*pow(t,n);

  double theta = n*(lambda_rad-lambda_0_rad);

  double Easting = rho*sin(theta);
  double Northing = rho_0 - rho*cos(theta);

  bool verbose = false;
  if(verbose)
  {
    cout.precision(12);  
    cout << "==================================" << endl;
    cout << "The parameters are: " << endl;
    cout << "phi_0: " << PP["phi_0"] << endl;
    cout << "lambda_0: " << PP["lambda_0"] << endl;
    cout << "phi_1: " << PP["phi_1"] << endl;
    cout << "phi_2: " << PP["phi_2"] << endl;
    cout << "Equatorial radius " << a << endl;
    cout << "Eccentricity: " << ecc << endl; 
    cout << endl;
    cout << "The location of the point is: " << endl;
    cout << "latitude: " << Lat_Long.first << endl;
    cout << "longitude: " << Lat_Long.second << endl;
    cout << "And the converted point is: " << endl;
    cout << "Easting: " << Easting << endl;
    cout << "Northing: " << Northing << endl;
    cout << "==================================" << endl << endl;
  }


  pair<double,double> Easting_Northing = make_pair(Easting,Northing);
  return Easting_Northing;  

}

// This does a lambert conical conformal inverse conversion. 
pair<double,double> LSDCoordinateConverter::LCC_inverse(pair<double,double> Easting_Northing)
{
  // get all the necessary parameters
  double ecc = RefEllipsoid.get_ecc();
  double eccSquared = RefEllipsoid.get_eccSquared();
  double a = RefEllipsoid.get_EquatorialRadius();

  double x = Easting_Northing.first;
  double y = Easting_Northing.second;

  // the project parameters are in a map. 
  // The checking of the map should happen before this stage when you call the EPSG code. 
  map<string, double> PP = ProjParams.get_ProjParameters();

  // get coordinates and point location, converting to radians
  double phi_0_rad = PP["phi_0"]*RADIANS_PER_DEGREE;
  double lambda_0_rad = PP["lambda_0"]*RADIANS_PER_DEGREE;
  double phi_1_rad = PP["phi_1"]*RADIANS_PER_DEGREE;
  double phi_2_rad = PP["phi_2"]*RADIANS_PER_DEGREE;


  // everything below this line comes from equations on page 107 and 108 of Snyder 1987.
  // The worked example is on page 296. 
  // This is a very inefficient way of doing this because there are constants being calculated
  // here that could be computed only one time rather than for every point. 
  // However that would require a more complicated object structure and we do not expect conversion
  // to these coordinates to be rate limiting since it only occurs either at the end of analysis during
  // file printing or at the beginning of the analysis. 
  // SMM 21/01/2019
  double m_1 = cos( phi_1_rad ) / (sqrt(1 - eccSquared*( sin( phi_1_rad ) )*( sin( phi_1_rad ) )));
  double m_2 = cos( phi_2_rad ) / (sqrt(1 - eccSquared*( sin( phi_2_rad ) )*( sin( phi_2_rad ) )));

  double t_1_top = tan(M_PI*0.25 - phi_1_rad*0.5);
  double t_1_bottom = pow(((1 - ecc*sin(phi_1_rad))/(1 + ecc*sin(phi_1_rad))),ecc*0.5);
  double t_1 = t_1_top/t_1_bottom;

  double t_2_top = tan(M_PI*0.25 - phi_2_rad*0.5);
  double t_2_bottom = pow(((1 - ecc*sin(phi_2_rad))/(1 + ecc*sin(phi_2_rad))),ecc*0.5);
  double t_2 = t_2_top/t_2_bottom;

  double t_0_top = tan(M_PI*0.25 - phi_0_rad*0.5);
  double t_0_bottom = pow(((1 - ecc*sin(phi_0_rad))/(1 + ecc*sin(phi_0_rad))),ecc*0.5);
  double t_0 = t_0_top/t_0_bottom;    
 

  double n = (log(m_1)- log(m_2))/(log(t_1)-log(t_2));
  double F = m_1/ ( n*pow(t_1,n) );
  
  double rho_0 = a*F*(pow(t_0,n));
  
  
  // Calculating rho diverges from the forward calculation
  // This is multiplied by the sign of n
  double rho = sqrt(x*x+ (rho_0-y)*(rho_0-y));
  if (n<0)
  {
    rho = -rho;
  }

  double theta = atan(x/(rho_0-y));

  // Now calculate the latitude and longitude
  // Longitude is easy
  double Longitude = (theta/n + lambda_0_rad)*DEGREES_PER_RADIAN;

  // Latitude needs a loop
  double threshold_error = 0.00000000001;
  double t = pow(  rho/(a*F), 1/n  );
  double Lat_guess = M_PI*0.5 - 2*atan(t);
  double diff = 100000;
  double new_phi;
  double old_phi = Lat_guess;
  cout << endl << endl << "--------------------------" << endl;
  do
  {
    new_phi = 0.5*M_PI - 2*atan( t*pow((1-ecc*sin(old_phi))/(1+ecc*sin(old_phi)),ecc/2) );
    diff = fabs(new_phi-old_phi);
    old_phi = new_phi;
    cout << "Diff is: " << diff << endl;
  } while (diff > threshold_error);
  double Latitude = new_phi*DEGREES_PER_RADIAN;

  bool verbose = false;
  if(verbose)
  {
    cout.precision(12);
    cout << "==================================" << endl;
    cout << "The parameters are: " << endl;
    cout << "phi_0: " << PP["phi_0"] << endl;
    cout << "lambda_0: " << PP["lambda_0"] << endl;
    cout << "phi_1: " << PP["phi_1"] << endl;
    cout << "phi_2: " << PP["phi_2"] << endl;
    cout << "Equatorial radius " << a << endl;
    cout << "Eccentricity: " << ecc << endl; 
    cout << endl;
    cout << "The location of the point is: " << endl;
    cout << "Easting: " << Easting_Northing.first << endl;
    cout << "Northing: " << Easting_Northing.second << endl;
    cout << "And the converted point is: " << endl;
    cout << "latitude: " << Latitude << endl;
    cout << "longitude: " << Longitude << endl;
    cout << "==================================" << endl << endl;
  }


  pair<double,double> Lat_Long = make_pair(Latitude,Longitude);
  return Lat_Long;  

}




//=================================================================================================
//=================================================================================================
// ALBERS
// This does a Albers conical equal area conversion. 
// See snyder 1987 page 101
pair<double,double> LSDCoordinateConverter::AEAC_forward(pair<double,double> Lat_Long)
{
  // get all the necessary parameters
  double ecc = RefEllipsoid.get_ecc();
  double eccSquared = RefEllipsoid.get_eccSquared();
  double a = RefEllipsoid.get_EquatorialRadius();

  // the project parameters are in a map. 
  // The checking of the map should happen before this stage when you call the EPSG code. 
  map<string, double> PP = ProjParams.get_ProjParameters();

  // get coordinates and point location, converting to radians
  double phi_0_rad = PP["phi_0"]*RADIANS_PER_DEGREE;
  double lambda_0_rad = PP["lambda_0"]*RADIANS_PER_DEGREE;
  double phi_1_rad = PP["phi_1"]*RADIANS_PER_DEGREE;
  double phi_2_rad = PP["phi_2"]*RADIANS_PER_DEGREE;


  //Make sure the longitude is between -180.00 .. 179.9
  double Long = Lat_Long.second;
  double LongTemp = (Long+180)-int((Long+180)/360)*360-180;

  double phi_rad = Lat_Long.first*RADIANS_PER_DEGREE;
  double lambda_rad = LongTemp*RADIANS_PER_DEGREE;

  // everything below this line comes from equations on page 101 and 102 of Snyder 1987.
  // The worked example is on page 293. 
  // This is a very inefficient way of doing this because there are constants being calculated
  // here that could be computed only one time rather than for every point. 
  // However that would require a more complicated object structure and we do not expect conversion
  // to these coordinates to be rate limiting since it only occurs either at the end of analysis during
  // file printing or at the beginning of the analysis. 
  // SMM 22/01/2019
  double m_1 = cos( phi_1_rad ) / (sqrt(1 - eccSquared*( sin( phi_1_rad ) )*( sin( phi_1_rad ) )));
  double m_2 = cos( phi_2_rad ) / (sqrt(1 - eccSquared*( sin( phi_2_rad ) )*( sin( phi_2_rad ) )));

  double q_0_log = log( (1 - ecc*sin(phi_0_rad))/(1 + ecc*sin(phi_0_rad)));
  double q_0_first = sin(phi_0_rad)/(1 - eccSquared*( sin( phi_0_rad ) )*( sin( phi_0_rad ) ));
  double q_0 = (1-eccSquared)*(q_0_first-(1/(2*ecc))*q_0_log);

  double q_1_log = log( (1 - ecc*sin(phi_1_rad))/(1 + ecc*sin(phi_1_rad)));
  double q_1_first = sin(phi_1_rad)/(1 - eccSquared*( sin( phi_1_rad ) )*( sin( phi_1_rad ) ));
  double q_1 = (1-eccSquared)*(q_1_first-(1/(2*ecc))*q_1_log);

  double q_2_log = log( (1 - ecc*sin(phi_2_rad))/(1 + ecc*sin(phi_2_rad)));
  double q_2_first = sin(phi_2_rad)/(1 - eccSquared*( sin( phi_2_rad ) )*( sin( phi_2_rad ) ));
  double q_2 = (1-eccSquared)*(q_2_first-(1/(2*ecc))*q_2_log);

  double q_log = log( (1 - ecc*sin(phi_rad))/(1 + ecc*sin(phi_rad)));
  double q_first = sin(phi_rad)/(1 - eccSquared*( sin( phi_rad ) )*( sin( phi_rad ) ));
  double q = (1-eccSquared)*(q_first-(1/(2*ecc))*q_log);
 
  // Now for the n values
  double n = (m_1*m_1-m_2*m_2)/(q_2-q_1);
  double C = m_1*m_1+n*q_1;
  
  double rho_0 = a*sqrt(C-n*q_0)/n;


  double rho = a*sqrt(C-n*q)/n;

  double theta = n*(lambda_rad-lambda_0_rad);

  double Easting = rho*sin(theta);
  double Northing = rho_0 - rho*cos(theta);

  bool verbose = false;
  if(verbose)
  {
    cout.precision(12);
    cout << "==================================" << endl;
    cout << "The parameters are: " << endl;
    cout << "phi_0: " << PP["phi_0"] << endl;
    cout << "lambda_0: " << PP["lambda_0"] << endl;
    cout << "phi_1: " << PP["phi_1"] << endl;
    cout << "phi_2: " << PP["phi_2"] << endl;
    cout << "Equatorial radius " << a << endl;
    cout << "Eccentricity: " << ecc << endl; 
    cout << endl;
    cout << "==================================" << endl;
    cout << "intermediate steps are:" << endl;
    cout << "m_1: " << m_1 << endl;
    cout << "m_2: " << m_2 << endl;
    cout << "q_1: " << q_1 << endl;
    cout << "q_2: " << q_2 << endl;
    cout << "q_0: " << q_0 << endl;
    cout << "n: " << n << endl;
    cout << "C: " << C << endl;
    cout << "rho_0: " << rho_0 << endl;
    cout << "q: " << q << endl;
    cout << "rho: " << rho << endl;
    cout << "theta: " << theta << endl;
    cout << endl<< endl << endl;
    cout << "The location of the point is: " << endl;
    cout << "latitude: " << Lat_Long.first << endl;
    cout << "longitude: " << Lat_Long.second << endl;
    cout << "And the converted point is: " << endl;
    cout << "Easting: " << Easting << endl;
    cout << "Northing: " << Northing << endl;
    cout << "==================================" << endl << endl;
  }


  pair<double,double> Easting_Northing = make_pair(Easting,Northing);
  return Easting_Northing;  

}

// This does a albers equal area conic inverse conversion. 
pair<double,double> LSDCoordinateConverter::AEAC_inverse(pair<double,double> Easting_Northing)
{
  // get all the necessary parameters
  double ecc = RefEllipsoid.get_ecc();
  double eccSquared = RefEllipsoid.get_eccSquared();
  double a = RefEllipsoid.get_EquatorialRadius();

  double x = Easting_Northing.first;
  double y = Easting_Northing.second;

  // the project parameters are in a map. 
  // The checking of the map should happen before this stage when you call the EPSG code. 
  map<string, double> PP = ProjParams.get_ProjParameters();

  // get coordinates and point location, converting to radians
  double phi_0_rad = PP["phi_0"]*RADIANS_PER_DEGREE;
  double lambda_0_rad = PP["lambda_0"]*RADIANS_PER_DEGREE;
  double phi_1_rad = PP["phi_1"]*RADIANS_PER_DEGREE;
  double phi_2_rad = PP["phi_2"]*RADIANS_PER_DEGREE;


  // everything below this line comes from equations on page 101 and 102 of Snyder 1987.
  // The worked example is on page 293. 
  // This is a very inefficient way of doing this because there are constants being calculated
  // here that could be computed only one time rather than for every point. 
  // However that would require a more complicated object structure and we do not expect conversion
  // to these coordinates to be rate limiting since it only occurs either at the end of analysis during
  // file printing or at the beginning of the analysis. 
  // SMM 22/01/2019
  double m_1 = cos( phi_1_rad ) / (sqrt(1 - eccSquared*( sin( phi_1_rad ) )*( sin( phi_1_rad ) )));
  double m_2 = cos( phi_2_rad ) / (sqrt(1 - eccSquared*( sin( phi_2_rad ) )*( sin( phi_2_rad ) )));

  double q_0_log = log( (1 - ecc*sin(phi_0_rad))/(1 + ecc*sin(phi_0_rad)));
  double q_0_first = sin(phi_0_rad)/(1 - eccSquared*( sin( phi_0_rad ) )*( sin( phi_0_rad ) ));
  double q_0 = (1-eccSquared)*(q_0_first-(1/(2*ecc))*q_0_log);

  double q_1_log = log( (1 - ecc*sin(phi_1_rad))/(1 + ecc*sin(phi_1_rad)));
  double q_1_first = sin(phi_1_rad)/(1 - eccSquared*( sin( phi_1_rad ) )*( sin( phi_1_rad ) ));
  double q_1 = (1-eccSquared)*(q_1_first-(1/(2*ecc))*q_1_log);

  double q_2_log = log( (1 - ecc*sin(phi_2_rad))/(1 + ecc*sin(phi_2_rad)));
  double q_2_first = sin(phi_2_rad)/(1 - eccSquared*( sin( phi_2_rad ) )*( sin( phi_2_rad ) ));
  double q_2 = (1-eccSquared)*(q_2_first-(1/(2*ecc))*q_2_log);
 
  // Now for the n values
  double n = (m_1*m_1-m_2*m_2)/(q_2-q_1);
  double C = m_1*m_1+n*q_1;
  
  double rho_0 = a*sqrt(C-n*q_0)/n;

  // Now we do the components that are specific to the inverse calculation
  double rho = sqrt(x*x+(rho_0-y)*(rho_0-y));
  double q = (C-(rho*rho*n*n)/(a*a))/n;
  double theta = atan(x/(rho_0-y));

  // Now calculate the latitude and longitude
  // Longitude is easy
  double Longitude = (theta/n + lambda_0_rad)*DEGREES_PER_RADIAN;

  // Latitude needs a loop
  double threshold_error = 0.00000000001;
  double Lat_guess = asin(q*0.5);
  double diff = 100000;
  double new_phi;
  double old_phi = Lat_guess;


  cout << endl << endl << "--------------------------" << endl;
  cout << "First guess is: " << Lat_guess*DEGREES_PER_RADIAN << endl;
  do
  {

    double phi_log = (1/(2*ecc))*log( (1 - ecc*sin(old_phi))/(1 + ecc*sin(old_phi)));
    double phi_second = sin(old_phi)/(1 - eccSquared*( sin( old_phi ) )*( sin( old_phi ) ));
    double phi_first = q/(1-eccSquared);
    double phi_front = (0.5*(1-eccSquared*( sin( old_phi ) )*( sin( old_phi ) ))*(1-eccSquared*( sin( old_phi ) )*( sin( old_phi ) )))/(cos(old_phi)); 
    new_phi = old_phi+phi_front*(phi_first-phi_second+phi_log);

    //cout << "New phi is: " << new_phi*DEGREES_PER_RADIAN << endl; 

    diff = fabs(new_phi-old_phi);
    old_phi = new_phi;
    cout << "Diff is: " << diff << endl;
  } while (diff > threshold_error);
  double Latitude = new_phi*DEGREES_PER_RADIAN;

  //cout << "Final latitude is: " << new_phi*DEGREES_PER_RADIAN << endl; 

  bool verbose = true;
  //cout << "YOYOYOYOYOY" << endl;
  if(verbose)
  {
    cout.precision(12);
    cout << "==================================" << endl;
    cout << "The parameters are: " << endl;
    cout << "phi_0: " << PP["phi_0"] << endl;
    cout << "lambda_0: " << PP["lambda_0"] << endl;
    cout << "phi_1: " << PP["phi_1"] << endl;
    cout << "phi_2: " << PP["phi_2"] << endl;
    cout << "Equatorial radius " << a << endl;
    cout << "Eccentricity: " << ecc << endl; 
    cout << endl;
    cout << "The location of the point is: " << endl;
    cout << "Easting: " << Easting_Northing.first << endl;
    cout << "Northing: " << Easting_Northing.second << endl;
    cout << "And the converted point is: " << endl;
    cout << "latitude: " << Latitude << endl;
    cout << "longitude: " << Longitude << endl;
    cout << "==================================" << endl << endl;
  }


  pair<double,double> Lat_Long = make_pair(Latitude,Longitude);
  return Lat_Long;  

}




// This does a universal transverse mercator forward conversion. 
// WARNING: this forces the UTM zone!!!
pair<double,double> LSDCoordinateConverter::UTM_forward(pair<double,double> Lat_Long)
{

  // get all the necessary parameters
  //double ecc = RefEllipsoid.get_ecc();
  double eccSquared = RefEllipsoid.get_eccSquared();
  double a = RefEllipsoid.get_EquatorialRadius();

  double Long = Lat_Long.second;
  double Lat = Lat_Long.first;

  // we need k0, eccPrimeSquared
  // eccPrimeSquared comes from snyder 1987 page 61 equation (8--12)
  double eccPrimeSquared = eccSquared/(1-eccSquared);  

  double UTMEasting, UTMNorthing;
  double LongOrigin,LongOriginRad;
  double N, T, C, A, M;

  //Make sure the longitude is between -180.00 .. 179.9
  double LongTemp = (Long+180)-int((Long+180)/360)*360-180;

  double LatRad = Lat*RADIANS_PER_DEGREE;
  double LongRad = LongTemp*RADIANS_PER_DEGREE;

  // Get the parameters from the projection
  map<string, double> PP = ProjParams.get_ProjParameters();
  map<string, int> int_PP = ProjParams.get_IntProjParameters();
  int    ZoneNumber = int_PP["zone"];
  double k0 = PP["k_0"];

  // It is possible to compute the zone based on the longitude but
  // we use the EPSG determined zone instead of the longitude
  // here so that points printed from a raster do not
  // appear in different zones
  /*
  ZoneNumber = int((LongTemp + 180)/6) + 1;
  if( Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 )
  {
    ZoneNumber = 32;
  }
        // Special zones for Svalbard
  if( Lat >= 72.0 && Lat < 84.0 )
  {
    if(      LongTemp >= 0.0  && LongTemp <  9.0 ) ZoneNumber = 31;
    else if( LongTemp >= 9.0  && LongTemp < 21.0 ) ZoneNumber = 33;
    else if( LongTemp >= 21.0 && LongTemp < 33.0 ) ZoneNumber = 35;
    else if( LongTemp >= 33.0 && LongTemp < 42.0 ) ZoneNumber = 37;
  }
  // +3 puts origin in middle of zone
  //compute the UTM Zone from the latitude and longitude
  cout << "Zone is  " << ZoneNumber << endl;
  Zone = ZoneNumber;
  */


  LongOrigin = (ZoneNumber - 1)*6 - 180 + 3;
  LongOriginRad = LongOrigin * RADIANS_PER_DEGREE;
  N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad));
  T = tan(LatRad)*tan(LatRad);
  C = eccPrimeSquared*cos(LatRad)*cos(LatRad);
  A = cos(LatRad)*(LongRad-LongOriginRad);

  M = a*((1 - eccSquared/4 - 3*eccSquared*eccSquared/64
                - 5*eccSquared*eccSquared*eccSquared/256) * LatRad
               - (3*eccSquared/8 + 3*eccSquared*eccSquared/32
                  + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
               + (15*eccSquared*eccSquared/256
                  + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad)
               - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad));

  UTMEasting = (double)
          (k0*N*(A+(1-T+C)*A*A*A/6
                 + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
           + 500000.0);

  UTMNorthing = (double)
          (k0*(M+N*tan(LatRad)
               *(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                 + (61-58*T+T*T+600*C-330*eccPrimeSquared)*A*A*A*A*A*A/720)));

  if(Lat < 0)
  {
    //10000000 meter offset for southern hemisphere
    UTMNorthing += 10000000.0;
  }

  cout << "=======================" << endl;
  cout << "UTM conversion" << endl;
  cout << "latitude is: " << Lat << endl;
  cout << "longitude is: " << Long << endl;
  cout << "Easting is: " << UTMEasting << endl;
  cout << "Northing is: " << UTMNorthing << endl;
  cout << "=======================" << endl;
  cout << endl << endl << endl;
  pair<double,double> Easting_Northing = make_pair(UTMEasting,UTMNorthing);
  return Easting_Northing;  


}


pair<double,double> LSDCoordinateConverter::UTM_inverse(pair<double,double> Easting_Northing)
{

  double UTMEasting = Easting_Northing.first;
  double UTMNorthing = Easting_Northing.second;
  
  // get all the necessary parameters
  //double ecc = RefEllipsoid.get_ecc();
  double eccSquared = RefEllipsoid.get_eccSquared();
  double a = RefEllipsoid.get_EquatorialRadius();

  // we need k0, eccPrimeSquared
  // eccPrimeSquared comes from snyder 1987 page 61 equation (8--12)
  double eccPrimeSquared = eccSquared/(1-eccSquared);  

  // Get the parameters from the projection
  map<string, double> PP = ProjParams.get_ProjParameters();
  map<string, int> int_PP = ProjParams.get_IntProjParameters();
  int    ZoneNumber = int_PP["zone"];
  double k0 = PP["k_0"];
  bool isNorth = false;
  if (int_PP["is_north"]== 1)
  {
    isNorth = true;
  }


  double e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared));
  double N1, T1, C1, R1, D, M;
  double LongOrigin;
  double mu, phi1Rad;
  double x, y;

  x = UTMEasting - 500000.0; //remove 500,000 meter offset for longitude
  y = UTMNorthing;

  if(isNorth == false)
  {
    //cout << "Line 1010, you are in the Southern hemisphere!"<< endl;
    //remove 10,000,000 meter offset used for southern hemisphere
    y -= 10000000.0;
  }

        //+3 puts origin in middle of zone
  LongOrigin = (ZoneNumber - 1)*6 - 180 + 3;
  eccPrimeSquared = (eccSquared)/(1-eccSquared);

  M = y / k0;
  mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64
                   -5*eccSquared*eccSquared*eccSquared/256));

  phi1Rad = mu + ((3*e1/2-27*e1*e1*e1/32)*sin(2*mu)
                        + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
                        + (151*e1*e1*e1/96)*sin(6*mu));

  N1 = a/sqrt(1-eccSquared*sin(phi1Rad)*sin(phi1Rad));
  T1 = tan(phi1Rad)*tan(phi1Rad);
  C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad);
  R1 = a*(1-eccSquared)/pow(1-eccSquared*sin(phi1Rad)*sin(phi1Rad), 1.5);
  D = x/(N1*k0);

  double Lat = phi1Rad - ((N1*tan(phi1Rad)/R1)
                         *(D*D/2
                           -(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24
                           +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared
                             -3*C1*C1)*D*D*D*D*D*D/720));

  Lat = Lat * DEGREES_PER_RADIAN;

  double Long = ((D-(1+2*T1+C1)*D*D*D/6
                 +(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)
                 *D*D*D*D*D/120)
                / cos(phi1Rad));
  Long = LongOrigin + Long * DEGREES_PER_RADIAN;

  pair<double,double> Lat_Long = make_pair(Lat,Long);

  cout << "=======================" << endl;
  cout << "UTM conversion" << endl;
  cout << "Easting is: " << UTMEasting << endl;
  cout << "Northing is: " << UTMNorthing << endl;
  cout << "latitude is: " << Lat << endl;
  cout << "longitude is: " << Long << endl;
  cout << "=======================" << endl;
  cout << endl << endl << endl;

  return Lat_Long;  
}




#endif
