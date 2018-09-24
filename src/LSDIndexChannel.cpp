//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexChannel
// Land Surface Dynamics IndexChannel
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for retaining information of channels
//
// Data is mostly pointers to indeces withing the LSDFlowInfo object
//  For data about the actuall channel details such as elevation,
//  drainage area, etc. once needs to use the derivative class
//  LSDChannel
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



// Source code for the LSDChannelIndex object
// this object contains the node indices as well as the
// row and col indices for individual channel segments
// these can be arranged arbitrailiy accoridng to channel
// junctions or simply nodes downstream of a given node and upstream
// of another arbitrary node EndNode




//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------


#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDShapeTools.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannel_CPP
#define LSDIndexChannel_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Empty create routine that throws error (you need to give it parameters)
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create()
{
	//cout << "LSDIndexChannel You need to initialize with some parameters" << endl;
	//exit(EXIT_FAILURE);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// first create routine
// creates an index channel with just the node index of the starting and ending nodes
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(int SJN, int EJN, LSDFlowInfo& FlowInfo)
{

	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
	GeoReferencingStrings =  FlowInfo.get_GeoReferencingStrings();

	StartJunction = -1;
	EndJunction = -1;
	StartNode = SJN;
	EndNode = EJN;

  //cout << "Start node: " << StartNode << " and end node: " << EndNode << endl;
  if (StartNode == EndNode)
  {
    cout << "This appears to be a channel made up of a single pixel!" << endl;
  }


	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;


		}
		else
		{
			curr_node = receive_node;
		}
	}

	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;

	//cout << "The number of nodes is: " << NodeSequence.size() << endl;


}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// second create routine
// creates an index channel with just the node index of the starting and ending nodes
// also includes junction information
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(int SJ, int SJN, int EJ, int EJN, LSDFlowInfo& FlowInfo)
{
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
	GeoReferencingStrings =  FlowInfo.get_GeoReferencingStrings();

	StartJunction = SJ;
	EndJunction = EJ;
	StartNode = SJN;
	EndNode = EJN;

  //cout << "Start node: " << StartNode << " and end node: " << EndNode << endl;
  if (StartNode == EndNode)
  {
    cout << "This appears to be a channel made up of a single pixel!" << endl;
  }

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		//cout << "receive_node: " << receive_node << " and Endnode: " << EndNode << endl;

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;

		}
		else
		{
			curr_node = receive_node;
		}
	}
	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;

 //cout << "The number of nodes is: " << NodeSequence.size() << endl;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// creates an index channel from two X and Y coordinate points. It will snap
// the upstream point to a threshold drainage area to ensure that you are on
// a channel.  Then moves downstream until you are within a reasonable distance
// from the downstream point.
// If this create function is used then the start and end junctions are not
// assigned in the object
// reads in a vector of X and Y coords where:
// upstream_X = X_coords[0]
// downstream_X = X_coords[1]
// upstream_Y = Y_coords[0]
// downstream_Y = Y_coords[1]
//
// FJC 17/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(vector<float>& X_coords, vector<float>& Y_coords, LSDFlowInfo& FlowInfo, float threshold_area, float threshold_distance)
{
	cout << "Threshold area is: " << threshold_area << " Threshold distance is: " << threshold_distance << endl;
	// populate the protected variables
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
	GeoReferencingStrings =  FlowInfo.get_GeoReferencingStrings();

	float upstream_X = X_coords[0];
	float downstream_X = X_coords[1];
	float upstream_Y = Y_coords[0];
	float downstream_Y = Y_coords[1];

	cout << "upstream X: " << upstream_X << " downstream X: " << downstream_X << " upstream Y: " << upstream_Y << " downstream Y: " << downstream_Y << endl;

	// don't set junctions for this create routine
	StartJunction = -1;
	EndJunction = -1;

	// find the nearest channel node to the upstream X and Y coords
	int ThisNode = FlowInfo.get_node_index_of_coordinate_point(upstream_X, upstream_Y);
	cout << "Node index: " << ThisNode << endl;
	int ThisRow, ThisCol, ChanNode;
	FlowInfo.retrieve_current_row_and_col(ThisNode, ThisRow, ThisCol);
	// get the drainage area of this node
	int ThisCP = FlowInfo.retrieve_contributing_pixels_of_node(ThisNode);
	float ThisArea = ThisCP*DataResolution*DataResolution;
	bool reached_chan = false;
	if (ThisArea >= threshold_area)
	{
		// you are already at the channel
		ChanNode = ThisNode;
		NodeSequence.push_back(ThisNode);
		RowSequence.push_back(ThisRow);
		ColSequence.push_back(ThisCol);
		cout << "At the channel!" << endl;
	}
	else
	{
		while(reached_chan == false)
		{
			int ReceiverNode, ReceiverRow, ReceiverCol;
			FlowInfo.retrieve_receiver_information(ThisNode,ReceiverNode, ReceiverRow, ReceiverCol);
			if (ThisNode == ReceiverNode)
			{
				cout << "WARNING: at base level. Can't get an index channel." << endl;
				break;
			}
			int ReceiverCP = FlowInfo.retrieve_contributing_pixels_of_node(ReceiverNode);
			float ReceiverArea = ReceiverCP*DataResolution*DataResolution;
			//cout << "ReceiverNode: " << ReceiverNode << " ReceiverArea: " << ReceiverArea << endl;
			if (ReceiverArea >= threshold_area)
			{
				// reached the channel
				reached_chan = true;
				ChanNode = ReceiverNode;
				NodeSequence.push_back(ReceiverNode);
				RowSequence.push_back(ReceiverRow);
				ColSequence.push_back(ReceiverCol);
				cout << "At the channel!" << endl;
			}
			else
			{
				// move to the next receiver
				ThisNode = ReceiverNode;
			}
		}
	}

	if (reached_chan == true)
	{
		cout << "Got the starting point, now moving downstream" << endl;
		// got the nearest channel to the starting coordinate point. Now move downstream and create the index channel.
		StartNode = ChanNode;
		// get the node index of the downstream coordinate point
		int DownstreamNode = FlowInfo.get_node_index_of_coordinate_point(downstream_X, downstream_Y);
		bool reached_downstream = false;
		ThisNode = StartNode;
		while (reached_downstream == false)
		{
			int ReceiverNode, ReceiverRow, ReceiverCol;
			FlowInfo.retrieve_receiver_information(ThisNode,ReceiverNode, ReceiverRow, ReceiverCol);
			NodeSequence.push_back(ReceiverNode);
			RowSequence.push_back(ReceiverRow);
			ColSequence.push_back(ReceiverCol);
			// check the receiver node and see how far it is from the end node
			float dist_to_end = FlowInfo.get_Euclidian_distance(ReceiverNode, DownstreamNode);
			if (dist_to_end < threshold_distance)
			{
				cout << "You are within the threshold distance of the end point!" << endl;
				reached_downstream = true;
				EndNode = ReceiverNode;
			}
			else
			{
				// move downstream
				ThisNode = ReceiverNode;
			}
		}
	}
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets the row and column of the end node
// row and col are replaced
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_row_column_of_end_node(LSDFlowInfo& FlowInfo, int& row, int& col)
{
  FlowInfo.retrieve_current_row_and_col(EndNode,row, col);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the n contributing pixels at one before the final pixel.
// useful for when you want to get the basin area just above the tributary
// junction, which is usually at EndNode
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_contributing_pixels_at_penultimate_node(LSDFlowInfo& FlowInfo)
{
	if(StartNode==EndNode)
	{
		return FlowInfo.retrieve_contributing_pixels_of_node(EndNode);
	}
	else
	{
		int n_nodes = RowSequence.size();
		return FlowInfo.retrieve_contributing_pixels_of_node( NodeSequence[n_nodes-2]);
	}
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function removes the final node in the index channel. This is used for
// channels that extend to a downstream junction and the user
// wants to truncate the channel before it encounters a junction
//
// SMM 05/11/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::truncate_final_node()
{
  EndJunction = NoDataValue;
  int n_nodes = NodeSequence.size();

  EndNode = NodeSequence[n_nodes-2];
  NodeSequence.pop_back();
  RowSequence.pop_back();
  ColSequence.pop_back();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the n contributing pixels at the node in the channel (not the node index)
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_contributing_pixels_at_node(int n_node, LSDFlowInfo& FlowInfo)
{
  int n_nodes_in_channel = int(NodeSequence.size());
  if (n_node < 0)
  {
    cout << "LINE 307 LSDIndexChannel Not in channel!" << endl;
    exit(EXIT_FAILURE);
  }
  if (n_node >= n_nodes_in_channel)
  {
    cout << "LINE 312 LSDIndexChannel Not in channel!" << endl;
    exit(EXIT_FAILURE);
  }
  return FlowInfo.retrieve_contributing_pixels_of_node( NodeSequence[n_node] );
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the node index (from flowInfo object
// of a node in the channel
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_node_in_channel(int n_node)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	//cout << "N nodes in channel: " << n_nodes_in_channel << endl;
	//cout << "Node: " << n_node << endl;
	if (n_node < 0)
	{
		cout << "LINE 330 LSDIndexChannel Not in channel!" << endl;
		cout << "Using the 0 node" << endl;
		n_node = 0;
		//exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 335 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	return NodeSequence[n_node];
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the node index (from flowInfo object) and row and column indices
// into an LSDRaster or IndexRaster
// of a node in the channel
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_node_row_col_in_channel(int n_node, int& node, int& row, int& col)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	if (n_node < 0)
	{
		cout << "LINE 355 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 360 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	node = NodeSequence[n_node];
	row = RowSequence[n_node];
	col = ColSequence[n_node];
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the channel onto an index raster. Used to test where the channel is.
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDIndexChannel::print_index_channel_to_index_raster()
{
	int n_nodes_in_channel = int(NodeSequence.size());
	cout << "NRows: " << NRows << " NCols: " << NCols << endl;

	Array2D<int> Channel_array(NRows,NCols,NoDataValue);
	for(int i = 0; i<n_nodes_in_channel; i++)
	{
		//cout << "row: " << RowSequence[i] << " col: " << ColSequence[i] << endl;
		Channel_array[RowSequence[i]][ColSequence[i]]= 1;
	}

	LSDIndexRaster Channel_loc(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Channel_array,GeoReferencingStrings);
	return Channel_loc;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this appends the channel onto an index raster. Used to test where the channel is.
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::append_index_channel_to_index_raster(LSDIndexRaster& old_raster)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	cout << "NRows: " << NRows << " NCols: " << NCols << endl;

	Array2D<int> Channel_array= old_raster.get_RasterData();
	for(int i = 0; i<n_nodes_in_channel; i++)
	{
		//cout << "row: " << RowSequence[i] << " col: " << ColSequence[i] << endl;
		Channel_array[RowSequence[i]][ColSequence[i]]= 1;
	}

	LSDIndexRaster Channel_loc(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Channel_array,GeoReferencingStrings);
	old_raster = Channel_loc;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the x and y coordinates of all points in the channel
//
// FJC 17/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_coordinates_of_channel_nodes(vector<double>& X_coordinates, vector<double>& Y_coordinates, LSDFlowInfo& FlowInfo)
{
	// these are for extracting element-wise data from the channel profiles.
	int this_node, row,col;
	double latitude,longitude;
	double x_loc,y_loc;
	LSDCoordinateConverterLLandUTM Converter;

	// find the number of nodes
	int n_nodes = int(NodeSequence.size());

	if (n_nodes <= 0)
	{
		cout << "Cannot print since you have not calculated channel properties yet." << endl;
	}

	for (int n = 0; n<n_nodes; n++)
	{
		this_node = NodeSequence[n];
		FlowInfo.retrieve_current_row_and_col(this_node,row,col);
		FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);
		FlowInfo.get_x_and_y_locations(row, col, x_loc, y_loc);

		X_coordinates.push_back(x_loc);
		Y_coordinates.push_back(y_loc);
	}
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_coordinates_of_channel_nodes(vector<float>& X_coordinates, vector<float>& Y_coordinates)
{
  vector<float> X_coords_temp;
  vector<float> Y_coords_temp;

  int n_nodes_in_channel = int(NodeSequence.size());

  for (int i = 0; i < n_nodes_in_channel; ++i)
  {
    // get the x and y coords for this row and col
    int row = RowSequence[i];
    int col = ColSequence[i];
    float x = XMinimum + float(col)*DataResolution + 0.5 * DataResolution;
    float y = YMinimum + float(NRows - row)*DataResolution - 0.5*DataResolution;
    X_coords_temp.push_back(x);
    Y_coords_temp.push_back(y);
  }

  // copy to vectors
  X_coordinates = X_coords_temp;
  Y_coordinates = Y_coords_temp;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function writes an index channel to a CSV file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDIndexChannel::write_channel_to_csv(string path, string filename)
{

  ofstream chan_out;
  string chan_fname = path+filename+"_index_chan.csv";
  chan_out.open(chan_fname.c_str());

  int NNodes = NodeSequence.size();

  int this_row,this_col;
  float x,y;

  int id = 0;

  chan_out << "id,x,y,row,column" << endl;
  for(int node = 0; node<NNodes; node++)
  {
    id++;
    this_row = RowSequence[node];
    this_col = ColSequence[node];
    x = XMinimum + float(this_col)*DataResolution + 0.5*DataResolution;

    // Slightly different logic for y because the DEM starts from the top corner
    y = YMinimum + float(NRows-this_row)*DataResolution - 0.5*DataResolution;

    chan_out << id << "," << x << "," << y << "," << this_row << ","
             << this_col << endl;
  }
  chan_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function writes an index channel to a CSV file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDIndexChannel::write_channel_to_csv(string path, string filename, LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& Elevation)
{

  // these are for extracting element-wise data from the channel profiles.
  int this_node, row,col;
  double latitude,longitude;
  double x_loc,y_loc;
  LSDCoordinateConverterLLandUTM Converter;

  // find the number of nodes
  int n_nodes = int(NodeSequence.size());

  // open the data file
  ofstream chan_out;
  string chan_fname = path+filename+"_index_chan.csv";
  chan_out.open(chan_fname.c_str());
  chan_out << "id,row,column,flow_distance,elevation,latitude,longitude,x,y" << endl;
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }

  for (int n = 0; n<n_nodes; n++)
  {
    this_node = NodeSequence[n];
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);
    FlowInfo.get_x_and_y_locations(row, col, x_loc, y_loc);

    chan_out << this_node << ","
                   << row << ","
                   << col << ","
									 << FlowDistance.get_data_element(row,col) << ","
									 << Elevation.get_data_element(row,col) << ",";
    chan_out.precision(9);
    chan_out << latitude << ","
             << longitude << ",";
    chan_out.precision(9);
    chan_out << x_loc << "," << y_loc << endl;
  }
}


#endif
