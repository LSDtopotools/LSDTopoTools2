//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDJunctionNetwork
// Land Surface Dynamics ChannelNetwork
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for organizing channel routing under the Fastscape algorithm
//  (see Braun and Willett, Geomorphology 2013, v180, p 170-179)
//  It uses the algorithm to create channel junction networks
//  that can be searched for network connectivity
//
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


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDJunctionNetwork.cpp
// LSDJunctionNetwork object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Stuart Grieve, University of Edinburgh
// <your name here>
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 1.0.0		15/07/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDJunctionNetwork_CPP
#define LSDJunctionNetwork_CPP

#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDRasterInfo.hpp"
#include "LSDSpatialCSVReader.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
using namespace std;
using namespace TNT;

/// @brief the copy constructor
LSDJunctionNetwork& LSDJunctionNetwork::operator=(const LSDJunctionNetwork& rhs)
{
  if (&rhs != this)
  {
    NRows = rhs.NRows;
    NCols = rhs.NCols;
    XMinimum = rhs.XMinimum;
    YMinimum = rhs.YMinimum;
    DataResolution = rhs.DataResolution;
    NoDataValue = rhs.NoDataValue;
    NJunctions = rhs.NJunctions;
    GeoReferencingStrings = rhs.GeoReferencingStrings;

    SourcesVector = rhs.SourcesVector;
    BaseLevelJunctions = rhs.BaseLevelJunctions;
    JunctionVector  = rhs.JunctionVector;
    StreamOrderVector  = rhs.StreamOrderVector;
    BLBasinVector  = rhs.BLBasinVector;
    NDonorsVector =  rhs.NDonorsVector;
    ReceiverVector  = rhs.ReceiverVector;
    DeltaVector = rhs.DeltaVector;
    DonorStackVector = rhs.DonorStackVector;
    SVector  = rhs.SVector;
    SVectorIndex  = rhs.SVectorIndex;
    NContributingJunctions  = rhs.NContributingJunctions;

    StreamOrderArray = rhs.StreamOrderArray.copy();
    JunctionArray = rhs.JunctionArray.copy();
    JunctionIndexArray = rhs.JunctionIndexArray.copy();

  }
  return *this;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This constructor does nothing but allows copying of these objects
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::create( void )
{
  cout << "I am an empty LSDJunctionNetwork" << endl;

  vector<int> emptyvec;

  SourcesVector = emptyvec;
  BaseLevelJunctions = emptyvec;
  JunctionVector  = emptyvec;
  StreamOrderVector  = emptyvec;
  BLBasinVector  = emptyvec;
  NDonorsVector =  emptyvec;
  ReceiverVector  = emptyvec;
  DeltaVector = emptyvec;
  DonorStackVector = emptyvec;
  SVector  = emptyvec;
  SVectorIndex  = emptyvec;
  NContributingJunctions  = emptyvec;

  Array2D<int> emptyarray(0,0);
  StreamOrderArray = emptyarray.copy();
  JunctionArray = emptyarray.copy();
  JunctionIndexArray = emptyarray.copy();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// create
// this defines a channel network based on a FlowInfo object
// and a list of source nodes
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::create(vector<int> Sources, LSDFlowInfo& FlowInfo)
{
  NRows = FlowInfo.NRows;
  NCols = FlowInfo.NCols;
  XMinimum = FlowInfo.XMinimum;
  YMinimum = FlowInfo.YMinimum;
  DataResolution = FlowInfo.DataResolution;
  NoDataValue = FlowInfo.NoDataValue;
  GeoReferencingStrings =  FlowInfo.GeoReferencingStrings;

  SourcesVector = Sources;

  // start arrays where the data all begins as nodata
  Array2D<int> TempLinkArray(NRows,NCols,NoDataValue);

  JunctionArray = TempLinkArray.copy();
  StreamOrderArray = TempLinkArray.copy();
  JunctionIndexArray = TempLinkArray.copy();

  vector<int> TempVector;
  JunctionVector = TempVector;
  BaseLevelJunctions = TempVector;

  // we loop through the sources file
  // for each source we burn down to either base level
  // or the next channel
  //
  // with link numbers by adding integers
  //
  int n_sources = SourcesVector.size();
  int current_node;
  int current_row,current_col;
  int receiver_node;
  int baselevel_switch;		// 0 if not a base level node, 1 if so
  int current_stream_order;
  int junction_switch;
  int donor_node, donor_row,donor_col;
  int n_current_stream_order_donors;

  // loop through sources.
  // this loop sets the stream orders and identifies the
  // junctions
  // it is the first of two loops through the stream network
  // the second will generate the stream link information
  for(int src = 0; src<n_sources; src++)
  {
    baselevel_switch =0;    // 0 == not base level
    junction_switch = 0;    // 0 == no junction so far
    current_node = SourcesVector[src];
    current_row = FlowInfo.RowIndex[current_node];
    current_col = FlowInfo.ColIndex[current_node];
    receiver_node = FlowInfo.ReceiverVector[current_node];

    current_stream_order = 1;

    if (current_node == receiver_node)
    {
      baselevel_switch = 1;
    }


    // follow node down through receivers until it hits a base level node
    // the switch is set to 2 because if there is a baselevel node
    // the algorithm has to go through the loop once to register the
    // 'downslope node' as the baselevel
    while ( baselevel_switch <2 )
    {

      // check if this node has already been added to the channel network
      // if it hasn't, then this becomes a channel of order of the current order
      if(StreamOrderArray[current_row][current_col] == NoDataValue)
      {
      	StreamOrderArray[current_row][current_col] = current_stream_order;
      }
      // if it isn't a nodata node:
      // note that the junction switch starts out as a zero.
      // the channel is followed looking at nodata values in the stream order array
      // once it hits the first junction, the StreamOrderArray has finite values,
      // so this logic is triggered.
      else if(StreamOrderArray[current_row][current_col] != NoDataValue)
      {

        // each source contributes a junction unless it is the
        // only source to go to a given baselevel node
        // check to see if this is the first time the channel
        // has hit another channel. If so, add the junction
        // and set the junction switch to 1 so that no further
        // junctions are added
        // also if it is a junction check to see if the stream order has been incremented
        // if it has not, the base level switch is turn to 2 and the
        // algorithm exits this loop and goes to the next source
        // if it has, this and all downstream nodes take on the new stream order
        // if junction switch is zero it means this is the first visit of a previously visited channel
        if (junction_switch == 0)
        {
          junction_switch = 1;
          JunctionArray[current_row][current_col]  = 1;

          // if it is the the first junction for this source, the current_stream_order
          // is one. Therefore any junction will result in a stream order
          // of at least 2.
          // If the junction is currently at a stream order of 1, then it
          // gets incremented
          if (StreamOrderArray[current_row][current_col] == current_stream_order)
          {
            current_stream_order ++;
            StreamOrderArray[current_row][current_col]= current_stream_order;
          }
          // if the junction is two or greater, the loop exits since there can
          // be no more incrementing of the stream order
          else
          {
            baselevel_switch = 2;
          }
        }
        else
        {
          // first, we check to see if it is not a junction. if not, we update the
          // stream order. If the stream order hasn't changed, then something
          // is amiss since there is no point moving downstream
          // when the current stream order is the same as the previous stream order,
          // since it can't increase stream order downstream
          // nodes following downstream will be at the current stream order
          if (JunctionArray[current_row][current_col]  != 1)
          {
            // THIS IS NOT A JUNCTION
            // if the current stream order is bigger than the existing stream order
            // at this point, then increase the existing stream order
            if ( current_stream_order > StreamOrderArray[current_row][current_col])
            {
              StreamOrderArray[current_row][current_col] = current_stream_order;
            }
            else
            {
              baselevel_switch = 2;
            }

          }
          // if it is a junction, see if there is an increment
          // in the current stream order
          // the node it has come from has an updated current stream order
          // so just look at the donor nodes to
          // see if there are two or more
          // donors at a higher stream order
          else if (JunctionArray[current_row][current_col]  == 1)
          {
            // THIS IS A JUNCTION
            // first, check to see if this has a higher stream order than the current stream
            // order
            if (StreamOrderArray[current_row][current_col] > current_stream_order)
            {
              // yes, the stream order at this junction is higher than the current stream order
              // even if the junction is incremented the downstream search would
              // stop here, so get out with the baselevel switch
              baselevel_switch = 2;
            }
            else if ( StreamOrderArray[current_row][current_col] == current_stream_order)
            {


              // this means that the current stream order is equal to or less than the streamorder
              // currently at this node this means you need to check and see if the thing is incremented
              n_current_stream_order_donors = 0;
              for(int dnode = 0; dnode<FlowInfo.NDonorsVector[current_node]; dnode++)
              {
                donor_node = FlowInfo.DonorStackVector[ FlowInfo.DeltaVector[current_node]+dnode];

                // ignore the base level donor
                if (donor_node != current_node)
                {
                  donor_row = FlowInfo.RowIndex[ donor_node ];
                  donor_col = FlowInfo.ColIndex[ donor_node ];

                  if (StreamOrderArray[donor_row][donor_col] == current_stream_order)
                  {
                    n_current_stream_order_donors++;
                  }
                }
              }

              // now check to see if the stream order has increased
              if (n_current_stream_order_donors >= 2)
              {
                current_stream_order++;
                StreamOrderArray[current_row][current_col] = current_stream_order;
              }
              else    // if it hasn't, the loop ends here.
              {
                baselevel_switch = 2;
              }

            }
            else if (StreamOrderArray[current_row][current_col] < current_stream_order)
            {

              // the current stream order is higher than the stream order at this point.
              // the node needs to update its stream order and keep going
              StreamOrderArray[current_row][current_col] = current_stream_order;
            }
            else
            {
              cout << "something about the logic has gone wrong. " <<endl;
              cout << "hhh node: " << current_node << ", current_stream_order " << current_stream_order
                   << " array: " << StreamOrderArray[current_row][current_col]<<" and src: " << SourcesVector[src] << endl;
            }

          }  // end logic for if this is a downstream junction
        }    // end logic for if this is not the first junction
      }      // end logic for if this is not a NoData node


      // get the next current node, which is this nodes receiver
      current_node = FlowInfo.ReceiverVector[current_node];

      // get the next receiver node, which is the next node
      receiver_node = FlowInfo.ReceiverVector[current_node];
      current_row = FlowInfo.RowIndex[current_node];
      current_col = FlowInfo.ColIndex[current_node];

      // if this is a baselevel node
      if (current_node == receiver_node)
      {
        baselevel_switch ++;
      }


    }    // end flow to baselevel loop

  }      // end sources loop


  // now you need to loop through the sources once more, creating links
  // each link has a starting node, and ending node
  // a stream order
  // a starting link and an ending link
  //
  // this should be arranged in an analagous way to the fastscape algorithm
  // all sources are on 1st order links
  int this_junction = -1;
  for(int src = 0; src<n_sources; src++)
  {
    this_junction++;				// increment the last junction
    baselevel_switch =0;			// 0 == not base level
    junction_switch = 0;			// 0 == no junction so far
    current_node = SourcesVector[src];
    current_row = FlowInfo.RowIndex[current_node];
    current_col = FlowInfo.ColIndex[current_node];
    receiver_node = FlowInfo.ReceiverVector[current_node];

    //cout << "LINE 257 ChNet, SOURCE: " << src <<  " n_src: " << n_sources << " current_node: " << current_node
    //     << " and rnode: " << receiver_node << endl;

    //each source is a junction node. Push back the junction vector
    JunctionVector.push_back(current_node);

    // set the junction Index Array
    JunctionIndexArray[current_row][current_col] = this_junction;

    // stream order only increases at junctions. So the junction node has a stream
    // order that remains the same until it gets to the next junction, where it possibly
    // could change
    StreamOrderVector.push_back( StreamOrderArray[current_row][current_col] );

    // check if this is a baselevel node
    if(receiver_node == current_node)
    {
      baselevel_switch = 1;			// turn the baselevel switch on
      ReceiverVector.push_back(this_junction);		// the Receiver node is iteself

      // this logic only applies to sources, which cannot lie downstream of another source
      // ***THIS MUST BE ENFORCED BY THE GET_SOURCES ALGORITHM***
      // this means that if a source is also a baselevel, then this is the one and only time
      // this baselevel junction is added, so add it to the baselevel vector
      BaseLevelJunctions.push_back(this_junction);
    }

    // the next element is the receiver junction, we need to follow the path to this receiver
    // the routine continues until the junction has been visited or until it hits
    // a baselevel node
    //cout << "LINE 280" << endl;
    while (baselevel_switch == 0 && junction_switch <2)
    {
      //cout << "Line 283" << endl;

      //cout << "Line 286, current node = " << current_node << " and rode: " << receiver_node << endl;
      current_node = receiver_node;
      //cout << "Line 288, current node = " << current_node << " and rode: " << receiver_node << endl;
      current_row = FlowInfo.RowIndex[current_node];
      current_col = FlowInfo.ColIndex[current_node];
      receiver_node = FlowInfo.ReceiverVector[current_node];

      // first we need logic for if this is a baselevel node
      if (current_node == receiver_node)
      {
        //cout << "source: " << src << " and BASELEVEL, node: " << current_node << " rnode: " << receiver_node << endl;
        // check to see if it has a junction index number.
        if(JunctionIndexArray[current_row][current_col] == NoDataValue)
        {
        	// it doens't have a JunctionIndexNumber. This is a new
        	// junction
        	this_junction++;

        	// this junction has the this_junction index. Set the JunctionIndexArray
        	JunctionIndexArray[current_row][current_col] = this_junction;

          // the receiver node of the previous junction is the new junction
          ReceiverVector.push_back( JunctionIndexArray[current_row][current_col] );

          //push back the junction vector
          JunctionVector.push_back(current_node);

          // because this is a baselevel node, the Receiver of this junction
          // is iteself
          ReceiverVector.push_back( JunctionIndexArray[current_row][current_col] );

          // the stream order of this node is also determined by the node
          StreamOrderVector.push_back( StreamOrderArray[current_row][current_col] );

          // finally, this is the first time we have visted this baselevel node.
          // So it gets added to the baselevel vector
          BaseLevelJunctions.push_back(this_junction);

        }
        else    // this junction does have an index number, no new junction is created
        {
          // the receiver node of the previous junction is the new junction
          ReceiverVector.push_back( JunctionIndexArray[current_row][current_col] );
        }
        junction_switch = 2;
        baselevel_switch = 1;      // this is a baselevel. It will exit the
                                    // loop and move to the next source
      }
      else                  // this is not a baselevel node
      {
        //cout << "LINE 330, not baselevel; src: " << src << " and node: " << current_node << " rnde: " << receiver_node << endl;
        // the node in the junction array is zero if it is not a
        // junction, 1 if it is an unvisited junction, and 2 or more if it
        // is a visited junction
        if(JunctionArray[current_row][current_col] != NoDataValue)
        {
          //cout << "LINE 338, found a junction at node: " << current_node
          //	 << " JArray: " << JunctionArray[current_row][current_col]  << endl;
          junction_switch = JunctionArray[current_row][current_col];
          JunctionArray[current_row][current_col] ++;		// increment the junction array
                        // it will be greater than 1 if
                        // the junction has been visited

          // if this junction has been visited, it will have a junction number
          // include the receiver vector
          if (JunctionIndexArray[current_row][current_col] != NoDataValue )
          {
            ReceiverVector.push_back( JunctionIndexArray[current_row][current_col] );

            // the loop will not continue; it will move onto the next
            // source since it has visited an already visited junction
          }
          else    // the junction has not been visited
          {
            // this is a new junction. Increment the 'last junction' int
            this_junction++;

            // this junction has the this_junction index. Set the JunctionIndexArray
            JunctionIndexArray[current_row][current_col] = this_junction;

            // the receiver node of the previous junction is the new junction
            ReceiverVector.push_back( JunctionIndexArray[current_row][current_col] );

            //push back the junction vector; this is a new junction
            JunctionVector.push_back(current_node);

            // get the stream order of this new junction
            StreamOrderVector.push_back( StreamOrderArray[current_row][current_col] );
          }
        }   // end logic for is this a junction
      }     // end logic for not a baselevel node
    }       // end baselevel logic
  }         // end sources loop

  //cout << "ChanNet; LINE 368; sz ReceiverVec: " << ReceiverVector.size() << " sz JuncVec: " << JunctionVector.size()
  //   << " sz SOVec: " << StreamOrderVector.size() << endl;

  //for(int i = 0; i< int(BaseLevelJunctions.size()); i++)
  //{
  //	cout << "LINE 382 bl node["<<i<<"]: " << BaseLevelJunctions[i] << endl;
  //}

  //cout << "LINE 385: links data " << endl;
  //for (int i = 0; i< int(StreamOrderVector.size()); i++)
  //{
  //	cout << "Junc: " << i << " node: " << JunctionVector[i]
  //	     << " receiv: " << ReceiverVector[i] << " Order: " << StreamOrderVector[i] << endl;
  //}


  // get the number of junctions
  NJunctions = int(JunctionVector.size());

  // now we implement the fastscape algorithm


  // set the sizes of the member vectors
  vector<int> ndn_vec(NJunctions,0);
  vector<int> ndn_nodata_vec(NJunctions,NoDataValue);
  vector<int> ndn_plusone_vec(NJunctions+1,0);
  vector<int> w_vector(NJunctions,0);

  NDonorsVector = ndn_vec;
  DonorStackVector = ndn_vec;
  DeltaVector = ndn_plusone_vec;

  SVector = ndn_nodata_vec;
  BLBasinVector = ndn_nodata_vec;

  // first create the number of donors vector
  // from braun and willett eq. 5
  for(int i = 0; i<NJunctions; i++)
  {
    NDonorsVector[ ReceiverVector[i] ]++;
  }

  // now create the delta vector
  // this starts on the last element and works its way backwards
  // from Braun and Willett eq 7 and 8
  DeltaVector[NJunctions] = NJunctions;
  for(int i = NJunctions; i>0; i--)
  {
    DeltaVector[i-1] = DeltaVector[i] -  NDonorsVector[i-1];
  }

  // now the DonorStack and the r vectors. These come from Braun and Willett
  // equation 9.
  // Note that in the manscript I have there is a typo in eqaution 9
  // (Jean Braun's code is correct)
  // it should be w_{r_i} = w_{r_i}+1
  int r_index;
  int w_index;
  int delta_index;
  for (int i = 0; i<NJunctions; i++)
  {
    r_index = ReceiverVector[i];
    delta_index = DeltaVector[ r_index ];
    w_index = w_vector[ r_index ];
    DonorStackVector[  delta_index+w_index ] = i;
    w_vector[r_index] += 1;
  }

  // now go through the base level node list, building the drainage tree for each of these nodes as one goes along
  int n_base_level_nodes;
  n_base_level_nodes = BaseLevelJunctions.size();

  int k;
  int j_index;
  int begin_delta_index, end_delta_index;
  int l_index;

  j_index = 0;
  for (int i = 0; i<n_base_level_nodes; i++)
  {
    k = BaseLevelJunctions[i];			// set k to the base level node

    // This doesn't seem to be in Braun and Willet but to get the ordering correct you
    // need to make sure that the base level node appears first in the donorstack
    // of nodes contributing to the baselevel node.
    // For example, if base level node is 4, with 4 donors
    // and the donor stack has 3 4 8 9
    // the code has to put the 4 first.
    if (DonorStackVector[ DeltaVector[k] ] != k)
    {
      int this_index = DonorStackVector[ DeltaVector[k] ];
      int bs_node = k;

      for(int ds_node = 1; ds_node < NDonorsVector[k]; ds_node++)
      {
        if( DonorStackVector[ DeltaVector[k] + ds_node ] == bs_node )
        {
          DonorStackVector[ DeltaVector[k] ] = k;
          DonorStackVector[ DeltaVector[k] + ds_node ] = this_index;
        }
      }
    }

    // now run recursive algorithm
    begin_delta_index = DeltaVector[k];
    end_delta_index = DeltaVector[k+1];

    for (int delta_index = begin_delta_index; delta_index<end_delta_index; delta_index++)
    {
      l_index = DonorStackVector[delta_index];
      add_to_stack(l_index, j_index, k);
    }
  }

  // now run the indexing and accumulation routine
  vector<int> vectorized_contributing_pixels(NJunctions,1);
  SVectorIndex = vectorized_contributing_pixels;

  int receiver_junction;
  int donor_junction;

  // loop through the s vector, adding pixels to receiver nodes
  for(int junc = NJunctions-1; junc>=0; junc--)
  {
    donor_junction = SVector[junc];
    receiver_junction = ReceiverVector[ donor_junction ];

    // every node is visited once and only once so we can map the
    // unique positions of the nodes to the SVector
    SVectorIndex[donor_junction] = junc;

    // add the upslope area (note no action is taken
    // for base level nodes since they donate to themselves and
    // we must avoid float counting
    if (donor_junction != receiver_junction)
    {
      vectorized_contributing_pixels[ receiver_junction ] +=  vectorized_contributing_pixels[ donor_junction ];
    }
  }
  //cout << "LINE 525 did area calcs " << endl;

  NContributingJunctions = vectorized_contributing_pixels;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_UTM_information(int& UTM_zone, bool& is_North)
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
// This function returns the x and y location of a row and column
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc)
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
void LSDJunctionNetwork::get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc)
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
void LSDJunctionNetwork::get_lat_and_long_locations(int row, int col, double& lat,
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

void LSDJunctionNetwork::get_x_and_y_from_latlong(vector<float> latitude, vector<float> longitude,
                                                   vector<float>& UTME,vector<float>& UTMN)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  int N_samples =  int(latitude.size());

  // set up some temporary vectors
  vector<float> this_UTMN(N_samples,0);
  vector<float> this_UTME(N_samples,0);

  double this_Northing;
  double this_Easting;

  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);


  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    cout << "Converting point " << i << " to UTM." << endl;
    Converter.LLtoUTM_ForceZone(eId, latitude[i], longitude[i],
                      this_Northing, this_Easting, UTM_zone);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
    cout << "Easting: " << this_Easting << " and northing: " << this_Northing << endl;
  }

  UTME = this_UTME;
  UTMN = this_UTMN;


}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// recursive add_to_stack routine, from Braun and Willett eq. 12 and 13
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::add_to_stack(int lm_index, int& j_index, int bl_node)
{
  //cout << "j_index: " << j_index << " and s_vec: " << lm_index << endl;

  SVector[j_index] = lm_index;
  BLBasinVector[j_index] = bl_node;
  j_index++;


  int begin_m,end_m;
  int l_index;
  // if donating to itself, need escape hatch
  if ( lm_index == bl_node)
  {
    begin_m = 0;
    end_m = 0;
  }
  else
  {
    begin_m = DeltaVector[lm_index];
    end_m =  DeltaVector[ lm_index+1];
  }
  //cout << "lm_index: " << lm_index << " begin_m: " << begin_m << " end m: " << end_m << endl;
  for( int m_index = begin_m; m_index<end_m; m_index++)
  {
    //cout << "recursion, begin_m: " << begin_m << " and end_m: " << end_m << endl;
    l_index = DonorStackVector[m_index];
    add_to_stack(l_index, j_index, bl_node);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function returns a integer vector containing all the junction numbers upslope
// of of the junction with number junction_number_outlet
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::get_upslope_junctions(int junction_number_outlet)
{
  vector<int> us_junctions;

  if(junction_number_outlet < 0 || junction_number_outlet > NJunctions-1)
  {
    cout << "Tried LSDJunctionNetwork::get_upslope_junctions but the"
         << "  junction number does not exist" << endl;
    exit(0);
  }

  int start_SVector_junction = SVectorIndex[junction_number_outlet];
  int end_SVector_junction = start_SVector_junction+NContributingJunctions[junction_number_outlet];

  for(int junction = start_SVector_junction; junction < end_SVector_junction; junction++)
  {
    us_junctions.push_back(SVector[junction]);
  }

  return us_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function returns a integer vector containing all the junction numbers upslope
// of of the junction with number junction_number_outlet of a specified order
//
// FJC 13/08/19
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::get_upslope_junctions_by_order(int junction_number_outlet, int stream_order)
{
  vector<int> us_junctions;

  if(junction_number_outlet < 0 || junction_number_outlet > NJunctions-1)
  {
    cout << "Tried LSDJunctionNetwork::get_upslope_junctions but the"
         << "  junction number does not exist" << endl;
    exit(0);
  }

  int start_SVector_junction = SVectorIndex[junction_number_outlet];
  int end_SVector_junction = start_SVector_junction+NContributingJunctions[junction_number_outlet];

  for(int junction = start_SVector_junction; junction < end_SVector_junction; junction++)
  {
    int this_jn = SVector[junction];
    int this_SO = get_StreamOrder_of_Junction(this_jn);
    if (this_SO == stream_order)
    {
      us_junctions.push_back(SVector[junction]);
    }
  }

  return us_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a junction and finds all the source junction upstream of the
// junction.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_all_source_junctions_of_an_outlet_junction(int junction_number_outlet)
{
  vector<int> us_junctions = get_upslope_junctions(junction_number_outlet);
  vector<int> source_junctions;

  int n_upslope_junctions = int(us_junctions.size());
  for (int j = 0; j<n_upslope_junctions; j++)
  {
    // if the junction has no donors, it is a source
    if (NDonorsVector[ us_junctions[j] ] == 0)
    {
      source_junctions.push_back(us_junctions[j]);
    }
  }
  return source_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a junction and finds all the source nodes
// from the flowinfo nodefile upstream of the
// outlet junction.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_all_source_nodes_of_an_outlet_junction(int junction_number_outlet)
{
  vector<int> us_junctions = get_upslope_junctions(junction_number_outlet);
  vector<int> source_nodes;

  int n_upslope_junctions = int(us_junctions.size());
  for (int j = 0; j<n_upslope_junctions; j++)
  {
    // if the junction has no donors, it is a source
    if (NDonorsVector[ us_junctions[j] ] == 0)
    {
      source_nodes.push_back(JunctionVector[ us_junctions[j] ]);
    }
  }
  return source_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns a vector of junction indices to all the donor
// junctions of a particular junction
//
// IMPORTANT: this has only retained the string "node" to keep equivalence
//  with the FlowInfo object. It takes junctions and returns junctions!!
// Also note that base level nodes have themselves as a donor
//
// SMM 16/6/2015
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_donor_nodes(int current_node)
{
  int start_D = DeltaVector[current_node];
  int end_D = DeltaVector[current_node+1];

  vector<int> donor_nodes;
  for(int this_node = start_D; this_node<end_D; this_node++)
  {
    //cout << "node " << current_node << " and donor: " << DonorStackVector[ this_node ] << endl;
    donor_nodes.push_back( DonorStackVector[ this_node ] );
  }

  return donor_nodes;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// If you get an upslope junction list the indexing starts at the furthest downslope
// junction. All of the junction pointing refers to the master junction list however.
//
// This function maps a junction onto the indexing of the upslope junction list
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::map_junction_to_upslope_junction_list(vector<int> upslope_junctions, int junction)
{
  // get the SVector location of the first junction
  int start_SVector_junction = SVectorIndex[ upslope_junctions[0] ];

  // now find how many elements away in the upslope_junction list
  // the current junction is from the start junction
  int mapped_us_junction_index = SVectorIndex[ junction ] - start_SVector_junction;

  if(mapped_us_junction_index < 0 || mapped_us_junction_index > int(upslope_junctions.size())-1)
  {
    cout << "Tried LSDJunctionNetwork::map_junction_to_upslope_junction_list"
         << "  junction number is not within the list of upslope junctions" << endl;
    exit(EXIT_FAILURE);
  }


	return mapped_us_junction_index;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the maximum stream order
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::get_maximum_stream_order()
{
  int max_stream_order = 0;
  for (int i = 0; i<NJunctions; i++)
  {
    if ( max_stream_order < StreamOrderVector[i])
    {
      max_stream_order = StreamOrderVector[i];
    }
  }
  return max_stream_order;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the number of streams of a given stream order
//
// FJC 15/03/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::get_number_of_streams(LSDFlowInfo& FlowInfo, int stream_order)
{
  int count = 0;
  for (int CurrentJN = 0; CurrentJN<NJunctions; ++CurrentJN)
  {
    int CurrentSO = get_StreamOrder_of_Junction(FlowInfo, CurrentJN);
    int ReceiverJN = get_Receiver_of_Junction(CurrentJN);
    if (CurrentSO == stream_order && CurrentJN != ReceiverJN)
    {
      int ReceiverSO = get_StreamOrder_of_Junction(FlowInfo, ReceiverJN);
      //check that you have increased SO downstream
      if (ReceiverSO > CurrentSO)
      {
        //reached end of stream segment, count
        count++;
      }
    }
  }
  return count;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the junction angles
// returns a map of float vectors. The float vector actually contains a float
// and 3 ints converted to float: the angle and the stream order of the junction
// and its 2 donors.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map<int, vector<float> > LSDJunctionNetwork::calculate_junction_angles(vector<int> JunctionList, LSDFlowInfo& FlowInfo)
{
  map<int, vector<float> > map_of_junction_angles;
  vector<float> temp_junctioninfo;
  vector<float> four_element(4,0);

  int NJuncs = int(JunctionList.size());
  if (NJuncs == 0)
  {
    NJuncs = get_NJunctions();
    cout << "You gave me an empty junction list, I am getting angles for all junctions." << endl;
    vector<int> JL;
    for(int i = 0; i<NJunctions; i++)
    {
      JL.push_back(i);
    }
    JunctionList = JL;
  }

  // now go through each junction getting the angles
  //vector<float> JunctionAngles;
  vector<int> donors;
  float this_angle;
  bool is_baselevel;
  bool channel_points_downstream = true;  // this is needed to get the correct
                                          // orientation of the channels
  int end_node, start_node1, start_node2;
  int junction_order, donor1_order, donor2_order;

  // vectors for holding the channel locations
  vector<float> x1, x2, y1, y2;
  int this_junc;

  // loop through junctions
  for(int junc = 0; junc < NJuncs; junc++ )
  {

    this_junc = JunctionList[junc];
    //cout << endl << "==================" << endl << "Junction is: " << this_junc << endl;

    // check to see if the junction exists
    if (this_junc >= NJunctions)
    {
      cout << "FATAL ERROR LSDJunctionNetwork::calculate_junction_angles" << endl;
      cout << "You have called a junction that doesn't exist." << endl;
      exit(EXIT_FAILURE);
    }

    // check if it is a baselevel node
    int ReceiverJN = get_Receiver_of_Junction(this_junc);
    is_baselevel = (this_junc == ReceiverJN) ? true : false;

    // if not a baselevel see if it has donors
    if (is_baselevel)
    {
      //cout << "This is a baselevel junction." << endl;
      //JunctionAngles.push_back(NoDataValue);
    }
    else
    {
      // check the junction to see if it has two or more donors
      donors = get_donor_nodes(this_junc);

      // it has donors
      if( int(donors.size()) >= 2)
      {
        if ( int(donors.size()) > 2)
        {
          cout << "Warning, this junction is a weirdo and has more than two donors. I am going to use the first two donors." << endl;
        }


        // reset the temp vector
        temp_junctioninfo = four_element;

        //cout << "The donor junctions are: " << donors[0] << ", " << donors[1] << endl;
        // now get the two segments.
        // The ending node is the current junction, the starting nodes are the
        // two donor junctions.
        end_node = get_Node_of_Junction(this_junc);
        start_node1 = get_Node_of_Junction(donors[0]);
        start_node2 = get_Node_of_Junction(donors[1]);

        // extract the channel information
        LSDIndexChannel c1(start_node1, end_node,FlowInfo);
        LSDIndexChannel c2(start_node2, end_node,FlowInfo);

        // now get the locations of the nodes in the channels in x,y coordinates
        c1.get_coordinates_of_channel_nodes(x1, y1);
        c2.get_coordinates_of_channel_nodes(x2, y2);

        junction_order = get_StreamOrder_of_Junction(FlowInfo,this_junc);
        donor1_order = get_StreamOrder_of_Junction(FlowInfo,donors[0]);
        donor2_order = get_StreamOrder_of_Junction(FlowInfo,donors[1]);

        temp_junctioninfo[1] = float(junction_order);
        temp_junctioninfo[2] = float(donor1_order);
        temp_junctioninfo[3] = float(donor2_order);

        // now calculate the angle
        this_angle = angle_between_two_vector_datasets(x1, y1,x2, y2,channel_points_downstream);
        this_angle = fabs(this_angle);
        //JunctionAngles.push_back(this_angle);
        //cout << "Angle is: " << this_angle << " radians, which is " << deg(this_angle) << " degrees." << endl;
        //cout << "Stream order is " << junction_order << " with donor 1: " << donor1_order << " and donor2: " << donor2_order << endl;
        temp_junctioninfo[0] = this_angle;

        map_of_junction_angles[this_junc] = temp_junctioninfo;
      }
      else
      {
        //cout << "This junction doesn't have 2 donors; it must be a source." << endl;
        //JunctionAngles.push_back(NoDataValue);
      }
    }
  }


  return map_of_junction_angles;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the junction angles
// It is a more complete version of the junction angle code
// The vector of floats actually contains a number of elements
// The int map has:
// [0] Stream order junction 1
// [1] Stream order junction 2
// [2] Stream order receiver junction
//
// The float map has
// [0] drainage area junction 1
// [1] drainage area junction 2
// [2] drainage area receiver junction
// [3] junction angle J1-J2
// [4] junction angle J1-R
// [5] junction angle J2-R
// [6] R^2 J1 segment
// [7] R^2 J2 segment
// [8] R^2 R segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::calculate_junction_angles_complete(vector<int> JunctionList,
                                                            LSDFlowInfo& FlowInfo,
                                                            map<int , vector<int> >& JA_int_info,
                                                            map<int, vector<float> >& JA_float_info )
{
  map<int, vector<float> > JA_float_map;
  map<int, vector<int> > JA_int_map;
  vector<float> temp_float_junctioninfo;
  vector<int> temp_int_junctioninfo;
  vector<float> nine_element_float(15,0);
  vector<int> four_element_int(4,0);

  int NJuncs = int(JunctionList.size());
  if (NJuncs == 0)
  {
    NJuncs = get_NJunctions();
    cout << "You gave me an empty junction list, I am getting angles for all junctions." << endl;
    vector<int> JL;
    for(int i = 0; i<NJunctions; i++)
    {
      JL.push_back(i);
    }
    JunctionList = JL;
  }

  // now go through each junction getting the angles
  //vector<float> JunctionAngles;
  vector<int> donors;
  float this_angle_D1_D2, this_angle_D1_R, this_angle_D2_R;
  bool is_baselevel;
  bool channel_points_downstream = true;  // this is needed to get the correct
                                          // orientation of the channels
  int end_node, start_node1, start_node2;
  int junction_order, donor1_order, donor2_order;
  int receiver_node, receiver_order;


  // vectors for holding the channel locations
  vector<float> x1, x2, y1, y2, xr, yr;


  // Some floats for holding vector information
  vector<float> D1_vec, D2_vec, R_vec;  // These are the vebor elements of the two donors and receivers.
                                        // They have an origin at 0,0, and the end point is at these coordinates.
  float intercept,gradient;
  float D1_R2, D2_R2, R_R2;

  // The hold drainage areas
  float DA1, DA2, DAR;

  float D1_bearing, D2_bearing, R_bearing;
  float D1_rotated_bearing, D2_rotated_bearing;
  float junction_sum;

  int this_junc;

  // loop through junctions
  for(int junc = 0; junc < NJuncs; junc++ )
  {

    this_junc = JunctionList[junc];
    //cout << endl << "==================" << endl << "Junction is: " << this_junc << endl;

    // check to see if the junction exists
    if (this_junc >= NJunctions)
    {
      cout << "FATAL ERROR LSDJunctionNetwork::calculate_junction_angles" << endl;
      cout << "You have called a junction that doesn't exist." << endl;
      exit(EXIT_FAILURE);
    }

    // check if it is a baselevel node
    int ReceiverJN = get_Receiver_of_Junction(this_junc);
    is_baselevel = (this_junc == ReceiverJN) ? true : false;

    // if not a baselevel see if it has donors
    if (is_baselevel)
    {
      //cout << "This is a baselevel junction." << endl;
      //JunctionAngles.push_back(NoDataValue);
    }
    else
    {
      // check the junction to see if it has two or more donors
      donors = get_donor_nodes(this_junc);

      // it has donors
      // It MUST have 2 donors. Junctions with 3  or more donors are rejected.
      if( int(donors.size()) == 2)
      {

        // reset the temp vectors
        temp_float_junctioninfo = nine_element_float;
        temp_int_junctioninfo = four_element_int;

        //cout << "The donor junctions are: " << donors[0] << ", " << donors[1] << endl;
        // now get the two segments.
        // The ending node is the current junction, the starting nodes are the
        // two donor junctions.
        end_node = get_Node_of_Junction(this_junc);
        start_node1 = get_Node_of_Junction(donors[0]);
        start_node2 = get_Node_of_Junction(donors[1]);
        receiver_node = get_Node_of_Junction(ReceiverJN);

        // extract the channel information
        LSDIndexChannel c1(start_node1, end_node,FlowInfo);
        LSDIndexChannel c2(start_node2, end_node,FlowInfo);
        LSDIndexChannel rchan(end_node, receiver_node,FlowInfo);

        // now get the locations of the nodes in the channels in x,y coordinates
        c1.get_coordinates_of_channel_nodes(x1, y1);
        c2.get_coordinates_of_channel_nodes(x2, y2);
        rchan.get_coordinates_of_channel_nodes(xr, yr);

        // now, for the angles, the channel segments need to both be pointing at the sam
        // "end node", which is the junction in question. Meaning that we
        // need to reverse the xr and yr vectors
        reverse(xr.begin(), xr.end());
        reverse(yr.begin(), yr.end());


        // Now make sure the links have at least 3 nodes
        int min_size = 500;
        if (int(x1.size()) < min_size)
        {
          min_size = int(x1.size());
        }
        if (int(x2.size()) < min_size)
        {
          min_size = int(x2.size());
        }
        if (int(xr.size()) < min_size)
        {
          min_size = int(xr.size());
        }

        if (min_size > 8)
        {
          junction_order = get_StreamOrder_of_Junction(FlowInfo,this_junc);
          donor1_order = get_StreamOrder_of_Junction(FlowInfo,donors[0]);
          donor2_order = get_StreamOrder_of_Junction(FlowInfo,donors[1]);
          receiver_order = get_StreamOrder_of_Junction(FlowInfo,ReceiverJN);

          // now calculate the angle
          // Get the R2 from the vectors. This is a bit stupid since it gets recalculated within
          // the angle code and thrown away, but I am trying to get the code written quickly and the
          // junction angle code is fast
          // The vecotr in the first line for each vector just gets discarded
          //cout << "This junction is: "  << this_junc << endl;
          D1_vec =  orthogonal_linear_regression( x1, y1, intercept, gradient, D1_R2);
          D1_vec =  get_directional_vector_coords_from_dataset(x1, y1, channel_points_downstream);
          D1_bearing = clockwise_angle_between_vector_and_north(0, 0, D1_vec[0], D1_vec[1]);

          D2_vec =  orthogonal_linear_regression( x2, y2, intercept, gradient, D2_R2);
          D2_vec =  get_directional_vector_coords_from_dataset(x2, y2, channel_points_downstream);
          D2_bearing = clockwise_angle_between_vector_and_north(0, 0, D2_vec[0], D2_vec[1]);

          R_vec =  orthogonal_linear_regression( xr, yr, intercept, gradient, R_R2);
          R_vec =  get_directional_vector_coords_from_dataset(xr, yr, channel_points_downstream);
          R_bearing = clockwise_angle_between_vector_and_north(0, 0, R_vec[0], R_vec[1]);

          // Now get the angles between junctions
          // We rotate all the bearings so that R is facing north
          D1_rotated_bearing = D1_bearing-R_bearing;
          if (D1_rotated_bearing < 0)
          {
            D1_rotated_bearing = D1_rotated_bearing+(2*M_PI);
          }

          D2_rotated_bearing = D2_bearing-R_bearing;
          if (D2_rotated_bearing < 0)
          {
            D2_rotated_bearing = D2_rotated_bearing+(2*M_PI);
          }

          // now figure out which one is bigger and calcualte the angles based on that
          if (D2_rotated_bearing > D1_rotated_bearing)
          {
            this_angle_D2_R = (2*M_PI)-D2_rotated_bearing;
            this_angle_D1_R = D1_rotated_bearing;
            this_angle_D1_D2 = D2_rotated_bearing-D1_rotated_bearing;
          }
          else
          {
            this_angle_D1_R = (2*M_PI)-D1_rotated_bearing;
            this_angle_D2_R = D2_rotated_bearing;
            this_angle_D1_D2 = D1_rotated_bearing-D2_rotated_bearing;
          }

          junction_sum = this_angle_D1_R+this_angle_D2_R+this_angle_D1_D2;

          // This attempts to remove junctions screwed up by nans in the bearing calculation.
          if (junction_sum < 2*M_PI+0.001  &&  junction_sum > 2*M_PI-0.001)
          {
            // The JunctionVector  in LSDJunctionNetwork is a vector of the
            // node indices indexed by the junction number
            int node_D1 = get_penultimate_node_from_stream_link(donors[0], FlowInfo);
            DA1 = FlowInfo.get_DrainageArea_square_m( node_D1);
            int node_D2 = get_penultimate_node_from_stream_link(donors[1], FlowInfo);
            DA2 = FlowInfo.get_DrainageArea_square_m( node_D2);

            DAR = FlowInfo.get_DrainageArea_square_m( end_node);

            temp_int_junctioninfo[0] = junction_order;
            temp_int_junctioninfo[1] = donor1_order;
            temp_int_junctioninfo[2] = donor2_order;
            temp_int_junctioninfo[3] = receiver_order;

            temp_float_junctioninfo[0] = DA1;
            temp_float_junctioninfo[1] = DA2;
            temp_float_junctioninfo[2] = DAR;
            temp_float_junctioninfo[3] = this_angle_D1_D2;
            temp_float_junctioninfo[4] = this_angle_D1_R;
            temp_float_junctioninfo[5] = this_angle_D2_R;
            temp_float_junctioninfo[6] = D1_R2;
            temp_float_junctioninfo[7] = D2_R2;
            temp_float_junctioninfo[8] = R_R2;
            temp_float_junctioninfo[9] = D1_bearing;
            temp_float_junctioninfo[10] = D2_bearing;
            temp_float_junctioninfo[11] = R_bearing;
            //temp_float_junctioninfo[12] = D2_vec[1];
            //temp_float_junctioninfo[13] = R_vec[0];
            //temp_float_junctioninfo[14] = R_vec[1];

            JA_int_map[this_junc] = temp_int_junctioninfo;
            JA_float_map[this_junc] = temp_float_junctioninfo;
          }
        }

      }
      else
      {
        //cout << "This junction doesn't have 2 donors; it must be a source." << endl;
        //JunctionAngles.push_back(NoDataValue);
      }
    }
  }

  JA_int_info = JA_int_map;
  JA_float_info = JA_float_map;
}






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This gets the junction angles
//
// It is a more complete version of the junction angle code
//
// This overloaded function includes elevation data and flow distance
//
// The vector of floats actually contains a number of elements
// The int map has:
// [0] Stream order junction 1
// [1] Stream order junction 2
// [2] Stream order receiver junction
//
// The float map has
// [0] drainage area junction 1
// [1] drainage area junction 2
// [2] drainage area receiver junction
// [3] junction angle J1-J2
// [4] junction angle J1-R
// [5] junction angle J2-R
// [6] R^2 J1 segment
// [7] R^2 J2 segment
// [8] R^2 R segment
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::calculate_junction_angles_complete(vector<int> JunctionList,
                                                            LSDFlowInfo& FlowInfo, LSDRaster& Elevation,
                                                            LSDRaster& FlowDistance, float vertical_interval,
                                                            map<int , vector<int> >& JA_int_info,
                                                            map<int, vector<float> >& JA_float_info )
{
  cout << "I'm calculating complete junction angles with elevation and flow distance." << endl;

  map<int, vector<float> > JA_float_map;
  map<int, vector<int> > JA_int_map;
  vector<float> temp_float_junctioninfo;
  vector<int> temp_int_junctioninfo;
  vector<float> nine_element_float(33,0);
  vector<int> four_element_int(4,0);

  int NJuncs = int(JunctionList.size());
  if (NJuncs == 0)
  {
    NJuncs = get_NJunctions();
    cout << "You gave me an empty junction list, I am getting angles for all junctions." << endl;
    vector<int> JL;
    for(int i = 0; i<NJunctions; i++)
    {
      JL.push_back(i);
    }
    JunctionList = JL;
  }

  // now go through each junction getting the angles
  //vector<float> JunctionAngles;
  vector<int> donors;
  float this_angle_D1_D2, this_angle_D1_R, this_angle_D2_R;
  bool is_baselevel;
  bool channel_points_downstream = true;  // this is needed to get the correct
                                          // orientation of the channels
  int end_node, start_node1, start_node2;
  int junction_order, donor1_order, donor2_order;
  int receiver_node, receiver_order;


  // vectors for holding the channel locations
  vector<float> x1, x2, y1, y2, xr, yr;


  // Some floats for holding vector information
  vector<float> D1_vec, D2_vec, R_vec;  // These are the vebor elements of the two donors and receivers.
                                        // They have an origin at 0,0, and the end point is at these coordinates.
  float intercept,gradient;
  float D1_R2, D2_R2, R_R2;

  // The hold drainage areas
  float DA1, DA2, DAR;

  float D1_bearing, D2_bearing, R_bearing;
  float D1_rotated_bearing, D2_rotated_bearing;
  float junction_sum;

  int this_junc;

  // Now we have variables for holding elevation and flow distance information
  float FD_J, FD_D1, FD_D2, FD_R;
  float z_J, z_D1, z_D2, z_R;

  // these are for the vertical intervals
  float FD_D1_vi, FD_D2_vi, FD_R_vi;
  float z_J_vi, z_D1_vi, z_D2_vi, z_R_vi;

  int row,col;
  float grad_D1, grad_D2, grad_R, grad_D1_vi, grad_D2_vi, grad_R_vi;

  // loop through junctions
  for(int junc = 0; junc < NJuncs; junc++ )
  {

    this_junc = JunctionList[junc];
    //cout << endl << "==================" << endl << "Junction is: " << this_junc << endl;

    // check to see if the junction exists
    if (this_junc >= NJunctions)
    {
      cout << "FATAL ERROR LSDJunctionNetwork::calculate_junction_angles" << endl;
      cout << "You have called a junction that doesn't exist." << endl;
      exit(EXIT_FAILURE);
    }

    // check if it is a baselevel node
    int ReceiverJN = get_Receiver_of_Junction(this_junc);
    is_baselevel = (this_junc == ReceiverJN) ? true : false;

    // if not a baselevel see if it has donors
    if (is_baselevel)
    {
      //cout << "This is a baselevel junction." << endl;
      //JunctionAngles.push_back(NoDataValue);
    }
    else
    {
      // check the junction to see if it has two or more donors
      donors = get_donor_nodes(this_junc);

      // it has donors
      // It MUST have 2 donors. Junctions with 3  or more donors are rejected.
      if( int(donors.size()) == 2)
      {

        // reset the temp vectors
        temp_float_junctioninfo = nine_element_float;
        temp_int_junctioninfo = four_element_int;

        //cout << "The donor junctions are: " << donors[0] << ", " << donors[1] << endl;
        // now get the two segments.
        // The ending node is the current junction, the starting nodes are the
        // two donor junctions.
        end_node = get_Node_of_Junction(this_junc);
        start_node1 = get_Node_of_Junction(donors[0]);
        start_node2 = get_Node_of_Junction(donors[1]);
        receiver_node = get_Node_of_Junction(ReceiverJN);


        //cout << "Node is: " << end_node << endl;
        // get all the elevations and flow distances
        FlowInfo.retrieve_current_row_and_col(end_node,row,col);
        //cout << "r,c junc: " << row <<","<< col << endl;
        FD_J = FlowDistance.get_data_element(row,col);
        z_J = Elevation.get_data_element(row,col);

        FlowInfo.retrieve_current_row_and_col(start_node1,row,col);
        //cout << "r,c D1: " << row <<","<< col << endl;
        FD_D1 = FlowDistance.get_data_element(row,col);
        z_D1 = Elevation.get_data_element(row,col);

        FlowInfo.retrieve_current_row_and_col(start_node2,row,col);
        //cout << "r,c D2: " << row <<","<< col << endl;
        FD_D2 = FlowDistance.get_data_element(row,col);
        z_D2 = Elevation.get_data_element(row,col);

        FlowInfo.retrieve_current_row_and_col(receiver_node,row,col);
        //cout << "r,c R: " << row <<","<< col << endl;
        FD_R = FlowDistance.get_data_element(row,col);
        z_R = Elevation.get_data_element(row,col);

        //cout << "FD_J: " << FD_J << " and z_J: " << z_J << endl;
        //cout << "FD_D1: " << FD_D1 << " and z_D1: " << z_D1 << endl;
        //cout << "FD_D2: " << FD_D2 << " and z_D2: " << z_D2 << endl;
        //cout << "FD_R: " << FD_R << " and z_R: " << z_R << endl;

        grad_D1 = (z_D1-z_J)/(FD_D1-FD_J);
        grad_D2 = (z_D2-z_J)/(FD_D2-FD_J);
        grad_R = (z_J-z_R)/(FD_J-FD_R);

        //cout << "grad_D1: " << grad_D1 << " grad_D2: " << grad_D2 << " grad_R: " << grad_R << endl;

        // WORKING HERE 07 DEC NEED TO ADD THE FIXED VEERTICAL INTERVAL
        // We now go through each  donor and the reciever to get the fixed drop
        float z_search;

        // First donor 1
        if ( (z_D1-z_J) <= vertical_interval)
        {
          z_D1_vi = z_D1;
          FD_D1_vi = FD_D1;
          grad_D1_vi = grad_D1;
        }
        else
        {
          // This is repeated in case the if statement for the vertical interval is
          // not triggered. Ensures this coponent returns some value.
          z_D1_vi = z_D1;
          FD_D1_vi = FD_D1;
          grad_D1_vi = grad_D1;

          // Logic if the fixed vertical interval is somwehere within the channel segment
          int this_node = start_node1;
          int search_node;

          while(this_node != end_node)
          {
            FlowInfo.retrieve_receiver_information(this_node,search_node, row, col);

            // Now check the elevation of the search node
            z_search = Elevation.get_data_element(row,col);

            if ((z_search-z_J) <= vertical_interval)
            {
              z_D1_vi = z_search;
              FD_D1_vi = FlowDistance.get_data_element(row,col);
              grad_D1_vi = (z_search-z_J)/(FD_D1_vi-FD_J);
              this_node = end_node;
            }
            else
            {
              this_node = search_node;
            }
          }
        }

        // Now donor 2
        if ( (z_D2-z_J) <= vertical_interval)
        {
          z_D2_vi = z_D2;
          FD_D2_vi = FD_D2;
          grad_D2_vi = grad_D2;
        }
        else
        {
          // This is repeated in case the if statement for the vertical interval is
          // not triggered. Ensures this coponent returns some value.
          z_D2_vi = z_D2;
          FD_D2_vi = FD_D2;
          grad_D2_vi = grad_D2;

          // Logic if the fixed vertical interval is somwehere within the channel segment
          int this_node = start_node2;
          int search_node;

          while(this_node != end_node)
          {
            FlowInfo.retrieve_receiver_information(this_node,search_node, row, col);

            // Now check the elevation of the search node
            z_search = Elevation.get_data_element(row,col);

            if ((z_search-z_J) <= vertical_interval)
            {
              z_D2_vi = z_search;
              FD_D2_vi = FlowDistance.get_data_element(row,col);
              grad_D2_vi = (z_search-z_J)/(FD_D2_vi-FD_J);
              this_node = end_node;
            }
            else
            {
              this_node = search_node;
            }
          }
        }

        // Now the receiver
        if ( (z_J-z_R) <= vertical_interval)
        {
          z_R_vi = z_R;
          FD_R_vi = FD_R;
          grad_R_vi = grad_R;
        }
        else
        {
          // This is repeated in case the if statement for the vertical interval is
          // not triggered. Ensures this coponent returns some value.
          z_R_vi = z_R;
          FD_R_vi = FD_R;
          grad_R_vi = grad_R;

          // Logic if the fixed vertical interval is somwehere within the channel segment
          int this_node = end_node;
          int search_node;

          while(this_node != receiver_node)
          {
            FlowInfo.retrieve_receiver_information(this_node,search_node, row, col);

            // Now check the elevation of the search node
            z_search = Elevation.get_data_element(row,col);

            if ((z_J- z_search) <= vertical_interval)
            {
              z_R_vi = z_search;
              FD_R_vi = FlowDistance.get_data_element(row,col);
              grad_R_vi = (z_J-z_search)/(FD_J-FD_R_vi);
              this_node = receiver_node;
            }
            else
            {
              this_node = search_node;
            }
          }
        }


        // extract the channel information
        LSDIndexChannel c1(start_node1, end_node,FlowInfo);
        LSDIndexChannel c2(start_node2, end_node,FlowInfo);
        LSDIndexChannel rchan(end_node, receiver_node,FlowInfo);

        // now get the locations of the nodes in the channels in x,y coordinates
        c1.get_coordinates_of_channel_nodes(x1, y1);
        c2.get_coordinates_of_channel_nodes(x2, y2);
        rchan.get_coordinates_of_channel_nodes(xr, yr);

        // now, for the angles, the channel segments need to both be pointing at the sam
        // "end node", which is the junction in question. Meaning that we
        // need to reverse the xr and yr vectors
        reverse(xr.begin(), xr.end());
        reverse(yr.begin(), yr.end());


        // Now make sure the links have at least 3 nodes
        int min_size = 500;
        if (int(x1.size()) < min_size)
        {
          min_size = int(x1.size());
        }
        if (int(x2.size()) < min_size)
        {
          min_size = int(x2.size());
        }
        if (int(xr.size()) < min_size)
        {
          min_size = int(xr.size());
        }

        if (min_size > 8)
        {
          junction_order = get_StreamOrder_of_Junction(FlowInfo,this_junc);
          donor1_order = get_StreamOrder_of_Junction(FlowInfo,donors[0]);
          donor2_order = get_StreamOrder_of_Junction(FlowInfo,donors[1]);
          receiver_order = get_StreamOrder_of_Junction(FlowInfo,ReceiverJN);

          // now calculate the angle
          // Get the R2 from the vectors. This is a bit stupid since it gets recalculated within
          // the angle code and thrown away, but I am trying to get the code written quickly and the
          // junction angle code is fast
          // The vector in the first line for each vector just gets discarded
          //cout << "This junction is: "  << this_junc << endl;
          D1_vec =  orthogonal_linear_regression( x1, y1, intercept, gradient, D1_R2);
          D1_vec =  get_directional_vector_coords_from_dataset(x1, y1, channel_points_downstream);
          D1_bearing = clockwise_angle_between_vector_and_north(0, 0, D1_vec[0], D1_vec[1]);

          D2_vec =  orthogonal_linear_regression( x2, y2, intercept, gradient, D2_R2);
          D2_vec =  get_directional_vector_coords_from_dataset(x2, y2, channel_points_downstream);
          D2_bearing = clockwise_angle_between_vector_and_north(0, 0, D2_vec[0], D2_vec[1]);

          R_vec =  orthogonal_linear_regression( xr, yr, intercept, gradient, R_R2);
          R_vec =  get_directional_vector_coords_from_dataset(xr, yr, channel_points_downstream);
          R_bearing = clockwise_angle_between_vector_and_north(0, 0, R_vec[0], R_vec[1]);

          // Now get the angles between junctions
          // We rotate all the bearings so that R is facing north
          D1_rotated_bearing = D1_bearing-R_bearing;
          if (D1_rotated_bearing < 0)
          {
            D1_rotated_bearing = D1_rotated_bearing+(2*M_PI);
          }

          D2_rotated_bearing = D2_bearing-R_bearing;
          if (D2_rotated_bearing < 0)
          {
            D2_rotated_bearing = D2_rotated_bearing+(2*M_PI);
          }

          // now figure out which one is bigger and calcualte the angles based on that
          if (D2_rotated_bearing > D1_rotated_bearing)
          {
            this_angle_D2_R = (2*M_PI)-D2_rotated_bearing;
            this_angle_D1_R = D1_rotated_bearing;
            this_angle_D1_D2 = D2_rotated_bearing-D1_rotated_bearing;
          }
          else
          {
            this_angle_D1_R = (2*M_PI)-D1_rotated_bearing;
            this_angle_D2_R = D2_rotated_bearing;
            this_angle_D1_D2 = D1_rotated_bearing-D2_rotated_bearing;
          }

          junction_sum = this_angle_D1_R+this_angle_D2_R+this_angle_D1_D2;

          // This attempts to remove junctions screwed up by nans in the bearing calculation.
          if (junction_sum < 2*M_PI+0.001  &&  junction_sum > 2*M_PI-0.001)
          {
            // The JunctionVector  in LSDJunctionNetwork is a vector of the
            // node indices indexed by the junction number
            int node_D1 = get_penultimate_node_from_stream_link(donors[0], FlowInfo);
            DA1 = FlowInfo.get_DrainageArea_square_m( node_D1);
            int node_D2 = get_penultimate_node_from_stream_link(donors[1], FlowInfo);
            DA2 = FlowInfo.get_DrainageArea_square_m( node_D2);

            DAR = FlowInfo.get_DrainageArea_square_m( end_node);

            temp_int_junctioninfo[0] = junction_order;
            temp_int_junctioninfo[1] = donor1_order;
            temp_int_junctioninfo[2] = donor2_order;
            temp_int_junctioninfo[3] = receiver_order;

            temp_float_junctioninfo[0] = DA1;
            temp_float_junctioninfo[1] = DA2;
            temp_float_junctioninfo[2] = DAR;
            temp_float_junctioninfo[3] = this_angle_D1_D2;
            temp_float_junctioninfo[4] = this_angle_D1_R;
            temp_float_junctioninfo[5] = this_angle_D2_R;
            temp_float_junctioninfo[6] = D1_R2;
            temp_float_junctioninfo[7] = D2_R2;
            temp_float_junctioninfo[8] = R_R2;
            temp_float_junctioninfo[9] = D1_bearing;
            temp_float_junctioninfo[10] = D2_bearing;
            temp_float_junctioninfo[11] = R_bearing;
            temp_float_junctioninfo[12] = FD_J;
            temp_float_junctioninfo[13] = z_J;
            temp_float_junctioninfo[14] = FD_D1;
            temp_float_junctioninfo[15] = z_D1;
            temp_float_junctioninfo[16] = FD_D2;
            temp_float_junctioninfo[17] = z_D2;
            temp_float_junctioninfo[18] = FD_R;
            temp_float_junctioninfo[19] = z_R;
            temp_float_junctioninfo[20] = grad_D1;
            temp_float_junctioninfo[21] = grad_D2;
            temp_float_junctioninfo[22] = grad_R;
            temp_float_junctioninfo[23] = FD_D1_vi;
            temp_float_junctioninfo[24] = z_D1_vi;
            temp_float_junctioninfo[25] = FD_D2_vi;
            temp_float_junctioninfo[26] = z_D2_vi;
            temp_float_junctioninfo[27] = FD_R_vi;
            temp_float_junctioninfo[28] = z_R_vi;
            temp_float_junctioninfo[29] = grad_D1_vi;
            temp_float_junctioninfo[30] = grad_D2_vi;
            temp_float_junctioninfo[31] = grad_R_vi;

            JA_int_map[this_junc] = temp_int_junctioninfo;
            JA_float_map[this_junc] = temp_float_junctioninfo;
          }
        }

      }
      else
      {
        //cout << "This junction doesn't have 2 donors; it must be a source." << endl;
        //JunctionAngles.push_back(NoDataValue);
      }
    }
  }

  JA_int_info = JA_int_map;
  JA_float_info = JA_float_map;
}








//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function gets the mean and standard error of every junction angle
// upslope of a given junction
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDJunctionNetwork::calculate_junction_angle_statistics_upstream_of_junction(int target_junction, LSDFlowInfo& FlowInfo)
{
  // get all the upslope junctions
  vector<int> JunctionList = get_upslope_junctions(target_junction);
  cout << "There are " << JunctionList.size() << " upslope junctions that I'll analyse" << endl;
  vector<float> JI_stats;

  // now get all the
  map<int, vector<float> >::iterator iter;
  map<int, vector<float> > JuncInfo = calculate_junction_angles(JunctionList,FlowInfo);

  // now get statistics from these
  vector<float> junc_angles;
  vector<float> this_JI;

  for(iter = JuncInfo.begin(); iter != JuncInfo.end(); ++iter)
  {
    this_JI = iter->second;
    if (isnan(this_JI[0]) == false)
    {
      junc_angles.push_back(this_JI[0]);
    }
  }
  cout << "N Junction angles: " << junc_angles.size() << endl;

  // now get the stats
  float mean = get_mean_ignore_ndv(junc_angles,NoDataValue);
  float stddev = get_standard_deviation(junc_angles,mean,NoDataValue);
  float stderr_ =  get_standard_error(junc_angles,stddev);
  // sort the data so we can get the median and percentiles
  vector<size_t> index_map;
  vector<float> junc_angles_sorted;
  matlab_float_sort(junc_angles, junc_angles_sorted, index_map);
  float median = get_median_sorted(junc_angles_sorted);
  float p25 = get_percentile(junc_angles_sorted, 25);
  float p75 = get_percentile(junc_angles_sorted, 75);
  float mad = get_median_absolute_deviation(junc_angles_sorted,median);

  JI_stats.push_back(mean);
  JI_stats.push_back(stderr_);
  JI_stats.push_back(float(junc_angles.size()));
  JI_stats.push_back(median);
  JI_stats.push_back(p25);
  JI_stats.push_back(p75);
  JI_stats.push_back(mad);

  return JI_stats;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints a csv of statistics of the junction angles from a series
// of basins.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::print_junction_angles_from_basin_list(vector<int> JunctionList, LSDFlowInfo& FlowInfo, string csv_outname)
{

  ofstream csv_out;
  csv_out.open(csv_outname.c_str());
  csv_out.precision(9);

  csv_out << "BasinJunction,StreamOrder,Median,25th_percentile,75th_percentile,MAD" << endl;

  // loop through each basin and get the stats
  for (int i = 0; i < int(JunctionList.size()); i++)
  {
    int outlet_jn = JunctionList[i];
    // get all the upslope junctions
    vector<int> upslope_junctions = get_upslope_junctions(outlet_jn);
    cout << "There are " << upslope_junctions.size() << " upslope junctions that I'll analyse" << endl;
    vector<float> JI_stats;

    // now get all the angles
    map<int, vector<float> >::iterator iter;
    map<int, vector<float> > JuncInfo = calculate_junction_angles(upslope_junctions,FlowInfo);

    // now get statistics from these
    vector<float> junc_angles;
    vector<int> stream_order;
    vector<float> this_JI;

    for(iter = JuncInfo.begin(); iter != JuncInfo.end(); ++iter)
    {
      this_JI = iter->second;
      if (isnan(this_JI[0]) == false)
      {
        junc_angles.push_back(this_JI[0]);
        stream_order.push_back(int(this_JI[1]));
      }
    }
    cout << "N Junction angles: " << junc_angles.size() << endl;

    // now sort the vector by stream order. We want to get individual stats
    // for each order
    vector<size_t> index_map;
    vector<int> stream_order_sorted;
    vector<float> junc_angles_sorted;
    matlab_int_sort(stream_order, stream_order_sorted, index_map);
    matlab_float_reorder(junc_angles, index_map, junc_angles_sorted);

    // get the max stream order for this basin = the SO of the junction you're at
    int max_order = get_StreamOrder_of_Junction(outlet_jn);

    // declare a map to store the results.  This has the format:
    // vector<float> medians, vector<float> 25th_percentiles, vector<float> 75th_percentiles, vector<float> median_absolute_deviation
    // then you can get the stream orders by indexing the individual vectors.  E.g. the median
    // of the 3rd stream order would be junction_angle_stats[0][2]
    //map<vector<float>, vector<float>, vector<float>, vector<float>> junction_angle_stats;

    for (int SO = 2; SO < max_order; SO++)
    {
      cout << "This stream order is: " << SO << endl;
      vector<float> these_angles;
      // find the angles of this stream order
      for (int j = 0; j < int(stream_order_sorted.size()); j++)
      {
        cout << "This SO is: " << stream_order_sorted[j] << endl;
        if (stream_order_sorted[j] == (SO))
        {
          these_angles.push_back(junc_angles_sorted[j]);
        }
      }

      // now get the stats
      matlab_float_sort(these_angles, these_angles, index_map);
      float median = get_median_sorted(these_angles);
      float p25 = get_percentile(these_angles, 25);
      float p75 = get_percentile(these_angles, 75);
      float mad = get_median_absolute_deviation(these_angles,median);

      // junction_angle_stats[0].push_back(median);
      // junction_angle_stats[1].push_back(p25);
      // junction_angle_stats[2].push_back(p75);
      // junction_angle_stats[3].push_back(mad);
      cout << "Got the stats, writing to csv" << endl;

      csv_out << outlet_jn << "," << SO << "," << deg(median) << "," << deg(p25) << "," << deg(p75) << "," << deg(mad) << endl;
    }
  }
  csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This gets the junction angle stats for all basins of a given size
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::calculate_junction_angle_statistics_for_order(LSDFlowInfo& FlowInfo, int BasinOrder,
                             vector<int>& junction_list,
                             vector<float>& junction_angle_averages,
                             vector<float>& junction_angle_stderr,
                             vector<int>& N_junctions)
{
  vector<float> jaavg;
  vector<float> jastderr;
  vector<int> N_j;

  // get all the junctions of a given order
  vector<int> OutletJunctions = ExtractBasinJunctionOrder(BasinOrder, FlowInfo);

  // the required information
  pair<vector<float>,vector<float> > JA_stats;

  // now loop through these, getting the statistics of the upstream junctions.
  vector<float> JA_info;
  int n_outlets = int(OutletJunctions.size());
  for(int i = 0; i<n_outlets; i++)
  {
    JA_info = calculate_junction_angle_statistics_upstream_of_junction(OutletJunctions[i], FlowInfo);
    jaavg.push_back(JA_info[0]);
    jastderr.push_back(JA_info[1]);
    N_j.push_back( int(JA_info[2]) );
  }

  // replace the data vectors
  junction_list = OutletJunctions;
  junction_angle_averages = jaavg;
  junction_angle_stderr = jastderr;
  N_junctions = N_j;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints junction angles to a csv file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::print_junction_angles_to_csv(vector<int> JunctionList,
                                                       LSDFlowInfo& FlowInfo,
                                                       string csv_name)
{
  ofstream csv_out;
  csv_out.open(csv_name.c_str());
  csv_out.precision(9);


  // get the junction information
  map<int, vector<float> > JuncInfo = calculate_junction_angles(JunctionList, FlowInfo);
  map<int, vector<float> >::iterator iter;
  vector<float> this_JI;
  int this_junc;
  int this_node,curr_row,curr_col;
  int jso, d1so, d2so;
  double latitude, longitude;

  csv_out << "latitude,longitude,junction_number,junction_stream_order,donor1_stream_order,donor2_stream_order,junction_angle" << endl;
  for(iter = JuncInfo.begin(); iter != JuncInfo.end(); ++iter)
  {
    this_junc = iter->first;
    this_JI = iter->second;

    jso = int(this_JI[1]);
    d1so = int(this_JI[2]);
    d2so = int(this_JI[3]);

    // get the row and column of the junction from the junction node
    this_node = JunctionVector[this_junc];
    LSDCoordinateConverterLLandUTM Converter;
    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);
    get_lat_and_long_locations(curr_row, curr_col, latitude,longitude, Converter);

    // print to the csv file
    csv_out << latitude <<"," << longitude <<"," << this_junc <<"," << jso << ","
            << d1so << "," << d2so << "," << deg(this_JI[0]) << endl;
  }

  csv_out.close();
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints junction angles to a csv file
// Uses the complete junction angle code so much more extensive statistics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::print_complete_junction_angles_to_csv(vector<int> JunctionList,
                                                       LSDFlowInfo& FlowInfo,
                                                       string csv_name)
{
  ofstream csv_out;
  csv_out.open(csv_name.c_str());
  csv_out.precision(9);


  // get the junction information
  map<int, vector<int> > JuncInfo_int;
  map<int, vector<float> > JuncInfo_float;
  calculate_junction_angles_complete(JunctionList, FlowInfo, JuncInfo_int,JuncInfo_float);
  map<int, vector<float> >::iterator iter;
  vector<float> this_JI_float;
  vector<int> this_JI_int;
  int this_junc;
  int this_node,curr_row,curr_col;
  double latitude, longitude;

  csv_out << "latitude,longitude,junction_number,";
  csv_out << "junction_stream_order,donor1_stream_order,donor2_stream_order,receiver_stream_order,";
  csv_out << "donor1_drainage_area,donor2_draiange_area,this_junction_drainage_area,";
  csv_out << "donors_junction_angle,donor1_receiver_junction_angle,donor2_receiver_junction_angle,";
  csv_out << "donor1_R2_on_channel,donor2_R2_on_channel,receiver_R2_on_channel,";
  csv_out << "D1_bearing,D2_bearing,R_bearing" << endl;
  for(iter = JuncInfo_float.begin(); iter != JuncInfo_float.end(); ++iter)
  {
    this_junc = iter->first;
    this_JI_float = iter->second;
    this_JI_int = JuncInfo_int[this_junc];

    // get the row and column of the junction from the junction node
    this_node = JunctionVector[this_junc];
    LSDCoordinateConverterLLandUTM Converter;
    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);
    get_lat_and_long_locations(curr_row, curr_col, latitude,longitude, Converter);

    // print to the csv file
    csv_out << latitude <<"," << longitude <<"," << this_junc <<","
            << this_JI_int[0] << "," << this_JI_int[1] << "," << this_JI_int[2] << "," << this_JI_int[3] << ","
            << this_JI_float[0] << "," << this_JI_float[1]  << "," << this_JI_float[2] << ","
            << deg(this_JI_float[3]) << "," << deg(this_JI_float[4])  << "," << deg(this_JI_float[5]) << ","
            << this_JI_float[6] << "," << this_JI_float[7]  << "," << this_JI_float[8] << ","
            << deg(this_JI_float[9]) << "," << deg(this_JI_float[10])  << "," << deg(this_JI_float[11]) <<endl;
  }

  csv_out.close();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints junction angles to a csv file
// Uses the complete junction angle code so much more extensive statistics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::print_complete_junction_angles_to_csv(vector<int> JunctionList,
                                                       LSDFlowInfo& FlowInfo, LSDRaster& Elevations,
                                                       LSDRaster& FlowDistance, float vertical_interval,
                                                       string csv_name)
{
  cout << "Let me fun the full junction angle code that includes elevations and flow distances." << endl;

  ofstream csv_out;
  csv_out.open(csv_name.c_str());
  csv_out.precision(9);


  // get the junction information
  map<int, vector<int> > JuncInfo_int;
  map<int, vector<float> > JuncInfo_float;

  cout << "Now I will calculate the complete junction angle information." << endl;
  calculate_junction_angles_complete(JunctionList, FlowInfo, Elevations,
                                     FlowDistance, vertical_interval,
                                     JuncInfo_int,JuncInfo_float);
  map<int, vector<float> >::iterator iter;
  vector<float> this_JI_float;
  vector<int> this_JI_int;
  int this_junc;
  int this_node,curr_row,curr_col;
  double latitude, longitude;

  csv_out << "latitude,longitude,junction_number,";
  csv_out << "junction_stream_order,donor1_stream_order,donor2_stream_order,receiver_stream_order,";
  csv_out << "donor1_drainage_area,donor2_drainage_area,this_junction_drainage_area,";
  csv_out << "donors_junction_angle,donor1_receiver_junction_angle,donor2_receiver_junction_angle,";
  csv_out << "donor1_R2_on_channel,donor2_R2_on_channel,receiver_R2_on_channel,";
  csv_out << "D1_bearing,D2_bearing,R_bearing,";
  csv_out << "FlowDistance_junction,Elevation_junction,FlowDistance_donor1,Elevation_donor1,";
  csv_out << "FlowDistance_donor2,Elevation_donor2,FlowDistance_receiver,Elevation_receiver,";
  csv_out << "gradient_donor1,gradient_donor2,gradient_receiver,";
  csv_out << "FlowDistance_d1_vertical_interval,Elevation_d1_vertical_interval,";
  csv_out << "FlowDistance_d2_vertical_interval,Elevation_d2_vertical_interval,";
  csv_out << "FlowDistance_r_vertical_interval,Elevation_r_vertical_interval,";
  csv_out << "gradient_d1_vertical_interval,gradient_d2_vertical_interval,gradient_r_vertical_interval" << endl;

  cout << "csv Header printed, moving on to looping through the junctions" << endl;
  for(iter = JuncInfo_float.begin(); iter != JuncInfo_float.end(); ++iter)
  {
    this_junc = iter->first;
    this_JI_float = iter->second;
    this_JI_int = JuncInfo_int[this_junc];

    // get the row and column of the junction from the junction node
    this_node = JunctionVector[this_junc];
    LSDCoordinateConverterLLandUTM Converter;
    FlowInfo.retrieve_current_row_and_col(this_node,curr_row,curr_col);
    get_lat_and_long_locations(curr_row, curr_col, latitude,longitude, Converter);

    // print to the csv file
    csv_out << latitude <<"," << longitude <<"," << this_junc <<","
            << this_JI_int[0] << "," << this_JI_int[1] << "," << this_JI_int[2] << "," << this_JI_int[3] << ","
            << this_JI_float[0] << "," << this_JI_float[1]  << "," << this_JI_float[2] << ","
            << deg(this_JI_float[3]) << "," << deg(this_JI_float[4])  << "," << deg(this_JI_float[5]) << ","
            << this_JI_float[6] << "," << this_JI_float[7]  << "," << this_JI_float[8] << ","
            << deg(this_JI_float[9]) << "," << deg(this_JI_float[10])  << "," << deg(this_JI_float[11]) << ","
            << this_JI_float[12] << "," << this_JI_float[13]  << "," << this_JI_float[14] << ","
            << this_JI_float[15] << "," << this_JI_float[16]  << "," << this_JI_float[17] << ","
            << this_JI_float[18] << "," << this_JI_float[19]  << "," << this_JI_float[20] << ","
            << this_JI_float[21] << "," << this_JI_float[22]  << ","
            << this_JI_float[23] << "," << this_JI_float[24]  << "," << this_JI_float[25] << ","
            << this_JI_float[26] << "," << this_JI_float[27]  << "," << this_JI_float[28] << ","
            << this_JI_float[29] << "," << this_JI_float[30]  << "," << this_JI_float[31] << endl;
  }

  csv_out.close();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function agregates junction level statistics for larger basins
// Is primarily intended for use with junction angle code but could
// also be used to get statistics about the different slopes or other
// junction level statistics
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
/*
int LSDJunctionNetwork::calculate_basin_averages_from_junction_list(vector<int>& BasinJunctions, LSDFlowInfo& FlowInfo)
{

  // loop through the basins getting the basin statistics.
  int n_basin_junctions = int(BasinJunctions.size());
  for (int jn = 0; jn<n_basin_junctions; j++)
  {

  }


}
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function gets the junction index of a node.  If there is not a junction at this node
// then it will give NoDataValue
//
// FC 31/10/13
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::get_Junction_of_Node(int Node, LSDFlowInfo& FlowInfo)
{
  int JunctionNumber, Row, Col;

  FlowInfo.retrieve_current_row_and_col(Node, Row, Col);
  JunctionNumber = JunctionIndexArray[Row][Col];

  return JunctionNumber;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function returns a vector giving the junction number of all the sources
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_Junctions_of_Sources(LSDFlowInfo& FlowInfo)
{
  vector<int> Sources_junctions;
  int this_junc;

  int n_sources = SourcesVector.size();
  for(int i = 0; i< n_sources; i++)
  {
    this_junc = get_Junction_of_Node(SourcesVector[i], FlowInfo);
    Sources_junctions.push_back(this_junc);
  }
  return Sources_junctions;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function returns the penultimate node of the stream link below a given
// junction
// DTM 04/06/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::get_penultimate_node_from_stream_link(int upstream_junction, LSDFlowInfo& FlowInfo)
{
  int current_junc = upstream_junction;
  int receiver_junc = ReceiverVector[current_junc];
  // First, get channel pixels draining from the current junction.
  LSDIndexChannel StreamLinkVector = LSDIndexChannel(current_junc, JunctionVector[current_junc],
                                                     receiver_junc, JunctionVector[receiver_junc], FlowInfo);
  // Find final nth order channel pixel, which is the penultimate pixel
  // in channel.
  int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  // penultimate node is given by the second to last node in the stream link.
  int penultimate_node;
  if(n_nodes_in_channel > 1) penultimate_node = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  else  penultimate_node = StreamLinkVector.get_node_in_channel(0);
  return penultimate_node;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function generates and LSDIndexChannel given a starting junction
// NOTE: Each junction is the UPSTREAM end of a channel
// this is because junctions can have one and only once receiver
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexChannel LSDJunctionNetwork::generate_link_index_channel_from_junction(int start_junction,
                                       LSDFlowInfo& FlowInfo)
{
  // check bounds on junction
  if(start_junction < 0 || start_junction > NJunctions-1)
  {
    cout << "Tried LSDJunctionNetwork::generate_link_index_channel_from_junction"
         << "  but the junction number does not exist" << endl;
    exit(0);
  }

  // get the starting and ending junctions and nodes
  int end_junction = ReceiverVector[start_junction];
  int start_node = JunctionVector[start_junction];
  int end_node = JunctionVector[end_junction];;

  // extract the channel segment
  LSDIndexChannel this_channel(start_junction, start_node,
                               end_junction, end_node, FlowInfo);

  return this_channel;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function extracts the longest channel originating from a junction number
// outlet_junction
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexChannel LSDJunctionNetwork::generate_longest_index_channel_from_junction(int outlet_junction, LSDFlowInfo& FInfo,
                                        LSDRaster& dist_from_outlet)
{

  if (outlet_junction >= int(JunctionVector.size()))
  {
    cout << "LSDJunctionNetwork::generate_longest_index_channel_from_junction junction not in list" << endl;
    exit(EXIT_FAILURE);
  }

  // first get the number of nodes within the junction
  vector<int> us_junctions = get_upslope_junctions(outlet_junction);

  int n_us_junctions = int(us_junctions.size());
  int farthest_junc = outlet_junction;
  float farthest_dist = 0;
  float current_dist;
  int current_junc, current_row, current_col, current_node;

  // loop through these junctions, looking for the one that is farthest from the outlet
  for (int j = 0; j<n_us_junctions; j++)
  {
    current_junc = us_junctions[j];
    current_node = JunctionVector[current_junc];
    current_row = FInfo.RowIndex[current_node];
    current_col = FInfo.ColIndex[current_node];
    current_dist = dist_from_outlet.get_data_element(current_row,current_col);
    if(current_dist > farthest_dist)
    {
      farthest_dist = current_dist;
      farthest_junc = current_junc;
    }
  }

  int start_junction_node = JunctionVector[farthest_junc];
  int end_junction_node = JunctionVector[outlet_junction];
  LSDIndexChannel Longest_channel(farthest_junc, start_junction_node,
                  outlet_junction, end_junction_node, FInfo);

  return Longest_channel;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function extracts the longest channel originating from a basin
// this differs from the extract_longest_channel_from_junction in that
// basins continue down to a reciever junction
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexChannel LSDJunctionNetwork::generate_longest_index_channel_in_basin(int basin_junction, LSDFlowInfo& FInfo,
                           LSDRaster& dist_from_outlet)
{
  if (basin_junction >= int(JunctionVector.size()))
  {
    cout << "LSDJunctionNetwork::generate_longest_index_channel_in_basin junction not in list" << endl;
    exit(EXIT_FAILURE);
  }

  // first get the number of junctions upslope of the junction
  vector<int> us_junctions = get_upslope_junctions(basin_junction);

  int n_us_junctions = int(us_junctions.size());
  int farthest_junc = basin_junction;
  float farthest_dist = 0;
  float current_dist;
  int current_junc, current_row, current_col, current_node;

  // loop through these junctions, looking for the one that is farthest from the outlet
  for (int j = 0; j<n_us_junctions; j++)
  {
    current_junc = us_junctions[j];
    current_node = JunctionVector[current_junc];
    current_row = FInfo.RowIndex[current_node];
    current_col = FInfo.ColIndex[current_node];
    current_dist = dist_from_outlet.get_data_element(current_row,current_col);
    if(current_dist > farthest_dist)
    {
      farthest_dist = current_dist;
      farthest_junc = current_junc;
    }
  }

  int start_junction_node = JunctionVector[farthest_junc];

  int outlet_junction = ReceiverVector[basin_junction];
  int end_junction_node = JunctionVector[outlet_junction];
  LSDIndexChannel Longest_channel(farthest_junc, start_junction_node,
                  outlet_junction, end_junction_node, FInfo);

  return Longest_channel;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes a vector of basin junctions and returns a vector of the farthest
// upstream source nodes
//
// FJC 21/03/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::get_basin_sources_from_outlet_vector(vector<int> basin_junctions, LSDFlowInfo& FlowInfo,
                           LSDRaster& dist_from_outlet)
{
  vector<int> basin_sources;

  for (int i = 0; i < int(basin_junctions.size()); i++)
  {
    int basin_junction = basin_junctions[i];
    cout << "This junction is: " << basin_junction << endl;
    // if (basin_junction >= int(JunctionVector.size()))
    // {
    //   cout << "LSDJunctionNetwork::generate_longest_index_channel_in_basin junction not in list" << endl;
    //   exit(EXIT_FAILURE);
    // }
    // first get the number of junctions upslope of the junction
    vector<int> us_junctions = get_upslope_junctions(basin_junction);

    int n_us_junctions = int(us_junctions.size());
    int farthest_junc = basin_junction;
    float farthest_dist = 0;
    float current_dist;
    int current_junc, current_row, current_col, current_node;

    // loop through these junctions, looking for the one that is farthest from the outlet
    for (int j = 0; j<n_us_junctions; j++)
    {
      current_junc = us_junctions[j];
      current_node = JunctionVector[current_junc];
      current_row = FlowInfo.RowIndex[current_node];
      current_col = FlowInfo.ColIndex[current_node];
      current_dist = dist_from_outlet.get_data_element(current_row,current_col);
      if(current_dist > farthest_dist)
      {
        farthest_dist = current_dist;
        farthest_junc = current_junc;
      }
    }
    int start_junction_node = JunctionVector[farthest_junc];
    basin_sources.push_back(start_junction_node);
  }

  return basin_sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function gets tributaries along a continous channel.
// What it does is goes down the index channel looking at the JunctionIndexArray
// to see if there is a junction. If it hits a junction then all the contributing junction
//
// it overwrites two vectors:
// tributary_junctions, which lists all junctions whose reciever is the main
// stem
// and
// nodes_on_main_stem_of_tributaries, which are the njodes on the main_stem LSDIndexChannel
// where the tributaries intersect the main stem
// this second vector is used to calcualte the chi values of the downstream node
// of the tributaries
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::extract_tributary_junctions_to_main_stem(LSDIndexChannel& MainStem, LSDFlowInfo& FlowInfo,
                              vector<int>& tributary_junctions,
                              vector<int>& nodes_on_main_stem_of_tributaries)
{
  int n_channel_nodes = MainStem.get_n_nodes_in_channel();
  int node, row, col;
  int curr_junc;
  int previous_junc;
  int n_donors;
  int donor_junc;

  vector<int> temp_tributary_junctions;
  vector<int> temp_nodes_on_main_stem_of_tributaries;

  // the channel starts at a junction, get this junction and set it to
  // be the previous_junction
  MainStem.get_node_row_col_in_channel(0, node, row, col);
  previous_junc = JunctionIndexArray[row][col];

  // now loop through channel, starting at the top
  for(int ch_node = 1; ch_node<n_channel_nodes; ch_node++)
  {
    // get the node index as well as the row and column of the current node in the channel
    MainStem.get_node_row_col_in_channel(ch_node, node, row, col);

    curr_junc = JunctionIndexArray[row][col];
    // if the current junction does not equal the no data value, look for
    // donor nodes
    if(curr_junc != NoDataValue)
    {
      // use the delta vector to get the donor nodes
      n_donors = NDonorsVector[curr_junc];
      for(int donor = 0; donor<n_donors; donor++)
      {
        donor_junc = DonorStackVector[ DeltaVector[curr_junc]+donor ];
        if (donor_junc!=previous_junc)
        {
          temp_tributary_junctions.push_back( donor_junc );
          temp_nodes_on_main_stem_of_tributaries.push_back( ch_node );
        }
      }

      // now set the current junction to the previous junction in order to adnvance further down
      // the channel
      previous_junc = curr_junc;
    }
  }

  tributary_junctions = temp_tributary_junctions;
  nodes_on_main_stem_of_tributaries = temp_nodes_on_main_stem_of_tributaries;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// vector<int> LSDJunctionNetwork::get_pruned_tributaries_from_main_stem
// this function extracts tributaries juncions to the main stem of the channel, then selects
// a sample based on various criteria set by an integer called pruning switch
// pruning_switch == 0  channels are only added if they exceed a threshold drainage area
// pruning_switch == 1  channels are only added if the ratio between them and the mainstem
//						exceeds a certain value (pruning_threshold)
// pruning_switch == 2	channels are only added if the ratio between them and the area of the
//						mainstem _at the junction_ exceeds a certain value
// pruning_switch == 3 channels are only added if the channel order is >= threshold
// DTM 30/04/2013
//--------------------------------------------------------------------------------------------
vector<int> LSDJunctionNetwork::get_pruned_tributaries_from_main_stem(LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChannelNetwork,
int starting_junction, LSDRaster& DistanceFromOutlet, int pruning_switch, float pruning_threshold)
{
  NRows = FlowInfo.get_NRows();
  NCols = FlowInfo.get_NCols();
  DataResolution = FlowInfo.get_DataResolution();
  XMinimum = FlowInfo.get_XMinimum();
  YMinimum = FlowInfo.get_YMinimum();
  NoDataValue = ChannelNetwork.get_NoDataValue();
  float pixel_area = DataResolution*DataResolution;

  //int nodes_in_channel;

  // initiate the main stem
  LSDIndexChannel main_stem = ChannelNetwork.generate_longest_index_channel_in_basin(starting_junction,
                                           FlowInfo, DistanceFromOutlet);
  //nodes_in_channel = main_stem.get_n_nodes_in_channel();

  // Get drainage area of main stem - for pruning tributaries later on
  float main_stem_drainage_area = float(main_stem.get_contributing_pixels_at_outlet(FlowInfo))*pixel_area;
  float contributing_pixel_area;
  float main_stem_junc_area;
  int tributary_order;
  // now get the tributaries for the main stem
  vector<int> tributary_junctions;
  vector<int> node_of_tributaries;
  // Get the tributary junctions
  ChannelNetwork.extract_tributary_junctions_to_main_stem(main_stem, FlowInfo, tributary_junctions, node_of_tributaries);

  // get the number of main stem tributaries
  int n_main_stem_tributaries = tributary_junctions.size();

  // now loop through the tributaries, pruning them so that only those large enough
  // or of the correct order are selected
  vector<int> target_tributary_junctions;
  for (int trib = 0; trib<n_main_stem_tributaries; trib++)
  {
    // Create LSDIndexChannel object for this tributary
    LSDIndexChannel main_stem_tributary = ChannelNetwork.generate_longest_index_channel_in_basin(tributary_junctions[trib],
    FlowInfo, DistanceFromOutlet);

    // get contributing area for tributary to junction with main stem
    contributing_pixel_area = pixel_area*float(main_stem_tributary.get_contributing_pixels_at_penultimate_node(FlowInfo));
    // get contributing area for main stem at the junction with the tributary
    main_stem_junc_area = float(main_stem.get_contributing_pixels_at_node(node_of_tributaries[trib], FlowInfo))*pixel_area;
    // get order of tributary channel
    tributary_order = StreamOrderVector[ tributary_junctions[trib] ];

    // PRUNE TRIBUTARY NETWORK TO ONLY CONSIDER TARGET CHANNELS - USING SPECIFIED PRUNING THRESHOLDS
    if (pruning_switch == 0)        // simple contributing area threshold
    {
      if (contributing_pixel_area > pruning_threshold) target_tributary_junctions.push_back(tributary_junctions[trib]);
    }
    else if (pruning_switch == 1)   // ratio of contributing area of tributary to that of main stem
    {
      if (contributing_pixel_area/main_stem_drainage_area >pruning_threshold) target_tributary_junctions.push_back(tributary_junctions[trib]);
    }
    else if (pruning_switch == 2)   // ratio of contributing area of tributary to that of main stem at junction
    {
      if ( contributing_pixel_area/main_stem_junc_area > pruning_threshold) target_tributary_junctions.push_back(tributary_junctions[trib]);
    }
    else if (pruning_switch == 3)   // channel order
    {
      if (tributary_order > pruning_threshold) target_tributary_junctions.push_back(tributary_junctions[trib]);
    }
    else                            // No pruning
    {
      target_tributary_junctions.push_back(tributary_junctions[trib]);
    }
  }
  return target_tributary_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts basins according to their accumulated drainage area.
// Added by DTM 07/05/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::extract_basin_nodes_by_drainage_area(float DrainageAreaThreshold, LSDFlowInfo& FlowInfo)
{
  // declare the vector containeing the ##node## numbers of the basin outlets
  vector<int> outlet_nodes_of_basins;
  float PixelArea = DataResolution*DataResolution;
  // Loop through junction network until you reach channel junction with required drainage area.
  //int BasinID = 0;
  for (int CurrentJunc=0; CurrentJunc<NJunctions; ++CurrentJunc)
  {
    // Current Junction (should be <= threshold area)
    int CurrentNode = get_Node_of_Junction(CurrentJunc);
    float CurrentJuncDrainageArea = FlowInfo.NContributingNodes[CurrentNode] * PixelArea;
    // ReceiverJunc (should be > threshold area)
    int ReceiverJunc = ReceiverVector[CurrentJunc];
    int ReceiverJunc_Node = get_Node_of_Junction(ReceiverJunc);
    float ReceiverJuncDrainageArea = FlowInfo.NContributingNodes[ReceiverJunc_Node] * PixelArea;

    // Loop through all stream junctions of the required drainage area.  Need to find the
    // stream link in which the drainage area threshold is crossed, before searching through
    // the nodes in this stream link to find the exact node for which to extract the basin
    if ((CurrentJuncDrainageArea <= DrainageAreaThreshold) && (ReceiverJuncDrainageArea > DrainageAreaThreshold))
    {
      cout << "JUNCTION " << CurrentJunc << "/" << NJunctions << " Drainage area = " << CurrentJuncDrainageArea << endl;
      cout << "\t -> found channel segment... upstream drainage area = " << CurrentJuncDrainageArea << " downstream = "
           << ReceiverJuncDrainageArea << "starting search to find the threshold" << endl;

      // First, get a vector containing channel link nodes
      LSDIndexChannel StreamLinkVector = LSDIndexChannel(CurrentJunc, JunctionVector[CurrentJunc],
                                                           ReceiverJunc, JunctionVector[ReceiverJunc], FlowInfo);
      // Find number of nodes in channel
      //int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();

      // Loop through nodes in stream link vector until we reach required drainage area threshold
      int flag = 0;
      //int ChannelSizeTest;
      int NodeCount = 0;
      while ((flag == 0))
      {
        int TargetNode = StreamLinkVector.get_node_in_channel(NodeCount);
        float TargetNodeDrainageArea = FlowInfo.NContributingNodes[TargetNode] * PixelArea;

        // Test for a minimum acceptable drainage area e.g. 80% of desired drainage area before
        // channel is accepted.
        if (TargetNodeDrainageArea >= DrainageAreaThreshold)
        {
          int OutletNode = StreamLinkVector.get_node_in_channel(NodeCount-1);
          float OutletNodeDrainageArea = FlowInfo.NContributingNodes[OutletNode] * PixelArea;

          if (OutletNodeDrainageArea >= (0.8*DrainageAreaThreshold))
          {
            cout << "\t\t\t got it - drainage area = " << FlowInfo.NContributingNodes[OutletNode] * PixelArea << endl;
            outlet_nodes_of_basins.push_back(OutletNode);
            flag = 1;
          }
          else
          {
            cout << "\t\t\t nope - drainage area = " << FlowInfo.NContributingNodes[OutletNode] * PixelArea
                 << " -> too small to be considered" << endl;
            flag = 2;
          }
        }
        ++ NodeCount;
      }
    }
  }
  // Return vector with the relevant junction numbers
  return outlet_nodes_of_basins;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function gets all basin junctions where BOTH upstream catchments are greater
// than a specified drainage area. It will continue to move dowmstream until a base
// level junction is reached, meaning that nested catchments can be selected.
// FJC 10/01/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::extract_basin_nodes_above_drainage_area_threshold(LSDFlowInfo& FlowInfo, float DrainageAreaThreshold)
{
	vector<int> basin_nodes;
	float PixelArea = DataResolution*DataResolution;
	vector<int> visited_before;

	//for each source junction, go downstream until you reach a drainage area > threshold
	//get the source junctions
	vector<int> SourceJunctions = get_Junctions_of_Sources(FlowInfo);
	int n_sources = SourceJunctions.size();
	for (int i = 0; i < n_sources; i++)
	{
		int CurrentJunc = SourceJunctions[i];
		int reached_baselevel = 0;
		int been_to_junction = 0;
		while (reached_baselevel == 0 && been_to_junction == 0)
		{
			//get receiver information
			int ReceiverJunc = get_Receiver_of_Junction(CurrentJunc);
			//logic to see if you've been here before
			vector<int>::iterator find_it;
			find_it = find(visited_before.begin(), visited_before.end(), ReceiverJunc);
			if (find_it == visited_before.end())	// you've never been to this junction before
			{
				// put the junction in the visited before vector
				visited_before.push_back(ReceiverJunc);
				// logic for base level
				if (CurrentJunc == ReceiverJunc)
				{
					reached_baselevel = 1;
					break;
				}
				else
				{
					//get penultimate node
					int OutletNode = get_penultimate_node_from_stream_link(CurrentJunc, FlowInfo);
					float OutletDrainageArea = FlowInfo.NContributingNodes[OutletNode]*PixelArea;
					if (OutletDrainageArea >= DrainageAreaThreshold)
					{
						//this basin has an area > threshold, check the other tributary basin
						vector<int> us_juncs = get_donor_nodes(ReceiverJunc);
						for (int j = 0; j < int(us_juncs.size()); j++)
						{
							//cout << "Donor juncs: " << us_juncs[j] << endl;
							int base_level = is_Junction_BaseLevel(us_juncs[j]);
							if (us_juncs[j] != CurrentJunc && base_level == 0)
							{
								// check the drainage area of the other basin
								int ThisOutletNode = get_penultimate_node_from_stream_link(us_juncs[j], FlowInfo);
								float ThisOutletDrainageArea = FlowInfo.NContributingNodes[ThisOutletNode]*PixelArea;
								if (ThisOutletDrainageArea >= DrainageAreaThreshold)
								{
									//both basins have drainage area > threshold, push back the penultimate nodes of both basins
									//cout << "Basin node: " << OutletNode << endl;
									//cout << "Basin node: " << ThisOutletNode << endl;
 									basin_nodes.push_back(OutletNode);
									basin_nodes.push_back(ThisOutletNode);
								}
							}
						}
					}
					CurrentJunc = ReceiverJunc;
					//cout << "Moving downstream" << endl;
				}
			}
			else	//You've visited this junction before, so you can move to the next source
			{
				been_to_junction = 1;
			}
		}
	}
	return basin_nodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function checks all of the basin nodes to check if they fall within a
// mask (input raster). If they fall within the mask raster then the first node
// upstream not in the mask is selected.
//
// FJC 31/01/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::modify_basin_nodes_from_mask(vector<int> basin_nodes, LSDFlowInfo& FlowInfo, LSDRaster& MaskRaster)
{
  vector<int> NewBasinNodes;

  for (int i = 0; i < int(basin_nodes.size()); i++)
  {
    // check mask to see if nodes fall within the mask
    int row, col;
    int this_node = basin_nodes[i];
    //cout << "This basin node: " << this_node << endl;
    FlowInfo.retrieve_current_row_and_col(this_node, row, col);
    float mask_value = MaskRaster.get_data_element(row,col);
    if (mask_value <= 0)
    {
      // this node is not masked, so just return as normal
      NewBasinNodes.push_back(this_node);
    }
    else
    {
      // this node is masked, need to return upstream node
      int reached_upstream = 0;
      while (reached_upstream == 0)
      {
        int UpstreamNode = get_upstream_node_max_stream_order(this_node, FlowInfo);
        if (UpstreamNode != NoDataValue)
        {
          // check the upstream node
          int upstream_row, upstream_col;
          FlowInfo.retrieve_current_row_and_col(UpstreamNode, upstream_row, upstream_col);
          //cout << " Upstream row: " << upstream_row << " Upstream col: " << upstream_col << endl;
          float UpstreamValue = MaskRaster.get_data_element(upstream_row, upstream_col);
          if (UpstreamValue == NoDataValue)
          {
            // node not masked, return the upstream node
            NewBasinNodes.push_back(UpstreamNode);
            cout << "New node: " << UpstreamNode << endl;
            reached_upstream = 1;
          }
          else
          {
            // node masked, move upstream
            this_node = UpstreamNode;
          }
        }
        else { break; }
      }
    }
  }
  return NewBasinNodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts basin junctions from a vector of the basin outlet nodes.
// Added by FJC 15/01/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::extract_basin_junctions_from_nodes(vector<int> basin_nodes, LSDFlowInfo& FlowInfo)
{
  vector<int> basin_junctions;

  for (unsigned int i =0; i < basin_nodes.size(); i++)
  {
    int node = basin_nodes[i];
    int junction = find_upstream_junction_from_channel_nodeindex(node, FlowInfo);
    //cout << junction << endl;
    basin_junctions.push_back(junction);
  }

  return basin_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts the base level junction to which a starting junction
// drains
// Added by SMM 21/02/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::find_base_level_node_of_junction(int StartingJunction)
{
  int ThisJunction = StartingJunction;

  if(ThisJunction >= int(ReceiverVector.size()))
  {
    cout << "Warning, you have selected a junction that doesn't exist!" << endl;
    cout << "Defaulting to Junction 0" << endl;
    ThisJunction = 0;
  }
  while(ReceiverVector[ThisJunction] != ThisJunction)
  {
       ThisJunction =  ReceiverVector[ThisJunction];
  }
  return ThisJunction;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the stream order of a node.
// FJC 29/09/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_StreamOrder_of_Node(LSDFlowInfo& FlowInfo, int node)
{
  int StreamOrder;
  int row,col;

  FlowInfo.retrieve_current_row_and_col(node,row,col);
  StreamOrder = StreamOrderArray[row][col];

  return StreamOrder;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the stream order of a junction.
// FJC 01/03/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_StreamOrder_of_Junction(LSDFlowInfo& FlowInfo, int junction)
{
  NJunctions = int(JunctionVector.size());
  if(junction <0 || junction >= NJunctions)
  {
    cout << "You have selected a junction that is not in the junction list." << endl
         << "defaulting to junction 0!" << endl;
    junction = 0;
  }

  int StreamOrder;
  int row,col;
  int node = get_Node_of_Junction(junction);

  FlowInfo.retrieve_current_row_and_col(node,row,col);
  StreamOrder = StreamOrderArray[row][col];

  return StreamOrder;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the stream order of a junction. It samples directly from
// the stream order vector
// SMM 26/10/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_StreamOrder_of_Junction(int junction)
{
  NJunctions = int(JunctionVector.size());
  if(junction <0 || junction >= NJunctions)
  {
    cout << "You have selected a junction that is not in the junction list." << endl
         << "defaulting to junction 0!" << endl;
    junction = 0;
  }
  int StreamOrder = StreamOrderVector[junction];

  return StreamOrder;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function finds the next junction that is of higher stream order
// SMM 26/10/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_Next_StreamOrder_Junction(int junction)
{
  // check bounds
  NJunctions = int(JunctionVector.size());
  if(junction <0 || junction >= NJunctions)
  {
    cout << "You have selected a junction that is not in the junction list." << endl
         << "defaulting to junction 0!" << endl;
    junction = 0;
  }

  // now follow the receiver junctions until you hit either a baselevel node
  // or a junction with a higher stream order.
  int current_SO = StreamOrderVector[junction];
  int next_junc = get_Receiver_of_Junction(junction);
  int next_SO = StreamOrderVector[next_junc];
  int this_junc = junction;
  //cout << "this junc: " << this_junc << " next junc: " << next_junc << " next_SO: " << next_SO << endl;
  while (next_junc != this_junc && next_SO == current_SO)
  {
    this_junc = next_junc;
    next_junc = get_Receiver_of_Junction(this_junc);
    next_SO = StreamOrderVector[next_junc];
    //cout << "this junc: " << this_junc << " next junc: " << next_junc << " next_SO: " << next_SO << endl;
  }

  //cout << "Starting junction: " << junction << " junction SO: " << current_SO << endl;
  //cout << "Finishing junction: " << next_junc << " finishing SO: " << next_SO << endl;
  if (next_junc == this_junc)
  {
    //cout << "You reached a baselevel node, returning nodata" << endl;
    next_junc = int(NoDataValue);
  }
  return next_junc;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function checks whether a junction is upsream of another junction
// SMM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
bool LSDJunctionNetwork::is_junction_upstream(int current_junction, int test_junction)
{
  bool i = false;

  int start_SVector_junction = SVectorIndex[current_junction];
  int end_SVector_junction = start_SVector_junction+NContributingJunctions[current_junction];

  int SVector_test_junction = SVectorIndex[test_junction];

  for(int junction = start_SVector_junction; junction < end_SVector_junction; junction++)
  {
    if (junction == SVector_test_junction)
    {
      i = true;
    }
  }

  return i;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function checks whether a junction is at base level
// FJC 11/01/2017
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::is_Junction_BaseLevel(int junction)
{
  int base_level = 0;
  int ReceiverJN = get_Receiver_of_Junction(junction);
  if (junction == ReceiverJN) { base_level = 1; }

  return base_level;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// Two getter functions that require bounds checking
// Added by SMM 21/02/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_Node_of_Junction(int junction) const
{
  if(junction >= int(JunctionVector.size()))
  {
    cout << "Warning, you have selected a junction that doesn't exist!" << endl;
    cout << "Defaulting to Junction 0" << endl;
    junction = 0;
  }
  return JunctionVector[junction];
}

int LSDJunctionNetwork::get_Receiver_of_Junction(int junction) const
{
  if(junction >= int(ReceiverVector.size()))
  {
    cout << "Warning, you have selected a junction that doesn't exist!" << endl;
    cout << "Defaulting to Junction 0" << endl;
    junction = 0;
  }
  return ReceiverVector[junction];
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// Converting a junction list into a node list
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::get_node_list_from_junction_list(vector<int> junction_list)
{
  int N_juncs = junction_list.size();
  vector<int> node_list;
  int this_junc;

  for (int i = 0; i< N_juncs; i++)
  {
    this_junc = junction_list[i];
    if (this_junc >= 0 && this_junc < int(JunctionVector.size()))
    {
      node_list.push_back(JunctionVector[this_junc]);
    }
  }

  return node_list;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// Converting a junction list into a node list
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::get_node_list_of_penultimate_node_from_junction_list(vector<int> junction_list, LSDFlowInfo& FlowInfo)
{
  int N_juncs = junction_list.size();
  vector<int> node_list;
  int this_junc;

  for (int i = 0; i< N_juncs; i++)
  {
    this_junc = junction_list[i];
    if (this_junc >= 0 && this_junc < int(JunctionVector.size()))
    {
      int outlet_node = get_penultimate_node_from_stream_link(this_junc,FlowInfo);
      node_list.push_back(outlet_node);
    }
  }

  return node_list;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// Function to get the junction downstream of the next
// Added by FJC 08/10/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_downstream_junction(int starting_junction, LSDFlowInfo& FlowInfo)
{
    int start_node = get_Node_of_Junction(starting_junction);
    int receiver_node, receiver_row, receiver_col, receiver_junction;
    int i = 0;
    while (i == 0)
    {
      FlowInfo.retrieve_receiver_information(start_node, receiver_node, receiver_row, receiver_col);
      receiver_junction = get_Junction_of_Node(receiver_node, FlowInfo);
      if (receiver_junction == NoDataValue)
      {
        start_node = receiver_node;
      }
      else
      {
        //cout << "Reached downstream junction" << endl;
        i=1;
      }
    }

    return receiver_junction;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts the junctions of a given basin order that are the lowermost
// junction in the basin
//
// IMPORTANT: the junctions always point downstream since they can have one and only
// one receiver. However, for a basin of given order, this starts just upstream of the
// confluence to the next basin order. So the basin ##INCLUDES## the channel flowing
// downstream to the penultimite node
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::extract_basins_order_outlet_junctions(int BasinOrder, LSDFlowInfo& FlowInfo)
{
  // declare the vector containeing the ##Junction## numbers of the basin outlets
  vector<int> outlet_junctions_of_basins_of_order_x;

  // Loop through junction network until you reach nth order channel junction.
  int BasinID = 0;
  int current_junc,receiver_junc,receiver_junc_SO;
  for (int junctionID=0; junctionID<NJunctions; ++junctionID)
  {
    // Loop through all stream junctions of the required basin order.
    // Note that the nth order basins are defined as capturing the full
    // drainage area for the nth order stream.  An nth order stream terminates
    // at a junction an order greater than order n.
    if (StreamOrderVector[junctionID] == BasinOrder)
    {
      // Get info from ChanelNetwork object regarding position of junction
      current_junc = junctionID;//JunctionVector[junctionID];
      receiver_junc = ReceiverVector[current_junc];
      receiver_junc_SO = StreamOrderVector[receiver_junc];
      // Identify outlet of nth order basin using the condition that the
      // receiver junction should be of higher order.
      if (receiver_junc_SO > BasinOrder)
      {
        bool IsTruncated = node_tester(FlowInfo,current_junc);
        if(IsTruncated == false)
        {
          outlet_junctions_of_basins_of_order_x.push_back(current_junc);
          // Increment BasinID to ensure that each basin is distinct.
          ++BasinID;
        }
      } // end of while logic - have searched through junction catchment
    }
  }
  // Return vector with the relevant junction numbers
  return outlet_junctions_of_basins_of_order_x;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the node indices of outlets of basins of a certain order
//
// IMPORTANT: the junctions always point downstream since they can have one and only
// one receiver. However, for a basin of given order, this starts just upstream of the
// confluence to the next basin order. So the basin ##INCLUDES## the channel flowing
// downstream to the penultimite node
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::extract_basins_order_outlet_nodes(vector<int>& BasinOutletJunctions, LSDFlowInfo& FlowInfo)
{
  // declare the vector containeing the ##Node## numbers of the basin outlets
  vector<int> outlet_nodes_of_basins_of_order_x;

  // find how many basins there are
  int n_basins = BasinOutletJunctions.size();

  int current_junc, receiver_junc, n_nodes_in_channel,basin_outlet;
  for (int basinID=0; basinID<n_basins; ++basinID)
  {
    current_junc = BasinOutletJunctions[basinID];
    receiver_junc = ReceiverVector[current_junc];

    // First, get channel pixels draining from the current junction.
    LSDIndexChannel StreamLinkVector = LSDIndexChannel(current_junc, JunctionVector[current_junc],
                     receiver_junc, JunctionVector[receiver_junc], FlowInfo);

    // Find final nth order channel pixel, which is the penultimate pixel
    // in channel.
    n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
    // ======================
    // gets crazy again here.
    basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
    outlet_nodes_of_basins_of_order_x.push_back(basin_outlet);
  }
  // Return vector with the relevant junction numbers
  return outlet_nodes_of_basins_of_order_x;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts the drainage basins of a given order, n
// DTM 17/10/2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::ExtractBasinsOrder(int BasinOrder, LSDFlowInfo& FlowInfo)
{
  // Declare drainage basins array
  Array2D<int> basins(NRows,NCols,NoDataValue);
  // Loop through junction network until you reach nth order channel junction.
  int BasinID = 0;
  int row,col,node,current_junc,receiver_junc,receiver_junc_SO,n_nodes_in_channel,basin_outlet;
  for (int junctionID=0; junctionID<NJunctions; ++junctionID)
  {
    // Loop through all stream junctions of the required basin order.
    // Note that the nth order basins are defined as capturing the full
    // drainage area for the nth order stream.  An nth order stream terminates
    // at a junction an order greater than order n.
    if (StreamOrderVector[junctionID] == BasinOrder)
    {
      // Get info from ChanelNetwork object regarding position of junction
      current_junc = junctionID;//JunctionVector[junctionID];
      receiver_junc = ReceiverVector[current_junc];
      receiver_junc_SO = StreamOrderVector[receiver_junc];
      // Identify outlet of nth order basin using the condition that the
      // receiver junction should be of higher order.
      if (receiver_junc_SO > BasinOrder)
      {
        bool IsTruncated = node_tester(FlowInfo,current_junc);

        if(IsTruncated == false)
        {
          // First, get channel pixels draining from the current junction.
          LSDIndexChannel StreamLinkVector = LSDIndexChannel(current_junc, JunctionVector[current_junc],
                                                             receiver_junc, JunctionVector[receiver_junc], FlowInfo);
          // Find final nth order channel pixel, which is the penultimate pixel
          // in channel.
          n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
          // ======================
          // gets crazy again here.
          basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
          // Get all contributing pixels and label with BasinID
          vector<int> BasinNodeVector = FlowInfo.get_upslope_nodes(basin_outlet);
          // Loop through basin to label basin pixels with basin ID
          for (int BasinIndex = 0; BasinIndex < int(BasinNodeVector.size()); ++BasinIndex)
          {
            node = BasinNodeVector[BasinIndex];
            FlowInfo.retrieve_current_row_and_col(node,row,col);
            basins[row][col] = BasinID;
          }
          // Increment BasinID to ensure that each basin is distinct.
          ++BasinID;
        }
      } // end of while logic - have searched through junction catchment
    }
  }
  // Return raster with all nth order drainage basins.
  LSDIndexRaster basin_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,basins,GeoReferencingStrings);
  return basin_raster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts the juctions of all non-beheaded drainage basins of a given order, n
// SWDG 23/10/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::ExtractBasinJunctionOrder(int BasinOrder, LSDFlowInfo& FlowInfo)
{
  vector<int> Junctions;
  // Loop through junction network until you reach nth order channel junction.

  int current_junc,receiver_junc,receiver_junc_SO;
  for (int junctionID=0; junctionID<NJunctions; ++junctionID)
  {
    // Loop through all stream junctions of the required basin order.
    // Note that the nth order basins are defined as capturing the full
    // drainage area for the nth order stream.  An nth order stream terminates
    // at a junction an order greater than order n.

    if (StreamOrderVector[junctionID] == BasinOrder)
    {
      //cout << "Found a junction that is of the correct order." << endl;
      //cout << "Junction is: " << junctionID << endl;
      //cout << "This junction order is: " <<  StreamOrderVector[junctionID] << endl;

      // Get info from ChanelNetwork object regarding position of junction
      current_junc = junctionID;//JunctionVector[junctionID];
      receiver_junc = ReceiverVector[current_junc];
      receiver_junc_SO = StreamOrderVector[receiver_junc];
      // Identify outlet of nth order basin using the condition that the
      // receiver junction should be of higher order.
      if (receiver_junc_SO > BasinOrder)
      {
        //use the node tester to get rid of any basins that are beheaded
        if (node_tester(FlowInfo, current_junc) == false)
        {
          Junctions.push_back(current_junc);
        }
      }
    }
  }

  return Junctions;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts the juctions of all non-beheaded drainage basins
// of a given order, n
// It keeps basins that abut nodata values

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::ExtractBasinJunctionOrderKeepEdgeBasins(int BasinOrder, LSDFlowInfo& FlowInfo)
{
  vector<int> Junctions;
  // Loop through junction network until you reach nth order channel junction.

  int current_junc,receiver_junc,receiver_junc_SO;
  for (int junctionID=0; junctionID<NJunctions; ++junctionID)
  {
    // Loop through all stream junctions of the required basin order.
    // Note that the nth order basins are defined as capturing the full
    // drainage area for the nth order stream.  An nth order stream terminates
    // at a junction an order greater than order n.

    if (StreamOrderVector[junctionID] == BasinOrder)
    {
      //cout << "Found a junction that is of the correct order." << endl;
      //cout << "Junction is: " << junctionID << endl;
      //cout << "This junction order is: " <<  StreamOrderVector[junctionID] << endl;

      // Get info from ChanelNetwork object regarding position of junction
      current_junc = junctionID;//JunctionVector[junctionID];
      receiver_junc = ReceiverVector[current_junc];
      receiver_junc_SO = StreamOrderVector[receiver_junc];
      // Identify outlet of nth order basin using the condition that the
      // receiver junction should be of higher order.
      if (receiver_junc_SO > BasinOrder)
      {
        Junctions.push_back(current_junc);
      }
    }
  }

  return Junctions;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function takes a junction number and then cycles through the upslope
// junctions, looking for sources. Once it finds a source it then traces
// the channel to the divide (finding the longest segment) and returns the
// node indices of these hilltop nodes.
//
// SMM 25/9/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
vector<int> LSDJunctionNetwork::FindFarthestUpslopeHilltopsFromSources(int JunctionNumber,
                                          LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance)
{
  // get the list of upslope junctions
  vector<int> upslope_juncs = get_upslope_junctions(JunctionNumber);

  vector<int> hilltops_from_sources;
  int upslope_node_index;
  int this_source_node;

  // loop through these junctions, looking for sources
  int n_upslope_juncs = upslope_juncs.size();
  for(int i = 0; i<n_upslope_juncs; i++)
  {
    // test to see if the junction has donor, if not it is a source
    if (NDonorsVector[ upslope_juncs[i] ] == 0)
    {
      //cout << "source found, junction: " << upslope_juncs[i] << endl;

      // get the node
      this_source_node =  JunctionVector[ upslope_juncs[i] ];

      // now find the upslope node of this souce
      upslope_node_index =  FlowInfo.find_farthest_upslope_node(this_source_node, FlowDistance);

      // add it to the vector
      hilltops_from_sources.push_back(upslope_node_index);
    }
  }

  return hilltops_from_sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns a vector of nodeindex values of potential channel heads
// above a given junction
//
// SMM 26/09/2013
//
// Modified by FC 26/09/2013 to extract the node of the predicted channel heads for each
// hilltop node; now returns a vector of integers.
// Modified by FC 18/11/13 to extract the channel head from the outlet junction rather than from
// the sources.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
int LSDJunctionNetwork::GetChannelHeadsChiMethodFromNode(int NodeNumber,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster)
{
  //vector<int> ChannelHeadNodes;
  float downslope_chi = 0;

  // get the second order junction from this node
  int node_at_junction = FlowInfo.ReceiverVector[NodeNumber];
  int first_order_junction = get_Junction_of_Node (node_at_junction, FlowInfo);
  int second_order_junction = get_Receiver_of_Junction(first_order_junction);
  int second_order_node = get_Node_of_Junction(second_order_junction);

  // get the hilltop node from this junction
  int hilltop_node = FlowInfo.find_farthest_upslope_node(NodeNumber, FlowDistance);

  //perform chi segment fitting
  LSDChannel new_channel(hilltop_node, second_order_node, downslope_chi, m_over_n, A_0, FlowInfo,  ElevationRaster);
  int channel_head_node = new_channel.calculate_channel_heads(MinSegLength, A_0, m_over_n, FlowInfo);

  // get the nodes of the hilltops.
  //vector<int> hilltop_nodes = FindFarthestUpslopeHilltopsFromSources(JunctionNumber,
  //                                        FlowInfo, FlowDistance);

  return channel_head_node;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function returns a vector of nodeindex values of potential channel heads
// above a given junction
//
// SMM 26/09/2013
//
// Modified by FC 26/09/2013 to extract the node of the predicted channel heads for each
// hilltop node; now returns a vector of integers.
// Modified by FC 18/11/13 to extract the channel head from the outlet junction rather than from
// the sources.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::GetChannelHeadsChiMethodFromSourceNode(int NodeNumber,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster, int NJunctions)
{
  int channel_head_node;
  //vector<int> ChannelHeadNodes;
  float downslope_chi = 0;

  // get the hilltop node from this source node
  int hilltop_node = FlowInfo.find_farthest_upslope_node(NodeNumber, FlowDistance);

  //get the junction at the source node
  int source_junction = get_Junction_of_Node(NodeNumber, FlowInfo);

  if (NJunctions == 0)
  {
    LSDChannel new_channel(hilltop_node, NodeNumber, downslope_chi, m_over_n, A_0, FlowInfo,  ElevationRaster);
    channel_head_node = new_channel.calculate_channel_heads(MinSegLength, A_0, m_over_n, FlowInfo);
  }
  else if (NJunctions > 0)
  {
    int count = 0;
    // get the nth junction downstream
    for (int i = 0; i < NJunctions; i++)
    {
      //cout << "Source junction: " << source_junction << endl;
      int downstream_junction = get_downstream_junction(source_junction, FlowInfo);
      //cout << "downstream junction: " << downstream_junction << endl;
      source_junction = downstream_junction;
      count++;
      //cout << "Moved " << count << " junctions downstream from source" << endl;
    }
    int final_node = get_Node_of_Junction(source_junction);
    //cout << "Node of downstream junction: " << final_node << endl;
    //cout << "Start node of channel: " << hilltop_node << endl;
    //perform chi segment fitting
    LSDChannel new_channel(hilltop_node, final_node, downslope_chi, m_over_n, A_0, FlowInfo,  ElevationRaster);
    channel_head_node = new_channel.calculate_channel_heads(MinSegLength, A_0, m_over_n, FlowInfo);
  }
  else
  {
    cout << "Something has gone wrong, you have a negative number of Junctions." << endl;
    channel_head_node = NoDataValue;
  }
  return channel_head_node;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function writes a csv file of the chi and elevation values for each valley to
// hilltop profile in the DEM
// csv file has the suffix of the source junction
// N junctions = number of junctions downstream of the source junction to get the profile for
// FJC 23/12/16
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
void LSDJunctionNetwork::write_valley_hilltop_chi_profiles_to_csv(vector<int> sources, float A_0, float m_over_n, LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance, LSDRaster& ElevationRaster, int NJunctions, string output_path, string DEM_ID)
{
	float downslope_chi = 0;

	//loop through all the sources and get the channel profile from the valley to the hilltop
	for (int i = 0; i < int(sources.size()); i++)
	{
		// get the farthest upslope hilltop node
		int hilltop_node = FlowInfo.find_farthest_upslope_node(sources[i], FlowDistance);
		// get the valley node
		int source_junction = get_Junction_of_Node(sources[i], FlowInfo);
		// move downstream the specified number of junctions
		for (int j = 0; j < NJunctions; j++)
		{
			int downstream_junction = get_downstream_junction(source_junction, FlowInfo);
			source_junction = downstream_junction;
		}
		int final_node = get_Node_of_Junction(source_junction);
		// get the LSDChannel
		LSDChannel new_channel(hilltop_node, final_node, downslope_chi, m_over_n, A_0, FlowInfo, ElevationRaster);
		// write to csv
		string jn_str = static_cast<ostringstream*>( &(ostringstream() << source_junction) )->str();
		string output_csv_filename = DEM_ID+"_chan_profile_"+jn_str;
		new_channel.write_channel_to_csv(output_path, output_csv_filename, FlowDistance);
	}
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns a vector of nodeindex values of potential channel heads
// above a given junction
//
// SMM 26/09/2013
//
// Modified by FC 26/09/2013 to extract the node of the predicted channel heads for each
// hilltop node; now returns a vector of integers.
// Modified by FC 18/11/13 to extract the channel head from the outlet junction rather than from
// the sources.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::GetChannelfromDreich(int NodeNumber,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster, string path_name, int NJunctions)
{
  //vector<int> ChannelHeadNodes;
  float downslope_chi = 0;

  //get the junction at the source node
  int Junction = get_Junction_of_Node(NodeNumber, FlowInfo);
  // get the hilltop node from this junction
  int hilltop_node = FlowInfo.find_farthest_upslope_node(NodeNumber, FlowDistance);

  // get the nth junction downstream
  for (int i = 0; i < NJunctions; i++)
  {
    int downstream_junction = get_Receiver_of_Junction(Junction);
    Junction = downstream_junction;
  }
  int final_node = get_Node_of_Junction(Junction);

  //perform chi segment fitting
  LSDChannel new_channel(hilltop_node, final_node, downslope_chi, m_over_n, A_0, FlowInfo,  ElevationRaster);
  //string node = static_cast<ostringstream*>( &(ostringstream() << NodeNumber) )->str();
  //string file_name = "full_profile_"+node;
  //int channel_head_node = new_channel.calculate_channel_heads_with_profile(MinSegLength, A_0, m_over_n, FlowInfo, node);
  //new_channel.write_channel_to_csv(path_name, file_name, FlowDistance);
  LSDIndexRaster ChannelRaster = new_channel.print_channel_to_IndexRaster(FlowInfo);

  return ChannelRaster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns all potential channel heads in a DEM. It looks for
// channel heads organized by a basin order which is fed to the code
// The basin order just determines how far downstream the algorithm looks for the 'fluvial'
// section.
// It returns a vector<int> of nodeindices where the channel heads are
//
// SMM 26/09/2013
//
//
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// vector<int> LSDJunctionNetwork::GetChannelHeadsChiMethodBasinOrder(int BasinOrder,
//                                       int MinSegLength, float A_0, float m_over_n,
// 									  LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
// 									  LSDRaster& ElevationRaster)
// {
// 	vector<int> ChannelHeadNodes;
// 	vector<int> ChannelHeadNodes_temp;
//
//   	vector<int> junction_list = extract_basins_order_outlet_junctions(BasinOrder, FlowInfo);
//   	int max_junctions = junction_list.size();
//   	cout << "No of junctions: " << max_junctions << endl;
//   	int junction_number = 0;
//
// 	 //loop through junctions collecting channel heads
//     for (int i = 0; i < max_junctions; i++)
//   	{
// 		  cout << "Junction " << i << " of " << max_junctions << endl;
//
// 		  // get the junction number
// 		  junction_number = junction_list[i];
//
// 		  // get a local list of channel heads
// 		  int channel_head_node = GetChannelHeadsChiMethodFromNode(junction_number,
//                                       			MinSegLength, A_0, m_over_n, FlowInfo,
//                                        			FlowDistance, ElevationRaster);
//
//       // now append these channel heads to the master list
// 		  //ChannelHeadNodes_temp.insert(ChannelHeadNodes_temp.end(), these_channel_heads.begin(), these_channel_heads.end());
// 		  ChannelHeadNodes_temp.push_back(channel_head_node);
// 	  }
//
// 	  // Removing any nodes that are not the furthest upstream
//     int upstream_test = 0;
//     vector<int>::iterator find_it;
//
//     for (unsigned int node =0; node < ChannelHeadNodes_temp.size(); node++)
//     {
//       vector<int> tests;
//       int current_node = ChannelHeadNodes_temp[node];
//       for (unsigned int i = 0; i < ChannelHeadNodes_temp.size(); i++)
//       {
//         if (ChannelHeadNodes_temp[i] != current_node)
//         {
//           int test_node = ChannelHeadNodes_temp[i];
//           upstream_test = FlowInfo.is_node_upstream(current_node, test_node);
//           tests.push_back(upstream_test);
//         }
//       }
//       find_it = find(tests.begin(), tests.end(), 1);
//       if (find_it == tests.end())
//       {
//         ChannelHeadNodes.push_back(current_node);
//       }
//     }
//
//     return ChannelHeadNodes;
// }
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns all potential channel heads in a DEM. It looks for
// channel heads based on the outlet junctions of the valleys (which are identified by looking
// for portions of the landscape with 10 or more nodes with a high curvature that are linked)
// It returns a vector<int> of nodeindices where the channel heads are
//
// FC 31/10/13
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
vector<int> LSDJunctionNetwork::GetChannelHeadsChiMethodFromValleys(vector<int> ValleyNodes,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster)
{
  vector<int> ChannelHeadNodes;
  vector<int> ChannelHeadNodes_temp;

  int max_nodes = ValleyNodes.size();
  int node_number = 0;

  //loop through junctions collecting channel heads
  for (int i = 0; i < max_nodes; i++)
  {
    cout << flush << "Node = " << i+1 << " of " << max_nodes << "\r";

    // get the junction number
    node_number = ValleyNodes[i];

    // get a local list of channel heads
	  int channel_head_node = GetChannelHeadsChiMethodFromNode(node_number,
                                      			MinSegLength, A_0, m_over_n, FlowInfo,
                                       			FlowDistance, ElevationRaster);

     // now append these channel heads to the master list
    ChannelHeadNodes_temp.push_back(channel_head_node);
  }

  //removing any nodes that are not the furthest upstream
  int upstream_test = 0;
  vector<int>::iterator find_it;

  for (int node = 0; node < int(ChannelHeadNodes_temp.size()); node++)
  {
    vector<int> tests;
    int current_node = ChannelHeadNodes_temp[node];
    for (int i = 0; i < int(ChannelHeadNodes_temp.size()); i++)
    {
      if (ChannelHeadNodes_temp[i] != current_node)
      {
        int test_node = ChannelHeadNodes_temp[i];
        upstream_test = FlowInfo.is_node_upstream(current_node, test_node);
        tests.push_back(upstream_test);
      }
    }
    find_it = find(tests.begin(), tests.end(), 1);
    if (find_it == tests.end())
    {
      ChannelHeadNodes.push_back(current_node);
    }

  }

  cout << "No of source nodes: " << ChannelHeadNodes.size() << endl;

  return ChannelHeadNodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns all potential channel heads in a DEM. It looks for
// channel heads based on the sources from the valley network which identifies concave portions
// of the landscape.  The channel profiles are run from the furthest upslope hilltop to a downstream
// junction: the length of the channel profile can be set by the number of junctions downstream (NJunctions)
//
// FC 10/09/15
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
vector<int> LSDJunctionNetwork::GetChannelHeadsChiMethodFromSources(vector<int> ValleySources,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster, int NJunctions)
{
  vector<int> ChannelHeadNodes;
  vector<int> ChannelHeadNodes_temp;

  int max_nodes = ValleySources.size();
  int node_number = 0;

  //loop through junctions collecting channel heads
  for (int i = 0; i < max_nodes; i++)
  {
    cout << flush << "Node = " << i+1 << " of " << max_nodes << "\r";

    // get the junction number
    node_number = ValleySources[i];

    // get a local list of channel heads
    int channel_head_node = GetChannelHeadsChiMethodFromSourceNode(node_number,
                                      MinSegLength, A_0, m_over_n, FlowInfo,
                                      FlowDistance, ElevationRaster, NJunctions);

     // now append these channel heads to the master list
    ChannelHeadNodes_temp.push_back(channel_head_node);
  }

  //removing any nodes that are not the furthest upstream
  int upstream_test = 0;
  vector<int>::iterator find_it;

  for (int node = 0; node < int(ChannelHeadNodes_temp.size()); node++)
  {
    vector<int> tests;
    int current_node = ChannelHeadNodes_temp[node];
    for (int i = 0; i < int(ChannelHeadNodes_temp.size()); i++)
    {
      if (ChannelHeadNodes_temp[i] != current_node)
      {
        int test_node = ChannelHeadNodes_temp[i];
        upstream_test = FlowInfo.is_node_upstream(current_node, test_node);
        tests.push_back(upstream_test);
      }
    }
    find_it = find(tests.begin(), tests.end(), 1);
    if (find_it == tests.end())
    {
      ChannelHeadNodes.push_back(current_node);
    }

  }

  cout << "No of source nodes: " << ChannelHeadNodes.size() << endl;

  return ChannelHeadNodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns all channels in the DEM that the DrEICH algorithm uses for segiment fitting.
// It looks for channels based on the outlet junctions of valleys.
// It returns a LSDIndexRaster with the channels.
//
// FC 21/08/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::GetChannelsDreich(vector<int> ValleySources,
                                      int MinSegLength, float A_0, float m_over_n,
                                      LSDFlowInfo& FlowInfo, LSDRaster& FlowDistance,
                                      LSDRaster& ElevationRaster, string path_name, int NJunctions)
{
  Array2D<int> channel_nodes(NRows,NCols,NoDataValue);

  int max_nodes = ValleySources.size();
  int node_number = 0;

  //loop through junctions collecting channels
  for (int i = 0; i < max_nodes; i++)
  {
    cout << flush << "Node = " << i+1 << " of " << max_nodes << "\r";

    // get the junction number
    node_number = ValleySources[i];

    // get an index raster with the channel data
    LSDIndexRaster Channel = GetChannelfromDreich(node_number, MinSegLength, A_0, m_over_n,
                         FlowInfo, FlowDistance, ElevationRaster, path_name, NJunctions);
    Array2D<int> ChannelData = Channel.get_RasterData();

    //copy this array to a master array
    for (int row = 0; row < NRows; row++)
    {
      for (int col = 0; col < NCols; col++)
      {
        if(ChannelData[row][col] != NoDataValue)
        {
          channel_nodes[row][col] = ChannelData[row][col];
        }
      }
    }
  }

  //now create the final LSDIndexRaster with the channels
  LSDIndexRaster AllChannels(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, channel_nodes, GeoReferencingStrings);
  return AllChannels;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns a 2D integer array containing the locations of all pixels identified
// as being part of the channel using chi profiles.  It calculates the chi and elevation value
// of every pixel upstream of the given junction, then bins this data and calculates the pixels
// in the 95th percentile of each bin.  Any pixels above the 95th percentile are considered part
// of the channel, and any below are considered to be hillslopes.  This is the first part of the
// channel head prediction using chi profiles.
//
// Parameters: Junction number, A_0, m over n, bin width (suggested value of 10), FlowInfo object,
// Elevation raster
// Returns: Array2D<int> with the channel pixel locations
//
// FC 01/10/13
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
Array2D<int> LSDJunctionNetwork::GetChannelHeadsChiMethodAllPixels(int JunctionNumber,
                                      float A_0, float m_over_n, float bin_width, LSDFlowInfo& FlowInfo,
                                      LSDRaster& ElevationRaster)
{
  Array2D<int> channel_pixels(NRows,NCols,NoDataValue);
  // get the node index of this junction
  int starting_node = JunctionVector[JunctionNumber];
  string jn_name = itoa(JunctionNumber);
  string uscore = "_";
  jn_name = uscore+jn_name;

  //get the chi and elevation values of each upslope node
  vector<int> upslope_nodes = FlowInfo.get_upslope_nodes(starting_node);
  vector<float> upslope_chi = FlowInfo.get_upslope_chi(starting_node, m_over_n, A_0);
  vector<float> elevation;
  int row,col;

  string string_filename_all;
  string filename_all = "chi_profile_all";
  string dot = ".";
  string extension = "txt";
  string_filename_all = filename_all+jn_name+dot+extension;
  ofstream chi_profile_all;
  chi_profile_all.open(string_filename_all.c_str());
  //cout << "The filename is " << string_filename_all << endl;

  for (int node=0; node < int(upslope_nodes.size()); node++)
  {
    FlowInfo.retrieve_current_row_and_col(upslope_nodes[node], row, col);
    float elev = ElevationRaster.get_data_element(row,col);
    elevation.push_back(elev);
    chi_profile_all << upslope_chi[node] << " " << elev << endl;
  }

  string string_filename;
  string filename = "chi_profile";
  string_filename = filename+jn_name+dot+extension;
  //cout << "The filename is " << string_filename << endl;

  ofstream chi_profile;
  chi_profile.open(string_filename.c_str());

  float lower_limit = 0;
  vector<float> mean_chi;
  vector<float> mean_elev;
  vector<float> midpoints;
  vector<float> medianY;
  vector<float> st_dev_chi;
  vector<float> st_dev_elev;
  vector<float> range_min;
  vector<float> range_max;
  vector<int> n_obs;

  float NoDataValue = FlowInfo.get_NoDataValue();

  //bin the data to find the range of the 95th percentile
  bin_data(upslope_chi, elevation, bin_width, mean_chi, mean_elev, midpoints, medianY,
           st_dev_chi, st_dev_elev, range_min, range_max, n_obs, lower_limit, NoDataValue);

  //find the most linear channel segment (highest r2 value)
  int n_bins = mean_chi.size();
  int min_seg_length = n_bins/10;
  vector<float> channel_chi;
  vector<float> channel_elev;
  vector<float>::iterator vec_iter_start;
  vector<float>::iterator vec_iter_end;
  float max_r2 = 0;
  float elev_limit = 0;

  for (int channel_segment = min_seg_length; channel_segment <= n_bins-min_seg_length; channel_segment++)
  {
    //assigning the chi values of the channel segment
    channel_chi.resize(channel_segment);
    vec_iter_start = mean_chi.begin();
    vec_iter_end = vec_iter_start+channel_segment;
    channel_chi.assign(vec_iter_start,vec_iter_end);

    // assigning the elevation values of the channel segment
    channel_elev.resize(channel_segment);
    vec_iter_start = range_min.begin();
    vec_iter_end = vec_iter_start+channel_segment;
    channel_elev.assign(vec_iter_start,vec_iter_end);

    //performing linear regression on channel segment: getting highest r2
    vector<float> residuals_chan;
    vector<float> results_chan = simple_linear_regression(channel_chi,channel_elev, residuals_chan);
    float r2 = results_chan[2];
    if (r2 > max_r2)
    {
      max_r2 = r2;
      elev_limit = channel_elev.back();
    }
  }

  //extend the linear channel segment up through the plot
  vector<float> mean_chi_regression;
  vector<float> range_min_regression;
  vector<float> elev_regression;
  mean_chi_regression.resize(range_min.size());
  range_min_regression.resize(range_min.size());
  elev_regression.resize(range_min.size());
  float regression_pointer = 0;

  for (int i=0; i<int(range_min.size()); i++)
  {
    if (range_min[i] <= elev_limit)
    {
      if (range_min[i] !=0)
      {
        mean_chi_regression[i] = mean_chi[i];
        range_min_regression[i] = range_min[i];
        regression_pointer = i;
      }
    }
  }

  float x1 = mean_chi_regression.front();
  float y1 = range_min_regression.front();
  float x2 = mean_chi_regression[regression_pointer];
  float y2 = range_min_regression[regression_pointer];
  float gradient = (y2 - y1)/(x2 - x1);
  float intercept = y2 - (gradient * x2);

  for (int i = 0; i < n_bins; i++)
  {
    if (range_min[i] <=elev_limit)
    {
      elev_regression[i] = range_min[i];
    }
    else
    {
      elev_regression[i] = mean_chi[i]*gradient + intercept;
    }
  }

  for(int i = 0 ; i< n_bins; i++)
  {
    if (mean_chi[i] != 0)
    {
      chi_profile << mean_chi[i] << " " << mean_elev[i] << " " << range_min[i] << " " << range_max[i] << " " << elev_regression[i] << endl;
    }
  }
  chi_profile.close();

  //classify any nodes to the left of the channel segment as the channel; any nodes to the
  // right are classified as hillslopes
  vector<int> channel_nodes;
  vector<int> upslope_nodes_temp;
  vector<int> source_nodes;
  vector<int>::iterator iterator_find;

  for (int i=0; i < int(upslope_nodes.size()); i++)
  {
    int bin_id = int((upslope_chi[i]-lower_limit)/bin_width);
    FlowInfo.retrieve_current_row_and_col(upslope_nodes[i], row, col);
    if (upslope_chi[i] <= mean_chi[bin_id] && elevation[i] >= elev_regression[bin_id])
    {
      channel_pixels[row][col] = 1;
      channel_nodes.push_back(upslope_nodes[i]);
    }
    else
    {
      channel_pixels[row][col] = 0;
    }
  }

  return channel_pixels;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
// This function returns a 2D integer array containing the locations of all pixels identified
// as being part of the channel using chi profiles.  It calculates the chi and elevation value
// of every pixel upstream of the given junction, then bins this data and calculates the pixels
// in the 95th percentile of each bin.  Any pixels above the 95th percentile are considered part
// of the channel, and any below are considered to be hillslopes.  This is the first part of the
// channel head prediction using chi profiles.
//
// Parameters: Junction number, A_0, m over n, bin width (suggested value of 10), FlowInfo object,
// Elevation raster
// Returns: Array2D<int> with the channel pixel locations
//
// FC 01/10/13
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-
vector<int> LSDJunctionNetwork::GetSourceNodesChiMethodAllPixels(int JunctionNumber,
                                      float A_0, float m_over_n, float bin_width, LSDFlowInfo& FlowInfo,
                                      LSDRaster& ElevationRaster)
{
  // get the node index of this junction
  int starting_node = JunctionVector[JunctionNumber];
  string jn_name = itoa(JunctionNumber);
  string uscore = "_";
  jn_name = uscore+jn_name;

  //get the chi and elevation values of each upslope node
  vector<int> upslope_nodes = FlowInfo.get_upslope_nodes(starting_node);
  vector<float> upslope_chi = FlowInfo.get_upslope_chi(starting_node, m_over_n, A_0);
  vector<float> elevation;
  int row,col;

  for (int node=0; node < int(upslope_nodes.size()); node++)
  {
    FlowInfo.retrieve_current_row_and_col(upslope_nodes[node], row, col);
    float elev = ElevationRaster.get_data_element(row,col);
    elevation.push_back(elev);
  }

  float lower_limit = 0;
  vector<float> mean_chi;
  vector<float> mean_elev;
  vector<float> midpoints;
  vector<float> medianY;
  vector<float> st_dev_chi;
  vector<float> st_dev_elev;
  vector<float> range_min;
  vector<float> range_max;
  vector<int> n_obs;

  float NoDataValue = FlowInfo.get_NoDataValue();

  //bin the data to find the range of the 95th percentile
  bin_data(upslope_chi, elevation, bin_width, mean_chi, mean_elev, midpoints, medianY,
           st_dev_chi, st_dev_elev, range_min, range_max, n_obs, lower_limit, NoDataValue);

  //find the most linear channel segment (highest r2 value)
  int n_bins = mean_chi.size();
  int min_seg_length = n_bins/10;
  vector<float> channel_chi;
  vector<float> channel_elev;
  vector<float>::iterator vec_iter_start;
  vector<float>::iterator vec_iter_end;
  float max_r2 = 0;
  float elev_limit = 0;

  for (int channel_segment = min_seg_length; channel_segment <= n_bins-min_seg_length; channel_segment++)
  {
    //assigning the chi values of the channel segment
    channel_chi.resize(channel_segment);
    vec_iter_start = mean_chi.begin();
    vec_iter_end = vec_iter_start+channel_segment;
    channel_chi.assign(vec_iter_start,vec_iter_end);

    // assigning the elevation values of the channel segment
    channel_elev.resize(channel_segment);
    vec_iter_start = range_min.begin();
    vec_iter_end = vec_iter_start+channel_segment;
    channel_elev.assign(vec_iter_start,vec_iter_end);

    //performing linear regression on channel segment: getting highest r2
    vector<float> residuals_chan;
    vector<float> results_chan = simple_linear_regression(channel_chi,channel_elev, residuals_chan);
    float r2 = results_chan[2];
    if (r2 > max_r2)
    {
      max_r2 = r2;
      elev_limit = channel_elev.back();
    }
  }

  //extend the linear channel segment up through the plot
  vector<float> mean_chi_regression;
  vector<float> range_min_regression;
  vector<float> elev_regression;
  mean_chi_regression.resize(range_min.size());
  range_min_regression.resize(range_min.size());
  elev_regression.resize(range_min.size());
  float regression_pointer = 0;

  for (int i=0; i< int(range_min.size()); i++)
  {
    if (range_min[i] <= elev_limit)
    {
      if (range_min[i] !=0)
      {
        mean_chi_regression[i] = mean_chi[i];
        range_min_regression[i] = range_min[i];
        regression_pointer = i;
      }
    }
  }

  float x1 = mean_chi_regression.front();
  float y1 = range_min_regression.front();
  float x2 = mean_chi_regression[regression_pointer];
  float y2 = range_min_regression[regression_pointer];
  float gradient = (y2 - y1)/(x2 - x1);
  float intercept = y2 - (gradient * x2);

  for (int i = 0; i < n_bins; i++)
  {
    if (range_min[i] <=elev_limit)
    {
      elev_regression[i] = range_min[i];
    }
    else
    {
      elev_regression[i] = mean_chi[i]*gradient + intercept;
    }
  }

  //classify any nodes to the left of the channel segment as the channel; any nodes to the
  // right are classified as hillslopes
  vector<int> channel_nodes;
  vector<int> source_nodes;

  for (unsigned int i=0; i < upslope_nodes.size(); i++)
  {
    int bin_id = int((upslope_chi[i]-lower_limit)/bin_width);
    FlowInfo.retrieve_current_row_and_col(upslope_nodes[i], row, col);
    if (upslope_chi[i] <= mean_chi[bin_id] && elevation[i] >= elev_regression[bin_id])
    {
      channel_nodes.push_back(upslope_nodes[i]);
    }
  }

  // find the furthest upslope nodes classified as being part of the channel network (use as sources for next
  // step of chi method)

  int upstream_test = 0;
  vector<int>::iterator find_it;

  for (unsigned int node =0; node < channel_nodes.size(); node++)
  {
     vector<int> tests;
     int current_node = channel_nodes[node];
     for (unsigned int i = 0; i < channel_nodes.size(); i++)
     {
      if (channel_nodes[i] != current_node)
      {
        int test_node = channel_nodes[i];
        upstream_test = FlowInfo.is_node_upstream(current_node, test_node);
        tests.push_back(upstream_test);
      }
     }
     find_it = find(tests.begin(), tests.end(), 1);
     if (find_it == tests.end())
     {
      source_nodes.push_back(current_node);
     }
  }
  cout << "No of channel nodes: " << channel_nodes.size() << endl;
  cout << "No of source nodes: " << source_nodes.size() << endl;
  return source_nodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// PREDICTING CHANNEL HEADS USING TANGENTIAL CURVATURE
//
// This function is used to predict channel head locations based on the method proposed by
// Pelletier (2013).  It creates a contour curvature map and identifies channel heads as pixels greater
// than a user defined contour curvature threshold value, set by default at 0.1.  The threshold curvature
// can also be defined as a multiple of the standard deviation of the curvature.  Before this function is called
// the DEM must be filtered using the wiener filter in the LSDRasterSpectral object in order to remove high frequency
// noise.
//
// Reference: Pelletier (2013) A robust, two-parameter method for the extraction of drainage
// networks from high-resolution digital elevation models (DEMs): Evaluation using synthetic and real-world
// DEMs, Water Resources Research 49: 1-15
//
// added by FC 16/07/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

vector<int> LSDJunctionNetwork::calculate_pelletier_channel_heads(float tan_curv_threshold, LSDFlowInfo& FlowInfo, Array2D<float>& tan_curv_array)
{
  cout << "Getting Pelletier channel heads" << endl;
  Array2D<float> chan_head_locations(NRows,NCols,NoDataValue);
  float total_curv = 0;
  int n_observations = 0;

  //get the mean of the tangential curvature
  for (int row = 0; row < NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if (tan_curv_array[row][col] != NoDataValue)
      {
        total_curv = total_curv + tan_curv_array[row][col];
        ++n_observations;
      }
    }
  }

  float mean_curv = total_curv/n_observations;
  float total_st_dev = 0;

  // get the standard deviation of the curvature and use 3*st dev as the threshold value
  for (int row = 0; row < NRows; row++)
	{
    for(int col = 0; col < NCols; col++)
    {
      if (tan_curv_array[row][col] != NoDataValue)
      {
        total_st_dev = ((tan_curv_array[row][col] - mean_curv)*(tan_curv_array[row][col] - mean_curv)) + total_st_dev;
      }
    }
  }

  float st_dev = sqrt(total_st_dev/n_observations);
  tan_curv_threshold = 3*st_dev;
  cout << "Got standard deviation" << endl;

  // Get all the locations where the tan curvature is greater than the user defined threshold
  vector<int> channel_nodes;
  int CurrentNodeIndex = 0;

  for (int row = 0; row < NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if (tan_curv_array[row][col] > tan_curv_threshold)
      {
        chan_head_locations[row][col] = tan_curv_array[row][col];
        CurrentNodeIndex = FlowInfo.NodeIndex[row][col];
        channel_nodes.push_back(CurrentNodeIndex);
      }
      else
      {
        chan_head_locations[row][col] = 0;
      }
    }
  }
  cout << "Got channel nodes" << endl;
  // STEP 3: Finding the furthest upstream channel node

  int upstream_test = 0;
  vector<int>::iterator find_it;
  vector<int> source_nodes;

  // identify whether there are any further upstream channel nodes, if not then
  // the node is a channel head
  for (unsigned int node =0; node < channel_nodes.size(); node++)
  {
     cout << flush << "\t\t node " << node << "/" << channel_nodes.size() << "\r";
     vector<int> tests;
     int current_node = channel_nodes[node];
     for (unsigned int i = 0; i < channel_nodes.size(); i++)
     {
      if (channel_nodes[i] != current_node)
      {
        int test_node = channel_nodes[i];
        upstream_test = FlowInfo.is_node_upstream(current_node, test_node);
        tests.push_back(upstream_test);
      }
     }
     find_it = find(tests.begin(), tests.end(), 1);
     if (find_it == tests.end())
     {
      source_nodes.push_back(current_node);
     }
  }
  cout << "Got source nodes" << endl;
  cout << "No of channel nodes: " << channel_nodes.size() << " No of source nodes: " << source_nodes.size() << endl;
  return source_nodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// PREDICTING CHANNEL HEADS USING TANGENTIAL CURVATURE
//
// This function is used to predict channel head locations based on the method
// proposed by Pelletier (2013).  It creates a contour curvature map and identifies
// channel heads as pixels greater than a user defined contour curvature threshold
// value, set by default at 0.1.  The threshold curvature can also be defined as a
// multiple of the standard deviation of the curvature.  Before this function is
// called the DEM must be filtered using the wiener filter in the LSDRasterSpectral
// object in order to remove high frequency noise.
//
// Reference: Pelletier (2013) A robust, two-parameter method for the extraction of
// drainage networks from high-resolution digital elevation models (DEMs): Evaluation
// using synthetic and real-world DEMs, Water Resources Research 49: 1-15
//
// DTM 03/06/2014
//
// Initial function gave a map of pixels with sufficient tangential curvature to
// be designated as a channel.  This map needed to be reduced to give the source
// pixels only.  This is done by i) sorting all the possible sources by
// elevation and ii) routing flow from each potential source using an adaption
// of Freeman MD flow.  Any potential sources that are located on ANY down-slope
// pathway within convergent part of the topography from previously visited
// source pixels are excluded.  This is held within the
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::calculate_pelletier_channel_heads_DTM(LSDFlowInfo& FlowInfo, Array2D<float> topography, float tan_curv_threshold, Array2D<float>& tan_curv_array, Array2D<float>& tan_curv_array_LW)
{
  Array2D<float> curv_array(NRows,NCols,NoDataValue);
  vector<int> possible_sources_row;
  vector<int> possible_sources_col;
  vector<float> possible_sources_elev;
  // Get all the locations where the tan curvature is greater than the user defined threshold
  for (int row = 0; row < NRows; row++)
  {
    for(int col = 0; col < NCols; col++)
    {
      if (tan_curv_array[row][col] > tan_curv_threshold && tan_curv_array_LW[row][col] > 0)
      {
        possible_sources_row.push_back(row);
        possible_sources_col.push_back(col);
        possible_sources_elev.push_back(topography[row][col]);
        curv_array[row][col] = tan_curv_array[row][col];
      }
    }
  }

  vector<size_t> index_map;
  // sort
  matlab_float_sort_descending(possible_sources_elev, possible_sources_elev, index_map);
  matlab_int_reorder(possible_sources_row, index_map, possible_sources_row);
  matlab_int_reorder(possible_sources_col, index_map, possible_sources_col);
  // remove all pixels on downstream pathway
  vector<int> Sources = identify_upstream_limits(FlowInfo, topography, possible_sources_row,possible_sources_col,tan_curv_array_LW);

  return Sources;
}
//------------------------------------------------------------------------------
// This is used to reduce a map of potential sources down to a simplified source
// network for channel extraction by removing potential sources that are on ANY
// downslope pathway from previous sources.  It uses a similar algorithm to the
// Freeman multi-directional flow routing algorithm in the LSDRaster object, but
// restricts flow to within convergent parts of the topography to prevent flow
// crossing noses etc.
// DTM 03/06/2014
vector<int> LSDJunctionNetwork::identify_upstream_limits(LSDFlowInfo& FlowInfo, Array2D<float>& topography, vector<int> source_row_vec,vector<int> source_col_vec, Array2D<float>& tan_curv)
{

  //create output array, populated with nodata
  Array2D<float> area(NRows, NCols, NoDataValue);
  Array2D<int> sources_array(NRows, NCols, int(NoDataValue));
  Array2D<int> times_visited(NRows, NCols, int(NoDataValue));
  //declare variables
  vector<float> flat;
  vector<float> sorted;
  vector<size_t> index_map;
  float one_ov_root_2 = 0.707106781187;
  float p = 1.1; //value avoids preferential flow to diagonals

  //loop through the dem cells creating a row major 1D vector, flat, and
  //setting the cell area to every npn ndv cell
  for (int i = 0; i < NRows; ++i)
  {
    for (int j = 0; j < NCols; ++j)
    {
      flat.push_back(topography[i][j]);
      if (topography[i][j] != NoDataValue)
      {
        area[i][j] = 0;
        times_visited[i][j] = 0;
      }
    }
  }
  int n_possible_sources = source_row_vec.size();
  int row,col;
  for(int i = 0; i<n_possible_sources; ++i)
  {
    row = source_row_vec[i];
    col = source_col_vec[i];
    area[row][col] = 1.0;
    times_visited[row][col] = 1;
  }
  Array2D<int> temp = times_visited.copy() ;
  //sort the 1D elevation vector and produce an index
  matlab_float_sort_descending(flat, sorted, index_map);

  //int i_source = 0;
  for(int q = 0 ;q < int(flat.size()); ++q)
  {

    if (sorted[q] != NoDataValue)
    {
      //use row major ordering to reconstruct each cell's i,j coordinates
      int i = index_map[q] / NCols;
       int j = index_map[q] % NCols;

      //skip edge cells and cells above the source pixel
      if (i != 0 && j != 0 && i != NRows-1 && j != NCols-1){

        //reset variables on each loop
        float total = 0;
        float slope1 = 0;
        float slope2 = 0;
        float slope3 = 0;
        float slope4 = 0;
        float slope5 = 0;
        float slope6 = 0;
        float slope7 = 0;
        float slope8 = 0;

        //Get sum of magnitude of downslope flow, total, and store the magnitude of
        //each of the 8 downslope cells as slope1->8 *Avoids NDVs*
        if (topography[i][j] > topography[i-1][j-1] && topography[i-1][j-1] != NoDataValue){
          slope1 = pow(((topography[i][j] - topography[i-1][j-1]) * one_ov_root_2),p);
          total += slope1;
          if(times_visited[i][j] > 0 && tan_curv[i-1][j-1]>0) times_visited[i-1][j-1]+=1;
        }
        if (topography[i][j] > topography[i-1][j] && topography[i-1][j] != NoDataValue){
          slope2 = pow((topography[i][j] - topography[i-1][j]),p);
          total += slope2;
          if(times_visited[i][j] > 0 && tan_curv[i-1][j]>0) times_visited[i-1][j]+=1;
        }
        if (topography[i][j] > topography[i-1][j+1] && topography[i-1][j+1] != NoDataValue){
          slope3 = pow(((topography[i][j] - topography[i-1][j+1]) * one_ov_root_2),p);
          total += slope3;
          if(times_visited[i][j] > 0 && tan_curv[i-1][j+1]>0) times_visited[i-1][j+1] += 1;
        }
        if (topography[i][j] > topography[i][j+1] && topography[i][j+1] != NoDataValue){
          slope4 = pow((topography[i][j] - topography[i][j+1]),p);
          total += slope4;
          if(times_visited[i][j] > 0 && tan_curv[i][j+1]>0) times_visited[i][j+1] += 1;
        }
        if (topography[i][j] > topography[i+1][j+1] && topography[i+1][j+1] != NoDataValue){
          slope5 = pow(((topography[i][j] - topography[i+1][j+1]) * one_ov_root_2),p);
          total += slope5;
          if(times_visited[i][j] > 0 && tan_curv[i+1][j+1]>0) times_visited[i+1][j+1]+=1;
        }
        if (topography[i][j] > topography[i+1][j] && topography[i+1][j] != NoDataValue){
          slope6 = pow((topography[i][j] - topography[i+1][j]),p);
          total += slope6;
          if(times_visited[i][j] > 0 && tan_curv[i+1][j]>0) times_visited[i+1][j]+=1;
        }
        if (topography[i][j] > topography[i+1][j-1] && topography[i+1][j-1] != NoDataValue){
          slope7 = pow(((topography[i][j] - topography[i+1][j-1]) * one_ov_root_2),p);
          total += slope7;
          if(times_visited[i][j] > 0 && tan_curv[i+1][j-1]>0) times_visited[i+1][j-1]+=1;
        }
        if (topography[i][j] > topography[i][j-1] && topography[i][j-1] != NoDataValue){
          slope8 = pow((topography[i][j] - topography[i][j-1]),p);
          total += slope8;
          if(times_visited[i][j] > 0 && tan_curv[i][j-1]>0) times_visited[i][j-1]+=1;
        }

        //divide slope by total to get the proportion of flow directed to each cell
        //and increment the downslope cells. If no downslope flow to a node, 0 is
        //added, so no change is seen.
        area[i-1][j-1] += (area[i][j] * (slope1)/total);
        area[i-1][j] += (area[i][j] * (slope2)/total);
        area[i-1][j+1] += (area[i][j] * (slope3)/total);
        area[i][j+1] += (area[i][j] * (slope4)/total);
        area[i+1][j+1] += (area[i][j] * (slope5)/total);
        area[i+1][j] += (area[i][j] * (slope6)/total);
        area[i+1][j-1] += (area[i][j] * (slope7)/total);
        area[i][j-1] += (area[i][j] * (slope8)/total);
      }
    }
  }
  vector<int> source_nodes;
  for(int i = 0; i<n_possible_sources; ++i)
  {
    row = source_row_vec[i];
    col = source_col_vec[i];
    if (times_visited[row][col] == 1)
    {
      source_nodes.push_back(FlowInfo.retrieve_node_from_row_and_column(row, col));
      //sources_array[row][col] = 1;
    }
  }
//   //write output LSDRaster object
//   LSDIndexRaster SourcesRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, int(NoDataValue), sources_array);
//   LSDRaster AreaRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, area);
//   AreaRaster.write_raster("area","flt");
  return source_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// IDENTIFYING VALLEYS USING TANGENTIAL CURVATURE
//
// This function is used to identify concave portions of the landscape using a tangential
// curvature threshold. It defines the threshold based on a multiple of the standard deviation
// of the curvature.  It then identifies valleys in which there are a linked series of pixels
// which have a curvature value greater than the threshold, and finds the outlet junction number
// of this valley.  This can be passed to the channel head prediction algorithm using the chi
// method.
// Reference: Peuker, T. K. and D. H. Douglas, (1975), "Detection of surface-specific points by local
// parallel processing of discrete terrain elevation data," Comput. Graphics Image Process., 4: 375-387
//
// added by FC 28/10/13
// edited by FC 18/11/13; put the user-defined parameter of the no connecting nodes into the arguments
// so it can be specified in the parameter file.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Array2D<int> LSDJunctionNetwork::find_valleys(LSDFlowInfo& FlowInfo, Array2D<float>& tan_curv_array, vector<int> sources, int no_connecting_nodes, float tan_curv_threshold)
{
  Array2D<int> NodesVisitedBefore(NRows,NCols,0);
  Array2D<int> NodesVisitedBeforeTemp(NRows,NCols,0);
  Array2D<int> valley_junctions(NRows,NCols,NoDataValue);
  vector<int> valley_nodes;
  //float tan_curv_threshold = 0.1;

  //Find valleys with linked pixels greater than the threshold
  int n_sources = sources.size();
  cout << "No of sources: " << n_sources << endl;

  // Loop through all the sources, moving downstream - keep a count of the number of connected
  // nodes that are above the threshold curvature.  If there are more than 10 nodes that are
  // connected then it is a valley - get the outlet junction of the valley and store in a vector

  for (int source = 0; source < n_sources; source++)
  {
    bool EndofReach = false;
    int max_no_connected_nodes =0;
    int CurrentNode = sources[source];
    int CurrentRow,CurrentCol,ReceiverNode,ReceiverRow,ReceiverCol;

    while (EndofReach == false)
    {
      FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
      FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
      if (tan_curv_array[CurrentRow][CurrentCol] != NoDataValue)
      {
        float node_curvature = tan_curv_array[CurrentRow][CurrentCol];
        //cout << node_curvature << endl;
        NodesVisitedBefore[CurrentRow][CurrentCol] = 1;

        if (node_curvature > tan_curv_threshold)
        {
          ++max_no_connected_nodes;
          //cout << "Max no of connecting nodes: " << max_no_connected_nodes << endl;
        }
        else
        {
          max_no_connected_nodes = 0;
        }

        //check whether the no of connected nodes has been reached; if it has then identify the
        // outlet junction of the valley
        if (max_no_connected_nodes > no_connecting_nodes)
        {
          EndofReach = true;
          int this_node = CurrentNode;
          int current_row,current_col,downslope_node,downslope_row,downslope_col,current_SO,downslope_SO;
          bool reached_outlet = false;
          while (reached_outlet == false)
          {
            FlowInfo.retrieve_current_row_and_col(this_node, current_row, current_col);
            FlowInfo.retrieve_receiver_information(this_node, downslope_node, downslope_row, downslope_col);
            current_SO = StreamOrderArray[current_row][current_col];
            downslope_SO = StreamOrderArray[downslope_row][downslope_col];
            NodesVisitedBeforeTemp[current_row][current_col] = 1;
            bool BeentoReceiver = false;
            if (downslope_SO > current_SO)
            {
              valley_nodes.push_back(this_node);
              int valley_junction = find_upstream_junction_from_channel_nodeindex(this_node, FlowInfo);
              valley_junctions[current_row][current_col] = valley_junction;
              reached_outlet = true;
             // cout << "valley junction: " << valley_junction << endl;
            }
            if (NodesVisitedBeforeTemp[downslope_row][downslope_col] ==1) BeentoReceiver = true;
            if(BeentoReceiver == false)
            {
              //Move downstream
              this_node = downslope_node;
            }
            else
            {
              //Push back the valley node
              reached_outlet = true;
            }
          }
        }

        bool ReceiverVisitedBefore = false;
        // test to see whether we have visited this node before
        if(NodesVisitedBefore[ReceiverRow][ReceiverCol]==1) ReceiverVisitedBefore = true;
        if(ReceiverVisitedBefore == false)
        {
          //Move downstream
          CurrentNode = ReceiverNode;
        }
        else
        {
          //Move to next source
          EndofReach = true;
        }
      }
      else
      {
        EndofReach = true;
      }
    }
  }

  return valley_junctions;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function is used to get the outlet node downstream of a series of source nodes
// It is used for identifying valleys to run the DrEICH algorithm on
//
// FJC 19/08/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDJunctionNetwork::get_outlet_nodes_from_sources(LSDFlowInfo& FlowInfo, vector<int> sources)
{
  Array2D<int> NodesVisitedBefore(NRows,NCols,0);
  Array2D<int> NodesVisitedBeforeTemp(NRows,NCols,0);
  vector<int> valley_nodes;

  for (int i = 0; i < int(sources.size()); i++)
  {
    int this_node = sources[i];
    int current_row,current_col,downslope_node,downslope_row,downslope_col,current_SO,downslope_SO;
    bool reached_outlet = false;
    while (reached_outlet == false)
    {
      FlowInfo.retrieve_current_row_and_col(this_node, current_row, current_col);
      FlowInfo.retrieve_receiver_information(this_node, downslope_node, downslope_row, downslope_col);
      current_SO = StreamOrderArray[current_row][current_col];
      downslope_SO = StreamOrderArray[downslope_row][downslope_col];
      NodesVisitedBeforeTemp[current_row][current_col] = 1;
      bool BeentoReceiver = false;
      int base_level = FlowInfo.is_node_base_level(downslope_node);
      if (base_level == 1)
      {
        valley_nodes.push_back(this_node);
      }
      if (downslope_SO > current_SO)
      {
        valley_nodes.push_back(this_node);
        //int valley_junction = find_upstream_junction_from_channel_nodeindex(this_node, FlowInfo);
        //valley_junctions[current_row][current_col] = valley_junction;
        //reached_outlet = true;
        // cout << "valley junction: " << valley_junction << endl;
      }
      if (NodesVisitedBeforeTemp[downslope_row][downslope_col] ==1) BeentoReceiver = true;
      if(BeentoReceiver == false)
      {
        //Move downstream
        this_node = downslope_node;
      }
      else
      {
        //Push back the valley node
        reached_outlet = true;
      }
    }
  }
  return valley_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function is used to identify concave portions of the landscape using a tangential
// curvature threshold. It defines the threshold based on the standard deviation
// of the curvature.  It then identifies valleys in which there are a linked series of pixels
// which have a curvature value greater than the threshold, and finds the outlet junction number
// of this valley.  This can be passed to the channel head prediction algorithm using the chi
// method.
// Reference: Peuker, T. K. and D. H. Douglas, (1975), "Detection of surface-specific points by local
// parallel processing of discrete terrain elevation data," Comput. Graphics Image Process., 4: 375-387
//
// FJC 20/07/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Array2D<int> LSDJunctionNetwork::find_valleys_adaptive_threshold(LSDFlowInfo& FlowInfo, Array2D<float>& tan_curv_array, vector<int> sources, int no_connecting_nodes, Array2D<float>& tan_curv_threshold)
{
  Array2D<int> NodesVisitedBefore(NRows,NCols,0);
  Array2D<int> NodesVisitedBeforeTemp(NRows,NCols,0);
  Array2D<int> valley_junctions(NRows,NCols,NoDataValue);
  vector<int> valley_nodes;

  //Find valleys with linked pixels greater than the threshold
  int n_sources = sources.size();
  cout << "No of sources: " << n_sources << endl;

  // Loop through all the sources, moving downstream - keep a count of the number of connected
  // nodes that are above the threshold curvature.  If there are more than 10 nodes that are
  // connected then it is a valley - get the outlet junction of the valley and store in a vector

  for (int source = 0; source < n_sources; source++)
  {
    bool EndofReach = false;
    int max_no_connected_nodes =0;
    int CurrentNode = sources[source];
    int CurrentRow,CurrentCol,ReceiverNode,ReceiverRow,ReceiverCol;

    while (EndofReach == false)
    {
      FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
      FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
      if (tan_curv_array[CurrentRow][CurrentCol] != NoDataValue)
      {
        float node_curvature = tan_curv_array[CurrentRow][CurrentCol];
        float threshold_curvature = tan_curv_threshold[CurrentRow][CurrentCol];
        //cout << node_curvature << endl;
        NodesVisitedBefore[CurrentRow][CurrentCol] = 1;

        if (node_curvature > threshold_curvature)
        {
          ++max_no_connected_nodes;
          //cout << "Max no of connecting nodes: " << max_no_connected_nodes << endl;
        }
        else
        {
          max_no_connected_nodes = 0;
        }

        //check whether the no of connected nodes has been reached; if it has then identify the
        // outlet junction of the valley
        if (max_no_connected_nodes > no_connecting_nodes)
        {
          EndofReach = true;
          int this_node = CurrentNode;
          int current_row,current_col,downslope_node,downslope_row,downslope_col,current_SO,downslope_SO;
          bool reached_outlet = false;
          while (reached_outlet == false)
          {
            FlowInfo.retrieve_current_row_and_col(this_node, current_row, current_col);
            FlowInfo.retrieve_receiver_information(this_node, downslope_node, downslope_row, downslope_col);
            current_SO = StreamOrderArray[current_row][current_col];
            downslope_SO = StreamOrderArray[downslope_row][downslope_col];
            NodesVisitedBeforeTemp[current_row][current_col] = 1;
            bool BeentoReceiver = false;
            if (downslope_SO > current_SO)
            {
              valley_nodes.push_back(this_node);
              int valley_junction = find_upstream_junction_from_channel_nodeindex(this_node, FlowInfo);
              valley_junctions[current_row][current_col] = valley_junction;
              reached_outlet = true;
             // cout << "valley junction: " << valley_junction << endl;
            }
            if (NodesVisitedBeforeTemp[downslope_row][downslope_col] ==1) BeentoReceiver = true;
            if(BeentoReceiver == false)
            {
              //Move downstream
              this_node = downslope_node;
            }
            else
            {
              //Push back the valley node
              reached_outlet = true;
            }
          }
        }

        bool ReceiverVisitedBefore = false;
        // test to see whether we have visited this node before
        if(NodesVisitedBefore[ReceiverRow][ReceiverCol]==1) ReceiverVisitedBefore = true;
        if(ReceiverVisitedBefore == false)
        {
          //Move downstream
          CurrentNode = ReceiverNode;
        }
        else
        {
          //Move to next source
          EndofReach = true;
        }
      }
      else
      {
        EndofReach = true;
      }
    }
  }

  return valley_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// DTM 16/07/2015
Array2D<int> LSDJunctionNetwork::find_valleys_using_channel_mask(LSDFlowInfo& FlowInfo, Array2D<int>& channel_mask, vector<int> sources, int no_connecting_nodes)
{
  Array2D<int> NodesVisitedBefore(NRows,NCols,0);
  Array2D<int> NodesVisitedBeforeTemp(NRows,NCols,0);
  Array2D<int> valley_junctions(NRows,NCols,NoDataValue);
  vector<int> valley_nodes;

  //Find valleys with linked pixels greater than the threshold
  int n_sources = sources.size();
  cout << "No of sources: " << n_sources << endl;

  // Loop through all the sources, moving downstream - keep a count of the number of connected
  // nodes that are above the threshold curvature.  If there are more than 10 nodes that are
  // connected then it is a valley - get the outlet junction of the valley and store in a vector

  for (int source = 0; source < n_sources; source++)
  {
    bool EndofReach = false;
    int max_no_connected_nodes =0;
    int CurrentNode = sources[source];
    int CurrentRow,CurrentCol,ReceiverNode,ReceiverRow,ReceiverCol;

    while (EndofReach == false)
    {
      FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
      FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
      if (channel_mask[CurrentRow][CurrentCol] != int(NoDataValue))
      {
        NodesVisitedBefore[CurrentRow][CurrentCol] = 1;

        if (channel_mask[CurrentRow][CurrentCol] == 1) ++max_no_connected_nodes;
        else max_no_connected_nodes = 0;

        //check whether the no of connected nodes has been reached; if it has then identify the
        // outlet junction of the valley
        if (max_no_connected_nodes > no_connecting_nodes)
        {
          EndofReach = true;
          int this_node = CurrentNode;
          int current_row,current_col,downslope_node,downslope_row,downslope_col,current_SO,downslope_SO;
          bool reached_outlet = false;
          while (reached_outlet == false)
          {
            FlowInfo.retrieve_current_row_and_col(this_node, current_row, current_col);
            FlowInfo.retrieve_receiver_information(this_node, downslope_node, downslope_row, downslope_col);
            current_SO = StreamOrderArray[current_row][current_col];
            downslope_SO = StreamOrderArray[downslope_row][downslope_col];
            NodesVisitedBeforeTemp[current_row][current_col] = 1;
            bool BeentoReceiver = false;
            if (downslope_SO > current_SO)
            {
              valley_nodes.push_back(this_node);
              int valley_junction = find_upstream_junction_from_channel_nodeindex(this_node, FlowInfo);
              valley_junctions[current_row][current_col] = valley_junction;
              reached_outlet = true;
             // cout << "valley junction: " << valley_junction << endl;
            }
            if (NodesVisitedBeforeTemp[downslope_row][downslope_col] ==1) BeentoReceiver = true;
            if(BeentoReceiver == false) this_node = downslope_node; //Move downstream
            else  reached_outlet = true; //Push back the valley node
          }
        }

        bool ReceiverVisitedBefore = false;
        // test to see whether we have visited this node before
        if(NodesVisitedBefore[ReceiverRow][ReceiverCol]==1) ReceiverVisitedBefore = true;
        if(ReceiverVisitedBefore == false)
        {
          //Move downstream
          CurrentNode = ReceiverNode;
        }
        else
        {
          //Move to next source
          EndofReach = true;
        }
      }
      else
      {
        EndofReach = true;
      }
    }
  }

  return valley_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function extracts the ridge network by defining it as the co-boundaries
// of basins of equivalent order, for all basin orders within the landscape.
// This is relatively trivial since in each array containing the basins of the
// same order, each basin is labelled with a unique identifier, thus co-
// boundaries are found by locating pixels that neighbour pixels from another
// basin of the same order.
//
// Updated to return an LSDRaster object as ridges can now be assigned CHT values,
// using LSDRaster::RidgeSample, which are not integers. - SWDG
//
// DTM 18/10/2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDJunctionNetwork::ExtractRidges(LSDFlowInfo& FlowInfo)
{
  //LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  Array2D<float> RidgeNetwork(NRows,NCols,NoDataValue);
  // Find maximum stream order
  int MaxStreamOrder = 0;
  for (int i = 0; i < int(StreamOrderVector.size()); ++i)
  {
    if (StreamOrderVector[i] != NoDataValue && StreamOrderVector[i] > MaxStreamOrder)
    {
      MaxStreamOrder = StreamOrderVector[i];
    } // end logic for updating maximum stream order
  }   // end loop through StreamOrderVector

  // Loop through basin orders getting basin and then finding adjacent basin
  // margins
  for (int order = 1; order < MaxStreamOrder + 1; ++order)
  {
    LSDIndexRaster basins = ExtractBasinsOrder(order, FlowInfo);
    for (int row = 0; row < NRows; ++row)
    {
      for (int col = 0; col < NCols; ++col)
      {
        if (basins.get_data_element(row, col) != NoDataValue)
        {
          // Deal with corners
          if (row == 0 && col == 0)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                && basins.get_data_element(row+1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                  && basins.get_data_element(row+1, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          else if (row == 0 && col == NCols - 1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                && basins.get_data_element(row+1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                  && basins.get_data_element(row+1, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          else if (row == NRows - 1  && col == 0)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                && basins.get_data_element(row-1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          else if (row == NRows - 1  && col == NCols - 1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                && basins.get_data_element(row-1, col) != NoDataValue) ||
              (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                 && basins.get_data_element(row, col-1) != NoDataValue) ||
              (basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          // Edge pixels
         else if (row == 0)
         {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                && basins.get_data_element(row+1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                  && basins.get_data_element(row+1, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                  && basins.get_data_element(row+1, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          else if (row == NRows -1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                && basins.get_data_element(row-1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          else if (col == 0)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                && basins.get_data_element(row+1, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                  && basins.get_data_element(row-1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                  && basins.get_data_element(row+1, col) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }
          else if (col == NCols - 1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                && basins.get_data_element(row+1, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                  && basins.get_data_element(row-1, col) != NoDataValue) ||
              	(basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                  && basins.get_data_element(row+1, col) != NoDataValue))
            {
              RidgeNetwork[row][col] = 1;
            }
          }

          // Non edge pixels
          else if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                	&& basins.get_data_element(row-1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                	&& basins.get_data_element(row-1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
               	 && basins.get_data_element(row, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                	&& basins.get_data_element(row, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                	&& basins.get_data_element(row+1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                	&& basins.get_data_element(row+1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                	&& basins.get_data_element(row+1, col+1) != NoDataValue))
          {
            RidgeNetwork[row][col] = 1;
          }
        } // end logic for locating basin margins
      }   // end loop through col
    }     // end loop through row
  }       // end loop through different basin orders.
  // Return raster with all nth order drainage basins.
  	LSDRaster ridge_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,RidgeNetwork,GeoReferencingStrings);
	return ridge_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This overloaded function extracts the ridge network for a defined stream
// order, passed in by the user.
//
// Updated to return an LSDRaster object as ridges can now be assigned CHT values,
// using LSDRaster::RidgeSample, which are not integers. - SWDG
//
// DTM 18/10/2012
// SWDG 28/03/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDJunctionNetwork::ExtractRidges(LSDFlowInfo& FlowInfo, int& min_order, int& max_order)
{
  //LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  Array2D<float> RidgeNetwork(NRows,NCols,NoDataValue);
  // Find maximum stream order
  int MaxStreamOrder = 0;
  for (int i = 0; i < int(StreamOrderVector.size()); ++i)
  {
    if (StreamOrderVector[i] != NoDataValue && StreamOrderVector[i] > MaxStreamOrder)
    {
      MaxStreamOrder = StreamOrderVector[i];
    } // end logic for updating maximum stream order
  }   // end loop through StreamOrderVector


  // Check min_order is lower than or equal to max_order, if test fails,
  // set min_order to 1
  if (min_order > max_order){
      min_order = 1;
      cout << "\tmin_order larger than max_order, min_order set to 1" << endl;
  }

  // test the upper limit against the maximum possible value and set
  // MaxStreamOrder to the user defined value. If user supplied value
  // is too high, use original MaxStreamOrder value.
  if (max_order < MaxStreamOrder){
      MaxStreamOrder = max_order;
  }
  else if(max_order > MaxStreamOrder){
      cout << "\tmax_order exceeds stream orders found in dem, max_order set to "
      << MaxStreamOrder << endl;
  }

  // Loop through basin orders getting basin and then finding adjacent basin
  // margins
  for (int order = min_order; order < MaxStreamOrder + 1; ++order)
  {
    LSDIndexRaster basins = ExtractBasinsOrder(order, FlowInfo);
    for (int row = 0; row < NRows; ++row)
    {
      for (int col = 0; col < NCols; ++col)
      {
        if (basins.get_data_element(row, col) != NoDataValue)
        {
          // Deal with corners
          if (row == 0 && col == 0)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                && basins.get_data_element(row+1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                  && basins.get_data_element(row+1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          else if (row == 0 && col == NCols - 1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                && basins.get_data_element(row+1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                  && basins.get_data_element(row+1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          else if (row == NRows - 1  && col == 0)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                && basins.get_data_element(row-1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          else if (row == NRows - 1  && col == NCols - 1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                && basins.get_data_element(row-1, col) != NoDataValue) ||
              (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                 && basins.get_data_element(row, col-1) != NoDataValue) ||
              (basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          // Edge pixels
         else if (row == 0)
         {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                && basins.get_data_element(row+1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                  && basins.get_data_element(row+1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                  && basins.get_data_element(row+1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          else if (row == NRows -1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                && basins.get_data_element(row-1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          else if (col == 0)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                && basins.get_data_element(row+1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                  && basins.get_data_element(row-1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                  && basins.get_data_element(row+1, col) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }
          else if (col == NCols - 1)
          {
            if ((basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                && basins.get_data_element(row+1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                  && basins.get_data_element(row-1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                  && basins.get_data_element(row+1, col) != NoDataValue))
            {
              RidgeNetwork[row][col] = order;
            }
          }

          // Non edge pixels
          else if ((basins.get_data_element(row, col) != basins.get_data_element(row-1, col-1)
                  && basins.get_data_element(row-1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col)
                  && basins.get_data_element(row-1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row-1, col+1)
                  && basins.get_data_element(row-1, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col-1)
                  && basins.get_data_element(row, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row, col+1)
                  && basins.get_data_element(row, col+1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col-1)
                  && basins.get_data_element(row+1, col-1) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col)
                  && basins.get_data_element(row+1, col) != NoDataValue) ||
                (basins.get_data_element(row, col) != basins.get_data_element(row+1, col+1)
                  && basins.get_data_element(row+1, col+1) != NoDataValue))
          {
            RidgeNetwork[row][col] = order;
          }
        } // end logic for locating basin margins
      }   // end loop through col
    }     // end loop through row
  }       // end loop through different basin orders.
  // Return raster with all nth order drainage basins.
  LSDRaster ridge_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,RidgeNetwork,GeoReferencingStrings);
	return ridge_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// ExtractHilltops
//------------------------------------------------------------------------------
// Resticts ridgeline to part of ridge network where the slope is less than a
// threshold value
//
// Now outputs an LSDRaster to fall in line with other hilltop tools - SWDG 29/8/13
//
// DTM 01/04/2013
//------------------------------------------------------------------------------
LSDRaster LSDJunctionNetwork::ExtractHilltops(LSDRaster& RidgeRaster, LSDRaster& SlopeRaster, float MaxSlope)
{
  Array2D<float> Hilltops(NRows,NCols,NoDataValue) ;
  for (int row = 0; row < NRows; ++row)
  {
    for (int col = 0; col < NCols; ++col)
    {
      if (SlopeRaster.get_data_element(row,col) < MaxSlope)
      {
        Hilltops[row][col] = RidgeRaster.get_data_element(row,col);
      }
    }
  }
  LSDRaster hilltop_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,Hilltops,GeoReferencingStrings);
	return hilltop_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function iterates through the junction nodes and assigns the unique
// junction ID (q) to every stream pixel. This can be used with the LSDRaster
// hilltop_flow_routing function to assign a unique ID to each hilltop
// section tying it to a specific section of the channel network.
//
// Takes a flowinfo object and returns an LSDIndexRaster of the indexed channel
// newtork.
//
// SWDG - 04/04/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::ChannelIndexer(LSDFlowInfo& flowinfo)
{
  Array2D<int> StreamOutput(NRows,NCols,NoDataValue);

  int g = 0;  //ints to store the row and col of the current px
  int h = 0;

  for (int q = 1; q < NJunctions; ++q){

    if (q % 100 == 0){
      cout << "\tJunction = " << q << " / " << NJunctions << "    \r";
    }

    int sourcenodeindex = JunctionVector[q]; //first cell of segment
    int recieverjunction = ReceiverVector[q];
    int recievernodeindex = JunctionVector[recieverjunction]; //last cell of segment

    //get row and col of last px in junction. This location should not be written,
    //as it is the start of a new junction.
    int a = 0;
    int b = 0;
    flowinfo.retrieve_current_row_and_col(recievernodeindex,a,b);

    //write first pixel
    flowinfo.retrieve_current_row_and_col(sourcenodeindex,g,h);
    StreamOutput[g][h] = q;

    bool Flag = false; //Flag used to indicate if end of stream segemnt has been reached
    int CurrentNodeIndex = 0;
    int next_receiver;
    while(Flag == false)
    {

      CurrentNodeIndex = flowinfo.NodeIndex[g][h]; //update node index to move 1 px downstream
      flowinfo.retrieve_receiver_information(CurrentNodeIndex, next_receiver, g, h);

      if (CurrentNodeIndex == next_receiver)
      {
        //need to stop 1 px before node
        //cout << "I found the base level" << endl;
        Flag = true;
      }
      else if(recievernodeindex== next_receiver)
      {
        //cout << "I found the receiver" << endl;
        Flag = true;
      }
      else
      {
          StreamOutput[g][h] = q;
      }
    }
  }
  cout << endl;
  LSDIndexRaster IndexedChannels(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, StreamOutput,GeoReferencingStrings);
  return IndexedChannels;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function writes three vector: a vector of node indices,
// a vector of junction indices and a vector of stream orders
// It is used to create an ordered channel vector that can be combined
// with other methods to produce a channel network
// csv file.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::GetChannelNodesAndJunctions(LSDFlowInfo& flowinfo, vector<int>& NIvec, vector<int>& JIvec, vector<int>& SOvec)
{
  vector<int> NI_vector;
  vector<int> JI_vector;
  vector<int> SO_vector;

  int row = 0;  //ints to store the row and col of the current px
  int col = 0;

  //cout << "I am going to go through " << NJunctions << " Junctions for you." << endl;

  for (int q = 0; q < NJunctions; ++q)
  {

    if (q % 100 == 0)
    {
      cout << "\tJunction = " << q << " / " << NJunctions << "    \r";
    }
    //cout << "Junction is: " << q << " ";

    int sourcenodeindex = JunctionVector[q]; //first cell of segment
    int recieverjunction = ReceiverVector[q];
    int recievernodeindex = JunctionVector[recieverjunction]; //last cell of segment

    //cout << "Source NI : " << sourcenodeindex << " and receiver NI: " << recievernodeindex << endl;

    //cout << "reciever is:"  <<  recieverjunction << " ";
    //get row and col of last px in junction. This location should not be written,
    //as it is the start of a new junction.
    int lp_row = 0;
    int lp_col = 0;
    flowinfo.retrieve_current_row_and_col(recievernodeindex,lp_row,lp_col);

    //write first pixel
    flowinfo.retrieve_current_row_and_col(sourcenodeindex,row,col);
    NI_vector.push_back(sourcenodeindex);
    JI_vector.push_back(q);
    SO_vector.push_back(StreamOrderArray[row][col]);

    bool Flag = false; //Flag used to indicate if end of stream segemnt has been reached
    int CurrentNodeIndex = 0;

    if(recieverjunction == q)
    {
      //cout << "You are on a baselevel junction" << endl;
      Flag = true;
    }

    int next_receiver;
    while(Flag == false)
    {

      CurrentNodeIndex = flowinfo.NodeIndex[row][col]; //update node index to move 1 px downstream
      flowinfo.retrieve_receiver_information(CurrentNodeIndex, next_receiver, row, col);

      //cout << "CNI: " <<  CurrentNodeIndex << " and RNI: " << recievernodeindex << endl;

      if (CurrentNodeIndex == next_receiver)
      {
        //need to stop 1 px before node
        //cout << "I found the base level" << endl;
        Flag = true;
      }
      else if(recievernodeindex== next_receiver)
      {
        //cout << "I found the receiver" << endl;
        Flag = true;
      }
      else
      {
        NI_vector.push_back(next_receiver);
        JI_vector.push_back(q);
        SO_vector.push_back(StreamOrderArray[row][col]);
      }
    }
  }
  //cout << "Okay, I've got the nodes" << endl;


  NIvec = NI_vector;
  JIvec = JI_vector;
  SOvec = SO_vector;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// Return a map of all nodes in the channel, useful to just chek if me node is a CNode or not
// B.G. 12/11/2018
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
map<int,bool> LSDJunctionNetwork::GetMapOfChannelNodes(LSDFlowInfo& flowinfo)
{

  map<int,bool> is_node_in_channel;
  int row = 0;  //ints to store the row and col of the current px
  int col = 0;

  //cout << "I am going to go through " << NJunctions << " Junctions for you." << endl;

  for (int q = 0; q < NJunctions; ++q)
  {

    if (q % 100 == 0)
    {
      cout << "\tJunction = " << q << " / " << NJunctions << "    \r";
    }
    //cout << "Junction is: " << q << " ";

    int sourcenodeindex = JunctionVector[q]; //first cell of segment
    int recieverjunction = ReceiverVector[q];
    int recievernodeindex = JunctionVector[recieverjunction]; //last cell of segment

    //cout << "Source NI : " << sourcenodeindex << " and receiver NI: " << recievernodeindex << endl;

    //cout << "reciever is:"  <<  recieverjunction << " ";
    //get row and col of last px in junction. This location should not be written,
    //as it is the start of a new junction.
    int lp_row = 0;
    int lp_col = 0;
    flowinfo.retrieve_current_row_and_col(recievernodeindex,lp_row,lp_col);

    //write first pixel
    flowinfo.retrieve_current_row_and_col(sourcenodeindex,row,col);
    is_node_in_channel[sourcenodeindex] = true;

    bool Flag = false; //Flag used to indicate if end of stream segemnt has been reached
    int CurrentNodeIndex = 0;

    if(recieverjunction == q)
    {
      //cout << "You are on a baselevel junction" << endl;
      Flag = true;
    }

    int next_receiver;
    while(Flag == false)
    {

      CurrentNodeIndex = flowinfo.NodeIndex[row][col]; //update node index to move 1 px downstream
      flowinfo.retrieve_receiver_information(CurrentNodeIndex, next_receiver, row, col);

      //cout << "CNI: " <<  CurrentNodeIndex << " and RNI: " << recievernodeindex << endl;

      if (CurrentNodeIndex == next_receiver)
      {
        //need to stop 1 px before node
        //cout << "I found the base level" << endl;
        Flag = true;
      }
      else if(recievernodeindex== next_receiver)
      {
        //cout << "I found the receiver" << endl;
        Flag = true;
      }
      else
      {
        is_node_in_channel[next_receiver] = true;
      }
    }
  }
  //cout << "Okay, I've got the nodes" << endl;


  return is_node_in_channel;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// SplitChannel
// This function splits the channel into a series of segments, providing a
// convenient unit with which to analyse landscapes.  The user provides the
// TargetSegmentLength, which specifies how many nodes should be in each
// segment, and a MinimumSegmentLength, which specifies the fewest permissable
// number of nodes.  Segments smaller than this are amalgamated into the
// upstream segment.
// The algorithm loops through the sources and traces downstream, stopping a
// segment after the target segment length, when the stream order increases (to
// preserve structure of drainage network), or when a channel pixel has already
// been visited.
//
// DTM 23/10/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::SplitChannel(LSDFlowInfo& FlowInfo, vector<int> Sources, int TargetSegmentLength)//, int MinimumSegmentLength)
{
  //LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  Array2D<int> ChannelSegments(NRows,NCols,int(NoDataValue));
  Array2D<int> NodesVisitedBefore(NRows,NCols,0);
  //----------------------------------------------------------------------------
  //
  int SegmentID = 0;
  int N_Sources = Sources.size();
  // Loop through sources
  for (int i_source = 0; i_source < N_Sources; ++i_source)
  {
    bool EndOfReach = false;
    int NodeCount = 0;
    vector<int> ChannelNodesInSegment;
    int CurrentNode = Sources[i_source];
    int CurrentRow,CurrentCol,ReceiverNode,ReceiverRow,ReceiverCol,CurrentStreamOrder,ReceiverStreamOrder;
    // Trace downstream until you rach the end of this channel reach
    while(EndOfReach == false)
    {
      FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
      FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);

      ChannelSegments[CurrentRow][CurrentCol] = SegmentID;
      NodesVisitedBefore[CurrentRow][CurrentCol] = 1;
      ++NodeCount;

      // Now check whether we have enough channel nodes
      if(NodeCount >= TargetSegmentLength)
      {
        ++SegmentID;
        NodeCount = 0;
      }
      // Now check to see whether stream order increases (want to start a new
      // segment if this is the case)
      ReceiverStreamOrder = StreamOrderArray[ReceiverRow][ReceiverCol];
      CurrentStreamOrder = StreamOrderArray[CurrentRow][CurrentCol];
      if (ReceiverStreamOrder > CurrentStreamOrder)
      {
        NodeCount = 0;
        ++SegmentID;
      }

      bool ReceiverVisitedBefore = false;
      // test to see whether we have visited this node before
      if(NodesVisitedBefore[ReceiverRow][ReceiverCol]==1) ReceiverVisitedBefore = true;

      if(ReceiverVisitedBefore == true)
      {
        EndOfReach = true;
        ++SegmentID;
        ++i_source;
      }
      else
      {
        // Move downstream
        CurrentNode = ReceiverNode;
      }
    }
  }
  LSDIndexRaster ChannelSegmentsRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ChannelSegments,GeoReferencingStrings);
  return ChannelSegmentsRaster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// TypologyModel
// This function splits the channel into a series of segments, providing a
// convenient unit with which to analyse landscapes.  Function modified from
// original SplitChannel function so that the segment length varies with the
// drainage area of the catchment.
// Length (m) is calculated based on: L = Min_reach_length * sqrt(Drainage area (km))
// User must pass in the minimum reach length in metres
// The algorithm loops through the sources and traces downstream, stopping a
// segment after the target segment length, when the stream order increases (to
// preserve structure of drainage network), or when a channel pixel has already
// been visited.
// User must pass in an empty IndexRaster which will be populated with the channel
// segments data, and two vector of vectors which will be populated:
// vector< vector<int> > SegmentInfoInts has the following layout:
// 0 - segment IDS
// 1 - start node of each segment (upstream)
// 2 - end nodes (downstream)
// vector< vector<float> > SegmentInfoFloats has the following layout:
// 0 - segment lengths
// 1 - flow length of each segment
// 2 - elevation of the start nodes
// 3 - slope of the segment
// 4 - slope of the segment based on flow length
// 5 - discharge of the segment
// 6 - transport capacity of the segment
// 7 - sediment supply of the segment
// 8 - algorithm value of the segment
//
// Modified FJC 06/02/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::TypologyModel(LSDFlowInfo& FlowInfo, vector<int> Sources, vector<int> BaselineSources, vector<int> CatchIDs, vector<int> HydroCodes, int MinReachLength, int search_radius, LSDRaster& ElevationRaster, LSDRaster& DischargeRaster, LSDIndexRaster& ChannelSegmentsRaster, vector< vector<int> >& SegmentInfoInts, vector< vector<float> >& SegmentInfoFloats)
{
  //vectors for storing information about the segments
  vector<int> SegmentIDs;       // ID of each segment
  vector<int> StartNodes;       // start node (upstream) of each segment
  vector<int> EndNodes;         // end node (downstream) of each segment
  vector<int> CatchIDs_final;   // catch ID of each DRN source
  vector<int> HydroCodes_final; // hydro code of each DRN source
  vector<float> SegmentLengths; // length of each segment
  vector<float> FlowLengths;    // flow length of each segment
  vector<float> Elevations;     // elevation of each start node
  vector<float> Slopes;         // slope of each segment
  vector<float> Slopes_Lf;      // slope based on flow length
  vector<float> Discharges;     // discharge of each segment
  vector<float> TransportCapacities; // transport capacity (Q_c) of each segment: Q*S
  vector<float> SedimentSupplies;    // sediment supply: Q_c of upstream reach * ratio of upstream length to current length
  vector<float> AVs;            // algorithm value, sqrt((Q_C)^2 + (Q_S)^2)

  //LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  Array2D<int> ChannelSegments(NRows,NCols,int(NoDataValue));
  Array2D<int> NodesVisitedBefore(NRows,NCols,0);
  Array2D<int> HydroCodes_array(NRows,NCols,int(NoDataValue));
  // need these arrays to deal with sediment supply
  Array2D<float> SedimentSupply(NRows,NCols,NoDataValue);
  Array2D<float> Distances(NRows,NCols,NoDataValue);
  //----------------------------------------------------------------------------
  //
  int SegmentID = 0;
  int N_Sources = Sources.size();
  // Loop through sources
  for (int i_source = 0; i_source < N_Sources; ++i_source)
  {
    int node_test = 0;
    int CatchID, HydroCode;
    // thin to the baseline channel network
    for (int j = 0; j < int(BaselineSources.size()); ++j)
    {
      // check if the source is the same as the baseline
      if (Sources[i_source] == BaselineSources[j]) {node_test = 1;}
      else
      {
        //check if node is downstream or upstream of the source
        int upstream_test = FlowInfo.is_node_upstream(Sources[i_source], BaselineSources[j]);
        int downstream_test = FlowInfo.is_node_upstream(BaselineSources[j], Sources[i_source]);
        if (upstream_test == 1 || downstream_test == 1) { node_test = 1; }

        //get stream order of baseline
        //int BaselineSO = FlowInfo.get_StreamOrder_of_Node(BaselineSources[j]);
      }
      if (node_test == 1)
      {
        //assign the catchID and hydrocode of this source
        cout << "Catch ID: " << CatchIDs[j] << " Hydro code: " << HydroCodes[j] << endl;
        CatchID = CatchIDs[j];
        HydroCode = HydroCodes[j];
      }
    }
    if (node_test == 1)
    {
      bool EndOfReach = false;
      //int NodeCount = 0;
      vector<int> ChannelNodesInSegment;
      int CurrentNode = Sources[i_source];
      int CurrentRow,CurrentCol,ReceiverNode,ReceiverRow,ReceiverCol,CurrentStreamOrder,ReceiverStreamOrder;
      FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
      // get the drainage area of the sources
      float ThisArea = FlowInfo.get_DrainageArea_square_km(CurrentNode);
      // calculate the segment length from the drainage area.
      float SegmentLength = MinReachLength*sqrt(ThisArea);
      // check that seg length is at least the minimum
      if (SegmentLength < MinReachLength) { SegmentLength = MinReachLength; }

      //push back the source to the vector of starting nodes
      int ThisStartNode = CurrentNode;
      int ThisStartRow = CurrentRow;
      int ThisStartCol = CurrentCol;
      StartNodes.push_back(ThisStartNode);
      SegmentIDs.push_back(SegmentID);
      float ThisElev = ElevationRaster.get_data_element(CurrentRow,CurrentCol);
      Elevations.push_back(ThisElev);
      //SegmentLengths.push_back(SegmentLength);

      float ThisDischarge;
      //find the nearest discharge value to this nodes
      if (DischargeRaster.get_data_element(CurrentRow,CurrentCol) != NoDataValue)
      {
        ThisDischarge = (DischargeRaster.get_data_element(CurrentRow,CurrentCol))*0.001;  // convert discharge to m^3/s
      }
      else { ThisDischarge = float(NoDataValue); }
      Discharges.push_back(ThisDischarge);

      //push back the catch IDs and HydroCodes
      CatchIDs_final.push_back(CatchID);
      HydroCodes_final.push_back(HydroCode);

      // Trace downstream until you rach the end of this channel reach
      while(EndOfReach == false)
      {
        FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
        FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
        // get the elevation information
        float ReceiverElev = ElevationRaster.get_data_element(ReceiverRow,ReceiverCol);

        ChannelSegments[CurrentRow][CurrentCol] = SegmentID;
        NodesVisitedBefore[CurrentRow][CurrentCol] = 1;

        // get the distance to check the segment length
        float ThisDistance = FlowInfo.get_Euclidian_distance(ThisStartNode,CurrentNode);

        ReceiverStreamOrder = StreamOrderArray[ReceiverRow][ReceiverCol];
        CurrentStreamOrder = StreamOrderArray[CurrentRow][CurrentCol];

        bool ReceiverVisitedBefore = false;
        // test to see whether we have visited this node before
        if(NodesVisitedBefore[ReceiverRow][ReceiverCol]==1) ReceiverVisitedBefore = true;

        // Now check whether we have a long enough segment or whether stream order
        // increases (want to start a new segment if this is the case)
        if((ThisDistance >= SegmentLength || ReceiverStreamOrder > CurrentStreamOrder) && ReceiverVisitedBefore == false)
        {
          ++SegmentID;
          //get the slope of this reach and push back to vector
          float ReachSlope = (ThisElev - ReceiverElev)/SegmentLength;

          // get the flow length of this reach
          float FlowLength = FlowInfo.get_flow_length_between_nodes(ThisStartNode, CurrentNode);

          float Slope_Lf = float(NoDataValue);
          // get the slope based on the flow length
          if (FlowLength > 0) { Slope_Lf = (ThisElev - ReceiverElev)/FlowLength; }

          // get discharge and push back to vector
          float ReceiverDischarge;
          if (DischargeRaster.get_data_element(ReceiverRow,ReceiverCol) != NoDataValue)
          {
            ReceiverDischarge = (DischargeRaster.get_data_element(ReceiverRow,ReceiverCol))*0.001;  // convert to m^3/s
          }
          else { ReceiverDischarge = float(NoDataValue); }

          // recalculate the segment length
          ThisArea = FlowInfo.get_DrainageArea_square_km(ReceiverNode);
          float NewSegmentLength = MinReachLength*sqrt(ThisArea);
          if (NewSegmentLength < MinReachLength) { NewSegmentLength = MinReachLength; }

          // get the transport capacity and sediment supply
          float TC = float(NoDataValue);
          if (ThisDischarge != NoDataValue && ReachSlope != NoDataValue) { TC = ThisDischarge*ReachSlope; }

          // get the sediment supply. We push this to an array because if we are at a junction
          // then we will need to add the sediment supply from upstream afterwards.
          //cout << "Getting upstream values..." << endl;
          if (ThisStartNode != Sources[i_source])
          {
            float UpstreamDistance = SegmentLengths.back();
            float UpstreamTC = TransportCapacities.back();
            //cout << "Got the upstream values" << endl;
            if (UpstreamTC != NoDataValue && ThisDistance > 0)
            {
              float ThisQs = UpstreamTC * (UpstreamDistance/ThisDistance);
              SedimentSupply[ThisStartRow][ThisStartCol] = ThisQs;
              Distances[ThisStartRow][ThisStartCol] = ThisDistance;
            }
          }

          // get the catchment id and hydrocodes. Need to check the source points
          // for hydrocodes and push back the one with the lowest value
          int LowestHC = 1000000;
          for (int i = 0; i < int(BaselineSources.size()); ++i)
          {
            int upstream_test = FlowInfo.is_node_upstream(ReceiverNode, BaselineSources[i]);
            if (upstream_test == 1)
            {
              int This_HC = HydroCodes[i];
              //cout << "THIS HC: " << This_HC << endl;
              if (This_HC < LowestHC)
              {
                LowestHC = This_HC;
              }
            }
          }
          if (LowestHC == 1000000) { LowestHC = int(NoDataValue); }
          CatchIDs_final.push_back(CatchID);
          HydroCodes_final.push_back(LowestHC);

          //push back info to vectors
          EndNodes.push_back(CurrentNode);
          StartNodes.push_back(ReceiverNode);
          SegmentIDs.push_back(SegmentID);

          SegmentLengths.push_back(ThisDistance);
          FlowLengths.push_back(FlowLength);
          Elevations.push_back(ReceiverElev);
          Slopes.push_back(ReachSlope);
          Slopes_Lf.push_back(Slope_Lf);
          Discharges.push_back(ReceiverDischarge);
          TransportCapacities.push_back(TC);

          //update variables
          ThisStartNode = ReceiverNode;
          ThisElev = ReceiverElev;
          ThisDischarge = ReceiverDischarge;
          SegmentLength = NewSegmentLength;
          ThisStartRow = ReceiverRow;
          ThisStartCol = ReceiverCol;

        }

        if(ReceiverVisitedBefore == true)
        {
          EndOfReach = true;
          ++SegmentID;

          //get the slope of this reach and push back to vector
          float ReachSlope = (ThisElev - ReceiverElev)/SegmentLength;

          // get the flow length of this reach
          float FlowLength = FlowInfo.get_flow_length_between_nodes(ThisStartNode, CurrentNode);

          float Slope_Lf = float(NoDataValue);
          // get the slope based on the flow length
          if (FlowLength > 0) { Slope_Lf = (ThisElev - ReceiverElev)/FlowLength; }

          // recalculate the segment length
          ThisArea = FlowInfo.get_DrainageArea_square_km(ReceiverNode);
          float NewSegmentLength = MinReachLength*sqrt(ThisArea);
          if (NewSegmentLength < MinReachLength) { NewSegmentLength = MinReachLength; }

          // get the transport capacity and sediment supply
          float TC = float(NoDataValue);
          if (ThisDischarge != NoDataValue && ReachSlope != NoDataValue) { TC = ThisDischarge*ReachSlope; }

          // get the sediment supply. We push this to an array because if we are at a junction
          // then we will need to add the sediment supply from upstream afterwards.
          if (SedimentSupply[ReceiverRow][ReceiverCol] != NoDataValue && TC != NoDataValue)
          {
            float DownstreamDistance = Distances[ReceiverRow][ReceiverCol];
            float ReceiverQS = TC * (ThisDistance/DownstreamDistance);
            if (ReceiverQS > 0) { SedimentSupply[ReceiverRow][ReceiverCol] += ReceiverQS; }
          }

          //push back info to vectors
          EndNodes.push_back(CurrentNode);
          SegmentLengths.push_back(ThisDistance);
          FlowLengths.push_back(FlowLength);
          Slopes.push_back(ReachSlope);
          Slopes_Lf.push_back(Slope_Lf);
          TransportCapacities.push_back(TC);
        }
        else
        {
          // Move downstream
          CurrentNode = ReceiverNode;
        }
      }
    }
  }


  // getting the sediment supply from array to vector. This is annoying but I'm too dumb to think of a more efficient way of doing this.
  for (int i = 0; i < int(StartNodes.size()); i++)
  {
    int row, col;
    FlowInfo.retrieve_current_row_and_col(StartNodes[i], row, col);
    float ThisQs = SedimentSupply[row][col];
    SedimentSupplies.push_back(ThisQs);

    //calculate the AV for each Q_s and Q_c
    float AV = float(NoDataValue);
    if (TransportCapacities[i] != NoDataValue && ThisQs != NoDataValue)
    {
      AV = sqrt((TransportCapacities[i]*TransportCapacities[i]) + (ThisQs*ThisQs));
    }
    AVs.push_back(AV);
  }


  cout << "Checking segment node finder, n_start nodes: " << StartNodes.size() << " N_end nodes: " << EndNodes.size() << " N_segment ids: " << SegmentIDs.size() << " Final segment ID: " << SegmentID << " N elevations: " << Elevations.size() << " N slopes: " << Slopes.size() << " N slopes lf: " << Slopes_Lf.size() << " N segment lengths: " << SegmentLengths.size() << " N flow lengths: " << FlowLengths.size() << " N TransportCapacities: " << TransportCapacities.size() << " N sediment supplies: " << SedimentSupplies.size() << " N algorithm values: " << AVs.size() << endl;

  //push back to master vectors
  SegmentInfoInts.push_back(SegmentIDs);
  SegmentInfoInts.push_back(StartNodes);
  SegmentInfoInts.push_back(EndNodes);
  SegmentInfoInts.push_back(CatchIDs_final);
  SegmentInfoInts.push_back(HydroCodes_final);

  SegmentInfoFloats.push_back(SegmentLengths);
  SegmentInfoFloats.push_back(FlowLengths);
  SegmentInfoFloats.push_back(Elevations);
  SegmentInfoFloats.push_back(Slopes);
  SegmentInfoFloats.push_back(Slopes_Lf);
  SegmentInfoFloats.push_back(Discharges);
  SegmentInfoFloats.push_back(TransportCapacities);
  SegmentInfoFloats.push_back(SedimentSupplies);
  SegmentInfoFloats.push_back(AVs);

  LSDIndexRaster SegmentsRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ChannelSegments,GeoReferencingStrings);
  ChannelSegmentsRaster = SegmentsRaster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Remove channel segments which are not downstream given a list of source
// nodes
//
// FJC 14/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::remove_tributary_segments(LSDFlowInfo& FlowInfo, vector<int> Sources, vector < vector<int> >& SegmentInfoInts, vector < vector<float> >& SegmentInfoFloats)
{
  vector < vector<int> > SegmentInfoInts_temp;
  vector < vector<float> > SegmentInfoFloats_temp;
  //vectors for storing information about the segments
  vector<int> SegmentIDs;       // ID of each segment
  vector<int> StartNodes;       // start node (upstream) of each segment
  vector<int> EndNodes;         // end node (downstream) of each segment
  vector<float> SegmentLengths; // length of each segment
  vector<float> FlowLengths;    // flow length of each segment
  vector<float> Elevations;     // elevation of each start node
  vector<float> Slopes;         // slope of each segment
  vector<float> Slopes_Lf;      // slopes based on flow length
  vector<float> Discharges;     // discharge of each segment
  vector<float> TransportCapacities; // transport capacity of each segment


  // loop through all the start nodes and check whether they are downstream of a source node
  for (int i = 0; i < int(SegmentInfoInts[0].size()); ++i)
  {
    int CurrentNode = SegmentInfoInts[1][i];
    int downstream_test = 0;
    //check whether this node is downstream of a source node
    for (int j = 0; j < int(Sources.size()); ++j)
    {
      int this_test = FlowInfo.is_node_upstream(CurrentNode, Sources[j]);
      if (this_test == 1) { downstream_test = 1; }
    }
    if (downstream_test == 1)
    {
      //cout << "This node is downstream of a source! Keeping..." << endl;
      SegmentIDs.push_back(SegmentInfoInts[0][i]);
      StartNodes.push_back(SegmentInfoInts[1][i]);
      EndNodes.push_back(SegmentInfoInts[2][i]);
      SegmentLengths.push_back(SegmentInfoFloats[0][i]);
      FlowLengths.push_back(SegmentInfoFloats[1][i]);
      Elevations.push_back(SegmentInfoFloats[2][i]);
      Slopes.push_back(SegmentInfoFloats[3][i]);
      Slopes_Lf.push_back(SegmentInfoFloats[4][i]);
      Discharges.push_back(SegmentInfoFloats[5][i]);
      TransportCapacities.push_back(SegmentInfoFloats[6][i]);
    }
  }
  // push back to master vectors
  SegmentInfoInts_temp.push_back(SegmentIDs);
  SegmentInfoInts_temp.push_back(StartNodes);
  SegmentInfoInts_temp.push_back(EndNodes);
  SegmentInfoFloats_temp.push_back(SegmentLengths);
  SegmentInfoFloats_temp.push_back(FlowLengths);
  SegmentInfoFloats_temp.push_back(Elevations);
  SegmentInfoFloats_temp.push_back(Slopes);
  SegmentInfoFloats_temp.push_back(Slopes_Lf);
  SegmentInfoFloats_temp.push_back(Discharges);
  SegmentInfoFloats_temp.push_back(TransportCapacities);

  //copy to output vecvecs
  SegmentInfoInts = SegmentInfoInts_temp;
  SegmentInfoFloats = SegmentInfoFloats_temp;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print channel segment information to csv file
//
// FJC 07/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::print_channel_segments_to_csv(LSDFlowInfo& FlowInfo, vector <vector <int> > SegmentInfoInts, vector <vector <float> > SegmentInfoFloats, string outfilename)
{
  // fid the last '.' in the filename to use in the csv filename
  unsigned dot = outfilename.find_last_of(".");

  string prefix = outfilename.substr(0,dot);
  //string suffix = str.substr(dot);
  string insert = "_sub_reaches.csv";
  string outfname = prefix+insert;
  ofstream csv_out;
  csv_out.open(outfname.c_str());
  csv_out.precision(8);

  csv_out << "SegmentID,NodeNumber,x,y,SegmentLength,FlowLength,Elevation,Slope,Slope_Lf,Qmed,Q_c,Q_s,t,CATCH_,HYDRO_CODE" << endl;
  int current_row, current_col;
  float x,y;
  cout << "N segments: " << SegmentInfoInts[0].size() << endl;

  // first print all the start nodes to csv file
  for (int i = 0; i < int(SegmentInfoInts[0].size()); i++)
  {
    int current_node = SegmentInfoInts[1][i];
    // get the row and col
    FlowInfo.retrieve_current_row_and_col(current_node, current_row, current_col);
    // get the x and y location of the node
    // the last 0.0001*DataResolution is to make sure there are no integer data points
    x = XMinimum + float(current_col)*DataResolution + 0.5*DataResolution + 0.0001*DataResolution;

    // the last 0.0001*DataResolution is to make sure there are no integer data points
    // y coord a bit different since the DEM starts from the top corner
    y = YMinimum + float(NRows-current_row)*DataResolution - 0.5*DataResolution + 0.0001*DataResolution;

    csv_out << SegmentInfoInts[0][i] << "," << current_node << "," << x << "," << y << "," << SegmentInfoFloats[0][i] << "," << SegmentInfoFloats[1][i] << "," << SegmentInfoFloats[2][i] << "," << SegmentInfoFloats[3][i] << "," << SegmentInfoFloats[4][i] << "," << SegmentInfoFloats[5][i] << "," << SegmentInfoFloats[6][i] << "," << SegmentInfoFloats[7][i] << "," << SegmentInfoFloats[8][i] << "," << SegmentInfoInts[3][i] << "," << SegmentInfoInts[4][i] << endl;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// SplitHillslopes
// This function is intended to follow the SplitChannel function.  It traces
// through the receiver nodes from every hillslope pixel and then assigns them
// an integer value that matches the index of the section of channel that is
// setting the base level of that hillslope.
//
// DTM 29/10/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::SplitHillslopes(LSDFlowInfo& FlowInfo, LSDIndexRaster& ChannelSegmentsRaster)
{
  Array2D<int> HillslopeSegmentArray(NRows,NCols,NoDataValue);
  Array2D<int> ChannelSegmentArray = ChannelSegmentsRaster.get_RasterData();
  vector<int> rows_visited,cols_visited;
  Array2D<int> VisitedBefore(NRows,NCols,0);
  int CurrentNode,ReceiverNode,ReceiverRow,ReceiverCol;
  // loop through the raster finding hillslope pixels
  for(int i = 0; i < NRows; ++i)
  {
    cout << flush << i+1 << "/" << NRows << "\r";
    for(int j = 0; j < NCols; ++j)
    {
      // Has node been visited before?
      bool VisitedBeforeTest;
      if(VisitedBefore[i][j]==1)VisitedBeforeTest = true;
      // If not visted before, then we can carry on, but mark as now visited
      else
      {
        VisitedBeforeTest = false;
        VisitedBefore[i][j]=1;
      }
      // Test that the node is a data node but not a channel node, and that it
      // hasn't been visited yet!
      if((FlowInfo.NodeIndex[i][j]!=NoDataValue) && (ChannelSegmentArray[i][j]==NoDataValue) && (VisitedBeforeTest == false))
      {
        bool finish_trace = false;
        CurrentNode = FlowInfo.NodeIndex[i][j];
        rows_visited.push_back(i);
        cols_visited.push_back(j);
        while(finish_trace == false)
        {
          FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
          // if the receiver is a stream pixel then read through the vector of
          // visited rows/columns and update hillslope segment array for each
          // using the index of the channel segment.
          if(ChannelSegmentArray[ReceiverRow][ReceiverCol] != NoDataValue)
          {
            finish_trace = true;
            int N_nodes = rows_visited.size();
            for (int i_vec = 0; i_vec < N_nodes; ++i_vec)
            {
              HillslopeSegmentArray[rows_visited[i_vec]][cols_visited[i_vec]] = ChannelSegmentArray[ReceiverRow][ReceiverCol];
            }
            rows_visited.clear();
            cols_visited.clear();
          }
          // else if the receiver is a base level node, in which case it will
          // never reach a channel -> set hillslope segment array for vector of
          // visited rows and columns as nodata.
          else if (ReceiverNode == CurrentNode)
          {
            finish_trace = true;
            int N_nodes = rows_visited.size();
            for (int i_vec = 0; i_vec < N_nodes; ++i_vec)
            {
              HillslopeSegmentArray[rows_visited[i_vec]][cols_visited[i_vec]] = NoDataValue;
            }
            rows_visited.clear();
            cols_visited.clear();
          }
          // else if the receiver has been visited before, then read through the
          // vector of visited rows/columns and update hillslope segment array
          // for each using the index of the receiver.
          else if(VisitedBeforeTest==true)
          {
            finish_trace = true;
            int N_nodes = rows_visited.size();
            for (int i_vec = 0; i_vec < N_nodes; ++i_vec)
            {
              HillslopeSegmentArray[rows_visited[i_vec]][cols_visited[i_vec]] = HillslopeSegmentArray[ReceiverRow][ReceiverCol];
            }
            rows_visited.clear();
            cols_visited.clear();
          }
          // otherwise the next pixel must be a hillslope pixel downslope that
          // has not yet been visited.  Add it to the vectors of visited points
          // and move downstream.
          else
          {
            rows_visited.push_back(ReceiverRow);
            cols_visited.push_back(ReceiverCol);
            VisitedBefore[ReceiverRow][ReceiverCol]=1;
            CurrentNode = ReceiverNode;
          }
        }
      }
    }
  }
  LSDIndexRaster HillslopeSegmentsRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, HillslopeSegmentArray,GeoReferencingStrings);
  return HillslopeSegmentsRaster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// SplitHillslopes
// This is an overloaded function doing the same as the previous version to
// segment hillslopes according to the channel index of the channel setting its
// base level.  However, this has been adapted to include an additional input
// raster - MultiThreadChannelRaster - which recognises that real channels may
// be multithreaded and/or have widths greater than or equal to one pixel.
// To be rigourous, these should be removed from analyses of hillslope
// properties.
//
// DTM 29/10/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::SplitHillslopes(LSDFlowInfo& FlowInfo, LSDIndexRaster& ChannelSegmentsRaster, LSDIndexRaster& MultiThreadChannelRaster)
{
  Array2D<int> HillslopeSegmentArray(NRows,NCols,NoDataValue);
  Array2D<int> ChannelSegmentArray = ChannelSegmentsRaster.get_RasterData();
  vector<int> rows_visited,cols_visited;
  Array2D<int> VisitedBefore(NRows,NCols,0);
  Array2D<int> MultiThreadChannelArray = MultiThreadChannelRaster.get_RasterData();
  int CurrentNode,ReceiverNode,ReceiverRow,ReceiverCol;

  // loop through the raster finding hillslope pixels
  for(int i = 0; i < NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      // Has node been visited before?
      bool VisitedBeforeTest;
      if(VisitedBefore[i][j]==1)VisitedBeforeTest = true;
      // If not visted before, then we can carry on, but mark as now visited
      else
      {
        VisitedBeforeTest = false;
        VisitedBefore[i][j]=1;
      }
      // Test that the node is a data node but not a channel node, and that it
      // hasn't been visited yet!
      if((FlowInfo.NodeIndex[i][j]!=NoDataValue) && (ChannelSegmentArray[i][j] == NoDataValue)
          && (MultiThreadChannelArray[i][j] == 0) && (VisitedBeforeTest == false))
      {
        bool finish_trace = false;
        CurrentNode = FlowInfo.NodeIndex[i][j];
        rows_visited.push_back(i);
        cols_visited.push_back(j);
        while(finish_trace == false)
        {
          FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
          if(VisitedBefore[ReceiverRow][ReceiverCol]==1)VisitedBeforeTest = true;
          // if the receiver is a stream pixel then read through the vector of
          // visited rows/columns and update hillslope segment array for each
          // using the index of the channel segment.
          if(ChannelSegmentArray[ReceiverRow][ReceiverCol] != NoDataValue)
          {
            finish_trace = true;
            int N_nodes = rows_visited.size();
            for (int i_vec = 0; i_vec < N_nodes; ++i_vec)
            {
              HillslopeSegmentArray[rows_visited[i_vec]][cols_visited[i_vec]] = ChannelSegmentArray[ReceiverRow][ReceiverCol];
            }
            rows_visited.clear();
            cols_visited.clear();
          }
          // else if the receiver is a base level node, in which case it will
          // never reach a channel -> set hillslope segment array for vector of
          // visited rows and columns as nodata.
          else if (ReceiverNode == CurrentNode)
          {
            finish_trace = true;
            int N_nodes = rows_visited.size();
            for (int i_vec = 0; i_vec < N_nodes; ++i_vec)
            {
              HillslopeSegmentArray[rows_visited[i_vec]][cols_visited[i_vec]] = NoDataValue;
            }
            rows_visited.clear();
            cols_visited.clear();
          }
          // else if the receiver has been visited before, then read through the
          // vector of visted rows/columns and update hillslope segment array
          // for each using the index of the receiver.
          else if(VisitedBeforeTest==true)
          {
            finish_trace = true;
            int N_nodes = rows_visited.size();
            for (int i_vec = 0; i_vec < N_nodes; ++i_vec)
            {
              HillslopeSegmentArray[rows_visited[i_vec]][cols_visited[i_vec]] = HillslopeSegmentArray[ReceiverRow][ReceiverCol];
            }
            rows_visited.clear();
            cols_visited.clear();
          }

          // otherwise the next pixel must be a hillslope pixel downslope that
          // has not yet been visited.  Add it to the vectors of visited points
          // and move downstream.
          else
          {
            if(MultiThreadChannelArray[ReceiverRow][ReceiverCol] != 1)
            {
              rows_visited.push_back(ReceiverRow);
              cols_visited.push_back(ReceiverCol);
              VisitedBefore[ReceiverRow][ReceiverCol]=1;
            }
            // Update CurrentNode to trace downstream
            CurrentNode = ReceiverNode;
          }
        }
      }
    }
  }
  LSDIndexRaster HillslopeSegmentsRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, HillslopeSegmentArray,GeoReferencingStrings);
  return HillslopeSegmentsRaster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the node indices of outlets of basins of a certain order
//
// IMPORTANT: the junctions always point downstream since they can have one and only
// one receiver. However, for a basin of given order, this starts just upstream of the
// confluence to the next basin order. So the basin ##INCLUDES## the channel flowing
// downstream to the penultamite node
//
// the basin_reference_number is just a reference number for printing to
// the IndexRaster
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::extract_basin_from_junction(int basin_junction, int basin_reference_number, LSDFlowInfo& FlowInfo)
{
  if (basin_junction >= int(JunctionVector.size()))
  {
    cout << "LSDJunctionNetwork::extract_basin_from_junction junction not in list" << endl;
    exit(EXIT_FAILURE);
  }

  int receiver_junc, n_nodes_in_channel,basin_outlet;
  Array2D<int> Basin(NRows,NCols,NoDataValue);
  // get the reciever junction
  receiver_junc = ReceiverVector[basin_junction];

  LSDIndexChannel StreamLinkVector = LSDIndexChannel(basin_junction, JunctionVector[basin_junction],
                                                           receiver_junc, JunctionVector[receiver_junc], FlowInfo);

  // Find final nth order channel pixel, which is the penultimate pixel in channel.
  n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  int node,row,col;

  basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  vector<int> BasinNodeVector = FlowInfo.get_upslope_nodes(basin_outlet);
  // Loop through basin to label basin pixels with basin ID
  for (int BasinIndex = 0; BasinIndex < int(BasinNodeVector.size()); ++BasinIndex)
  {
    node = BasinNodeVector[BasinIndex];
    FlowInfo.retrieve_current_row_and_col(node,row,col);
    Basin[row][col] = basin_reference_number;
  }
  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Basin,GeoReferencingStrings);
  return IR;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function converts a list of sources used to generate the initial channel network
// into a list of junction indexes of channel heads which can be used to extract hollows.
//
// SWDG 05/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Get_Channel_Head_Junctions(vector<int> Sources, LSDFlowInfo& FlowInfo){

  vector<int> Channel_Head_Junctions;
  int i;
  int j;

  for (int q = 0; q < int(Sources.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(Sources[q],i,j);
    Channel_Head_Junctions.push_back(JunctionIndexArray[i][j]);
  }

  return Channel_Head_Junctions;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts a single hollow from a given channel head junction.
//
// SWDG 05/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::extract_hollow(int CH_junction, LSDFlowInfo& FlowInfo){

  if (CH_junction >= int(JunctionVector.size()))
  {
    cout << "LSDJunctionNetwork::extract_hollow junction not in list" << endl;
    exit(EXIT_FAILURE);
  }

  // declare variables used during the extraction of the hollows
  int channel_head_row;
  int channel_head_col;
  int node_of_ch;
  int receiver_junc;
  int hollow_outlet;
  int node;
  int row;
  int col;
  Array2D<int> Hollow(NRows,NCols,NoDataValue); //final output array

  //get the coordinates of the channel head
  node_of_ch = get_Node_of_Junction(CH_junction);
  FlowInfo.retrieve_current_row_and_col(node_of_ch, channel_head_row, channel_head_col);

  // get the reciever junction
  receiver_junc = ReceiverVector[CH_junction];

  LSDIndexChannel StreamLinkVector = LSDIndexChannel(CH_junction, JunctionVector[CH_junction],
                                                           receiver_junc, JunctionVector[receiver_junc], FlowInfo);

  hollow_outlet = StreamLinkVector.get_node_in_channel(0); //get the tip of the channel
  vector<int> HollowNodeVector = FlowInfo.get_upslope_nodes(hollow_outlet);

  // Loop through basin to label basin pixels with basin ID
  for (int HollowIndex = 0; HollowIndex < int(HollowNodeVector.size()); ++HollowIndex)
  {
    node = HollowNodeVector[HollowIndex];
    FlowInfo.retrieve_current_row_and_col(node,row,col);

    Hollow[row][col] = CH_junction;

  }

  //remove channel head pixel from hollow
  Hollow[channel_head_row][channel_head_col] = NoDataValue;

  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Hollow,GeoReferencingStrings);
  return IR;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This function extracts a series of hollows from a vector of channel head junctions.
//
// SWDG 05/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::extract_hollow(vector<int> CH_junctions, LSDFlowInfo& FlowInfo){

  // declare variables used during the extraction of the hollows
  int channel_head_row;
  int channel_head_col;
  int node_of_ch;
  int receiver_junc;
  int hollow_outlet;
  int node;
  int row;
  int col;
  Array2D<int> Hollows(NRows,NCols,NoDataValue); //final output array

  for (int q = 0; q < int(CH_junctions.size()); ++q){

    if (CH_junctions[q] >= int(JunctionVector.size())){
      cout << "LSDJunctionNetwork::extract_hollow junction not in list" << endl;
      exit(EXIT_FAILURE);
    }

    //get the coordinates of the channel head
    node_of_ch = get_Node_of_Junction(CH_junctions[q]);
    FlowInfo.retrieve_current_row_and_col(node_of_ch, channel_head_row, channel_head_col);

    // get the reciever junction
    receiver_junc = ReceiverVector[CH_junctions[q]];

    LSDIndexChannel StreamLinkVector = LSDIndexChannel(CH_junctions[q], JunctionVector[CH_junctions[q]],
                                                             receiver_junc, JunctionVector[receiver_junc], FlowInfo);

    hollow_outlet = StreamLinkVector.get_node_in_channel(0); //get the tip of the channel
    vector<int> HollowNodeVector = FlowInfo.get_upslope_nodes(hollow_outlet);

    // Loop through basin to label basin pixels with basin ID
    for (int HollowIndex = 0; HollowIndex < int(HollowNodeVector.size()); ++HollowIndex){
      node = HollowNodeVector[HollowIndex];
      FlowInfo.retrieve_current_row_and_col(node,row,col);

      Hollows[row][col] = CH_junctions[q];

    }

  //remove channel head pixel from hollow
  Hollows[channel_head_row][channel_head_col] = NoDataValue;

  }

  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Hollows,GeoReferencingStrings);
  return IR;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
//This function gets the an LSDIndexRaster of basins draining from a vector of junctions.
//
// IMPORTANT: the junctions always point downstream since they can have one and only
// one receiver. However, for a basin of given order, this starts just upstream of the
// confluence to the next basin order. So the basin ##INCLUDES## the channel flowing
// downstream to the penultamite node
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::extract_basins_from_junction_vector(vector<int> basin_junctions, LSDFlowInfo& FlowInfo)
{

  Array2D<int> Basin(NRows,NCols,NoDataValue);

  for (vector<int>::iterator it = basin_junctions.begin(); it !=  basin_junctions.end(); ++it)
  {

    int basin_junction = *it;

    if (basin_junction >= int(JunctionVector.size()))
    {
      cout << "LSDJunctionNetwork::extract_basin_from_junction junction not in list" << endl;
      exit(EXIT_FAILURE);
    }

    int receiver_junc, n_nodes_in_channel, basin_outlet;

    // get the reciever junction
    receiver_junc = ReceiverVector[basin_junction];

    LSDIndexChannel StreamLinkVector = LSDIndexChannel(basin_junction, JunctionVector[basin_junction],
                                                             receiver_junc, JunctionVector[receiver_junc], FlowInfo);

    // Find final nth order channel pixel, which is the penultimate pixel
    // in channel.
    n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
    int node,row,col;

    basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
    vector<int> BasinNodeVector = FlowInfo.get_upslope_nodes(basin_outlet);
    // Loop through basin to label basin pixels with basin ID
    for (int BasinIndex = 0; BasinIndex < int(BasinNodeVector.size()); ++BasinIndex)
    {
      node = BasinNodeVector[BasinIndex];
      FlowInfo.retrieve_current_row_and_col(node,row,col);
      Basin[row][col] = basin_junction;
    }

  }

  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Basin,GeoReferencingStrings);
  return IR;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
//This function gets the an LSDIndexRaster of basins draining from a vector of junctions.
//
// IMPORTANT: the junctions always point downstream since they can have one and only
// one receiver. However, for a basin of given order, this starts just upstream of the
// confluence to the next basin order. So the basin ##INCLUDES## the channel flowing
// downstream to the penultamite node
//
// SMM 01/09/2012
//
// UPDATED so that later basins don't overwrite smaller basins. Input junction
// vector is first sorted by upslope drainage area - do the nested basins first,
// then larger basins won't overwrite these.  FJC 10/01/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::extract_basins_from_junction_vector_nested(vector<int> basin_junctions, LSDFlowInfo& FlowInfo)
{
  Array2D<int> Basin(NRows,NCols,NoDataValue);

	// sort basin junctions by contributing nodes
	vector<size_t> index_map;
	vector<int> basin_junctions_CP;
	for (int i = 0; i < int(basin_junctions.size()); i++)
	{
		// get contributing pixels of each junction
		int this_node = get_Node_of_Junction(basin_junctions[i]);
		basin_junctions_CP.push_back(FlowInfo.NContributingNodes[this_node]);
	}
	matlab_int_sort(basin_junctions_CP, basin_junctions_CP, index_map);
	matlab_int_reorder(basin_junctions, index_map, basin_junctions);

  for (vector<int>::iterator it = basin_junctions.begin(); it !=  basin_junctions.end(); ++it)
  {
		cout << "Basin junction: " << *it << endl;

    int basin_junction = *it;

    if (basin_junction >= int(JunctionVector.size()))
    {
      cout << "LSDJunctionNetwork::extract_basin_from_junction junction not in list" << endl;
      exit(EXIT_FAILURE);
    }

    int receiver_junc, n_nodes_in_channel, basin_outlet;

    // get the reciever junction
    receiver_junc = ReceiverVector[basin_junction];

    LSDIndexChannel StreamLinkVector = LSDIndexChannel(basin_junction, JunctionVector[basin_junction],
                                                             receiver_junc, JunctionVector[receiver_junc], FlowInfo);

    // Find final nth order channel pixel, which is the penultimate pixel
    // in channel.
    n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
    int node,row,col;

    basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
    vector<int> BasinNodeVector = FlowInfo.get_upslope_nodes(basin_outlet);
    // Loop through basin to label basin pixels with basin ID
    for (int BasinIndex = 0; BasinIndex < int(BasinNodeVector.size()); ++BasinIndex)
    {
      node = BasinNodeVector[BasinIndex];
      FlowInfo.retrieve_current_row_and_col(node,row,col);
			if (Basin[row][col] == NoDataValue)
			{
      	Basin[row][col] = basin_junction;
			}
    }

  }

  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Basin,GeoReferencingStrings);
  return IR;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
//This function gets basins in a rudimentary way: it just takes a list of nodes
// numbers the basins.
// Basins later in the list overwrite basins earlier in the list
//
// SMM 07/05/2015
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDJunctionNetwork::extract_basins_from_junctions_rudimentary(vector<int> junctions, LSDFlowInfo& FlowInfo)
{

  Array2D<int> Basin(NRows,NCols,NoDataValue);
  int node,row,col;
  int n_juncs = junctions.size();
  for (int i = 0; i<n_juncs; i++)
  {
    int basin_outlet = get_Node_of_Junction(junctions[i]);
    vector<int> BasinNodeVector = FlowInfo.get_upslope_nodes(basin_outlet);
    for (int BasinIndex = 0; BasinIndex < int(BasinNodeVector.size()); ++BasinIndex)
    {
      node = BasinNodeVector[BasinIndex];
      FlowInfo.retrieve_current_row_and_col(node,row,col);
      Basin[row][col] = junctions[i];
    }
  }

  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Basin,GeoReferencingStrings);
  return IR;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this sends the StreamOrderArray to a LSDIndexRaster
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::StreamOrderArray_to_LSDIndexRaster()
{
  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, StreamOrderArray,GeoReferencingStrings);
  return IR;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this sends the StreamOrderArray to a WGS84 csv
//
// SMM 12/11/2016
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::StreamOrderArray_to_WGS84CSV(string FileName_prefix)
{
  // append csv to the filename
  string FileName = FileName_prefix+".csv";

  //open a file to write
  ofstream WriteData;
  WriteData.open(FileName.c_str());

  WriteData.precision(8);
  WriteData << "latitude,longitude,Stream Order" << endl;

  // the x and y locations
  double latitude,longitude;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  //loop over each cell and if there is a value, write it to the file
  for(int i = 0; i < NRows; ++i)
  {
    for(int j = 0; j < NCols; ++j)
    {
      if (StreamOrderArray[i][j] != NoDataValue)
      {
        get_lat_and_long_locations(i, j, latitude, longitude, Converter);

        WriteData << latitude << "," << longitude << "," << StreamOrderArray[i][j] << endl;
      }
    }
  }

  WriteData.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This prints a channel network to csv in WGS84
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::PrintChannelNetworkToCSV(LSDFlowInfo& flowinfo, string fname_prefix)
{

  // first get the vectors
  vector<int> NIvec;
  vector<int> SOvec;
  vector<int> JIvec;
  GetChannelNodesAndJunctions(flowinfo, NIvec, JIvec, SOvec);

  // Deal with setting up the file
  // append csv to the filename
  string FileName = fname_prefix+".csv";

  //open a file to write
  ofstream WriteData;
  WriteData.open(FileName.c_str());

  WriteData.precision(8);
  WriteData << "latitude,longitude,Junction Index,Stream Order,NI,receiver_NI" << endl;

  // the x and y locations
  double latitude,longitude;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;


  // now get the number of channel nodes
  int this_NI, receiver_NI;
  int row,col;
  int NNodes = int(NIvec.size());
  //cout << "The number of nodes is: " << NNodes << endl;
  for(int node = 0; node<NNodes; node++)
  {
    this_NI = NIvec[node];
    flowinfo.retrieve_current_row_and_col(this_NI,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);
    flowinfo.retrieve_receiver_information(this_NI, receiver_NI);

    WriteData << latitude << "," << longitude << "," << JIvec[node] << "," << SOvec[node] << "," << NIvec[node] << "," << receiver_NI << endl;

  }

  WriteData.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This prints a channel network to csv in WGS84
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::PrintChannelNetworkToCSV_nolatlon(LSDFlowInfo& flowinfo, LSDRaster& elevation, string fname_prefix)
{

  // first get the vectors
  vector<int> NIvec;
  vector<int> SOvec;
  vector<int> JIvec;
  GetChannelNodesAndJunctions(flowinfo, NIvec, JIvec, SOvec);

  // Deal with setting up the file
  // append csv to the filename
  string FileName = fname_prefix+".csv";

  //open a file to write
  ofstream WriteData;
  WriteData.open(FileName.c_str());

  WriteData.precision(8);
  WriteData << "X,Y,Junction Index,Stream Order,NI,elevation,receiver_NI" << endl;

  // the x and y locations
  double X,Y;



  // now get the number of channel nodes
  int this_NI, receiver_NI;
  int row,col;
  int NNodes = int(NIvec.size());
  float this_elev;
  //cout << "The number of nodes is: " << NNodes << endl;
  for(int node = 0; node<NNodes; node++)
  {
    this_NI = NIvec[node];
    flowinfo.retrieve_current_row_and_col(this_NI,row,col);
    flowinfo.retrieve_receiver_information(this_NI, receiver_NI);
    this_elev = elevation.get_data_element(row,col);
    flowinfo.get_x_and_y_locations(row,col,X,Y);
    WriteData << X << "," << Y << "," << JIvec[node] << "," << SOvec[node] << "," << NIvec[node] << "," << this_elev << "," << receiver_NI << endl;

  }

  WriteData.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
// This prints a channel network to csv in WGS84
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::PrintChannelNetworkToCSV_WithElevation(LSDFlowInfo& flowinfo, string fname_prefix, LSDRaster& Elevation)
{

  // first get the vectors
  vector<int> NIvec;
  vector<int> SOvec;
  vector<int> JIvec;
  GetChannelNodesAndJunctions(flowinfo, NIvec, JIvec, SOvec);

  // Deal with setting up the file
  // append csv to the filename
  string FileName = fname_prefix+".csv";

  //open a file to write
  ofstream WriteData;
  WriteData.open(FileName.c_str());

  WriteData.precision(8);
  WriteData << "latitude,longitude,Junction Index,Stream Order,NI,receiver_NI,elevation(m)" << endl;

  // the x and y locations
  double latitude,longitude;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;


  // now get the number of channel nodes
  int this_NI, receiver_NI;
  int row,col;
  int NNodes = int(NIvec.size());
  float this_elev;
  //cout << "The number of nodes is: " << NNodes << endl;
  for(int node = 0; node<NNodes; node++)
  {
    this_NI = NIvec[node];
    flowinfo.retrieve_current_row_and_col(this_NI,row,col);
    get_lat_and_long_locations(row, col, latitude, longitude, Converter);
    flowinfo.retrieve_receiver_information(this_NI, receiver_NI);
    this_elev = Elevation.get_data_element(row,col);

    WriteData << latitude << "," << longitude << "," << JIvec[node] << "," << SOvec[node] << "," << NIvec[node] << "," << receiver_NI << "," << this_elev << endl;

  }

  WriteData.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this sends the JunctionArray to a LSDIndexRaster
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::JunctionArray_to_LSDIndexRaster()
{
  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, JunctionArray,GeoReferencingStrings);
  return IR;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this sends the JunctionIndexArray to a LSDIndexRaster
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::JunctionIndexArray_to_LSDIndexRaster()
{
  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, JunctionIndexArray,GeoReferencingStrings);
  return IR;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this turns the StreamOrderArray into a binary rastser
// where
// 1 == channel
// 0 == hillslope
//
// SMM 01/09/2012
// fixed bug where output was georef to Xmin,Xmin instead of Xmin,Ymin --SWDG 2/7/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::StreamOrderArray_to_BinaryNetwork_LSDIndexRaster()
{
  Array2D<int> BinaryNetwork(NRows,NCols,0);
  for(int row = 0; row<NRows; row++)
  {
    for (int col = 0; col<NCols; col++)
    {
      if(StreamOrderArray[row][col] == NoDataValue)
      {
        BinaryNetwork[row][col] = NoDataValue;
      }
      else if (StreamOrderArray[row][col] >= 1)
      {
        BinaryNetwork[row][col] = 1;
      }
    }
  }

  LSDIndexRaster IR(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, BinaryNetwork,GeoReferencingStrings);
  return IR;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the junction information
// it does it in binary information
// there is a seperate 'pickle' function that puts everyting into binary format
//
// SMM 01/09/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::print_junction_info_vectors(string filename)
{
  string string_filename;
  string dot = ".";
  string extension = "txt";
  string_filename = filename+dot+extension;
  cout << "The filename is " << string_filename << endl;

  // print out all the donor, reciever and stack info
  ofstream donor_info_out;
  donor_info_out.open(string_filename.c_str());
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << i << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << JunctionVector[i] << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << StreamOrderVector[i] << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << ReceiverVector[i] << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << NDonorsVector[i] << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions+1; i++)
  {
    donor_info_out << DeltaVector[i] << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << DonorStackVector[i] << " ";
  }
  donor_info_out << endl;
  for(int i = 0; i<NJunctions; i++)
  {
    donor_info_out << SVector[i] << " ";
  }
  donor_info_out << endl;

  if( int(SVectorIndex.size()) == NJunctions)
  {
    for(int i = 0; i<NJunctions; i++)
    {
      donor_info_out << SVectorIndex[i] << " ";
    }
    donor_info_out << endl;
    for(int i = 0; i<NJunctions; i++)
    {
      donor_info_out << NContributingJunctions[i] << " ";
    }
    donor_info_out << endl;
  }

  cout << "LINE 746 " << endl;
  donor_info_out.close();
  cout << "LINE 749" << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Quick and dirty way to get channels of a defined stream order, no input error handling,
// will return an LSDIndexRaster of NoDataValues if an erroneous order is passed in.
//
// Input an integer of the required stream order, returns an LSDIndexRaster of the desired
// channels. - SWDG 04/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::GetStreams(int order)
{
  Array2D<int> SingleStream(NRows,NCols,NoDataValue);

  for (int i = 0; i < NRows; ++i){
    for (int j = 0; j < NCols; ++j){
      if (StreamOrderArray[i][j] == order){
        SingleStream[i][j] = StreamOrderArray[i][j];
      }
    }
  }

  LSDIndexRaster Stream(NRows,NCols, XMinimum, YMinimum, DataResolution,
                        NoDataValue, SingleStream,GeoReferencingStrings);
	return Stream;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Overloaded version of the GetStreams function wich operates in the same way but takes a
// minimum and a maximum order integer and will return an LSDIndexRaster of the relevant
// channels. No error handling - could give odd results if strange values are passed in.
// SWDG 04/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDJunctionNetwork::GetStreams(int min_order, int max_order)
{
  Array2D<int> SelectedStreams(NRows,NCols,NoDataValue);

  for (int i = 0; i < NRows; ++i){
    for (int j = 0; j < NCols; ++j){
      if (StreamOrderArray[i][j] >= min_order && StreamOrderArray[i][j] <= max_order){
        SelectedStreams[i][j] = StreamOrderArray[i][j];
      }
    }
  }

  LSDIndexRaster Stream(NRows,NCols, XMinimum, YMinimum, DataResolution,
                        NoDataValue, SelectedStreams,GeoReferencingStrings);
  return Stream;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to test whether a junction's upstream nodes border areas of No Data
// important to ensure basins are not being artificially truncated
//
// Pass in a flowinfo object and the junction index you want to test
//
// returns a boolean indicating if no data values are present or not:
// false (0) = only good data values
// true (1) = no data values present
//
// Updated 24/10/13 to handle junction numbers in the same way that the basin extraction code does,
// by searching one junction downstream of the given junction and then going back up by one node. SWDG
//
// SWDG 27/06/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
bool LSDJunctionNetwork::node_tester(LSDFlowInfo& FlowInfo, int input_junction)
{

  Array2D<int> FlowDirection = FlowInfo.get_FlowDirection();  //used as a proxy of the elevation data
  bool flag = false;

  //get reciever junction of the input junction
  int receiver_junc = ReceiverVector[input_junction];

  // This is the node where I will check all the upstream nodes
  int basin_outlet;


  // Now see if it is the reciever junction
  if(input_junction == receiver_junc)
  {
    basin_outlet = JunctionVector[input_junction];
  }
  else   // this is a bit more complux but saves masses of computational time
  {
    //cout << "input junction:" << input_junction << " and reciever junction " << receiver_junc << endl;
    // Create channel segement from input junction down to receiver junction
    LSDIndexChannel StreamLinkVector = LSDIndexChannel(input_junction, JunctionVector[input_junction],
                                                     receiver_junc, JunctionVector[receiver_junc], FlowInfo);

    // Get the number of nodes (DEM Cells) that make up the channel segment
    int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();

    // get the penultimate node in the channel. Eg one pixel upstream from the outlet node of a basin.
    // -2 is used due to zero indexing.
    //   if(n_nodes_in_channel == 1) basin_outlet = StreamLinkVector.get_node_in_channel(0); // test for 1 pixel tributary
    //   else basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
    basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  }

  //Get all cells upslope of a junction - eg every cell of the drainage basin of interest
  vector<int> upslope_nodes = FlowInfo.get_upslope_nodes(basin_outlet);

  //loop over each cell in the basin and test for No Data values
  for(vector<int>::iterator it = upslope_nodes.begin(); it != upslope_nodes.end(); ++it)
  {
    int i;
    int j;
    FlowInfo.retrieve_current_row_and_col(*it,i,j);

    //check for edges of the file
    if (i == 0 || i == (NRows - 1) || j == 0 || j == (NCols - 1)){
    flag = true;
    return flag;}

    // check surrounding cells for NoDataValue
    else if (FlowDirection[i+1][j+1] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i+1][j] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i+1][j-1] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i][j+1] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i][j-1] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i-1][j+1] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i-1][j] == NoDataValue){flag = true;
    return flag;}
    else if (FlowDirection[i-1][j-1] == NoDataValue){flag = true;
    return flag;}
  }
  return flag;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns a integer vector containing the junction number of the largest
// donor catchment (i.e. donor junction with greatest drainage area) upslope of each
// baselevel node. These can then be used as the starting locations for performing chi
// analysis.
//
// // IMPORTANT: the junctions always point downstream since they can have one and only
// one receiver. However, for a basin of given order, this starts just upstream of the
// confluence to the next basin order. So the basin ##INCLUDES## the channel flowing
// downstream to the penultamite node
//
// MDH 19/6/13
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_BaseLevel_DonorJunctions()
{
  vector<int> BL_Donor_junctions;

  //loope through baselevel junctions and get junction immediately upstrea,
  for (int junc = 0; junc< int(BaseLevelJunctions.size()); ++junc)
  {
    if(BaseLevelJunctions[junc] < 0 || BaseLevelJunctions[junc] > NJunctions-1)
    {
      cout << " Tried LSDJunctionNetwork::get_BaseLevel_DonorJunctions but the"
           << " junction number does not exist" << endl;
      exit(0);
    }

    int SVector_junction = SVectorIndex[BaseLevelJunctions[junc]];
    BL_Donor_junctions.push_back(SVector[SVector_junction+1]);
  }
  return BL_Donor_junctions;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prunes a list of baselevel junctions by removing juntions that have
// basins bouned by NoData.
// It seeks to remove basins draining from the edge.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_Edge(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo)
{
  vector<int> BL_Donor_junctions_pruned;
  int N_BaseLevelJuncs = int(BaseLevelJunctions_Initial.size());
  cout << endl << endl << "I am going to remove any basins draining to the edge." << endl;

  for(int i = 0; i < N_BaseLevelJuncs; ++i)
  {
    //cout << "I'm checking node " << i << " to see if it is truncated." << endl;
    bool keep_base_level_node = true;

    // get donor nodes to base level nodes node
    vector<int> DonorJunctions = get_donor_nodes(BaseLevelJunctions_Initial[i]);
    int N_Donors = DonorJunctions.size();
    //cout << "It has " << N_Donors << " donor junctions." << endl;

    if (N_Donors == 1)    // this is a tiny base level node with no donors, we won't keep it.
    {
      //cout << "This is a base level junction! Junction is: " << BaseLevelJunctions_Initial[i] << endl
      //     << "and donor is: " << DonorJunctions[0] << endl;
      keep_base_level_node = false;
    }
    else
    {
      // check to see if either donor nodes are truncated - basically want to keep
      // basins that flow onto edge of DEM

      for(int i_donor = 0; i_donor < N_Donors; ++ i_donor)
      {
        if (DonorJunctions[i_donor] == BaseLevelJunctions_Initial[i])
        {
          //cout << "This is a baselevel junction." << endl;
        }
        else
        {
          bool IsTruncated = node_tester(FlowInfo,DonorJunctions[i_donor]);
          if(IsTruncated == true)
          {
            keep_base_level_node = false;
          }
        }
      }
    }
    if(keep_base_level_node == true)
    {
      BL_Donor_junctions_pruned.push_back(BaseLevelJunctions_Initial[i]);
    }
  }

  cout << "I have removed the channels that are draining from the edge of the DEM." << endl;
  cout << "I now have " << BL_Donor_junctions_pruned.size() << " base level junctions" << endl;
  return BL_Donor_junctions_pruned;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prunes a list of baselevel junctions by removing juntions that have
// basins bouned by NoData.
// It seeks to remove basins draining from the edge.
// Similar to the previous function but this one ignores the outlet reach
// because in DEMs draining to a cut edge it is common to have a nodata
// node near the outlet.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_Edge_Ignore_Outlet_Reach(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& TestRaster)
{
  vector<int> BL_Donor_junctions_pruned;
  int N_BaseLevelJuncs = int(BaseLevelJunctions_Initial.size());
  //cout << endl << endl << "I am going to remove any basins draining to the edge, ignoring the outlet reach." << endl;

  for(int i = 0; i < N_BaseLevelJuncs; ++i)
  {
    //cout << "I'm checking node " << i << " to see if it is truncated." << endl;
    bool keep_base_level_node = true;

    // get donor nodes to base level nodes node
    vector<int> DonorJunctions = get_donor_nodes(BaseLevelJunctions_Initial[i]);
    int N_Donors = int(DonorJunctions.size());
    //cout << "Looking at baselevel junction: " << BaseLevelJunctions_Initial[i] << " It has " << N_Donors << " donor junctions." << endl;

    if (N_Donors == 1)    // this is a tiny base level node with no donors, we won't keep it.
    {
      //cout << "This is a base level junction! Junction is: " << BaseLevelJunctions_Initial[i] << endl
      //     << "and donor is: " << DonorJunctions[0] << endl;
      keep_base_level_node = false;
    }
    else
    {
      // get the node indices of the donors
      for (int i_donor = 0; i_donor<N_Donors; i_donor++)
      {
        //cout << "Donor " << i_donor << " of " << N_Donors << " and junction is: " << DonorJunctions[i_donor] << endl;

        // make sure you are not getting the baselevl node since it will certainly
        // have nodata bounding it.
        if(DonorJunctions[i_donor] != BaseLevelJunctions_Initial[i])
        {
          int this_NI = JunctionVector[ DonorJunctions[i_donor] ];
          bool is_influenced_by_nodata = FlowInfo.is_upstream_influenced_by_nodata(this_NI, TestRaster);
          if (is_influenced_by_nodata)
          {
            //cout << "This node has a NoData influence upslope." << endl;
            keep_base_level_node = false;
          }
        }
      }
    }

    // Keep this baselelvel node if it is true
    if(keep_base_level_node == true)
    {
      BL_Donor_junctions_pruned.push_back(BaseLevelJunctions_Initial[i]);
    }
  }
  //cout << "I have removed the channels that are draining from the edge of the DEM." << endl;
  //cout << "I now have " << BL_Donor_junctions_pruned.size() << " base level junctions" << endl;
  return BL_Donor_junctions_pruned;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loops through all basins and takes the largest basin
// in each baselevel basin that is not influenced by nodata
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_To_Largest_Complete_Basins(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& TestRaster, LSDIndexRaster& FlowAcc)
{
  // We will loop through each basin. We get all the junctions in the basin
  // we then test if each is influenced by nodata: we then keep the largest one
  // influence by nodata
  vector<int> pruned_basin_list;
  int this_pruned_basin;
  int N_BL_Nodes = int(BaseLevelJunctions_Initial.size());
  for(int BLJ = 0; BLJ<N_BL_Nodes; BLJ++)
  {
    // get all the donor junctions. We need to have the baselevel junction
    // first so we have to append the upslope junctions after that
    vector<int> upslope_juncs;
    upslope_juncs.push_back(BaseLevelJunctions_Initial[BLJ]);
    vector<int> new_upslope_juncs = get_upslope_junctions( BaseLevelJunctions_Initial[BLJ] );
    upslope_juncs.insert(upslope_juncs.end(), new_upslope_juncs.begin(), new_upslope_juncs.end());

    // get the flow accumulation from each of these basins
    vector<int> contributing_pixels_junctions = get_contributing_pixels_from_specified_junctions(upslope_juncs,
                                                     FlowInfo, FlowAcc);

    int NDValue = -99;
    this_pruned_basin = NDValue;

    int N_total_juncs = int(upslope_juncs.size());
    int max_contributing_pixels = 0;

    // loop thyrough all the upslope junctions, testing if they are
    // influenced by the edge and how many contributing pixels they have
    for(int this_junc_index = 0; this_junc_index< N_total_juncs; this_junc_index++)
    {
      // get the current node index
      int this_NI = JunctionVector[ upslope_juncs[this_junc_index] ];
      bool is_influenced_by_nodata = FlowInfo.is_upstream_influenced_by_nodata(this_NI, TestRaster);

      // only record data if it is not influenced by nodata
      if (is_influenced_by_nodata == false)
      {
        // only record data if it is bigger than the previous biggest node
        if( contributing_pixels_junctions[this_junc_index] > max_contributing_pixels)
        {
          max_contributing_pixels =  contributing_pixels_junctions[this_junc_index];
          this_pruned_basin = upslope_juncs[this_junc_index];
        }
      }
    }

    // only keep the basin if it has a sensible junction number
    if(this_pruned_basin != NDValue)
    {
      pruned_basin_list.push_back(this_pruned_basin);
    }
  }

  return pruned_basin_list;

}







//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loops through all basins and takes the largest basin
// in each baselevel basin that is not influenced by nodata
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_By_Contributing_Pixel_Window(vector<int>& Junctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc,
                                              int lower_limit, int upper_limit)
{

  // get the flow accumulation (in pixels) from each of these basins
  vector<int> contributing_pixels_junctions = get_contributing_pixels_from_specified_junctions(Junctions_Initial,
                                                     FlowInfo, FlowAcc);

  // We will loop through each basin. We get all the junctions in the basin
  // we then test if each is influenced by nodata: we then keep the largest one
  // influence by nodata
  vector<int> pruned_basin_list;
  int this_CP;
  int N_J = int(contributing_pixels_junctions.size());

  // loop through all the basins getting the ones in the size window
  int min_basin_size = 1000000000;
  int max_basin_size = 0;
  for(int J = 0; J<N_J; J++)
  {
    this_CP = contributing_pixels_junctions[J];

    if(this_CP > max_basin_size)
    {
      max_basin_size = this_CP;
    }
    if(this_CP < min_basin_size)
    {
      min_basin_size = this_CP;
    }

    //cout << "Junction in BL list["<< J << "] has " << this_CP << " contributing pixels" << endl;

    // if the junction is within the contributing pixel window, keep it.
    if( this_CP >= lower_limit && this_CP < upper_limit)
    {
      pruned_basin_list.push_back( Junctions_Initial[J] );
    }
  }

  if(  pruned_basin_list.size() == 0)
  {
    cout << "LSDJunctionNetwork::Prune_Junctions_By_Contributing_Pixel_Window" << endl;
    cout << "I have failed to find any basins meeting your criteria. " << endl;
    cout << "The maximum basin in this DEM is: " << max_basin_size
         << " and the minimum is: " << min_basin_size << " pixels." << endl;
  }

  return pruned_basin_list;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loops through all basins and takes the largest basin
// in each baselevel basin that is not influenced by nodata
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_By_Contributing_Pixel_Window_Remove_Nested_And_Nodata(LSDFlowInfo& FlowInfo,
                                              LSDRaster& TestRaster, LSDIndexRaster& FlowAcc,
                                              int lower_limit, int upper_limit)
{

  // We use ALL the junctions for this operation
  vector<int> Junctions_Initial;
  for (int i = 0; i< NJunctions; i++)
  {
    Junctions_Initial.push_back(i);
  }

  // get the flow accumulation (in pixels) from each of these basins
  cout << "First let me prune within the contributing area window." << endl;
  cout << " I'm starting with " << Junctions_Initial.size() << " junctions." << endl;
  cout << "The lower limit is: " << lower_limit << " and upper limit "<< upper_limit << endl;
  vector<int> first_pruning = Prune_Junctions_By_Contributing_Pixel_Window(Junctions_Initial,FlowInfo,
                                              FlowAcc, lower_limit, upper_limit);
  cout << "Right, I've pruned those and have " << first_pruning.size() << " junctions left." << endl;

  // So now we need to prune the basins bounded by nodata, and prune the nested
  // basins. Which to do first? The nodata pruning is computationally expensive:
  // it requires a search of all the pixels in a basin. But the problem is
  // that if we prune by nesting we might remove a load of basins in a large
  // basin that are nested, only for that large basin to be removed later by the
  // nodata pruning. So even though it will be slow we need to prune by
  // nodata first.
  cout << "Now I am going to see if any are draining to the edge. " << endl;
  int N_total_juncs = int(first_pruning.size());
  vector<int> second_pruning;
  for(int this_junc_index = 0; this_junc_index< N_total_juncs; this_junc_index++)
  {
    // get the current node index
    int this_NI = JunctionVector[ first_pruning[this_junc_index] ];
    bool is_influenced_by_nodata = FlowInfo.is_upstream_influenced_by_nodata(this_NI, TestRaster);

    // only record data if it is not influenced by nodata
    if (is_influenced_by_nodata == false)
    {
      second_pruning.push_back( first_pruning[this_junc_index] );
    }
  }
  cout << "I now have " << second_pruning.size() << " Junctions left." << endl;

  // Now prune based on nesting
  cout << "Now I'm pruning out any nested junctions." << endl;
  vector<int> third_pruning = Prune_Junctions_If_Nested(second_pruning,FlowInfo, FlowAcc);
  cout << "Finished with pruning, I have " << third_pruning.size() << " junctions left." << endl;
  return third_pruning;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of junctions and then prunes them on the basis
// of whether they are nested. Nested junctions are eliminated.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_If_Nested(vector<int>& Junctions_Initial,
                                      LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc)
{
  // Find out how many junctions there are
  int N_Juncs = int(Junctions_Initial.size());
  //cout << "There are " << N_Juncs << " Junctions, from which I'll compare" << endl;

  // get the contributing pixels of these junctions
  // get the flow accumulation (in pixels) from each of these basins
  vector<int> contributing_pixels_junctions = get_contributing_pixels_from_specified_junctions(Junctions_Initial,
                                                     FlowInfo, FlowAcc);


  // now get all the possible two pair combinations of these junctions
  bool zero_indexed = true;   // this is just because the junction indices are numbered from zero
  int k = 2;                  // We want combinations of 2 junctions

  // the vec vec holds a vector of each possible combination of channels
  // each vector has two elements in it: the first and second channel in the comibination
  vector< vector<int> > combo_vecvec = combinations(N_Juncs, k, zero_indexed);

  //cout << "The combo vecvec has: " << combo_vecvec.size() << " elements" << endl;

  //vector<int> combo_first = combo_vecvec[0];
  //vector<int> combo_second = combo_vecvec[1];
  //cout << "Running through combinations, of which there are: " << combo_first.size() << endl;


  // get the number of combinations
  int N_combinations = int(combo_vecvec.size());

  // loop through all the combinations. There is always a bigger and smaller junction
  // in the combination. We store these in bigger smaller
  vector<int> bigger;
  vector<int> smaller;

  int first_index;
  int second_index;

  for (int i = 0; i<N_combinations; i++)
  {

    // get the indices into the junction vector of all the combinations
    first_index = combo_vecvec[i][0];
    second_index = combo_vecvec[i][1];
    if ( contributing_pixels_junctions[ first_index ] > contributing_pixels_junctions[ second_index ])
    {
      bigger.push_back(first_index);
      smaller.push_back(second_index);
    }
    else
    {
      bigger.push_back(second_index);
      smaller.push_back(first_index);
    }
  }

  map<int,int> Nested_Junctions;
  bool is_upstream;
  // now we loop through all combos, looking for  nested basins
  for (int i = 0; i<N_combinations; i++)
  {
    // Only check if the junction has not already been found to be nested
    if ( Nested_Junctions.find( smaller[i] ) == Nested_Junctions.end() )
    {

      // If it is upstream, add it to the nested map
      is_upstream = is_junction_upstream(Junctions_Initial[ bigger[i] ], Junctions_Initial[ smaller[i] ]);
      if(is_upstream)
      {
        Nested_Junctions[ smaller[i] ] = 1;
      }
    }
  }

  // Now loop thrugh all the junctions, checking to see if they are nested
  vector<int> non_nested_junctions;
  for (int i = 0; i<N_Juncs; i++)
  {
    // keep the junctions that are not in the nested junctions map
    if ( Nested_Junctions.find( i ) == Nested_Junctions.end() )
    {
      non_nested_junctions.push_back( Junctions_Initial[i]);
    }
  }
  return non_nested_junctions;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prunes a list of baselevel junctions by removing junctions whose
// contributing pixels are less than a threshold
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_Area(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc, int Threshold)
{
  vector<int> BL_Donor_junctions_pruned;
  int N_BaseLevelJuncs = int(BaseLevelJunctions_Initial.size());
  cout << endl << endl << "I am going to remove any basins smaller than " << Threshold << " pixels." << endl;
  cout << "We are starting with: " << N_BaseLevelJuncs << " juntions." << endl;

  int row,col, current_node;

  for(int i = 0; i < N_BaseLevelJuncs; ++i)
  {
    current_node = JunctionVector[BaseLevelJunctions_Initial[i]];
    FlowInfo.retrieve_current_row_and_col(current_node,row,col);

    // get the flow accumulation
    int Acc =  FlowAcc.get_data_element(row,col);

    cout << "The flow accumulation for this baselevel node is: " << Acc << endl;

    if(Acc >= Threshold)
    {
      BL_Donor_junctions_pruned.push_back(BaseLevelJunctions_Initial[i]);
    }
  }

  cout << "I have removed the channels smaller than a threshold area." << endl;
  cout << "I now have " << BL_Donor_junctions_pruned.size() << " base level junctions" << endl;
  return BL_Donor_junctions_pruned;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prunes a list of baselevel junctions by retaining ONLY the largest
//  basin in the list of junctions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_Largest(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc)
{
  vector<int> BL_Donor_junctions_pruned;
  int N_BaseLevelJuncs = int(BaseLevelJunctions_Initial.size());

  int row,col, current_node;
  if(BaseLevelJunctions_Initial.size() <= 0)
  {
    cout << "I am afraid you have no junctions in your junction list. Exiting." << endl;
    exit(0);
  }

  int largest_junc = BaseLevelJunctions_Initial[0];
  int largest_ncontrib = 0;

  for(int i = 0; i < N_BaseLevelJuncs; ++i)
  {
    current_node = JunctionVector[BaseLevelJunctions_Initial[i]];
    FlowInfo.retrieve_current_row_and_col(current_node,row,col);

    // get the flow accumulation
    int Acc =  FlowAcc.get_data_element(row,col);

    if (Acc > largest_ncontrib)
    {
      largest_ncontrib = Acc;
      largest_junc = BaseLevelJunctions_Initial[i];
    }
  }

  BL_Donor_junctions_pruned.push_back(largest_junc);

  return BL_Donor_junctions_pruned;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prunes a list of baselevel junctions by selecting only those junctions
// that are above or below (depending on the bool keep_junctions_below_threshold)
// a threshold elevation
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_Threshold_Elevation(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& Elev,
                                              float threshold_elevation, bool keep_junctions_below_threshold)
{
  vector<int> BL_Donor_junctions_pruned;
  int N_BaseLevelJuncs = int(BaseLevelJunctions_Initial.size());

  int row,col, current_node;
  if(BaseLevelJunctions_Initial.size() <= 0)
  {
    cout << "I am afraid you have no junctions in your junction list. Exiting." << endl;
    exit(0);
  }


  for(int i = 0; i < N_BaseLevelJuncs; ++i)
  {
    current_node = JunctionVector[BaseLevelJunctions_Initial[i]];
    FlowInfo.retrieve_current_row_and_col(current_node,row,col);

    float this_elevation = Elev.get_data_element(row,col);

    cout << "Junction: " <<  BaseLevelJunctions_Initial[i] << " and elevation is: " << this_elevation << endl;

    if (keep_junctions_below_threshold)
    {
      if (this_elevation <= threshold_elevation)
      {
        BL_Donor_junctions_pruned.push_back(BaseLevelJunctions_Initial[i]);
      }
    }
    else
    {
      if (this_elevation >= threshold_elevation)
      {
        BL_Donor_junctions_pruned.push_back(BaseLevelJunctions_Initial[i]);
      }
    }
  }

  return BL_Donor_junctions_pruned;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prunes a list of baselevel junctions by selecting only those junctions
// that are within an elevation window
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::Prune_Junctions_Elevation_Window(vector<int>& BaseLevelJunctions_Initial,
                                              LSDFlowInfo& FlowInfo, LSDRaster& Elev,
                                              float lower_threshold, float upper_threshold)
{
  vector<int> BL_Donor_junctions_pruned;
  int N_BaseLevelJuncs = int(BaseLevelJunctions_Initial.size());

  int row,col, current_node;
  if(BaseLevelJunctions_Initial.size() <= 0)
  {
    cout << "I am afraid you have no junctions in your junction list. Exiting." << endl;
    exit(0);
  }


  for(int i = 0; i < N_BaseLevelJuncs; ++i)
  {
    current_node = JunctionVector[BaseLevelJunctions_Initial[i]];
    FlowInfo.retrieve_current_row_and_col(current_node,row,col);

    float this_elevation = Elev.get_data_element(row,col);

    cout << "Junction: " <<  BaseLevelJunctions_Initial[i] << " and elevation is: " << this_elevation;

    if (this_elevation >= lower_threshold && this_elevation <= upper_threshold)
    {
      cout << " KEEPING this one." << endl;
      BL_Donor_junctions_pruned.push_back(BaseLevelJunctions_Initial[i]);
    }
    else
    {
      cout << " NOT keeping that one. " << endl;
    }
  }

  return BL_Donor_junctions_pruned;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function returns a vector of the contributing pixels from a list
// of junctions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_contributing_pixels_from_specified_junctions(vector<int>& JunctionList,
                                              LSDFlowInfo& FlowInfo, LSDIndexRaster& FlowAcc)
{
  int N_Juncs = int(JunctionList.size());
  int row,col, current_node;
  vector<int> N_contributing_pixels;
  if(JunctionList.size() <= 0)
  {
    cout << "I am afraid you have no junctions in your junction list. Exiting." << endl;
  }
  else
  {
    for(int i = 0; i < N_Juncs; ++i)
    {
      current_node = JunctionVector[JunctionList[i]];
      FlowInfo.retrieve_current_row_and_col(current_node,row,col);

      // get the flow accumulation
      N_contributing_pixels.push_back(FlowAcc.get_data_element(row,col));
    }
  }

  return N_contributing_pixels;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// GET RECEIVER JUNCTION FOR SPECIFIED COORDINATES
// For input coordinates (e.g. the location of a cosmogenic radionuclide sample), get the
// closest downslope node of the catchment.  This enables easy extraction of catchment
// for analysis.
// DTM 17/10/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_receiver_junction_for_specified_coordinates(float X_coordinate,
                                  float Y_coordinate, LSDFlowInfo& FlowInfo)
{
  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

  // Find first downstream junction by running through receiver nodes until you
  // find a junction.
  int at_junction = 0;
  int CurrentNode = FlowInfo.NodeIndex[row_point][col_point];
  int ReceiverRow, ReceiverCol, ReceiverNode, junction;
  while(at_junction<1)
  {
    FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
    CurrentNode = ReceiverNode;
    junction = retrieve_junction_number_at_row_and_column(ReceiverRow,ReceiverCol);
    //test to see if receiver node is in channel
    if(junction != NoDataValue) ++at_junction;
  }
  return junction;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Get node index of nearest channel for specified coordinates
// for lat long
// FJC 20/11/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_junction_of_nearest_channel_from_lat_long(double latitude, double longitude, LSDFlowInfo& FlowInfo, LSDCoordinateConverterLLandUTM Converter)
{
  double this_Northing, this_Easting;
  int UTM_zone;
  bool is_North;
  FlowInfo.get_UTM_information(UTM_zone, is_North);
  int eId = 22;             // defines the ellipsiod. This is WGS
  Converter.LLtoUTM_ForceZone(eId, latitude, longitude,
                    this_Northing, this_Easting, UTM_zone);

  // now get the nearest channel node
  int this_junc = get_receiver_junction_for_specified_coordinates(this_Easting, this_Northing, FlowInfo);

  return this_junc;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Get upstream junction from lat long
// FJC 27/09/18
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_upstream_junction_from_lat_long(double latitude, double longitude, LSDFlowInfo& FlowInfo, LSDCoordinateConverterLLandUTM Converter)
{
  double this_Northing, this_Easting;
  int UTM_zone;
  bool is_North;
  FlowInfo.get_UTM_information(UTM_zone, is_North);
  int eId = 22;             // defines the ellipsiod. This is WGS
  Converter.LLtoUTM_ForceZone(eId, latitude, longitude,
                    this_Northing, this_Easting, UTM_zone);

  int this_node = get_nodeindex_of_nearest_channel_for_specified_coordinates(this_Easting, this_Northing, 1, 5, FlowInfo);
  int this_junc = find_upstream_junction_from_channel_nodeindex(this_node, FlowInfo);

  return this_junc;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// GET node index of nearest channel for specificed coordinates
// For input coordinates (e.g. the location of a cosmogenic radionuclide sample), get the
// node index of the nearest channel.  This enables easy extraction of catchment
// for analysis.
//
// The X_coordinate and Y_coordinate should be in the
// same spatial reference as the DEM, typically in UTM
//
// The threshold stream order is the stream order that qualifies as a 'channel'
//
// SMM 21/10/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDJunctionNetwork::get_nodeindex_of_nearest_channel_for_specified_coordinates(float X_coordinate,
             float Y_coordinate, int search_radius_nodes, int threshold_stream_order, LSDFlowInfo& FlowInfo)
{

  // variables neighbor search
  int kernal_size = search_radius_nodes*2+1;
  int this_krow, this_kcol;
  int largest_SO_in_kernal;
  int this_SO;
  int largest_SO_row, largest_SO_col;

  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum - DataResolution*0.5;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum - DataResolution*0.5;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));
  int CurrentNode = FlowInfo.NodeIndex[row_point][col_point];

  bool is_in_raster = true;
  int NearestChannel;

  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1 || CurrentNode == NoDataValue)
  {
    is_in_raster = false;
  }

  if (is_in_raster)
  {
    // Find first downstream junction by running through receiver nodes until you
    // find a junction.
    int CurrentNode = FlowInfo.NodeIndex[row_point][col_point];
    int ReceiverRow, ReceiverCol, ReceiverNode, CurrentCol, CurrentRow;
    cout << "Current node: " << CurrentNode << endl;

    // get the current row and column
    FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);

    cout << "line 5937" << endl;
    // get the first receiver
    FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);

    cout << "CurrentNode: " << CurrentNode << " ReceiverNode: " << ReceiverNode << endl;
    // make sure you are not at base level
    NearestChannel = NoDataValue;

    // check to see if this node has a stream order >= 1
    if(StreamOrderArray[CurrentRow][CurrentCol] >= threshold_stream_order)
    {
      NearestChannel = CurrentNode;
    }

    // loop until you find a channel
    while(NearestChannel == NoDataValue && CurrentNode != ReceiverNode)
    {
      //cout << "Found a channel" << endl;
      // now move down one node
      CurrentNode = ReceiverNode;

      // get the current row and column
      FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);

      // now search the kernal
      largest_SO_in_kernal = NoDataValue;
      largest_SO_row = NoDataValue;
      largest_SO_col = NoDataValue;
      for (int krow = 0; krow<kernal_size; krow++)
      {
        for (int kcol = 0; kcol<kernal_size; kcol++)
        {
          // get the row and column
          this_krow = CurrentRow-search_radius_nodes+krow;
          this_kcol = CurrentCol-search_radius_nodes+kcol;

          // only test if it within size of the Stream Order array
          if(this_krow >= 0 && this_krow < NRows-1 && this_kcol >= 0 && this_kcol < NCols-1)
          {
            this_SO = StreamOrderArray[this_krow][this_kcol];
            if (this_SO >= threshold_stream_order && this_SO > largest_SO_in_kernal)
            {
              largest_SO_in_kernal = this_SO;
              largest_SO_row = this_krow;
              largest_SO_col = this_kcol;
            }
          }
        }
      }


      // check to if the kernal returned a channel node
      if(largest_SO_in_kernal != NoDataValue)
      {
        NearestChannel = FlowInfo.NodeIndex[largest_SO_row][largest_SO_col];
      }
      else    // get the next node
      {
        FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
        //cout << "CurrentNode: " << CurrentNode << " RN: " << ReceiverNode << endl;
      }
    }
  }
  else
  {
    NearestChannel = NoDataValue;
  }

  //cout << "LSDJunctionNetwork 5789 Nearest_channel node is: " << NearestChannel << endl;
  return NearestChannel;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function just flows downslope until it finds a channel of stream order greater
// than the threshold stream order
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::find_nearest_downslope_channel(float X_coordinate, float Y_coordinate, int threshold_SO, LSDFlowInfo& FlowInfo)
{
  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum - DataResolution*0.5;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum - DataResolution*0.5;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));
  int CurrentNode = FlowInfo.NodeIndex[row_point][col_point];
  int channel_node = find_nearest_downslope_channel(CurrentNode,threshold_SO, FlowInfo);
  return channel_node;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function just flows downslope until it finds a channel of stream order greater
// than the threshold stream order
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::find_nearest_downslope_channel(int StartingNode, int threshold_SO, LSDFlowInfo& FlowInfo)
{
	int ChannelNode = -9999;
  int row, col;
	int CurrentNode = StartingNode;
	int BaseLevel = FlowInfo.is_node_base_level(CurrentNode);
	FlowInfo.retrieve_current_row_and_col(StartingNode, row, col);
	//check if you are already at a channel
	if (StreamOrderArray[row][col] != NoDataValue && StreamOrderArray[row][col] >= threshold_SO
	&& BaseLevel == 0)
	{
		ChannelNode = FlowInfo.NodeIndex[row][col];
	}
	//if not at a channel, move downstream
	else
	{
		bool ReachedChannel = false;
		while (ReachedChannel == false)
		{
			//get receiver information
			int ReceiverNode, ReceiverRow, ReceiverCol;
			FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
			//if node is at baselevel then exit
			if (CurrentNode == ReceiverNode)
			{
				ChannelNode = FlowInfo.NodeIndex[ReceiverRow][ReceiverCol];
				ReachedChannel = true;
				//cout << "You reached a baselevel node, returning baselevel" << endl;
			}
			//if receiver is a channel > threshold then get the stream order
			if (StreamOrderArray[ReceiverRow][ReceiverCol] != NoDataValue &&
			StreamOrderArray[ReceiverRow][ReceiverCol] >= threshold_SO)
			{
				ChannelNode = FlowInfo.NodeIndex[ReceiverRow][ReceiverCol];
				ReachedChannel = true;
			}
			else
			{
				//move downstream
				CurrentNode = ReceiverNode;
			}
		}
	}
  return ChannelNode;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns information on the nearest channel pixel to a node
//
// FJC 08/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_info_nearest_channel_to_node(int& StartingNode, int& threshold_SO, LSDFlowInfo& FlowInfo, LSDRaster& DistFromOutlet, int& ChannelNode, float& FlowLength, float& DistanceUpstream)
{
	int row, col;
	FlowLength = 0;
	float root_2 = 1.4142135623;
	int CurrentNode = StartingNode;
	int BaseLevel = FlowInfo.is_node_base_level(CurrentNode);
	FlowInfo.retrieve_current_row_and_col(StartingNode, row, col);
	//check if you are already at a channel
	if (StreamOrderArray[row][col] != NoDataValue && StreamOrderArray[row][col] >= threshold_SO
	&& BaseLevel == 0)
	{
		ChannelNode = FlowInfo.NodeIndex[row][col];
		//cout << "You are already at a channel" << endl;
		// get the upstream distance
		DistanceUpstream = DistFromOutlet.get_data_element(row,col);
	}
	//if not at a channel, move downstream
	else
	{
		bool ReachedChannel = false;
		while (ReachedChannel == false)
		{
			//get receiver information
			int ReceiverNode, ReceiverRow, ReceiverCol;
			FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
			//if node is at baselevel then exit
			if (CurrentNode == ReceiverNode)
			{
				ChannelNode = FlowInfo.NodeIndex[ReceiverRow][ReceiverCol];
				ReachedChannel = true;
				//cout << "You reached a baselevel node, returning baselevel" << endl;
			}
			//if receiver is a channel > threshold then get the stream order
			if (StreamOrderArray[ReceiverRow][ReceiverCol] != NoDataValue &&
			StreamOrderArray[ReceiverRow][ReceiverCol] >= threshold_SO)
			{
				ChannelNode = FlowInfo.NodeIndex[ReceiverRow][ReceiverCol];
				// get the upstream distance of the nearest channel node
				DistanceUpstream = DistFromOutlet.get_data_element(ReceiverRow,ReceiverCol);
				//cout << "You've reached a channel!" << endl;
				ReachedChannel = true;
			}
			else
			{
				//move downstream
				CurrentNode = ReceiverNode;
				// update length
				if (FlowInfo.retrieve_flow_length_code_of_node(ReceiverNode) == 1){ FlowLength += DataResolution; }
				else if (FlowInfo.retrieve_flow_length_code_of_node(ReceiverNode) == 2){ FlowLength += (DataResolution * root_2); }
			}
		}
	}
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns information on the nearest channel pixel to a node. It
// checks to see if the channel pixel is in the main stem channel - if not, it
// keeps moving downstream.
//
// FJC 08/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_info_nearest_channel_to_node_main_stem(int& StartingNode, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, LSDIndexChannel& MainStem, int& ChannelNode, float& FlowLength, float& DistanceUpstream, float& Relief)
{
	int row, col;
	FlowLength = 0;
	float root_2 = 1.4142135623;
	int CurrentNode = StartingNode;
	FlowInfo.retrieve_current_row_and_col(StartingNode, row, col);

  // get the nodes in the main stem
	vector<int> main_stem_nodes = MainStem.get_NodeSequence();
	vector<int>::iterator find_it;
	find_it = find(main_stem_nodes.begin(), main_stem_nodes.end(), StartingNode);
	//check if you are already at a channel
	if (find_it != main_stem_nodes.end())
	{
		ChannelNode = StartingNode;
		//cout << "You are already at a channel" << endl;
		// get the upstream distance
		DistanceUpstream = DistFromOutlet.get_data_element(row,col);
		FlowLength=0;
		Relief=0;
	}
	//if not at a channel, move downstream
	else
	{
		bool ReachedChannel = false;
		while (ReachedChannel == false)
		{
			//get receiver information
			int ReceiverNode, ReceiverRow, ReceiverCol;
			FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
			//if node is at baselevel then exit
			if (CurrentNode == ReceiverNode)
			{
				ChannelNode = FlowInfo.NodeIndex[ReceiverRow][ReceiverCol];
				Relief = ElevationRaster.get_data_element(row,col) - ElevationRaster.get_data_element(ReceiverRow,ReceiverCol);
				ReachedChannel = true;
				//cout << "You reached a baselevel node, returning baselevel" << endl;
			}
			//if receiver is a channel node in the main stem
			vector<int>::iterator find_main_stem;
			find_main_stem = find(main_stem_nodes.begin(), main_stem_nodes.end(), ReceiverNode);
			if (find_main_stem != main_stem_nodes.end())
			{
				ChannelNode = ReceiverNode;
				// get the upstream distance of the nearest channel node
				DistanceUpstream = DistFromOutlet.get_data_element(ReceiverRow,ReceiverCol);
				Relief = ElevationRaster.get_data_element(row,col) - ElevationRaster.get_data_element(ReceiverRow,ReceiverCol);
				//cout << "You've reached a channel!" << endl;
				ReachedChannel = true;
			}
			else
			{
				//move downstream
				CurrentNode = ReceiverNode;
				// update length
				if (FlowInfo.retrieve_flow_length_code_of_node(ReceiverNode) == 1){ FlowLength += DataResolution; }
				else if (FlowInfo.retrieve_flow_length_code_of_node(ReceiverNode) == 2){ FlowLength += (DataResolution * root_2); }
			}
		}
	}

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Find upstream junction from channel nodeindex
//
// This function takes a nodeindex, checks to see if it is a channel, and if so
// it marches upstream until it finds a junction
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::find_upstream_junction_from_channel_nodeindex(int ChannelNodeIndex,
                                                   LSDFlowInfo& FlowInfo)
{
  int UpstreamJunction = NoDataValue;
  int CurrentNode = ChannelNodeIndex;
  int CurrentRow, CurrentCol;
  //int DonorRow, DonorCol;
  int junction;
  int donor_channel_order;

  // get the current row and column
  FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);
  if(CurrentRow == NoDataValue || CurrentCol == NoDataValue)
  {
    cout << "LSDJunctionNetwork::find_upstream_junction_from_channel_nodeindex::FATALERROR" << endl;
    cout << "Cannot find node on raster" << endl;
    exit(EXIT_FAILURE);
  }

  // get the stream order
  int this_channel_order = StreamOrderArray[CurrentRow][CurrentCol];

  if (this_channel_order != NoDataValue)
  {
    // first test if this is a junction
    junction = retrieve_junction_number_at_row_and_column(CurrentRow,CurrentCol);
    if(junction != NoDataValue)
    {
      UpstreamJunction = junction;
    }



    while (UpstreamJunction == NoDataValue)
    {
      // get the donors
      vector<int> donors = FlowInfo.get_donor_nodes(CurrentNode);

      // now loop through the donors looking for a channel of the same order
      int n_donors = donors.size();
      int this_donor = 0;
      do
      {
        //cout << "this donor: " << this_donor << " and the donor NI: " << donors[this_donor] << endl;
        FlowInfo.retrieve_current_row_and_col(donors[this_donor],CurrentRow,CurrentCol);
        donor_channel_order = StreamOrderArray[CurrentRow][CurrentCol];
        //cout << "donor_channel_order: " << donor_channel_order << " and tcho: " << this_channel_order << endl;

        this_donor++;
      } while( donor_channel_order != this_channel_order && this_donor<n_donors);

      // now check if the donor is a junction
      junction = retrieve_junction_number_at_row_and_column(CurrentRow,CurrentCol);
      //cout << "Junction is: " << junction << endl;

      if(junction != NoDataValue)		// if it is, set the junction
      {
        UpstreamJunction = junction;
      }
      else      // if it isn't, go upslope
      {
        CurrentNode = donors[this_donor-1];
      }

      //cout << "Current node yo1: " << CurrentNode << endl;
    }
  }

  return UpstreamJunction;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function checks the upslope nodes of a junction to test if any of them
// are the same stream order as the junction
// FJC and MAH 18/03/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::check_stream_order_of_upstream_nodes(int junction, LSDFlowInfo& FlowInfo)
{
  int same_SO = 0;
  //get the node of the current junction
  int CurrentNode = get_Node_of_Junction(junction);
  int CurrentRow, CurrentCol, UpstreamRow, UpstreamCol;
  // get the current row and column
  FlowInfo.retrieve_current_row_and_col(CurrentNode,CurrentRow,CurrentCol);

  // get the stream order
  int CurrentSO = StreamOrderArray[CurrentRow][CurrentCol];
  //cout << "Current SO: " << CurrentSO << endl;

  //loop through all the donor nodes and check the stream order
  vector<int> donors = FlowInfo.get_donor_nodes(CurrentNode);
  for(int i = 0; i < int(donors.size()); i++)
  {
    // get the upstream row and column
    FlowInfo.retrieve_current_row_and_col(donors[i],UpstreamRow,UpstreamCol);
    // get the stream order
    int UpstreamSO = StreamOrderArray[UpstreamRow][UpstreamCol];
    //cout << "Upstream stream order: " << UpstreamSO << endl;
    if(UpstreamSO == CurrentSO)
    {
      same_SO = 1;
    }
  }
  return same_SO;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the node index of the node upstream of a given node
// with the highest stream order
// FJC
// 31/01/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDJunctionNetwork::get_upstream_node_max_stream_order(int current_node, LSDFlowInfo& FlowInfo)
{
  vector<int> DonorNodes = FlowInfo.get_donor_nodes(current_node);
  int max_SO = 0;
  int max_NI = NoDataValue;
  for (int i = 0; i < int(DonorNodes.size()); i++)
  {
    int this_SO = get_StreamOrder_of_Node(FlowInfo, DonorNodes[i]);
    if (this_SO > max_SO)
    {
      max_SO = this_SO;
      max_NI = DonorNodes[i];
    }
  }
  if (max_NI == NoDataValue)
  {
    cout << "Couldn't find a donor node with a valid stream order, returning NoDataValue" << endl;
  }

  return max_NI;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a CRN file and returns the node index and the junction
// index of the nearest channels
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::snap_point_locations_to_channels(vector<float> x_locs,
                vector<float> y_locs,
                int search_radius_nodes, int threshold_stream_order,
                LSDFlowInfo& FlowInfo, vector<int>& valid_cosmo_points,
                vector<int>& snapped_node_indices, vector<int>& snapped_junction_indices)
{
  // First reset the vectors that will be copied
  vector<int> empty_vec;
  valid_cosmo_points = empty_vec;
  snapped_node_indices = empty_vec;
  snapped_junction_indices = empty_vec;


  vector<int> reciever_junc;
  int this_junc, this_chan_node;

  float x_loc,y_loc;

  //cout << "JN LINE 4287, XMinimum is: " << XMinimum << " YMinimum is: "
  //     << YMinimum << endl;

  // now loop through cosmo points recording the junctions
  int n_cosmo_points = int(x_locs.size());
  for (int samp = 0; samp<n_cosmo_points; samp++)
  //for (int samp = 0; samp<1; samp++)
  {
    x_loc = x_locs[samp];
    y_loc = y_locs[samp];

    cout << "JN 6322 x_loc: " << x_loc << " y_loc: " << y_loc << endl;

    // check to see if the node is in the raster
    bool is_in_raster = FlowInfo.check_if_point_is_in_raster(x_loc,y_loc);

    // check that the point is not nodata and set is_in_raster to false if it is
    int tmpNode = FlowInfo.get_node_index_of_coordinate_point(x_loc,y_loc);
    int tmpRow;
    int tmpCol;
    //cout << "Node Index: " << tmpNode << endl;
    if (tmpNode != NoDataValue)
    {

      FlowInfo.retrieve_current_row_and_col(tmpNode, tmpRow, tmpCol);

      if (FlowInfo.get_LocalFlowDirection(tmpRow, tmpCol) == NoDataValue){
        is_in_raster = false;
      }
    }
    else { is_in_raster = false; }

    if(is_in_raster)
    {
      cout << "Snapping: This point is in the raster!" << endl;
      this_chan_node = get_nodeindex_of_nearest_channel_for_specified_coordinates(x_loc, y_loc,
                       search_radius_nodes, threshold_stream_order,
                       FlowInfo);
      cout << "Snapping: Got channel!, channel node is: " << this_chan_node << endl;
      if(this_chan_node != NoDataValue)
      {
        this_junc = find_upstream_junction_from_channel_nodeindex(this_chan_node, FlowInfo);
        cout << "Snapping, got_this_junc with the node index " << this_chan_node <<  endl;
        snapped_node_indices.push_back(this_chan_node);
        snapped_junction_indices.push_back(this_junc);
        valid_cosmo_points.push_back(samp);
      }
      else
      {
        cout << endl << "+++" << endl;
        cout << "WARNING LSDJunctionNetwork::snap_point_locations_to_channels." << endl;
        cout << "This node is in the DEM but I have not found a nearby channel." << endl;
        cout << "+++" << endl << endl;
      }

      //cout << "channel node index is: " << this_chan_node << " and receiver_junc is: "
      //     << this_junc << endl;
    }
    else
    {
      cout << "WARNING LSDJunctionNetwork::snap_point_locations_to_channels." << endl;
      cout << "This point at location: " << x_loc << " " << y_loc << endl
           << "does not seem to be in the raster, or is in a NoData region." << endl;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a CRN file and returns the node index
// of the nearest channels.
// It does not move to nearest junction index
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::snap_point_locations_to_nearest_channel_node_index(vector<float> x_locs,
                vector<float> y_locs,
                int search_radius_nodes, int threshold_stream_order,
                LSDFlowInfo& FlowInfo, vector<int>& valid_cosmo_points,
                vector<int>& snapped_node_indices)
{
  // First reset the vectors that will be copied
  vector<int> empty_vec;
  valid_cosmo_points = empty_vec;
  snapped_node_indices = empty_vec;

  vector<int> reciever_junc;
  int this_junc, this_chan_node;

  float x_loc,y_loc;

  //cout << "JN LINE 4287, XMinimum is: " << XMinimum << " YMinimum is: "
  //     << YMinimum << endl;

  // now loop through cosmo points recording the junctions
  int n_cosmo_points = int(x_locs.size());
  for (int samp = 0; samp<n_cosmo_points; samp++)
  //for (int samp = 0; samp<1; samp++)
  {
    x_loc = x_locs[samp];
    y_loc = y_locs[samp];

    cout << "Snapping to nearest channel node index. x_loc: " << x_loc << " y_loc: " << y_loc << endl;

    // check to see if the node is in the raster
    bool is_in_raster = FlowInfo.check_if_point_is_in_raster(x_loc,y_loc);

    // check that the point is not nodata and set is_in_raster to false if it is
    int tmpNode = FlowInfo.get_node_index_of_coordinate_point(x_loc,y_loc);
    int tmpRow;
    int tmpCol;
    //cout << "Node Index: " << tmpNode << endl;
    if (tmpNode != NoDataValue)
    {

      FlowInfo.retrieve_current_row_and_col(tmpNode, tmpRow, tmpCol);

      if (FlowInfo.get_LocalFlowDirection(tmpRow, tmpCol) == NoDataValue){
        is_in_raster = false;
      }
    }
    else { is_in_raster = false; }

    if(is_in_raster)
    {
      cout << "Snapping: This point is in the raster!" << endl;
      this_chan_node = find_nearest_downslope_channel(x_loc, y_loc,threshold_stream_order,FlowInfo);
      cout << "Snapping: Got channel!, channel node is: " << this_chan_node << endl;
      if(this_chan_node != NoDataValue)
      {
        cout << "Snapping, got the node index!" << endl;
        snapped_node_indices.push_back(this_chan_node);
        valid_cosmo_points.push_back(samp);
      }
      else
      {
        cout << endl << "+++" << endl;
        cout << "WARNING LSDJunctionNetwork::snap_point_locations_to_channels." << endl;
        cout << "This node is in the DEM but I have not found a nearby channel." << endl;
        cout << "+++" << endl << endl;
      }

      //cout << "channel node index is: " << this_chan_node << " and receiver_junc is: "
      //     << this_junc << endl;
    }
    else
    {
      cout << "WARNING LSDJunctionNetwork::snap_point_locations_to_channels." << endl;
      cout << "This point at location: " << x_loc << " " << y_loc << endl
           << "does not seem to be in the raster, or is in a NoData region." << endl;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



vector<int> LSDJunctionNetwork::snap_point_locations_to_upstream_junctions_from_latlong_csv(string csv_filename,
                int search_radius_nodes, int threshold_stream_order,
                LSDFlowInfo& FlowInfo, LSDRasterInfo& RI)
{
  // read the csv file 
  LSDSpatialCSVReader outlets_data( RI, csv_filename );

  // Now get the lat long points
  // First get the lat-long
  int UTM_zone;
  bool is_North;
  outlets_data.get_UTM_information(UTM_zone, is_North);
  cout << "The UTM zone is: " << UTM_zone << " and it ";
  if(is_North)
  {
    cout << "is north." << endl;
  }
  else
  {
    cout << "is south." << endl;
  }

  // Get the local coordinates
  vector<float> fUTM_easting,fUTM_northing;
  outlets_data.get_x_and_y_from_latlong(fUTM_easting,fUTM_northing);
  cout << "The x and y data are: " << endl;
  for (int i = 0; i< int(fUTM_easting.size()); i++)
  {
    cout << fUTM_easting[i] << "," << fUTM_northing[i] << endl;
  }

  // now loop through cosmo points recording the junctions
  int n_outlet_points = int(fUTM_easting.size());
  vector<int> snapped_junction_indices;
  float x_loc,y_loc;
  int this_chan_node;
  for (int samp = 0; samp<n_outlet_points; samp++)
  //for (int samp = 0; samp<1; samp++)
  {
    x_loc = fUTM_easting[samp];
    y_loc = fUTM_northing[samp];

    cout << "Snapping to nearest channel node index. x_loc: " << x_loc << " y_loc: " << y_loc << endl;

    // check to see if the node is in the raster
    bool is_in_raster = FlowInfo.check_if_point_is_in_raster(x_loc,y_loc);

    // check that the point is not nodata and set is_in_raster to false if it is
    int tmpNode = FlowInfo.get_node_index_of_coordinate_point(x_loc,y_loc);
    int tmpRow;
    int tmpCol;
    //cout << "Node Index: " << tmpNode << endl;
    if (tmpNode != NoDataValue)
    {

      FlowInfo.retrieve_current_row_and_col(tmpNode, tmpRow, tmpCol);

      if (FlowInfo.get_LocalFlowDirection(tmpRow, tmpCol) == NoDataValue){
        is_in_raster = false;
      }
    }
    else { is_in_raster = false; }

    if(is_in_raster)
    {
      cout << "Snapping: This point is in the raster!" << endl;
      this_chan_node = find_nearest_downslope_channel(x_loc, y_loc,threshold_stream_order,FlowInfo);
      cout << "Snapping: Got channel!, channel node is: " << this_chan_node << endl;
      if(this_chan_node != NoDataValue)
      {
        cout << "Snapping, got the node index and the junction" << endl;       
        // now get the upslope junction
        int temp_junc = find_upstream_junction_from_channel_nodeindex(this_chan_node, FlowInfo);
        snapped_junction_indices.push_back(temp_junc);
      }
      else
      {
        cout << endl << "+++" << endl;
        cout << "WARNING LSDJunctionNetwork::snap_point_locations_to_channels." << endl;
        cout << "This node is in the DEM but I have not found a nearby channel." << endl;
        cout << "+++" << endl << endl;
      }

    }
    else
    {
      cout << "WARNING LSDJunctionNetwork::snap_point_locations_to_channels." << endl;
      cout << "This point at location: " << x_loc << " " << y_loc << endl
           << "does not seem to be in the raster, or is in a NoData region." << endl;
    }



  }  

  return snapped_junction_indices;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// THESE FUNCTIONS ARE WRITTEN TO COUPLE EACH HILLSLOPE PIXEL TO THE CHANNEL NODE THAT
// SETS THEIR LOWER BOUNDARY CONDITION.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// couple_hillslope_nodes_to_channel_nodes  => Should be in LSDFlowInfo
//----------------------------------------------------------------------------------------
// This function couples all hillslope pixels within a given basin to the channel node
// that sets the baselevel for that hillslope pixel.  The node on the channel network for
// which this occurs is determined using Stuart's rather wonderful hillslope flow routing
// method (REFERENCE TO GO HERE), which exploits the D-Infinity flow routing algorithm.
void LSDJunctionNetwork::couple_hillslope_nodes_to_channel_nodes(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& D_inf_Flowdir, LSDIndexRaster& ChannelNodeNetwork, int OutletJunction, vector<int>& hillslope_nodes, vector<int>& baselevel_channel_nodes)
{
  LSDIndexChannel StreamLinkVector = LSDIndexChannel(OutletJunction, JunctionVector[OutletJunction],ReceiverVector[OutletJunction],
                                                     JunctionVector[ReceiverVector[OutletJunction]], FlowInfo);
  // Step 1: get basin
  int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  int basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  vector<int> BasinNodes = FlowInfo.get_upslope_nodes(basin_outlet);
  int N_BasinNodes = BasinNodes.size();
  //-----------------------------------------//
  // option to clip all rasters to basin here//
  //-----------------------------------------//
  // Step 2: sort basin nodes by elevation
  vector<float> ElevationValues(N_BasinNodes,float(NoDataValue));
  vector<size_t> index_map;
  for(int i = 0; i<N_BasinNodes; ++i)
  {
    int row,col;
    FlowInfo.retrieve_current_row_and_col(BasinNodes[i],row,col);
    ElevationValues[i]=Elevation.get_data_element(row,col);
  }
  cout << "\t\t sorting values" << endl;
  matlab_float_sort_descending(ElevationValues,ElevationValues,index_map);
  matlab_int_reorder(BasinNodes,index_map,BasinNodes);
  // Step 3: For each node, route flow to channel
  bool skip_trace,hillslope_node_test;
  vector< vector<float> > vv_temp;
  vector<float> v_temp;
  int output_channel_node;
  vector<int> ChannelNodes,HillslopeNodes;
  int target_row,target_col;
  int channel_node_count = 0;
  cout << "\t\t routing flow from each pixel" << endl;
  for(int i = 0; i<N_BasinNodes; ++i)
  {
    cout << flush << i+1 << "/" << N_BasinNodes << "\r";
    FlowInfo.retrieve_current_row_and_col(BasinNodes[i],target_row,target_col);
    if(ChannelNodeNetwork.get_data_element(target_row,target_col) == NoDataValue)
    {
      hillslope_node_test = true;
      FlowInfo.D_Inf_single_trace_to_channel(Elevation, BasinNodes[i], ChannelNodeNetwork, D_inf_Flowdir, vv_temp, v_temp, output_channel_node, skip_trace);
    }
    else hillslope_node_test = false;
    if(skip_trace == false)
    {
      if(hillslope_node_test == true)
      {
        ChannelNodes.push_back(output_channel_node);
        HillslopeNodes.push_back(BasinNodes[i]);
      }
      else ++channel_node_count;
    }
    else cout << "\n\t trace failed - skipping!" << endl;
  }
  cout << "\t\t hillslope flow routing complete - note within basin there were " << channel_node_count << " channel pixels from " << N_BasinNodes << "; moving on..." << endl;
  // return output vectors
  baselevel_channel_nodes = ChannelNodes;
  hillslope_nodes = HillslopeNodes;
}

//----------------------------------------------------------------------------------------
// get_channel_characteristics_for_nodes => should be in LSDChiNetwork or LSDChannel
//----------------------------------------------------------------------------------------
// This function gets the values from a given channel characteristic corresponding to the
// base level channel node for given hillslope nodes.  For example channel longitudinal
// coordinate, chi coordinate or channel number.
// vector<float> LSDRaster::get_channel_characteristics_for_nodes(LSDFlowInfo& FlowInfo, vector<float> channel_characteristics, vector<int>& node_indices)
// {
//
// }
// vector<int> LSDRaster::get_channel_characteristics_for_nodes(LSDFlowInfo& FlowInfo, vector<int> channel_characteristics, vector<int>& node_indices)
// {
//
// }

//----------------------------------------------------------------------------------------
//
//  .----..-.    .----.  .----. .----. .----. .-.     .--.  .-..-. .-. .----.
//  | {_  | |   /  {}  \/  {}  \| {}  \| {}  }| |    / {} \ | ||  `| |{ {__
//  | |   | `--.\      /\      /|     /| .--' | `--./  /\  \| || |\  |.-._} }
//  `-'   `----' `----'  `----' `----' `-'    `----'`-'  `-'`-'`-' `-'`----'
//
//----------------------------------------------------------------------------------------

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Calculate relief relative to channel
// This calculates relief of each pixel compared to the nearest channel pixel
// Uses a threshold stream order to avoid small tributaries
//
// FJC 17/11/15
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDJunctionNetwork::calculate_relief_from_channel(LSDRaster& ElevationRaster, LSDFlowInfo& FlowInfo, int threshold_SO)
{
  Array2D<float> ReliefArray(NRows, NCols, NoDataValue);

  for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
      float this_elevation = ElevationRaster.get_data_element(row,col);
      if (this_elevation != NoDataValue)
      {
        //get the nearest channel pixel
        int CurrentNode = FlowInfo.retrieve_node_from_row_and_column(row,col);
        int BaseLevel = FlowInfo.is_node_base_level(CurrentNode);
        //if already at a channel then set relief to 0
        if (StreamOrderArray[row][col] != NoDataValue && StreamOrderArray[row][col] >= threshold_SO
        && BaseLevel == 0)
        {
          ReliefArray[row][col] = 0;
        }
        //if not at a channel, move downstream
        else
        {
          bool ReachedChannel = false;
          while (ReachedChannel == false)
          {
            //get receiver information
            int ReceiverNode, ReceiverRow, ReceiverCol;
            FlowInfo.retrieve_receiver_information(CurrentNode, ReceiverNode, ReceiverRow, ReceiverCol);
            //if node is at baselevel then exit
            if (CurrentNode == ReceiverNode)
            {
              ReachedChannel = true;
            }
            //if receiver is a channel > threshold then get the relief
            if (StreamOrderArray[ReceiverRow][ReceiverCol] != NoDataValue &&
            StreamOrderArray[ReceiverRow][ReceiverCol] >= threshold_SO)
            {
              ReachedChannel = true;
              float channel_elevation = ElevationRaster.get_data_element(ReceiverRow, ReceiverCol);
              //get the relief of the pixel (Pixel Elevation - Channel Elevation)
              ReliefArray[row][col] = (this_elevation - channel_elevation);
            }
            else
            {
              //move downstream
              CurrentNode = ReceiverNode;
            }
          }
        }
      }
    }
  }

  LSDRaster Relief(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ReliefArray, GeoReferencingStrings);
  return Relief;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Calculate relief relative to channel from connected components
// This calculates relief of each pixel compared to the nearest channel pixel
// Uses a threshold stream order for each connected components patch so that the
// whole patch is connected to the same channel.
//
// FJC 29/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDJunctionNetwork::calculate_relief_from_channel_connected_components(LSDRaster& ElevationRaster, LSDIndexRaster& ConnectedComponents, LSDRaster& DistFromOutlet, LSDFlowInfo& FlowInfo, int threshold_SO, int search_distance)
{
  Array2D<float> ReliefArray(NRows, NCols, NoDataValue);
  Array2D<int> Elevations = Get_Elevation_of_Nearest_Channel_for_Connected_Components(ConnectedComponents, ElevationRaster, DistFromOutlet, FlowInfo, threshold_SO, search_distance);
	cout << "Got the elevations of channel reaches" << endl;

  //calculate the relief (elevation of node - elevation of channel reach)
	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			int patch_id = ConnectedComponents.get_data_element(row,col);
			if (patch_id != NoDataValue)
			{
				// get the elevation of this pixel
				float this_elev = ElevationRaster.get_data_element(row,col);
				// get the elevation of the channel reach
				float channel_elev = Elevations[row][col];
				float relief = this_elev - channel_elev;
				if (relief < 0) { relief = 0; }
				ReliefArray[row][col] = relief;
			}
		}
	}

  LSDRaster Relief(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, ReliefArray, GeoReferencingStrings);
  return Relief;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes in a raster of connected component patches. It finds
// the stream order of the nearest channel for the patch.  Returns an array
// with the node index of the nearest channel for the patch.
//
// FJC 29/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Array2D<int> LSDJunctionNetwork::Get_Elevation_of_Nearest_Channel_for_Connected_Components(LSDIndexRaster& ConnectedComponents, LSDRaster& ElevationRaster, LSDRaster& DistFromOutlet, LSDFlowInfo& FlowInfo, int threshold_SO, int search_distance)
{
	Array2D<int> Elevations(NRows,NCols,NoDataValue);
	Array2D<float> FlowLengths(NRows,NCols,NoDataValue);
	Array2D<int> ChannelNodes(NRows,NCols,NoDataValue);
	Array2D<int> PatchIDs = ConnectedComponents.get_RasterData();

  // Get an array with the nearest channels for each pixel in the patch
	for (int row = 0; row < NRows; row++)
  {
    for (int col = 0; col < NCols; col++)
    {
			// find the patch ID of this node
			int patch_id = PatchIDs[row][col];
			if (patch_id != NoDataValue)
			{
				// get the flow length of the nearest channel
				int CurrentNode = FlowInfo.retrieve_node_from_row_and_column(row,col);
				int ChannelNode;
				float FlowLength, DistanceUpstream;
				get_info_nearest_channel_to_node(CurrentNode, threshold_SO, FlowInfo, DistFromOutlet, ChannelNode, FlowLength, DistanceUpstream);
				//cout << "Channel node: " << ChannelNode << " Flow length: " << FlowLength << endl;
				FlowLengths[row][col] = FlowLength;
				ChannelNodes[row][col] = ChannelNode;
			}
		}
	}

	cout << "Got the flow lengths and nodes for each patch, now finding nearest channel..." << endl;

	// Find the nearest channel node for each patch ID

	vector<int> PatchIDs_vector = Flatten_Without_Nodata(PatchIDs, NoDataValue);
	vector<float> FlowLengths_vector = Flatten_Without_Nodata(FlowLengths, NoDataValue);
	vector<int> ChannelNodes_vector = Flatten_Without_Nodata(ChannelNodes, NoDataValue);
	 //get unique patch IDs
  vector<int> Unique_Patches = Unique(PatchIDs_vector);
	vector<int> Elevation_vector;

	for (int i =0; i < int(Unique_Patches.size()); i++)
	{
		float ShortestLength = 100000000000;
		int NearestChannel = 0;
		for (int j = 0; j < int (PatchIDs_vector.size()); j++)
		{
			// find the nearest channel node
			if (PatchIDs_vector[j] == Unique_Patches[i])
			{
				if (FlowLengths_vector[j] < ShortestLength)
				{
					//update the flow length and node
					ShortestLength = FlowLengths_vector[j];
					NearestChannel = ChannelNodes_vector[j];
				}
			}
		}
		//cout << "Length: " << ShortestLength << " Channel node: " << NearestChannel << endl;

		// Get the average elevation of the reach for this patch ID
		float MeanElev = find_mean_elevation_of_channel_reach(NearestChannel,search_distance,ElevationRaster,FlowInfo);
		Elevation_vector.push_back(MeanElev);
		//cout << "Patch ID: " << Unique_Patches[i] << " Elevation of channel: " << MeanElev << endl;
	}

	// Update the array with the nearest channel node for each pixel
	cout << "Updating array with nearest channels" << endl;

	vector<int>::iterator it;
	for (int row = 0; row < NRows; row++)
	{
		for (int col = 0; col < NCols; col++)
		{
			// find the patch ID of this node
			int patch_id = PatchIDs[row][col];
			if (patch_id != NoDataValue)
			{
				// find the elevation for this patch ID
				it = find(Unique_Patches.begin(), Unique_Patches.end(), patch_id);
				int index = it - Unique_Patches.begin();
				//scout << "Index: " << index << endl;
				// update the vector with the elevation for this patch ID
				Elevations[row][col] = Elevation_vector[index];
			}
		}
	}
	cout << "Got the nearest channel elevations" << endl;

	return Elevations;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a node on the channel network and looks a specified distance
// upstream and downstream - it then calculates the average elevation of the
// reach.
//
// FJC 29/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDJunctionNetwork::find_mean_elevation_of_channel_reach(int StartingNode, int search_distance, LSDRaster& ElevationRaster, LSDFlowInfo& FlowInfo)
{
	int upstream_dist = 0;
	int downstream_dist = 0;
	int row, col, this_node;

	this_node = StartingNode;
	vector<float> elevations;
	//go upstream
	while (upstream_dist <= search_distance)
	{
		int SO_test = 0;
		//get the current stream order
		FlowInfo.retrieve_current_row_and_col(this_node, row, col);
		int this_SO = StreamOrderArray[row][col];

		//look through the donor nodes for the same stream order
		vector<int> donor_nodes = FlowInfo.get_donor_nodes(this_node);
		for (int i = 0; i < int(donor_nodes.size()); i++)
		{
			int donor_row, donor_col;
			FlowInfo.retrieve_current_row_and_col(donor_nodes[i], donor_row, donor_col);
			int DonorSO = StreamOrderArray[donor_row][donor_col];
			if (DonorSO == this_SO)
			{
				SO_test = 1;
				//push back the current elevation to the vector
				elevations.push_back(ElevationRaster.get_data_element(row,col));
				//move upstream
				this_node = donor_nodes[i];
				upstream_dist++;
			}
		}
		if (SO_test == 0)
		{
			cout << "You have reached a tributary junction, I won't check any further upstream" << endl;
			break;
		}
	}

	//go downstream
	while (downstream_dist <= search_distance)
	{
		//get current stream order
		FlowInfo.retrieve_current_row_and_col(this_node, row, col);
		int this_SO = StreamOrderArray[row][col];

		//get receiver info
		int receiver_node, receiver_row, receiver_col;
		FlowInfo.retrieve_receiver_information(this_node, receiver_node, receiver_row, receiver_col);
		int receiver_SO = StreamOrderArray[receiver_row][receiver_col];

		//push back receiver elevation to the vector
		elevations.push_back(ElevationRaster.get_data_element(receiver_row, receiver_col));

		//check if you've reached a junction
		if (this_SO == receiver_SO)
		{
			// move downstream
			this_node = receiver_node;
			downstream_dist++;
		}
		else
		{
			cout << "You've reached a junction, stopping" << endl;
			break;
		}
	}

	//calculate the mean elevation for this reach
	float mean_elev = get_mean(elevations);
	return mean_elev;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function takes a point on the channel network and finds the distance to the
// nearest predicted floodplain initiation point from the FIRTH method.  It requires the
// node index of the point on the channel network, and the connected components
// raster from the FIRTH method. It also requires the search distance (number of pixels // to search upstream and downstream).
//
// FJC 08/09/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDJunctionNetwork::find_distance_to_nearest_floodplain_pixel(int point_node, int search_distance, LSDRaster& FloodplainRaster, LSDFlowInfo& FlowInfo)
{
	// find out if you are already in the floodplain
	int row, col, this_node;
	float root_2 = 1.4142135623;
	bool reached_FIP = false;
	float distance;
	int upstream_dist = 0;
	int downstream_dist = 0;
	this_node = point_node;

	//check the upstream direction
	while (reached_FIP == false && upstream_dist <= search_distance)
	{
		int SO_test = 0;
		FlowInfo.retrieve_current_row_and_col(this_node, row, col);
		int this_SO = StreamOrderArray[row][col];
		int this_FP = FloodplainRaster.get_data_element(row, col);
		vector<int> donor_nodes = FlowInfo.get_donor_nodes(this_node);
		for (int i = 0; i < int(donor_nodes.size()); i++)
		{
			int donor_row, donor_col;
			FlowInfo.retrieve_current_row_and_col(donor_nodes[i], donor_row, donor_col);
			int DonorSO = StreamOrderArray[donor_row][donor_col];
			if (DonorSO == this_SO)
			{
				SO_test = 1;
				//at the FIP you are in the FIP but the donor is not
				int DonorFP = FloodplainRaster.get_data_element(donor_row, donor_col);
				if (DonorFP == NoDataValue && this_FP != NoDataValue)
				{
					// you've reached the FIP!
					cout << "Reached FIP" << endl;
					reached_FIP = true;
				}
				else
				{
					this_node = donor_nodes[i];
					//update length
					if (FlowInfo.retrieve_flow_length_code_of_node(donor_nodes[i]) == 1){ upstream_dist += DataResolution; }
      		else if (FlowInfo.retrieve_flow_length_code_of_node(donor_nodes[i]) == 2){ upstream_dist += (DataResolution * root_2); }
					else if (FlowInfo.retrieve_flow_length_code_of_node(donor_nodes[i]) == 0){ break; }
				}
			}
		}
		if (SO_test == 0)
		{
			cout << "You have reached a tributary junction, I won't check any further upstream" << endl;
			break;
		}
	}
	cout << "Now checking the downstream direction" << endl;

	//check the downstream direction
	reached_FIP = false;
	this_node = point_node;
	while (reached_FIP == false && downstream_dist <= search_distance)
	{
		FlowInfo.retrieve_current_row_and_col(this_node, row, col);
		int this_FP = FloodplainRaster.get_data_element(row,col);
		if (this_FP == NoDataValue)
		{
			int receiver_node, receiver_row, receiver_col;
			FlowInfo.retrieve_receiver_information(this_node, receiver_node, receiver_row, receiver_col);
			int receiver_FP = FloodplainRaster.get_data_element(receiver_row, receiver_col);
			if (receiver_FP != NoDataValue)
			{
				reached_FIP = true;
			}
			else
			{
				this_node = receiver_node;
				//update length
				if (FlowInfo.retrieve_flow_length_code_of_node(receiver_node) == 1){ downstream_dist += DataResolution; }
      	else if (FlowInfo.retrieve_flow_length_code_of_node(receiver_node) == 2){ downstream_dist += (DataResolution * root_2); }
			}
		}
		else
		{
			downstream_dist = 1000000;
		}
	}

	//find the nearest node
	if (upstream_dist < downstream_dist)
	{
		distance = upstream_dist*-1;
	}
	else if (downstream_dist < upstream_dist)
	{
		distance = downstream_dist;
	}
	else if (upstream_dist == downstream_dist && upstream_dist < search_distance)
	{
		distance = upstream_dist*-1;
	}
  else
	{
		cout << "I couldn't find a FIP within the search radius, returning NDV" << endl;
		distance = NoDataValue;
	}

	return distance;
	cout << "Got the distance for this FIP" << endl;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// END OF FLOODPLAIN FUNCTIONS

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function prints all junctions to a csv file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::print_junctions_to_csv(LSDFlowInfo& FlowInfo, string fname)
{

  int this_node;
  int row,col;
  double x_loc,y_loc;
  double latitude,longitude;

  vector<int> JunctionList;
  for (int i = 0; i<NJunctions; i++)
  {
    JunctionList.push_back(i);
    //cout << "The stream order of junction " << i << " is " << StreamOrderVector[i] << endl;
  }

  // open the outfile
  ofstream sources_out;
  sources_out.open(fname.c_str());
  sources_out.precision(9);

  sources_out << "junction,node,x,y,latitude,longitude, stream_order" << endl;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  for (int i = 0; i<NJunctions; i++)
  {
    this_node = get_Node_of_Junction(JunctionList[i]);

    // get the row and column
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);

    // get the x and y locations
    FlowInfo.get_x_and_y_locations(row, col, x_loc, y_loc);

    // get the lat and long locations
    FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    // print to file
    sources_out << JunctionList[i] << "," << this_node << "," << x_loc << ","
                << y_loc << "," << latitude << "," << longitude << "," << StreamOrderVector[ JunctionList[i] ] << endl;

  }
  sources_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function takes a list of junctions and prints them to a csv file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::print_junctions_to_csv(LSDFlowInfo& FlowInfo, vector<int> JunctionList, string fname)
{
  int n_junctions = int(JunctionList.size());
  int this_node;
  int row,col;
  double x_loc,y_loc;
  double latitude,longitude;

  int NJunctions = int(JunctionVector.size());
  if (n_junctions == 0)
  {
    cout << "You passed me an empty junction list. Printing all the junctions!" << endl;
    for (int i = 0; i<NJunctions; i++)
    {
      JunctionList.push_back(i);
    }
    n_junctions  = NJunctions;
  }

  // open the outfile
  ofstream sources_out;
  sources_out.open(fname.c_str());
  sources_out.precision(9);

  sources_out << "junction,node,x,y,latitude,longitude,stream_order" << endl;

  // this is for latitude and longitude
  LSDCoordinateConverterLLandUTM Converter;

  for (int i = 0; i<n_junctions; i++)
  {
    this_node = get_Node_of_Junction(JunctionList[i]);

    // get the row and column
    FlowInfo.retrieve_current_row_and_col(this_node,row,col);

    // get the x and y locations
    FlowInfo.get_x_and_y_locations(row, col, x_loc, y_loc);

    // get the lat and long locations
    FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);

    // print to file
    sources_out << JunctionList[i] << "," << this_node << "," << x_loc << ","
                << y_loc << "," << latitude << "," << longitude << "," << StreamOrderVector[ JunctionList[i] ] << endl;

  }

  sources_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating a bonehead version of the chi slope
// and the chi intercept
//
// What this does is it takes a list of baselevel junctions
// It then finds all sources in the basin
// It then sorts these according to length
// Starting with the longest channel it moves down marking channels as read
// For all but the longest channel, it will eventually come upon a channel that has
// already been visited. So it will stop there.
// BUT we want to keep a few chi nodes downslope (to make sure there are not
// big jumps in chi at tributary junctions)
// so it follows the visited channel down a few nodes. It does this using the
// flowinfo function   get_downslope_node_after_fixed_visited_nodes
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_overlapping_channels(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    int n_nodes_to_visit)
{
  // Get the number of baselevel nodes
  int N_baselevel_nodes = int(BaseLevel_Junctions.size());

  // create the visited array
  int not_visited = 0;
  LSDIndexRaster VisitedRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, GeoReferencingStrings,not_visited);

  vector<int> NewSources;
  vector<int> NewOutlets;
  int thisOutlet;

  // loop through these nodes
  for (int BL = 0; BL < N_baselevel_nodes; BL++)
  {
    int outlet_node = JunctionVector[BaseLevel_Junctions[BL] ];
    //cout << "The outlet node is: " << outlet_node << endl;

    // get all the source nodes of the base level
    vector<int> source_nodes = get_all_source_nodes_of_an_outlet_junction(BaseLevel_Junctions[BL]);
    //cout << "The number of sources is: " << source_nodes.size() << endl;


    // sort the nodes by flow distance in ascending order
    vector<int> SortedSources = FlowInfo.sort_node_list_based_on_raster(source_nodes, DistanceFromOutlet);

    // get them in descending order
    reverse(SortedSources.begin(),SortedSources.end());

    // now loop through the sorted sources
    int n_sources = int(SortedSources.size());

    for(int s = 0; s<n_sources; s++)
    {
      // get the channel from this source and mark up the covered raster
      thisOutlet = FlowInfo.get_downslope_node_after_fixed_visited_nodes(SortedSources[s],
                  outlet_node, n_nodes_to_visit, VisitedRaster);

      //cout << "Source number " << s << " source node is: " << SortedSources[s] << " BL node: " << outlet_node
      //     << " and new outlet: " << thisOutlet << endl;


      NewSources.push_back(SortedSources[s]);
      NewOutlets.push_back(thisOutlet);
    }

  }

  outlet_nodes = NewOutlets;
  source_nodes = NewSources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating a bonehead version of the chi slope
// and the chi intercept
// overlaoded, overwrite the baselevel nodes vector (used for visualisation).
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_overlapping_channels(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    vector<int>& baselevel_nodes,
                                    int n_nodes_to_visit)
{
  // Get the number of baselevel nodes
  int N_baselevel_nodes = int(BaseLevel_Junctions.size());

  // create the visited array
  int not_visited = 0;
  LSDIndexRaster VisitedRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, GeoReferencingStrings,not_visited);

  vector<int> NewSources;
  vector<int> NewOutlets;
  vector<int> NewBaselevelNodes;
  int thisOutlet;

  // loop through these nodes
  for (int BL = 0; BL < N_baselevel_nodes; BL++)
  {
    int outlet_node = JunctionVector[BaseLevel_Junctions[BL] ];
    //cout << "The outlet node is: " << outlet_node << endl;

    // get all the source nodes of the base level
    vector<int> source_nodes = get_all_source_nodes_of_an_outlet_junction(BaseLevel_Junctions[BL]);
    //cout << "The number of sources is: " << source_nodes.size() << endl;


    // sort the nodes by flow distance in ascending order
    vector<int> SortedSources = FlowInfo.sort_node_list_based_on_raster(source_nodes, DistanceFromOutlet);

    // get them in descending order
    reverse(SortedSources.begin(),SortedSources.end());

    // now loop through the sorted sources
    int n_sources = int(SortedSources.size());

    for(int s = 0; s<n_sources; s++)
    {
      // get the channel from this source and mark up the covered raster
      thisOutlet = FlowInfo.get_downslope_node_after_fixed_visited_nodes(SortedSources[s],
                  outlet_node, n_nodes_to_visit, VisitedRaster);

      //cout << "Source number " << s << " source node is: " << SortedSources[s] << " BL node: " << outlet_node
      //     << " and new outlet: " << thisOutlet << endl;


      NewSources.push_back(SortedSources[s]);
      NewOutlets.push_back(thisOutlet);
      NewBaselevelNodes.push_back(outlet_node);
    }

  }

  outlet_nodes = NewOutlets;
  source_nodes = NewSources;
  baselevel_nodes = NewBaselevelNodes;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating a bonehead version of the chi slope
// and the chi intercept and is written for the condition that baselevel junctions
// have been specified by the user
// MDH 19/6/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_overlapping_channels_to_downstream_outlets(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    int n_nodes_to_visit)
{
  // Get the number of baselevel nodes
  int N_baselevel_nodes = int(BaseLevel_Junctions.size());

  // create the visited array
  int not_visited = 0;
  LSDIndexRaster VisitedRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, GeoReferencingStrings,not_visited);

  vector<int> NewSources;
  vector<int> NewOutlets;
  int thisOutlet;

  // loop through these nodes
  for (int BL = 0; BL < N_baselevel_nodes; BL++)
  {
    int outlet_node = get_penultimate_node_from_stream_link(BaseLevel_Junctions[BL],FlowInfo);
    //cout << "The outlet node is: " << outlet_node << endl;

    // get all the source nodes of the base level
    vector<int> source_nodes = get_all_source_nodes_of_an_outlet_junction(BaseLevel_Junctions[BL]);
    //cout << "The number of sources is: " << source_nodes.size() << endl;


    // sort the nodes by flow distance in ascending order
    vector<int> SortedSources = FlowInfo.sort_node_list_based_on_raster(source_nodes, DistanceFromOutlet);

    // get them in descending order
    reverse(SortedSources.begin(),SortedSources.end());

    // now loop through the sorted sources
    int n_sources = int(SortedSources.size());

    for(int s = 0; s<n_sources; s++)
    {
      // get the channel from this source and mark up the covered raster
      thisOutlet = FlowInfo.get_downslope_node_after_fixed_visited_nodes(SortedSources[s],
                  outlet_node, n_nodes_to_visit, VisitedRaster);

      //cout << "Source number " << s << " source node is: " << SortedSources[s] << " BL node: " << outlet_node
      //     << " and new outlet: " << thisOutlet << endl;


      NewSources.push_back(SortedSources[s]);
      NewOutlets.push_back(thisOutlet);
    }

  }

  outlet_nodes = NewOutlets;
  source_nodes = NewSources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating a bonehead version of the chi slope
// and the chi intercept and is written for the condition that baselevel junctions
// have been specified by the user
// Overloaded version that spits out the baselevel nodes as well.
// This is a bit screwy as the baslevel node will be the next junction down
//
// This is working now but need to check if going upstream on junction leads
// to funny basin shapes
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_overlapping_channels_to_downstream_outlets(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    vector<int>& baselevel_nodes,
                                    int n_nodes_to_visit)
{
  // Get the number of baselevel nodes
  int N_baselevel_nodes = int(BaseLevel_Junctions.size());

  // create the visited array
  int not_visited = 0;
  LSDIndexRaster VisitedRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, GeoReferencingStrings,not_visited);

  vector<int> NewSources;
  vector<int> NewOutlets;
  vector<int> NewBaselevelNodes;
  int thisOutlet;
  //int outlet_junction;

  // loop through these nodes
  for (int BL = 0; BL < N_baselevel_nodes; BL++)
  {
    int outlet_node = get_penultimate_node_from_stream_link(BaseLevel_Junctions[BL],FlowInfo);
    int outlet_junction_node = JunctionVector[BaseLevel_Junctions[BL] ];
    //cout << "The outlet node is: " << outlet_node << endl;

    // get all the source nodes of the base level
    vector<int> source_nodes = get_all_source_nodes_of_an_outlet_junction(BaseLevel_Junctions[BL]);
    //cout << "The number of sources is: " << source_nodes.size() << endl;


    // sort the nodes by flow distance in ascending order
    vector<int> SortedSources = FlowInfo.sort_node_list_based_on_raster(source_nodes, DistanceFromOutlet);

    // get them in descending order
    reverse(SortedSources.begin(),SortedSources.end());

    // now loop through the sorted sources
    int n_sources = int(SortedSources.size());

    for(int s = 0; s<n_sources; s++)
    {
      // get the channel from this source and mark up the covered raster
      thisOutlet = FlowInfo.get_downslope_node_after_fixed_visited_nodes(SortedSources[s],
                  outlet_node, n_nodes_to_visit, VisitedRaster);

      //cout << "Source number " << s << " source node is: " << SortedSources[s] << " BL node: " << outlet_node
      //     << " and new outlet: " << thisOutlet << endl;


      NewSources.push_back(SortedSources[s]);
      NewOutlets.push_back(thisOutlet);


      //cout << "Search for Heybubba: This will probably cause a segfault. It needs to be a junction!" << endl;
      NewBaselevelNodes.push_back(outlet_junction_node);

    }

  }

  outlet_nodes = NewOutlets;
  source_nodes = NewSources;
  baselevel_nodes = NewBaselevelNodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating a bonehead version of the chi slope
// and the chi intercept and is written for the condition that baselevel junctions
// have been specified by the user
// Overloaded version that spits out the baselevel nodes as well.
// This is a bit screwy as the baslevel node will be the next junction down
//
// This is working now but need to check if going upstream on junction leads
// to funny basin shapes
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::get_overlapping_channels_to_downstream_outlets(LSDFlowInfo& FlowInfo,
                                    vector<int> BaseLevel_Junctions,
                                    LSDRaster& DistanceFromOutlet,
                                    vector<int>& source_nodes,
                                    vector<int>& outlet_nodes,
                                    vector<int>& baselevel_nodes,
                                    int n_nodes_to_visit,
                                    vector<int>& input_nodes)
{
  // Get the number of baselevel nodes
  int N_baselevel_nodes = int(BaseLevel_Junctions.size());

  // create the visited array
  int not_visited = 0;
  LSDIndexRaster VisitedRaster(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, GeoReferencingStrings,not_visited);

  vector<int> NewSources;
  vector<int> NewOutlets;
  vector<int> NewBaselevelNodes;
  int thisOutlet;
  //int outlet_junction;

  // loop through these nodes
  for (int BL = 0; BL < N_baselevel_nodes; BL++)
  {
    // int outlet_node = get_penultimate_node_from_stream_link(BaseLevel_Junctions[BL],FlowInfo);
    int outlet_node = input_nodes[BL];
    int outlet_junction_node = JunctionVector[BaseLevel_Junctions[BL] ];
    //cout << "The outlet node is: " << outlet_node << endl;

    // get all the source nodes of the base level
    vector<int> source_nodes = get_all_source_nodes_of_an_outlet_junction(BaseLevel_Junctions[BL]);
    //cout << "The number of sources is: " << source_nodes.size() << endl;


    // sort the nodes by flow distance in ascending order
    vector<int> SortedSources = FlowInfo.sort_node_list_based_on_raster(source_nodes, DistanceFromOutlet);

    // get them in descending order
    reverse(SortedSources.begin(),SortedSources.end());

    // now loop through the sorted sources
    int n_sources = int(SortedSources.size());
    // int this_outlet;


    for(int s = 0; s<n_sources; s++)
    {
      // get the channel from this source and mark up the covered raster

      int this_node = SortedSources[s]; int row,col; FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      VisitedRaster.set_data_element(row,col,1);
      int recnode = -9999;
      bool keepgoing = true;
      do
      {
        FlowInfo.retrieve_receiver_information(this_node,recnode);
        int rrow,rcol; FlowInfo.retrieve_current_row_and_col(recnode,rrow,rcol);
        if(this_node == recnode || VisitedRaster.get_data_element(rrow,rcol) == 1 || this_node == input_nodes[BL] )
        {
          thisOutlet = this_node;
          recnode = this_node;
          keepgoing = false;
        }
        else
        {
          VisitedRaster.set_data_element(rrow,rcol,1);
          this_node = recnode;
        }
      }while(keepgoing);


      NewSources.push_back(SortedSources[s]);
      NewOutlets.push_back(thisOutlet);


      //cout << "Search for Heybubba: This will probably cause a segfault. It needs to be a junction!" << endl;
      // NewBaselevelNodes.push_back(outlet_junction_node);
      NewBaselevelNodes.push_back(input_nodes[BL]);

    }

  }

  outlet_nodes = NewOutlets;
  source_nodes = NewSources;
  baselevel_nodes = NewBaselevelNodes;

}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// THis function aims to select basins from a node rather than a junction.
// it takes a list of outlet nodes and OVERWRITES 4 vectors: baselevel_nodes, baselevel_junctions, outlet_nodes and source_nodes
// It requires the nodes to be the EXACT location of the required outlets, and may require a previous function to select them
// B.G. - 09/12/2019 <- NEARLY CHRISTMAS
// (             )
//  `--(_   _)--'
//       Y-Y
//      /@@ \
//     /     \
//     `--'.  \             ,
//         |   `.__________/)
void LSDJunctionNetwork::select_basin_from_nodes(vector<int>& input_nodes, vector<int>& sources, vector<int>& baselevel_nodes, vector<int>& baselevel_junctions,vector<int>& outlet_nodes,LSDFlowInfo& FlowInfo, 
  LSDRaster& DistanceFromOutlet, bool check_edges)
{

  // First I want to get the baselevel junctions
  vector<int>new_baselevel_junctions;
  cout << "DEBUG::STARTING BASIN SELECTION FROM NODE" << endl;

  // preprocessing the stack for checking the edges
  vector<bool> checker;
  if(check_edges)
    checker = this->check_nodata_influence(FlowInfo, DistanceFromOutlet);
  
  vector<int> newput_doe(input_nodes.size());

  for(size_t i=0; i<input_nodes.size(); i++)
  {
    // first checking if basin cutted:
    bool is_influenced_by_nodata = false;

    if(check_edges)
      is_influenced_by_nodata = checker[input_nodes[i]];

    if(is_influenced_by_nodata)
      continue;


    int this_row, this_col;
    FlowInfo.retrieve_current_row_and_col(input_nodes[i], this_row, this_col);

    cout << "DEBUG::TRYING TO FIND::" << input_nodes[i] << " row: " << this_row << " col: " << this_col << endl;

    int this_junction = this->find_upstream_junction_from_channel_nodeindex(input_nodes[i], FlowInfo);

    cout << "DEBUG::FOUND JUNCTION::" << this_junction << endl;
    // vector<int> source_nodes = get_all_source_nodes_of_an_outlet_junction(this_junction);
    // int outlet_node = get_penultimate_node_from_stream_link(this_junction,FlowInfo);
    // int outlet_junction_node = JunctionVector[BaseLevel_Junctions[BL]];


    // // int this_outlet = input_nodes[i];
    // // int this_baselevel = input_nodes[i];
    // // new_baselevel_nodes.push_back(this_baselevel);
    new_baselevel_junctions.push_back(this_junction);
    newput_doe[i] = input_nodes[i];
    // // new_outlet_nodes.push_back(this_outlet);
    // // for(size_t Ngazul =0; Ngazul<source_nodes.size();Ngazul++)
    // //   new_sources.push_back(source_nodes[Ngazul]);
  }

  newput_doe.shrink_to_fit();
  baselevel_junctions = new_baselevel_junctions;

  this->get_overlapping_channels_to_downstream_outlets(FlowInfo,
                                    baselevel_junctions,
                                    DistanceFromOutlet,
                                    sources,
                                    outlet_nodes,
                                    baselevel_nodes,
                                    -9999,
                                    input_nodes);


  // for(auto& v:baselevel_junctions)
  // {
  //   if(v == -9999)
  //   {
  //     v = 
  //   }
  // }




  // sources = new_sources;
  // baselevel_nodes = new_baselevel_nodes;
  // outlet_nodes = new_outlet_nodes;
  // SHould be done

}



vector<bool> LSDJunctionNetwork::check_nodata_influence(LSDFlowInfo& FlowInfo, LSDRaster& testrast)
{
  vector<bool> is_blurp;
  int cpt =0, cpt2 = 0;;
  // Getting the stack
  const vector<int>& stack = FlowInfo.get_SVector();
  // assuming none of the data is affected
  is_blurp = std::vector<bool>(stack.size(),false);
  // Nothing is processed
  std::vector<bool> is_processed(stack.size(),false);

  // reverse iterations
  for(int row =0; row< NRows; row++)
  for(int col =0; col< NCols; col++)
  {
    // If already done: skip
    int this_node = FlowInfo.retrieve_node_from_row_and_column(row,col);

    // MAh node
    if(row == 0 || row == NRows - 1 || col == 0 || col == NCols-1 || testrast.get_data_element(row,col) == NoDataValue)
    {
      // is_processed[this_node] = true;
      if(this_node != -9999)
        is_blurp[this_node] = true;
      cpt++;
      //OUH
      for(int i = -1; i <= 1;i++)
      {
        for(int j = -1; j <= 1;j++)
        {
          if(i==0 && j==0)
          {
            continue;
          }

          if(i + row < 0 || i + row >= NRows || j + col < 0 || j + col >= NCols)
          {
            continue;
          }

          int neighnode = FlowInfo.retrieve_node_from_row_and_column(i+row, j+col);
          if(neighnode == -9999)
            continue;

          // is_processed[neighnode] = true;
          is_blurp[neighnode] = true;
          cpt++;
        }
      }
    }
  }

  // std::cout << "Salub  " << cpt << " || " << stack.size() << endl ;

  for(int ri = stack.size()-1; ri>=0; ri--)
  {
    int this_node = stack[ri];
    bool is_affected = is_blurp[this_node];
    if(is_affected)
    {
      is_blurp[this_node] = true;
      int tested_node = this_node;
      int recnode;FlowInfo.retrieve_receiver_information(tested_node,recnode);
      while(recnode != tested_node)
      {
        tested_node = recnode;

        is_blurp[tested_node] = true;
        cpt++;
        // is_processed[tested_node] = true;
        FlowInfo.retrieve_receiver_information(tested_node,recnode);
      }
    }
  }
  // std::cout << "Salub2  " << cpt << endl ;

  // LSDRaster this_rast;
  // Array2D<float> thisarr(NRows,NCols,float(0));

  // for( int ri = stack.size()-1; ri>=0; ri--)
  // {
  //   int row,col; FlowInfo.retrieve_current_row_and_col(stack[ri],row,col);
  //   if(is_blurp[stack[ri]])
  //     thisarr[row][col] = 1;

  // }

  // Old check
  // LSDRaster this_rast(NRows, NCols, XMinimum, YMinimum,
  //     DataResolution,NoDataValue, thisarr, GeoReferencingStrings);
  // this_rast.write_raster("FLUB","bil");

  return is_blurp;
}  


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function aims to select basins from a list of node approximatively located around the outlet of their main basins
// It snaps the inputed node to the closest LARGEST basin.
// B.G. - 10/12/2019 <- NEARLY CHRISTMAS
// (             )
//  `--(_   _)--'
//       Y-Y
//      /@@ \
//     /     \
//     `--'.  \             ,
//         |   `.__________/)
void LSDJunctionNetwork::basin_from_node_snap_to_largest_surrounding_DA(vector<int>& input_nodes, vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  int n_pixels, bool check_edges)
{

  vector<int> corrected_nodes;
  // for each of my input nodes
  for(size_t i=0; i<input_nodes.size(); i++)
  {
    // Getting various noe informations
    int this_node = input_nodes[i];
    int this_row, this_col; FlowInfo.retrieve_current_row_and_col(this_node,this_row,this_col);
    vector<int> vabul; vabul.push_back(-1); vabul.push_back(0); vabul.push_back(1);
    double max_DA = -9999;
    int target_row, target_col, target_node;
    // Actually checking around my node to find the largest DA
    DrainageArea.snap_to_row_col_with_greatest_value_in_window(this_row, this_col, target_row, target_col, n_pixels);
    target_node = FlowInfo.retrieve_node_from_row_and_column(target_row, target_col); 

    // Retrieving the targeted nodes
    corrected_nodes.push_back(target_node);
  }

  // now I can actually get the basins
  this->select_basin_from_nodes( corrected_nodes,  sources,  baselevel_nodes,  baselevel_junctions, outlet_nodes, FlowInfo,  DistanceFromOutlet, check_edges);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function aims to select basins By minimum size
// B.G. - 10/12/2019 <- NEARLY CHRISTMAS
// (             )
//  `--(_   _)--'
//       Y-Y
//      /@@ \
//     /     \
//     `--'.  \             ,
//         |   `.__________/)
vector<int> LSDJunctionNetwork::basin_from_node_minimum_DA(vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  double min_DA, bool check_edges)
{
  vector<int> target_nodes;
  vector<bool> checker;
  if(check_edges)
    checker = this->check_nodata_influence(FlowInfo, DistanceFromOutlet);

  const vector<int>& stack = FlowInfo.get_SVector();
  for(size_t i =0; i<stack.size();i++)
  { 
    int this_node = stack[i];
    int row,col; FlowInfo.retrieve_current_row_and_col(this_node,row,col);

    if(DrainageArea.get_data_element(row,col) == NoDataValue )
      continue;

    if(check_edges)
    {
      if(checker[this_node])
        continue;
    }

    if(DrainageArea.get_data_element(row,col) >= min_DA)
    {
      target_nodes.push_back(this_node);
      int ncontributing_nodes = FlowInfo.retrieve_contributing_pixels_of_node(this_node);
      i += size_t(ncontributing_nodes);
    }
  }


  if(target_nodes.size() == 0)
  {

    cout << endl << "FATALERROR::LSDJunctionNetwork::basin_from_node_snap_to_largest_surrounding_DA" << endl;
    cout << "No basin below your minimum size!" << endl;
    exit(EXIT_FAILURE);
  }


  this->select_basin_from_nodes( target_nodes,  sources,  baselevel_nodes,  baselevel_junctions, outlet_nodes, FlowInfo,  DistanceFromOutlet, false);
  return target_nodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function aims to select basins By minimum size
// B.G. - 10/12/2019 <- NEARLY CHRISTMAS
// (             )
//  `--(_   _)--'
//       Y-Y
//      /@@ \
//     /     \
//     `--'.  \             ,
//         |   `.__________/)
vector<int> LSDJunctionNetwork::basin_from_node_range_DA(vector<int>& sources, vector<int>& baselevel_nodes,
  vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
  double min_DA, double max_DA, bool check_edges)
{
  // Getting the stack
  const vector<int>& stack = FlowInfo.get_SVector(); 
  // Initialising the target nodes, reserve memory for worst case scenario ()
  vector<int> target_nodes;target_nodes.reserve(stack.size());

  // std::cout << "endel" << std::endl;
  vector<bool> checker;
  if(check_edges)
    checker = this->check_nodata_influence(FlowInfo, DistanceFromOutlet);

  // std::cout << "flark" << std::endl;

  // going through the stack, down to top
  for( size_t i=0; i < stack.size(); i++)
  {
    // std::cout << "FLFLFLFL" << std::endl;
    // Getting the row,col for DA
    int row,col; FlowInfo.retrieve_current_row_and_col(stack[i],row,col);
    if(DrainageArea.get_data_element(row,col) >= min_DA && DrainageArea.get_data_element(row,col) <= max_DA)
    {
      // influenced by edges
      if(check_edges)
      {  
        if(checker[stack[i]])
        {
          continue;
        }
      }
      // This node is a base-level, saving it (emplace_back is better than push_back when used with reserve)
      target_nodes.emplace_back(stack[i]);
      // Jumping to the next basin in the stack to avoid multiple nested basins
      int ncontributing_nodes = FlowInfo.retrieve_contributing_pixels_of_node(stack[i]);
      i += size_t(ncontributing_nodes);
    }
    // Next test
  }
  // Fixing my vector memory
  target_nodes.shrink_to_fit();

  // Double checking few things
  if(target_nodes.size() == 0)
  {

    cout << endl << "FATALERROR::LSDJunctionNetwork::basin_from_node_snap_to_largest_surrounding_DA" << endl;
    cout << "No basin below your minimum size!" << endl;
    exit(EXIT_FAILURE);
  }

  this->select_basin_from_nodes( target_nodes,  sources,  baselevel_nodes,  baselevel_junctions, outlet_nodes, FlowInfo,  DistanceFromOutlet, false);
  
  return target_nodes;
}

std::vector<int> LSDJunctionNetwork::basin_from_node_minimum_DA_draining_to_list_of_nodes(vector<int>& input_nodes, vector<int>& sources, vector<int>& baselevel_nodes,
vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
double min_DA)
{

  // First I need to make sure I ignore my nodes!
  map<int,bool> dontdo;
  for(size_t no=0; no<input_nodes.size(); no++)
  {
    dontdo[input_nodes[no]] = true;
  }

  vector<int> target_nodes;
  // cout << "Strug" << endl;
  for(size_t no=0; no<input_nodes.size(); no++)
  {
    // cout << "Geft::" << no << "/" << input_nodes.size() << "||" << input_nodes[no] << endl;

    int this_node = input_nodes[no];
    if(this_node == NoDataValue)
      continue;
    // cout << "gust" << endl;
    vector<int> nenodes = FlowInfo.get_donor_nodes(this_node);
    // cout << "dolik" << endl;
    for(size_t on=0; on<nenodes.size(); on++)
    {
      // cout << "dolik::nono::" <<on << "/" <<  nenodes.size()<< endl;
      if(nenodes[on] == NoDataValue)
        continue;

      int this_row,this_col; FlowInfo.retrieve_current_row_and_col(nenodes[on], this_row,this_col);

      if(DrainageArea.get_data_element(this_row,this_col) == NoDataValue)
      {
        // cout << "FLEKT?!" << endl;
        continue;
      }

      if( dontdo.find(nenodes[on]) != dontdo.end() || dontdo[nenodes[on]] == true )
        continue;

      if(DrainageArea.get_data_element(this_row,this_col)>=min_DA )
      {
        target_nodes.push_back(nenodes[on]);
        // cout << dontdo[nenodes[on]] << endl;
      }

    }

  }
  // cout << "GAbul::" << target_nodes.size() << endl;
  // exit(EXIT_FAILURE);


  if(target_nodes.size() == 0)
  {

    cout << endl << "FATALERROR::LSDJunctionNetwork::basin_from_node_snap_to_largest_surrounding_DA" << endl;
    cout << "No basin below your minimum size!" << endl;
    exit(EXIT_FAILURE);
  }
  this->select_basin_from_nodes( target_nodes,  sources,  baselevel_nodes,  baselevel_junctions, outlet_nodes, FlowInfo,  DistanceFromOutlet, false);

  return target_nodes;
}




std::vector<int> LSDJunctionNetwork::basin_from_node_all_minimum_DA_for_one_watershed(int outlet_node_of_the_watershed, vector<int>& sources, vector<int>& baselevel_nodes,
vector<int>& baselevel_junctions,vector<int>& outlet_nodes, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& DrainageArea,
double min_DA, double max_DA)
{

  // I first need to get the river stack from the baselevel
  vector<int> SV = FlowInfo.get_SVector(), new_stack;
  std::cout << SV.size() << endl;
  map<int,bool> visited;

  // First, I am building a new stack and initialising a visited map.
  // gather all the nodes in my basin and now shich one I need to do or not
  int index_baselevel_node = -9999;
  bool is_visited = false;
  cout << "1" << endl;
  for(size_t i=0; i < SV.size() && is_visited == false; i++)
  {
    // If I haven't found my baselevel yet: ignore that node
    if(SV[i] != outlet_node_of_the_watershed && is_visited == false)
    {
      visited[SV[i]] = true;
    }
    // If this is my vbaselevel: start to gather the node
    else if(is_visited == false)
    {
      index_baselevel_node = i;
      is_visited = true;
      visited[SV[i]] = false;
      new_stack = FlowInfo.get_upslope_nodes_include_outlet(SV[i]);
      // std::cout << "FOUND IT" << std::endl;
    }
  }
  cout << "2" << endl;

  for(size_t i=0; i < new_stack.size(); i++)
    visited[new_stack[i]] = false;

  // this will host all the baselevel nodes for the next basins
  vector<int> target_nodes;

  // I have my stack, let's go through it
  for(size_t i=0; i < new_stack.size(); i++)
  {
    int this_node = new_stack[i], this_row, this_col;
    // If my node has already been visited -> ignore
    if(visited[this_node])
      continue;

    FlowInfo.retrieve_current_row_and_col(this_node,this_row,this_col);

    if(DrainageArea.get_data_element(this_row,this_col) <= max_DA && DrainageArea.get_data_element(this_row,this_col) >= min_DA)
    {
      // std::cout << "FOUND ONE" << std::endl;
      // this node will be a base level
      target_nodes.push_back(this_node);
      // gathering all the nodes draining
      vector<int> these_nodes = FlowInfo.get_upslope_nodes_include_outlet(this_node);
      for(size_t u=0; u<these_nodes.size(); u++)
      {
        visited[these_nodes[u]] = true;
      }
    }
    else
      visited[this_node] = true;
  }
  cout << "3" << endl;

  // I have all my baselevel nodes hopefully
  this->select_basin_from_nodes( target_nodes,  sources,  baselevel_nodes,  baselevel_junctions, outlet_nodes, FlowInfo,  DistanceFromOutlet, false);
  return target_nodes;

}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes in two vectors with the rows and cols of a line and finds
//the pixels along the line >= a threshold stream order
// FJC 16/04/17
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDJunctionNetwork::get_channel_pixels_along_line(vector<int> line_rows, vector<int> line_cols, int threshold_SO, LSDFlowInfo& FlowInfo)
{
  vector<int> outlet_nodes;
  int n_pixels = line_rows.size();
  int this_NI = 0;
  cout << "Returning the points on the line with a SO >= " << threshold_SO << endl;

  //for each pixel along the line, check if the SO is greater than the threshold
  for (int i = 0; i < n_pixels; i++)
  {
    int NI = FlowInfo.retrieve_node_from_row_and_column(line_rows[i], line_cols[i]);
    if (this_NI != NI)  // check for duplicates
    {
      int this_SO = get_StreamOrder_of_Node(FlowInfo, NI);
      if (this_SO >= threshold_SO)
      {
        //push back the node to the outlet
        outlet_nodes.push_back(NI);
      }
      this_NI = NI;
    }
  }

  return outlet_nodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a vector of basin outlet junctions and writes data about
// the longest channel in each to a csv.
// FJC 06/05/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::write_river_profiles_to_csv(vector<int>& BasinJunctions, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& Elevation, string csv_filename, int window_size)
{
  int this_node, row, col, stream_order;
  double latitude, longitude;
  float dist_from_outlet, drainage_area;
  LSDCoordinateConverterLLandUTM Converter;

  // open the csv
  ofstream chan_out;
  string this_fname = csv_filename+".csv";
  chan_out.open(this_fname.c_str());

  chan_out << "basin_id,id,node,distance_from_outlet,elevation,drainage_area,stream_order,slope,latitude,longitude" << endl;

  // for each basin, get the profile
  for (int i = 0; i < int(BasinJunctions.size()); i++)
  {
    // get the longest channel in this basin
    LSDIndexChannel ThisChannel = generate_longest_index_channel_in_basin(BasinJunctions[i],FlowInfo,DistanceFromOutlet);
    vector<int> NodeSequence = ThisChannel.get_NodeSequence();
    int UpstreamNode = NodeSequence.front();
    int DownstreamNode = NodeSequence.back();

    // make an LSDChannel object. this is the same as the index channel but has elevation data, etc.
    float downslope_chi=1;
    float m_over_n=0.5;
    LSDChannel ThisChan(UpstreamNode, DownstreamNode, downslope_chi, m_over_n, downslope_chi, FlowInfo, Elevation);
    // get channel slopes from LSDChannel
    vector<float> channel_slopes = ThisChan.calculate_channel_slopes(window_size, DistanceFromOutlet);

    for (int n = 0; n < int(NodeSequence.size()); n++)
    {
      this_node = NodeSequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);
      drainage_area = FlowInfo.get_DrainageArea_square_m(this_node);
      stream_order = get_StreamOrder_of_Node(FlowInfo, this_node);
      dist_from_outlet = FlowInfo.get_flow_length_between_nodes(this_node, DownstreamNode);

      chan_out << BasinJunctions[i] << ","
               << UpstreamNode << ","
               << this_node << ","
               << dist_from_outlet << ","
               << Elevation.get_data_element(row,col) << ","
               << drainage_area << ","
               << stream_order << ","
               << channel_slopes[n] << ",";
      chan_out.precision(9);
      chan_out << latitude << "," << longitude << endl;
    }
  }

  chan_out.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a vector of basin outlet junctions and writes data about
// all the tributaries to a csv
// FJC 06/05/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDJunctionNetwork::write_river_profiles_to_csv_all_tributaries(vector<int>& BasinJunctions, LSDFlowInfo& FlowInfo, LSDRaster& DistanceFromOutlet, LSDRaster& Elevation, string csv_filename)
{
  int this_node, row, col, stream_order;
  double latitude, longitude;
  float dist_from_outlet, drainage_area;
  LSDCoordinateConverterLLandUTM Converter;

  // open the csv
  ofstream chan_out;
  string this_fname = csv_filename+".csv";
  chan_out.open(this_fname.c_str());

  chan_out << "basin_id,id,node,distance_from_outlet,elevation,drainage_area,stream_order,latitude,longitude" << endl;

  // for each basin, get the profile
  for (int i = 0; i < int(BasinJunctions.size()); i++)
  {
    // get all the channel heads upstream of this junction
    vector<int> SourceNodes, SourceJunctions;
    for (int src = 0; src < int(SourcesVector.size()); src++)
    {
      int source_jn = get_Junction_of_Node(SourcesVector[src], FlowInfo);
      bool us = is_junction_upstream(BasinJunctions[i], source_jn);
      if (us)
      {
        SourceNodes.push_back(SourcesVector[src]);
        SourceJunctions.push_back(source_jn);
      }
    }

    // get the node of the basin junction
    int outlet_node = get_Node_of_Junction(BasinJunctions[i]);

    // now get the index channel between each source and the outlet junction
    for (int j = 0; j < int(SourceNodes.size()); j++)
    {

      LSDIndexChannel ThisChannel(SourceJunctions[j], SourceNodes[j], BasinJunctions[i], outlet_node, FlowInfo);
      vector<int> NodeSequence = ThisChannel.get_NodeSequence();
      int DownstreamNode = NodeSequence.back();
      for (int n = 0; n < int(NodeSequence.size()); n++)
      {
        this_node = NodeSequence[n];
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        dist_from_outlet = FlowInfo.get_flow_length_between_nodes(this_node, DownstreamNode);
        drainage_area = FlowInfo.get_DrainageArea_square_m(this_node);
        stream_order = get_StreamOrder_of_Node(FlowInfo, this_node);

        chan_out << BasinJunctions[i] << ","
                 << SourceJunctions[j] << ","
                 << this_node << ","
                 << dist_from_outlet << ","
                 << Elevation.get_data_element(row,col) << ","
                 << drainage_area << ","
                 << stream_order << ",";
        chan_out.precision(9);
        chan_out << latitude << "," << longitude << endl;
      }
    }
  }
  chan_out.close();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to get total length of channels upstream of a node
// FJC 30/04/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
float LSDJunctionNetwork::GetTotalChannelLengthUpstream(int this_node, LSDFlowInfo& FlowInfo)
{
  // get the nodes upslope of this node
  float TotalLength = 0;
  float root2 = 1.41421356;
  vector<int> UpslopeNodes = FlowInfo.get_upslope_nodes(this_node);

  for (int i = 0; i < int(UpslopeNodes.size()); i++)
  {
    // check if this upslope node is part of the stream network
    int this_SO = get_StreamOrder_of_Node(FlowInfo, UpslopeNodes[i]);
    if (this_SO != NoDataValue)
    {
      // if it's a channel, get the flow length code
      int FlowLengthCode = FlowInfo.retrieve_flow_length_code_of_node(UpslopeNodes[i]);
      if (FlowLengthCode == 1)
      {
        TotalLength += DataResolution; // cardinal
      }
      else if (FlowLengthCode == 2)
      {
        TotalLength += DataResolution * root2; // Diagonal
      }
    }
  }

  return TotalLength;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Function to write profiles from the channel head to a certain flow distance downstream
// This writes profiles for every source in the network!
// FJC  02/05/18
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDJunctionNetwork::write_river_profiles_to_csv_all_sources(float channel_length, int slope_window_size, LSDFlowInfo& FlowInfo, LSDRaster& Elevation, string csv_filename)
{
  int this_node, row, col;
  double latitude, longitude, x_loc, y_loc;
  float ThisLength, drainage_area;
  LSDCoordinateConverterLLandUTM Converter;

  // open the csv
  ofstream chan_out;
  string this_fname = csv_filename+".csv";
  chan_out.open(this_fname.c_str());

  chan_out << "id,node,row,column,distance_from_source,elevation,drainage_area,latitude,longitude,easting,northing" << endl;

  for (int i = 0; i < int(SourcesVector.size()); i++)
  {
      vector<int> channel_nodes;
      bool reached_end = false;
      int start_node = SourcesVector[i];
      int start_jn = get_Junction_of_Node(start_node, FlowInfo);
      this_node = start_node;
      // now go downstream until you are at the threshold channel length
      while (reached_end == false)
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        FlowInfo.get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        FlowInfo.get_x_and_y_locations(row, col, x_loc, y_loc);
        ThisLength = FlowInfo.get_flow_length_between_nodes(start_node,this_node);
        drainage_area = FlowInfo.get_DrainageArea_square_m(this_node);

        if (ThisLength <= channel_length)
        {
          // find the receiver node
          int receiver_node;
          FlowInfo.retrieve_receiver_information(this_node, receiver_node);

          // write to csv
          chan_out << start_jn << ","
                   << this_node << ","
                   << row << ","
                   << col << ","
                   << ThisLength << ","
                   << Elevation.get_data_element(row,col) << ","
                   << drainage_area << ",";
          chan_out.precision(9);
          chan_out << latitude << ","
                   << longitude << ",";
          chan_out.precision(9);
          chan_out << x_loc << "," << y_loc << endl;


          if (receiver_node == this_node)
          {
            cout << "I've reached a base level before the defined channel length, exiting" << endl;
            reached_end = true;
          }
          else
          {
            //cout << "This node: " << this_node << " receiver node: " << receiver_node << endl;
            channel_nodes.push_back(this_node);
            this_node = receiver_node;
          }
        }
        else
        {
          reached_end = true;
        }
      }
    }

    chan_out.close();
}

#endif
