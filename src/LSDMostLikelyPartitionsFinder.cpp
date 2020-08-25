//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDMostLikeleyPartitionsFinder
// Land Surface Dynamics MostLikeleyPartitionsFinder
//
// An object for extracting segments from x,y data
//  developed for the University of Edinburgh
//  Land Surface Dynamics group topographic toolbox.
//
// This object is mainly used in the analysis of channel profiles
//  transformed using the integral method of channel analysis.
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
// LSDMostLikelyPartitionsFinder.cpp
// source code for the LSDMostLikelyPartitionsFinder object
// this object looks for the most likeley partitions or segments
// of 2D data, and is principally used to identify segments of
// differing channel steepness in chi-zeta space
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
//
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 1.0    03/01/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


// TODO:
// 1. Peel out the partitioning so that you can save and load partition files.
// 3. p, t and ANCOVA test to see if regressions are significant and different

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
#include "LSDMostLikelyPartitionsFinder.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDMostLikelyPartitionsFinder_CPP
#define LSDMostLikelyPartitionsFinder_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create function, makes an partitions finder object with some x and y data.
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::create(int this_min_seg_length, vector<float> this_x_data, vector<float> this_y_data)
{
  minimum_segment_length = this_min_seg_length;
  x_data = this_x_data;
  y_data = this_y_data;

  base_sigma = 100.0;      // this is arbitrary

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function resets all the derived data members.
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::reset_derived_data_members()
{

  Array2D<float> empty_array;
  like_array =empty_array;
  m_array =empty_array;
  b_array = empty_array;
  rsquared_array = empty_array;
  DW_array = empty_array;

  vector<float> empty_dv;
  MLE_of_segments = empty_dv;

  vector< vector<int> > vv_int;
  segments_for_each_n_segments = vv_int;

  vector< vector < vector<int> > > vvvi;
  partitions = vvvi;

  vector<int> empty_vec;
  best_fit_AIC = empty_vec;
  best_fit_AICc = empty_vec;

  vector< vector<float> > empty_vecvec;
  AIC_for_each_n_segments= empty_vecvec;
  AICc_for_each_n_segments = empty_vecvec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Various algorithms are implemented for thinning the data.
// Several options are available.
//
// this one collects data as close as possible to some target dx, but does
// not modify invidivual data points.
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::thin_data_target_dx_preserve_data(float dx)
{
  vector<float> thinned_x;
  vector<float> thinned_y;

  //cout << "x_data size: " << x_data.size() << endl;

  thinned_x.push_back(x_data[0]);
  thinned_y.push_back(y_data[0]);

  float start_x =  x_data[0];
  float next_x = start_x+dx;

  int n_nodes = x_data.size();
  int last_picked = 0;

  for (int i = 1; i<n_nodes; i++)
  {
    //cout << "next x is: " << next_x << " and x is: " << x_data[i] << endl;
    if(x_data[i] >=next_x)
    {
      thinned_x.push_back(x_data[i]);
      thinned_y.push_back(y_data[i]);

      next_x += dx;
      last_picked = i;
    }

  }

  // make sure the last data element is included
  if (last_picked != n_nodes-1)
  {
    thinned_x.push_back(x_data[n_nodes-1]);
    thinned_y.push_back(y_data[n_nodes-1]);
  }

  x_data = thinned_x;
  y_data = thinned_y;
  reset_derived_data_members();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// same as above but returns an index vector into the data points that were selected
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::thin_data_target_dx_preserve_data(float dx, vector<int>& node_ref)
{
  vector<float> thinned_x;
  vector<float> thinned_y;
  vector<int> node_reference;

  //cout << "x_data size: " << x_data.size() << endl;

  thinned_x.push_back(x_data[0]);
  thinned_y.push_back(y_data[0]);
  node_reference.push_back(0);

  float start_x =  x_data[0];
  float next_x = start_x+dx;

  int n_nodes = x_data.size();
  int last_picked = 0;

  for (int i = 1; i<n_nodes; i++)
  {
    //cout << "next x is: " << next_x << " and x is: " << x_data[i] << endl;
    if(x_data[i] >=next_x)
    {
      thinned_x.push_back(x_data[i]);
      thinned_y.push_back(y_data[i]);
      node_reference.push_back(i);

      next_x += dx;
      last_picked = i;
    }
  }

  // make sure the last data element is included
  if (last_picked != n_nodes-1)
  {
    thinned_x.push_back(x_data[n_nodes-1]);
    thinned_y.push_back(y_data[n_nodes-1]);
    node_reference.push_back(n_nodes-1);
  }

  x_data = thinned_x;
  y_data = thinned_y;
  reset_derived_data_members();
  node_ref = node_reference;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// same as above but spawns a new object
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDMostLikelyPartitionsFinder LSDMostLikelyPartitionsFinder::spawn_thinned_data_target_dx_preserve_data(float dx)
{
  vector<float> thinned_x;
  vector<float> thinned_y;

  thinned_x.push_back(x_data[0]);
  thinned_y.push_back(y_data[0]);

  float start_x =  x_data[0];
  float next_x = start_x+dx;

  int n_nodes = x_data.size();
  int last_picked = 0;

  for (int i = 1; i<n_nodes; i++)
  {
    if(x_data[i] >=next_x)
    {
      thinned_x.push_back(x_data[i]);
      thinned_y.push_back(y_data[i]);

      next_x += dx;
      last_picked = i;
    }
  }

  // make sure the last data element is included
  if (last_picked != n_nodes-1)
  {
    thinned_x.push_back(x_data[n_nodes-1]);
    thinned_y.push_back(y_data[n_nodes-1]);
  }

  LSDMostLikelyPartitionsFinder new_object(minimum_segment_length,thinned_x,thinned_y);
  return new_object;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this thins the data by skipping elements.
// a positive number N means it skips N elements after each element
// a negative number -N means that after N elements it skips one element
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::thin_data_skip(int N, vector<int>& node_ref)
{
  vector<float> thinned_x;
  vector<float> thinned_y;
  vector<int> node_reference;
  int n_nodes = x_data.size();

  int node = 0;
  int n_counted = 0;

  if (N == 0)
  {
    thinned_x = x_data;
    thinned_y = y_data;

    for (int i = 0; i<n_nodes; i++)
    {
      node_reference.push_back(i);
    }
  }
  else
  {
    while (node < n_nodes)
    {
      if (N < 0)
      {
        thinned_x.push_back(x_data[node]);
        thinned_y.push_back(y_data[node]);
        node_reference.push_back(node);

        n_counted++;
        if (n_counted == -N)    // this triggers the skip
        {
          node+=2;
          n_counted = 0;
        }
        else
        {
          node++;
        }
      }
      else    // this is for skipping after each node
      {
        thinned_x.push_back(x_data[node]);
        thinned_y.push_back(y_data[node]);
        node_reference.push_back(node);
        node+=N+1;
      }
    }
  }

  // ensure the final daat element is selected
  if (node != n_nodes-1)
  {
    node = n_nodes-1;
    thinned_x.push_back(x_data[node]);
    thinned_y.push_back(y_data[node]);
    node_reference.push_back(node);
  }

  x_data = thinned_x;
  y_data = thinned_y;
  reset_derived_data_members();
  node_ref = node_reference;

  //cout << "Did the skipping (LSDMLPF line 280), n_nodes: " << n_nodes
  //     << " after skip: " << x_data.size() << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this one recasts data at fixed values of dx, using linear interpoloation
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::thin_data_target_dx_linear_interpolation(float dx)
{
  vector<float> thinned_x;
  vector<float> thinned_y;

  thinned_x.push_back(x_data[0]);
  thinned_y.push_back(y_data[0]);

  float start_x =  x_data[0];
  float next_x = start_x+dx;
  float data_slope, next_y;

  int n_nodes = x_data.size();

  for (int i = 1; i<n_nodes; i++)
  {
    if(x_data[i] >=next_x)
    {
      thinned_x.push_back(next_x);


      data_slope = (y_data[i]-y_data[i-1])/(x_data[i]-x_data[i-1]);
      next_y = data_slope*(next_x-x_data[i-1])+y_data[i-1];

      thinned_y.push_back(next_y );

      next_x += dx;
    }
  }
  x_data = thinned_x;
  y_data = thinned_y;
  reset_derived_data_members();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// same as above but spawns a new object.
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDMostLikelyPartitionsFinder LSDMostLikelyPartitionsFinder::spawn_thinned_data_target_dx_linear_interpolation(float dx)
{
  vector<float> thinned_x;
  vector<float> thinned_y;

  thinned_x.push_back(x_data[0]);
  thinned_y.push_back(y_data[0]);

  float start_x =  x_data[0];
  float next_x = start_x+dx;
  float data_slope, next_y;

  int n_nodes = x_data.size();

  for (int i = 1; i<n_nodes; i++)
  {
    if(x_data[i] >=next_x)
    {
      thinned_x.push_back(next_x);


      data_slope = (y_data[i]-y_data[i-1])/(x_data[i]-x_data[i-1]);
      next_y = data_slope*(next_x-x_data[i-1])+y_data[i-1];

      thinned_y.push_back(next_y );

      next_x += dx;
    }
  }
  LSDMostLikelyPartitionsFinder new_object(minimum_segment_length,thinned_x, thinned_y);
  return new_object;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this skips nodes but using a Monte Carlo scheme that samples random points along the channel profile
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::thin_data_monte_carlo_skip(int Mean_skip,int skip_range, vector<int>& node_ref)
{



  int minimum_skip = Mean_skip - 0.5*skip_range;
  long seed = time(NULL);

  int N = int((float(skip_range))*(ran3(&seed))+0.5)+minimum_skip;
  vector<float> thinned_x;
  vector<float> thinned_y;
  vector<int> node_reference;
  int n_nodes = x_data.size();

  int node = 0;      // the node along the channel that will be used for the next data element
  int n_counted = 0;
  int new_N_switch;    // this is used if N is negative and multiple nodes are selelcted in a row

  // if N is greater than zero, new_N_switch is true (== 1)
  if(N >= 0)
  {
    new_N_switch = 1;
  }
  else
  {
    new_N_switch = 0;
  }

  int last_node = 0;      // the last node used

  while (node < n_nodes)
  {
    last_node = node;
    if (N == 0)
    {
      thinned_x.push_back(x_data[node]);
      thinned_y.push_back(y_data[node]);
      node_reference.push_back(node);

      node++;
    }
    else if (N < 0)
    {
      thinned_x.push_back(x_data[node]);
      thinned_y.push_back(y_data[node]);
      node_reference.push_back(node);

      n_counted++;
      if (n_counted == -N)    // this triggers the skip
      {
        node+=2;
        n_counted = 0;
        new_N_switch = 1;    // once the skip is triggered, you will need a new N
      }
      else
      {
        node++;
      }
    }
    else    // this is for skipping after each node
    {
      thinned_x.push_back(x_data[node]);
      thinned_y.push_back(y_data[node]);
      node_reference.push_back(node);
      node+=N+1;
    }

    if (new_N_switch == 1)
    {
      float random_N = ran3(&seed);
      float skippy = (float(skip_range));
      N = int(skippy*(random_N)+0.5)+minimum_skip;
      //cout << "N is: " << N << " and random: " << random_N << " and skppy: " << skippy
      //     << " and skip_range: " << skip_range << " and minimum_skip: " << minimum_skip << endl;
      if(N >= 0)
      {
        new_N_switch = 1;
      }
      else
      {
        new_N_switch = 0;
      }
    }
  }

  //cout << "node: " << node << " and n_nodes: " << n_nodes << endl;

  // ensure the final data element is selected
  if (last_node != n_nodes-1)
  {
    //cout << "yo, node before = " << node;
    node = n_nodes-1;
    //cout << " and after: " << node << endl;
    thinned_x.push_back(x_data[node]);
    thinned_y.push_back(y_data[node]);
    node_reference.push_back(node);
  }

  x_data = thinned_x;
  y_data = thinned_y;
  reset_derived_data_members();
  node_ref = node_reference;

  //cout << "Did the skipping (LSDMLPF line 463), n_nodes: " << n_nodes
  //     << " after skip: " << x_data.size() << endl;
  //for (int i = 0; i< int(node_ref.size()); i++)
  //{
  //  cout << "node ref["<<i<< "]: " << node_ref[i] << endl;
  //}


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Thins object based on a monte carlo approach using a mean, max and minimum dchi
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::thin_data_monte_carlo_dchi(float mean_dchi, float variation_dchi, vector<int>& node_ref)
{

  //cout << "LSDMostLikelyPartitionsFinder, LINE 391, mean dchi: " << mean_dchi << endl;

  vector<float> thinned_x;
  vector<float> thinned_y;
  vector<int> node_reference;

  // set the max standard deviation to mean dchi so that you don't get
  // negative dchi values
  if (variation_dchi > mean_dchi)
  {
    variation_dchi = mean_dchi;
  }

  float min_dchi = mean_dchi-variation_dchi;
  float range_chi = 2*variation_dchi;

  long seed = time(NULL);

  // get dx using a random seed
  float dx = ran3(&seed)*range_chi+min_dchi;

  thinned_x.push_back(x_data[0]);
  thinned_y.push_back(y_data[0]);
  node_reference.push_back(0);

  float start_x =  x_data[0];
  float next_x = start_x+dx;

  int n_nodes = x_data.size();
  int last_picked = 0;

  for (int i = 1; i<n_nodes; i++)
  {
    //cout << "next x is: " << next_x << " and x is: " << x_data[i] << endl;
    if(x_data[i] >=next_x)
    {
      thinned_x.push_back(x_data[i]);
      thinned_y.push_back(y_data[i]);
      node_reference.push_back(i);

      dx = ran3(&seed)*range_chi+min_dchi;
      next_x += dx;
      last_picked = i;
    }
  }

  // make sure the last data element is included
  if (last_picked != n_nodes-1)
  {
    thinned_x.push_back(x_data[n_nodes-1]);
    thinned_y.push_back(y_data[n_nodes-1]);
    node_reference.push_back(n_nodes-1);
  }

  x_data = thinned_x;
  y_data = thinned_y;
  reset_derived_data_members();
  node_ref = node_reference;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Prints the x and y data to screen (for bug checking)
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::print_x_y_data_to_screen()
{
  int n_nodes = x_data.size();
  for (int i = 0; i<n_nodes; i++)
  {
    cout << "i: " << i << " \t x: " << x_data[i] << " \t y:  " << y_data[i] << endl;
  }
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This is a wrapper function to find the likelihood of all the segments
// you can enter a vector of sigmas such that each node has a different sigma value
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::best_fit_driver_AIC_for_linear_segments(vector<float> sigma_values)
{
  // cout << "best_fit_driver_AIC_for_linear_segments, getting like data" <<endl;
  calculate_segment_matrices(base_sigma);
  // cout << "best_fit_driver_AIC_for_linear_segments, got like data" <<endl;

  // get the maximum liklihood of segments
  find_max_like_of_segments();

  // cout << "||" << endl;

  get_n_segments_for_various_sigma(sigma_values);
  // cout << "!!" << endl;

  //print_AIC_and_AICc_to_screen(sigma_values);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This is a wrapper function to find the likelihood of all the segments
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::best_fit_driver_AIC_for_linear_segments(float sigma)
{
  //cout << "best_fit_driver_AIC_for_linear_segments, getting like data" <<endl;
  calculate_segment_matrices(base_sigma);
  //cout << "best_fit_driver_AIC_for_linear_segments, got like data" <<endl;

  // get the maximum liklihood of segments
  find_max_like_of_segments();

  vector<float> sigma_values;
  sigma_values.push_back(sigma);

  get_n_segments_for_various_sigma(sigma_values);

  //print_AIC_and_AICc_to_screen(sigma_values);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function returns data for a given sigma value
// the 'node' int is the index into the sigma vector
// this function replaces a large number of data elements.
// the vectors b, m, r^2 and DW are from the regressions of each segment
// the fit y is a vector with data in the x positions that is derived from
// the best fit lines
// the seg lengths are the individual segment lengths,
// this_MLE, this_n_segments and this_n_nodes are all returned
// so the user can combine two or more segments and get an AIC or AICc
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::get_data_from_best_fit_lines(int node, vector<float> sigma_values,
                      vector<float>& b_values, vector<float>& m_values,
                      vector<float>& r2_values, vector<float>& DW_values,
                       vector<float>& fitted_y,vector<int>& seg_lengths,
                 float& this_MLE, int& this_n_segments, int& this_n_nodes,
                 float& this_AIC, float& this_AICc)
{
  //cout <<"LSDMLPF line 613, getting data" << endl;

  int n_sigma_values = best_fit_AICc.size();
  if (node >= n_sigma_values)
  {
    cout << "LSDMostLikelyPartitionsFinder::get_data_from_best_fit_lines " << endl
         << "you have not calculated AIC with this many nodes" << endl;
  }


  // get the m, b, etc from the data
  int n_sigma_for_printing = node;
  int bestfit_segments_node = best_fit_AICc[n_sigma_for_printing];

  //cout << "LINE 627 getting properties, bf segments node: " << bestfit_segments_node <<  endl;
  get_properties_of_best_fit_segments(bestfit_segments_node,
                     m_values,b_values, r2_values,DW_values);

  //cout << "LINE 631 got properties" << endl;


  vector<int> segment_length;
  segment_length = segments_for_each_n_segments[ best_fit_AICc[n_sigma_for_printing] ];


  vector<float> new_sig_MLE = transform_like_from_sigma1_to_sigma2(base_sigma,
                  MLE_of_segments, sigma_values[n_sigma_for_printing]);

  vector<float> AICc_values = AICc_for_each_n_segments[node];
  float AICc_value = AICc_values[ best_fit_AICc[n_sigma_for_printing] ];

  // now print this data
  //cout << endl << endl << endl << "The data from the best fit: " << endl;
  //cout << "sigma is: " << sigma_values[n_sigma_for_printing]
  //   << " MLE: " << new_sig_MLE [ best_fit_AICc[n_sigma_for_printing] ]
  //   << " and the number of segments is: " << best_fit_AICc[n_sigma_for_printing]+1
  //   << " AICc: " << AICc_value << endl;
  //for (int i = 0; i< best_fit_AICc[n_sigma_for_printing]+1; i++)
  //{
  //  cout << "seg_length: " << segment_length[i] << " " << m_values[i] << " " << b_values[i] << " "
  //       << r2_values[i] << " " << DW_values[i] << endl;
  //}


  // create a vector of y values from the best fit segments
  int n_nodes = x_data.size();
  vector<float> fit_y(n_nodes);
  int this_node = 0;
  for (int seg = 0; seg< best_fit_AICc[n_sigma_for_printing]+1; seg++)
  {
    for(int seg_node =0; seg_node< segment_length[seg]; seg_node++)
    {
      fit_y[this_node] = m_values[seg]*x_data[this_node]+ b_values[seg];
      this_node++;
    }
  }

  // now replace data for the purposes of getting it to the calling function
  fitted_y = fit_y;
  seg_lengths = segment_length;
  this_MLE = new_sig_MLE [ best_fit_AICc[n_sigma_for_printing] ];
  this_n_segments =  best_fit_AICc[n_sigma_for_printing]+1;
  this_n_nodes = n_nodes;
  this_AICc = AICc_value;
  vector<float> AIC_values = AIC_for_each_n_segments[node];
  this_AIC =  AICc_values[ best_fit_AICc[n_sigma_for_printing] ];
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this replaces two vectors, which have the starting and ending position of the best fit segments
// for a given sigma. This gets data from (most likeley) the get_data_from_best_fit_lines function
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::get_start_and_end_x_for_segments(vector<float>& start_x,
                        vector<float>& end_x, vector<int> seg_lengths)
{
  vector<float> empty_vec;
  start_x = empty_vec;
  end_x = empty_vec;
  int n_segs = seg_lengths.size();

  start_x.push_back(x_data[0]);
  end_x.push_back(seg_lengths[0]-1);

  for (int i = 1; i<n_segs; i++)
  {
    start_x.push_back( x_data[ seg_lengths[i-1] ]);
    end_x.push_back( x_data[ seg_lengths[i]-1 ] );
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function is used to calculate the slope, intercept, and likelihood of
// all possible linear segments along a series of data points.
// the function requires the data in x and y vectors, the maximum segment length
// and sigma, the standard deviation of the measured data. This will be approxamately
// the error in the surface elevation, although it might have to be increased simply because
// likelihood will tend to zero if this is too small. sigma should also be considered to
// contain the 'noise' inherent in channel incision so perhaps 1-5 metres is appropriate
// the maximum segment length is an integer: it is the number of data points used.
// these data points from raw chi data are irregularly spaced so two segments of the same
// 'length' can have different lengths in chi space. One remedey for this is a preprocessor that
// places the zeta vs chi data along evenly spaced points.
//
// The routine generates three matrices. The row of the matrix is the starting node of the segment.
// The column of the matrix is the ending node of the segment. Thus the routine will generate a
// matrix that is dimension n x n where n is the number of data points.
//
// One potential future development is to implement this using a sparse matrix from the boost mtl
// library to reduce the memory usage.
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::calculate_segment_matrices(float sigma)
{
  int n_data_points = x_data.size();
  if (minimum_segment_length>n_data_points)
  {
    //cout << "LSDMostLikelyPartitionsFinder::calculate_segment_matrices: your segment length is greater than half the number of data points" << endl;
    //cout << "This means that there can only be overlapping segments. Changing segment length to maximum segment length "<< endl;
    minimum_segment_length = n_data_points;
  }

  // set up the arrays
  // in the future I might consider using sparse arrays but for now we'll just populate
  // the empty spots with placeholders
  float no_data_value = -9999;
  Array2D<float> temp_array(n_data_points,n_data_points,no_data_value);
  like_array = temp_array.copy();
  m_array = temp_array.copy();
  b_array = temp_array.copy();
  rsquared_array = temp_array.copy();
  DW_array = temp_array.copy();

  int start_node = 0;
  int end_node = n_data_points-1;

  // populate the matrix.
  // the get segment row function is recursive so it moves down through all the possible
  // starting nodes
  //cout << "LINE 518, sigma is: " << sigma << endl;
  populate_segment_matrix(start_node, end_node, no_data_value, sigma);

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function populates the matrices of liklihood, m and b values
// it is a recursive algorithm so in fact it doesn't just get one row
// but drills down through all the possible starting nodes to complete the
// matrix
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::populate_segment_matrix(int start_node, int end_node, float no_data_value,float sigma)
{

  if (like_array[start_node][end_node] == no_data_value)
  {
    // create the two segments
    vector<float> segment_x;
    vector<float> segment_y;
    vector<float> residuals;
    vector<float> regression_results;
    float this_MLE;

    // now create iterators to deal with these segments
    vector<float>::iterator vec_iter_start;
    vector<float>::iterator vec_iter_end;

    // the first step is to get the segment starting on the
    // first node and ending on the last node
    // find out how many nodes are in the segment
    int n_nodes_in_segment = end_node - start_node+1;

    // resize the vectors accordingly
    segment_x.resize(n_nodes_in_segment);
    segment_y.resize(n_nodes_in_segment);

    // get the iterators for the start and end of the vectors
    vec_iter_start = x_data.begin()+start_node;
    vec_iter_end = vec_iter_start+n_nodes_in_segment;
    segment_x.assign(vec_iter_start,vec_iter_end);

    vec_iter_start = y_data.begin()+start_node;
    vec_iter_end = vec_iter_start+n_nodes_in_segment;
    segment_y.assign(vec_iter_start,vec_iter_end);

    // do the least squares regression on this segment
    // cout << "LINE 568, sigma is: " << sigma << endl;
    regression_results = simple_linear_regression(segment_x, segment_y, residuals);
    this_MLE = calculate_MLE_from_residuals( residuals, sigma);

    //cout << "LINE 584 doing start: " << start_node << " end: " << end_node << endl;

    like_array[start_node][end_node] = this_MLE;
    m_array[start_node][end_node] = regression_results[0];
    b_array[start_node][end_node] = regression_results[1];
    rsquared_array[start_node][end_node] = regression_results[2];
    DW_array[start_node][end_node] = regression_results[3];

    // now loop through all the end nodes that are allowed that are not the final node.
    // that is the first end node is first plus the maximum length -1 , and then
    // the final end node before the end of the data is the last node minus the
    // maximum length of the segment
    for (int loop_end = start_node+minimum_segment_length-1; loop_end< end_node-minimum_segment_length+1; loop_end++)
    {
      if (like_array[start_node][loop_end] == no_data_value)
      {
        // get this segment and resize the vectors
        n_nodes_in_segment = loop_end - start_node+1;
        segment_x.resize(n_nodes_in_segment);
        segment_y.resize(n_nodes_in_segment);

        // get the iterators for the start and end of the vectors
        vec_iter_start = x_data.begin()+start_node;
        vec_iter_end = vec_iter_start+n_nodes_in_segment;
        segment_x.assign(vec_iter_start,vec_iter_end);

        vec_iter_start = y_data.begin()+start_node;
        vec_iter_end = vec_iter_start+n_nodes_in_segment;
        segment_y.assign(vec_iter_start,vec_iter_end);

        // do the least squares regression on this segment
        regression_results = simple_linear_regression(segment_x, segment_y, residuals);
        this_MLE = calculate_MLE_from_residuals( residuals, sigma);

        // fill in the matrices
        like_array[start_node][loop_end] = this_MLE;
        m_array[start_node][loop_end] = regression_results[0];
        b_array[start_node][loop_end] = regression_results[1];
        rsquared_array[start_node][loop_end] = regression_results[2];
        DW_array[start_node][loop_end] = regression_results[3];
        //cout << "LINE 612 doing start: " << start_node << " end: " << loop_end << endl;

        // now get the row from the next segment
        populate_segment_matrix(loop_end+1, end_node, no_data_value,sigma);
      }
    }

  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function calculates the most likeley combination of segments given the liklihood
// of individual segments calculated by the calculate_segment_matrices function
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::find_max_like_of_segments()
{
  // first get the number of nodes
  int n_data_points = like_array.dim1();
  if (minimum_segment_length>n_data_points)
  {
    cout << "LSDStatsTools find_max_AIC_of_segments: your segment length is greater than the number of data points" << endl;
    cout << "This means that there can only be overlapping segments. Changing segment length to minimum segment length "<< endl;
    minimum_segment_length = n_data_points;
  }

  // get the maximum number of segments
  int max_n_segments = n_data_points/minimum_segment_length;

  // initialize a vector for holding the MLE of each n_segments
  vector<float> MLE_for_segments(max_n_segments,-0.00000000001);

  // initialize a vecvec for holding the MLE segment partition
  vector< vector <int> > most_likely_segments(max_n_segments);

  // create the partition data element
  partition_driver_to_vecvecvec(n_data_points);

  // now loop through the number of segments, calucalting the maximum likelihood each time
  vector< vector <int> > partition_vecvec;
  float this_MLE;
  int start_node,end_node;
  for (int n_elem = 0; n_elem< int(partitions.size()); n_elem++)
  {

    partition_vecvec = partitions[n_elem];
    int n_partitions_this_nsegments = partition_vecvec.size();
    // cout << "n_segments: " << n_elem+1 << " and number of partitions of this n segments: " << n_partitions_this_nsegments << endl;
    for (int n_partition = 0; n_partition< n_partitions_this_nsegments; n_partition++)
    {
      // cout << "number of partitions: " << n_partition+1 << " of " << n_partitions_this_nsegments << endl;
      vector<int> individual_partition = partition_vecvec[n_partition];
      int n_elements = individual_partition.size();

      do
      {
        // calcualte the MLE for this particular permutation
        this_MLE = 1;
        start_node = 0;
        for (int i = 0; i<n_elements; i++)
        {
          end_node = start_node+individual_partition[i]-1;
          //cout << "start node: " << start_node << " " << " end node: " << end_node
          //     << " and like: " << like_array[start_node][end_node] << endl;
          this_MLE = this_MLE*like_array[start_node][end_node];
          start_node = end_node+1;
        }
        //cout << "This MLE: " << this_MLE << " seg MLE: " << MLE_for_segments[n_elem] << endl;

        // now test if this MLE is better than the the best MLE so far for this number of
        // partitions
        if( this_MLE > MLE_for_segments[n_elem] )
        {
          //cout << "element is: " << n_elem << " new MLE" << endl;
          MLE_for_segments[n_elem] = this_MLE;
          most_likely_segments[n_elem] = individual_partition;
        }
      } while ( prev_permutation(individual_partition.begin(),individual_partition.end()) );

    }
  }
  segments_for_each_n_segments = most_likely_segments;
  MLE_of_segments = MLE_for_segments;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function drives the partitioning algorithms
// k is the number of elements in the partition
// minimum length is the minimum length of the segment
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::partition_driver_to_vecvecvec(int k)
{
  int n = 2*k;
  int t = 0;
  vector<int> p(k,0);

  int max_segments = k/minimum_segment_length;
  vector< vector < vector<int> > > this_partition(max_segments);
  partitions = this_partition;

  // run the partitioning code
  //cout << "partition_driver_to_vecvecvec, doing partitions" << endl;
  partitions_with_minimum_length( n, k, t, p);
  //cout << "partition_driver_to_vecvecvec, finished partitions" << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// an integer partition algorithm
// Algorithm and original Pascal implementation: Frank Ruskey, 1995.
// Translation to C: Joe Sawada, 1997
// grabbed from http://theory.cs.uvic.ca/inf/nump/NumPartition.html
// adapted smm 21/12/2012
// algorith described in
// http://mathworld.wolfram.com/PartitionFunctionP.html
// and
// Skiena, S. Implementing Discrete Mathematics: Combinatorics and Graph Theory with Mathematica. Reading, MA: Addison-Wesley, 1990.
// this is a further adaptation that only presents solution to the partition
// with segments of a minimum length
// it stores all the partitions
// http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5123/Permutations-in-C.htm
// http://www.cplusplus.com/reference/algorithm/next_permutation/
// http://mdm4u1.wetpaint.com/page/4.3+Permutations+with+Some+Identical+Elements
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::partitions_with_minimum_length(int n, int k, int t, vector<int>& p)
{

  int j;

  p[t] = k;
  if (n==k)
  {
    partition_assign(t, p);
  }
  for (j=LSDpartitions_min(k,n-k); j>=minimum_segment_length; j--)
  {
    partitions_with_minimum_length(n-k,j,t+1,p);
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// a function for use with the permutations
// this assigns values into the vecvecvec that contains all the partitioning information
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::partition_assign(int t, vector<int>& p)
{

  vector< vector<int> > this_nsegments_vecvec = partitions[t-1];
  vector<int> this_partitions_partitions;

  for(int i=1; i<=t; i++)
  {
    this_partitions_partitions.push_back( p[i] );
  }
  this_nsegments_vecvec.push_back(this_partitions_partitions);
  partitions[t-1] = this_nsegments_vecvec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// a function for use with the permutations
// gets the mininum of two values
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDMostLikelyPartitionsFinder::LSDpartitions_min( int x, int y)
{
  if (x<y)
  {
    return x;
  }
  else
  {
    return y;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this takes a likelihood array that has been calcualted with a given sigma value and
// normalizes the sigma values as though sigma was equal to 1.
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Array2D<float> LSDMostLikelyPartitionsFinder::normalize_like_matrix_to_sigma_one(float sigma)
{
  // get dimensions of the like array
  int nRows = like_array.dim1();
  int nCols = like_array.dim2();
  float no_data_value = -9999;
  float sigsquared = sigma*sigma;
  Array2D<float> sig1_like_array = like_array.copy();

  for (int row = 0; row<nRows; row++)
  {
    for(int col = 0; col<nCols; col++)
    {
      if(sig1_like_array[row][col] != no_data_value)
      {
        sig1_like_array[row][col] = pow(sig1_like_array[row][col],sigsquared);
      }
    }
  }

  return sig1_like_array;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this normalizes but with vector data, for use with MLE vector for segments
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDMostLikelyPartitionsFinder::normalize_like_vector_to_sigma_one(float sigma, vector<float> like_vector)
{
  // get dimensions of the like vector
  int ndata = like_vector.size();
  float sigsquared = sigma*sigma;
  vector<float> sig1_like_vector(ndata);

  for (int i = 0; i<ndata; i++)
  {

    sig1_like_vector[i] = pow(like_vector[i],sigsquared);

  }

  return sig1_like_vector;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this takes a normalize likelihood array and updates the values to a new sigma value
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::change_normalized_like_matrix_to_new_sigma(float sigma, Array2D<float>& sig1_like_array)
{
  // get dimensions of the like array
  int nRows = sig1_like_array.dim1();
  int nCols = sig1_like_array.dim2();
  float no_data_value = -9999;
  float one_over_sigsquared = 1/(sigma*sigma);
  like_array = sig1_like_array.copy();

  for (int row = 0; row<nRows; row++)
  {
    for(int col = 0; col<nCols; col++)
    {
      if(like_array[row][col] != no_data_value)
      {
        like_array[row][col] = pow(like_array[row][col],one_over_sigsquared);
      }
    }
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this does the same as above except for vector data
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDMostLikelyPartitionsFinder::change_normalized_like_vector_to_new_sigma(float sigma, vector<float> sig1_like_vector)
{
  // get dimensions of the like vector
  int ndata = sig1_like_vector.size();
  float one_over_sigsquared = 1/(sigma*sigma);
  vector<float> like_vector(ndata);

  for (int i = 0; i<ndata; i++)
  {
    like_vector[i] = pow(sig1_like_vector[i],one_over_sigsquared);
  }

  return like_vector;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this does the same as above except for vector data
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> LSDMostLikelyPartitionsFinder::transform_like_from_sigma1_to_sigma2(float sigma1,
                          vector<float> sig1_like_vector, float sigma2)
{
  // get dimensions of the like vector
  int ndata = sig1_like_vector.size();
  float sigsq_transform = (sigma1*sigma1)/(sigma2*sigma2);
  vector<float> like_vector(ndata);

  for (int i = 0; i<ndata; i++)
  {
    like_vector[i] = pow(sig1_like_vector[i],sigsq_transform);
  }

  return like_vector;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes the normalized MLE values (normalized with sigma = 1) and returns the
// best fit number of segments from both the AIC and the AICc measures. It also returns
// two vector of vectors which are the AIC values for the various values of sigma
// passed to the function in the sigma values vector
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::get_n_segments_for_various_sigma(vector<float> sigma_values)
{
  int n_sigma = sigma_values.size();
  vector<float> empty_vec;
  vector<float> AIC_of_segments;
  vector<float> AICc_of_segments;
  vector< vector<float> > AIC_for_each(n_sigma);
  vector< vector<float> > AICc_for_each(n_sigma);
  vector<int> bf_AIC(n_sigma);
  vector<int> bf_AICc(n_sigma);

  // loop through the sigma values collecting the AIC and AICc values
  for (int i = 0; i< n_sigma; i++)
  {
    // calculate the AIC values for this value of sigma
    calculate_AIC_of_segments_with_variable_sigma(sigma_values[i],
                AIC_of_segments,AICc_of_segments);
    AIC_for_each[i] = AIC_of_segments;
    AICc_for_each[i] = AICc_of_segments;

    // now find the minimum AIC and AICc
    float minimum_AIC = 10000;
    int min_AIC_segments = 0;
    float minimum_AICc = 10000;
    int min_AICc_segments = 0;

    int n_AIC = AIC_of_segments.size();
    for (int n_seg = 0; n_seg<n_AIC; n_seg++)
    {
      if(AIC_of_segments[n_seg] < minimum_AIC)
      {
        minimum_AIC = AIC_of_segments[n_seg];
        min_AIC_segments = n_seg;
      }

      if(AICc_of_segments[n_seg] < minimum_AICc)
      {
        minimum_AICc = AICc_of_segments[n_seg];
        min_AICc_segments = n_seg;
      }
    }

    bf_AIC[i] = min_AIC_segments;
    bf_AICc[i] = min_AICc_segments;
  }

  // replace vectors
  AIC_for_each_n_segments = AIC_for_each;
  AICc_for_each_n_segments = AICc_for_each;
  best_fit_AIC = bf_AIC;
  best_fit_AICc = bf_AICc;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function calculates AIC and AICc of segments taking the maximum_MLE based on a sigma of one
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::calculate_AIC_of_segments_with_variable_sigma(float sigma,
                vector<float>& AIC_of_segments,vector<float>& AICc_of_segments)
{

  // recast the MLE_vector
  vector<float> new_sig_MLE = transform_like_from_sigma1_to_sigma2(base_sigma, MLE_of_segments, sigma);

  // initialize the vector holding the Aikake Information Criterion
  // and then calcualte the AIC and the AICc
  vector<float> AIC(MLE_of_segments.size(),0.0);
  vector<float> AICc(MLE_of_segments.size(),0.0);
  for (int n_elem = 0; n_elem< int(MLE_of_segments.size()); n_elem++)
  {
    float AICk = (float(n_elem)+1);           // this is the number of segments
    float AICn = float(x_data.size());        // this is the number of data elements
    // if the MLE is 0 or less, this will throw an error. This only happens if the fit
    // is terrible so in this case set AIC and AICc to a large positive number
    if(new_sig_MLE[n_elem]<= 0)
    {
      AIC[n_elem] = 9999;
      AICc[n_elem] = 9999;
    }
    else
    {
      //cout << "n_segs: " << n_elem+1
      //   << " MLE: " <<  new_sig_MLE[n_elem] << " log: "
      //   << log( new_sig_MLE[n_elem]) << " 2nd term: "
      //   << -2*log( new_sig_MLE[n_elem]) << endl;
      AIC[n_elem] = 4*AICk-2*log( new_sig_MLE[n_elem]);    // the 4 comes from the fact that
                                                             // for each segment there are 2 parameters
      AICc[n_elem] = AIC[n_elem] + 2*AICk*(AICk+1)/(AICn-AICk-1);
    }
    //cout << "AIC: " << AIC[n_elem] << " and AICc: " << AICc[n_elem] << endl << endl;

  }
  AIC_of_segments = AIC;
  AICc_of_segments = AICc;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function returns the m, b, r^2 and D-W values for the best fit segments
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::get_properties_of_best_fit_segments(int bestfit_segments_node,
                     vector<float>& m_values, vector<float>& b_values,
                     vector<float>& r2_values, vector<float>& DW_values)
{


  int n_segments = bestfit_segments_node+1;     // this converts from 0-based indexing to
                          // the actual number of segments

  //cout << "Line 1249, n_segments: " << n_segments << endl;

  // initialize temp_vectors
  vector<float> m(n_segments);
  vector<float> b(n_segments);
  vector<float> r2(n_segments);
  vector<float> DW(n_segments);

  // the start node and end node, used to index into the arrays
  int start_node,end_node;

  // get the segment lengths
  vector<int> individual_partition = segments_for_each_n_segments[bestfit_segments_node];

  //cout << "Line 1263 segments_for_each_n_segments.size(): "  << segments_for_each_n_segments.size() << endl;
  //cout << "individual_partition.size(): " << individual_partition.size() << endl;

  // now loop through the segments
  start_node = 0;
  for (int i = 0; i<n_segments; i++)
  {
    end_node = start_node+individual_partition[i]-1;
    //cout << "start node: " << start_node << " " << " end node: " << end_node
    //     << " m: " << m_array[start_node][end_node] << " b: " << b_array[start_node][end_node]
    //     << " r^2: " << rsquared_array[start_node][end_node]
    //     << " DW: " << DW_array[start_node][end_node] << endl;
    m[i] = m_array[start_node][end_node];
    b[i] = b_array[start_node][end_node];
    r2[i] = rsquared_array[start_node][end_node];
    DW[i] = DW_array[start_node][end_node];
    start_node = end_node+1;
  }

  m_values = m;
  b_values = b;
  r2_values = r2;
  DW_values = DW;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function prints the most likeley segment lengths to screen
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::print_to_screen_most_likeley_segment_lengths()
{

  // loop through the number of segments, printing the
  // most likeley segment lengths to screen
  cout << endl << "printing most likeley segment lenghts: " << endl;
  for (int n_elem = 0; n_elem< int(segments_for_each_n_segments.size()); n_elem++)
  {
    cout << "n elements: " << n_elem << " and MLE: " << MLE_of_segments[n_elem] << endl;
    vector<int> individual_partition = segments_for_each_n_segments[n_elem];
    cout << "segment lengths: " << individual_partition.size() << endl;
    if ( int(individual_partition.size()) != n_elem+1)
    {
      cout << "LINE 707 statstools something is wrong n partitions is incorrect" << endl;
    }

    for (int i = 0; i<=n_elem; i++)
    {
      cout << individual_partition[i] << " ";
    }
  }
  cout << endl << "finished printing most likeley segment lengths" << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function prints the AIC and AICc values to screen
//
// SMM 01/02/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDMostLikelyPartitionsFinder::print_AIC_and_AICc_to_screen(vector<float> sigma_values)
{
  // get the number of sigma values
  int n_sigs = sigma_values.size();
  vector<float> AI_values;
  vector<int> this_minimum;
  int AI_sz;
  int n_mins;


  cout << endl << "LSDMostLikelyPartitionsFinder::print_AIC_and_AICc_to_screen" << endl;

  // loop through sigma values printing the best fits
  for (int i = 0; i<n_sigs; i++)
  {
    cout << endl << endl << "sigma is: " << sigma_values[i] << endl;
    AI_values = AIC_for_each_n_segments[i];
    AI_sz = int(AI_values.size());
    cout << "Min AIC node is: " << best_fit_AIC[i] << " and AIC values: ";
    for (int j = 0; j<AI_sz; j++)
    {
      cout << AI_values[j] << " ";
    }
    cout << endl;

    cout << "the segments lengths are: ";
    this_minimum = segments_for_each_n_segments[ best_fit_AIC[i] ];
    n_mins = best_fit_AIC[i]+1;
    for (int seg = 0; seg<n_mins; seg++)
    {
      cout << this_minimum[seg] << " ";
    }
    cout << endl;

    cout << "Min AICc node is: " << best_fit_AICc[i] << " and AICc values: ";
    AI_values = AICc_for_each_n_segments[i];
    for (int j = 0; j<AI_sz; j++)
    {
      cout << AI_values[j] << " ";
    }
    cout << endl;

    cout << "the segments lengths are: ";
    this_minimum = segments_for_each_n_segments[ best_fit_AICc[i] ];
    n_mins = best_fit_AICc[i]+1;
    for (int seg = 0; seg<n_mins; seg++)
    {
      cout << this_minimum[seg] << " ";
    }
    cout << endl;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#endif
