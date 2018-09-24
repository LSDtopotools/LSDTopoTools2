//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDPorewaterColumn
// Land Surface Dynamics PorewterColumn object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic tools
//  This object calculates porewater pressure based on the Iverson 2000 WRR
// model
//
// Developed by:
//  Simon M. Mudd
//  Stuart W.D. Grieve
//
// Copyright (C) 2016 Simon M. Mudd 2013 6
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


/** @file LSDPorewaterColumn.hpp
@author Simon M. Mudd, University of Edinburgh
@author Stuart W. D. Grieve, University of Edinburgh

**/



#ifndef LSDPorewaterColumn_CPP
#define LSDPorewaterColumn_CPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "LSDPorewaterColumn.hpp"
#include "LSDPorewaterParams.hpp"
#include "LSDStatsTools.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Empty create function
// Starts with some defaults. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterColumn::create()
{
  row = 0;
  col = 0;
  node_index = 0;
  cout << "I am an empty LSDPorewaterColumn object." <<endl;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function just uses an initial Psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterColumn::create(vector<float> Initial_Psi)
{
  row = 0;
  col = 0;
  node_index = 0;
  Psi = Initial_Psi;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function takes a porewater parameter object and uses the
// steady infiltration rate to set psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void  LSDPorewaterColumn::create(LSDPorewaterParams LSDPP)
{
  row = 0;
  col = 0;
  node_index = 0;
  Psi = LSDPP.calculate_steady_psi();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the Psi values to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterColumn::PrintPsiToScreen()
{
  for(int i = 0; i<int(Psi.size()); i++)
  {
    cout << "Psi["<<i<<"]: " << Psi[i] << endl;
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This caluclates the response function
// THis comes from iverson's equation 27e
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDPorewaterColumn::CalculateResponseFunction(float t_star)
{
    float R;
    
    float sqrt_term = sqrt(t_star/M_PI);
    float exp_term = exp(-1/t_star);
    
    float multiple_bit = sqrt_term*exp_term;
    
    if (t_star != 0)
    {
      R = multiple_bit- erfcf(1/ (sqrt(t_star)));
    }
    else   // If t_star is 0, then 1/sqrt(t_star) is infinity, meaning erfc is 0)
    {
      R = 0;
    }

    return R;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This caluclates the Psi value based on iverson's equation 27
// Has only the transient component of psi
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::CalcualtePsiTransient(LSDPorewaterParams& LSDPP, float Iz_over_Kz, 
                           float t_star, float T_star)
{
    vector<float> Depths = LSDPP.get_Depths();
    vector<float> transient_Psi(Depths.size());
    
    float R;
    if (t_star < T_star)
    {
      R = CalculateResponseFunction(t_star);
    }
    else
    {
      R = CalculateResponseFunction(t_star-T_star);
    }
    
    // This solves the equation, based on the response function (R_fn),
    // which is equation 27e
    for (int i = 0; i< int(Depths.size()) ; i++)
    {
      transient_Psi[i] = Depths[i]*Iz_over_Kz*R;
    }

    return transient_Psi;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Compute psi from equation 27a and b, but using dimensional time
// A bit slow since I haven't vectorised the calculations.
// Only calculates the transient component of psi for use with 
// time series of rainfall
// times need to be in seconds
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::CalculatePsiDimensionalTimeTransient(LSDPorewaterParams& LSDPP, float t, float T, float Iz_over_Kz)
{

  vector<float> Depths = LSDPP.get_Depths();
  vector<float> transient_Psi(Depths.size());
  float t_star;
  float T_star;
  float zsquare;
  
  float D_hat = LSDPP.get_D_hat();
  float R;
  
  // loop through depths: each depth has a different t_star and T_star since
  // these depend on depth
  for(int i = 0; i< int(Depths.size()) ; i++)
  {
    // first get the nondimensional time. Note that according to
    // equations 27c,d the dimensionless time is a function of depth,
    // so each point below the surface has a different t_star and T_star
    zsquare = Depths[i]*Depths[i];
    t_star = t * D_hat / zsquare;
    T_star = T * D_hat / zsquare;
    
    if (t_star < T_star)
    {
      R = CalculateResponseFunction(t_star);
    }
    else
    {
      R = CalculateResponseFunction(t_star)-CalculateResponseFunction(t_star-T_star);
    }
    transient_Psi[i] =Depths[i]*Iz_over_Kz*R;
    
    //cout << "depth: " << Depths[i] << " t_star: " << t_star << " T_star: " << T_star << endl;
    
    
  }
  return transient_Psi;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This caluclates the Psi value based on iverson's equation 27
// It parses a time series
// The durations and time are in seconds. 
// The intensities are in Iz_over_Kz
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterColumn::CalculatePsiFromTimeSeries(vector<float> durations, vector<float> intensities, 
                                LSDPorewaterParams& LSDPP, float t)
{
  // Get the steady state time
  vector<float> steady_psi = LSDPP.calculate_steady_psi();
  vector<float> cumulative_psi = steady_psi;
  
  // see what the result is:
  //cout << endl << "======================" << endl << "Steady psi: " << endl;
  //for(int i = 0; i<int(cumulative_psi.size()); i++)
  //{
  //  cout << "Psi["<<i<<"]: " << cumulative_psi[i] << endl;
  //}
  //cout << "======================" << endl << endl << endl;

  // Now we try to construct the transient pressure. 
  // loop through the record getting cumulative times
  vector<float> starting_times;
  starting_times.push_back(0);
  float cumulative_time = 0;
  int count = 0; 
  bool end_count_found = false;
  int end_count = 0;
  
  for (int i = 0; i< int(durations.size()); i++)
  {
    cumulative_time += durations[i];
    
    //cout << "t: " << t << " cumulative time: " << cumulative_time << endl;
    
    // the cumulative time is the time at the end of this timestep. 
    // if the cumulative time  is less than the time of simulation, 
    // then we need to acount for this pulse of rainfall        
    if (t < cumulative_time)
    {
      if (end_count_found == false)
      {
        end_count_found = true;
        end_count = count;
      }
    }
    count++;
    starting_times.push_back(cumulative_time);
  }
  
  // we don't need the last element
  starting_times.pop_back();
  
  //cout << "N starting times: " <<  starting_times.size() << endl;

  // If we didn't find the end count it means the rainfall records have ended and we need
  // all of the data        
  if (end_count_found == false)
  {
    // The minus one is needed since we have counted past the end of the index
    end_count = count-1;
  }
  
  //cout << "end count is: " << end_count << endl;


  // check starting times, etc
  //for(int i = 0; i< int(starting_times.size()); i++)
  //{
  //  cout << "st: " << starting_times[i] << " i: " << intensities[i] << " d: " << durations[i] << endl;
  //
  //}


  // okay, now get the transients from superposition 
  // First we need to figure out how many of these we will need
  float eff_t, this_intensity, this_duration;
  vector<float> this_transient_Psi;
  for(int i = 0; i< int(starting_times.size()); i++)
  {
    if(i<= end_count)
    {
      eff_t = t-starting_times[i];
      this_intensity = intensities[i];
      this_duration = durations[i];
      
      //cout << "Eff t: " << eff_t << " and intensity: " << this_intensity << " dur: " << this_duration << endl;
      
      // get this steps Psi value
      this_transient_Psi = CalculatePsiDimensionalTimeTransient(LSDPP, eff_t, this_duration, this_intensity);
      
      // check values
      //cout << "Transient psi is:"<< endl;
      //for(int i = 0; i<int(cumulative_psi.size()); i++)
      //{
      //  cout << this_transient_Psi[i] << endl;
      //}

      // add this step's transient Psi values.
      for(int i = 0; i<int(cumulative_psi.size()); i++)
      {
        cumulative_psi[i]+=this_transient_Psi[i];
      }
    }
  }
  
  // see what the result is:
  for(int i = 0; i<int(cumulative_psi.size()); i++)
  {
    cout << "Psi["<<i<<"]: " << cumulative_psi[i] << endl;
  }
  
  Psi = cumulative_psi;
  
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Factor of safety calculations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is the friction
float LSDPorewaterColumn::F_f(LSDPorewaterParams& LSDPP)
{
  float alpha =  LSDPP.get_alpha();
  float friction_angle = LSDPP.get_friction_angle();
  
  float tan_alpha = tan(alpha);
  float tan_friction_angle = tan(friction_angle);
  
  return tan_friction_angle/tan_alpha;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is from the cohesion
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::F_c(LSDPorewaterParams& LSDPP)
{
  float cohesion = LSDPP.get_cohesion();
  float alpha =  LSDPP.get_alpha();
  float weight_of_soil = LSDPP.get_weight_of_soil();
  vector<float> Depths = LSDPP.get_Depths();
  
  float denom;
  float denom2 = sin(alpha)*cos(alpha);
  
  vector<float> F_c_vec;
  
  for(int i = 0; i< int(Depths.size()); i++)
  {
    denom = Depths[i]*weight_of_soil;
    
    F_c_vec.push_back( cohesion/(denom*denom2) );
  }
  
  return F_c_vec;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is from the water
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::F_w(LSDPorewaterParams& LSDPP)
{
  vector<float> Depths = LSDPP.get_Depths();
  float alpha =  LSDPP.get_alpha();
  float friction_angle = LSDPP.get_friction_angle();
  float weight_of_soil = LSDPP.get_weight_of_soil();
  float weight_of_water = LSDPP.get_weight_of_water();
    
  float denom, num1, num2;
  float denom2 = sin(alpha)*cos(alpha);
  float denom_tot;
  
  vector<float> F_w_vec;
  
  for(int i = 0; i< int(Depths.size()); i++)
  {
    
    num1 = Psi[i]*weight_of_water;
    num2 = -num1*tan(friction_angle);
    
    denom = Depths[i]*weight_of_soil;
    denom_tot = denom*denom2;
    
    F_w_vec.push_back( num2/denom_tot );
  }
  
  return F_w_vec;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is the total factor of safety
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDPorewaterColumn::FS(LSDPorewaterParams& LSDPP)
{

  // get the components of the factor of safety
  float F_f_float = F_f(LSDPP);
  vector<float> F_c_vec = F_c(LSDPP);
  vector<float> F_w_vec = F_w(LSDPP);

  // get the 
  vector<float> FS = F_c_vec;
  for(int i = 0; i< int(F_c_vec.size()); i++)
  {
    FS[i] = F_f_float+F_c_vec[i]+F_w_vec[i];
    cout << "FS["<<i<<"]: " << FS[i] << endl;
  }
  
  return FS;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks to see a failure depth
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDPorewaterColumn::DepthOfFailure(LSDPorewaterParams& LSDPP, float minimum_depth)
{
  float depth_of_failure = -9999;
  
  // get the factor of safety vector
  vector<float> FoS = FS(LSDPP);
  
  // get the depth vector
  vector<float> Depths = LSDPP.get_Depths();
  
  int N_depths = int(Depths.size());
  for(int i = 0; i< N_depths; i++)
  {
    // only check FS if above minimum depth
    if(Depths[i]>= minimum_depth)
    {
      if(FoS[i] < 1.0)
      {
        depth_of_failure = Depths[i];
        i = N_depths;
      }
    }
  }
  return depth_of_failure;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This checks to see a failure depth
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterColumn::GetMinFS(LSDPorewaterParams& LSDPP, float minimum_depth, float& depth_of_minFS, float& minFS)
{
  depth_of_minFS = minimum_depth;
  float min_FS = 9999;
  
  // get the factor of safety vector
  vector<float> FoS = FS(LSDPP);
  
  // get the depth vector
  vector<float> Depths = LSDPP.get_Depths();
  
  int N_depths = int(Depths.size());
  for(int i = 0; i< N_depths; i++)
  {
    // only check FS if above minimum depth
    if(Depths[i]>= minimum_depth)
    {
      if(FoS[i] < min_FS)
      {
        depth_of_minFS = Depths[i];
        min_FS  = FoS[i];
      }
    }
  }
  
  minFS = min_FS;
  
  
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This scans a timeseries for a failure
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDPorewaterColumn::ScanTimeseriesForFailure(vector<float> durations, vector<float> intensities,
                                   LSDPorewaterParams& LSDPP, float minimum_depth, 
                                   vector<float> times)
{
  // loop through times
  int n_times = int(times.size());
  float depth_of_minFS;
  float min_FS;
  for(int i = 0; i< n_times; i++)
  {
    // get the pore pressure
    CalculatePsiFromTimeSeries(durations, intensities, LSDPP, times[i]);
    
    // get the min
    GetMinFS(LSDPP, minimum_depth, depth_of_minFS, min_FS);
    
    cout << "Time in weeks is: " << LSDPP.seconds_to_weeks(times[i]) << " d min FS: " 
         << depth_of_minFS << " min FS: " << min_FS << endl;

  }
}

#endif