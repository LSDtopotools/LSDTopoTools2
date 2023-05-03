//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDParticleColumn.cpp
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for holding particles in a column
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

#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include "LSDStatsTools.hpp"
#include "LSDCRNParameters.hpp"
#include "LSDParticle.hpp"
#include "LSDParticleColumn.hpp"
using namespace std;

#ifndef LSDParticleColumn_CPP
#define LSDParticleColumn_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  create, the default constructor
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::create()
{
   Row = 0;
   Col = 0;
   NodeIndex = 0;
   SoilThickness = 0;
   
   RockDensity = 2650;
   SoilDensity = 1300;
   
   UseDenstyProfile = false;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This create function uses all data members
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::create(int tRow, int tCol, int tNodeIndex, 
        double tSoilDensity, double tRockDensity,   
        float tSoilThickness, float tDataResolution, bool tUseDenstyProfile, 
        vector<double> tDensityDepths, vector<double> tDensityDensities, 
        list<LSDCRNParticle> tCRNParticleList)
{
  Row = tRow;
  Col = tCol;
  NodeIndex = tNodeIndex; 
  SoilDensity = tSoilDensity; 
  RockDensity = tRockDensity;
  SoilThickness = tSoilThickness;
  DataResolution = tDataResolution;
  UseDenstyProfile = tUseDenstyProfile; 
  DensityDepths = tDensityDepths;
  DensityDensities = tDensityDensities;
  CRNParticleList = tCRNParticleList;
}        
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The const copy constructor for a const object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDParticleColumn& LSDParticleColumn::operator=(const LSDParticleColumn& rhs)
{
	if (&rhs != this)
    {
		  create(rhs.getRow(), rhs.getCol(), rhs.getNodeIndex(), 
               rhs.getSoilDensity(),rhs.getRockDensity(),
               rhs.getSoilThickness(),rhs.getDataResolution(),
               rhs.getUseDenstyProfile(),rhs.getDensityDepths(),
               rhs.getDensityDensities(), rhs.getCRNParticleList());
    }
    return *this;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The const copy constructor  for a non const object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDParticleColumn& LSDParticleColumn::operator=(LSDParticleColumn& rhs)
{
	if (&rhs != this)
    {
		  create(rhs.getRow(), rhs.getCol(), rhs.getNodeIndex(), 
               rhs.getSoilDensity(),rhs.getRockDensity(),
               rhs.getSoilThickness(),rhs.getDataResolution(),
               rhs.getUseDenstyProfile(),rhs.getDensityDepths(),
               rhs.getDensityDensities(), rhs.getCRNParticleList());
    }
    return *this;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-






//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This resets the scaling of the CRN particle which is required since the 
// snow and self shielding alter the scaling values
// The scaling is updated within the object which is passed by reference
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::reset_scaling(LSDCRNParameters& LSDCRNP, string Muon_scaling)
{
  if (Muon_scaling == "Schaller" )
  {
    LSDCRNP.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    LSDCRNP.set_Braucher_parameters();
  }
  else if (Muon_scaling == "BraucherBorchers" )
  {
    LSDCRNP.set_BraucherBorchers_parameters();
  }
  else if (Muon_scaling == "Granger" )
  {
    LSDCRNP.set_Granger_parameters();
  }
  else if (Muon_scaling == "newCRONUS" )
  {
    LSDCRNP.set_newCRONUS_parameters();
  }
  else
  {
    cout << "You didn't set the muon scaling." << endl
         << "Options are Schaller, Braucher, newCRONUS, BraucherBorchers, and Granger." << endl
         << "You chose: " << Muon_scaling << endl
         << "Defaulting to Braucher et al (2009) scaling" << endl;
    LSDCRNP.set_Braucher_parameters();
  }
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This creates a column with particles at constant depth spacing and
// wit SS cosmo concentration
// SMM 25 July 2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::initiate_SS_cosmo_column_3CRN(int start_type, 
          double startxLoc, double startyLoc,
		      double start_depth, double particle_spacing, 
		      double zeta, double eff_eros_rate,
		      LSDCRNParameters& CRN_param)
{
  // get number of particles
  int N_particles = double(start_depth/particle_spacing)+1;

  // create the list
  //list<LSDCRNParticle> CRN_plist;

  double this_depth;

  // now loop over the depths, inserting particles and setting them to steady state
  double top_depth = start_depth-double(N_particles-1)*particle_spacing;
  //cout << "Inserting particles, top depth is: " << top_depth << endl;
  //cout << "N_particles are: "  << N_particles << endl;
  for (int p = 0; p< N_particles; p++)
  {
    // first the depth
    this_depth = top_depth+double(p)*particle_spacing;

    // now insert the particle into the list
    insert_particle_into_column(start_type, startxLoc, startyLoc,
			      this_depth,zeta);
  }

  // now loop through the particles, setting to steady
  list<LSDCRNParticle>::iterator part_iter;
  part_iter = CRNParticleList.begin();
  while(part_iter != CRNParticleList.end())
  {
    // update the CRN_concntrations
    ( *part_iter ).update_10Be_SSfull(eff_eros_rate,CRN_param);
    ( *part_iter ).update_14C_SSfull(eff_eros_rate,CRN_param);
    ( *part_iter ).update_21Ne_SSfull(eff_eros_rate,CRN_param);
 
    part_iter++;   
  }
  
  // check the apparent erosion
  //cout << "eff eros: " << eff_eros_rate << " and eros is: " << 10*eff_eros_rate/RockDensity << endl;
  //part_iter = CRNParticleList.begin();
  //while(part_iter != CRNParticleList.end())
  //{
  //  // update the CRN_concntrations
  //  cout << " app eros of part: " << (*part_iter).apparent_erosion_10Be_neutron_only(RockDensity, CRN_param); 
  //  part_iter++;   
  //}  

  //CRNParticleList = CRN_plist;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 




void LSDParticleColumn::read_erosion_time_series(string filename)
{
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  vector<float> this_time_vec;
  vector<float> this_rate_vec;
  vector<float> this_removal_vec;
  
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv data file, but the file" << filename
         << " doesn't exist;  LSDSpatialCSVReader::load_csv_data" << endl;
    exit(EXIT_FAILURE);
  }

  // get the headers from the first line
  getline(ifs, line_from_file);

  // now loop through the rest of the lines, getting the data.
  while( getline(ifs, line_from_file))
  {
    //cout << "Getting line, it is: " << line_from_file << endl;
    // reset the string vec
    this_string_vec = empty_string_vec;

    // create a stringstream
    stringstream ss(line_from_file);

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

    //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
    if ( int(this_string_vec.size()) <= 0)
    {
      cout << "Hey there, I am trying to load your csv data but you seem not to have" << endl;
      cout << "enough columns in your file. I am ignoring a line" << endl;
    }
    else if (int(this_string_vec.size()) == 3)
    {
      this_time_vec.push_back( atof(this_string_vec[0].c_str()) );
      this_rate_vec.push_back( atof(this_string_vec[1].c_str()) );
      this_removal_vec.push_back( atof(this_string_vec[2].c_str()) );

    }
    else
    {
      cout << "I am trying to read an erosion/deposition history but you have the wrong number of columns." << endl;
      cout << "There should be three columns: time,rate,removal" << endl;
      exit(0);
    }

  }

  rate_vec = this_rate_vec;
  time_vec = this_time_vec;
  removal_vec = this_removal_vec;

  check_erosion_time_series();

}


void LSDParticleColumn::check_erosion_time_series()
{
  int n_inputs = int(rate_vec.size());
  if (n_inputs == 0)
  {
    cout << "There does not seem to be any time series information. Exiting. " << endl;
    exit(0);
  }
  else if (time_vec[n_inputs-1] != -9999)
  {
    cout << "Error, the last time must be -9999" << endl;
    cout << "This indicates the long-term erosion rate." << endl;
    cout << "Check your time series and ensure you have a long term erosion rate." << endl;
  }
  else if (rate_vec[n_inputs-1] <= 0)
  {
    cout << "The final row must contain an erosion rate greater than zero which is the background erosion rate." << endl;
    cout << "This is not the case here so exiting." << endl;
    exit(0);
  }
  else
  {
    // first do the first step. Check if this is a removal or erosion step
    if (rate_vec[0] == 0)
    {
      if (time_vec[0] != 0)
      {
        cout << "You cannot have the first removal at a nonzero time, exiting." << endl;
        exit(0);
      }
    }
    else
    {
      if (removal_vec[0] != 0)
      {
        cout << "You cannot have a nonzero removal and a nonzero erosion at the same time, exiting." << endl;
        exit(0);
      }
    }

    // now loop through the rest
    for (int i = 1; i<n_inputs-1; i++)
    {
      if (removal_vec[i] != 0 && rate_vec[i] != 0)
      {
        cout << "You cannot have a nonzero removal and a nonzero erosion at the same time, exiting." << endl;
        exit(0);
      }
      if (removal_vec[i] != 0)
      {
        if (time_vec[i] != time_vec[i-1])
        {
          cout << "I have a removal but it must be at the same time the previous step ended. Exiting." << endl;
          exit(0);
        }
      }
    }
  }
}

// This just reads the erosion timeseries and calculates how much depth must be added
// to te sampling location to get the start of the time series. 
float LSDParticleColumn::calculate_increase_in_depth_from_erosion_timeseries()
{
  int n_inputs = int(time_vec.size());

  float d_change = 0;

  // the first input must be a rate 
  d_change+= rate_vec[0]*time_vec[0];

  // Loop through the inputs, removing/depositing and eroding as you go along. 
  // we ignore the last input since that is the background rate
  for (int i = 1; i<n_inputs-1; i++)
  {
    if (rate_vec[i] != 0)
    {
      d_change+= rate_vec[i]*(time_vec[i]-time_vec[i-1]);
    }
    else
    {
      d_change+= removal_vec[i];
    }
  }

  return d_change;
}



void LSDParticleColumn::calculate_CRN_conc_from_erosion_timeseries(vector<float> sampling_depths, float d_change, 
                                  double latitude, double longitude, double elevation,
                                  string Muon_scaling, string path_to_atmospheric_data)
{

  // Deal with some cosmogenic parameters first
  LSDCRNParameters LSDCRNP;  
  // set up the scaling
  if (Muon_scaling == "Schaller" )
  {
    LSDCRNP.set_Schaller_parameters();
  }
  else if (Muon_scaling == "Braucher" )
  {
    LSDCRNP.set_Braucher_parameters();
  }
  else if (Muon_scaling == "BraucherBorchers" )
  {
    LSDCRNP.set_BraucherBorchers_parameters();
  }
  else if (Muon_scaling == "Granger" )
  {
    LSDCRNP.set_Granger_parameters();
  }
  else if (Muon_scaling == "newCRONUS" )
  {
    LSDCRNP.set_newCRONUS_parameters();
  }
  else
  {
    cout << "You didn't set the muon scaling." << endl
        << "Options are Schaller, Braucher, newCRONUS, BraucherBorchers, and Granger." << endl
        << "You chose: " << Muon_scaling << endl
        << "Defaulting to Braucher et al (2009) scaling" << endl;
    LSDCRNP.set_Braucher_parameters();
  }

  // get the atmospheric parameters
  LSDCRNP.load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  LSDCRNP.set_CRONUS_data_maps();
  
  // a function for scaling stone production, defaults to 1
  double Fsp = 1.0;

    // now the pressure
  double this_pressure = LSDCRNP.NCEPatm_2(latitude, longitude,elevation);

  // get the production
  double prod_scaling = LSDCRNP.stone2000sp(latitude,this_pressure, Fsp); 


  // ONLY DO THIS FOR TESTING
  //cout << "Prod scaling of this test is: " << prod_scaling << endl;
  //cout << "For testing I am changing this to 1" << endl;
  //prod_scaling = 1;

  // First we get the initial depths and convert these to effective depths
  vector<double> initial_depths;
  vector<double> initial_effective_depths;
  int n_samples = int(sampling_depths.size());
  for (int s = 0; s< n_samples; s++)
  {
    initial_depths.push_back(sampling_depths[s]+d_change);
    initial_effective_depths.push_back( (sampling_depths[s]+d_change)*RockDensity*0.1 );
    //cout << "Insert depth is: " << sampling_depths[s]+d_change << " and effective depth: " << (sampling_depths[s]+d_change)*RockDensity*0.1 << endl;
  }

  // now insert these particles. 
  // for the time being we have default variables for x,y,z but this will need to change if we are to get
  // actual cosmo concentrations.
  double startx = 0;
  double starty = 0;
  double startz = elevation;
  int starttype = 0;

  vector<bool> nuclide_scaling_switches(4,false);
  nuclide_scaling_switches[0] = true;   // 10Be
  nuclide_scaling_switches[1] = true;   // 26Al
  nuclide_scaling_switches[3] = true;   // 14C

  // reset the scaling for this particular location
  // Then calculate the F values. 
  reset_scaling(LSDCRNP, Muon_scaling);
  LSDCRNP.scale_F_values(prod_scaling,nuclide_scaling_switches);


  int n_time_steps = int(time_vec.size());
  double ss_rate = rate_vec[n_time_steps-1];
  double eff_ss_rate = ss_rate*RockDensity*0.1;

  float dt;
  double removal_depth;
  double effective_erate;

  // now loop through the particles and calculate the concentrations
  // We are going to do the concentrations of 10Be, 26Al, and 14C
  for (int s = 0; s< n_samples; s++)
  {
    
    // Initiate the particle. 
    // It has a start depth and initial starting depth that is the sampling depth plus the change
    // induced by the erosion history. 
    LSDCRNParticle CRN_tp(starttype,startx,starty,initial_depths[s],initial_effective_depths[s],startz);


    CRN_tp.update_10Be_SSfull(eff_ss_rate, LSDCRNP);
    CRN_tp.update_26Al_SSfull(eff_ss_rate, LSDCRNP);
    CRN_tp.update_14C_SSfull(eff_ss_rate, LSDCRNP);

    //cout << "Getting initial concentration, eff depth is: " <<  initial_effective_depths[s] << " initial eff rate is: " 
    //     <<   eff_ss_rate << " initial conc is: " << CRN_tp.getConc_10Be() << endl;

    // now we loop through the erosion record.
    // We do the penultimate time last, and then work our
    // way forward in time
    // you need a different logic for the last step
    for (int ts = n_time_steps-2; ts >= 1; ts --)
    {

      // check if this is a removal or erosion step
      if (rate_vec[ts]==0)
      {
        // this is a removal step
        removal_depth = removal_vec[ts];

        CRN_tp.remove_mass_one_layer(removal_depth , RockDensity);
        // we don't need to update time or anything since the removal depth is instantaneous
      }
      else
      {
        
        
        // get the time interval
        dt = time_vec[ts]-time_vec[ts-1];

        // get the effective erosion rate
        effective_erate = rate_vec[ts]*RockDensity*0.1;

        CRN_tp.update_10Be_conc(dt,effective_erate, LSDCRNP);
        CRN_tp.update_26Al_conc(dt,effective_erate, LSDCRNP);
        CRN_tp.update_14C_conc(dt,effective_erate, LSDCRNP);

        // now update the depths by removing the mass
        CRN_tp.remove_mass_one_layer(rate_vec[ts]*dt , RockDensity);

        //cout << "This is a change in erosion rate step" << endl;
        //cout << "dt: " << dt << " and eff erosion rate is: " << effective_erate << endl; 
        //cout << "removed " << rate_vec[ts]*dt << " m of material" << endl;       
      }     
    }
    // now we need to do the last step
    dt = time_vec[0];
    // get the effective erosion rate
    effective_erate = rate_vec[0]*RockDensity*0.1;

    CRN_tp.update_10Be_conc(dt,effective_erate, LSDCRNP);
    CRN_tp.update_26Al_conc(dt,effective_erate, LSDCRNP);
    CRN_tp.update_14C_conc(dt,effective_erate, LSDCRNP);
    CRN_tp.remove_mass_one_layer(rate_vec[0]*dt , RockDensity);

    //cout << "This is the last step" << endl;
    //cout << "dt: " << dt << " and eff erosion rate is: " << effective_erate << endl; 
    //cout << "removed " << rate_vec[0]*dt << " m of material" << endl;    

    CRNParticleList.push_back(CRN_tp);
  }   // end particle, return loop to next particle.


  cout << "I've finished your transient history." << endl;  

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// update_list
// this is the function that inserts a particle into a list
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::insert_particle_into_column(int start_type, 
   double startxLoc, double startyLoc,
	 double start_depth, double zeta)
{
	double d = start_depth;
	double eff_d;                   // effective depth in g/cm^2
	double z_p = zeta-start_depth;  // starting elevation of the particle
	
	// first check to see if we use depth column
	if (UseDenstyProfile)
	{
    cout << "LSDParticleColumn::insert_part_into_column, you asked to use density\n";
    cout << "profile but this is not implemented yet. Fatal error";
    exit(EXIT_FAILURE);
  }
  else
  {
    if (d > SoilThickness)
    {
       // get effective depth if the particle is below the soil layer
       eff_d = SoilThickness*SoilDensity*0.1+ (d-SoilThickness)*RockDensity*0.1;
    }
    else
    {
      // gets the effective depth if the particle is in the soil
      eff_d =  SoilThickness*SoilDensity*0.1;
    }
  }
		
	LSDCRNParticle CRN_tp(start_type,startxLoc,startyLoc,d,eff_d,z_p);
	CRNParticleList.push_back(CRN_tp);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function updates the cosmo concentrations 
// It also calculates 'effective' concentration for particles
// leaving the surface so they only gain nuclide concetration
// while still in the regolith. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDParticleColumn LSDParticleColumn::update_CRN_list_rock_only_eros_limit_3CRN(
	double dt, double uplift_rate,
	int start_type, double start_depth,
	double startxLoc, double startyLoc,
	double zeta_old,double zeta_new,
  double particle_spacing, LSDCRNParameters& CRN_param)
{
  double d;                // depth of particle (m)
  double eff_d;            // effective depth in g/cm^2
  double z_p;              // elevation of particle (m)
  double eff_eros_rate;    // effective erosion rate in g/cm^2/yr
  int eroded_test;         // used to see if particle has eroded
  double effective_dt;     // the time of exposure for particles that have eroded
  double d_frac;           // fraction of depth eroded particles have spent in soil

  // the list and iterators for eroded particles
  list<LSDCRNParticle> eroded_list;
  list<LSDCRNParticle>::iterator part_iter;
  list<LSDCRNParticle>::iterator remove_iter;

  // so first, determine the depth of material lost
  double depth_lost = uplift_rate*dt- (zeta_new-zeta_old);
  //cout << "YoYoMa, line 267, uplift rate is: " << uplift_rate << endl;
  
  // the elevation of the uplifted old surface, for determining d_frac
  double uplift_surf = zeta_old+uplift_rate*dt;
  
  // now get the effective erosion rate
  eff_eros_rate = RockDensity*0.1*depth_lost/dt;

  // print the effective erosion rate for debugging
  //cout << "Heyjabbajsbba, line 275, effective erosion rate in LSDParticleColumn is: "
  //     << eff_eros_rate << " and eros is: " << depth_lost/dt << endl;

  // go through and update the CRN concentrations
  // in the list
  part_iter = CRNParticleList.begin();
  while(part_iter != CRNParticleList.end())
  {
    // get the old zeta location
    z_p = ( *part_iter ).get_zetaLoc();
    
    // now add the uplift to that location
    z_p =  z_p +  uplift_rate*dt; 
    
    d = zeta_new-z_p;
    eff_d = d*RockDensity*0.1;

    // check to see if the particle is above the surface, if
    // it is use an abbreviated exposure
    if (z_p >= zeta_new)
    {
      d_frac = (uplift_surf-z_p)/(depth_lost);
      effective_dt = d_frac*dt;

      // it has zero depth (sampled form surface)
      d = 0;
      eff_d = 0;
    }
    else
    {
      effective_dt = dt;
    }

    // update the CRN_concentrations
    ( *part_iter ).update_10Be_conc(effective_dt,eff_eros_rate, CRN_param);
    ( *part_iter ).update_14C_conc(effective_dt,eff_eros_rate, CRN_param);
    ( *part_iter ).update_21Ne_conc(effective_dt,eff_eros_rate, CRN_param);

    // update the depths
    ( *part_iter ).update_depths(d, eff_d);

    // update the zeta locations (these have been advected by uplift
    ( *part_iter ).update_zetaLoc(z_p);

    part_iter++;
  }

  // now go through the list and see if the particles
  // are either eroded or a new particle needs to be added
  // particles are added to the back of the list and eroded from the front
  // so first check for erosion
  eroded_test = 0;
  do
  {
    part_iter = CRNParticleList.begin();
    z_p = ( *part_iter ).get_zetaLoc();

    // if the elevation of the particle is greater than the elevation
    // of the surface, it is eroded from the particle list
    // and eroded_test remains unchanged
    // if the elevation of the particle is not greater than the surface
    // elevation, then eroded_test goes to 1 and the loop is exited
    if (z_p > zeta_new)
    {
      //cout << "LINE 82 popping out!" << endl;
      eroded_list.push_back(*part_iter);
      CRNParticleList.pop_front();
    }
    else
    {
      eroded_test = 1;
    }
  } while (eroded_test == 0);

  // now see if we insert a particle
  d = (CRNParticleList.back()).getdLoc();
  //cout << "bottom particle depth is: " << d << endl;
  if ( (start_depth - d) >= particle_spacing)
  {
    //cout << "LINE 95 inserting_particle!" << endl;
    insert_particle_into_column( start_type,startxLoc, startyLoc,
                                start_depth, zeta_new);
  }

  // create a column in the same place but with the eroded nodes only
  LSDParticleColumn eroded_column(Row, Col, NodeIndex, SoilDensity, RockDensity,   
        SoilThickness, DataResolution, UseDenstyProfile, 
        DensityDepths, DensityDensities, eroded_list);

  return eroded_column;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function collects particles so that they can be used by 
// to estimate erosion rates or used for other aggregating purposes
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDParticleColumn::collect_particles(vector<LSDParticleColumn>& ColList_vec)
{
  // get the number of columns 
  int N_cols = ColList_vec.size();
  
  // loop through them, aggregating the particles. 
  list<LSDCRNParticle>::iterator part_iter;

  for (int i = 0; i<N_cols; i++)
  {
    list<LSDCRNParticle> list_from_col = ColList_vec[i].getCRNParticleList();
    
    // loop through the columns, collecting particles
    part_iter =  list_from_col.begin();
    while(part_iter != list_from_col.end())
    {  
      CRNParticleList.push_back(*part_iter);
      part_iter++;
    }
  }  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function calculates the apparent erosion rate from 3CRNs based on neutron
// only assumption from the first particle in the list
// without mixing, this will be the top particle
// rho_r is the rock density
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDParticleColumn::calculate_app_erosion_3CRN_neutron_rock_only(
                                                   LSDCRNParameters& CRN_param )
{
  vector<double> apparent_erosion(3,0.0);
  list<LSDCRNParticle>::iterator part_iter;

  // first go through and update the CRN concentrations
  // in the list
  part_iter = CRNParticleList.begin();
  apparent_erosion[0] = (*part_iter).apparent_erosion_10Be_neutron_only(RockDensity, CRN_param);
  apparent_erosion[1] =(*part_iter).apparent_erosion_14C_neutron_only(RockDensity, CRN_param);
  apparent_erosion[2] =(*part_iter).apparent_erosion_21Ne(RockDensity, CRN_param);
  
  //cout << "app erosion from particle at "<< (*part_iter).getdLoc() << " depth" << endl;
    
  return apparent_erosion; 

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for bug checking. It prints out the particles, and their
// properties to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParticleColumn::print_particle_properties_3CRN(string csv_outname)
{

  ofstream outfile;

  // make sure the file has the csv extension
  string ext_str = ".csv";
  if (csv_outname.find(ext_str) == std::string::npos) 
  {
    csv_outname = csv_outname + ext_str;
  }
  outfile.open(csv_outname.c_str());

  outfile << "depth,effective_depth, C_10Be, C_26Al, C_14C" << endl;

  list<LSDCRNParticle>::iterator part_iter;
  part_iter = CRNParticleList.begin();
  while(part_iter != CRNParticleList.end())
  {
    // get the zeta location
    double d_loc =   ( *part_iter ).getdLoc();
    double effD =  ( *part_iter ).geteffective_dLoc();
    double conc10Be = ( *part_iter ).getConc_10Be();
    double conc26Al = ( *part_iter ).getConc_26Al();    
    double conc14C = ( *part_iter ).getConc_14C();

    outfile << d_loc << "," << effD << "," << conc10Be << "," << conc26Al << "," << conc14C << endl;

    part_iter++;
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for bug checking. It prints out the particles, and their
// properties to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDParticleColumn::print_particle_properties_to_screen(LSDCRNParameters& CRNparam)
{

  cout << "\n\nColumn at ["<<Row<<"]["<<Col<<"]; density is: " << RockDensity << endl;
  list<LSDCRNParticle>::iterator part_iter;
  part_iter = CRNParticleList.begin();
  while(part_iter != CRNParticleList.end())
  {
    // get the zeta location
    double z_p = ( *part_iter ).get_zetaLoc();
    double d_loc =   ( *part_iter ).getdLoc();
    double effD =  ( *part_iter ).geteffective_dLoc();
    double conc10Be = ( *part_iter ).getConc_10Be();
          //double conc14C = ( *part_iter ).getConc_14C();
          //double conc21Ne = ( *part_iter ).getConc_21Ne();
    double appEros10Be = ( *part_iter ).apparent_erosion_10Be_neutron_only(RockDensity, CRNparam);
          //double appEros14C = ( *part_iter ).apparent_erosion_14C_neutron_only(RockDensity, CRNparam);
          //double appEros21Ne = ( *part_iter ).apparent_erosion_21Ne_neutron_only(RockDensity, CRNparam);


    cout << z_p << "\t " << d_loc << "\t" << effD << "\t" << conc10Be << "\t" << appEros10Be << endl;


    part_iter++;
  }
}

#endif
