//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDCRNParameters.cpp
//
// Land Surface Dynamics Cosmogenic Radionuclide Parameters Object
//
// This keeps track of paramters used to calculate the evolution of 
// in situ cosmogenic nuclides. 
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
#include <map>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
#include "LSDCRNParameters.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDCRNParameters_CPP
#define LSDCRNParameters_CPP


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// the LSDCRNParameters object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// function to set CRN parameters
void LSDCRNParameters::create()
{
  version = "1.0";
  
  S_t = 1;
  neutron_S_t = 1;

  // from Vermeesh 2007
  lambda_10Be = 456e-9;		// in yr-1
  lambda_26Al = 980e-9;		// in yr-1
  lambda_14C = 121e-6;		// in yr-1
  lambda_36Cl = 230e-8;		// in yr-1

  // from Vermeesh 2007
  P0_10Be = 5.11;					// in a/g/yr
  P0_26Al = 30.31;				// in a/g/yr
  P0_14C = 5.86;					// in a/g/yr
  P0_36Cl = 55.45;				// in a/g/yr
  P0_21Ne = 20.29;				// in a/g/yr
  P0_3He = 97.40;					// in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 0.9724;
  F_10Be[1] = 0.0186;
  F_10Be[2] = 0.004;
  F_10Be[3] = 0.005;

  // dimensionless
  F_26Al[0] = 0.9655;
  F_26Al[1] = 0.0233;
  F_26Al[2] = 0.005;
  F_26Al[3] = 0.0062;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.0691;
  F_14C[2] = 0.0809;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.903;
  F_36Cl[1] = 0.0447;
  F_36Cl[2] = 0.05023;
  F_36Cl[3] = 0.0;
}

// this function gets the parameters used to convert elevation to 
// pressure
void LSDCRNParameters::load_parameters_for_atmospheric_scaling(string path_to_data)
{
  cout.precision(8);
  
  // first load the levels
  levels.push_back(1000);
  levels.push_back(925);
  levels.push_back(850);
  levels.push_back(700);
  levels.push_back(600);
  levels.push_back(500);
  levels.push_back(400);
  levels.push_back(300);
  
  // the dimensions of the data
  int n_levels = 8;
  int NRows = 73;
  int NCols = 145;
  Array2D<double> new_slp(NRows,NCols,0.0);
  Array2D<double> new_meant(NRows,NCols,0.0);
    
  // now load the mean sea level pressure
  string filename = "NCEP2.bin";
  filename = path_to_data+filename;
  //cout << "Loading mean sea level, file is: " << endl << filename << endl;

  ifstream ifs_data(filename.c_str(), ios::in | ios::binary);
  if( ifs_data.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  
  

  double temp;
  //cout << "The size of a double is: " << sizeof(temp) << endl;
  for (int i=0; i<NCols; ++i)
  {
    for (int j=0; j<NRows; ++j)
    {
      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      new_slp[j][i] = temp;
      //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    }
  }
  
  for (int i=0; i<NCols; ++i)
  {
    for (int j=0; j<NRows; ++j)
    {
      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      new_meant[j][i] = temp;
      //cout << "new_meant100["<<j+1<<"]["<<i+1<<"]: " << new_meant[j][i] << endl;
    }
  }  
  
  // now get the indices
  vector<double> temp_lat(NRows,0.0);
  for (int i=0; i<NRows; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    temp_lat[i] = temp;
    //cout << "Lat["<<i+1<<"]: " << temp_lat[i] << endl;
  }
  vector<double> temp_long(NCols,0.0);
  for (int i=0; i<NCols; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    temp_long[i] = temp;
    //cout << "Long["<<i+1<<"]: " << temp_long[i] << endl;
  }  
  
  ifs_data.close();
  
  
  // now the data with levels
  filename = "NCEP_hgt.bin";
  filename = path_to_data+filename;
  //cout << "Loading hgt, file is: " << endl << filename << endl;

  ifstream ifs_data2(filename.c_str(), ios::in | ios::binary);
  if( ifs_data2.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  
  

  // get the gm heights
  vector< Array2D<double> > vec_hgt_gm_array;
  for (int lvl = 0; lvl < n_levels; lvl++)
  {
    Array2D<double> current_hgt_array(NRows,NCols,0.0);
    for (int i=0; i<NCols; ++i)
    {
      for (int j=0; j<NRows; ++j)
      {
        ifs_data2.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        current_hgt_array[j][i] = temp;
        //cout << "new_slp["<<lvl<<"]["<<j+1<<"]["<<i+1<<"]: " << current_hgt_array[j][i] << endl;
      }
    }
    //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    vec_hgt_gm_array.push_back(current_hgt_array.copy());
  }

  // now the gp heights
  vector< Array2D<double> > vec_hgt_gp_array;
  for (int lvl = 0; lvl < n_levels; lvl++)
  {
    Array2D<double> current_hgt_array(NRows,NCols,0.0);
    for (int i=0; i<NCols; ++i)
    {
      for (int j=0; j<NRows; ++j)            
      {
        ifs_data2.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        current_hgt_array[j][i] = temp;
        //cout << "new_slp["<<lvl<<"]["<<j+1<<"]["<<i+1<<"]: " << current_hgt_array[j][i] << endl;
      }
    }
    //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    vec_hgt_gp_array.push_back(current_hgt_array.copy());
  }
  ifs_data2.close();
  
  // now update the data elements
  NCEPlat = temp_lat;
  NCEPlon = temp_long;
  
  meanslp = new_slp.copy();
  meant1000 = new_meant.copy();
  gp_hgt = vec_hgt_gp_array;
  gm_hgt = vec_hgt_gp_array;  
  
  //cout << "Size lat: " << NCEPlat.size() << " size long: " << NCEPlon.size() << endl;
  //cout << "size slp:" << meanslp.dim1() << " " << meanslp.dim2() << endl;
  //cout << "size t1000: " << meant1000.dim1() << " " << meant1000.dim2() << endl;
  

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function sets a number of paramters that are used for acalucaltion
// of scaling and error propigation
//
// They are constants used in the CRONUS caluclator, and have been ported
// from make_al_be_consts_v22.m written by Greg Balco
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_CRONUS_data_maps()
{
  //cout << "Line 278, creating the CRONUS data maps" << endl;
  
  map<string,double> temp_map;
  
  // 10Be decay: this is specific to the CRONUS calculator
  temp_map["l10"] = -log(0.5)/1.387e6; // Chmeleff/Korschinek value
  // Al-26 decay constant -- value compatible with Nishiizumi standards
  temp_map["l26"] = 9.83e-7;
  
  
  //cout << "Creating CRONUS data maps, lambda 10: " << temp_map["l10"]
  //     << " and lambda 26: " << temp_map["l26"] << endl;

  // Note that the uncertainty is not used in exposure-age or erosion-rate
  // calculators. Here only for development purposes
  double dldt = -log(0.5)*(1/(1.387e6*1.387e6));
  temp_map["dell10"] = sqrt((dldt*1/(0.012e6*0.012e6)));
  temp_map["dell26"] = 2.5e-8;

  //Effective attenuation length for spallation in rock
  // Commonly accepted value: 160 g/cm2
  // For discussion see Gosse and Phillips (2000)
  temp_map["Lsp"] = 160.0;

  // Fsp - fraction of total production by spallation rather than muons
  // For use with Lal/Stone scaling scheme in exposure age calculation
  // For details, see Stone (2000)
  // This aspect of Stone(2000) is de-emphasized in version 2. These constants
  // are only in use for historical comparisons and quick initial guesses for 
  // the exposure age and erosion rate solvers. 

  // Note that they are probably incorrect WRT Be-10 restandardization. Don't
  // use these for anything important. 
  temp_map["Fsp10"] = 0.978;
  temp_map["Fsp26"] = 0.974;

  // Be-10 standardization info.
  // Standards comparison/conversion lookup table
  temp_map["Be_std_07KNSTD"] = 1.0000;
  temp_map["Be_std_KNSTD"] = 0.9042;
  temp_map["Be_std_NIST_Certified"] = 1.0425;
  temp_map["Be_std_LLNL31000"] =  0.8761;
  temp_map["Be_std_LLNL10000"] = 0.9042;
  temp_map["Be_std_LLNL3000"] = 0.8644;
  temp_map["Be_std_LLNL1000"] = 0.9313;
  temp_map["Be_std_LLNL300"] = 0.8562;
  temp_map["Be_std_NIST_30000"] =  0.9313;
  temp_map["Be_std_NIST_30200"] = 0.9251;
  temp_map["Be_std_NIST_30300"] = 0.9221;
  temp_map["Be_std_NIST_30600"] = 0.9130;
  temp_map["Be_std_NIST_27900"] = 1;
  temp_map["Be_std_S555"] = 0.9124;
  temp_map["Be_std_S2007"] = 0.9124;
  temp_map["Be_std_BEST433"] = 0.9124;
  temp_map["Be_std_BEST433N"] = 1;
  temp_map["Be_std_S555N"] = 1;
  temp_map["Be_std_S2007N"] = 1;;

  // Same for Al-26. A zero placeholder is
  // also allowed. 
  temp_map["Al_std_KNSTD"] = 1.0;
  temp_map["Al_std_ZAL94"] = 0.9134;
  temp_map["Al_std_SMAL11"] = 1.021;
  temp_map["Al_std_0"] = 1.0;
  temp_map["Al_std_ZAL94N"] = 1.0;
  temp_map["Al_std_ASTER"] = 1.021;
  temp_map["Al_std_Z92-0222"] = 1;

  // Reference production rates at SLHL for spallation according to various
  // scaling schemes. Letter codes De, Du, Li, and St refer to
  // scaling schemes of Desilets et al. (2006), Dunai (2001), Lifton (2006),
  // and Stone (2000) respectively. The final letter code Lm refers to the
  // paleomagnetically corrected version of the Lal 1991/Stone 2000 scaling
  // factors.

  // These reference production rates are derived using the current version of
  // get_al_be_age.m. 

  // Al-26 and Be-10 production rates are linked by the production ratio,
  // not determined independently. This is because there are not very
  // many geological-calibration sites for Al-26. See the calibration data
  // sets and the accompanying paper for more information. 

  // This reflects the v 2.1 air pressure upgrade.
  // Also the update to 07KNSTD in version 2.2. 

  // Be-10 production rates. Brent's values. Greg agrees. 
  temp_map["P10_ref_St"] = 4.49;
  temp_map["delP10_ref_St"] = 0.39;

  temp_map["P10_ref_Lm"] = 4.39; 
  temp_map["delP10_ref_Lm"] = 0.37;

  temp_map["P10_ref_De"] = 4.41;
  temp_map["delP10_ref_De"] = 0.52;

  temp_map["P10_ref_Du"] = 4.43;
  temp_map["delP10_ref_Du"] = 0.52;

  temp_map["P10_ref_Li"] = 4.87;
  temp_map["delP10_ref_Li"] = 0.48;

  // Al-26 production rates are derived from Be-10 production rates 
  double R2610 = 6.1*1.106; // Update assumed production ratio
  temp_map["P26_ref_St"] = temp_map["P10_ref_St"]*R2610;
  temp_map["delP26_ref_St"] = temp_map["delP10_ref_St"]*R2610;
  temp_map["P26_ref_Lm"] = temp_map["P10_ref_Lm"]*R2610; 
  temp_map["delP26_ref_Lm"] = temp_map["delP10_ref_Lm"]*R2610;
  temp_map["P26_ref_De"] = temp_map["P10_ref_De"]*R2610;
  temp_map["delP26_ref_De"] = temp_map["delP10_ref_De"]*R2610;
  temp_map["P26_ref_Du"] = temp_map["P10_ref_Du"]*R2610;
  temp_map["delP26_ref_Du"] = temp_map["delP10_ref_Du"]*R2610;
  temp_map["P26_ref_Li"] = temp_map["P10_ref_Li"]*R2610;
  temp_map["delP26_ref_Li"] = temp_map["delP10_ref_Li"]*R2610;

  // Muon interaction cross-sections. All follow Heisinger (2002a,b).
  // Note that the energy-dependence-of-muon-interaction-cross-section
  // exponent alpha is treated as model-dependent -- it's internal to 
  // P_mu_total.m and can't be passed.  

  temp_map["Natoms10"] = 2.006e22;
  temp_map["Natoms26"] = 1.003e22;

  // Be-10 interaction cross-sections
  // Restandardized by ETH-07KNSTD factor for version 2.2.1. 
  temp_map["k_neg10"] = (0.704 * 0.1828 * 0.0043)/1.096;
  temp_map["delk_neg10"] = (0.704 * 0.1828 * 0.0003)/1.096;
  temp_map["sigma190_10"] = (0.094e-27)/1.096;
  temp_map["delsigma190_10"] = (0.013e-27)/1.096;
  temp_map["Be10_Natoms_times_sigma190"] =  1.7205e-6;
  temp_map["Be10_Natoms_times_delsigma190"] =  2.37938e-7;
  temp_map["Be10_Natoms_sigma_modified_for_fast"] = 3.361886e-8;
  //cout << "LINE 405, be10 prod: " <<  temp_map["Be10_Natoms_sigma_modified_for_fast"] << endl;

  // Al-26 interaction cross-sections
  temp_map["k_neg26"] = 0.296 * 0.6559 * 0.022;
  temp_map["delk_neg26"] = 0.296 * 0.6559 * 0.002;
  temp_map["sigma190_26"] = 1.41e-27;
  temp_map["delsigma190_26"] = 0.17e-27;
  temp_map["Al26_Natoms_times_sigma190"] =  1.41423e-5;
  temp_map["Al26_Natoms_times_delsigma190"] =  1.7051e-6;
  temp_map["Al26_Natoms_sigma_modified_for_fast"] =  2.76347e-8;
  

  // Paleomagnetic records for use in time-dependent production rate schemes
  // Derived from Nat Lifton's compilation of paleomagnetic data from
  // various sources. See Lifton et al. (2006) and Pigati and Lifton (2005).

  // Load the magnetic field data
  // load PMag_Mar07
  // Relative dipole moment and time vector
  // al_be_consts.M = MM0; 
  // al_be_consts.t_M = t_M; 
  // These start at 7500 yr -- time slices are 7500,8500,9500,10500,11500
  // in order to use data from Yang et al; subsequent time slices are 
  // 12000:1000:800000 for SINT800 data; final two time points are 801000
  // and Inf. 

  // Cutoff rigidity blocks for past 6900 yr. 
  // TTRc and IHRC are lon x lat x time blocks of Rc values for the past 
  // 6900 years.
  // Both are derived by Nat Lifton from the magnetic field reconstructions of
  // Korte and Constable. 
  // TTRC has cutoff rigidity obtained by trajectory tracing -- these are for
  // the Lifton and Desilets scaling factors. IHRc has cutoff rigidity
  // obtained by finding magnetic inclination and horizontal field strength
  // from the field model, then applying Equation 2 of Dunai(2001). 
  // al_be_consts.TTRc = TTRc; % data block
  // al_be_consts.IHRc = IHRc; % data block
  // al_be_consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
  // al_be_consts.lon_Rc = lon_Rc;
  // al_be_consts.t_Rc = t_Rc; % time vector for Rc data block

  // Effective pole positions and field strengths inferred from K and C field
  // reconstructions for last 7000 yr. These are used in the
  // paleomagnetically-corrected implementation of the Lal SF. They are for
  // the same times as the RC slices in the data block above. Again,
  // generated by Nat Lifton -- hence KCL = Korte-Constable-Lifton. 

  // al_be_consts.MM0_KCL = MM0_KCL;
  // al_be_consts.lat_pp_KCL = lat_pp_KCL;
  // al_be_consts.lon_pp_KCL = lon_pp_KCL;

  // Solar variability from Lifton et al. 2005
  // Has been averaged and resampled to the same time slices as everything
  // else. 

  // al_be_consts.S = S; 
  // al_be_consts.SInf = 0.95; % Long-term mean S value;

  // set the data member map
  CRONUS_data_map = temp_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function retrieves the Stone production reference values for 10Be and 
// 26Al
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParameters::get_Stone_Pref()
{
  // first check to see if CRONUS data maps are set
  if(CRONUS_data_map.find("l10") == CRONUS_data_map.end())
  {
    cout << "You haven't set the CRONUS data map. I'm doing that for you now!" << endl;
    set_CRONUS_data_maps();
  }
  
  vector<double> Stone_pref(2,0.0);
  Stone_pref[0] = CRONUS_data_map["P10_ref_St"];
  Stone_pref[1] = CRONUS_data_map["P26_ref_St"];  
  
  return Stone_pref;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function wraps the CRONUS muon production function
// It returns a vector with elements
// Muon_production[0] = 10Be fast
// Muon_production[1] = 26Al fast
// Muon_production[2] = 10Be neg
// Muon_production[3] = 26Al neg
// 14/10/2014
// SMM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParameters::calculate_muon_production_CRONUS(double z, double h)
{
  vector<double> Muon_production(4,0.0);
  P_mu_total(z,h);
  
  
  Muon_production[0] = CRONUS_muon_data["P_fast_10Be"];
  Muon_production[1] = CRONUS_muon_data["P_fast_26Al"];
  Muon_production[2] = CRONUS_muon_data["P_neg_10Be"];
  Muon_production[3] = CRONUS_muon_data["P_neg_26Al"];

  return Muon_production;
  
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  
// This function sets CRONUS muon production paramteters
//
// Calculates the production rate of Al-26 or Be-10 by muons
// as a function of depth below the surface z (g/cm2) and
// site atmospheric pressure h (hPa).
//
// out.phi_vert_slhl muons/cm2/s/sr
// out.R_vert_slhl muons/g/s/sr
// out.R_vert_site muons/g/s/sr
// out.phi_vert_site muons/cm2/s/sr
// out.phi muons/cm2/yr
// out.R muons/g/yr
// out.P_fast atoms/g/yr
// out.P_neg atoms/g/yr
// out.Beta nondimensional
// out.Ebar GeV
// out.H g/cm2
// out.LZ g/cm2
//
// This uses the scheme in Heisinger and others (2002, 2 papers). The
// vertically traveling muon flux is scaled to the site elevation using
// energy-dependent attenuation lengths from Boezio et al. (2000). See the 
// hard-copy documentation for detailed citations and a full discussion of
// the calculation. 
//
// Note that some constants are internal to the function. The only ones that
// get passed from upstream are the ones that a) are nuclide-specific, or b) 
// actually have quoted uncertainties in Heisinger's papers. 
// The fraction of muons that are negative is internal; so is the
// energy-dependence exponent alpha.
//
// Original Written by Greg Balco -- UW Cosmogenic Nuclide Lab
// balcs@u.washington.edu
// March, 2006
// Part of the CRONUS-Earth online calculators: 
//      http://hess.ess.washington.edu/math
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::P_mu_total(double z,double h)
{
  map<string,double> temp_data_map;

  //cout << "CHECKING MUON FLUX, Line 504 " << endl;

  // first check to see if CRONUS data maps are set
  if(CRONUS_data_map.find("l10") == CRONUS_data_map.end())
  {
    cout << "You haven't set the CRONUS data map. I'm doing that for you now!" << endl;
    set_CRONUS_data_maps();
  }


  // calculator the atmospheric depth in g/cm2
  double H = (1013.25 - h)*1.019716;
  //cout << "Atmospheric depth is: " << H << " g/cm^2" << endl;
  

  // find the vertical flux at SLHL
  double a = 258.5*(pow(100,2.66));
  double b = 75*(pow(100,1.66));
  
  // this is equation (1) From Heisinger 2002a and (3) from Heisinger 2002b
  double phi_vert_slhl = (a/((z+21000.0)*((pow((z+1000),1.66)) + b)))
                            *exp(-5.5e-6* z);

  // The above expression is only good to 2e5 g/cm2. We don't ever consider production
  // below that depth. The full-depth scheme appears in the comments below.
  // ------ begin full-depth flux equations -------
  //phiz_1 = (a./((z+21000).*(((z+1000).^1.66) + b))).*exp(-5.5e-6 .* z);
  //phiz_2 = 1.82e-6.*((121100./z).^2).*exp(-z./121100) + 2.84e-13;
  //out(find(z<200000)) = phiz_1(find(z<200000));
  //out(find(z>=200000)) = phiz_2(find(z>=200000));
  // ------ end full-depth flux equations -------

  // find the stopping rate of vertical muons at SLHL
  // this is done in a subfunction Rv0, because it gets integrated later.
  double R_vert_slhl = Rv0(z);

  // find the stopping rate of vertical muons at site
  double R_vert_site = R_vert_slhl*exp(H/LZ(z));
  //cout << "LZ(" << z << "): " << LZ(z) << endl;
  //cout << "R_vert_slhl: " << R_vert_slhl << " R_vert_site: " << R_vert_site << endl;

  // find the flux of vertical muons at site
  // integrate
  // ends at 200,001 g/cm2 to avoid being asked for an zero
  // range of integration -- 
  // get integration tolerance -- want relative tolerance around
  // 1 part in 10^4. 
  double tol = phi_vert_slhl*1e-4;
  double phi_vert_site = integrate_muon_flux(z, H, tol);
  
  //=====================
  // I THINK below here is an error in Balco's code
  // The below equation uses a which is calculated above, 
  // but in balco's code a is then used as an index, so the below
  // equation takes the index value rather than the precalculated 
  // value of a
  //=======================
  
  // invariant flux at 2e5 g/cm2 depth - constant of integration
  // calculated using commented-out formula above
  double phi_200k = (a/((2.0e5+21000.0)*((pow((2.0e5+1000.0),1.66)) + b)))
                      *exp(-5.5e-6 * 2.0e5);
  //double test_balco_error = (1.0/((2.0e5+21000.0)*((pow((2.0e5+1000.0),1.66)) + b)))
  //                    *exp(-5.5e-6 * 2.0e5);
                      
  phi_vert_site = phi_vert_site + phi_200k;
  
  // find the total flux of muons at site
  // angular distribution exponent
  double nofz = 3.21 - 0.297*log((z+H)/100.0 + 42.0) + 1.21e-5*(z+H);
  // derivative of same
  double dndz = (-0.297/100.0)/((z+H)/100.0 + 42.0) + 1.21e-5;

  // caluculate phi in muons/cm2/s
  double phi_temp = (phi_vert_site*2* M_PI) / (nofz+1.0);

  // convert to muons/cm2/yr
  double phi = phi_temp*60.0*60.0*24.0*365.0;

  // find the total stopping rate of muons at site in muons/g/s
  double R_temp = (2*M_PI/(nofz+1.0))*R_vert_site 
                  - phi_vert_site*(-2*M_PI*(1/((nofz+1.0)*(nofz+1.0))))*dndz;
    
  // convert to negative muons/g/yr
  double this_R = R_temp*0.44*60.0*60.0*24.0*365.0;

  // Now calculate the production rates. 
  // Depth-dependent parts of the fast muon reaction cross-section
  double Beta = 0.846 - 0.015 * log((z/100.0)+1.0) 
                      + 0.003139 * (log((z/100.0)+1.0)*log((z/100.0)+1.0));
  double Ebar = 7.6 + 321.7*(1 - exp(-8.059e-6*z)) 
                    + 50.7*(1-exp(-5.05e-7*z));

  // internally defined constants
  double aalpha = 0.75;
  
  double sigma0_Be10 = CRONUS_data_map["sigma190_10"]/(pow(190.0,aalpha));
  double sigma0_Al26 = CRONUS_data_map["sigma190_26"]/(pow(190.0,aalpha));
  
  // fast muon production
  double P_fast_Be10 = phi*Beta*(pow(Ebar,aalpha))
                          *sigma0_Be10*CRONUS_data_map["Natoms10"];
  double P_fast_Al26 = phi*Beta*(pow(Ebar,aalpha))
                          *sigma0_Al26*CRONUS_data_map["Natoms26"];
  
  //cout << "Phi: " << phi << " Beta " << Beta << " Ebar: " << Ebar << endl
  //     << "aalpha: " << aalpha << endl
  //     << " prod10: " << CRONUS_data_map["Be10_Natoms_sigma_modified_for_fast"] << endl
  //     << " prod26: " << CRONUS_data_map["Al26_Natoms_sigma_modified_for_fast"] << endl;
       
         
  //double P2_10Be = phi*Beta*(pow(Ebar,aalpha))*
  //                  CRONUS_data_map["Be10_Natoms_sigma_modified_for_fast"];
  //double P2_26Al = phi*Beta*(pow(Ebar,aalpha))*
  //                  CRONUS_data_map["Al26_Natoms_sigma_modified_for_fast"];                    
  
  //cout << "Pfast10be: " << P_fast_Be10 << " P2: " << P2_10Be << endl;
  //cout << "Pfast26Al: " << P_fast_Al26 << " P2: " << P2_26Al << endl;
  
  // negative muon capture
  double P_neg_Be10 = this_R*CRONUS_data_map["k_neg10"];
  double P_neg_Al26 = this_R*CRONUS_data_map["k_neg26"];

  //cout << "Sig0: " << sigma0_Be10 << " Pfast: " << P_fast_Be10 << " P_neg: " << P_neg_Be10 << endl;

  temp_data_map["phi_vert_slhl"] = phi_vert_slhl;
  temp_data_map["R_vert_slhl"] = R_vert_slhl;
  temp_data_map["phi_vert_site"] = phi_vert_site;
  temp_data_map["R_vert_site"] = R_vert_site;
  temp_data_map["phi"] = phi;
  temp_data_map["R"] = this_R;
  temp_data_map["Beta"] = Beta;
  temp_data_map["Ebar"] = Ebar;
  temp_data_map["P_fast_10Be"] = P_fast_Be10;
  temp_data_map["P_fast_26Al"] = P_fast_Al26;
  temp_data_map["P_neg_10Be"] = P_neg_Be10;
  temp_data_map["P_neg_26Al"] = P_neg_Al26;
  temp_data_map["H"] = H;
  temp_data_map["LZ"] = LZ(z);

  CRONUS_muon_data = temp_data_map;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function wraps the P_mu_tot function and replaces the 10Be and
// 26Al production rates
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::P_mu_total_return_nuclides(double z,double h, 
                                   double& Be10_total_mu, double& Al26_total_mu)
{
  // get the muon roduction
  P_mu_total(z,h);
  
  double p10 = CRONUS_muon_data["P_fast_10Be"]+CRONUS_muon_data["P_neg_10Be"];
  double p26 = CRONUS_muon_data["P_fast_26Al"]+CRONUS_muon_data["P_neg_26Al"];
  
  Be10_total_mu = p10;
  Al26_total_mu = p26;
}                                   


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this subfunction does the integration of muon sotpping to get the total
// muon flux
//
// NOTE: THIS IS RATE LIMITING
// Will probably need to come back and try to speed up. 
// The best way is to retain previously calculated balues, so at new
// node spacings only the intermediate values are recalculated. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::integrate_muon_flux(double z, double H, double tolerance)
{
  // we just use a simple simpsons rule but we keep fining the grid until there 
  // is no difference in the solution
  int n_nodes = 101;
  double start_z = z;
  double end_z = 2.0e5+1.0;
  double spacing = (end_z-z)/double(n_nodes-1);
  double a,b, intermediate;
  double fa,fb,fi,sum, last_sum;
  double integral_difference;
  
  // get the inital guess
  b = start_z;
  fb = Rv0(b)*exp(H/LZ(b));
  sum = 0;
  for(int i = 1; i< n_nodes; i++)
  {
    // locations of spacings
    a = b;
    b = a+spacing; 
    intermediate = (b+a)/2.0;
    
    //cout << "z["<<i<<"]: " << a << " i+1/2: " << intermediate << " i+1: " << b << endl;
    
    // functions evaluated at spacings
    fa = fb;
    fi = Rv0(intermediate)*exp(H/LZ(intermediate));
    fb = Rv0(b)*exp(H/LZ(b));
    
    sum+= ((b-a)/6.0)*(fa+4.0*fi+fb);
  }
  last_sum = sum;
  
  //cout << "LINE 660, integrating muon flux, initial guess: " << last_sum << endl;
  double log_error_ratio = 0.5;
  double node_multiplier;
  // now loop until the error tolerance is reached
  do
  {
    // this implements and adaptive spacing between nodes, determined by the 
    // ratio between the error and the tolerance in an attempt to speed up
    // the integration
    node_multiplier = 1.0+log_error_ratio;
    
    // increase the density of the nodes
    n_nodes = int(double(n_nodes)*node_multiplier);
    //cout << "Looping for simpsons, n_nodes: " << n_nodes << endl;
    spacing = (end_z-z)/double(n_nodes-1);
    
    //cout << "Nodes: " << n_nodes << " and spacing: " << spacing << endl;
    
    b = start_z;
    fb = Rv0(b)*exp(H/LZ(b));
    sum = 0;
    for(int i = 1; i< n_nodes; i++)
    {
      // locations of spacings
      a = b;
      b = a+spacing; 
      intermediate = (b+a)/2.0;
    
      // functions evaluated at spacings
      fa = fb;
      fi = Rv0(intermediate)*exp(H/LZ(intermediate));
      fb = Rv0(b)*exp(H/LZ(b));
    
      sum+= ((b-a)/6.0)*(fa+4.0*fi+fb);
      
      
    }
    // compare this integral with the last one
    integral_difference = fabs(last_sum-sum);
    
    log_error_ratio = log(integral_difference/tolerance);
    //cout << "log error_ratio: " << log_error_ratio << endl;
    
    // reset the last sum
    last_sum = sum;
    
    //cout << "this sum: " << last_sum << endl;
  
  } while(integral_difference>tolerance);
  
  // now return the integral
  return last_sum;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this subfunction returns the stopping rate of vertically traveling muons
// as a function of depth z at sea level and high latitude.
// Modified from Greg Balco's CRONUS calculator
// z is the depth below the surface in g/cm^2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::Rv0(double z)
{
  double a = exp(-5.5e-6*z);
  double b = z + 21000.0;
  double c = pow((z + 1000.0),1.66) + 1.567e5;
  double dadz = -5.5e-6 * exp(-5.5e-6*z);
  double dbdz = 1.0;
  double dcdz = 1.66*(pow((z + 1000),0.66));
  
  double out = -5.401e7*(b*c*dadz-a*(c*dbdz+b*dcdz))/(b*b*c*c);
  return out;

  // full depth calculation appears in comments below
  // testing indicates this isn't really necessary
  //R_1 = -5.401e7 .* (b.*c.*dadz - a.*(c.*dbdz + b.*dcdz))./(b.^2 .* c.^2);
  //f = (121100./z).^2;
  //g = exp(-z./121100);
  //dfdz = (-2.*(121100.^2))./(z.^3);
  //dgdz = -exp(-z./121100)./121100;
  //R_2 = -1.82e-6.*(g.*dfdz + f.*dgdz);
  //out(find(z<200000)) = R_1(find(z<200000));
  //out(find(z>=200000)) = R_2(find(z>=200000));
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// this subfunction returns the effective atmospheric attenuation length for
// muons of range Z
// z is the depth in g/cm^2
//
// Original by Greg Balco as part of the CRONUS calculator
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::LZ(double z)
{

  //define range/momentum relation
  // table for muons in standard rock in Groom and others 2001
  // units are range in g cm-2 (column 2)
  //momentum in MeV/c (column 1)
  vector<double> data_for_LZ_range;
  vector<double> data_for_LZ_momentum;
  data_for_LZ_momentum.push_back(4.704e1);
  data_for_LZ_range.push_back(8.516e-1);
  data_for_LZ_momentum.push_back(5.616e1); 
  data_for_LZ_range.push_back(1.542e0);
  data_for_LZ_momentum.push_back(6.802e1); 
  data_for_LZ_range.push_back(2.866e0);
  data_for_LZ_momentum.push_back(8.509e1);
  data_for_LZ_range.push_back(5.698e0);
  data_for_LZ_momentum.push_back(1.003e2);
  data_for_LZ_range.push_back(9.145e0);
  data_for_LZ_momentum.push_back(1.527e2);
  data_for_LZ_range.push_back(2.676e1);
  data_for_LZ_momentum.push_back(1.764e2); 
  data_for_LZ_range.push_back(3.696e1);
  data_for_LZ_momentum.push_back(2.218e2);
  data_for_LZ_range.push_back(5.879e1);
  data_for_LZ_momentum.push_back(2.868e2);
  data_for_LZ_range.push_back(9.332e1);
  data_for_LZ_momentum.push_back(3.917e2);
  data_for_LZ_range.push_back(1.524e2);
  data_for_LZ_momentum.push_back(0.945e2);
  data_for_LZ_range.push_back(2.115e2);
  data_for_LZ_momentum.push_back(8.995e2);
  data_for_LZ_range.push_back(4.418e2);
  data_for_LZ_momentum.push_back(1.101e3);
  data_for_LZ_range.push_back(5.534e2);
  data_for_LZ_momentum.push_back(1.502e3);
  data_for_LZ_range.push_back(7.712e2);
  data_for_LZ_momentum.push_back(2.103e3);
  data_for_LZ_range.push_back(1.088e3);
  data_for_LZ_momentum.push_back(3.104e3);
  data_for_LZ_range.push_back(1.599e3);
  data_for_LZ_momentum.push_back(4.104e3);
  data_for_LZ_range.push_back(2.095e3);
  data_for_LZ_momentum.push_back(8.105e3);
  data_for_LZ_range.push_back(3.998e3);
  data_for_LZ_momentum.push_back(1.011e4);
  data_for_LZ_range.push_back(4.920e3);
  data_for_LZ_momentum.push_back(1.411e4);
  data_for_LZ_range.push_back(6.724e3);
  data_for_LZ_momentum.push_back(2.011e4);
  data_for_LZ_range.push_back(9.360e3);
  data_for_LZ_momentum.push_back(3.011e4);
  data_for_LZ_range.push_back(1.362e4);
  data_for_LZ_momentum.push_back(4.011e4);
  data_for_LZ_range.push_back(1.776e4);
  data_for_LZ_momentum.push_back(8.011e4);
  data_for_LZ_range.push_back(3.343e4);
  data_for_LZ_momentum.push_back(1.001e5);
  data_for_LZ_range.push_back(4.084e4);
  data_for_LZ_momentum.push_back(1.401e5);
  data_for_LZ_range.push_back(5.495e4);
  data_for_LZ_momentum.push_back(2.001e5);
  data_for_LZ_range.push_back(7.459e4);
  data_for_LZ_momentum.push_back(3.001e5); 
  data_for_LZ_range.push_back(1.040e5);
  data_for_LZ_momentum.push_back(4.001e5);
  data_for_LZ_range.push_back(1.302e5);
  data_for_LZ_momentum.push_back(8.001e5);
  data_for_LZ_range.push_back(2.129e5);

  // deal with zero situation
  if(z < 1)
  {
    z = 1.0;
  }

  double log_z = log(z);
  //cout << "z is:" << z << " and log z is: " << log_z << endl;


  // obtain momenta
  // use log-linear interpolation
  int n_momentum_dpoints = int(data_for_LZ_momentum.size());
  vector<double> log_momentum;
  vector<double> log_range;
  for(int i = 0; i<n_momentum_dpoints; i++)
  {
    log_momentum.push_back(log(data_for_LZ_momentum[i]));
    log_range.push_back(log(data_for_LZ_range[i]));
    //cout << "Momentum: " << log_momentum[i] << " range: " << log_range[i] << endl;
  }
  double P_MeVc = exp(interp1D_ordered(log_range,log_momentum,log_z));
  //cout << "log_z: " <<  log_z << " interp: " 
  //     << interp1D_ordered(log_range,log_momentum,log_z) << " P_MeVc: " << P_MeVc << endl;

  // obtain attenuation lengths
  double out = 263.0 + 150*(P_MeVc/1000.0);
  return out;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets the parameters to those used by Granger and Smith 2000
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_Granger_parameters()
{
  //S_t = 1;

  // from Vermeesh 2007
  // 10Be from Chmeleff/Korschinek 10Be decay constant;
  lambda_10Be = 500e-9;    // in yr-1
  lambda_26Al = 980e-9;    // in yr-1
  lambda_14C = 121e-6;     // in yr-1
  lambda_36Cl = 230e-8;    // in yr-1

  // from Vermeesh 2007
  // These are calibrated to the Stone scaling
  // Also linke to the nishizumii standards
  // These come with Cosmocalc version 2.0
  // http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
  P0_10Be = 4.30;          // in a/g/yr
  P0_26Al = 31.10;         // in a/g/yr
  P0_14C = 15.21;          // in a/g/yr
  P0_36Cl = 58.95;         // in a/g/yr
  P0_21Ne = 18.23;         // in a/g/yr
  P0_3He = 121.59;         // in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 0.9724;
  F_10Be[1] = 0.0186;
  F_10Be[2] = 0.004;
  F_10Be[3] = 0.005;

  // dimensionless
  F_26Al[0] = 0.9655;
  F_26Al[1] = 0.0233;
  F_26Al[2] = 0.005;
  F_26Al[3] = 0.0062;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.0691;
  F_14C[2] = 0.0809;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.903;
  F_36Cl[1] = 0.0447;
  F_36Cl[2] = 0.05023;
  F_36Cl[3] = 0.0;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// function to set CRN parameters
// based on the Vermeesch approximation of the Schaller et al (2000)
// formulation
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_Schaller_parameters()
{
  //S_t =1;

  // from Vermeesh 2007
  // 10Be from Chmeleff/Korschinek 10Be decay constant;
  lambda_10Be = 500e-9;    // in yr-1
  lambda_26Al = 980e-9;    // in yr-1
  lambda_14C = 121e-6;     // in yr-1
  lambda_36Cl = 230e-8;    // in yr-1

  // from Vermeesh 2007
  // These are calibrated to the Stone scaling
  // Also linke to the nishizumii standards
  // These come with Cosmocalc version 2.0
  // http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
  P0_10Be = 4.30;          // in a/g/yr
  P0_26Al = 31.10;         // in a/g/yr
  P0_14C = 15.21;          // in a/g/yr
  P0_36Cl = 58.95;         // in a/g/yr
  P0_21Ne = 18.23;         // in a/g/yr
  P0_3He = 121.59;         // in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 0.964;
  F_10Be[1] = 0.0266;
  F_10Be[2] = -0.0074;
  F_10Be[3] = 0.0168;

  // dimensionless
  F_26Al[0] = 0.9575;
  F_26Al[1] = 0.0315;
  F_26Al[2] = -0.009;
  F_26Al[3] = 0.02;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.1363;
  F_14C[2] = 0.0137;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.903;
  F_36Cl[1] = 0.0793;
  F_36Cl[2] = 0.0177;
  F_36Cl[3] = 0.0;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets the parameters to those used by Braucher et al 2009
// as implemented in cosmocalc v2.0
// http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_Braucher_parameters()
{
  //S_t = 1;

  // from Vermeesh 2007
  // 10Be from Chmeleff/Korschinek 10Be decay constant;
  lambda_10Be = 500e-9;    // in yr-1
  lambda_26Al = 980e-9;    // in yr-1
  lambda_14C = 121e-6;     // in yr-1
  lambda_36Cl = 230e-8;    // in yr-1          

  // from Vermeesh 2007
  // These are calibrated to the Stone scaling
  // Also linke to the nishizumii standards
  // These come with Cosmocalc version 2.0
  // http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
  P0_10Be = 4.30;          // in a/g/yr
  P0_26Al = 31.10;         // in a/g/yr
  P0_14C = 15.21;          // in a/g/yr
  P0_36Cl = 58.95;         // in a/g/yr
  P0_21Ne = 18.23;         // in a/g/yr
  P0_3He = 121.59;         // in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 1500;
  Gamma[2] = 1500;
  Gamma[3] = 4320;

  // dimensionless
  F_10Be[0] = 0.9887;
  F_10Be[1] = 0.0027;
  F_10Be[2] = 0.0;
  F_10Be[3] = 0.0086;

  // dimensionless
  F_26Al[0] = 0.9699;
  F_26Al[1] = 0.0275;
  F_26Al[2] = 0.000;
  F_26Al[3] = 0.0026;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.15;
  F_14C[2] = 0.0;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.9456;
  F_36Cl[1] = 0.0324;
  F_36Cl[2] = 0.00;
  F_36Cl[3] = 0.022;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// 10Be is set to a new production curve provided by Shasta Marrero
// All others: sets the parameters to those used by Braucher et al 2009
// as implemented in cosmocalc v2.0
// http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_newCRONUS_parameters()
{
  //S_t = 1;

  // from Vermeesh 2007
  // 10Be from Chmeleff/Korschinek 10Be decay constant;
  lambda_10Be = 500e-9;    // in yr-1
  lambda_26Al = 980e-9;    // in yr-1
  lambda_14C = 121e-6;     // in yr-1
  lambda_36Cl = 230e-8;    // in yr-1          

  // from Data provided by Shasta mararreo for 10Be, everyting
  // else is from the Braucher constants
  //
  // All but 10Be are calibrated to the Stone scaling
  // Also linke to the nishizumii standards
  // These come with Cosmocalc version 2.0
  // http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
  
  // UPDATE The Al, Ne, He and 14C data come from the newcronus production files:
  // https://bitbucket.org/cronusearth/cronus-calc/src/bfdbded294a6432b9f4262a3859aee9fb9df9b18/de/production/physpars.m?at=master&fileviewer=file-view-default
  P0_10Be = 4.075213;          // in a/g/yr
  P0_26Al = 27.9;         // in a/g/yr
  P0_14C = 12.2;          // in a/g/yr
  P0_36Cl = 58.95;         // in a/g/yr
  P0_21Ne = 16.63;         // in a/g/yr
  P0_3He = 118.0;         // in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 1459.76761923;
  Gamma[2] = 11039.2402217;
  Gamma[3] = 4320;

  // ONLY 10Be reflects new CRONUS!!!!!!!
  // dimensionless
  F_10Be[0] = 0.98374;
  F_10Be[1] = 0.0137188126531;
  F_10Be[2] = 0.00252519164093;
  F_10Be[3] = 0.0;

  // dimensionless
  F_26Al[0] = 0.9699;
  F_26Al[1] = 0.0275;
  F_26Al[2] = 0.000;
  F_26Al[3] = 0.0026;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.15;
  F_14C[2] = 0.0;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.9456;
  F_36Cl[1] = 0.0324;
  F_36Cl[2] = 0.00;
  F_36Cl[3] = 0.022;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets the production rates to mimic those of CRONUS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_CRONUS_stone_parameters()
{
  //S_t =1;

  // from CRONUS calculator
  lambda_10Be = 499.75e-9;     // in yr-1
  lambda_26Al = 983e-9;     // in yr-1


  // from the CRONUS calculator
  // the numbers are divided by the scaling factor that determines neutron
  // spallation. The intended effect is to replicate the neutron production of 
  // the CRONUS calculator
  // IMPORTANT: these are scaled by the Granger scaling factors and not the 
  // scaller factors
  P0_10Be = 4.49/0.9724;          // in a/g/yr
  P0_26Al = 30.2922/0.9655;       // in a/g/yr
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this forces a neutron only calculation
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_Neutron_only_parameters()
{
  //S_t =1;

  // 10Be from Chmeleff/Korschinek 10Be decay constant;
  lambda_10Be = 500e-9;    // in yr-1
  lambda_26Al = 980e-9;    // in yr-1
  lambda_14C = 121e-6;     // in yr-1
  lambda_36Cl = 230e-8;    // in yr-1

  // from Vermeesh 2007
  P0_10Be = 4.30;          // in a/g/yr
  P0_26Al = 31.10;         // in a/g/yr
  P0_14C = 15.21;          // in a/g/yr
  P0_36Cl = 58.95;         // in a/g/yr
  P0_21Ne = 18.23;         // in a/g/yr
  P0_3He = 121.59;         // in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 738.6;
  Gamma[2] = 2688;
  Gamma[3] = 4360;

  // dimensionless
  F_10Be[0] = 1;
  F_10Be[1] = 0;
  F_10Be[2] = 0;
  F_10Be[3] = 0;

  // dimensionless
  F_26Al[0] = 1;
  F_26Al[1] = 0;
  F_26Al[2] = 0;
  F_26Al[3] = 0;

  // dimensionless
  F_14C[0] = 1;
  F_14C[1] = 0;
  F_14C[2] = 0;
  F_14C[3] = 0;

  // dimensionless
  F_36Cl[0] = 1;
  F_36Cl[1] = 0;
  F_36Cl[2] = 0;
  F_36Cl[3] = 0;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this gets difference in the fraction of spallation for different
// paris of muon production schemes
// key for pairs:
// 0 Braucher-Schaller
// 1 Braucher-Granger
// 2 Braucher-Neutron Only
// 3 Schaller-Granger
// 4 Schaller-Neutron Only
// 5 Granger-Neutron Only
vector<double> LSDCRNParameters::get_uncertainty_scaling_pair(int pair)
{
  double Be10_difference;
  double Al26_difference;
  
  vector<double> error_vec(2,0.0);
  
  if(pair == 0)
  {
    Be10_difference = 0.9887-0.9640;
    Al26_difference = 0.9669-0.9675;
  }
  else if(pair == 1)
  {
    Be10_difference = 0.9887-0.9744;
    Al26_difference = 0.9669-0.9625;
  }
  else if(pair == 2)
  {
    Be10_difference = 0.9887-1;
    Al26_difference = 0.9669-1;
  }
  else if(pair == 3)
  {
    Be10_difference = 0.9640-0.9744;
    Al26_difference = 0.9675-0.9625;
  }
  else if(pair == 4)
  {
    Be10_difference = 0.9640-1;
    Al26_difference = 0.9675-1;
  }  
  else if(pair == 4)
  {
    Be10_difference = 0.9744-1;
    Al26_difference = 0.9625-1;
  }
  else
  {
    cout << "LSDCRNP, line 1268, you did not supply a valud scaling pair" << endl
         << "Defaulting to Braucher-Schaller" << endl;
    Be10_difference = 0.9887-0.9640;
    Al26_difference = 0.9669-0.9675;
  }  
  
  error_vec[0] = Be10_difference;
  error_vec[1] = Al26_difference;  
  
  return error_vec;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function modifies the P0 values of Al and Be so that they
// mirror the uncertainty reported in CRONUS v2.2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParameters::set_P0_CRONUS_uncertainty_plus()
{
  // get the uncertainty percentage from the numbers in CRONUS
  // these numbers are from the stone scaling of CRONUS routine
  // make_al_be_consts_v22
  double CRONUS_P10 = 4.49;       // in atoms/g/year
  double CRONUS_delP10 = 0.39;
  
  double fraction_uncert = CRONUS_delP10/CRONUS_P10;
  
  // now update the 10Be and 26Al uncertainty values
  P0_10Be = P0_10Be*(1+fraction_uncert);
  P0_26Al = P0_26Al*(1+fraction_uncert);
  
  // get the uncertainty value. This gets passed on to gaussian error propagation
  vector<double> uncert_change;
  uncert_change.push_back(P0_10Be*fraction_uncert);
  uncert_change.push_back(P0_26Al*fraction_uncert);
  return uncert_change;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function modifies the P0 values of Al and Be so that they
// mirror the uncertainty reported in CRONUS v2.2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParameters::set_P0_CRONUS_uncertainty_minus()
{
  // get the uncertainty percentage from the numbers in CRONUS
  // these numbers are from the stone scaling of CRONUS routine
  // make_al_be_consts_v22
  double CRONUS_P10 = 4.49;       // in atoms/g/year
  double CRONUS_delP10 = 0.39;
  
  double fraction_uncert = CRONUS_delP10/CRONUS_P10;
  
  // now update the 10Be and 26Al uncertainty values
  P0_10Be = P0_10Be*(1-fraction_uncert);
  P0_26Al = P0_26Al*(1-fraction_uncert);
  
  // get the uncertainty value. This gets passed on to gaussian error propagation
  vector<double> uncert_change;
  uncert_change.push_back(P0_10Be*fraction_uncert);
  uncert_change.push_back(P0_26Al*fraction_uncert);
  return uncert_change;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Scaling from the Stone 2000 paper
// Units:
// latitude in decimal degrees
// pressure in hPa
// fsp is the fraction (between 0 and 1) of production at sea level
// and high latitude due to spallation (as opposed to muons).
// This argument is optional and defaults to 0.978, which is the value
// used by Stone (2000) for Be-10. The corresponding value for Al-26
// is 0.974. Note that using 0.844 for Be-10 and 0.826 for Al-26 will
// closely reproduce the Lal, 1991 scaling factors as long as the standard
// atmosphere is used to convert sample elevation to atmospheric pressure.
// Also note that this function will yield the scaling factor for spallation
// only when fsp=1, and that for muons only when fsp=0.
//
// IMPORTANT: This (and the Rc version) is probably the best scaling method!
// See https://cosmognosis.wordpress.com/2014/01/07/high-altitude-low-latitude-calibration-sites-i/
//
// Elevation can be converted to pressure with the functions
// stdatm.m (general use) and antatm.m (Antarctica).
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::stone2000sp(double lat,double P, double Fsp)
{
  if (Fsp > 1)
  {
    Fsp = 0.978;
  }
  
  if (fabs(lat) > 90)
  {
    cout << "Your latitude is > 90! Defaulting to 45 degrees" << endl;
    lat = 45;
  }

  // Spallogenic production at index latitudes;
  // enter constants from Table 1
  vector<double> a;
  a.push_back(31.8518);
  a.push_back(34.3699);
  a.push_back(40.3153);
  a.push_back(42.0983);
  a.push_back(56.7733);
  a.push_back(69.0720);
  a.push_back(71.8733);

  vector<double> b;
  b.push_back(250.3193);
  b.push_back(258.4759);
  b.push_back(308.9894);
  b.push_back(512.6857);
  b.push_back(649.1343);
  b.push_back(832.4566);
  b.push_back(863.1927);

  vector<double> c;
  c.push_back(-0.083393);
  c.push_back(-0.089807);
  c.push_back(-0.106248);
  c.push_back(-0.120551);
  c.push_back(-0.160859);
  c.push_back(-0.199252);
  c.push_back(-0.207069);

  vector<double> d;
  d.push_back(7.4260e-5);
  d.push_back(7.9457e-5);
  d.push_back(9.4508e-5);
  d.push_back(1.1752e-4);
  d.push_back(1.5463e-4);
  d.push_back(1.9391e-4);
  d.push_back(2.0127e-4);

  vector<double> e;
  e.push_back(-2.2397e-8);
  e.push_back(-2.3697e-8);
  e.push_back(-2.8234e-8);
  e.push_back(-3.8809e-8);
  e.push_back(-5.0330e-8);
  e.push_back(-6.3653e-8);
  e.push_back(-6.6043e-8);

  vector<double> ilats;
  ilats.push_back(0);
  ilats.push_back(10);
  ilats.push_back(20);
  ilats.push_back(30);
  ilats.push_back(40);
  ilats.push_back(50);
  ilats.push_back(60);

  // calculate index latitudes at given P's
  double lat0  = a[0] + (b[0] * exp(P/(-150.0))) + (c[0]*P) + (d[0]*(P*P)) + (e[0]*(P*P*P));
  double lat10 = a[1] + (b[1] * exp(P/(-150.0))) + (c[1]*P) + (d[1]*(P*P)) + (e[1]*(P*P*P));
  double lat20 = a[2] + (b[2] * exp(P/(-150.0))) + (c[2]*P) + (d[2]*(P*P)) + (e[2]*(P*P*P));
  double lat30 = a[3] + (b[3] * exp(P/(-150.0))) + (c[3]*P) + (d[3]*(P*P)) + (e[3]*(P*P*P));
  double lat40 = a[4] + (b[4] * exp(P/(-150.0))) + (c[4]*P) + (d[4]*(P*P)) + (e[4]*(P*P*P));
  double lat50 = a[5] + (b[5] * exp(P/(-150.0))) + (c[5]*P) + (d[5]*(P*P)) + (e[5]*(P*P*P));
  double lat60 = a[6] + (b[6] * exp(P/(-150.0))) + (c[6]*P) + (d[6]*(P*P)) + (e[6]*(P*P*P));

  vector<double> lat_at_specifics(7,0.0);
  lat_at_specifics[0] = lat0;
  lat_at_specifics[1] = lat10;
  lat_at_specifics[2] = lat20;
  lat_at_specifics[3] = lat30;
  lat_at_specifics[4] = lat40;
  lat_at_specifics[5] = lat50;
  lat_at_specifics[6] = lat60;

  //northernize southern-hemisphere inputs
  lat = fabs(lat);

  //set high lats to 60;
  if(lat > 60)
  {
    lat = 60.0;
  }

  // interpoloate elevation
  double S = interp1D_ordered(ilats,lat_at_specifics, lat);

  // Production by muons

  //constants
  vector<double> mk;
  mk.push_back(0.587);
  mk.push_back(0.600);
  mk.push_back(0.678);
  mk.push_back(0.833);
  mk.push_back(0.933);
  mk.push_back(1.000);
  mk.push_back(1.000);

  // index latitudes at given P's
  vector<double> m_index_at_given_P;
  m_index_at_given_P.push_back(mk[0]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[1]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[2]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[3]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[4]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[5]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[6]*exp( (1013.25-P)/242.0));
   
  // interpolate for actual elevation
  double M = interp1D_ordered(ilats, m_index_at_given_P, lat);

  //cout << "S: " << S << " M: " << M << " Fsp: " << Fsp << endl;

  // Combine spallogenic and muogenic production; return
  double Fm = 1 - Fsp;
  double out = ((S * Fsp) + (M * Fm));
  
  //cout << "Stone 2000 scaling is: "  << out << endl;

  return out;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function sets the total scaling value
// It is a product of the scaling, topographic shielding and snow shielding. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_scaling(double scaling, double topo_shield, double snow_shield)
{
  S_t = scaling*topo_shield*snow_shield;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function sets the total scaling value for use with neutron only 
// It is a product of the scaling, topographic shielding and snow shielding. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::set_neutron_scaling(double scaling, double topo_shield, 
                                           double snow_shield)
{
  neutron_S_t = scaling*topo_shield*snow_shield;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function takes a single scaling factor for
// elevation scaling, self shielding, snow shielding,
// and latitude scaling and produces scaling factors
// for each production mechamism.
// the scaling follows the approach of vermeesch 2008
// it uses a 'virstual' shielding depth to calcualte
// the updated scaling factors
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::scale_F_values(double single_scaling)
{

  cout << "Line 1411 LSDCRNParameters THIS scaling is very slow! I suggest " 
       << "using the Newton-Raphson version!" << endl;

  double tol = 1e-7;
  double x = 0;
  double new_x = 0;
  double test_scaling = 1e8;
  double dx =-10;

  // first do 10Be
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }


  //int iterations = 0;
  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    //iterations ++;
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_10Be[0]+
                   exp(-new_x/Gamma[1])*F_10Be[1]+
                   exp(-new_x/Gamma[2])*F_10Be[2]+
                   exp(-new_x/Gamma[3])*F_10Be[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }
  //cout << "Iterations were: " << iterations << endl;

  // now reset the F_values
  F_10Be[0] = exp(-new_x/Gamma[0])*F_10Be[0];
  F_10Be[1] = exp(-new_x/Gamma[1])*F_10Be[1];
  F_10Be[2] = exp(-new_x/Gamma[2])*F_10Be[2];
  F_10Be[3] = exp(-new_x/Gamma[3])*F_10Be[3];

  //cout << "============================================================" << endl;
  //cout << "LSDCRNP OLD SCALING" << endl;
  //cout << "FINISHED 10Be x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_10Be[0] << endl << F_10Be[1] << endl << F_10Be[2] << endl << F_10Be[3] << endl;
  //cout << "============================================================" << endl;
  //cout << "Total scaling is: " << single_scaling << " and sum is: " 
  //     << F_10Be[0]+F_10Be[1]+F_10Be[2]+F_10Be[3] << endl;

  // now do 26Al
  x = 0;
  new_x = 0;
  test_scaling = 1e8;
  dx =-10;
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_26Al[0]+
                   exp(-new_x/Gamma[1])*F_26Al[1]+
                   exp(-new_x/Gamma[2])*F_26Al[2]+
                   exp(-new_x/Gamma[3])*F_26Al[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_26Al[0] = exp(-new_x/Gamma[0])*F_26Al[0];
  F_26Al[1] = exp(-new_x/Gamma[1])*F_26Al[1];
  F_26Al[2] = exp(-new_x/Gamma[2])*F_26Al[2];
  F_26Al[3] = exp(-new_x/Gamma[3])*F_26Al[3];

  //cout << "FINISHED 26Al x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_26Al[0] << endl << F_26Al[1] << endl << F_26Al[2] << endl << F_26Al[3] << endl;

  // now do 36Cl
  x = 0;
  new_x = 0;
  test_scaling = 1e8;
  dx =-10;
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_36Cl[0]+
                   exp(-new_x/Gamma[1])*F_36Cl[1]+
                   exp(-new_x/Gamma[2])*F_36Cl[2]+
                   exp(-new_x/Gamma[3])*F_36Cl[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_36Cl[0] = exp(-new_x/Gamma[0])*F_36Cl[0];
  F_36Cl[1] = exp(-new_x/Gamma[1])*F_36Cl[1];
  F_36Cl[2] = exp(-new_x/Gamma[2])*F_36Cl[2];
  F_36Cl[3] = exp(-new_x/Gamma[3])*F_36Cl[3];

  //cout << "FINISHED 36Cl x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_36Cl[0] << endl << F_36Cl[1] << endl << F_36Cl[2] << endl << F_36Cl[3] << endl;

  // now do 14C
  x = 0;
  new_x = 0;
  test_scaling = 1e8;
  dx =-10;
  if (single_scaling > 1)
  {
    dx = -10;
  }
  else if (single_scaling < 1)
  {
    dx = 10;
  }
  else if (single_scaling == 1)
  {
    dx = 10;
    test_scaling = 1;
  }

  while (fabs(test_scaling - single_scaling) > tol)
  //for (int i = 0; i< 100; i++)
  {
    x = new_x;
    new_x = x+dx;

    // calcualte the new scaling
    test_scaling = exp(-new_x/Gamma[0])*F_14C[0]+
                   exp(-new_x/Gamma[1])*F_14C[1]+
                   exp(-new_x/Gamma[2])*F_14C[2]+
                   exp(-new_x/Gamma[3])*F_14C[3];

    // if you have overshot, halve dx and reset new_x
    if (single_scaling > 1)
    {
      if (test_scaling > single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
    else if (single_scaling < 1)
    {
      if (test_scaling < single_scaling)
      {
        dx = 0.5*dx;
        new_x = x;
      }
    }
  }

  // now reset the F_values
  F_14C[0] = exp(-new_x/Gamma[0])*F_14C[0];
  F_14C[1] = exp(-new_x/Gamma[1])*F_14C[1];
  F_14C[2] = exp(-new_x/Gamma[2])*F_14C[2];
  F_14C[3] = exp(-new_x/Gamma[3])*F_14C[3];

  //cout << "FINISHED 14C x is: " << x << " and test_scaling is: " << test_scaling << endl;
  //cout << F_14C[0] << endl << F_14C[1] << endl << F_14C[2] << endl << F_14C[3] << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function is similar to the scaling function for total nuclides but
// it allows the user to pick the nuclides they want scaled, and
// also uses a faster newton raphson iterator to get the correct scaling
//
// The bool vector nuclides_for_scaling has four elements
// nuclides_for_scaling[0] = true: calculate 10Be
// nuclides_for_scaling[1] = true: calculate 26Al
// nuclides_for_scaling[2] = true: calculate 36Cl
// nuclides_for_scaling[3] = true: calculate 14C
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::scale_F_values(double single_scaling, vector<bool> nuclides_for_scaling)
{
  // first check the boolean vector. If it is an incorrect size, print a warning
  // and then default to all true
  if(nuclides_for_scaling.size() != 4)
  {
    cout << "nulides for scaling vector is the wrong size! Defaulting to a " << endl;
    cout << "vector that calculates all nuclides" << endl;
    vector<bool> temp_vec(4,true);
    nuclides_for_scaling = temp_vec;
  }

  // set up the parameters for the newton-raphson iteration
  double initial_guess = 0;         // an initial test depth
  double new_x;                 // the new test depth
  double displace_x = 1e-6;     // how far you diplace the test depth
  double displace_scaling;      // the scaling after displacement
  double scaling_this_step;     // the scaling at the current step
  double fx;                    // function for Newton Raphson to find root
  double fx_displace;           // function after displacement
  double fx_derivative;         // derivative of the root finding function
  double x_change;              // the change in the scaling
  vector<double> F(4,0.0);      // holds the F values, cycles between nuclides
  double tolerance = 1e-7;      // the tolerance over which the scaling can change

  // go through the nuclides, selecting each based on the booleans fed to 
  // the routine  
  // first 10Be  
  if(nuclides_for_scaling[0])
  {
    // replace the F values
    F[0] =  F_10Be[0];
    F[1] =  F_10Be[1];
    F[2] =  F_10Be[2];
    F[3] =  F_10Be[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_10Be[0] = exp(-new_x/Gamma[0])*F_10Be[0];
    F_10Be[1] = exp(-new_x/Gamma[1])*F_10Be[1];
    F_10Be[2] = exp(-new_x/Gamma[2])*F_10Be[2];
    F_10Be[3] = exp(-new_x/Gamma[3])*F_10Be[3];  
  }

  // now 26Al
  if(nuclides_for_scaling[1])
  {
    // replace the F values
    F[0] =  F_26Al[0];
    F[1] =  F_26Al[1];
    F[2] =  F_26Al[2];
    F[3] =  F_26Al[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_26Al[0] = exp(-new_x/Gamma[0])*F_26Al[0];
    F_26Al[1] = exp(-new_x/Gamma[1])*F_26Al[1];
    F_26Al[2] = exp(-new_x/Gamma[2])*F_26Al[2];
    F_26Al[3] = exp(-new_x/Gamma[3])*F_26Al[3];  
  }

  // now 36Cl
  if(nuclides_for_scaling[2])
  {
    // replace the F values
    F[0] =  F_36Cl[0];
    F[1] =  F_36Cl[1];
    F[2] =  F_36Cl[2];
    F[3] =  F_36Cl[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_36Cl[0] = exp(-new_x/Gamma[0])*F_36Cl[0];
    F_36Cl[1] = exp(-new_x/Gamma[1])*F_36Cl[1];
    F_36Cl[2] = exp(-new_x/Gamma[2])*F_36Cl[2];
    F_36Cl[3] = exp(-new_x/Gamma[3])*F_36Cl[3];  
  }

  // now 14C
  if(nuclides_for_scaling[3])
  {
    // replace the F values
    F[0] =  F_14C[0];
    F[1] =  F_14C[1];
    F[2] =  F_14C[2];
    F[3] =  F_14C[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_14C[0] = exp(-new_x/Gamma[0])*F_14C[0];
    F_14C[1] = exp(-new_x/Gamma[1])*F_14C[1];
    F_14C[2] = exp(-new_x/Gamma[2])*F_14C[2];
    F_14C[3] = exp(-new_x/Gamma[3])*F_14C[3];  
  }


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function prints F values to screen for bug checking
// The bool vector nuclides_for_scaling has four elements
// nuclides_for_scaling[0] = true: calculate 10Be
// nuclides_for_scaling[1] = true: calculate 26Al
// nuclides_for_scaling[2] = true: calculate 36Cl
// nuclides_for_scaling[3] = true: calculate 14C
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::print_F_values_to_screen(vector<bool> nuclides_for_scaling)
{

  if(nuclides_for_scaling[0])
  {
    cout << "10Be F\t"<<F_10Be[0]<<"\t"<<F_10Be[1]<<"\t"<<F_10Be[2]<<"\t"<<F_10Be[3]<<endl;
  }
  if(nuclides_for_scaling[1])
  {
    cout << "26Al F\t"<<F_26Al[0]<<"\t"<<F_26Al[1]<<"\t"<<F_26Al[2]<<"\t"<<F_26Al[3]<<endl;
  }
  if(nuclides_for_scaling[2])
  {
    cout << "36Cl F\t"<<F_36Cl[0]<<"\t"<<F_36Cl[1]<<"\t"<<F_36Cl[2]<<"\t"<<F_36Cl[3]<<endl;
  }
  if(nuclides_for_scaling[3])
  {
    cout << "14C F\t"<<F_14C[0]<<"\t"<<F_14C[1]<<"\t"<<F_14C[2]<<"\t"<<F_14C[3]<<endl;
  }


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints the CRN parameters to screen for bug checking
// nuclides_for_scaling[0] = true: calculate 10Be
// nuclides_for_scaling[1] = true: calculate 26Al
// nuclides_for_scaling[2] = true: calculate 36Cl
// nuclides_for_scaling[3] = true: calculate 14C
void LSDCRNParameters::print_parameters_to_screen(vector<bool> nuclides_for_scaling)
{
  cout << "=======================================================" << endl;
  cout << "Cosmogenic parameters" << endl;
  cout << "Gammas: \t" << Gamma[0] << "\t"<< Gamma[1]<< "\t"<<Gamma[2]<<"\t"<<Gamma[3]<<endl;
  

  if(nuclides_for_scaling[0])
  {
    cout << "+++++++++++" << endl;
    cout << "P0 10Be: " << P0_10Be << " lambda 10Be: " << lambda_10Be << endl;
    cout << "10Be F\t"<<F_10Be[0]<<"\t"<<F_10Be[1]<<"\t"<<F_10Be[2]<<"\t"<<F_10Be[3]<<endl;
  }
  if(nuclides_for_scaling[1])
  {
    cout << "+++++++++++" << endl;
    cout << "P0 26Al: " << P0_26Al << " lambda 26Al: " << lambda_26Al << endl;
    cout << "26Al F\t"<<F_26Al[0]<<"\t"<<F_26Al[1]<<"\t"<<F_26Al[2]<<"\t"<<F_26Al[3]<<endl;
  }
  if(nuclides_for_scaling[2])
  {
    cout << "+++++++++++" << endl;
    cout << "P0 36Cl: " << P0_36Cl << " lambda 36Cl: " << lambda_36Cl << endl;
    cout << "36Cl F\t"<<F_36Cl[0]<<"\t"<<F_36Cl[1]<<"\t"<<F_36Cl[2]<<"\t"<<F_36Cl[3]<<endl;
  }
  if(nuclides_for_scaling[3])
  {
    cout << "+++++++++++" << endl;
    cout << "P0 14C: " << P0_14C << " lambda 14C: " << lambda_14C << endl;
    cout << "14C F\t"<<F_14C[0]<<"\t"<<F_14C[1]<<"\t"<<F_14C[2]<<"\t"<<F_14C[3]<<endl;
  }
  cout << "==========================================================="  << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function converts elevation to atmospheric pressure
// It follows the cronus calculator NCEPatm_2 function
// written by Greg Balco
// SMM
// 4/12/2014
//
// Looks up surface pressure and 1000 mb temp from NCEP reanalysis
// and calculates site atmospheric pressures using these as inputs to the
// standard atmosphere equation. 
//
// Syntax: pressure = NCEPatm_2(site_lat,site_lon,site_elv);
// 
// Requires:
///       site_lat: latitude (DD). Southern hemisphere is negative.
//       site_lon: longitude (DD). Western hemisphere is negative.
//           Tries to deal with 0-360 longitudes gracefully.
//       site_elv: elevation (m).
//
// Returns site pressure in hPa.
//
// Vectorized. Send vectors of equal length.
//
// Note: this must load the data file NCEP2.mat whenever called. 
// Repeated calls to this function will be slow for this reason. 
//
// Also: This function is OK but not great for Antarctica.
// Use antatm.m instead. 
//
// Remember: it is always better to estimate the average pressure at your 
// site using a pressure-altitude relation obtained from nearby station
// data.
//
// Original m code Written by Greg Balco -- UW Cosmogenic Nuclide Lab
// balcs@u.washington.edu
// October, 2007
// Part of the CRONUS-Earth online calculators: 
//      http://hess.ess.washington.edu/math
//
// Copyright 2001-2007, University of Washington
// All rights reserved
// Developed in part with funding from the National Science Foundation.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double LSDCRNParameters::NCEPatm_2(double site_lat, double site_lon, double site_elev)
{
  // deal with negative longitudes
  if(site_lon < 0)
  {
    site_lon = site_lon+360.0;
  }
  
  //cout << "LSDCRNP, line 1618, Site lat: " << site_lat << " and long: " << site_lon << endl;
  
  // check to see if data is loaded:
  if (int(gm_hgt.size()) != 8)
  {
    string path_to_data;
    cout << "You didn't load the NCEP data. Doing that now. " << endl;
    cout << "Enter path to data files: " << endl;
    cin >> path_to_data;
    load_parameters_for_atmospheric_scaling(path_to_data);
  }
  
  // now, interpolate sea level pressure and temperature
  //cout << "interpolating pressure " << endl;
  double site_slp = interp2D_bilinear(NCEPlat, NCEPlon, meanslp, 
                        site_lat, site_lon);
  //cout << "Did pressure, now temp:" << endl;                      
  double site_T = interp2D_bilinear(NCEPlat, NCEPlon, meant1000, 
                        site_lat, site_lon);
  //cout << "Did temp" << endl;                      
  
  
  double site_T_degK = site_T + 273.15;

  // Some More parameters
  double gmr = -0.03417; // Assorted constants (this has come from Greg Balco's code)
  double dtdz = 0.0065;  // Lapse rate from standard atmosphere 
  
  // Calculate site pressure using the site-specific SLP and T1000 with the
  // standard atmosphere equation.

  //cout << site_T_degK << endl;
  //cout << "Log term: " << log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) << endl;
  //cout << "Exp term: " << exp( (gmr/dtdz)*( log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) ) ) << endl;

  double out = site_slp*exp( (gmr/dtdz)*( log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) ) );
  
  //cout << endl;
  //cout << "Site sea level pressure: " << site_slp << " and site Temp: "<< site_T 
  //     << " and pressure: " << out << endl << endl <<endl;
  return out;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function gets the spallation attenuation lenth in g/cm^2
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double LSDCRNParameters::get_spallation_attenuation_length(bool use_CRONUS)
{
  
  double att_length = 160;
  if(use_CRONUS == true)
  {
    // first check to see if CRONUS data maps are set
    if(CRONUS_data_map.find("l10") == CRONUS_data_map.end())
    {
      //cout << "You haven't set the CRONUS data map. I'm doing that for you now!" << endl;
      set_CRONUS_data_maps();
    }
    att_length = CRONUS_data_map["Lsp"];    
  
  }
  else
  {
    att_length = Gamma[0];
  }

  return att_length;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This functioon gets the cosmogenic decay coefficients
// for 10Be and 26Al
// It has a boolean to choose the CRONUS or Cosmocalc values
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> LSDCRNParameters::get_decay_coefficients(bool use_CRONUS)
{
  
  vector<double> decay_constants(2,0.0);
  if(use_CRONUS == true)
  {
    // first check to see if CRONUS data maps are set
    if(CRONUS_data_map.find("l10") == CRONUS_data_map.end())
    {
      //cout << "You haven't set the CRONUS data map. I'm doing that for you now!" << endl;
      set_CRONUS_data_maps();
    }  
    decay_constants[0] = CRONUS_data_map["l10"];
    decay_constants[1] = CRONUS_data_map["l26"];
      
    //cout << "Getting lambda, it is for 10Be: " << CRONUS_data_map["l10"] << endl;
        
  
  }
  else
  {
    decay_constants[0] = lambda_10Be;
    decay_constants[1] = lambda_26Al;
  }

  return decay_constants;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function calculates the CRONUS version of the mu production vector
// takes atmospheric pressure in HPa
// and the effective depth in g/cm^2
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDCRNParameters::get_CRONUS_P_mu_vectors(double pressure, double sample_effective_depth,
                               vector<double>& z_mu, vector<double>& P_mu_z_10Be,
                               vector<double>&  P_mu_z_26Al)
{
  // parameters for setting up the log vector
  int n_z = 101;
  double end_log = 5.3;
  double spacing = end_log/(double(n_z-2));
  double this_log10;
  
  // set up vectors to hold production rates
  vector<double> zz_mu;
  vector<double> zP_mu_z_10Be;
  vector<double> zP_mu_z_26Al;
  
  // get the top values of the vectors
  double this_z = 0.5*sample_effective_depth;
  double this_P_mu_10Be =0;
  double this_P_mu_26Al =0;
  P_mu_total_return_nuclides(this_z,pressure, this_P_mu_10Be, this_P_mu_26Al);
  
  zz_mu.push_back(this_z);
  zP_mu_z_10Be.push_back(this_P_mu_10Be);
  zP_mu_z_26Al.push_back(this_P_mu_26Al);
  
  // first build the vector of depths
  for(int i = 0; i<n_z-1; i++)
  {
    this_log10 = double(i)*spacing;
    this_z = pow(10.0,this_log10)+0.5*sample_effective_depth;
    P_mu_total_return_nuclides(this_z,pressure, this_P_mu_10Be, this_P_mu_26Al);
    
    zz_mu.push_back(this_z);
    zP_mu_z_10Be.push_back(this_P_mu_10Be);
    zP_mu_z_26Al.push_back(this_P_mu_26Al); 
  }

  // we then need to take away some component of z_mu
  // so this is matched in the ET objective function
  // (see section 4 in get_al_be_erosion.m)
  int n_nodes = int(zz_mu.size());
  for(int i = 0; i<n_nodes; i++)
  {
    zz_mu[i] =  zz_mu[i] - sample_effective_depth*0.5;
  }
  
  z_mu=zz_mu;
  P_mu_z_10Be = zP_mu_z_10Be;
  P_mu_z_26Al = zP_mu_z_26Al;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Integrate muons flux as as a function of erosion rate
// Note this is different from CRONUS calculator since this uses Simpson's
// Rule whereas CRONUS uses trapezoid rule
//
// This seems to produce much more muon activity than the vermeesh/grange/schaller
// equations!!
// Will need to work out why. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::integrate_muon_flux_for_erosion(double E, 
                           vector<double> z_mu, vector<double> P_mu_10Be,
                           vector<double> P_mu_26Al,
                           double& Be10_mu_N, double& Al26_mu_N)
{                          
  
  if(E == 0)
  {
    cout << "Warning, you are about to divide by zero in" << endl
         << "LSDCRNParameters::integrate_muon_flux_for_erosion" << endl;
  }
  
  // get the decay coefficients
  bool use_CRONUS = true;
  vector<double> decay_coeff = get_decay_coefficients(use_CRONUS);
  
  // get the time vector, given by dividing depth by erosion. 
  // !!IMPORTANT Erosion is in g/cm^2/yr!!
  int n_z_nodes = int(z_mu.size());
  vector<double> t_mu(n_z_nodes,0.0);
  for(int i = 0; i<n_z_nodes; i++)
  {
    t_mu[i] = z_mu[i]/E;
    //cout << "t_mu["<<i+1<<"]: " << t_mu[i] << endl;
  
  }

  // now integrate over this vector using simpsons rule
  double a,b;
  double fa10,fb10;
  double fa26,fb26;
  
  
  /*
  // set up the locations to evaluate the integral
  b = t_mu[0];
  
  double intermediate;
  double fi10,fi26;
  double sum10,sum26;
  
  fb10 = P_mu_10Be[0]*(exp(-decay_coeff[0]*t_mu[0]));
  fb26 = P_mu_26Al[0]*(exp(-decay_coeff[1]*t_mu[0]));
  sum10 = 0;
  sum26 = 0;
  for(int i = 1; i< n_z_nodes; i++)
  {
    // locations of spacings
    a = b;
    b = t_mu[i]; 
    intermediate = (b+a)/2.0;
    
    // functions evaluated at spacings
    fa10 = fb10;
    fa26 = fb26;
    
    fi10 = P_mu_10Be[i]*(exp(-decay_coeff[0]*intermediate));
    fi26 = P_mu_26Al[i]*(exp(-decay_coeff[1]*intermediate));
    
    fb10 = P_mu_10Be[i]*(exp(-decay_coeff[0]*b));
    fb26 = P_mu_26Al[i]*(exp(-decay_coeff[1]*b));
    
    sum10+= ((b-a)/6.0)*(fa10+4.0*fi10+fb10);
    sum26+= ((b-a)/6.0)*(fa26+4.0*fi26+fb26);
  }
  */
  
  // for error checking, try trapezoid rule: this is not as accurate as simpsons
  // but is used by CRONUS calculator
  double sum_trap10 = 0;
  double sum_trap26 = 0;
  b = t_mu[0];
  
  fb10 = P_mu_10Be[0]*(exp(-decay_coeff[0]*t_mu[0]));
  fb26 = P_mu_26Al[0]*(exp(-decay_coeff[1]*t_mu[0]));
  for(int i = 1; i< n_z_nodes; i++)
  {
    // locations of spacings
    a = b;
    b = t_mu[i]; 
    
    // functions evaluated at spacings
    fa10 = fb10;
    fa26 = fb26;
    
    fb10 = P_mu_10Be[i]*(exp(-decay_coeff[0]*b));
    fb26 = P_mu_26Al[i]*(exp(-decay_coeff[1]*b));
    
    sum_trap10+= (b-a)*0.5*(fb10+fa10);
    sum_trap26+= (b-a)*0.5*(fb26+fa26);
  }
  //cout << "Simpsons  10: " << sum10 << " 26: " << sum26 << endl;
  //cout << "Trapezoid 10: " << sum_trap10 << " 26: " << sum_trap26 << endl;
  
  
  // now return the integral
  Be10_mu_N = sum_trap10;
  Al26_mu_N = sum_trap26;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function gets the spallation flux for a given erosion rate
// for non-time dependant production
// It is used in the CRONUS emulator
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::integrate_nonTD_spallation_flux_for_erosion(double E, 
                                   double thick_SF,
                                   double P_sp_10Be, double P_sp_26Al, 
                                   double& Be10_sp_N,double& Al26_sp_N)
{
  // get analytical solution of spallation production
  bool use_CRONUS = true;
  
  vector<double> decay = get_decay_coefficients(use_CRONUS);
  double A_10Be =  decay[0]+E/get_spallation_attenuation_length(use_CRONUS);
  double A_26Al =  decay[1]+E/get_spallation_attenuation_length(use_CRONUS);

  //cout << "Integrating spallation, LINE 1719" << endl;
  //cout << "ThickSF: " << thick_SF << endl;
  //cout << "Psp10: " << P_sp_10Be << " Psp26: " << P_sp_26Al << endl;
  //cout << "decay10: " << decay[0] << " decay26: " << decay[1] << endl;
  //cout << "A10: " << A_10Be << " A26: " << A_26Al << endl;
  //cout << "thick/A: " <<  thick_SF/A_10Be << endl;
  
  Be10_sp_N = thick_SF*P_sp_10Be/A_10Be;
  Al26_sp_N = thick_SF*P_sp_26Al/A_26Al;

} 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//  This calculates some muogenic paramters for the uncertainty analysis
// returns the uncertainty parameters in the form of a vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParameters::CRONUS_get_muon_uncertainty_params(double pressure)
{
  // first check to see if CRONUS data maps are set
  if(CRONUS_data_map.find("l10") == CRONUS_data_map.end())
  {
    cout << "You haven't set the CRONUS data map. I'm doing that for you now!" << endl;
    set_CRONUS_data_maps();
  }
  
  // now get the production parameters
  double test_elev = 0.0;
  vector<double> muon_prod 
             = calculate_muon_production_CRONUS(test_elev, pressure);
  
  // calculate the parameters
  double delPfast_10 = muon_prod[0]*(CRONUS_data_map["delsigma190_10"]/
                                     CRONUS_data_map["sigma190_10"]);
  double delPfast_26 = muon_prod[1]*(CRONUS_data_map["delsigma190_26"]/
                                     CRONUS_data_map["sigma190_26"]);
  double delPneg_10 = muon_prod[2]*(CRONUS_data_map["delk_neg10"]/
                                     CRONUS_data_map["k_neg10"]);
  double delPneg_26 = muon_prod[3]*(CRONUS_data_map["delk_neg26"]/
                                     CRONUS_data_map["k_neg26"]);;
  double delPmu0_10 = sqrt(delPfast_10*delPfast_10 + delPneg_10*delPneg_10);
  double delPmu0_26 = sqrt(delPfast_26*delPfast_26 + delPneg_26*delPneg_26);
  
  // return the vector with the paramters
  vector<double> uncert_params(6,0.0);
  uncert_params[0]=delPfast_10;
  uncert_params[1]=delPfast_26;
  uncert_params[2]=delPneg_10;
  uncert_params[3]=delPneg_26;
  uncert_params[4]=delPmu0_10;
  uncert_params[5]=delPmu0_26;
  
  return uncert_params;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// 
// This function gets relative scaling errors
// Used in the CRONUS uncertanty analysis
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<double> LSDCRNParameters::CRONUS_get_uncert_production_ratios(string scaling_name)
{
  // first check to see if CRONUS data maps are set
  if(CRONUS_data_map.find("l10") == CRONUS_data_map.end())
  {
    cout << "You haven't set the CRONUS data map. I'm doing that for you now!" << endl;
    set_CRONUS_data_maps();
  }

  // vector to hold the uncertainty of relative production
  vector<double> rel_delP(2,0.0);

  // get the strings for entry into the data map
  string P10_map_string = "P10_ref_"+scaling_name;
  string delP10_map_string = "delP10_ref_"+scaling_name;
  string P26_map_string = "P26_ref_"+scaling_name;
  string delP26_map_string = "delP26_ref_"+scaling_name;
  if (scaling_name == "St" || scaling_name == "Lm" || scaling_name == "Du" ||
      scaling_name == "Li" || scaling_name == "De")
  {
    rel_delP[0] = CRONUS_data_map[delP10_map_string]/CRONUS_data_map[P10_map_string];
    rel_delP[1] = CRONUS_data_map[delP26_map_string]/CRONUS_data_map[P26_map_string];
  }
  else
  {
    cout << "Warning: you have not selected a valid scaling name, options are: " << endl;
    cout << "St  for Stone scaling" << endl;
    cout << "Li  for Lifton scaling" << endl;
    cout << "De  for Deselits scaling" << endl;
    cout << "Du  for Dunai scaling" << endl;
    cout << "Lm  for Lal magnetic scaling" << endl;
    cout << "You selected: " << scaling_name << endl;
  }
  return rel_delP;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints parameters to allow both checking and tracking of the
// parameters used in calculations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::print_parameters_to_file(string fname, string muon_scaling)
{
  ofstream CRNparams_out;
  CRNparams_out.open(fname.c_str());
  
  CRNparams_out << "CRNParameters version: " << version << endl;
  
  
  if(muon_scaling == "Braucher")
  {
    set_Braucher_parameters();
    
  }
  else if (muon_scaling == "Schaller")
  {
    set_Schaller_parameters();
  }
  else if (muon_scaling == "Granger")
  {
    set_Granger_parameters();
  }
  else
  {
    muon_scaling = "Braucher_as_default";
    set_Braucher_parameters();
  }
  CRNparams_out << "Muon scaling: " << muon_scaling;
  
  CRNparams_out << "lambda_10Be: " << lambda_10Be << " yr^-1" << endl;
  CRNparams_out << "lambda_26Al: " << lambda_26Al << " yr^-1" << endl;
  CRNparams_out << "lambda_14C: " << lambda_14C << " yr^-1" << endl;
  CRNparams_out << "lambda_36Cl: " << lambda_36Cl << " yr^-1" << endl;
  set_CRONUS_data_maps();

  CRNparams_out << "P0_10Be: " << P0_10Be << " a/g/yr, delP0_10Be: " << CRONUS_data_map["delP10_ref_St"] << endl;
  CRNparams_out << "P0_26Al: " << P0_26Al << " a/g/yr, delP0_26Al: " << CRONUS_data_map["del26Al_ref_St"] << endl;
  CRNparams_out << "P0_14C: " << P0_14C << " a/g/yr, delP0_14C: " << CRONUS_data_map["del14C_ref_St"] << endl;
  CRNparams_out << "P0_36Cl: " << P0_36Cl << " a/g/yr, delP0_36Cl: " << CRONUS_data_map["del36Cl_ref_St"] << endl;
  CRNparams_out << "P0_21Ne: " << P0_21Ne << " a/g/yr, delP0_21Ne: " << CRONUS_data_map["del21Ne_ref_St"] << endl;
  CRNparams_out << "P0_3He: " << P0_3He << " a/g/yr, delP0_3He: " << CRONUS_data_map["del3He_ref_St"] << endl;
  
  CRNparams_out << "Gamma[0]: " << Gamma[0] << " g/cm^2" << endl;
  CRNparams_out << "Gamma[1]: " << Gamma[1] << " g/cm^2" << endl;
  CRNparams_out << "Gamma[2]: " << Gamma[2] << " g/cm^2" << endl;
  CRNparams_out << "Gamma[3]: " << Gamma[3] << " g/cm^2" << endl;

  CRNparams_out << "10Be F\t"<<F_10Be[0]<<"\t"<<F_10Be[1]<<"\t"<<F_10Be[2]<<"\t"<<F_10Be[3]<<endl;
  CRNparams_out << "26Al F\t"<<F_26Al[0]<<"\t"<<F_26Al[1]<<"\t"<<F_26Al[2]<<"\t"<<F_26Al[3]<<endl;
  CRNparams_out << "36Cl F\t"<<F_36Cl[0]<<"\t"<<F_36Cl[1]<<"\t"<<F_36Cl[2]<<"\t"<<F_36Cl[3]<<endl;
  CRNparams_out << "14C F\t"<<F_14C[0]<<"\t"<<F_14C[1]<<"\t"<<F_14C[2]<<"\t"<<F_14C[3]<<endl;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This  function prints out the production rates of a vector of depths
// it is used to compare production rates amongst different production schemes
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDCRNParameters::Print_10Beproduction_csv(string filename, string path_to_atmospheric_data)
{
  // open the file
  ofstream prod_file_out;
  prod_file_out.open(filename.c_str());
  prod_file_out << "Production rates from muons, at 70N, 0W, 0 elevation."  << endl;
  prod_file_out << "All prduction rates in atoms/g/yr" << endl;
  prod_file_out << "eff_depth(g/cm^2),mu_CRONUS_prod,mu_Schaller_prod,mu_Granger_prod,"
                << "mu_Braucher_prod,mu_newCRONUS_prod," 
                << "total_CRONUS_prod,total_Schaller_prod,total_Granger_prod,"
                << "total_Braucher_prod,total_newCRONUS_prod" << endl;
  
  // we are going to print these at a fixed location
  double site_lat = 70;
  double site_lon = 0;
  double site_elev = 0;

  // get the atmospheric parameters
  load_parameters_for_atmospheric_scaling(path_to_atmospheric_data);
  set_CRONUS_data_maps();
  double pressure = NCEPatm_2(site_lat, site_lon, site_elev);
  
  // calculate the scaling. We assume no topographic, self or snow shielding
  double Fsp = 1.0;
  double this_P = stone2000sp(site_lat,pressure, Fsp);

  // get some parameters for Stone production
  vector<double> Prefs_st = get_Stone_Pref();

  // retrieve the pre-scaling factors
  double P_ref_St_10 = Prefs_st[0];
  double Fsp_CRONUS = CRONUS_data_map["Fsp10"];

  // also get the P0 for the CRONUS spallation
  double CRONUS_P = this_P*P_ref_St_10*Fsp_CRONUS;
  bool use_CRONUS = true;
  double Gamma_CRONUS = get_spallation_attenuation_length(use_CRONUS);
  
  cout << "The pressure is: " << pressure << endl;
  
  // this tells the scaling routine to on;y scale for 10Be;
  vector<bool> nuclides_for_scaling(4,false);
  nuclides_for_scaling[0]= true;
  
  vector<double> z_mu;
  vector<double> P_mu_z_10Be;
  vector<double> P_total_z_10Be;
  vector<double> P_mu_z_26Al;
  
  vector<double> P_mu_schaller;
  vector<double> P_mu_granger;
  vector<double> P_mu_braucher;
  vector<double> P_mu_newCRONUS;
  
  vector<double> P_total_schaller;
  vector<double> P_total_granger;
  vector<double> P_total_braucher;
  vector<double> P_total_newCRONUS;
  
  // get the production vectors from the CRONUS scheme
  double sample_effective_depth = 0;
  get_CRONUS_P_mu_vectors(pressure, sample_effective_depth, z_mu, P_mu_z_10Be,P_mu_z_26Al);
  
  // now the production vectors from the other schemes
  int n_z = int(z_mu.size());
  for(int i = 0; i<n_z; i++)
  {
    double z = z_mu[i];
    
    // get the CRONUS total production rates
    double this_total_CRONUS = P_mu_z_10Be[i] + CRONUS_P*exp(-z/Gamma_CRONUS);
    P_total_z_10Be.push_back(this_total_CRONUS);
    
    set_Braucher_parameters();
    scale_F_values(this_P, nuclides_for_scaling);
    P_mu_braucher.push_back( P0_10Be*(F_10Be[1]*exp(-z/Gamma[1]) + F_10Be[2]*exp(-z/Gamma[2]) +
                             F_10Be[3]*exp(-z/Gamma[3])));
    P_total_braucher.push_back( P0_10Be*(F_10Be[0]*exp(-z/Gamma[0])+F_10Be[1]*exp(-z/Gamma[1]) + 
                             F_10Be[2]*exp(-z/Gamma[2]) + F_10Be[3]*exp(-z/Gamma[3])));
                             
    set_Schaller_parameters();
    scale_F_values(this_P, nuclides_for_scaling);
    P_mu_schaller.push_back( P0_10Be*(F_10Be[1]*exp(-z/Gamma[1]) + F_10Be[2]*exp(-z/Gamma[2]) +
                             F_10Be[3]*exp(-z/Gamma[3])));
    P_total_schaller.push_back( P0_10Be*(F_10Be[0]*exp(-z/Gamma[0])+F_10Be[1]*exp(-z/Gamma[1]) + 
                             F_10Be[2]*exp(-z/Gamma[2]) + F_10Be[3]*exp(-z/Gamma[3])));
                                                          
    set_Granger_parameters();
    scale_F_values(this_P, nuclides_for_scaling);
    P_mu_granger.push_back( P0_10Be*(F_10Be[1]*exp(-z/Gamma[1]) + F_10Be[2]*exp(-z/Gamma[2]) +
                             F_10Be[3]*exp(-z/Gamma[3])));
    P_total_granger.push_back( P0_10Be*(F_10Be[0]*exp(-z/Gamma[0])+F_10Be[1]*exp(-z/Gamma[1]) + 
                             F_10Be[2]*exp(-z/Gamma[2]) + F_10Be[3]*exp(-z/Gamma[3])));
                                                          
    set_newCRONUS_parameters();
    scale_F_values(this_P, nuclides_for_scaling);
    P_mu_newCRONUS.push_back( P0_10Be*(F_10Be[1]*exp(-z/Gamma[1]) + F_10Be[2]*exp(-z/Gamma[2]) +
                             F_10Be[3]*exp(-z/Gamma[3])));
    P_total_newCRONUS.push_back( P0_10Be*(F_10Be[0]*exp(-z/Gamma[0])+F_10Be[1]*exp(-z/Gamma[1]) + 
                             F_10Be[2]*exp(-z/Gamma[2]) + F_10Be[3]*exp(-z/Gamma[3])));
  }
  
  for(int i = 0; i<n_z; i++)
  {
    prod_file_out  << z_mu[i] << "," << P_mu_z_10Be[i] << "," << P_mu_schaller[i] << ","
                   << P_mu_granger[i] << ","  << P_mu_braucher[i] << "," 
                   << P_mu_newCRONUS[i] << "," << P_total_z_10Be[i] << "," 
                   << P_total_schaller[i] << ","
                   << P_total_granger[i] << ","  << P_total_braucher[i] << "," 
                   << P_total_newCRONUS[i] << endl;
  }
  
  prod_file_out.close();

}

#endif
