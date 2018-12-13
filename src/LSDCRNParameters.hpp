//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDCRNParameters.hpp
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
#include <iostream>
#include <vector>
#include <map>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

// Sorting compiling problems with MSVC
#ifdef _WIN32
#ifndef M_PI
extern double M_PI;
#endif
#endif


#ifndef LSDCRNParameters_H
#define LSDCRNParameters_H

/// @brief This class contains parameters used in cosmogenic nuclide calculations
/// It sits seperately from the particle object since it applies to an
/// entire environment and not just an individual particle. 
/// Seperating the object in this way reduces memory redundancy 
class LSDCRNParameters
{
  public:
  
  /// @brief The default constructor. It is the only possible constructor
  LSDCRNParameters()         { create(); }
  
  /// This is a friend class so that it can be called from the particle 
  friend class LSDCRNParticle;

  /// @brief function for loading parameters that allow pressure calculation
  /// from elevation
  /// @author SMM
  /// @date 02/12/2014
  void load_parameters_for_atmospheric_scaling(string path_to_params);
  
  /// @brief This function sets a numer of parameters that are used
  /// to replicate the CRONUS calculator.
  /// @details the original parameters are derived from the 
  /// make_al_be_consts_v22
  /// Written by Greg Balco -- Berkeley Geochronology Center
  ///  balcs@u.washington.edu -- balcs@bgc.org
  ///  February, 2008
  ///  Part of the CRONUS-Earth online calculators: 
  ///     http://hess.ess.washington.edu/math
  /// @author SMM
  /// @date 06/12/2014
  void set_CRONUS_data_maps();

  /// @param This function returns the stone production prescalings
  ///  for 10Be and 26Al
  /// @return Prefs a vector<double> that holds:
  ///  Prefs[0] = stone prescaling of 10Be preduction 
  ///  Prefs[1] = stone prescaling of 26Be preduction 
  /// @author SMM
  /// @date 14/10/2014
  vector<double> get_Stone_Pref();

  /// @brief This function wraps the CRONUS muon production function
  ///  It returns a vector with elements
  ///  Muon_production[0] = 10Be fast
  ///  Muon_production[1] = 26Al fast
  ///  Muon_production[2] = 10Be neg
  ///  Muon_production[3] = 26Al neg
  /// @param z depth below the surface z (g/cm2)
  /// @param h atmospheric pressure (hPa)
  /// @return a four element vector containing:
  ///   Muon_production[0] = 10Be fast
  ///   Muon_production[1] = 26Al fast
  ///   Muon_production[2] = 10Be neg
  ///   Muon_production[3] = 26Al neg
  /// @author SMM
  /// @date 14/10/2014
  vector<double> calculate_muon_production_CRONUS(double z, double h);


  /// @brief calculates the production rate of Al-26 or Be-10 by muons
  /// @detail This uses the scheme in Heisinger and others (2002, 2 papers). The
  ///  vertically traveling muon flux is scaled to the site elevation using
  ///  energy-dependent attenuation lengths from Boezio et al. (2000). See the 
  ///  hard-copy documentation for detailed citations and a full discussion of
  ///  the calculation. 
  ///  Note that some constants are internal to the function. The only ones that
  ///  get passed from upstream are the ones that a) are nuclide-specific, or b) 
  ///  actually have quoted uncertainties in Heisinger's papers. 
  ///  The fraction of muons that are negative is internal; so is the
  ///  energy-dependence exponent alpha.
  ///  Original Written by Greg Balco -- UW Cosmogenic Nuclide Lab
  ///  balcs@u.washington.edu
  ///  March, 2006
  ///  Part of the CRONUS-Earth online calculators: 
  ///      http://hess.ess.washington.edu/math
  /// @param z depth below the surface z (g/cm2)
  /// @param h atmospheric pressure (hPa)
  /// @author SMM
  /// @date 06/12/2014
  void P_mu_total(double z,double h);
 
  /// @brief A wrapper for the P_mu_total function that replaces 
  ///  total production for 10Be and 26Al due to muons
  /// @param z depth below the surface z (g/cm2)
  /// @param h atmospheric pressure (hPa)
  /// @param Be10_total_mu the total muon production for 10Be. This
  ///  is replaced by the function. 
  /// @param 26Al_total_mu the total muon production for 26Al. This
  ///  is replaced by the function.  
  /// @author SMM
  /// @date 15/12/2014
  void P_mu_total_return_nuclides(double z,double h, double& Be10_total_mu,
                                     double& Al26_total_mu);
                                     
                                     
                                     
 
  /// @brief this subfunction returns the stopping rate of vertically traveling muons
  /// as a function of depth z at sea level and high latitude.
  /// @detail Modified from Greg Balco's CRONUS calculator
  /// @param z is the depth below the surface in g/cm^2
  /// @return Rv0 the muon stopping rate
  /// @author SMM
  /// @date 06/12/2014
  double Rv0(double z); 

  /// @brief this subfunction returns the effective atmospheric attenuation length for
  /// muons of range Z
  /// @detail Original by Greg Balco as part of the CRONUS calculator
  /// @param z is the depth in g/cm^2
  /// @return effective atmospheric attenuation length in g/cm^2
  /// @author SMM
  /// @date 06/12/2014
  double LZ(double z);

  /// @brief subroutine for integrating the muon flux
  /// @detail uses simplsons rule, keeps refining nodes until a tolerance is 
  ///  reached. 
  /// @param z the depth of the sample in g/cm^2
  /// @param H the atmospheric depth in g/cm^2
  /// @param the tolerance; successive refined meshes must exceed this tolerance
  ///  in order for the iteration to be sucessfull
  /// @author SMM
  /// @date 07/12/2014
  double integrate_muon_flux(double z, double H, double tolerance);
 
  // functions for altering the parameter values
  
  /// @brief This resets the F, Gamma and P0 values so that they conform to 
  /// Granger and Smith 2000 scaling. Adopted from from Vermeesh 2007
  /// @author SMM
  /// @date 01/01/2010
  void set_Granger_parameters();
  
  /// @brief This resets the F, Gamma and P0 values so that they conform to 
  /// Schaller (2009) scaling. Adopted from from Vermeesh 2007
  /// @author SMM
  /// @date 01/01/2010
  void set_Schaller_parameters();

  /// @brief This resets the F, Gamma and P0 values so that they conform to 
  /// Braucher et al (2009) scaling. Adopted from from Vermeesh 2007, 
  /// @detail From version 2.0 of cosmocalc
  /// @author SMM
  /// @date 27/01/2015
  void set_Braucher_parameters();
  
  /// @brief This resets the F, Gamma and P0 values 
  /// For 10Be, these correspond to new production curves provided by Shasta Marerro
  //  For the rest they conform to 
  /// Braucher et al (2009) scaling. Adopted from from Vermeesh 2007, 
  /// @detail From version 2.0 of cosmocalc
  /// @author SMM
  /// @date 28/01/2016
  void set_newCRONUS_parameters();
  
    
  /// @brief this resets the production and decay coefficients of 10Be and 26Al
  ///  to mimic the parameters for stone scaling in CRONUS calculator
  /// @detail IMPORTANT the F and Gamma numbers are not changed so you will 
  ///  need to set granger or schaller parameters beforehand. 
  ///  ALSO this is no longer necessary due to changes in Vermeesch's reported
  ///  values which now correspond to the CRONUS values
  /// @author SMM
  /// @date 17/12/2014
  void set_CRONUS_stone_parameters();

  /// @brief This sets the F values to use neutron only production
  /// @details F0 == 1, all other F values == 0
  /// @author SMM
  /// @date 14/07/2014	
  void set_Neutron_only_parameters();

  /// @brief This function resets the P0 using the error from the CRONUS 
  ///  calculator. It allows one to test the uncertainty in the
  ///  calculated erosion rates. 
  /// @detail this version adds to the production. 
  /// @return a vector continaing the change in the 10Be and 26Al production rates
  ///  this is used in the gaussian error propigation
  /// @author SMM
  /// @date 03/05/2014
  vector<double> set_P0_CRONUS_uncertainty_plus();

  /// @brief This function resets the P0 using the error from the CRONUS 
  ///  calculator. It allows one to test the uncertainty in the
  ///  calculated erosion rates. 
  /// @detail this version adds to the production
  /// @return a vector continaing the change in the 10Be and 26Al production rates
  ///  this is used in the gaussian error propigation
  /// @author SMM
  /// @date 03/05/2014
  vector<double> set_P0_CRONUS_uncertainty_minus();

  /// @breif this gets difference in the fraction of spallation for different
  /// paris of muon production schemes
  /// @detail key for pairs:
  ///  0 Braucher-Schaller
  ///  1 Braucher-Granger
  ///  2 Braucher-Neutron Only
  ///  3 Schaller-Granger
  ///  4 Schaller-Neutron Only
  ///  5 Granger-Neutron Only
  /// @param pair key for the pair (see key above)
  /// @return vector containing the error from the pair
  /// @author SMM
  /// @date 03/02/2015
  vector<double> get_uncertainty_scaling_pair(int pair);


  /// @brief This sets the internal scaling for the particle. It includes
  ///  topographic shielding, snow shielding and scaling from latitude and
  ///  other factors
  /// @param scaling the lat-magnetic scaling 
  /// @param topo_shield the topograpgic shielding
  /// @param snow_shield shielding from snow
  /// @author SMM
  /// @date 17/12/2014
  void set_scaling(double scaling, double topo_shield, double snow_shield);

  /// @brief This sets the internal scaling for the particle. It includes
  ///  topographic shielding, snow shielding and scaling from latitude and
  ///  other factors. It is for use with neutron only calculations
  /// @param scaling the lat-magnetic scaling 
  /// @param topo_shield the topograpgic shielding
  /// @param snow_shield shielding from snow
  /// @author SMM
  /// @date 17/12/2014
  void set_neutron_scaling(double scaling, double topo_shield, double snow_shield);




  /// @brief This gets the Lifton Scaling
  /// Modified from Greg Balco's code:
  ///  http://hess.ess.washington.edu/math
  /// @param latitude in decimal degrees
  /// @param pressure in hPa
  /// @param fsp is the fraction (between 0 and 1) of production at sea level
  ///  and high latitude due to spallation (as opposed to muons).
  ///  This argument is optional and defaults to 0.978, which is the value
  ///  used by Stone (2000) for Be-10. The corresponding value for Al-26
  ///  is 0.974. Note that using 0.844 for Be-10 and 0.826 for Al-26 will
  ///  closely reproduce the Lal, 1991 scaling factors as long as the standard
  ///  atmosphere is used to convert sample elevation to atmospheric pressure.
  ///  Also note that this function will yield the scaling factor for spallation
  ///  only when fsp=1, and that for muons only when fsp=0.
  ///  @return the scaling factor
  /// @author SMM
  /// @date 5/12/2014
  double stone2000sp(double lat,double P, double Fsp);

  /// @brief this function takes a single scaling factor for
  /// elevation scaling, self shielding, snow shielding,
  /// and latitude scaling and produces scaling factors
  /// for each production mechamism.
  /// the scaling follows the approach of vermeesch 2008
  /// it uses a 'virtual' shielding depth to calculate
  /// the updated scaling factors
  /// @param single_scaling a lumped scaling factor	
  /// @author SMM
  /// @date 01/01/2010
  void scale_F_values(double single_scaling);

  /// @brief this function takes a single scaling factor for
  /// elevation scaling, self shielding, snow shielding,
  /// and latitude scaling and produces scaling factors
  /// for each production mechamism.
  /// the scaling follows the approach of vermeesch 2008
  /// it uses a 'virtual' shielding depth to calculate
  /// the updated scaling factors
  /// @param single_scaling a lumped scaling factor
  /// @param nuclides_for_scaling this is a vector of bool telling the code 
  ///  which nuclides to calculate. The values are:
  ///  nuclides_for_scaling[0] = true: calculate 10Be
  ///  nuclides_for_scaling[1] = true: calculate 26Al
  ///  nuclides_for_scaling[2] = true: calculate 36Cl
  ///  nuclides_for_scaling[3] = true: calculate 14C
  /// @author SMM
  /// @date 01/02/2015
  void scale_F_values(double single_scaling, vector<bool> nuclides_for_scaling);

  /// @brief Prints F calues to screen for bug checking
  /// @param nuclides_for_scaling this is a vector of bool telling the code 
  ///  which nuclides to calculate. The values are:
  ///  nuclides_for_scaling[0] = true: calculate 10Be
  ///  nuclides_for_scaling[1] = true: calculate 26Al
  ///  nuclides_for_scaling[2] = true: calculate 36Cl
  ///  nuclides_for_scaling[3] = true: calculate 14C
  /// @author SMM
  /// @date 23/02/2015
  void print_F_values_to_screen(vector<bool> nuclides_for_scaling);

  /// @brief Prints parameters to screen for bug checking
  /// @param nuclides_for_scaling this is a vector of bool telling the code 
  ///  which nuclides to calculate. The values are:
  ///  nuclides_for_scaling[0] = true: calculate 10Be
  ///  nuclides_for_scaling[1] = true: calculate 26Al
  ///  nuclides_for_scaling[2] = true: calculate 36Cl
  ///  nuclides_for_scaling[3] = true: calculate 14C
  /// @author SMM
  /// @date 23/02/2015
  void print_parameters_to_screen(vector<bool> nuclides_for_scaling);

  /// @brief this changes the 10Be decay. It is here because
  /// 10Be decay rates reported in the literature have changed
  /// @author SMM
  /// @date 01/01/2010
  void update_10Be_decay(double new_decay)	{ lambda_10Be = new_decay; }
  
  /// @brief this changes the 10Be P0 value. It is here because
  /// 10Be decay rates reported in the literature have changed
  /// @author SMM
  /// @date 01/01/2010	
  void update_10Be_P0(double new_P0)			{ P0_10Be = new_P0; }

  /// @brief This calcualtes the atmospheric pressure given latidude, longitude
  /// and elevation
  /// @details Looks up surface pressure and 1000 mb temp from NCEP reanalysis
  /// and calculates site atmospheric pressures using these as inputs to the
  /// standard atmosphere equation. 
  /// Also: This function is OK but not great for Antarctica.
  /// Use antatm.m instead. 
  /// Remember: it is always better to estimate the average pressure at your 
  /// site using a pressure-altitude relation obtained from nearby station
  /// data.
  ///
  /// Original m code Written by Greg Balco -- UW Cosmogenic Nuclide Lab
  /// @param site_lat latitude (DD). Southern hemisphere is negative.
  /// @param site_lon longitude (DD). Western hemisphere is negative.
  ///       Tries to deal with 0-360 longitudes gracefully.
  /// @param site_elv elevation (m).
  /// @return site pressure in hPa.
  /// @author SMM
  /// @date 04/12/2014
  double NCEPatm_2(double site_lat, double site_lon, double site_elev);
  
  /// @brief This gets the attenuation depth in g/cm^2
  ///  You tell it if you want the CRONUS values
  /// @param use_CRONUS a bool that allows you to select the CRONUS
  ///  value
  /// @return the spallation attenuation length in g/cm^2
  /// @author SMM
  /// @date 14/12/2014
  double get_spallation_attenuation_length(bool use_CRONUS);
  
  /// @brief This gets the decay coefficnets for both 10Be and 26Al
  ///  You tell it if you want the CRONUS values
  /// @param use_CRONUS a bool that allows you to select the CRONUS
  ///  value
  /// @return the decay coefficients decay_coeff in vector<double>
  ///   decay_coeff[0] = 10Be decay in yr^-1
  ///   decay_coeff[1] = 26Al decay in yr^-1
  /// @author SMM
  /// @date 14/12/2014
  vector<double> get_decay_coefficients(bool use_CRONUS);
  
  /// @brief This function calculates the CRONUS version 
  ///  of the mu production vector
  /// @param pressure takes atmospheric pressure in HPa
  /// @param effective_depth and the effective depth in g/cm^2
  /// @param z_mu a vector<double> that is replaced in the function with the 
  ///  muon production depths
  /// @param P_mu_z_10Be a vector<double> that is replaced in the function with 
  ///  muon production at the depths in z_mu. This one for 10Be.
  /// @param P_mu_z_10Be a vector<double> that is replaced in the function with 
  ///  muon production at the depths in z_mu. This one for 26Al.  
  /// @author SMM
  /// @date 15/12/2014
  void get_CRONUS_P_mu_vectors(double pressure, double sample_effective_depth,
                               vector<double>& z_mu, vector<double>& P_mu_z_10Be,
                               vector<double>& P_mu_z_26Al);
  
  /// @brief subroutine for integrating the muon flux for a given erosion rate
  /// @detail uses simplsons rule, keeps refining nodes until a tolerance is 
  ///  reached. 
  /// @detail Note this is a little different from the CRONUS calculator since
  ///  CRONUS uses trapezoid rule and this uses Simpson's rule so this function
  ///  is a little bit more accurate.
  /// @param E the target erosion rate
  /// @param z_mu a vector of effective depths over which to integrate. 
  ///  comes from the get_CRONUS_P_mu_vectors
  /// @param P_mu_z_10Be a vector<double> that is muon production at the depths 
  ///   in z_mu. This one for 10Be.
  /// @param P_mu_z_10Be a vector<double> that is muon production at the depths 
  ///   in z_mu. This one for 26Al.
  /// @param Be10_mu_N atoms producted of 10Be. Is replaced in the function.
  /// @param Al26_mu_N atoms produced of 26Al. Is replaced in the function
  /// @author SMM
  /// @date 15/12/2014
  void integrate_muon_flux_for_erosion(double E, vector<double> z_mu,
                           vector<double> P_mu_10Be, vector<double> P_mu_26Al,
                           double& Be10_mu_N, double& Al26_mu_N);
  
  /// @brief this function calculates the total atoms from spallation
  /// for a given erosion rate that replicates the CRONUS calculator
  /// @param E the erosion rate in g/cm^2/yr
  /// @param thick_SF the thickness scaling factor (between 0 and 1)
  /// @param P_sp_10Be The production rate of 10Be. Includes shielding corrections
  /// @param P_sp_26Al The production rate of 26Al. Includes shielding corrections
  /// @param Be10_sp_N The number of atoms from spallation of Be10
  ///  replaced within this function
  /// @param Al26_sp_N The number of atoms from spallation of Al26
  ///  replaced within this function
  /// @author SMM
  /// @date 15/12/2014
  void integrate_nonTD_spallation_flux_for_erosion(double E, 
                                   double thick_SF,
                                   double P_sp_10Be, double P_sp_26Al, 
                                   double& Be10_sp_N,double& Al26_sp_N);
  
  /// @brief THis gets parameters for uncertainty analysis for muogenic 
  ///  production. It is part of the CRONUS calculator replication
  /// @param pressure the atmospheric pressure in HPa
  /// @return a vector<double> with the uncertanty parameters. These are:
  ///  uncert_params[0]=delPfast_10;
  ///  uncert_params[1]=delPfast_26;
  ///  uncert_params[2]=delPneg_10;
  ///  uncert_params[3]=delPneg_26;
  ///  uncert_params[4]=delPmu0_10;
  ///  uncert_params[5]=delPmu0_26;
  /// @author SMM
  /// @date 15/12/2014
  vector<double> CRONUS_get_muon_uncertainty_params(double pressure);                         
  
  /// @brief Gets the uncertanty of spallation and muon production as
  ///  a ratio of the total production
  /// @detail used in error propagation for the CRONUS calculator
  /// @param scaling_name a string containing the scaling scheme you want
  ///  choices are:
  ///  St  for Stone scaling
  ///  Li  for Lifton scaling
  ///  De  for Deselits scaling
  ///  Du  for Dunai scaling
  ///  Lm  for Lal magnetic scaling
  /// @author SMM
  /// @date 18/12/2014
  vector<double> CRONUS_get_uncert_production_ratios(string scaling_name);

  /// @brief this Prints the parameters to a file that is used for 
  ///  testing and reconstruction of results
  /// @param fname the filename with full file path
  /// @param muon_scaling the scaling scheme, either Granger, Schaller or Braucher
  ///  the default is Braucher
  /// @author SMM
  /// @date 15/03/2015
  void print_parameters_to_file(string fname, string muon_scaling);

  /// @brief this Prints the prodcution rates for muons for various schemes to 
  /// a csv file
  /// @detail The parameter file contains a column for effective depth, and then
  ///  columns for the production of muons from different production schemes
  /// @param fname the filename with full file path
  /// @param path_to_atmospheric_data the path to the folder containing atmospheric data
  /// @author SMM
  /// @date 18/02/2016
  void Print_10Beproduction_csv(string filename, string path_to_atmospheric_data);
  
  private:
  /// @brief This is called by the default constructor. 
  /// It is the only possible constructor
  void create();
  
  /// the version number of this CRNParameters object
  string version;

  /// The attenudation lengths in g/cm^2. Each element in the array
  /// refers to a different production mechanism
  double Gamma[4];			// attenuation legths in g/cm^2

  /// F values allocating production to spallation and muon production for 10Be
  double F_10Be[4];
  
  /// F values allocating production to spallation and muon production for 14C
  double F_14C[4];
  
  /// F values allocating production to spallation and muon production for 26Al
  double F_26Al[4];
  
  /// F values allocating production to spallation and muon production for 36Cl
  double F_36Cl[4];
  
  /// F values allocating production to spallation and muon production for 21Ne
  double F_21Ne[4];
  
  /// F values allocating production to spallation and muon production for 3He
  double F_3He[4];

  /// topographic shielding
  /// other shielding and scaling calucalted using the scale_F_values function
  double S_t;
  
  /// other shielding and scaling calucalted using the scale_F_values function
  /// for neutrons only
  double neutron_S_t;

  /// decay rate for 10Be in yr-1
  double lambda_10Be;
  
  /// decay rate for 26Al in yr-1	
  double lambda_26Al;
  
  /// decay rate for 14C in yr-1	
  double lambda_14C;	
  
  /// decay rate for 36Cl in yr-1	
  double lambda_36Cl;		

  /// production rate for 10Be in a/g/yr
  double P0_10Be;	
  
  /// production rate for 26Al in a/g/yr		
  double P0_26Al;		
  
  /// production rate for 14C in a/g/yr
  double P0_14C;			
  
  /// production rate for 36Cl in a/g/yr
  double P0_36Cl;		
  
  /// production rate for 21Ne in a/g/yr
  double P0_21Ne;			
  
  /// production rate for 3He in a/g/yr
  double P0_3He;			
  
  /// This is a data map used for storing CRONUS calcluator parameters
  map<string,double> CRONUS_data_map;
  
  /// This is a data map used for storing CRONUS muon parameters
  map<string,double> CRONUS_muon_data;
  
  /// levels: the levels for the atmospheric scaling of pressure
  vector<double> levels;
  
  /// This is an index for the latitudes for atmospheric scaling
  vector<double> NCEPlat;
  
  /// This is an index for the longitudes for atmospheric scaling
  vector<double> NCEPlon;
  
  /// This is an array holding sea level perssures
  Array2D<double> meanslp;
  
  /// This is an array healing mean temperatures
  Array2D<double> meant1000;
  
  /// This is a vector of arrays holding something called gp_hgt;
  vector< Array2D<double> > gp_hgt;
  
  /// This is a vector of arrays holding something called gp_hgt;
  vector< Array2D<double> > gm_hgt;  
  
  
};

#endif
