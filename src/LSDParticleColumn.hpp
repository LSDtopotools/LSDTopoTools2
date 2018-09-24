//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDParticleColumn
// Land Surface Dynamics Particle Column
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for tracing particles and retaining geochemical information
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
#include <list>
#include "LSDCRNParameters.hpp"
#include "LSDParticle.hpp"
using namespace std;

#ifndef LSDParticleColumn_H
#define LSDParticleColumn_H


/// @brief This is a class for a particle that can be tracked through simulations
/// and retains data about position and chemical content
class LSDParticleColumn
{
	public:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Constructors
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    /// @brief The default constructor. It is the only possible constructor
	  LSDParticleColumn()			{ create(); }
	
	  /// @brief Constructor where all data member are assigned. Used for copy 
	  /// constructors
    LSDParticleColumn(int tRow, int tCol, int tNodeIndex, 
        double tSoilDensity, double tRockDensity,   
        float tSoilThickness, float tDataResolution, bool tUseDenstyProfile, 
        vector<double> tDensityDepths, vector<double> tDensityDensities, 
        list<LSDCRNParticle> tCRNParticleList)
        { create(tRow, tCol, tNodeIndex, tSoilDensity, tRockDensity,   
                 tSoilThickness,tDataResolution, tUseDenstyProfile, 
                 tDensityDepths, tDensityDensities, tCRNParticleList); }
 
    /// @brief const reference operator
    LSDParticleColumn(const LSDParticleColumn& tPC)
    	{ create(tPC.getRow(), tPC.getCol(), tPC.getNodeIndex(), 
               tPC.getSoilDensity(),tPC.getRockDensity(),
               tPC.getSoilThickness(),tPC.getDataResolution(),
               tPC.getUseDenstyProfile(),tPC.getDensityDepths(),
               tPC.getDensityDensities(), tPC.getCRNParticleList()); }
    
    /// @brief const reference operator
    LSDParticleColumn(LSDParticleColumn& tPC)
    	{ create(tPC.getRow(), tPC.getCol(), tPC.getNodeIndex(), 
               tPC.getSoilDensity(),tPC.getRockDensity(),
               tPC.getSoilThickness(),tPC.getDataResolution(),
               tPC.getUseDenstyProfile(),tPC.getDensityDepths(),
               tPC.getDensityDensities(), tPC.getCRNParticleList()); }

    /// @brief const copy constructor
    LSDParticleColumn& operator=(const LSDParticleColumn& tPC); 

    /// @brief copy constructor
    LSDParticleColumn& operator=(LSDParticleColumn& tPC); 
         
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // GETTER FUNCTIONS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-         
    /// @brief function for the row
    /// @return Row the row
    int getRow() const { return Row;}
      
    /// @brief function for the Col
    /// @return Col the column of the particle column (links to LSDRaster)
    int getCol() const { return Col;}
           
    /// @brief function for the NodeIndex
    /// @return NodeIndex of the particle column (links to LSDRaster)
    int getNodeIndex() const { return NodeIndex;}
    
    /// @brief function for the SoilDensity
    /// @return SoilDensity of the particle column
    double  getSoilDensity() const { return SoilDensity;}    

    /// @brief function for the RockDensity
    /// @return RockDensity of the particle column
    double  getRockDensity() const { return RockDensity;}        

    /// @brief function for the SoilThickness
    /// @return SoilThickness of the particle column
    float getSoilThickness() const { return SoilThickness;}   
    
    /// @brief function for the DataResolution
    /// @return DataResolution of the particle column (links to LSDRaster)
    float getDataResolution() const { return DataResolution;}       

    /// @brief function for the UseDenstyProfile;
    /// @return UseDenstyProfile tells object to use density profile
    bool getUseDenstyProfile() const { return UseDenstyProfile;} 

    /// @brief function for the DensityDepths;
    /// @return DensityDepths depth density measurments in the profile
    vector<double> getDensityDepths() const { return DensityDepths;} 

    /// @brief function for the DensityDensities;
    /// @return DensityDensities density measurments in the profile (kg/m^3)
    vector<double> getDensityDensities() const { return DensityDensities;} 

    /// @brief function for the CRNParticleList
    /// @return CRNParticleList the underlying particle list
     list<LSDCRNParticle> getCRNParticleList() const { return CRNParticleList;} 


    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Some setter functions
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    /// @brief set the rock density 
    void set_RockDensity( double new_rhoR) {RockDensity = new_rhoR;}
    
    /// @brief set the soil density
    void set_SoilDensity( double new_rhoS) {SoilDensity = new_rhoS;}
    
    /// @brief set the soil thickness
    void set_SoilThickness( double new_h)  {SoilThickness = new_h; }

    /// @brief this sets the row and column data members
    void set_Row_and_Col(int thisRow, int thisCol) { Row = thisRow; Col = thisCol; }

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // Computations 
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-      
    /// @brief This function initiates a list with particles evenly spaced
    /// @param start_type the starting type of the particle
    /// @param startxLoc the starting x location of the particle
    /// @param startxLoc the starting y location of the particle
    /// @param start_depth the starting depth of the particle in metres
    /// @param particle_spacing the spacing in metres between particles (in depth)
    /// @param zeta the surface elevation (m)
    /// @param eff_eros_rate the erosion rate in g/cm^2/yr
    /// @param CRN_param a CRNParameters object
    /// @author SMM
    /// @date 25/07/2014
    void initiate_SS_cosmo_column_3CRN(int start_type, 
          double startxLoc, double startyLoc,
		      double start_depth, double particle_spacing, 
		      double zeta, double eff_eros_rate,
		      LSDCRNParameters& CRN_param);

    /// @brief This inserts particles into the column
    /// @param start_type the starting type of the particle
    /// @param startxLoc the starting x location of the particle
    /// @param startxLoc the starting y location of the particle
    /// @param start_depth the starting depth of the particle in metres
    /// @param zeta the elevation of the surface
    void insert_particle_into_column(int start_type, 
       double startxLoc, double startyLoc,
	     double start_depth, double zeta);

    /// @brief This function updates the CRN concentrations, updates the
    /// zeta values in the face of uplift, and gives partial exposure for
    /// particles near the surface so these particles should represent the 
    /// concentration at the time of erosion. 
    /// @detail The returned particle column contains the eroded particles. 
    /// The function uses rock only since at the moment there is not time
    /// to include soil layers and all the computation required to sort out
    /// the density expansion
    /// @param dt the time step
    /// @param uplift_rate the rate of uplift in m/yr
    /// @param start_type the starting type of the particle
    /// @param start_depth the starting depth of the particle in metres
    /// @param startxLoc the starting x location of the particle
    /// @param startxLoc the starting y location of the particle
    /// @param zeta_old the elevation of the previous timestep
    /// @param zeta_new the elevation at this timestep
    /// @param particle spacing the distance between particles
    /// @param CRN_param a LSDCRNParameters object 
    /// @return a particle column that contains eroded particles. These particles
    /// have been eroded in the last timestep but they are only exposed 
    /// during the time they are in the ground. This is to replicate collection
    /// of particles that are in streams and have only just emerged from the ground
    /// @author SMM
    /// @date 25/07/2014
    LSDParticleColumn update_CRN_list_rock_only_eros_limit_3CRN(	double dt, 
                     double uplift_rate, int start_type, double start_depth,
	                   double startxLoc, double startyLoc,
	                   double zeta_old,double zeta_new,
                     double particle_spacing, LSDCRNParameters& CRN_param);

 
    /// @brief the purpose of this object is to collect particles into one
    /// column that can be used, for example, to explore the CRN concetrations
    /// of particles being eroded from a landscape
    /// @param ColList_vec a vector of particle columns
    /// @author SMM
    /// @date 25/07/2014
    void collect_particles(vector<LSDParticleColumn>& ColList_vec);

    /// @brief this function calculates the apparent erosion from the top
    /// particle in the column using neutron only assumption, and only rock density
    /// @param CRN_param and LSDCRNParameters object
    /// @return a vector with the apprent erosion rates in m/yr
    /// [0] from 10Be
    /// [1] from 14C
    /// [2] from 21Ne
    /// @author SMM
    /// @date 25/7/2014 
    vector<double> calculate_app_erosion_3CRN_neutron_rock_only(LSDCRNParameters& CRN_param ); 
 
  /// @brief this prints particles to screen for bug checking
  /// @param CRNParams the CRN parameter file (for apparent erosion
  /// @author SMM
  /// @date 02/08/2014
  void print_particle_properties_to_screen(LSDCRNParameters& CRNParam);  
  
	protected:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // DATA MEMBERS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    /// The row of the particle column. 
    /// This allows the obect to be linked to an LSDRaster	
    int Row;

    /// The row of the particle column. 
    /// This allows the obect to be linked to an LSDRaster	    
    int Col;

    /// The row of the particle column. 
    /// This allows the obect to be linked to an LSDFlowInfo object	    
    int NodeIndex;
    
    /// the soil density in kg/m^3
    double SoilDensity;
    
    /// the rock density in kh/m^3
    double RockDensity;
    
    /// the soil thickness. It is a float since LSDRasters are floats
    float SoilThickness;

    /// the cellsize of the raster if the column is linked to a raster.
    /// It is a float since LSDRasters are floats
    /// used to track if the particle is within the cell
    float DataResolution;
    
    /// This tells the model if it needs to use the density profile
    bool UseDenstyProfile;
    
    /// This is the density profile. The depths are the points where
    /// density has been sampled
    vector<double> DensityDepths;
    
    /// This is the density profile. The densities (in kg/m^3) are the density 
    /// measurments at depths set by DenstiyDepths    
    vector<double> DensityDensities;
    
    /// This is the list for holding the particles
    list<LSDCRNParticle> CRNParticleList;
    

    
    
  private:

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //
    // CREATE FUNCTIONS
    //
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    /// @brief default create function
    void create();		
    
    /// @brief create function using all data members
    void create(int tRow, int tCol, int tNodeIndex, 
        double tSoilDensity, double tRockDensity,   
        float tSoilThickness, float tDataResolution, bool tUseDenstyProfile, 
        vector<double> tDensityDepths, vector<double> tDensityDensities, 
        list<LSDCRNParticle> tCRNParticleList);
};

#endif
