//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDParticle
// Land Surface Dynamics Particle
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
#include "LSDCRNParameters.hpp"
using namespace std;

#ifndef LSDParticle_H
#define LSDParticle_H

// empty class definition so that we can use friend functions
class LSDCRNParameters;

/// @brief This is a class for a particle that can be tracked through simulations
/// and retains data about position and chemical content
class LSDParticle
{
 public:
 
    /// @brief default constructor. Makes a particle with type 0, age 0, and 
    /// no data for other parameters
    LSDParticle()					{ create(); }
    
    /// @brief Constructor. Similar to default, but assignes type. 
    /// @param StartType the type of the particle
    LSDParticle( int StartType)			{ create(StartType); }
    
    /// @brief Constructor. Similar to default, but assignes type and age. 
    /// @param StartType the type of the particle    
    /// @param StartAge the starting age of the particle
    LSDParticle( int StartType, double StartAge)	{ create(StartType, StartAge); }
    
    /// @brief Constructor. Assignes all data members except yLoc
    /// @param StartType the type of the particle 
    /// @param StartCI the starting cell index   
    /// @param StartAge the starting age of the particle   
    /// @param StartOSLage the starting OSL age
    /// @param StartxLoc the starting x location
    /// @param StartdLoc the starting depth
    LSDParticle( int StartType, int StartCI, double StartAge, double StartOSLage, 
                double StartxLoc, double StartdLoc)
                { create(StartType, StartCI, StartAge, StartOSLage, StartxLoc, StartdLoc); }

    /// @brief Constructor. Assignes all data members
    /// @param StartType the type of the particle 
    /// @param StartCI the starting cell index   
    /// @param StartAge the starting age of the particle   
    /// @param StartOSLage the starting OSL age
    /// @param StartxLoc the starting x location
    /// @param StartyLoc the starting y location
    /// @param StartdLoc the starting depth
    LSDParticle( int StartType, int StartCI, double StartAge, double StartOSLage, 
                double StartxLoc, double StartyLoc, double StartdLoc)
                { create(StartType, StartCI, StartAge, StartOSLage, StartxLoc, StartyLoc,StartdLoc); }

    /// @brief Constructor. Assignes type, x location and d location. 
    /// other parameters are default
    /// @param StartType the type of the particle 
    /// @param StartxLoc the starting x location
    /// @param StartdLoc the starting depth    						
    LSDParticle( int StartType, double StartxLoc, double StartdLoc)
                { create(StartType, StartxLoc, StartdLoc); }

    /// @brief Get the type
    /// @return Type the type of the particle
    int    getType() const			{ return Type; }
    
    /// @brief Get the CellIndex
    /// @return CellIndex the CellIndex of the particle
    int    getCellIndex() const		{ return CellIndex; }
    
    /// @brief Get the Age
    /// @return Age the Age of the particle. This is the time since it has
    /// entered the soil layer  
    double getAge() const			{ return Age; }
    
    /// @brief Get the OSLage
    /// @return OSLage the Optically stimulated luminesce age of the particle      
    double getOSLage() const		{ return OSLage; }
    
    /// @brief Get the xLoc
    /// @return  xLoc the  x location of the particle         
    double getxLoc() const			{ return xLoc; }
    
    /// @brief Get the yLoc
    /// @return  yLoc the  y location of the particle         
    double getyLoc() const			{ return yLoc; }    
    
    /// @brief Get the dLoc
    /// @return  dLoc the  depth of the particle       
    double getdLoc() const			{ return dLoc; }

    /// @brief const reference operator
    LSDParticle(const LSDParticle& tP)
    	{ create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),
               tP.getxLoc(),tP.getyLoc(),tP.getdLoc()); }
    
    /// @brief const reference operator
    LSDParticle(LSDParticle& tP)
    	{ create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),
               tP.getxLoc(),tP.getyLoc(),tP.getdLoc()); }

    /// @brief copy constructor
    LSDParticle& operator=(const LSDParticle& tP);
    
    /// @brief outstream operator
    std::ostream& operator<< (std::ostream&);

    /// @brief Increase the age of the particle
    /// @param dt the time increment
    void incrementAge(double dt);
    
    /// @brief Change the cell index. Used for referencing where the particle is
    /// this could be changed in the future to an octree structure
    /// @param CI the cell index
    void setCellIndex(int CI);
    
    /// @brief This 'exposes' the particle to light and resets the OSL age
    void OSLexpose();
    
    /// @brief This 'exposes' the particle and resets its soil age to zero
    void SoilAgeExpose();
    
    /// @brief Allows the user to update the horizontal location (dangerous!)    
    void update_xLoc(double new_xLoc);

    /// @brief Allows the user to update the horizontal location (dangerous!)    
    void update_yLoc(double new_yLoc);
	  
    /// @brief this changes a particles horizontal and vertical position. If the
    /// particle reaches the surface it 'reflects' back into the soil layer
    /// @detail IMPORTANT: this increments the OSL age but not the soil age!!
    /// @param dx the distance moved horizontally. 
    /// @param dd the distance moved in depth
    /// @param h the thickness of the soil
    /// @param dt the time over which the particle moves.
    /// @author SMM
    /// @date 01/01/2010 
    void displaceReflect(double dx,double dd,double h,double dt);   
    
    /// @brief this function test to see if the particle is within a hillslope, 
    /// and if not sets the data to no data values
    /// @param lambda the length of the hillslope
    /// @author SMM
    /// @date 01/01/2010   
    int  test_domain(double lambda);

 protected:
 
    /// The type: used to identify, generally, what mineral the particle is
    int Type;
    
    /// The cell index, used for referencing the location of a particle in 
    /// a grid 
    int CellIndex;
    
    /// Age the particle has spent in the soil
    double Age;
    
    /// Optically stimulated luminescence age
    double OSLage;
    
    /// Horizontal location of the particle
    double xLoc;
    
    /// Other coordinate, used in 3D simulations
    double yLoc;
    
    /// Depth of the particle
    double dLoc;

 private:
 
    /// @brief creates a default particle of type 0, 
    /// cell index of -1, age of 0 and OSLage of -9999
    /// @author SMM
    /// @date 01/01/2008
    void create();
    
    /// @brief creates a default particle 
    /// cell index of -1, age of 0 and OSLage of -9999
    /// @param StartType the starting type of the particle
    /// @author SMM
    /// @date 01/01/2008
    void create(int);
    
    /// @brief creates a default particle of 
    /// cell index of -1, and OSLage of -9999
    /// @param StartType the starting type of the particle
    /// @param StartAge the starting age of the particle
    /// @author SMM
    /// @date 01/01/2008    
    void create(int, double);
    
    /// @brief create function where user assigns all data members excepth
    /// yLoc (yLoc is only used in 3D simulations)
    void create(int, int, double, double, double, double);
    
    /// @brief create function where user assigns all data members except 
    void create(int, int, double, double, double, double, double);
    
    /// @brief create function where user assigns type, dLoc and xLoc
    void create(int, double, double);
};

/// @brief CRN tracer particle object, a particle that contains nuclide information
class LSDCRNParticle: public LSDParticle
{
  public:
  /// default constructor
  LSDCRNParticle()          { create(); }

  /// @ brief constructor with starting x location, depth and z location
  LSDCRNParticle(int startType, double startxLoc,double startdLoc,
          double start_effdloc, double startzloc)
              { create(startType,startxLoc,startdLoc,
                start_effdloc,startzloc); }

  /// @ brief constructor with starting location information
  LSDCRNParticle(int startType, double startxLoc, double startyLoc,
               double startdLoc, double start_effdloc, double startzloc)
            { create(startType,startxLoc,startyLoc,startdLoc,
              start_effdloc,startzloc); }
  
  LSDCRNParticle(int startType, double startxLoc,double startzeta_Loc)
            { create(startType,startxLoc,startzeta_Loc); }

  /// @brief a create function for a volume particle
  LSDCRNParticle(int startType, int startGSDType,
              double startxLoc,
              double startdLoc, double start_effdloc,
              double startzLoc, double startMass,
              double startSurfaceArea)
              { create(startType, startGSDType,startxLoc,
                       startdLoc, start_effdloc, startzLoc, startMass,
                       startSurfaceArea); }
              
  /// @brief a create function for a volume particle  with y loc          
  LSDCRNParticle(int startType, int startGSDType, double startxLoc, double startyLoc,
              double startdLoc, double start_effdloc,
              double startzLoc, double startMass,
              double startSurfaceArea)
              { create(startType, startGSDType, startxLoc, startyLoc,
                       startdLoc, start_effdloc, startzLoc, startMass,
                       startSurfaceArea); } 

  /// @brief constructor that includes some CRN information
  LSDCRNParticle(int startType, double startxLoc,double startzeta_Loc,
              double startdLoc, double start_effdLoc,
              double start_C10Be,double start_C14C)
              { create(startType,startxLoc,startzeta_Loc,
                startdLoc, start_effdLoc, start_C10Be,start_C14C); }
                
  
  /// @brief constructor that includes CRN information with many nuclides
  LSDCRNParticle(int startType, double startxLoc,double startzeta_Loc,
              double startdLoc, double start_effdLoc,
              double start_C10Be,double start_C14C, double start_21Ne)
              { create(startType,startxLoc,startzeta_Loc,
                startdLoc, start_effdLoc,
                start_C10Be, start_C14C, start_21Ne); }
  
  /// @brief constructor that retains all data 
  LSDCRNParticle(int startType, int startGSDType, int startCellIndex, double startAge, double startOSLAge,
              double startxLoc,double startyLoc, double startdLoc, double startefdLoc,
              double startzLoc, double start_C10Be, double start_C26Al,
              double start_C36Cl, double start_C14C,
              double start_C21Ne, double start_C3He,
              double start_Cf7Be, double start_Cf10Be,
              double start_Cf210Pb, double start_Cf137Cs, 
              double start_Mass, double start_StartingMass,
              double start_SurfaceArea)
              { create(startType, startGSDType, startCellIndex, startAge, startOSLAge,
                startxLoc, startyLoc,startdLoc, startefdLoc,
                startzLoc, start_C10Be, start_C26Al,
                start_C36Cl, start_C14C,
                start_C21Ne, start_C3He,
                start_Cf7Be, start_Cf10Be,
                start_Cf210Pb, start_Cf137Cs,
                start_Mass, start_StartingMass,start_SurfaceArea); }


  /// @brief The copy constructor
  /// @param tP an LSDCRNParticle object
  /// @author SMM
  /// @date 01/01/2010
  LSDCRNParticle(const LSDCRNParticle& tP)
      { create(tP.getType(), tP.getGSDType(), tP.getCellIndex(), tP.getAge(),tP.getOSLage(),
               tP.getxLoc(),tP.getyLoc(),tP.getdLoc(),
               tP.geteffective_dLoc(),tP.get_zetaLoc(),tP.getConc_10Be(),
               tP.getConc_26Al(), tP.getConc_36Cl(), tP.getConc_14C(),
               tP.getConc_21Ne(), tP.getConc_3He(),
               tP.getConc_f7Be(), tP.getConc_f10Be(),
               tP.getConc_f210Pb(), tP.getConc_f137Cs(),
               tP.getMass(), tP.getStartingMass(),
               tP.getSurfaceArea()); }

  /// @brief The copy constructor
  /// @param tP a constant LSDCRNParticle object
  /// @author SMM
  /// @date 01/01/2010
  LSDCRNParticle(LSDCRNParticle& tP)
      { create(tP.getType(), tP.getGSDType(), tP.getCellIndex(), tP.getAge(),tP.getOSLage(),
               tP.getxLoc(),tP.getyLoc(),tP.getdLoc(),
               tP.geteffective_dLoc(),tP.get_zetaLoc(),tP.getConc_10Be(),
               tP.getConc_26Al(), tP.getConc_36Cl(), tP.getConc_14C(),
               tP.getConc_21Ne(), tP.getConc_3He(),
               tP.getConc_f7Be(), tP.getConc_f10Be(),
               tP.getConc_f210Pb(), tP.getConc_f137Cs(),
               tP.getMass(), tP.getStartingMass(),
               tP.getSurfaceArea());   }

  /// @brief the copy constructor for constant LSDCRNParticle objects  
  LSDCRNParticle& operator=(const LSDCRNParticle& tP);

  /// @brief the copy constructor 
  LSDCRNParticle& operator=(LSDCRNParticle& tP);

  /// @brief this sets 10Be conc
  /// @param the new 10Be conc in atoms per gram
  void setConc_10Be(double new_10BeConc) { Conc_10Be = new_10BeConc; }

  /// @brief this sets 26Al conc
  /// @param the new 26Al conc in atoms per gram
  void setConc_26Al(double new_26AlConc) { Conc_26Al = new_26AlConc; }

  /// @brief Retrieves the concentration of 10Be
  /// @return concentration of 10Be
  double getConc_10Be() const      { return Conc_10Be; }

  /// @brief Retrieves the concentration of 26Al
  /// @return concentration of 26Al
  double getConc_26Al() const      { return Conc_26Al; }

  /// @brief Retrieves the concentration of 36Cl
  /// @return concentration of 36Cl
  double getConc_36Cl() const      { return Conc_36Cl; }

  /// @brief Retrieves the concentration of 14C
  /// @return concentration of 14C
  double getConc_14C() const			{ return Conc_14C; }

  /// @brief Retrieves the concentration of 21Ne
  /// @return concentration of 21Ne
  double getConc_21Ne() const			{ return Conc_21Ne; }

  /// @brief Retrieves the concentration of 3He
  /// @return concentration of 3He
  double getConc_3He() const			{ return Conc_3He; }

  /// @brief Retrieves the concentration of fallout 7Be
  /// @return concentration of fallout 7Be
  double getConc_f7Be() const			{ return Conc_f7Be; }

  /// @brief Retrieves the concentration of fallout 10Be
  /// @return concentration of fallout 10Be
  double getConc_f10Be() const		{ return Conc_f10Be; }
	
  /// @brief Retrieves the concentration of 201Pb
  /// @return concentration of 210Pb
  double getConc_f210Pb() const		{ return Conc_f210Pb; }

  /// @brief Retrieves the concentration of 137Cs
  /// @return concentration of 137Cs
  double getConc_f137Cs() const		{ return Conc_f137Cs; }

  /// @brief Retrieves the effective depth
  /// @return the effective depth
  double geteffective_dLoc() const	{ return effective_dLoc; }

  /// @brief Retrieves the zeta location
  /// @return zetaLoc, which the the elevation relative to some arbitrary datum
  double get_zetaLoc() const			{ return zetaLoc; }

  /// @brief returns the mass of the particle
  /// @return the mass
  double getMass() const  				{ return Mass; }

  /// @brief returns the starting mass of the particle
  /// @return the starting mass
  double getStartingMass() const			{ return StartingMass; }

  /// @brief returns the surface area of the particle
  /// @return the surface area	
  double getSurfaceArea()	const			{ return SurfaceArea; }

  /// @brief returns the GSDType of the particle, which is an index into a 
  /// LSDVolumeParticleInfo object
  /// @return the surface GSDType	
  int getGSDType()	const				{ return GSDType; }

  /// @brief a utility function to convert between metres and g/cm^2
  /// @param l_in_m the length in metres
  /// @param rho the density in kg/m^3
  /// @author SMM
  /// @date 17/12/2014
  double convert_m_to_gpercm2(double l_in_m, double rho);

  /// @brief a utility function to convert between metres and g/cm^2
  /// @param l_in_gpercm2 the length in g per cm^2
  /// @param rho the density in kg/m^3
  /// @author SMM
  /// @date 17/12/2014
  double convert_gpercm2_to_m(double l_in_gpercm2, double rho);
  
  
  /// @brief update the 10Be concentration based on a constant erosion rate
  /// using the full production range, including muons. 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_10Be_conc(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 26Al concentration based on a constant erosion rate
  /// using the full production range, including muons. 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_26Al_conc(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 14C concentration based on a constant erosion rate
  /// using the full production range, including muons. 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_14C_conc(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 36Cl concentration based on a constant erosion rate
  /// using the full production range, including muons. 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_36Cl_conc(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 21Ne concentration based on a constant erosion rate
  /// using the full production range, including muons. 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_21Ne_conc(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 3He concentration based on a constant erosion rate
  /// using the full production range, including muons. 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_3He_conc(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief This updates 10Be nuclide concentration if erosion is increasing (or decreasing) linearly
  /// alpha is the change in erosion rate--do not set alpha to zero!
  /// @details This solves an analytical solution for cosmo concertration with a linear
  /// increase in cosmo concentration
  /// @param dt time step to calculate next cosmo concetration
  /// @param erosion_rate the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
  /// @param alpha the increase in erosion rate (in g/cm^2/yr^2)
  /// @param CRNp a CRN parameters object
  /// @author SMM
  /// @date 01/01/2010
  void update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, 
       LSDCRNParameters& CRNp);

  /// @brief update the 10Be concentration based on a constant erosion rate
  /// using the ONLY NEUTRON porduction. It is less computationally expensive
  /// than calcualting the full production and is a good approximation of the
  /// erosion rates are slow 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr      POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_10Be_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 26Al concentration based on a constant erosion rate
  /// using the ONLY NEUTRON porduction. It is less computationally expensive
  /// than calcualting the full production and is a good approximation of the
  /// erosion rates are slow 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr       POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_26Al_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 14C concentration based on a constant erosion rate
  /// using the ONLY NEUTRON porduction. It is less computationally expensive
  /// than calcualting the full production and is a good approximation of the
  /// erosion rates are slow 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_14C_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief update the 36Cl concentration based on a constant erosion rate
  /// using the ONLY NEUTRON porduction. It is less computationally expensive
  /// than calcualting the full production and is a good approximation of the
  /// erosion rates are slow 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
  /// @param dt the timestep over which erosion occurs
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_36Cl_conc_neutron_only(double dt,double erosion, LSDCRNParameters& CRNp);

  /// @brief This assigns nuclide concentrations with constant values
  /// @param C_10Be the 10Be concentration
  /// @param C_26Al the 26Al concentration
  /// @param C_36Cl the 36Cl concentration
  /// @param C_14C the 14C concentration
  /// @param C_21Ne the 21Ne concentration
  /// @param C_3He the 3He concentration
  void update_cosmo_conc_const(double C_10Be, double C_26Al, double C_36Cl,
                          double C_14C, double C_21Ne, double C_3He);

  /// @brief Bring the 10Be concentration to steady state based
  /// on a constant erosion rate using full muogenic production.  
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_10Be_SSfull(double erosion, LSDCRNParameters& CRNp);

  /// @brief Bring the 10Be concentration to steady state based
  ///  on a constant erosion rate using full muogenic production. This version
  ///  is for a DEPTH INTEGRATED concentration and as such should only be used
  ///  for basinwide calculations 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @param top_eff_depth the effective depth g/cm^2 at the top of the
  ///  section being depth-integrated
  /// @param bottom_eff_depth the effective depth g/cm^2 at the bottom of the
  ///  section being depth-integrated
  /// @author SMM
  /// @date 20/02/2015
void update_10Be_SSfull_depth_integrated(double erosion_rate, LSDCRNParameters& CRNp,
                                  double top_eff_depth, double bottom_eff_depth);


  /// @brief Bring the 26Al concentration to steady state based
  /// on a constant erosion rate using full muogenic production.  
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_26Al_SSfull(double erosion, LSDCRNParameters& CRNp);

  /// @brief Bring the 26Al concentration to steady state based
  ///  on a constant erosion rate using full muogenic production. This version
  ///  is for a DEPTH INTEGRATED concentration and as such should only be used
  ///  for basinwide calculations 
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @param top_eff_depth the effective depth g/cm^2 at the top of the
  ///  section being depth-integrated
  /// @param bottom_eff_depth the effective depth g/cm^2 at the bottom of the
  ///  section being depth-integrated
  /// @author SMM
  /// @date 20/02/2015
  void update_26Al_SSfull_depth_integrated(double erosion_rate, LSDCRNParameters& CRNp,
                                  double top_eff_depth, double bottom_eff_depth);


  /// @brief Bring the 14C concentration to steady state based
  /// on a constant erosion rate using full muogenic production.  
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_14C_SSfull(double erosion, LSDCRNParameters& CRNp);

  /// @brief Bring the 36Cl concentration to steady state based
  /// on a constant erosion rate using full muogenic production.  
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_36Cl_SSfull(double erosion, LSDCRNParameters& CRNp);

  /// @brief Bring the 21Ne concentration to steady state based
  /// on a constant erosion rate using full muogenic production.  
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_21Ne_SSfull(double erosion, LSDCRNParameters& CRNp);

  /// @brief Bring the 3He concentration to steady state based
  /// on a constant erosion rate using full muogenic production.  
  /// @details This function solves for the updated concentration assuming
  /// a constant erosion rate (in g/cm^2/yr). It is an analytical
  /// steady state solution
  /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// to approximate production from the different production mechanisms
  /// @author SMM
  /// @date 01/01/2010
  void update_3He_SSfull(double erosion, LSDCRNParameters& CRNp);

  /// @brief A wrapper function to update all the nuclide concentrations
  /// in one go. It uses full muogenic production
  /// @param dt the timestep over which the concetrations are updated
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @author SMM
  /// @date 01/01/2014
  void update_all_CRN(double dt, double erosion, LSDCRNParameters& CRNp);
  
  /// @brief A wrapper function to update all the nuclide concentrations
  /// in one go. It uses NEUTRON PRODUCTION ONLY to save computational
  /// expense. This is a reasonable approximation in slowly eroding landscapes
  /// @param dt the timestep over which the concetrations are updated
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @author SMM
  /// @date 01/01/2014
  void update_all_CRN_neutron_only(double dt, double erosion, LSDCRNParameters& CRNp);

  /// @brief A wrapper function to update all the nuclide concentrations
  /// to steady state using full muon production
  /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @author SMM
  /// @date 01/01/2014
  void update_all_CRN_SSfull(double erosion_rate, LSDCRNParameters& CRNp);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density above the sampling depth and for full muon 
  /// production based on COSMOCALC scaling
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @param muon_string a string containting the muon scaling
  ///  can be Schaller, Granger or Braucher. Default is Braucher
  /// @param top_eff_depth the effective depth g/cm^2 at the top of the
  ///  section being depth-integrated
  /// @param bottom_eff_depth the effective depth g/cm^2 at the bottom of the
  ///  section being depth-integrated
  /// @return the apparent erosion rate in g/cm^2/yr and m/yr
  /// @author SMM
  /// @date 13/03/2015 
  vector<double> apparent_erosion_10Be_COSMOCALC(double rho, LSDCRNParameters& CRNp,
                                    double scaling_factor, string Muon_scaling,
                                    double top_eff_depth, double bottom_eff_depth);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density above the sampling depth and if the
  /// production was assumed to be from neutrons only for 10Be
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @return the apparent erosion rate in m/yr
  /// @author SMM
  /// @date 01/01/2010 
  double apparent_erosion_10Be_neutron_only(double rho, LSDCRNParameters& CRNp);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density above the sampling depth and for full muon 
  /// production based on COSMOCALC scaling
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @param muon_string a string containting the muon scaling
  ///  can be Schaller, Granger or Braucher. Default is Braucher
  /// @param top_eff_depth the effective depth g/cm^2 at the top of the
  ///  section being depth-integrated
  /// @param bottom_eff_depth the effective depth g/cm^2 at the bottom of the
  ///  section being depth-integrated
  /// @return the apparent erosion rate in g/cm^2/yr and m/yr
  /// @author SMM
  /// @date 13/03/2015 
  vector<double> apparent_erosion_26Al_COSMOCALC(double rho, LSDCRNParameters& CRNp,
                                    double scaling_factor, string Muon_scaling,
                                    double top_eff_depth, double bottom_eff_depth);


  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density above the sampling depth and if the
  /// production was assumed to be from neutrons only for 26Al
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @return the apparent erosion rate in m/yr
  /// @author SMM
  /// @date 01/01/2010 
  double apparent_erosion_26Al_neutron_only(double rho, LSDCRNParameters& CRNp);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density above the sampling depth and if the
  /// production was assumed to be from neutrons only for 14C
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @return the apparent erosion rate in m/yr
  /// @author SMM
  /// @date 01/01/2010 
  double apparent_erosion_14C_neutron_only(double rho, LSDCRNParameters& CRNp);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density abouve the sampling depth and if the
  /// production was assumed to be from neutrons only
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @return the apparent erosion rate in m/yr
  /// @author SMM
  /// @date 01/01/2010 
  double apparent_erosion_36Cl_neutron_only(double rho, LSDCRNParameters& CRNp);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density above the sampling depth and if the
  /// production was assumed to be from neutrons only for 21Ne
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @return the apparent erosion rate in m/yr
  /// @author SMM
  /// @date 01/01/2010 
  double apparent_erosion_21Ne(double rho, LSDCRNParameters& CRNp);

  /// @brief This returns the 'apparent' erosion rate that one would calcualte
  /// based on an assumed density aboue the sampling depth and if the
  /// production was assumed to be from neutrons only for 3He
  /// @param rho the density in kg/m^3 above the 'sampling' point
  /// @param CRNp a CRN parameters object that stores the coefficients
  /// @return the apparent erosion rate in m/yr
  /// @author SMM
  /// @date 01/01/2010 
  double apparent_erosion_3He(double rho, LSDCRNParameters& CRNp);

  /// @brief  functions for dealing with fallout nuclides
  void update_fallout10Be_simple_density(double dt, double M_supply_surface,
      double rho_s, double k_f10Be, double deltad, LSDCRNParameters& CRNp);


  void update_fallout10Be_simple_density_2exp(double dt, double M_supply_surface,
    double rho_skg, double k1_f10Be, double k2_f10Be, double chi_f10Be,
    double deltad_m, LSDCRNParameters& CRNp);


  /// @brief This assignes depths. Note: it does not check if the
  /// effective depth is compatible with the depth!
  /// @param delta_d the new depth (in metres)
  /// @param delta_ed the new effective depth in g/cm^2
  void update_depths(double delta_d, double delta_ed);

  /// @brief Updates the depth and calculates and updates the effective depth
  ///  for a one layer density model (that is, all material has the same density)
  /// @param d the depth of the particle in metres
  /// @param rho the density of the material in kg/m^3
  /// @author SMM
  /// @date 16/12/2014
  void calculate_effective_depth_one_layer(double d, double rho);


  /// @brief This assigns a new value for zeta
  /// @param new_zeta the new zeta, which is elevation above an arbitrary datum
  void update_zetaLoc(double new_zeta);

  void update_zetaLoc_with_new_surface(double new_zeta);

  void erode_mass_only(double dt, double mass_erosion_rate);
  
  void erode_mass_only_linear_increase(double dt, double mass_erosion_rate, double alpha);

  //##################
  // !!!!!! SCALING
  //##################

  /// @brief This integrates over beta for Deselits scaling (see the paper)
  /// Modified from Greg Balco's code:
  ///  http://hess.ess.washington.edu/math
  /// @param x integrating variable
  /// @param Rc is cutoff rigidity (GV)
  ///  @return the scaling factor
  /// @author SMM
  /// @date 1/12/2014
  double intOfBeta(double x, double Rc);

  /// @brief This gets the deselits 2006 scaling
  /// Modified from Greg Balco's code:
  ///  http://hess.ess.washington.edu/math
  /// @param h is atmospheric pressure (hPa)
  /// @param Rc is cutoff rigidity (GV)
  ///  @return the scaling factor
  /// @author SMM
  /// @date 1/12/2014
  double desilets2006sp(double h,double Rc);  

  /// @brief This gets the Dunai Scaling
  /// Modified from Greg Balco's code:
  ///  http://hess.ess.washington.edu/math
  /// @param h is atmospheric pressure (hPa)
  /// @param Rc is cutoff rigidity (GV)
  ///  @return the scaling factor
  /// @author SMM
  /// @date 1/12/2014
  double dunai2001sp(double h, double Rc);     

  /// @brief This gets the Stone scaling
  /// Modified from Greg Balco's code:
  ///  http://hess.ess.washington.edu/math
  /// @param h is atmospheric pressure (hPa)
  /// @param Rc is cutoff rigidity (GV)
  /// @param S solar modulation factor (nondimensional, see source paper)
  ///  @return the scaling factor
  /// @author SMM
  /// @date 2/12/2014
  double lifton2006sp(double h, double Rc, double S);     

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

  /// @brief This implements a cutoff-rigidity rather than latitude based scaling
  /// scheme based on the Lal spallation polynomials. For use in
  /// paleomagnetically-corrected exposure-age calculations. 
  /// @detail Original version Written by Greg Balco -- UW Cosmogenic Nuclide Lab
  ///  balcs@u.washington.edu
  ///  March, 2007
  ///  Part of the CRONUS-Earth online calculators: 
  ///  http://hess.ess.washington.edu/math
  /// @param h = tmospheric pressure (hPa)
  /// @param Rc cutoff rigidity (GV)
  /// @return scaling factor
  /// @author SMM
  /// @date 5/12/2014
  double stone2000Rcsp(double h, double Rc);

  /// @brief this calculates a thickness scaling factor that is used in the
  ///  CRONUS calculator calculations
  /// @param LSDCRNP an LSDCRNParameters object
  /// @param use_CRONUS a boolean that if true uses the attenuation depth
  ///  of the CRONUS calculator (at the moment this makes no difference)
  /// @return The thickness scaling factor (between 0 and 1)
  /// @author SMM
  /// @date 14/12/2014
  double thickness_scaling_factor(LSDCRNParameters& LSDCRNP, bool use_CRONUS);

  /// @brief This function is the initial guess of erosion rate
  ///  that repicates the initial guesss component of the CRONUS 
  ///  erosion rate function
  /// @detail This replicates the 'initial guess' stage in Greg Balco's 
  ///  CRONUS calculator code get_al_be_erosion.m
  /// @param LSDCRNP and LSDCRNParameters object
  /// @param pressure the atmospheric pressure in hPa
  /// @param lat the latitude  
  /// @param N_10Be the number of 10Be atoms
  /// @param N_26 the number of 26Al atoms
  /// @param topo_scale the topographic scaling (between 0 and 1)
  /// @param snow_scale the snow scaling (between 0 and 1)
  /// @return Erosion_rates a vector<double> Eguess
  ///  Eguess[0] = 10Be
  ///  Eguess[1] = 26Al
  /// @author SMM
  /// @date 14/12/2014
  vector<double> CRONUS_initial_guess(LSDCRNParameters& LSDCRNP, double pressure,
                                 double lat, double N_10Be, double N_26Al, 
                                 double topo_scale, double snow_scale);

  /// @brief This function wraps the functions for getting the erosion rate 
  ///  from Al and Be data. 
  /// @detail The function emulates the get_al_be_erosion.m from the
  ///  CRONUS calculator, written by Greg Balco
  /// @param LSDCRNP and LSDCRNParameters object
  /// @param pressure the atmospheric pressure in hPa
  /// @param lat the latitude  
  /// @param rho the density in kg/m^3
  /// @param N_10Be the number of 10Be atoms
  /// @param N_26 the number of 26Al atoms
  /// @param sample_del10 The accelerator uncertanty in the number of 10Be atoms
  /// @param sample_del26 The accelerator uncertanty in the number of 26Al atoms
  /// @param topo_scale the topographic scaling (between 0 and 1)
  /// @param snow_scale the snow scaling (between 0 and 1)
  /// @return erate_consts A vector containing information about erosion rates and 
  ///   uncertainties
  ///  erate_consts[0] = erate_10Be in g/cm^2/yr
  ///  erate_consts[1] = erate_10Be in M/Myrs
  ///  erate_consts[2] = erate 10Be internal uncertainty in M/Myrs
  ///  erate_consts[3] = erate 10Be external uncertainty in M/Myrs
  ///  erate_consts[4] = 10Be muon production rate atoms/g/yr
  ///  erate_consts[5] = 10Be spallation production rate atoms/g/yr
  ///  erate_consts[6] = erate_26Al in g/cm^2/yr
  ///  erate_consts[7] = erate_26Al in M/Myrs
  ///  erate_consts[9] = erate 26Al internal uncertainty in M/Myrs
  ///  erate_consts[9] = erate 26Al external uncertainty in M/Myrs
  ///  erate_consts[10] = 26Al muon production rate atoms/g/yr
  ///  erate_consts[11] = 26Al spallation production rate atoms/g/yr  
  /// @author SMM
  /// @date 16/12/2014
  vector<double> CRONUS_get_Al_Be_erosion(LSDCRNParameters& LSDCRNP, double pressure,
                      double lat, double rho, double N_10Be, double N_26Al,
                      double sample_del10, double sample_del26,
                      double topo_scale, double snow_scale);

  /// @brief This function wraps the functions for getting the erosion rate 
  ///  from Al and Be data. Similar to above but modifies production rates 
  /// @detail The function emulates the get_al_be_erosion.m from the
  ///  CRONUS calculator, written by Greg Balco
  /// @param LSDCRNP and LSDCRNParameters object
  /// @param pressure the atmospheric pressure in hPa
  /// @param lat the latitude  
  /// @param rho the density in kg/m^3
  /// @param N_10Be the number of 10Be atoms
  /// @param N_26 the number of 26Al atoms
  /// @param sample_del10 The accelerator uncertanty in the number of 10Be atoms
  /// @param sample_del26 The accelerator uncertanty in the number of 26Al atoms
  /// @param topo_scale the topographic scaling (between 0 and 1)
  /// @param snow_scale the snow scaling (between 0 and 1)
  /// @param Spallation_frac fraction change in the default spallation
  /// @param P_mu_frac fraction change in the default muon production
  /// @return erate_consts A vector containing information about erosion rates and 
  ///   uncertainties
  ///  erate_consts[0] = erate_10Be in g/cm^2/yr
  ///  erate_consts[1] = erate_10Be in M/Myrs
  ///  erate_consts[2] = erate 10Be internal uncertainty in M/Myrs
  ///  erate_consts[3] = erate 10Be external uncertainty in M/Myrs
  ///  erate_consts[4] = 10Be muon production rate atoms/g/yr
  ///  erate_consts[5] = 10Be spallation production rate atoms/g/yr
  ///  erate_consts[6] = erate_26Al in g/cm^2/yr
  ///  erate_consts[7] = erate_26Al in M/Myrs
  ///  erate_consts[9] = erate 26Al internal uncertainty in M/Myrs
  ///  erate_consts[9] = erate 26Al external uncertainty in M/Myrs
  ///  erate_consts[10] = 26Al muon production rate atoms/g/yr
  ///  erate_consts[11] = 26Al spallation production rate atoms/g/yr  
  /// @author SMM
  /// @date 16/12/2014
  vector<double> CRONUS_get_Al_Be_erosion_modified_production(LSDCRNParameters& LSDCRNP, double pressure,
                      double lat, double rho, double N_10Be, double N_26Al,
                      double sample_del10, double sample_del26,
                      double topo_scale, double snow_scale,
                      double Spallation_frac, double P_mu_frac);


  /// @brief This produces a screen report that replicates the CRONUS calculator
  ///  screen report
  /// @param erate_consts the vector that is returned by the function
  ///  LSDCRNParticle::CRONUS_get_Al_Be_erosion
  /// @author SMM
  /// @date 19/12/2014
  void CRONUS_screen_report(vector<double> erate_consts);


  /// @brief This function returns the number of atoms given an effective erosion rate
  ///   (this is the erosion rate in g/cm^2/yr)
  ///   It is used in an optimisation loop (equivalent to fzero in matlab)
  /// @param effective_erosion_rate erosion rate in g/cm^2/yr
  /// @param LSDCRNP and LSDCRNParameters object
  /// @param z_mu a vector of depths for muon flux. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_mu_z_10Be a vector of production of muons for 10Be. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_mu_z_26Al a vector of production of muons for 26Al. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param thickSF the thickness scaling factor, generated by 
  ///   LSDCRNParticle::thickness_scaling_factor
  /// @param P_sp_10Be an adjusted production rate that takes into account
  ///   scaling, snow and topo shielding (stoneP*toposhield*snowshield)
  /// @param P_sp_26Al an adjusted production rate that takes into account
  ///   scaling, snow and topo shielding (stoneP*toposhield*snowshield)
  /// @param N_10Be the number of 10Be atoms. Replaced in this function.
  /// @param N_26 the number of 26Al atoms. Replaced in this function.
  /// @author SMM
  /// @date 17/12/2014
  void CRONUS_calculate_N_forward(double effective_erosion_rate,
                            LSDCRNParameters& LSDCRNP,
                            vector<double>& z_mu, vector<double>& P_mu_z_10Be, 
                            vector<double>& P_mu_z_26Al, double thickSF, 
                            double P_sp_10Be, double P_sp_26Al,
                            double& N_Be10, double& N_Al26);

  /// @brief This function returns the number of atoms given an effective erosion rate
  ///   (this is the erosion rate in g/cm^2/yr)
  ///   It is used in an optimisation loop (equivalent to fzero in matlab)
  /// @detail this is an overloaded version that returns the atoms from the
  ///   two production mechanisms
  /// @param effective_erosion_rate erosion rate in g/cm^2/yr
  /// @param LSDCRNP and LSDCRNParameters object
  /// @param z_mu a vector of depths for muon flux. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_mu_z_10Be a vector of production of muons for 10Be. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_mu_z_26Al a vector of production of muons for 26Al. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param thickSF the thickness scaling factor, generated by 
  ///   LSDCRNParticle::thickness_scaling_factor
  /// @param P_sp_10Be an adjusted production rate that takes into account
  ///   scaling, snow and topo shielding (stoneP*toposhield*snowshield)
  /// @param P_sp_26Al an adjusted production rate that takes into account
  ///   scaling, snow and topo shielding (stoneP*toposhield*snowshield)
  /// @param N_10Be the number of 10Be atoms. Replaced in this function.
  /// @param N_26 the number of 26Al atoms. Replaced in this function.
  /// @param Be10_mu_N atoms of 10Be from muons
  /// @param Al26_mu_N atoms of 26Al from muons
  /// @param Be10_sp_N atoms of 10Be from spallation
  /// @param Al26_sp_N atoms of 26Al from spallation
  /// @author SMM
  /// @date 17/12/2014
  void CRONUS_calculate_N_forward(double effective_erosion_rate, 
                            LSDCRNParameters& LSDCRNP,
                            vector<double>& z_mu, vector<double>& P_mu_z_10Be, 
                            vector<double>& P_mu_z_26Al, double thickSF, 
                            double P_sp_10Be, double P_sp_26Al,
                            double& N_Be10, double& N_Al26, 
                            double& Be10_mu_N, double& Al26_mu_N,
                            double& Be10_sp_N, double& Al26_sp_N);

  /// @brief This function returns the number of atoms given an effective erosion rate
  ///   (this is the erosion rate in g/cm^2/yr)
  ///   It is used in an optimisation loop (equivalent to fzero in matlab)
  /// @param pressure atmospgeric pressure in HPa
  /// @param LSDCRNP and LSDCRNParameters object
  /// @param thickSF the thickness scaling factor, get from
  ///   LSDCRNParticle::thickness_scaling_factor
  /// @param rho the sample density in kg/m^3
  /// @param sample_N10 The number of 10Be atoms in the sample
  /// @param sample_N26 The number of 26Al atoms in the sample
  /// @param sample_del10 The accelerator uncertanty in the number of 10Be atoms
  /// @param sample_del26 The accelerator uncertanty in the number of 26Al atoms
  /// @param z_mu a vector of depths for muon flux. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_mu_z_10Be a vector of production of muons for 10Be. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_mu_z_26Al a vector of production of muons for 26Al. Generated by 
  ///   LSDCRNParameters::get_CRONUS_P_mu_vectors
  /// @param P_sp_10Be an adjusted production rate that takes into account
  ///   scaling, snow and topo shielding (stoneP*toposhield*snowshield)
  /// @param P_sp_26Al an adjusted production rate that takes into account
  ///   scaling, snow and topo shielding (stoneP*toposhield*snowshield)  
  /// @param eff_e_10 effective erosion rate from 10Be in g/cm^2/yr
  /// @param eff_e_26 effective erosion rate from 26Al in g/cm^2/yr
  /// @return errors a vector<double> of the errors in the erosion rate
  ///   errors[0] = external 10Be error
  ///   errors[1] = internal 10Be error (i.e., from the ams)
  ///   errors[2] = Production from spallation at surface for 10Be in atoms/g/yr
  ///   errors[3] = Production from muons at surface for 10Be in atoms/g/yr
  ///   errors[4] = external 10Be error
  ///   errors[5] = internal 10Be error (i.e., from the ams)
  ///   errors[6] = Production from spallation at surface for 26Al in atoms/g/yr
  ///   errors[7] = Production from muons at surface for 26Al in atoms/g/yr
  ///   errors[8] = Atoms from spallation of 10Be in atoms/g
  ///   errors[9] = Atoms from muons of 10Be in atoms/g
  ///   errors[10] = Atoms from spallation of 26Al in atoms/g
  ///   errors[11] = Atoms from muons of 26Al in atoms/g
  /// @author SMM
  /// @date 17/12/2014
  vector<double> CRONUS_error_propagation(double pressure, 
                            LSDCRNParameters& LSDCRNP, double thickSF, double rho,
                            double sample_N10, double sample_N26,
                            double sample_del10, double sample_del26,
                            vector<double>& z_mu, vector<double>& P_mu_z_10Be, 
                            vector<double>& P_mu_z_26Al, 
                            double P_sp_10Be, double P_sp_26Al,
                            double eff_e_10, double eff_e_26);

  /// @brief This finds the roots of a simple version of the number of atoms
  ///  given an effective erosion rate
  /// @detail This function is used by the CRONUS calculator to determine
  ///  uncertanties. 
  /// @param intial_erate the first guess of erosion rate in g/cm^2/yr.
  /// @param target the target number of atoms
  /// @param Psp production from spallation
  /// @param Pmu production from muons
  /// @param GammaSp effective attenuation length for spallation (g/cm^2)
  /// @param GammaMu effective attenuation length for muons (g/cm^2)
  /// @param GammaMu effective attenuation length for muons (g/cm^2)
  /// @param decay The decay coefficient of this nuclide
  /// @return the number of atoms per gram of the nuclide
  /// @author SMM
  /// @date 18/12/2014
  double CRONUS_simple_N_findroots(double inital_erate, double target,
                                           double Psp, double Pmu, double GammaSp,  
                                           double GammaMu,double decay);

  /// @brief This calculates a simple version of the number of atoms
  ///  given an effective erosion rate
  /// @detail This function is used by the CRONUS calculator to determine
  ///  uncertanties. 
  /// @param eff_eros erosion rate in g/cm^2/yr.
  /// @param Psp production from spallation
  /// @param Pmu production from muons
  /// @param GammaSp effective attenuation length for spallation (g/cm^2)
  /// @param GammaMu effective attenuation length for muons (g/cm^2)
  /// @param GammaMu effective attenuation length for muons (g/cm^2)
  /// @param decay The decay coefficient of this nuclide
  /// @return the number of atoms per gram of the nuclide
  /// @author SMM
  /// @date 18/12/2014
  double CRONUS_simple_N_forward(double eff_eros, double Psp, double Pmu,
                                double GammaSp, double GammaMu, double decay);

  protected:

  /// The effective depth of the particle in g/cm^2
  double effective_dLoc;

  /// the elevation (in metres) relative to an arbitrary datum (used to interface with 
  /// landscape evolution)
  double zetaLoc;

  /// concentration of 10Be in atoms/g
  double Conc_10Be;

  /// concentration of 26Al in atoms/g	
  double Conc_26Al;      // a/g

  /// concentration of 36Cl in atoms/g
  double Conc_36Cl;      // a/g

  /// concentration of 14C in atoms/g
  double Conc_14C;      // a/g

  /// concentration of 21Ne in atoms/g
  double Conc_21Ne;      // a/g

  /// concentration of 3He in atoms/g
  double Conc_3He;			// a/g

  /// concentration of fallout 7Be in units tba
  double Conc_f7Be;			// fallout units tba

  /// concentration of fallout 10Be in units a/g
  double Conc_f10Be; 			// fallout, units a/g

  /// concentration of fallout 210Pb in units tba
  double Conc_f210Pb;			// fallout, units tba

  /// concentration of fallout 137Cs in units tba
  double Conc_f137Cs;			// fallout, units tba

  /// the mass of the particle. Default units kg
  double Mass;
  
  /// Starting mass of the particle. Default units kg
  double StartingMass;	
  
  /// Surface area in m^2
  double SurfaceArea;
  
  /// an integer used to denote the GSD type
  /// this is an INDEX into a vector in the LSDVolumeParticleInfo
  /// object  
  int GSDType;
  
  private:
  /// @breif default create function
  void create();
  
  /// @brief this just sets simple location information
  void create(int startType, double startxLoc,
            double startzLoc);
  
  /// @brief a create function setting location data
  void create(int startType, double startxLoc,
             double startdLoc, double start_effdloc,
	            double startzloc);
  
  /// @brief a create function for a setting location data with y loc
  void create(int startType, double startxLoc, double startyLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc);	    
  
  /// @brief a create function for a volume particle
  void create(int startType, int startGSDType,
				double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc, double startMass,
	            double startSurfaceArea);
	            
	/// @brief a create function for a volume particle  with y loc          
  void create(int startType, int startGSDType, double startxLoc, double startyLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc, double startMass,
	            double startSurfaceArea);	 
              
  /// @brief a create function if you just want to work with cosmo                                 
  void create(int startType, double startxLoc,double startzeta_Loc,
			double startdLoc, double start_effdLoc,
			double start_C10Be,double start_C14C);
	
  /// @brief a create function if you just want to work with cosmo, but with more nuclides
  /// than previous version		
  void create(int startType, double startxLoc,double startzeta_Loc,
			double startdLoc, double start_effdLoc,
			double start_C10Be,double start_C14C, double start_21Ne);
	
  /// @brief this is the create function used in the copy constructors 		
  void create(int startType, int startGSDType, int startCellIndex, 
          double startAge, double startOSLAge,
          double startxLoc,double startyLoc,double startdLoc, double startefdLoc,
          double startzLoc, double start_C10Be, double start_C26Al,
          double start_C36Cl, double start_C14C,
          double start_C21Ne, double start_C3He,
          double start_Cf7Be, double start_Cf10Be,
          double start_Cf210Pb, double start_Cf137Cs, 
          double start_Mass, double start_StartingMass,double startSurfaceArea);
          



   

};


#endif
